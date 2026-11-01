"""Calibrating predictive claims for expression scores in randomized trials.

A score predicts benefit from a treatment only if its association with outcome differs between the arms of a
randomized comparison: the score x arm interaction. This module estimates that interaction (Cox or logistic model),
pools it over trials, and calibrates it against two references:

- size-matched random gene sets (is the interaction specific to the score's genes?), and
- single genes, by a paired within-arm bootstrap (does the score do better than one gene?),

with the difference detectable at the trials' precision (``mdd``). Interactions are benefit-oriented: positive means
more benefit from the experimental arm per SD of the score (-log HR for survival, log OR for response).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import norm, t as tdist


def standardize(x) -> np.ndarray:
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std(ddof=1)


# ------------------------------------------------------------------ models
def cox_ph(X, time, event, max_iter: int = 60, tol: float = 1e-10) -> Tuple[np.ndarray, np.ndarray]:
    """Cox proportional-hazards fit (Efron ties) by Newton-Raphson; returns coefficients and standard errors."""
    X = np.asarray(X, float); X = X[:, None] if X.ndim == 1 else X
    time = np.asarray(time, float); event = np.asarray(event).astype(bool)
    o = np.argsort(-time, kind="stable"); X, time, event = X[o], time[o], event[o]
    ut = np.unique(time[event])[::-1]                                       # distinct event times, descending
    end = np.searchsorted(-time, -ut, side="right")                         # risk set = rows [0, end)
    grp = [np.where(event & (time == tj))[0] for tj in ut]

    def parts(b):
        eta = X @ b; c = eta.max(); w = np.exp(eta - c); p = X.shape[1]
        S0 = np.cumsum(w); S1 = np.cumsum(w[:, None] * X, 0); S2 = np.cumsum(w[:, None, None] * X[:, :, None] * X[:, None, :], 0)
        ll, g, H = 0.0, np.zeros(p), np.zeros((p, p))
        for e_, D in zip(end, grp):
            d = len(D); s0, s1, s2 = S0[e_ - 1], S1[e_ - 1], S2[e_ - 1]
            t0 = w[D].sum(); t1 = (w[D, None] * X[D]).sum(0); t2 = (w[D, None, None] * X[D, :, None] * X[D, None, :]).sum(0)
            ll += (eta[D] - c).sum(); g += X[D].sum(0)
            for l in range(d):
                f = l / d; A = s0 - f * t0; B = s1 - f * t1; C = s2 - f * t2
                ll -= np.log(A); g -= B / A; H += C / A - np.outer(B, B) / A ** 2
        return ll, g, H

    b = np.zeros(X.shape[1]); ll, g, H = parts(b)
    for _ in range(max_iter):
        step = np.linalg.solve(H, g); s = 1.0
        while True:
            ll2, g2, H2 = parts(b + s * step)
            if ll2 >= ll - 1e-12 or s < 1e-8:
                break
            s /= 2
        b = b + s * step; ll, g, H = ll2, g2, H2
        if np.max(np.abs(s * step)) < tol:
            break
    return b, np.sqrt(np.diag(np.linalg.inv(H)))


def _logit_ml(X, y, max_iter=100, tol=1e-10):
    b = np.zeros(X.shape[1])
    for _ in range(max_iter):
        p = 1 / (1 + np.exp(-np.clip(X @ b, -30, 30))); I = X.T @ (X * (p * (1 - p))[:, None])
        step = np.linalg.solve(I, X.T @ (y - p)); b = b + step
        if np.max(np.abs(step)) < tol:
            p = 1 / (1 + np.exp(-X @ b)); return b, np.linalg.inv(X.T @ (X * (p * (1 - p))[:, None])), True
    return b, None, False


def _logit_firth(X, y, max_iter=200, tol=1e-9):
    b = np.zeros(X.shape[1])

    def pll(b):
        p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
        return np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)) + 0.5 * np.linalg.slogdet(X.T @ (X * W[:, None]))[1]
    for _ in range(max_iter):
        p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p); I = X.T @ (X * W[:, None]); Ii = np.linalg.inv(I)
        h = np.einsum("ij,jk,ik->i", X * np.sqrt(W)[:, None], Ii, X * np.sqrt(W)[:, None])
        step = np.linalg.solve(I, X.T @ (y - p + h * (0.5 - p))); f0, s = pll(b), 1.0
        while s > 1e-6 and pll(b + s * step) < f0 - 1e-12:
            s /= 2
        b = b + s * step
        if np.max(np.abs(s * step)) < tol:
            break
    p = 1 / (1 + np.exp(-X @ b)); return b, np.linalg.inv(X.T @ (X * (p * (1 - p))[:, None]))


def logistic(X, y, firth: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Logistic regression (X includes the intercept column); Firth's penalized fit when ML separates."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    try:
        b, V, ok = _logit_ml(X, y)
    except np.linalg.LinAlgError:
        b, V, ok = None, None, False
    sep = (not ok) or np.any(np.abs(b) > 10) or np.any(~np.isfinite(np.diag(V))) or np.any(np.sqrt(np.diag(V)) > 10)
    if sep:
        if not firth:
            return np.full(X.shape[1], np.nan), np.full(X.shape[1], np.nan)
        b, V = _logit_firth(X, y)
    return b, np.sqrt(np.diag(V))


@dataclass
class Trial:
    """One randomized comparison. ``Z``: genes x samples, z-scored across the samples of this comparison (both arms).
    ``arm``: 1 = experimental (checkpoint-containing), 0 = control. Survival: ``time`` and ``event``; binary: ``response``
    (1 = benefit, e.g. pCR). ``covariates``: samples x k array kept in the logistic model (e.g. stratification factors)."""
    name: str
    Z: pd.DataFrame
    arm: np.ndarray
    time: Optional[np.ndarray] = None
    event: Optional[np.ndarray] = None
    response: Optional[np.ndarray] = None
    covariates: Optional[np.ndarray] = None
    firth: bool = True
    kind: str = field(init=False)

    def __post_init__(self):
        self.arm = np.asarray(self.arm, float)
        self.kind = "cox" if self.time is not None else "logit"
        if self.kind == "logit" and self.response is None:
            raise ValueError(f"{self.name}: give time and event, or response")

    def take(self, idx) -> "Trial":
        g = lambda v: None if v is None else np.asarray(v)[idx]
        return Trial(self.name, self.Z.iloc[:, idx], self.arm[idx], g(self.time), g(self.event), g(self.response),
                     g(self.covariates), self.firth)


def interaction(trial: Trial, score) -> Tuple[float, float]:
    """Benefit-oriented score x arm interaction per SD of the score, and its standard error."""
    s = standardize(score); a = trial.arm
    if trial.kind == "cox":
        b, se = cox_ph(np.column_stack([s, a, s * a]), trial.time, trial.event); return -b[2], se[2]
    cov = [] if trial.covariates is None else [np.asarray(trial.covariates, float).reshape(len(s), -1)]
    b, se = logistic(np.column_stack([np.ones_like(s), s, a, s * a] + cov), trial.response, trial.firth)
    return b[3], se[3]


# ------------------------------------------------------------------ pooling
def _reml_tau2(y, v):
    def nll(t2):
        w = 1 / (v + t2); mu = np.sum(w * y) / w.sum()
        return 0.5 * (np.sum(np.log(v + t2)) + np.log(w.sum()) + np.sum(w * (y - mu) ** 2))
    hi = max(10 * np.var(y, ddof=1), 10 * v.max(), 1e-4)
    r = minimize_scalar(nll, bounds=(0, hi), method="bounded", options=dict(xatol=1e-12))
    return r.x if nll(r.x) < nll(0.0) else 0.0


def pool(estimates: Sequence[float], ses: Sequence[float]) -> Dict[str, float]:
    """Random-effects pooling, conservative rule: the wider of the Hartung-Knapp (REML) and DerSimonian-Laird
    intervals and the larger of their P values. Non-finite estimates are dropped."""
    y, se = np.asarray(estimates, float), np.asarray(ses, float)
    ok = np.isfinite(y) & np.isfinite(se) & (se > 0); y, se = y[ok], se[ok]; k = len(y)
    if k < 2:
        return dict(est=np.nan, lo=np.nan, hi=np.nan, p=np.nan, k=k)
    v = se ** 2; t2 = _reml_tau2(y, v); w = 1 / (v + t2); mu = np.sum(w * y) / w.sum()
    se_hk = np.sqrt(np.sum(w * (y - mu) ** 2) / ((k - 1) * w.sum())); q = tdist.ppf(0.975, k - 1)
    hk = (mu - q * se_hk, mu + q * se_hk, 2 * tdist.sf(abs(mu / se_hk), k - 1))
    w0 = 1 / v; m0 = np.sum(w0 * y) / w0.sum(); Q = np.sum(w0 * (y - m0) ** 2)
    t2d = max(0.0, (Q - (k - 1)) / (w0.sum() - np.sum(w0 ** 2) / w0.sum())); wd = 1 / (v + t2d)
    md = np.sum(wd * y) / wd.sum(); sd = np.sqrt(1 / wd.sum()); z = norm.ppf(0.975)
    dl = (md - z * sd, md + z * sd, 2 * norm.sf(abs(md / sd)))
    lo, hi = hk[:2] if (hk[1] - hk[0]) >= (dl[1] - dl[0]) else dl[:2]
    return dict(est=mu, lo=lo, hi=hi, p=max(hk[2], dl[2]), k=k)


def mdd(se: float, power: float = 0.8, alpha: float = 0.05) -> float:
    """Smallest difference detectable with the given power by a two-sided test at the given standard error."""
    return (norm.ppf(1 - alpha / 2) + norm.ppf(power)) * se


def events_needed(delta: float, rho: float, allocation: float = 0.5, power: float = 0.8, alpha: float = 0.05,
                  inflation: float = 1.0) -> float:
    """Events a randomized survival comparison needs to detect a difference ``delta`` (per SD) between the treatment
    interactions of two standardized scores correlated at ``rho``, by normal theory:
    var(difference) = 2 (1 - rho) / (D p (1 - p)), with p the fraction allocated to the experimental arm.
    ``inflation`` multiplies the result (in the manuscript's six comparisons the bootstrap variance was 1.17 times the
    theoretical value in the survival comparisons)."""
    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    return inflation * 2 * (1 - rho) * z ** 2 / (delta ** 2 * allocation * (1 - allocation))


# ------------------------------------------------------------------ references
def weighted_score(Z: pd.DataFrame, weights: Dict[str, float], min_genes: int = 3, min_frac: float = 0.5):
    """sum(w z) / sum(|w|) over the genes present (finite, non-constant); None when fewer than ``min_genes`` or than
    ``min_frac`` of the weighted genes are present. Returns (score, genes used, weights used)."""
    sd = Z.std(axis=1)
    genes = [g for g in weights if g in Z.index and np.isfinite(Z.loc[g].to_numpy(float)).all() and sd[g] > 0]
    if len(genes) < min_genes or len(genes) < min_frac * len(weights):
        return None, genes, None
    w = np.array([weights[g] for g in genes], float)
    return Z.loc[genes].to_numpy(float).T @ w / np.abs(w).sum(), genes, w


def random_reference(trials: Sequence[Trial], sizes_weights: Dict[str, np.ndarray], n_draws: int = 1000,
                     seed: int = 0) -> np.ndarray:
    """Pooled interactions of size-matched random gene sets: in each trial, draw b takes as many genes as the score
    uses there (from that trial's complete, non-constant genes) with the score's weights in random order; the draws
    are pooled trial by trial with ``pool``. Returns the n_draws pooled estimates."""
    per = []
    for i, tr in enumerate(trials):
        w = np.asarray(sizes_weights[tr.name], float); n = len(w)
        Z = tr.Z[(tr.Z.std(axis=1) > 0) & np.isfinite(tr.Z).all(axis=1)].to_numpy(float)
        rng = np.random.default_rng([seed, i, n]); out = []
        for _ in range(n_draws):
            idx = rng.choice(Z.shape[0], n, replace=False); ww = w[rng.permutation(n)]
            try:
                out.append(interaction(tr, Z[idx].T @ ww / np.abs(ww).sum()))
            except Exception:
                out.append((np.nan, np.nan))
        per.append(np.array(out))
    D = np.stack(per, 1)
    return np.array([pool(D[b, :, 0], D[b, :, 1])["est"] for b in range(n_draws)])


def bootstrap_indices(trial: Trial, n_boot: int = 300, seed: int = 0) -> List[np.ndarray]:
    """Within-arm patient resamples (each arm resampled with replacement at its own size)."""
    rng = np.random.default_rng(seed); i1, i0 = np.where(trial.arm == 1)[0], np.where(trial.arm == 0)[0]
    return [np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]) for _ in range(n_boot)]


def bootstrap_interactions(trial: Trial, score, idx: List[np.ndarray]) -> np.ndarray:
    score = np.asarray(score, float); out = []
    for ix in idx:
        try:
            out.append(interaction(trial.take(ix), score[ix])[0])
        except Exception:
            out.append(np.nan)
    return np.array(out)


def paired_difference(diffs: Sequence[float], boot_a: Sequence[np.ndarray], boot_b: Sequence[np.ndarray]) -> Dict[str, float]:
    """Inverse-variance pooled difference between two scores' interactions over trials; each trial's variance is
    that of the bootstrap differences on common resamples."""
    v = []
    for a, b in zip(boot_a, boot_b):
        db = np.asarray(a) - np.asarray(b); v.append(np.var(db[np.isfinite(db)], ddof=1))
    d, v = np.asarray(diffs, float), np.asarray(v, float); w = 1 / v
    D = np.sum(w * d) / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(diff=D, se=se, lo=D - 1.96 * se, hi=D + 1.96 * se, p=2 * norm.sf(abs(D / se)), mdd=mdd(se), k=len(d))


# ------------------------------------------------------------------ the whole calibration
def calibrate(trials: Sequence[Trial], scores: Dict[str, Dict[str, float]], genes: Sequence[str],
              n_random: int = 1000, n_boot: int = 300, seed: int = 0) -> Dict[str, pd.DataFrame]:
    """Calibrate gene-weighted scores (``scores``: name -> {gene: weight}) over randomized ``trials``.

    Returns DataFrames: ``per_trial`` (interaction per score, gene and trial), ``pooled`` (pooled interaction, random
    reference median and 95th percentile, empirical P, whether the score exceeds the reference) and ``paired``
    (difference of each score from each single gene, with its P and the detectable difference)."""
    S, used = {}, {}
    for nm, wts in scores.items():
        for tr in trials:
            sc, g, w = weighted_score(tr.Z, wts)
            if sc is not None:
                S[(nm, tr.name)], used[(nm, tr.name)] = sc, w
    for g in genes:
        for tr in trials:
            if g in tr.Z.index and np.isfinite(tr.Z.loc[g].to_numpy(float)).all() and tr.Z.loc[g].std() > 0:
                S[(g, tr.name)] = tr.Z.loc[g].to_numpy(float)
    idx = {tr.name: bootstrap_indices(tr, n_boot, seed + i) for i, tr in enumerate(trials)}
    TR = {tr.name: tr for tr in trials}
    R = {key: interaction(TR[key[1]], sc) for key, sc in S.items()}
    BT = {key: bootstrap_interactions(TR[key[1]], sc, idx[key[1]]) for key, sc in S.items()}
    per = pd.DataFrame([dict(score=k[0], trial=k[1], est=v[0], se=v[1]) for k, v in R.items()])
    rows = []
    for nm in list(scores) + list(genes):
        ks = [tr.name for tr in trials if (nm, tr.name) in R]
        q = pool([R[(nm, t)][0] for t in ks], [R[(nm, t)][1] for t in ks]); row = dict(score=nm, **q)
        if nm in scores and len(ks) >= 2:
            sub = [TR[t] for t in ks]
            ref = random_reference(sub, {t: used[(nm, t)] for t in ks}, n_random, seed)
            row.update(ref_median=np.nanmedian(ref), ref_p95=np.nanpercentile(ref, 95),
                       emp_p=(1 + np.sum(ref >= q["est"])) / (n_random + 1), exceeds=bool(q["est"] > np.nanpercentile(ref, 95)))
        rows.append(row)
    prs = []
    for nm in scores:
        for g in genes:
            ks = [tr.name for tr in trials if (nm, tr.name) in R and (g, tr.name) in R]
            if len(ks) < 2:
                continue
            pdf = paired_difference([R[(nm, t)][0] - R[(g, t)][0] for t in ks], [BT[(nm, t)] for t in ks], [BT[(g, t)] for t in ks])
            prs.append(dict(score=nm, gene=g, **pdf))
    return dict(per_trial=per, pooled=pd.DataFrame(rows), paired=pd.DataFrame(prs))
