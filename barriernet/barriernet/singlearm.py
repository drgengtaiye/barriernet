"""Single-arm checkpoint cohorts: log-odds of response per SD of a score, pooled over cohorts, and the agreement of
that ranking with the randomized interaction (the manuscript's Fig. 2h). Single-arm associations cannot separate a
predictive from a prognostic score; they are reported for comparison with the randomized evidence."""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Optional, Sequence
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr, t as tdist
from .calibrate import logistic, standardize
from .library import score_signature


@dataclass
class Cohort:
    """One single-arm cohort. ``Z``: genes x samples, z-scored within the cohort; ``response``: 1 = responder."""
    name: str
    Z: pd.DataFrame
    response: np.ndarray
    covariates: Optional[np.ndarray] = None


def pool_dl_hk(y, se) -> dict:
    """DerSimonian-Laird estimate; CI and P from the wider of the Hartung-Knapp (DL tau2) and DL intervals."""
    y, se = np.asarray(y, float), np.asarray(se, float); ok = np.isfinite(y) & np.isfinite(se) & (se > 0); y, v = y[ok], se[ok] ** 2
    k = len(y)
    if k < 2:
        return dict(est=np.nan, lo=np.nan, hi=np.nan, p=np.nan, k=k, tau2=np.nan)
    w0 = 1 / v; m0 = np.sum(w0 * y) / w0.sum(); Q = np.sum(w0 * (y - m0) ** 2)
    t2 = max(0.0, (Q - (k - 1)) / (w0.sum() - np.sum(w0 ** 2) / w0.sum())); w = 1 / (v + t2); mu = np.sum(w * y) / w.sum()
    sd = np.sqrt(1 / w.sum()); z = norm.ppf(0.975); dl = (mu - z * sd, mu + z * sd, 2 * norm.sf(abs(mu / sd)))
    sh = np.sqrt(np.sum(w * (y - mu) ** 2) / ((k - 1) * w.sum())); q = tdist.ppf(0.975, k - 1)
    hk = (mu - q * sh, mu + q * sh, 2 * tdist.sf(abs(mu / sh), k - 1)) if sh > 0 else dl
    lo, hi = (hk[0], hk[1]) if (hk[1] - hk[0]) >= (dl[1] - dl[0]) else (dl[0], dl[1])
    return dict(est=mu, lo=lo, hi=hi, p=max(dl[2], hk[2]), k=k, tau2=t2)


def pool_fixed(y, se) -> dict:
    y, se = np.asarray(y, float), np.asarray(se, float); ok = np.isfinite(y) & np.isfinite(se) & (se > 0); y, se = y[ok], se[ok]
    w = 1 / se ** 2; mu = np.sum(w * y) / w.sum(); s = np.sqrt(1 / w.sum()); z = norm.ppf(0.975)
    return dict(est=mu, lo=mu - z * s, hi=mu + z * s, p=2 * norm.sf(abs(mu / s)), k=len(y))


def single_arm(cohorts: Sequence[Cohort], scores: Dict[str, Dict[str, float]], min_cohorts: Optional[int] = None):
    """Per cohort, Firth-fallback logistic regression of response on the standardized score; pooled by the
    DL/Hartung-Knapp rule (and by fixed effect) when scored in at least ``min_cohorts`` (default half the cohorts)."""
    rows = []
    for c in cohorts:
        y = np.asarray(c.response, float)
        for nm, w in scores.items():
            sc, used, _ = score_signature(c.Z, w, min_genes=1 if len(w) == 1 else 3)
            if sc is None or np.std(sc) == 0:
                continue
            cov = [] if c.covariates is None else [np.asarray(c.covariates, float).reshape(len(y), -1)]
            b, se = logistic(np.column_stack([np.ones(len(y)), standardize(sc)] + cov), y, True)
            rows.append(dict(score=nm, cohort=c.name, n=len(y), responders=int(y.sum()), genes=len(used), logOR=b[1], se=se[1]))
    PC = pd.DataFrame(rows)
    need = min_cohorts or max(2, int(np.ceil(len(cohorts) / 2)))
    out = []
    for nm, g in PC.groupby("score", sort=False) if len(PC) else []:
        if len(g) >= need:
            r = pool_dl_hk(g.logOR, g.se); f = pool_fixed(g.logOR, g.se)
            out.append(dict(score=nm, **r, fixed_est=f["est"], fixed_lo=f["lo"], fixed_hi=f["hi"], fixed_p=f["p"]))
    return dict(per_cohort=PC, pooled=pd.DataFrame(out))


def concordance(single_pooled: pd.DataFrame, randomized_pooled: pd.DataFrame) -> dict:
    """Spearman correlation between pooled single-arm log-odds and pooled randomized interactions, over the scores
    pooled in both."""
    T = randomized_pooled[["score", "est", "p"]].rename(columns={"est": "interaction", "p": "interaction_p"}).merge(
        single_pooled[["score", "est", "p"]].rename(columns={"est": "single_arm", "p": "single_arm_p"}), on="score")
    T = T.dropna(subset=["interaction", "single_arm"])
    r = spearmanr(T["interaction"], T["single_arm"]) if len(T) > 2 else (np.nan, np.nan)
    return dict(rho=float(r[0]), p=float(r[1]), n=len(T), table=T)
