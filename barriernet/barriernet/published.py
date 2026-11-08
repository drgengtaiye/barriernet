"""Pooling published subgroup results of randomized trials (the manuscript's Fig. 5f).

Each trial reports the hazard ratio, checkpoint arm against control, with its 95% CI separately for patients with a
high and a low signature. The benefit-oriented log ratio of the two hazard ratios (positive: more benefit when the
signature is high) is pooled by random effects (REML with the Knapp-Hartung adjustment, as R metafor's
rma(method = "REML", test = "knha")) and converted approximately to a per-SD scale from the split reported."""
from __future__ import annotations
from importlib import resources
from typing import Dict, Optional, Sequence
import numpy as np
import pandas as pd
from scipy import optimize, stats

Z = stats.norm.ppf(0.975)
REQUIRED = ("trial", "endpoint", "hr_high", "lo_high", "hi_high", "hr_low", "lo_low", "hi_low")


def logse(lo, hi):
    return (np.log(hi) - np.log(lo)) / (2 * Z)


def pair_y(h, hl, hh, l, ll, lh):
    return -(np.log(h) - np.log(l)), np.sqrt(logse(hl, hh) ** 2 + logse(ll, lh) ** 2)


def per_sd_divisor(frac_high, frac_low):
    """Expected standard-normal score difference between the high and low groups of the split reported: top
    fraction p_hi against bottom fraction p_lo (p_lo = 1 - p_hi unless a middle group was set aside)."""
    frac_high = np.asarray(frac_high, float); frac_low = np.asarray(frac_low, float)
    p_hi = np.where(np.isnan(frac_high), 0.5, frac_high); p_lo = np.where(np.isnan(frac_low), 1 - p_hi, frac_low)
    return stats.norm.pdf(stats.norm.ppf(1 - p_hi)) / p_hi + stats.norm.pdf(stats.norm.ppf(p_lo)) / p_lo


def _reml_tau2(y, v, X):
    def nll(t2):
        w = 1.0 / (v + t2); XtWX = X.T @ (w[:, None] * X); b = np.linalg.solve(XtWX, X.T @ (w * y)); r = y - X @ b
        return 0.5 * (np.sum(np.log(v + t2)) + np.linalg.slogdet(XtWX)[1] + np.sum(w * r * r))
    hi = max(100.0, 10 * np.var(y) + 10 * np.max(v)); grid = np.concatenate([[0.0], np.geomspace(1e-8, hi, 400)])
    vals = [nll(g) for g in grid]; i = int(np.argmin(vals))
    if i == 0:
        r = optimize.minimize_scalar(nll, bounds=(0.0, grid[1]), method="bounded", options={"xatol": 1e-14})
        return float(r.x) if r.fun < vals[0] else 0.0
    r = optimize.minimize_scalar(nll, bounds=(grid[i - 1], grid[min(i + 1, len(grid) - 1)]), method="bounded", options={"xatol": 1e-14})
    return max(float(r.x), 0.0)


def rma(y, se, X=None) -> dict:
    """REML random-effects model with the Knapp-Hartung adjustment."""
    y = np.asarray(y, float); v = np.asarray(se, float) ** 2; k = len(y)
    X = np.ones((k, 1)) if X is None else np.column_stack([np.ones(k), np.asarray(X, float)]); p = X.shape[1]
    t2 = _reml_tau2(y, v, X); w = 1.0 / (v + t2); XtWX = X.T @ (w[:, None] * X); b = np.linalg.solve(XtWX, X.T @ (w * y))
    r = y - X @ b; s2 = float(np.sum(w * r * r) / (k - p)); vb = s2 * np.linalg.inv(XtWX); seb = np.sqrt(np.diag(vb))
    tc = stats.t.ppf(0.975, k - p); tv = b / seb; pv = 2 * stats.t.sf(np.abs(tv), k - p)
    wi = 1.0 / v; s2t = (k - 1) * wi.sum() / (wi.sum() ** 2 - np.sum(wi ** 2))
    out = dict(k=k, est=float(b[0]), se=float(seb[0]), lo=float(b[0] - tc * seb[0]), hi=float(b[0] + tc * seb[0]), p=float(pv[0]),
               tau2=t2, I2=float(100 * t2 / (t2 + s2t)), t=tv, pvals=pv)
    if p == 1:
        ps = np.sqrt(t2 + vb[0, 0]); out["pi_lo"], out["pi_hi"] = float(b[0] - tc * ps), float(b[0] + tc * ps)
    return out


def leave_one_out(y, se):
    y = np.asarray(y, float); se = np.asarray(se, float)
    return [rma(np.delete(y, i), np.delete(se, i)) for i in range(len(y))]


def egger(y, se):
    """Regression test of the estimates on their standard errors (metafor regtest for a knha fit)."""
    m = rma(y, se, X=np.asarray(se, float)); return float(m["t"][1]), float(m["pvals"][1])


def _col(df, c):
    return df[c].astype(float).values if c in df else np.full(len(df), np.nan)


def trial_inputs(df: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in REQUIRED if c not in df]
    if missing:
        raise ValueError("published-trials table lacks column(s): " + ", ".join(missing))
    E = df.copy()
    y, se = pair_y(*(_col(E, c) for c in ("hr_high", "lo_high", "hi_high", "hr_low", "lo_low", "hi_low")))
    if "int_est" in E:
        ui = np.isnan(y) & E["int_est"].notna().values & (E.get("int_scale", pd.Series("", index=E.index)).values == "ratio_high_over_low")
        y[ui] = -np.log(_col(E, "int_est")[ui]); se[ui] = logse(_col(E, "int_lo")[ui], _col(E, "int_hi")[ui])
    E["y"], E["se"] = y, se
    d = per_sd_divisor(_col(E, "frac_high"), _col(E, "frac_low")); E["y_sd"], E["se_sd"] = y / d, se / d
    E["pdl1_y"], E["pdl1_se"] = pair_y(*(_col(E, c) for c in ("pdl1_hr_pos", "pdl1_lo_pos", "pdl1_hi_pos", "pdl1_hr_neg", "pdl1_lo_neg", "pdl1_hi_neg")))
    E["tmb_y"], E["tmb_se"] = pair_y(*(_col(E, c) for c in ("tmb_hr_high", "tmb_lo_high", "tmb_hi_high", "tmb_hr_low", "tmb_lo_low", "tmb_hi_low")))
    E["survival"] = E["endpoint"].isin(["OS", "PFS"])
    E["exclude"] = E["exclude"].astype(str).str.lower().isin(("yes", "true", "1")) if "exclude" in E else False
    E["figure"] = E["from_figure"].astype(str).str.lower().str.startswith(("yes", "estimated")) if "from_figure" in E else False
    return E


def pool_published(df: pd.DataFrame) -> Dict[str, object]:
    """The analyses of the manuscript: primary (survival endpoints, rows not marked ``exclude``), per SD, values given
    in text or tables only, leave one out, Egger test, and PD-L1 and TMB over the same trials where reported."""
    I = trial_inputs(df); P = I[I["survival"] & ~I["exclude"] & np.isfinite(I["y"])]
    out: Dict[str, object] = {}

    def fit(name, d, y="y", se="se"):
        d = d[np.isfinite(d[y]) & np.isfinite(d[se])]
        if len(d) >= 2:
            m = rma(d[y].values, d[se].values); m["t"] = float(m["t"][0]); m.pop("pvals")
            out[name] = m | {"trials": "; ".join(d["trial"])}
    fit("primary", P); fit("per_sd", P, "y_sd", "se_sd"); fit("text_or_table_only", P[~P["figure"]]); fit("all_rows", I[I["survival"] & np.isfinite(I["y"])])
    fit("pdl1", P, "pdl1_y", "pdl1_se"); fit("tmb", P, "tmb_y", "tmb_se"); fit("signature_where_pdl1", P[np.isfinite(P["pdl1_y"])])
    if len(P) >= 3:
        l1 = leave_one_out(P["y"].values, P["se"].values)
        out["leave_one_out"] = dict(k=len(P) - 1, lo=min(m["est"] for m in l1), hi=max(m["est"] for m in l1), p=max(m["p"] for m in l1))
    if len(P) >= 10:
        t, p = egger(P["y"].values, P["se"].values); out["egger"] = dict(k=len(P), t=t, p=p)
    out["per_trial"] = I
    return out


def example_table() -> pd.DataFrame:
    """The 19 published randomized comparisons extracted for the manuscript (Supplementary Table S28); IMmotion150,
    which is also among the six comparisons with patient-level data, is marked ``exclude``."""
    with resources.files("barriernet").joinpath("data/published_trials.tsv").open("r", encoding="utf-8") as fh:
        return pd.read_csv(fh, sep="\t")
