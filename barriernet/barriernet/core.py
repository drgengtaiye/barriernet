from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple

def zscore_series(x: pd.Series) -> pd.Series:
    """Per-feature z-score across samples (population SD, ddof = 0)."""
    x = pd.to_numeric(x, errors="coerce")
    mu = np.nanmean(x.values)
    sd = np.nanstd(x.values, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x*0.0
    return (x - mu) / (sd + 1e-9)

def zscore_df_cols(df: pd.DataFrame) -> pd.DataFrame:
    return df.apply(zscore_series, axis=0)

def score_modules(expr_samples_by_genes: pd.DataFrame,
                  modules: Dict[str, List[str]],
                  min_genes: int = 3,
                  zscore: bool = True) -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    """Module score = mean of per-gene z-scores (z across samples) over the module genes present
    in `expr_samples_by_genes`; NaN when fewer than `min_genes` are present."""
    X = expr_samples_by_genes.copy()
    X.columns = [str(c).upper() for c in X.columns]
    if zscore:
        X = zscore_df_cols(X)

    scores = pd.DataFrame(index=X.index)
    used: Dict[str, List[str]] = {}
    for m, genes in modules.items():
        g = [str(x).upper() for x in genes]
        g2 = [x for x in g if x in X.columns]
        used[m] = g2
        scores[m] = np.nan if len(g2) < min_genes else X[g2].mean(axis=1, skipna=True)
    return scores, used

def score_axes(module_scores: pd.DataFrame, axes: Dict[str, dict]) -> pd.DataFrame:
    """Axis = mean(positive modules) - mean(negative modules). Modules that could not be scored
    (too few genes) are skipped with a warning, so the axis is built from the rest."""
    out = pd.DataFrame(index=module_scores.index)
    for ax, spec in axes.items():
        pos = [x for x in spec.get("positive_modules", []) if x in module_scores.columns]
        neg = [x for x in spec.get("negative_modules", []) if x in module_scores.columns]
        missing = [m for m in spec.get("positive_modules", []) + spec.get("negative_modules", [])
                   if m not in module_scores.columns or module_scores[m].isna().all()]
        if missing:
            warnings.warn(f"Axis {ax}: module(s) {missing} unavailable; axis built from the remaining modules.")
        pos_v = module_scores[pos].mean(axis=1, skipna=True) if len(pos) else 0.0
        neg_v = module_scores[neg].mean(axis=1, skipna=True) if len(neg) else 0.0
        out[ax] = pos_v - neg_v
    return out

def summarise_used_genes(used: Dict[str, List[str]]) -> pd.DataFrame:
    rows = [{"module": m, "n_used": len(g), "genes_used": ",".join(g)} for m, g in used.items()]
    return pd.DataFrame(rows)
