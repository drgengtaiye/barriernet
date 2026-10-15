from __future__ import annotations
import numpy as np
import pandas as pd
from typing import Tuple

def _marker_genes():
    from .genesets import load_yaml, default_yaml
    gs = load_yaml(default_yaml("rna"))
    return [g for v in gs.modules.values() for g in v]


def looks_linear(expr: pd.DataFrame) -> bool:
    """True when values look like un-logged TPM/counts (large maximum, strong right skew)."""
    v = expr.to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    return bool(v.size and np.nanmax(v) > 100 and np.nanmedian(v) >= 0)


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")

def load_expr_matrix(path: str, gene_col: str | None=None, orientation: str="auto") -> Tuple[pd.DataFrame, str]:
    """Return (expr, orientation) with expr as samples×genes."""
    df = read_table(path)

    # Set gene index if provided
    if gene_col is not None and gene_col in df.columns:
        df = df.set_index(gene_col)

    # Try autodetect gene column in first column
    if gene_col is None and df.index.dtype == int and df.shape[1] > 1:
        c0 = df.columns[0]
        if df[c0].dtype == object:
            tmp = df.set_index(c0)
            numfrac = tmp.apply(lambda s: pd.to_numeric(s, errors='coerce')).notna().mean().mean()
            if numfrac > 0.6:
                df = tmp

    if orientation not in ("auto","samples_by_genes","genes_by_samples"):
        raise ValueError("orientation must be auto|samples_by_genes|genes_by_samples")

    if orientation == "auto":
        # Decide which axis holds gene symbols by overlap with the bundled marker genes; fall
        # back to shape (genes usually outnumber samples) only when neither axis matches.
        # (<= 0.2.0 used shape alone, which misreads e.g. 500 samples x 200 genes.)
        markers = set(_marker_genes())
        in_rows = len(markers & {str(x).upper() for x in df.index})
        in_cols = len(markers & {str(x).upper() for x in df.columns})
        if in_rows != in_cols:
            orientation = "genes_by_samples" if in_rows > in_cols else "samples_by_genes"
        else:
            orientation = "genes_by_samples" if df.shape[0] > df.shape[1] * 2 else "samples_by_genes"

    # Coerce numeric
    df = df.apply(lambda s: pd.to_numeric(s, errors="coerce"))

    expr = df.T if orientation == "genes_by_samples" else df
    expr.index = expr.index.astype(str)
    expr.columns = expr.columns.astype(str)
    return expr, orientation
