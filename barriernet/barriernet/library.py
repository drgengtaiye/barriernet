"""The published immune-signature library and the T-cell references benchmarked in the manuscript.

The 52 immuno-oncology entries come from the SignatureSets collection (BHKLAB, https://github.com/bhklab/SignatureSets,
MIT License; see data/SignatureSets_LICENSE.txt), with their published gene weights (weight 1 where none is given),
declared direction (sensitivity +1, resistance -1) and the family used in the manuscript (Inflamed/IFN, Other,
Confounder). TIS_Damotte and Tcell-inflamed_Ayers hold the same 18 genes, so the distinct library has 51 signatures.
"""
from __future__ import annotations
from importlib import resources
from typing import Dict, Optional, Tuple
import numpy as np
import pandas as pd

T_CELL_PANEL = ["CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "CD2", "CD247", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7"]
IFN_CHECKPOINT_PANEL = ["IFNG", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCR5", "STAT1", "IDO1", "CD274", "PDCD1", "LAG3",
                        "TIGIT", "CTLA4", "HAVCR2"]
PANEL_GENES = T_CELL_PANEL + IFN_CHECKPOINT_PANEL
REFERENCES = {"T-cell mean (13 genes)": T_CELL_PANEL, "T-cell mean (CD8A, PRF1, GZMA)": ["CD8A", "PRF1", "GZMA"]}
MAIN_REFERENCE = "T-cell mean (13 genes)"


def load_library(distinct: bool = True) -> pd.DataFrame:
    """Long table: signature, family, direction, published_engine, pmid, gene, weight, duplicate_of."""
    with resources.files("barriernet").joinpath("data/signatures_io.tsv").open("r", encoding="utf-8") as fh:
        L = pd.read_csv(fh, sep="\t", dtype={"pmid": str}, keep_default_na=False, na_values=[""])
    L["duplicate_of"] = L["duplicate_of"].fillna("")
    return L[L["duplicate_of"] == ""].reset_index(drop=True) if distinct else L


def read_signature_table(path) -> pd.DataFrame:
    """A user's signatures in the same long format; needs signature and gene, optional weight (default 1), family
    (default 'User') and direction (default +1)."""
    L = pd.read_csv(path, sep=None, engine="python")
    missing = {"signature", "gene"} - set(L.columns)
    if missing:
        raise ValueError(f"signature table lacks column(s): {', '.join(sorted(missing))}")
    L = L.copy(); L["gene"] = L["gene"].astype(str).str.upper()
    for c, v in (("weight", 1.0), ("family", "User"), ("direction", 1), ("published_engine", ""), ("pmid", ""), ("duplicate_of", "")):
        if c not in L:
            L[c] = v
    L["weight"] = pd.to_numeric(L["weight"], errors="coerce").fillna(1.0)
    return L.drop_duplicates(["signature", "gene"]).reset_index(drop=True)


def combine(library: pd.DataFrame, user: pd.DataFrame) -> pd.DataFrame:
    """The library plus a user's signatures; a user signature may not reuse a library name."""
    clash = sorted(set(user.signature) & set(library.signature))
    if clash:
        raise ValueError("signature name(s) already in the library, rename them: " + ", ".join(clash))
    return pd.concat([library, user], ignore_index=True)


def signature_weights(L: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    return {s: dict(zip(g["gene"], g["weight"].astype(float))) for s, g in L.groupby("signature", sort=False)}


def signature_meta(L: pd.DataFrame) -> pd.DataFrame:
    M = L.groupby("signature", sort=False).agg(family=("family", "first"), direction=("direction", "first"),
                                               published_engine=("published_engine", "first"), pmid=("pmid", "first"),
                                               n_genes=("gene", "size"))
    return M.reset_index()


def reference_weights(references: Dict[str, list] | None = None) -> Dict[str, Dict[str, float]]:
    return {k: {g: 1.0 for g in v} for k, v in (references or REFERENCES).items()}


def score_signature(Z: pd.DataFrame, weights: Dict[str, float], min_genes: int = 3, min_frac: float = 0.5
                    ) -> Tuple[Optional[np.ndarray], list, Optional[np.ndarray]]:
    """The manuscript's common engine: sum(w z) / sum(|w|) over the listed genes that are present and non-constant,
    missing values skipped sample by sample; scored when at least ``min_genes`` genes and ``min_frac`` of the list
    are present. Returns (score or None, genes used, weights used)."""
    sd = Z.std(axis=1, skipna=True)
    genes = [g for g in weights if g in Z.index and np.isfinite(sd[g]) and sd[g] > 0]
    if len(genes) < min_genes or len(genes) < min_frac * len(weights):
        return None, genes, None
    w = np.array([weights[g] for g in genes], float); M = Z.loc[genes].to_numpy(float)
    W = np.where(np.isfinite(M), 1.0, 0.0) * np.abs(w)[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        sc = np.nansum(M * w[:, None], 0) / W.sum(0)
    return (sc if np.isfinite(sc).all() else None), genes, w
