#!/usr/bin/env python3
"""E17 - independent single-cell attribution of the barriernet marker modules.

Pre-registered in analysis/E17_E22_prereg.md (section E17 and the 2026-10-04 "E17 data" deviation).

    python E17_singlecell_attribution.py --stage extract   # cell tables, 48-gene matrices, class mapping (no scores)
    python E17_singlecell_attribution.py --stage score     # attribution, sensitivity, expectations
    python E17_singlecell_attribution.py --stage figure    # figures/Fig_E17_attribution.{png,pdf}
    python E17_singlecell_attribution.py                   # all three, in that order

Method = the discovery-atlas code, GC_TME_Reprogrammer/validation/attribution/cell_attribution.py
(snapshot: code/analysis/00_upstream/cell_attribution.py), which produced the published
z = 2.6 / 2.8 and 73% / 91% / 28%:
  * per-cell module score = mean log-normalised expression over the module genes present (>= 3);
  * mean score per cell class;
  * z = (class mean - mean over classes) / (SD over classes, pandas ddof=1, + 1e-9);
  * share of signal = (cell fraction x class mean) / sum over classes of |cell fraction x class mean| x 100.
The discovery atlas (GSE183904, processed h5ad) is re-run through this file to show the definitions
reproduce the published numbers.

Data (tumour samples only): OMIX001073 (Sun 2022; primary), GSE206785 (Kang 2022; secondary).
Memory: only the 48 module genes are kept after one pass over each matrix (no cell subsampling).
Environment: conda env r45 (Python 3.14, numpy 2.3, pandas 2.3, scipy 1.17, pyarrow 25, h5py 3.16).
"""
import argparse, functools, gzip, json, os, subprocess, sys, time, warnings
from pathlib import Path

warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np
import pandas as pd

ROOT = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2")
RES = ROOT / "analysis/results"
CACHE = ROOT / "analysis/cache_E17"
FIGDIR = ROOT / "figures"
VAL = ROOT / "data_download/validation"
OMIX = VAL / "scrna_OMIX001073"
KANG = VAL / "scrna_GSE206785"
GC = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer")
DISC_H5 = GC / "out_scanpy_strict_GSE183904/GSE183904_TumorOnly_scanpy_strict_final_scored.h5ad"
DISC_PUB = GC / "validation/attribution"          # outputs of cell_attribution.py (published numbers)
SCRATCH = Path(os.environ.get("E17_SCRATCH", "/groups/adv2105_gp/yichen/tmp/E17_scratch"))

# ---------------------------------------------------------------- fixed definitions
MODULES = {   # barriernet 0.2.1, identical to cell_attribution.py
    "CAF_ECM": ["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
    "TAM_SUP": ["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
    "APM": ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    "IFN": ["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
    "TLS_B": ["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
    "T_CELL": ["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"],
}
MODS = list(MODULES)
GENES48 = [g for m in MODS for g in MODULES[m]]
GENE_MODULE = {g: m for m in MODS for g in MODULES[m]}
MIN_GENES = 3

CLASSES = ["T cell", "NK cell", "B cell", "plasma cell", "macrophage/monocyte", "dendritic cell",
           "mast cell", "fibroblast", "endothelial", "pericyte/smooth muscle", "epithelial/tumour"]
SHORT = {"T cell": "T", "NK cell": "NK", "B cell": "B", "plasma cell": "Plasma",
         "macrophage/monocyte": "Mac/Mono", "dendritic cell": "DC", "mast cell": "Mast",
         "fibroblast": "Fibro", "endothelial": "Endo", "pericyte/smooth muscle": "Peri/SMC",
         "epithelial/tumour": "Epi"}
LYMPHOID = ["T cell", "NK cell", "B cell", "plasma cell"]     # "lymphoid cells" of the 91% (Tcell+NK+Bcell+Plasma)
TNK = ["T cell", "NK cell"]
INTENDED = {"CAF_ECM": ["fibroblast"], "TAM_SUP": ["macrophage/monocyte"], "TLS_B": ["B cell", "plasma cell"],
            "T_CELL": ["T cell", "NK cell"], "APM": None, "IFN": None}      # prereg E17 expectations
LINEAGE = {"CAF_ECM": ["fibroblast", "pericyte/smooth muscle"],          # looser gene-level reading
           "TAM_SUP": ["macrophage/monocyte", "dendritic cell"],
           "TLS_B": ["B cell", "plasma cell"], "T_CELL": ["T cell", "NK cell"], "APM": None, "IFN": None}
DATASETS = ["GSE183904", "OMIX001073", "GSE206785"]
DSLAB = {"GSE183904": "GSE183904 (discovery)", "OMIX001073": "OMIX001073 (primary)",
         "GSE206785": "GSE206785 (secondary)"}
MIN_CELLS_PT, MIN_TYPES_PT = 20, 6      # per-patient top-class rule, from icb/attrib_per_patient.py
TENRICHED = ["GC02T-Tcell", "GC10T-Tcell"]
OMIX_DC = {"cDC1_XCR1", "cDC2_CD1C", "DC_LAMP3", "pDC_LILRA4"}
OMIX_MAC = {"M\u03c6_APOE", "M\u03c6_THBS1", "Mono_CD14", "Mono_FCGR3A"}

LOG = None


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    if LOG is not None:
        with open(LOG, "a") as fh:
            fh.write(s + "\n")


# ===================================================================== stage 1: extract
def _harmonise(lin_mod, total):
    """ln(1 + CP10k) from linear values: the discovery atlas's normalize_total(1e4) + log1p."""
    return np.log1p(lin_mod / total[:, None] * 1e4).astype(np.float32)


def extract_discovery():
    import h5py
    from anndata.io import read_elem
    t0 = time.time()
    with h5py.File(DISC_H5, "r") as f:
        obs = read_elem(f["obs"])
        rvar = read_elem(f["raw/var"]).index.astype(str)
        xvar = read_elem(f["var"]).index.astype(str)
        assert (rvar == xvar).all()
        genes_up = pd.Index([g.upper() for g in rvar])          # cell_attribution.py upper-cases names
        pres = [g for g in GENES48 if g in genes_up]
        gi = np.array([genes_up.get_loc(g) for g in pres])
        out = {}
        for key in ("raw/X", "X"):
            g = f[key]
            data, ind, ptr = g["data"][:], g["indices"][:], g["indptr"][:]
            look = np.full(len(rvar), -1, np.int64); look[gi] = np.arange(len(gi))
            pos = look[ind]; keep = pos >= 0
            rows = np.repeat(np.arange(len(ptr) - 1), np.diff(ptr))
            M = np.zeros((len(ptr) - 1, len(gi)), np.float32)
            M[rows[keep], pos[keep]] = data[keep]
            out[key] = M
            if key == "raw/X":
                tot = np.bincount(rows, weights=np.expm1(data.astype(np.float64)), minlength=len(ptr) - 1)
            del data, ind, ptr, rows, pos, keep
    same = np.array_equal(out["raw/X"], out["X"])
    log(f"[discovery] {len(obs):,} cells, {len(rvar):,} genes; module genes present {len(pres)}/48; "
        f"missing {[g for g in GENES48 if g not in pres]}; raw.X == X on module genes: {same}; "
        f"per-cell sum of expm1 (all genes): median {np.median(tot):.1f}, range {tot.min():.1f}-{tot.max():.1f} "
        f"-> ln(1+CP10k), used as is; {time.time()-t0:.0f}s")
    M = out["raw/X"]                       # cell_attribution.py uses adata.raw when present
    cmap = {"Epi": "epithelial/tumour", "Fibro": "fibroblast", "Endo": "endothelial", "Peri": "pericyte/smooth muscle",
            "Macro": "macrophage/monocyte", "DC": "dendritic cell", "Tcell": "T cell", "NK": "NK cell",
            "Bcell": "B cell", "Plasma": "plasma cell"}
    o = pd.DataFrame({"cell": obs.index.astype(str), "patient": obs["sample_id"].astype(str).values,
                      "library": obs["sample_id"].astype(str).values, "platform": "10x",
                      "author_major": obs["anno_level1"].astype(str).values, "author_sub": ""})
    o["class"] = o.author_major.map(cmap)
    o["fallback"] = False; o["tenriched"] = False
    return o, pres, M, M.copy(), tot


def extract_omix():
    import scipy.io
    t0 = time.time()
    ann = pd.read_csv(OMIX / "OMIX001073-20-04.tsv", sep="\t", index_col=0)
    bc = pd.read_csv(OMIX / "all_cell_matrix/barcodes.tsv.gz", header=None)[0].values
    assert (bc == ann.index.values).all(), "barcode order differs from annotation order"
    feat = pd.read_csv(OMIX / "all_cell_matrix/features.tsv.gz", sep="\t", header=None)
    sym = feat[1].astype(str).values
    gi = np.array([np.where(sym == g)[0][0] for g in GENES48])
    assert all((sym == g).sum() == 1 for g in GENES48)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    plain = SCRATCH / "omix_matrix.mtx"
    if not plain.exists():
        subprocess.run(f"pigz -dc -p 16 {OMIX/'all_cell_matrix/matrix.mtx.gz'} > {plain}.tmp && mv {plain}.tmp {plain}",
                       shell=True, check=True)
    A = scipy.io.mmread(str(plain))           # genes x cells, COO
    r, c, v = A.row, A.col, A.data
    del A
    ncell = len(bc)
    tot = np.bincount(c, weights=np.exp2(v) - 1.0, minlength=ncell)     # deposited values are log2(1 + x)
    look = np.full(len(sym), -1, np.int64); look[gi] = np.arange(len(gi))
    pos = look[r]; keep = pos >= 0
    dep = np.zeros((ncell, len(gi)), np.float32)
    dep[c[keep], pos[keep]] = v[keep]
    nnz_cell = np.bincount(c, minlength=ncell)
    del r, c, v, pos, keep
    log(f"[OMIX001073] matrix {len(sym):,} genes x {ncell:,} cells; nnz/cell == n_genes: "
        f"{(nnz_cell == ann.n_genes.values).all()}; sum(2^v-1) median {np.median(tot):.0f}; {time.time()-t0:.0f}s")
    tum = (ann.tissue == "Tumor").values
    ann = ann[tum]; dep = dep[tum]; tot = tot[tum]
    lin = np.exp2(dep.astype(np.float64)) - 1.0
    harm = _harmonise(lin, tot)

    # author sub-annotation for the two major types that span two of our classes
    T = pd.read_csv(OMIX / "OMIX001073-17.tsv", sep="\t", index_col=0, usecols=[0, 5])
    NK = pd.read_csv(OMIX / "OMIX001073-20-09.tsv", sep="\t", index_col=0, usecols=[0, 7])
    MY = pd.read_csv(OMIX / "OMIX001073-20-15.tsv", sep="\t", index_col=0, usecols=[0, 7])
    assert list(T.columns) == ["T_cluster"] and list(NK.columns) == ["NK_cluster"] and list(MY.columns) == ["myeloid_cluster"]
    sub = pd.Series("", index=ann.index, dtype=object)
    sub.loc[sub.index.intersection(T.index)] = T.loc[sub.index.intersection(T.index), "T_cluster"].astype(str)
    sub.loc[sub.index.intersection(NK.index)] = NK.loc[sub.index.intersection(NK.index), "NK_cluster"].astype(str)
    sub.loc[sub.index.intersection(MY.index)] = MY.loc[sub.index.intersection(MY.index), "myeloid_cluster"].astype(str)
    major_map = {"Epithelial cells": "epithelial/tumour", "Endocrine cells": "epithelial/tumour",
                 "Fibroblasts": "fibroblast", "Endothelial cells": "endothelial",
                 "Smooth muscle cells": "pericyte/smooth muscle", "Mast cells": "mast cell", "B cells": "B cell",
                 "B cells(Plasma cells)": "plasma cell", "Erythrocytes": None}
    cls, fb = [], []
    for mj, s in zip(ann.cluster.values, sub.values):
        if mj == "T cells & NK cells":
            if s == "":
                cls.append("T cell"); fb.append(True)                 # no author subcluster -> majority class
            elif s.startswith("NK_"):
                cls.append("NK cell"); fb.append(False)               # NK_FCGR3A/GNLY/IL7R/XCL1
            else:
                cls.append("T cell"); fb.append(False)                # CD4/CD8/Treg/cycling, NKT_*, ab-T, gd-T
        elif mj == "Myeloid cells":
            if s == "":
                cls.append("macrophage/monocyte"); fb.append(True)    # no author subcluster -> majority class
            elif s in OMIX_DC:
                cls.append("dendritic cell"); fb.append(False)
            else:
                assert s in OMIX_MAC, s
                cls.append("macrophage/monocyte"); fb.append(False)
        else:
            cls.append(major_map[mj]); fb.append(False)
    o = pd.DataFrame({"cell": ann.index.values, "patient": ann.patient.astype(str).values,
                      "library": ann.batch.astype(str).values, "platform": "10x 5'",
                      "author_major": ann.cluster.astype(str).values,
                      "author_sub": np.where(np.isin(ann.cluster.values, ["T cells & NK cells", "Myeloid cells"]),
                                             np.where(sub.values == "", "(no subcluster)", sub.values), ""),
                      "class": cls, "fallback": fb})
    o["tenriched"] = o.library.isin(TENRICHED)
    log(f"[OMIX001073] tumour cells {len(o):,}; T-cell-enriched libraries {o.tenriched.sum():,}")
    return o, list(GENES48), harm, dep, tot


def extract_kang():
    import pyarrow as pa, pyarrow.csv as pcsv
    t0 = time.time()
    meta = pd.read_csv(KANG / "GSE206785_metadata.txt.gz", sep="\t")
    path = KANG / "GSE206785_scgex.txt.gz"
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split(",")
    genes = hdr[1:]
    gpos = {g: i for i, g in enumerate(genes)}
    gi = [gpos[g] for g in GENES48]
    types = {h: pa.float32() for h in genes}; types[hdr[0]] = pa.string()
    rd = pcsv.open_csv(pa.input_stream(str(path), compression="gzip"),
                       read_options=pcsv.ReadOptions(block_size=1 << 26, use_threads=True),
                       convert_options=pcsv.ConvertOptions(column_types=types))
    cells, deps, tots = [], [], []
    for b in rd:
        cells.append(np.asarray(b.column(0).to_pylist(), dtype=object))
        cols = [b.column(j + 1).to_numpy(zero_copy_only=False) for j in range(len(genes))]
        tot = np.zeros(b.num_rows)
        for x in cols:
            tot += np.expm1(x.astype(np.float64))
        tots.append(tot)
        deps.append(np.column_stack([cols[j] for j in gi]).astype(np.float32))
    cell = np.concatenate(cells); dep = np.vstack(deps); tot = np.concatenate(tots)
    assert len(cell) == len(meta)
    suffix = np.array([c.rsplit("-", 1)[1] for c in cell])
    assert (suffix == meta.Sample.astype(str).values).all(), "matrix rows and metadata rows are not aligned"
    lin = np.expm1(dep.astype(np.float64))
    log(f"[GSE206785] {len(cell):,} cells x {len(genes):,} genes; row order matches metadata (sample suffix of "
        f"every cell ID equals metadata Sample); expm1 on module genes integer to within "
        f"{np.abs(lin - np.round(lin)).max():.4f}; sum(expm1) median {np.median(tot):.0f}; {time.time()-t0:.0f}s")
    tum = (meta.Tissue == "Tumor").values
    meta = meta[tum].reset_index(drop=True); dep = dep[tum]; tot = tot[tum]; cell = cell[tum]
    harm = _harmonise(np.expm1(dep.astype(np.float64)), tot)
    tmap = {"Epithelial": "epithelial/tumour", "Fibroblast": "fibroblast", "Endothelial": "endothelial",
            "Mural": "pericyte/smooth muscle", "Mast": "mast cell", "CD4+ T": "T cell", "CD8+ T": "T cell",
            "Innate lymphoid": "NK cell", "B": "B cell", "Plasma": "plasma cell", "Glial": None}
    cls, fb = [], []
    for ty, an in zip(meta.Type.values, meta.Annotation.astype(str).values):
        if ty == "Myeloid":
            if "_Mono_" in an or "_Macro_" in an:
                cls.append("macrophage/monocyte"); fb.append(False)
            elif "_cDC_" in an or "_pDC_" in an:
                cls.append("dendritic cell"); fb.append(False)
            elif "_Neutro_" in an:
                cls.append(None); fb.append(False)                     # neutrophils: none of the 11 classes
            else:
                cls.append("macrophage/monocyte"); fb.append(True)    # lowQual / UnIdent -> majority class
        else:
            cls.append(tmap[ty]); fb.append(False)
    o = pd.DataFrame({"cell": cell, "patient": meta.Patient.astype(str).values, "library": meta.Sample.astype(str).values,
                      "platform": meta.Platform.astype(str).values, "author_major": meta.Type.astype(str).values,
                      "author_sub": np.where(meta.Type.values == "Myeloid", meta.Annotation.astype(str).values, ""),
                      "class": cls, "fallback": fb})
    o["tenriched"] = False
    log(f"[GSE206785] tumour cells {len(o):,}")
    return o, list(GENES48), harm, dep, tot


def stage_extract():
    CACHE.mkdir(parents=True, exist_ok=True)
    maps = []
    for ds, fn in [("GSE183904", extract_discovery), ("OMIX001073", extract_omix), ("GSE206785", extract_kang)]:
        o, genes, harm, dep, tot = fn()
        o["class"] = o["class"].astype(object).where(o["class"].notna(), "(excluded)")
        o.to_csv(CACHE / f"{ds}_obs.csv.gz", index=False)
        np.savez_compressed(CACHE / f"{ds}_expr.npz", genes=np.array(genes), harm=harm, dep=dep, total_lin=tot)
        g = o.groupby(["author_major", "author_sub", "class", "fallback"], dropna=False)
        m = g.size().rename("n_tumour_cells").reset_index()
        if ds == "OMIX001073":
            m = m.merge(o[~o.tenriched].groupby(["author_major", "author_sub", "class", "fallback"]).size()
                        .rename("n_excl_Tenriched").reset_index(), how="left").fillna({"n_excl_Tenriched": 0})
        m.insert(0, "dataset", ds)
        maps.append(m)
    M = pd.concat(maps, ignore_index=True)
    M["rule"] = np.select(
        [M["class"] == "(excluded)", M.fallback, M.author_sub.astype(str) != ""],
        ["not one of the 11 classes: excluded before scoring",
         "author major type spans two classes and the cell has no informative author subtype: majority class "
         "(excluded in the 'no_fallback' sensitivity)",
         "author major type spans two classes: author subtype decides"],
        "author major type")
    M.to_csv(RES / "E17_celltype_mapping.csv", index=False)
    log(f"\nwrote {RES/'E17_celltype_mapping.csv'} ({len(M)} rows) before any module score was computed")
    for ds in DATASETS:
        s = M[M.dataset == ds].groupby("class").n_tumour_cells.sum()
        log(f"  {ds}: " + ", ".join(f"{SHORT.get(k, k)} {v:,}" for k, v in s.reindex(CLASSES + ['(excluded)']).dropna()
                                    .astype(int).items()))


# ===================================================================== stage 2: score
def load(ds):
    o = pd.read_csv(CACHE / f"{ds}_obs.csv.gz", keep_default_na=False, low_memory=False)
    z = np.load(CACHE / f"{ds}_expr.npz", allow_pickle=False)
    o["fallback"] = o.fallback.astype(str) == "True"
    o["tenriched"] = o.tenriched.astype(str) == "True"
    return o, list(z["genes"]), z["harm"], z["dep"]


def cell_scores(X, genes, keep_genes=None, gene_z=False):
    """Per-cell module score as in cell_attribution.py: mean expression of the module genes present (>= 3).
    gene_z=True: the Supplementary-Methods wording (mean of per-gene z-scores over the cells analysed)."""
    gl = list(genes)
    if gene_z:
        mu = X.mean(0, dtype=np.float64); sd = X.std(0, dtype=np.float64); sd[sd == 0] = 1
        X = ((X - mu) / sd).astype(np.float32)
    S, used = {}, {}
    for m in MODS:
        av = [g for g in MODULES[m] if g in gl and (keep_genes is None or g in keep_genes)]
        used[m] = av
        S[m] = X[:, [gl.index(g) for g in av]].mean(1, dtype=np.float64) if len(av) >= MIN_GENES else np.full(len(X), np.nan)
    return pd.DataFrame(S), used


def attribution(S, cls):
    """cell_attribution.py, sections 4-5: class means, z across classes, fraction-weighted share."""
    present = [c for c in CLASSES if (cls == c).any()]
    attr = S.groupby(cls.values).mean().reindex(present)
    z = attr.copy()
    for col in z.columns:
        mu, sd = z[col].mean(), z[col].std()
        z[col] = (z[col] - mu) / (sd + 1e-9) if sd > 0 else 0
    cnt = cls.value_counts().reindex(present)
    frac = cnt / cnt.sum()
    contrib = attr.mul(frac, axis=0)
    share = contrib / contrib.abs().sum() * 100
    return attr, z, share, cnt, frac


def top_rows(ds, an, z, share, attr):
    rows = []
    for m in MODS:
        if z[m].isna().all():
            continue
        zs = z[m].sort_values(ascending=False)
        it = INTENDED[m]
        row = dict(dataset=ds, analysis=an, module=m, top_class=zs.index[0], top_z=zs.iloc[0],
                   second_class=zs.index[1], second_z=zs.iloc[1],
                   intended=";".join(it) if it else "(not prespecified)",
                   intended_z=max(z.loc[c, m] for c in it if c in z.index) if it else np.nan,
                   top_is_intended=(zs.index[0] in it) if it else np.nan)
        if share is not None:
            row.update(share_fibroblast=share.reindex(["fibroblast"])[m].sum(),
                       share_macrophage=share.reindex(["macrophage/monocyte"])[m].sum(),
                       share_T_NK=share.reindex(TNK)[m].sum(), share_lymphoid=share.reindex(LYMPHOID)[m].sum(),
                       share_intended=share.reindex(it)[m].sum() if it else np.nan)
        rows.append(row)
    return rows


def run_variant(ds, an, o, S, mask=None, with_share=True):
    sel = (o["class"] != "(excluded)").values
    if mask is not None:
        sel &= mask
    attr, z, share, cnt, frac = attribution(S[sel].reset_index(drop=True), o.loc[sel, "class"].reset_index(drop=True))
    long = []
    for c in attr.index:
        for m in MODS:
            long.append(dict(dataset=ds, analysis=an, class_=c, n_cells=int(cnt[c]), cell_fraction=frac[c], module=m,
                             mean_score=attr.loc[c, m], z=z.loc[c, m],
                             share_pct=share.loc[c, m] if with_share else np.nan))
    return pd.DataFrame(long).rename(columns={"class_": "class"}), top_rows(ds, an, z, share if with_share else None, attr), sel


def patient_level(ds, an, o, S, sel):
    oo = o[sel].reset_index(drop=True); SS = S[sel].reset_index(drop=True)
    sums = SS.groupby([oo.patient, oo["class"]]).sum()
    n = oo.groupby(["patient", "class"]).size()
    rows, tops = [], []
    for p in sorted(oo.patient.unique()):
        sp = sums.loc[p]; tot = sp.sum()
        for c in sp.index:
            for m in MODS:
                rows.append(dict(dataset=ds, analysis=an, patient=p, n_cells_patient=int(n.loc[p].sum()), **{"class": c},
                                 n_cells_class=int(n.loc[(p, c)]), module=m, share_pct=sp.loc[c, m] / tot[m] * 100))
        # per-patient top class: classes with >= 20 cells, >= 6 classes (icb/attrib_per_patient.py rule)
        ok = [c for c in CLASSES if c in n.loc[p].index and n.loc[(p, c)] >= MIN_CELLS_PT]
        if len(ok) < MIN_TYPES_PT:
            continue
        m_ = (oo.patient == p) & oo["class"].isin(ok)
        attr, z, *_ = attribution(SS[m_].reset_index(drop=True), oo.loc[m_, "class"].reset_index(drop=True))
        for m in MODS:
            it = INTENDED[m]
            tops.append(dict(dataset=ds, analysis=an, patient=p, module=m, n_classes=len(ok),
                             intended_evaluable=bool(it) and any(c in ok for c in it),
                             top_class=z[m].idxmax(), top_z=z[m].max(),
                             top_is_intended=(z[m].idxmax() in it) if it else np.nan))
    return pd.DataFrame(rows), pd.DataFrame(tops)


def summarise_patients(P):
    out = []
    keyset = {"fibroblast": ["fibroblast"], "macrophage/monocyte": ["macrophage/monocyte"], "T+NK": TNK,
              "lymphoid": LYMPHOID, "B+plasma": ["B cell", "plasma cell"]}
    for (ds, an), d in P.groupby(["dataset", "analysis"]):
        for m in MODS:
            dm = d[d.module == m]
            for k, cs in list(keyset.items()) + [(c, [c]) for c in CLASSES if c not in ("fibroblast", "macrophage/monocyte")]:
                v = dm[dm["class"].isin(cs)].groupby("patient").share_pct.sum().reindex(dm.patient.unique()).fillna(0)
                out.append(dict(dataset=ds, analysis=an, module=m, class_group=k, n_patients=len(v),
                                median=v.median(), q1=v.quantile(.25), q3=v.quantile(.75), min=v.min(), max=v.max(),
                                n_ge50=int((v >= 50).sum())))
    return pd.DataFrame(out)


def gene_expression(ds, o, genes, X, sel):
    oo = o[sel].reset_index(drop=True); XX = X[sel]
    rows = []
    for c in [c for c in CLASSES if (oo["class"] == c).any()]:
        mc = (oo["class"] == c).values
        mu = XX[mc].mean(0, dtype=np.float64); pe = (XX[mc] > 0).mean(0)
        for j, g in enumerate(genes):
            rows.append(dict(dataset=ds, gene=g, module=GENE_MODULE[g], **{"class": c}, n_cells=int(mc.sum()),
                             mean_expr=mu[j], pct_expressing=pe[j] * 100))
    E = pd.DataFrame(rows)
    chk = []
    for g in genes:
        e = E[E.gene == g].set_index("class"); m = GENE_MODULE[g]
        it, lin = INTENDED[m], LINEAGE[m]
        top = e.mean_expr.idxmax()
        row = dict(dataset=ds, module=m, gene=g, top_class=top, top_mean=e.mean_expr.max(),
                   top_pct=e.loc[top, "pct_expressing"],
                   intended=";".join(it) if it else "(not prespecified)")
        if it:
            ib = e.loc[[c for c in it if c in e.index], "mean_expr"].idxmax()
            other = e.drop(index=[c for c in it if c in e.index])
            row.update(intended_best=ib, intended_mean=e.loc[ib, "mean_expr"], intended_pct=e.loc[ib, "pct_expressing"],
                       best_other=other.mean_expr.idxmax(), best_other_mean=other.mean_expr.max(),
                       ratio_intended_to_best_other=e.loc[ib, "mean_expr"] / max(other.mean_expr.max(), 1e-9),
                       top_is_intended=top in it, top_in_lineage=top in lin)
        chk.append(row)
    return E, pd.DataFrame(chk)


def reproduction(S_disc, o_disc):
    attr, z, share, cnt, frac = attribution(S_disc, o_disc["class"])
    inv = {"Epi": "epithelial/tumour", "Fibro": "fibroblast", "Endo": "endothelial", "Peri": "pericyte/smooth muscle",
           "Macro": "macrophage/monocyte", "DC": "dendritic cell", "Tcell": "T cell", "NK": "NK cell",
           "Bcell": "B cell", "Plasma": "plasma cell"}
    pub_s = pd.read_csv(DISC_PUB / "cell_attribution_scores.csv", index_col=0).rename(index=inv)[MODS]
    pub_z = pd.read_csv(DISC_PUB / "cell_attribution_zscore.csv", index_col=0).rename(index=inv)[MODS]
    pub_c = pd.read_csv(DISC_PUB / "cell_attribution_contribution_pct.csv", index_col=0).rename(index=inv)[MODS]
    rows = [
        ("z, CAF_ECM in fibroblasts", "2.6", pub_z.loc["fibroblast", "CAF_ECM"], z.loc["fibroblast", "CAF_ECM"]),
        ("z, TAM_SUP in macrophages", "2.8", pub_z.loc["macrophage/monocyte", "TAM_SUP"], z.loc["macrophage/monocyte", "TAM_SUP"]),
        ("share %, fibroblasts of CAF_ECM", "73%", pub_c.loc["fibroblast", "CAF_ECM"], share.loc["fibroblast", "CAF_ECM"]),
        ("share %, lymphoid (T, NK, B, plasma) of T_CELL", "91%", pub_c.reindex(LYMPHOID)["T_CELL"].sum(),
         share.reindex(LYMPHOID)["T_CELL"].sum()),
        ("share %, macrophages of TAM_SUP", "28%", pub_c.loc["macrophage/monocyte", "TAM_SUP"],
         share.loc["macrophage/monocyte", "TAM_SUP"]),
        ("share %, T + NK of T_CELL (prereg grouping)", "-", pub_c.reindex(TNK)["T_CELL"].sum(), share.reindex(TNK)["T_CELL"].sum()),
        ("max |diff|, class mean score table (10 x 6)", "-", 0.0, (attr[MODS] - pub_s.reindex(attr.index)).abs().max().max()),
        ("max |diff|, z table (10 x 6)", "-", 0.0, (z[MODS] - pub_z.reindex(z.index)).abs().max().max()),
        ("max |diff|, share table (10 x 6), %", "-", 0.0, (share[MODS] - pub_c.reindex(share.index)).abs().max().max()),
    ]
    R = pd.DataFrame(rows, columns=["quantity", "published_text", "cell_attribution_py_output", "E17_code"])
    R["abs_diff"] = (R.cell_attribution_py_output - R.E17_code).abs()
    R.loc[R.quantity.str.startswith("max"), "abs_diff"] = R.loc[R.quantity.str.startswith("max"), "E17_code"]
    R["reproduced"] = R.abs_diff < 1e-4
    top = pd.DataFrame(top_rows("GSE183904", "primary", z, share, attr))
    return R, top


def stage_score():
    allA, allT, allP, allPT, allE, allG = [], [], [], [], [], []
    data = {ds: load(ds) for ds in DATASETS}
    disc_genes = data["GSE183904"][1]
    for ds in DATASETS:
        o, genes, harm, dep = data[ds]
        S, used = cell_scores(harm, genes)
        log(f"\n[{ds}] genes used per module: " + "; ".join(f"{m} {len(used[m])}/{len(MODULES[m])}" for m in MODS))
        variants = [("primary", S, None, True)]
        if ds != "GSE183904":
            variants += [("deposited_scale", cell_scores(dep, genes)[0], None, True),
                         ("no_fallback", S, ~o.fallback.values, True),
                         ("genes42_discovery", cell_scores(harm, genes, keep_genes=set(disc_genes))[0], None, True)]
        if ds == "OMIX001073":
            variants += [("drop_Tenriched", S, ~o.tenriched.values, True)]
        if ds == "GSE206785":
            variants += [("SC3P_only", S, (o.platform == "SC3P").values, True),
                         ("SC5P_only", S, (o.platform == "SC5P").values, True)]
        sel0 = (o["class"] != "(excluded)").values
        Sz, _ = cell_scores(harm[sel0], genes, gene_z=True)      # z over the analysed tumour cells
        Sz_full = pd.DataFrame(np.nan, index=range(len(o)), columns=MODS); Sz_full.loc[sel0] = Sz.values
        variants += [("gene_z_scoring", Sz_full, None, False)]
        for an, SS, mask, ws in variants:
            A, T, sel = run_variant(ds, an, o, SS, mask, ws)
            allA.append(A); allT += T
            if an in ("primary", "drop_Tenriched"):
                P, PT = patient_level(ds, an, o, SS, sel)
                allP.append(P); allPT.append(PT)
            if an == "primary":
                E, G = gene_expression(ds, o, genes, harm, sel)
                allE.append(E); allG.append(G)
                if ds == "GSE183904":
                    R, _ = reproduction(SS[sel].reset_index(drop=True), o[sel].reset_index(drop=True))
                    R.to_csv(RES / "E17_reproduction.csv", index=False)
                    log("\nREPRODUCTION of the discovery numbers (cell_attribution.py outputs vs this code):")
                    log(R.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
    A = pd.concat(allA, ignore_index=True); T = pd.DataFrame(allT)
    # post hoc decomposition (not prespecified): new-atlas class means weighted by the DISCOVERY atlas's captured
    # composition, to separate per-cell expression from capture composition. Mast cells have no discovery class: weight 0.
    dfrac = A[(A.dataset == "GSE183904") & (A.analysis == "primary") & (A.module == "CAF_ECM")].set_index("class").cell_fraction
    extra = []
    for ds in ("OMIX001073", "GSE206785"):
        a = A[(A.dataset == ds) & (A.analysis == "primary")]
        idx = [c for c in CLASSES if c in set(a["class"])]
        attr = a.pivot(index="class", columns="module", values="mean_score").reindex(idx)[MODS]
        z = a.pivot(index="class", columns="module", values="z").reindex(idx)[MODS]
        w = dfrac.reindex(idx).fillna(0); w = w / w.sum()
        contrib = attr.mul(w, axis=0); share = contrib / contrib.abs().sum() * 100
        n = a[a.module == "CAF_ECM"].set_index("class").n_cells
        for c in idx:
            for m in MODS:
                extra.append(dict(dataset=ds, analysis="discovery_composition_weights", **{"class": c}, n_cells=int(n[c]),
                                  cell_fraction=w[c], module=m, mean_score=attr.loc[c, m], z=z.loc[c, m],
                                  share_pct=share.loc[c, m]))
        allT += top_rows(ds, "discovery_composition_weights", z, share, attr)
    A = pd.concat([A, pd.DataFrame(extra)], ignore_index=True); T = pd.DataFrame(allT)
    P = pd.concat(allP, ignore_index=True); PT = pd.concat(allPT, ignore_index=True)
    E = pd.concat(allE, ignore_index=True); G = pd.concat(allG, ignore_index=True)
    PS = summarise_patients(P)
    A.to_csv(RES / "E17_attribution_by_class.csv", index=False)
    T.to_csv(RES / "E17_top_class.csv", index=False)
    P.to_csv(RES / "E17_share_by_patient.csv", index=False)
    PS.to_csv(RES / "E17_share_by_patient_summary.csv", index=False)
    PT.to_csv(RES / "E17_top_class_by_patient.csv", index=False)
    E.to_csv(RES / "E17_gene_expression_by_class.csv", index=False)
    G.to_csv(RES / "E17_gene_check.csv", index=False)

    # counts
    C = A[(A.module == "CAF_ECM")][["dataset", "analysis", "class", "n_cells", "cell_fraction"]]
    C.to_csv(RES / "E17_cell_counts.csv", index=False)

    # sensitivity table (prespecified + additional), wide per module
    sens = T.copy()
    sens["kind"] = sens.analysis.map({"primary": "primary", "drop_Tenriched": "prespecified sensitivity",
                                      "deposited_scale": "additional", "no_fallback": "additional",
                                      "genes42_discovery": "additional", "SC3P_only": "additional",
                                      "SC5P_only": "additional", "gene_z_scoring": "additional (ranking only)",
                                      "discovery_composition_weights": "additional (post hoc decomposition)"})
    skip = pd.DataFrame([dict(dataset=d, analysis="bulk_composition_weights", kind="prespecified sensitivity",
                              module="(all)", top_class="NOT RUN",
                              note="no TCGA-STAD cell-type deconvolution (proportions) on disk; not fetched, per instructions")
                         for d in ("OMIX001073", "GSE206785")])
    pd.concat([sens, skip], ignore_index=True).to_csv(RES / "E17_sensitivity.csv", index=False)

    # prespecified expectations
    ex = []
    def tget(ds, an, m):
        r = T[(T.dataset == ds) & (T.analysis == an) & (T.module == m)]
        return r.iloc[0] if len(r) else None
    def pget(ds, an, m, grp):
        r = PS[(PS.dataset == ds) & (PS.analysis == an) & (PS.module == m) & (PS.class_group == grp)]
        return r.iloc[0] if len(r) else None
    def ptop(ds, an, m):
        d = PT[(PT.dataset == ds) & (PT.analysis == an) & (PT.module == m) & PT.intended_evaluable]
        return f"{int(d.top_is_intended.sum())}/{len(d)}" if len(d) else "n/a"
    specs = [
        ("CAF_ECM highest in fibroblasts", "CAF_ECM", "rank", None),
        ("TAM_SUP highest in macrophages/monocytes", "TAM_SUP", "rank", None),
        ("TLS_B highest in B or plasma cells", "TLS_B", "rank", None),
        ("T_CELL highest in T or NK cells", "T_CELL", "rank", None),
        ("fibroblasts contribute >= 50% of CAF_ECM signal", "CAF_ECM", "share_fibroblast", ("ge", "fibroblast")),
        ("T + NK cells contribute >= 50% of T_CELL signal", "T_CELL", "share_T_NK", ("ge", "T+NK")),
        ("macrophages contribute < 50% of TAM_SUP signal", "TAM_SUP", "share_macrophage", ("lt", "macrophage/monocyte")),
    ]
    for ds, an, role in [("OMIX001073", "primary", "primary"), ("GSE206785", "primary", "secondary"),
                         ("OMIX001073", "drop_Tenriched", "prespecified sensitivity")]:
        for lab, m, key, thr in specs:
            r = tget(ds, an, m); d0 = tget("GSE183904", "primary", m)
            if key == "rank":
                held = bool(r.top_is_intended)
                val = f"top {SHORT[r.top_class]} (z = {r.top_z:.2f}); intended z = {r.intended_z:.2f}; next {SHORT[r.second_class]} (z = {r.second_z:.2f})"
                disc = f"top {SHORT[d0.top_class]} (z = {d0.top_z:.2f})"
                pat = f"intended class on top in {ptop(ds, an, m)} patients"
            else:
                v = r[key]; held = bool(v >= 50) if thr[0] == "ge" else bool(v < 50)
                ps = pget(ds, an, m, thr[1])
                val = f"{v:.1f}%"
                disc = f"{d0[key]:.1f}%"
                pat = (f"median {ps['median']:.1f}% (IQR {ps.q1:.1f}-{ps.q3:.1f}); n = {int(ps.n_patients)}; "
                       f"{int(ps.n_ge50)} patients >= 50%")
            ex.append(dict(dataset=ds, analysis=an, role=role, expectation=lab, result=val,
                           held="held" if held else "not held", patient_level=pat, discovery_GSE183904=disc))
    X = pd.DataFrame(ex)
    X.to_csv(RES / "E17_expectations.csv", index=False)
    log("\nTOP CLASS (primary and sensitivity):")
    log(T.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    log("\nEXPECTATIONS:")
    log(X.to_string(index=False))


# ===================================================================== stage 3: figure
def stage_figure():
    sys.path.insert(0, str(FIGDIR))
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.patches import Rectangle
    import bnstyle as bs
    import figstyle as fs
    from bnstyle import INK, GREY, LIGHT, MM, W2
    A = pd.read_csv(RES / "E17_attribution_by_class.csv")
    T = pd.read_csv(RES / "E17_top_class.csv")
    E = pd.read_csv(RES / "E17_gene_expression_by_class.csv")
    PS = pd.read_csv(RES / "E17_share_by_patient_summary.csv")
    FS = 5.2
    DOT = 24
    fig = plt.figure(figsize=(W2, 205 * MM))
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.0, 0.95, 1.15], hspace=0.42, left=0.075, right=0.985,
                  top=0.965, bottom=0.07)

    # ---- (a) z heatmaps
    ga = GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[0], width_ratios=[1, 1, 1, 0.06], wspace=0.12)
    for k, ds in enumerate(DATASETS):
        ax = fig.add_subplot(ga[0, k])
        d = A[(A.dataset == ds) & (A.analysis == "primary")].pivot(index="class", columns="module", values="z")
        d = d.reindex(index=CLASSES, columns=MODS)
        im = ax.imshow(d.values, cmap=fs.DIV, vmin=-3, vmax=3, aspect="auto", interpolation="nearest")
        for i in range(d.shape[0]):
            for j in range(d.shape[1]):
                v = d.values[i, j]
                if np.isnan(v):
                    ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor="white", edgecolor=LIGHT, lw=0.4, hatch="////"))
                    if j == 2:
                        ax.text(2.5, i, "no mast class in this atlas", ha="center", va="center", fontsize=FS - 0.8,
                                color=fs.SLATE, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
                    continue
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=FS - 0.9,
                        color="white" if abs(v) > 1.8 else INK)
        for j, m in enumerate(MODS):
            col = d[m]
            if col.notna().any():
                i = list(d.index).index(col.idxmax())
                ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False, edgecolor=INK, lw=0.9))
        ax.set_xticks(range(len(MODS))); ax.set_xticklabels(MODS, rotation=40, ha="right", fontsize=FS)
        ax.set_yticks(range(len(CLASSES)))
        ax.set_yticklabels([SHORT[c] for c in CLASSES] if k == 0 else [], fontsize=FS)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        n = int(A[(A.dataset == ds) & (A.analysis == "primary") & (A.module == "CAF_ECM")].n_cells.sum())
        npt = int(PS[(PS.dataset == ds) & (PS.analysis == "primary")].n_patients.iloc[0])
        ax.set_title(f"{DSLAB[ds]}\n{n:,} cells, {npt} patients", fontsize=6, loc="left", pad=3)
        if k == 0:
            ax.text(-0.33, 1.13, "a", transform=ax.transAxes, fontsize=8.2, fontweight="bold", va="bottom")
    cax = fig.add_subplot(ga[0, 3])
    cb = fig.colorbar(im, cax=cax); cb.outline.set_visible(False)
    cb.set_label("module score, z across classes", fontsize=FS); cb.ax.tick_params(labelsize=FS - 0.4, length=1.5)

    # ---- (b) share-of-signal stacked bars (11 classes folded into 7 stacked groups; full table in CSV)
    ax = fig.add_subplot(gs[1])
    yt, yl = [], []
    y = 0
    for m in MODS:
        y0 = y
        for ds in DATASETS:
            d = A[(A.dataset == ds) & (A.analysis == "primary") & (A.module == m)].set_index("class").share_pct
            left = 0
            for gname, cs, col in GROUPS:
                w = d.reindex(cs).fillna(0).sum()
                if w > 0:
                    ax.barh(y, w, left=left, height=0.8, color=col, edgecolor="white", lw=0.5)
                left += w
            key = KEYGROUP.get(m)
            if key:
                v = d.reindex(key[1]).fillna(0).sum()
                ps = PS[(PS.dataset == ds) & (PS.analysis == "primary") & (PS.module == m) & (PS.class_group == key[2])]
                txt = f"{key[0]} {v:.0f}%"
                if len(ps):
                    txt += f";  per patient {ps['median'].iloc[0]:.0f}% (IQR {ps.q1.iloc[0]:.0f}–{ps.q3.iloc[0]:.0f})"
                ax.text(101.5, y, txt, fontsize=FS - 0.8, va="center", ha="left", color=fs.SLATE)
            yt.append(y); yl.append(ds)
            y += 1
        ax.annotate(m, xy=(0, y0 + 1), xycoords=("axes fraction", "data"), xytext=(-46, 0),
                    textcoords="offset points", fontsize=FS + 0.2, fontweight="bold", va="center", ha="right")
        y += 0.6
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=FS - 0.7)
    ax.set_ylim(y - 0.6 + 0.0, -0.6)
    ax.set_xlim(0, 100); ax.set_xticks([0, 25, 50, 75, 100])
    ax.axvline(50, color=INK, lw=0.45, ls=(0, (2, 2)), zorder=5)
    ax.set_xlabel("share of module signal, %  (cell fraction × class mean score, normalised over classes)", fontsize=FS + 0.2)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    pos = ax.get_position()
    ax.set_position([0.155, pos.y0, 0.45, pos.height])
    ax.text(-0.3, 1.01, "b", transform=ax.transAxes, fontsize=8.2, fontweight="bold", va="bottom")
    fs.swatch_legend(ax, [(col, g) for g, cs, col in GROUPS], x=1.62, y=0.97, dy=0.07, size=FS - 0.2)
    ax.text(1.62, 0.97 - 7 * 0.07 - 0.02, "dashed line: 50%", transform=ax.transAxes, fontsize=FS - 0.8,
            color=fs.SLATE, va="top")

    # ---- (c) dot plot, OMIX001073 primary
    ax = fig.add_subplot(gs[2])
    ds = "OMIX001073"
    e = E[E.dataset == ds]
    mean = e.pivot(index="class", columns="gene", values="mean_expr").reindex(index=CLASSES, columns=GENES48)
    pct = e.pivot(index="class", columns="gene", values="pct_expressing").reindex(index=CLASSES, columns=GENES48)
    sc = (mean - mean.min()) / (mean.max() - mean.min())             # scaled per gene (0-1)
    xs, ys = np.meshgrid(np.arange(len(GENES48)), np.arange(len(CLASSES)))
    ax.scatter(xs.ravel(), ys.ravel(), s=(pct.values.ravel() / 100) * DOT, c=sc.values.ravel(), cmap=fs.SEQ,
               vmin=0, vmax=1, lw=0.2, edgecolors=LIGHT)
    for m in MODS:
        js = [GENES48.index(g) for g in MODULES[m]]
        for c in INTENDED[m] or []:
            i = CLASSES.index(c)
            ax.add_patch(Rectangle((min(js) - .5, i - .5), len(js), 1, fill=False, edgecolor=bs.STROMA if m == "CAF_ECM"
                                   else bs.MACRO if m == "TAM_SUP" else bs.IMMUNE, lw=0.6))
        ax.text(np.mean(js), -1.25, m, ha="center", va="bottom", fontsize=FS, fontweight="bold")
        if m != MODS[-1]:
            ax.axvline(max(js) + .5, color=LIGHT, lw=0.5)
    ax.set_xticks(range(len(GENES48))); ax.set_xticklabels(GENES48, rotation=90, fontsize=FS - 0.8)
    ax.set_yticks(range(len(CLASSES))); ax.set_yticklabels([SHORT[c] for c in CLASSES], fontsize=FS)
    ax.set_xlim(-.7, len(GENES48) - .3); ax.set_ylim(len(CLASSES) - .5, -1.6)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.text(-0.065, 1.11, "c", transform=ax.transAxes, fontsize=8.2, fontweight="bold", va="bottom")
    ax.text(0, 1.11, "Module-gene expression by class, OMIX001073 tumour cells (boxes: prespecified intended class)",
            transform=ax.transAxes, fontsize=6, va="bottom")
    # size and colour keys
    for k, p in enumerate([10, 50, 100]):
        ax.scatter([1.005], [0.80 - k * 0.12], s=(p / 100) * DOT, c="#9AA5B1", transform=ax.transAxes, clip_on=False)
        ax.text(1.017, 0.80 - k * 0.12, f"{p}%", transform=ax.transAxes, fontsize=FS - 0.8, va="center")
    ax.text(1.0, 0.93, "% expressing", transform=ax.transAxes, fontsize=FS - 0.6, va="center")
    cax2 = fig.add_axes([ax.get_position().x1 + 0.004, ax.get_position().y0 + 0.02, 0.008, 0.07])
    import matplotlib as mpl
    cb2 = fig.colorbar(mpl.cm.ScalarMappable(cmap=fs.SEQ, norm=mpl.colors.Normalize(0, 1)), cax=cax2)
    cb2.outline.set_visible(False); cb2.ax.tick_params(labelsize=FS - 1, length=1.2)
    cb2.set_label("mean, scaled\nper gene", fontsize=FS - 1)
    bs.save(fig, "Fig_E17_attribution")
    log(f"wrote {FIGDIR/'Fig_E17_attribution'}.{{png,pdf}}")


# stacked groups for panel b (adjacent pairs pass the dataviz CVD / normal-vision checks, figures/check_palette.py)
GROUPS = [("T / NK", ["T cell", "NK cell"], "#0F6FC5"), ("B / plasma", ["B cell", "plasma cell"], "#86B8E8"),
          ("Mac/Mono", ["macrophage/monocyte"], "#C96F12"), ("DC / mast", ["dendritic cell", "mast cell"], "#EDC07A"),
          ("Fibro", ["fibroblast"], "#138A5E"), ("Endo / Peri/SMC", ["endothelial", "pericyte/smooth muscle"], "#93CDB0"),
          ("Epi", ["epithelial/tumour"], "#5E6B78")]
KEYGROUP = {"CAF_ECM": ("Fibro", ["fibroblast"], "fibroblast"),
            "TAM_SUP": ("Mac/Mono", ["macrophage/monocyte"], "macrophage/monocyte"),
            "TLS_B": ("B/plasma", ["B cell", "plasma cell"], "B+plasma"),
            "T_CELL": ("T/NK", ["T cell", "NK cell"], "T+NK")}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["extract", "score", "figure", "all"], default="all")
    a = ap.parse_args()
    LOG = RES / "E17_log.txt"
    log(f"\n===== E17 stage={a.stage}  {time.strftime('%Y-%m-%d %H:%M:%S')} =====")
    if a.stage in ("extract", "all"):
        stage_extract()
    if a.stage in ("score", "all"):
        stage_score()
    if a.stage in ("figure", "all"):
        stage_figure()
