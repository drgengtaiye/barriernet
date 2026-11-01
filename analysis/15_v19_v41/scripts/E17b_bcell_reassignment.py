"""Sensitivity analysis for the discovery-atlas attribution (Fig. 2, v19).

The independent figure review found cells annotated 'B cell' inside the macrophage cluster of
GSE183904. Many express myeloid markers and no B-cell marker, which points to an annotation error
in the processed atlas. This script reassigns B-labelled cells that express at least one myeloid
marker and neither B-cell marker to macrophage/monocyte, and recomputes each module's share of
signal and class z-score.

  myeloid marker: CD68, LYZ, C1QA, FCER1G or LST1 > 0
  B-cell marker:  MS4A1 or CD79A > 0

Descriptive only; the rule was written after seeing the UMAP. Per-cell module scores are those
plotted in Fig. 2b (figures/cns/data/fig2_cell_module_scores.parquet, computed as in E17).
Writes results/E17b_bcell_reassignment.csv.
"""
import pathlib
import numpy as np
import pandas as pd

G = pathlib.Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2")
obs = pd.read_parquet(G / "figures/data/sc/cells_obs.parquet")
expr = pd.read_parquet(G / "figures/data/sc/cells_expr.parquet")
sc = pd.read_parquet(G / "figures/cns/data/fig2_cell_module_scores.parquet")
assert (sc.index == obs.index).all() and (expr.index == obs.index).all()

lab = obs.anno_level1.astype(str)
myeloid = (expr[["CD68", "LYZ", "C1QA", "FCER1G", "LST1"]] > 0).any(axis=1)
bmark = (expr[["MS4A1", "CD79A"]] > 0).any(axis=1)
move = (lab == "Bcell") & myeloid & ~bmark
lab2 = lab.where(~move, "Macro")

rows = []
for mod in sc.columns:
    for name, L in [("as annotated", lab), ("reassigned", lab2)]:
        tot = sc[mod].groupby(L.values).sum()
        mean = sc[mod].groupby(L.values).mean()
        z = (mean - mean.mean()) / mean.std()
        for cls in tot.index:
            rows.append({"module": mod, "labels": name, "class": cls, "n_cells": int((L == cls).sum()),
                         "share_pct": 100 * tot[cls] / tot.sum(), "class_mean": mean[cls], "z": z[cls]})
out = pd.DataFrame(rows)
out.to_csv(G / "analysis/results/E17b_bcell_reassignment.csv", index=False)

print(f"B-labelled cells: {(lab == 'Bcell').sum():,}; reassigned (myeloid marker, no B marker): {int(move.sum()):,}")
print("marker detection in reassigned cells:",
      {g: round(float((expr.loc[move, g] > 0).mean()), 2) for g in ["MS4A1", "CD79A", "CD68", "LYZ", "C1QA"]})
for mod, cls in [("TAM_SUP", ["Macro"]), ("TLS_B", ["Bcell", "Plasma"]), ("CAF_ECM", ["Fibro"]), ("T_CELL", ["Tcell", "NK"])]:
    d = out[(out.module == mod) & out["class"].isin(cls)].groupby("labels").share_pct.sum()
    top = out[out.module == mod].sort_values("z").groupby("labels").tail(1).set_index("labels")["class"]
    print(f"{mod}: share in {'+'.join(cls)} {d['as annotated']:.1f}% -> {d['reassigned']:.1f}%; "
          f"top class {top['as annotated']} -> {top['reassigned']}")
