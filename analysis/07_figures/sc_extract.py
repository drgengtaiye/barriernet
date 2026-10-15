"""Pull what Figure 2 and Supplementary Figs S1-S4 need out of the 5 GB GSE183904 object:
per-cell obs (cell type, sample, programme scores), UMAP, and log-normalised expression of
the canonical markers and the 48 barriernet module genes; plus the fibroblast and
macrophage sub-objects (UMAP, subtype/state, programme scores) and per-sample clinical data."""
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/data/sc"
MOD = dict(CAF_ECM=["COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"],
           TAM_SUP=["SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"],
           APM=["B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"], IFN=["STAT1","IRF1","CXCL9","CXCL10","IFNG"],
           TLS_B=["CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"], T_CELL=["CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"])
CANON = ["CD3D","CD3E","NKG7","GNLY","MS4A1","CD79A","MZB1","JCHAIN","LYZ","CD68","C1QA","LST1","FCER1G",
         "CLEC10A","COL1A1","DCN","PECAM1","VWF","RGS5","MCAM","EPCAM","KRT19","KRT8"]
a = ad.read_h5ad(f"{R}/out_scanpy_strict_GSE183904/GSE183904_TumorOnly_scanpy_strict_final_scored.h5ad", backed="r")
obs = a.obs[["sample_id", "tissue", "anno_level1", "score_TLS", "score_CYTO", "score_EXHAUST", "score_IFN",
             "score_CAF_ECM", "score_TGFb", "score_TAM_SUP"]].copy()
obs["umap1"], obs["umap2"] = a.obsm["X_umap"][:, 0], a.obsm["X_umap"][:, 1]
genes = [g for g in dict.fromkeys(CANON + sum(MOD.values(), [])) if g in a.var_names]
missing = [g for g in dict.fromkeys(CANON + sum(MOD.values(), [])) if g not in a.var_names]
print("missing genes:", missing)
idx = [int(np.where(a.var_names == g)[0][0]) for g in genes]
X = a.X[:, sorted(idx)]
X = X.toarray() if sp.issparse(X) else np.asarray(X)
order = np.argsort(np.argsort(idx)); cols = [genes[i] for i in np.argsort(idx)]
E = pd.DataFrame(X, index=a.obs_names, columns=cols)[genes]
print("X max per gene (log-normalised check):", E.max().round(2).head(5).to_dict())
obs.to_parquet(f"{O}/cells_obs.parquet"); E.astype("float32").to_parquet(f"{O}/cells_expr.parquet")
for name, f, lab, progs in [("fibro", "results/Fig5_Fibro/Fibro_only_annotated.h5ad", "fibro_subtype",
                             ["prog_CAF_ECM", "prog_iCAF", "prog_STF3_FAP_INHBA", "prog_myCAF", "prog_CAF_TGFB"]),
                            ("macro", "results/Fig5_Macro/Macro_only_annotated.h5ad", "macro_state",
                             ["prog_SPP1_TAM", "prog_Inflam", "prog_HLAII_APC", "prog_IL10_STAT3", "prog_IFN"])]:
    s = ad.read_h5ad(f"{R}/{f}", backed="r")
    o = s.obs[["sample_id", "Patient_ID", lab] + progs].copy()
    o["umap1"], o["umap2"] = s.obsm["X_umap"][:, 0], s.obsm["X_umap"][:, 1]
    o.to_parquet(f"{O}/{name}_obs.parquet"); print(name, o.shape, o[lab].value_counts().to_dict())
    if name == "fibro":
        cl = s.obs[["sample_id", "Patient_ID", "sample_type", "Age", "Gender", "Stage", "Laurens", "TCGA", "MMR", "EBV", "StageGroup"]]
        cl.astype(str).drop_duplicates("sample_id").to_csv(f"{O}/sample_clinical_fibro.csv", index=False)
    else:
        cl = s.obs[["sample_id", "Patient_ID", "sample_type", "Age", "Gender", "Stage", "Laurens", "TCGA", "MMR", "EBV", "StageGroup"]]
        cl.astype(str).drop_duplicates("sample_id").to_csv(f"{O}/sample_clinical_macro.csv", index=False)
print("done")
