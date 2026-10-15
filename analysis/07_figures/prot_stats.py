"""Protein-level statistics for Fig. 7 (PDC000214 iTRAQ tumour/adjacent log ratios, 130 paired
samples, 80 patients), using the barriernet 48-gene module definitions."""
import numpy as np, pandas as pd
R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
MOD = dict(CAF_ECM=["COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"],
 TAM_SUP=["SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"],
 APM=["B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"], IFN=["STAT1","IRF1","CXCL9","CXCL10","IFNG"],
 TLS_B=["CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"], T_CELL=["CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"])
E = pd.read_csv(f"{R}/expr.tsv", sep="\t", index_col=0)
E = E.drop(columns=[c for c in ["Mean", "Median", "StdDev"] if c in E.columns])
det = {g: int(E[g].notna().sum()) if g in E.columns else 0 for G in MOD.values() for g in G}
Z = (E - E.mean()) / E.std()
used = {m: [g for g in G if det[g] >= 30] for m, G in MOD.items()}
score = pd.DataFrame({m: Z[G].mean(1, skipna=True) for m, G in used.items()})

def coherence(genes):
    c = Z[genes].corr(method="spearman", min_periods=30).values
    return np.nanmedian(c[np.triu_indices(len(genes), 1)])

rng = np.random.default_rng(2026)
pool = [g for g in E.columns if E[g].notna().sum() >= 100]
rows = []
for m, G in used.items():
    obs = coherence(G)
    null = np.array([coherence(list(rng.choice(pool, len(G), replace=False))) for _ in range(500)])
    rows.append(dict(module=m, n_used=len(G), n_defined=len(MOD[m]), rho_obs=obs, null_median=np.median(null),
                     null_p95=np.quantile(null, .95), p_emp=((null >= obs).sum() + 1) / 501))
coh = pd.DataFrame(rows)
mm = score.corr(method="spearman")
if __name__ == "__main__":
    print(pd.Series(det).to_string())
    print(coh.round(3).to_string(index=False))
    print(mm.round(2))
    coh.to_csv("Fig7_module_coherence.csv", index=False); mm.to_csv("Fig7_module_correlation.csv")
