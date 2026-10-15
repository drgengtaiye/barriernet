"""E7a - PRJEB25780 (pembrolizumab, Kim 2018) scores on log2(TPM + 1).

SCALE (corrected after v15): mGC_Pembrolizumab_Kim2018.TPM is already log2(TPM + 1) despite its
name (max 16.8; sum of 2^x - 1 per sample ~ 1e6). v15 applied log2(x + 1) again; this version
uses the matrix as provided. The locked script also used it as provided, so its module scores
were on the right scale; its comparators, however, were not the published methods:
  GEP-18  unweighted z-mean of the 18 genes;
  IPS     z-mean(MHC 8) + z-mean(EC 6) - z-mean(SC 5) - z-mean(CP 2), a crude
          stand-in for the 162-gene Charoentong algorithm;
  TIDE    tidepy was not installed, so the script fell back to an 'exclusion proxy'
          (mean of two 8-gene CAF and TAM z-means). The published 'TIDE' column equals
          that proxy to 1e-9. TIDE itself cannot be run offline here.
This script scores on log2(TPM + 1) (the matrix as deposited):
  barriernet modules/axes  row-z mean over the 44 evaluable patients, >= 3 genes (as E3);
  GEP-18   mean log2 expression of the 18 Ayers genes with equal weights, as in the
           PredictIO/SignatureSets implementation (Tcell-inflamed_Ayers, weight 1); the
           Cristescu 2018 weights are not available offline. A z-mean variant is kept too;
  IPS      Charoentong 2017 algorithm (IPS_Charoentong.csv from SignatureSets: 162 genes,
           26 groups): per-sample z over all genes, group means x group weight, MHC/CP/EC/SC
           class means, raw aggregate AZ = sum; IPS = 0..10 map. AZ is used for AUC;
  PD-L1    CD274 on the deposited log2(TPM + 1) scale;
  exclusion proxy  the locked script's TIDE fallback, recomputed on log2 and labelled as
           what it is.
With --twice, everything is recomputed on log2(x + 1) of the matrix (the v15 error) and
written to E7_scores_log2x2.csv, as a labelled sensitivity only.
"""
import numpy as np, pandas as pd

V = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
SIG = "/groups/adv2105_gp/yichen/Yi/com/gastric/icb/SignatureSets-main/data-raw/IO"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"

MOD = dict(
    CAF_ECM=["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
    TAM_SUP=["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
    APM=["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    IFN=["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
    TLS_B=["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
    T_CELL=["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"])
GEP18 = pd.read_csv(f"{SIG}/Tcell-inflamed_Ayers.csv").gene_name.tolist()
EXCL = (["COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "FAP", "ACTA2", "PDGFRB"],
        ["CD163", "MSR1", "MRC1", "IL10", "TGFB1", "SPP1", "APOE", "C1QA"])

b = pd.read_csv(f"{V}/benchmark_scores.tsv", sep="\t", index_col=0)
tpm = pd.read_csv(f"{V}/mGC_Pembrolizumab_Kim2018.TPM", sep="\t", index_col=0)[b.index]
import sys
TWICE = "--twice" in sys.argv
L = np.log2(tpm + 1) if TWICE else tpm                # file is already log2(TPM + 1)
L = L[np.isfinite(L).all(1)]
Z = ((L.T - L.mean(1)) / L.std(1, ddof=0).replace(0, 1)).T      # gene z across patients

def zmean(genes, need=3):
    g = [x for x in genes if x in Z.index]
    miss = [x for x in genes if x not in Z.index]
    return (Z.loc[g].mean() if len(g) >= need else pd.Series(np.nan, index=Z.columns)), g, miss

out = pd.DataFrame(index=b.index)
log = []
for m, G in MOD.items():
    out[m], g, miss = zmean(G); log.append(f"{m}: {len(g)}/{len(G)} {('missing ' + ','.join(miss)) if miss else ''}")
out["Barrier"] = out[["CAF_ECM", "TAM_SUP"]].mean(1)
out["Immune"] = out[["APM", "IFN", "TLS_B", "T_CELL"]].mean(1)
out["NetImmune"] = out.Immune - out.Barrier

g = [x for x in GEP18 if x in L.index]
log.append(f"GEP18: {len(g)}/{len(GEP18)} missing {[x for x in GEP18 if x not in L.index]}")
out["GEP18"] = L.loc[g].mean()                       # equal-weight mean of log2 expression
out["GEP18_zmean"] = zmean(GEP18)[0]

ips = pd.read_csv(f"{SIG}/IPS_Charoentong.csv")
log.append(f"IPS: {ips.gene_name.isin(L.index).sum()}/{len(ips)} genes present; missing "
           f"{ips.gene_name[~ips.gene_name.isin(L.index)].tolist()}")
ips = ips[ips.gene_name.isin(L.index)]
az, ipsd = {}, {}
for s in L.columns:
    ge = L[s]
    z1 = (ge.loc[ips.gene_name].values - ge.mean()) / ge.std()      # per-sample z over all genes
    t = ips.assign(z=z1).groupby(["class", "name"]).agg(z=("z", "mean"), w=("weight", "mean"))
    wg = (t.z * t.w).groupby(level="class").mean()
    a = wg[["MHC", "CP", "EC", "SC"]].sum()
    az[s] = a
    ipsd[s] = 0 if a <= 0 else (10 if a >= 3 else round(a * 10 / 3))
out["IPS_AZ"] = pd.Series(az); out["IPS"] = pd.Series(ipsd)

out["PDL1"] = L.loc["CD274"]                        # log2(TPM + 1)
out["ExclusionProxy"] = (zmean(EXCL[0], 2)[0] + zmean(EXCL[1], 2)[0]) / 2
out["MSI_H"] = b.MSI_H.astype(str).isin(["True", "TRUE"]).astype(int)
out["EBV"] = b.EBV_pos.astype(str).isin(["True", "TRUE"]).astype(int)
out["MSI_EBV"] = ((out.MSI_H + out.EBV) > 0).astype(int)
out["Responder"] = b.Responder.astype(int)
tag = "log2x2" if TWICE else "log2"
out.to_csv(f"{O}/E7_scores_{tag}.csv")
open(f"{O}/E7_scores_{tag}.log", "w").write("\n".join(log) + "\n")
print("\n".join(log))
from sklearn.metrics import roc_auc_score
print(f"check NetImmune AUC ({tag}; single log 0.797, double log 0.815):", round(roc_auc_score(out.Responder, out.NetImmune), 3))
