"""E5 (PRJEB25780) - response AUC of the TAM_SUP arms and the CXCL9:SPP1 polarity score in the
gastric pembrolizumab cohort (44 patients, 12 responders). The deposited matrix is already
log2(TPM + 1) and is used as provided ("log2", primary); "log2x2" re-applies log2(x + 1), the v15
error, kept as a labelled sensitivity. Peer-review addition;
expectations as in E5_trials.py (C1Q arm and polarity higher in responders; SPP1 arm not).
AUC > 0.5 means higher in responders; 2,000-resample stratified bootstrap 95% CI."""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
MOD = dict(CAF_ECM=["COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"],
           TAM_SUP=["SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"],
           C1Q_arm=["C1QA","C1QB","C1QC"], SPP1_arm=["SPP1","TREM2","MARCO","APOE","LGALS3"])
e0 = pd.read_csv(f"{V}/mGC_Pembrolizumab_Kim2018.TPM", sep="\t", index_col=0)
b = pd.read_csv(f"{V}/benchmark_scores.tsv", sep="\t", index_col=0); y = b.Responder.values.astype(int)
rng = np.random.default_rng(2026); rows = []
for scale in ["log2", "log2x2"]:
    e = e0[b.index] if scale == "log2" else np.log2(e0[b.index] + 1)
    Z = ((e.T - e.mean(1)) / e.std(1, ddof=0).replace(0, 1)).T
    S = {k: Z.loc[G].mean() for k, G in MOD.items()}; S["CXCL9_SPP1"] = Z.loc["CXCL9"] - Z.loc["SPP1"]
    for k, v in S.items():
        v = v.values; a = roc_auc_score(y, v); bs = []
        i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
        for _ in range(2000):
            j = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]; bs.append(roc_auc_score(y[j], v[j]))
        rows.append(dict(scale=scale, score=k, AUC=a, lo=np.percentile(bs, 2.5), hi=np.percentile(bs, 97.5)))
R = pd.DataFrame(rows); R.to_csv(f"{O}/E5_prjeb_arms.csv", index=False); print(R.round(3).to_string(index=False))
