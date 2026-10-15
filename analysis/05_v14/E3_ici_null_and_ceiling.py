"""E3 - PRJEB25780 (pembrolizumab, n=44, 12 responders).

(a) Structure-matched random-axis null for the ICI AUCs of NetImmune, Immune and IFN:
    random modules of the same sizes drawn from the 14,364-gene bulk pool, scored by the
    identical row-z mean. Random axes are credited with their better direction
    (AUC -> max(AUC, 1-AUC)); the observed scores keep their a-priori direction.
    SCALE (corrected after v15): the deposited matrix mGC_Pembrolizumab_Kim2018.TPM is
    already log2(TPM + 1) (max 16.8; sum of 2^x - 1 per sample ~ 1e6), so it is used as
    provided, which is also what the locked script (benchmark_barriernet_official.py) did.
    v15 wrongly took the file for linear TPM and applied log2(x + 1) a second time; that
    double-logged version is kept as a labelled sensitivity ("log2x2") for transparency.
(b) Re-run of the full-transcriptome ceiling (revision2/scripts/S12b.py), whose log was
    not kept: L2 logistic regression, 40 x repeated stratified 5-fold CV.
"""
import warnings, functools; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

V = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
NRAND = 2000
rng = np.random.default_rng(2026)
MOD = dict(
    CAF_ECM=["COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"],
    TAM_SUP=["SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"],
    APM=["B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"],
    IFN=["STAT1","IRF1","CXCL9","CXCL10","IFNG"],
    TLS_B=["CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"],
    T_CELL=["CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"])
IMM, BAR = ["APM","IFN","TLS_B","T_CELL"], ["CAF_ECM","TAM_SUP"]

e = pd.read_csv(f"{V}/mGC_Pembrolizumab_Kim2018.TPM", sep="\t", index_col=0)
b = pd.read_csv(f"{V}/benchmark_scores.tsv", sep="\t", index_col=0)
e0 = e[b.index]
y = b["Responder"].values.astype(int)
print(f"PRJEB25780: {e0.shape[1]} patients, {y.sum()} responders")
pool_all = open(f"{O}/bulk_gene_pool.txt").read().split()
si = [len(MOD[m]) for m in IMM]; sb = [len(MOD[m]) for m in BAR]
null_rows = []
for scale in ["log2", "log2x2"]:        # log2 = as deposited (primary); log2x2 = v15 error, sensitivity
    e = e0 if scale == "log2" else np.log2(e0 + 1)
    Z = ((e.T - e.mean(1)) / e.std(1, ddof=0).replace(0, 1)).T

    # reproduce the published scores from genes, as a check on the scoring path
    mod = pd.DataFrame({m: Z.loc[[g for g in G if g in Z.index]].mean() for m, G in MOD.items()})
    net = mod[IMM].mean(1) - mod[BAR].mean(1)
    print(scale, "max |recomputed - published| NetImmune:", float(np.abs(net - b["NetImmune"]).max()))
    obs = {"NetImmune": roc_auc_score(y, net), "Immune": roc_auc_score(y, mod[IMM].mean(1)),
           "IFN": roc_auc_score(y, mod["IFN"])}
    print("observed AUC:", {k: round(v, 3) for k, v in obs.items()})

    pool = [g for g in pool_all if g in Z.index]
    Zp = Z.loc[pool].values
    print("null gene pool (bulk pool present in PRJEB25780):", len(pool))

    def score(idx, sizes):
        out, k = [], 0
        for s in sizes:
            out.append(Zp[idx[k:k + s]].mean(0)); k += s
        return np.mean(out, 0)

    rows = []
    for i in range(NRAND):
        g = rng.choice(len(pool), sum(si) + sum(sb), replace=False)
        imm = score(g[:sum(si)], si); bar = score(g[sum(si):], sb)
        ifn = Zp[rng.choice(len(pool), 5, replace=False)].mean(0)
        a = [roc_auc_score(y, imm - bar), roc_auc_score(y, imm), roc_auc_score(y, ifn)]
        rows.append([max(x, 1 - x) for x in a])
    R = pd.DataFrame(rows, columns=["NetImmune", "Immune", "IFN"])
    R.to_csv(f"{O}/E3_ici_null_draws_{scale}.csv", index=False)
    null_rows += [{"scale": scale, "score": k, "AUC_obs": obs[k], "null_median": R[k].median(),
                          "null_p95": R[k].quantile(0.95),
                          "p_emp": ((R[k] >= obs[k]).sum() + 1) / (NRAND + 1)} for k in obs]
null = pd.DataFrame(null_rows)
null.to_csv(f"{O}/E3_ici_matched_null.csv", index=False)
print("\n=== matched random-axis null (direction credited to random axes) ===")
print(null.round(4).to_string(index=False))

e = e0[np.isfinite(e0).all(1)]                      # already log2(TPM + 1); do not log again
Z = ((e.T - e.mean(1)) / e.std(1).replace(0, 1)).T
X = np.ascontiguousarray(Z.T.values)
def rep_cv(C, rep=40, k=5):
    a = []
    for i in range(rep):
        p = np.zeros(len(y))
        for tr, te in StratifiedKFold(k, shuffle=True, random_state=i).split(X, y):
            m = LogisticRegression(solver="liblinear", dual=True, C=C, max_iter=2000).fit(X[tr], y[tr])
            p[te] = m.predict_proba(X[te])[:, 1]
        a.append(roc_auc_score(y, p))
    return np.array(a)
print(f"\n=== ICI ceiling: full transcriptome ({X.shape[1]} genes, log2 TPM as deposited), 40 x 5-fold CV ===")
cl = []
for C in [1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
    s = rep_cv(C)
    cl.append({"C": C, "cv_auc": s.mean(), "lo": np.percentile(s, 2.5), "hi": np.percentile(s, 97.5)})
    print(f"  C={C:<7g} CV AUC = {s.mean():.3f} [{np.percentile(s, 2.5):.3f}-{np.percentile(s, 97.5):.3f}]")
pd.DataFrame(cl).to_csv(f"{O}/E3_ici_ceiling.csv", index=False)
