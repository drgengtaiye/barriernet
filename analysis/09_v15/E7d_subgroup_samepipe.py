"""E7d - (i) responders by molecular subgroup and NetImmune AUC within MSS/EBV-negative
patients; (ii) fair comparison with the full-transcriptome classifier: single scores passed
through exactly the E3 pipeline (L2 logistic regression, liblinear dual, C in 1e-4..1,
40 repetitions of stratified 5-fold CV with random_state = repetition, pooled out-of-fold AUC
per repetition, mean across repetitions; best C reported, as for the classifier)."""
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
d = pd.read_csv(f"{O}/E7_scores_log2.csv", index_col=0)
y = d.Responder.values.astype(int)
rng = np.random.default_rng(2026)

# ---- (i) subgroup
d["group"] = np.where(d.MSI_H == 1, "MSI-H", np.where(d.EBV == 1, "EBV+", "MSS/EBV-negative"))
ct = pd.crosstab(d.group, d.Responder.map({1: "responder", 0: "non-responder"}))
print(ct)
sub = d[d.group == "MSS/EBV-negative"]
def boot_auc(yv, s, n=2000):
    a = roc_auc_score(yv, s); pos, neg = np.where(yv == 1)[0], np.where(yv == 0)[0]; bs = []
    for _ in range(n):
        i = np.r_[rng.choice(pos, len(pos)), rng.choice(neg, len(neg))]
        bs.append(roc_auc_score(yv[i], s[i]))
    return a, np.percentile(bs, 2.5), np.percentile(bs, 97.5)
srows = []
for k in ["NetImmune", "Immune", "IFN", "GEP18", "PDL1"]:
    if sub.Responder.nunique() == 2:
        a, lo, hi = boot_auc(sub.Responder.values, sub[k].values)
        srows.append(dict(score=k, n=len(sub), responders=int(sub.Responder.sum()), AUC=a, lo=lo, hi=hi))
S = pd.DataFrame(srows); print(S.round(3).to_string(index=False))
ct.to_csv(f"{O}/E7_subgroup_counts.csv"); S.to_csv(f"{O}/E7_subgroup_auc.csv", index=False)

# ---- (ii) same pipeline as the E3 full-transcriptome classifier
def rep_cv(X, C, rep=40, k=5):
    a = []
    for i in range(rep):
        p = np.zeros(len(y))
        for tr, te in StratifiedKFold(k, shuffle=True, random_state=i).split(X, y):
            m = LogisticRegression(solver="liblinear", dual=True, C=C, max_iter=2000).fit(X[tr], y[tr])
            p[te] = m.predict_proba(X[te])[:, 1]
        a.append(roc_auc_score(y, p))
    return np.array(a)
rows = []
for k in ["NetImmune", "Immune", "IFN"]:
    x = d[k].values; X = ((x - x.mean()) / x.std())[:, None]    # z across the 44, as genes in E3
    for C in [1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
        s = rep_cv(X, C)
        rows.append(dict(score=k, C=C, cv_auc=s.mean(), lo=np.percentile(s, 2.5), hi=np.percentile(s, 97.5),
                         fixed_auc=roc_auc_score(y, x)))
R = pd.DataFrame(rows); R.to_csv(f"{O}/E7_samepipeline_cv.csv", index=False)
best = R.loc[R.groupby("score").cv_auc.idxmax()]
print(R.round(3).to_string(index=False)); print("\nbest C per score:"); print(best.round(3).to_string(index=False))
ceil = pd.read_csv(f"{O}/E3_ici_ceiling.csv"); print("\nfull-transcriptome classifier (E3):"); print(ceil.round(3).to_string(index=False))
