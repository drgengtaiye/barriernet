import warnings, functools; warnings.filterwarnings('ignore')
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

V = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
e = pd.read_csv(f"{V}/mGC_Pembrolizumab_Kim2018.TPM", sep="\t", index_col=0)
e = np.log2(e + 1); e = e[np.isfinite(e).all(1)]
Z = ((e.T - e.mean(1)) / e.std(1).replace(0, 1)).T
b = pd.read_csv(f"{V}/benchmark_scores.tsv", sep="\t", index_col=0)
X = Z[b.index].T.values; y = b['Responder'].values.astype(int)
X = np.ascontiguousarray(X[:, np.isfinite(X).all(0)])
print("ICI full transcriptome:", X.shape, "| responders", y.sum())

def rep_cv(C, rep=40, k=5):
    a = []
    for i in range(rep):
        p = np.zeros(len(y))
        for tr, te in StratifiedKFold(k, shuffle=True, random_state=i).split(X, y):
            m = LogisticRegression(solver='liblinear', dual=True, C=C, max_iter=2000).fit(X[tr], y[tr])
            p[te] = m.predict_proba(X[te])[:, 1]
        a.append(roc_auc_score(y, p))
    return np.array(a)

print("\n=== ICI CEILING: full transcriptome (L2 logistic, 40x repeated 5-fold CV) ===")
best = None
for C in [1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
    s = rep_cv(C)
    print(f"  full transcriptome C={C:<7g} CV AUC = {s.mean():.3f} [{np.percentile(s,2.5):.3f}-{np.percentile(s,97.5):.3f}]")
    if best is None or s.mean() > best[0]: best = (s.mean(), C)
print(f"\n  BEST full-transcriptome CV AUC = {best[0]:.3f}  (C={best[1]})")
for nm, c in [("IFN (5 genes, fixed)", "IFN"), ("Immune (24 genes, fixed)", "Immune"),
              ("NetImmune (48 genes, fixed)", "NetImmune")]:
    print(f"  {nm:30s} AUC = {roc_auc_score(y, b[c]):.3f}  (no fitting)")
