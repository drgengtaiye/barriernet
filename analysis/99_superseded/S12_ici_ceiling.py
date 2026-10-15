import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
V="/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
O="results"
e=pd.read_csv(f"{V}/mGC_Pembrolizumab_Kim2018.TPM",sep="\t",index_col=0)
e=np.log2(e+1); e=e[np.isfinite(e).all(1)]
Z=((e.T-e.mean(1))/e.std(1).replace(0,1)).T          # gene-wise z
b=pd.read_csv(f"{V}/benchmark_scores.tsv",sep="\t",index_col=0)
X=Z[b.index].T.values; y=b['Responder'].values
X=X[:,np.isfinite(X).all(0)]
print("ICI full transcriptome:",X.shape,"responders",y.sum())

def rep_cv(Xm,C,rep=100,k=5):
    a=[]
    for i in range(rep):
        p=np.zeros(len(y))
        for tr,te in StratifiedKFold(k,shuffle=True,random_state=i).split(Xm,y):
            m=LogisticRegression(penalty='l2',C=C,max_iter=3000).fit(Xm[tr],y[tr])
            p[te]=m.predict_proba(Xm[te])[:,1]
        a.append(roc_auc_score(y,p))
    return np.array(a)

print("\n=== ICI CEILING: full transcriptome (L2 logistic, repeated 5-fold CV) ===")
best=None
for C in [1e-4,1e-3,1e-2,1e-1]:
    s=rep_cv(X,C); print(f"  full transcriptome  C={C:<8g} CV AUC = {s.mean():.3f} [{np.percentile(s,2.5):.3f}-{np.percentile(s,97.5):.3f}]")
    if best is None or s.mean()>best[0]: best=(s.mean(),C)
print(f"\n  best full-transcriptome CV AUC = {best[0]:.3f}")
for nm,col in [("Immune (24 genes, fixed)","Immune"),("NetImmune (48 genes, fixed)","NetImmune"),("IFN (5 genes, fixed)","IFN")]:
    print(f"  {nm:28s} AUC = {roc_auc_score(y,b[col]):.3f}   (no fitting, no CV needed)")
