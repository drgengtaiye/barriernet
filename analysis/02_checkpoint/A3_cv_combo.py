import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import statsmodels.api as sm
rng=np.random.default_rng(42)

V="/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
O="/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
d=pd.read_csv(f"{V}/benchmark_scores.tsv",sep="\t",index_col=0)
y=d['Responder'].values.astype(int)
msi=d['MSI_H'].astype(str).isin(['True','TRUE']).values
ebv=d['EBV_pos'].astype(str).isin(['True','TRUE']).values
me=(msi|ebv).astype(float)
net=d['NetImmune'].values.astype(float)

def mk(): return make_pipeline(StandardScaler(),LogisticRegression(C=1e9,max_iter=5000))

MODELS={'MSI/EBV only':np.c_[me],'NetImmune only':np.c_[net],'NetImmune + MSI/EBV':np.c_[net,me]}

def loo_auc(X,y):
    p=np.zeros(len(y))
    for tr,te in LeaveOneOut().split(X):
        m=mk().fit(X[tr],y[tr]); p[te]=m.predict_proba(X[te])[:,1]
    return roc_auc_score(y,p)

def rep_cv(X,y,rep=200,k=5):
    a=[]
    for i in range(rep):
        skf=StratifiedKFold(k,shuffle=True,random_state=i)
        p=np.zeros(len(y))
        for tr,te in skf.split(X,y):
            m=mk().fit(X[tr],y[tr]); p[te]=m.predict_proba(X[te])[:,1]
        a.append(roc_auc_score(y,p))
    return np.array(a)

rows=[]
for name,X in MODELS.items():
    app=roc_auc_score(y,mk().fit(X,y).predict_proba(X)[:,1])
    lo=loo_auc(X,y); cv=rep_cv(X,y)
    rows.append(dict(Model=name,Apparent=app,LOOCV=lo,CV5_mean=cv.mean(),
                     CV5_lo=np.percentile(cv,2.5),CV5_hi=np.percentile(cv,97.5),
                     Optimism=app-cv.mean()))
r=pd.DataFrame(rows)
pd.set_option('display.width',200)
print("=== Combination model: apparent vs cross-validated AUC (n=44, 12 responders) ===")
print(r.round(3).to_string(index=False))

# likelihood-ratio test: does NetImmune add over MSI/EBV?
X0=sm.add_constant(np.c_[me]); X1=sm.add_constant(np.c_[me,(net-net.mean())/net.std()])
m0=sm.Logit(y,X0).fit(disp=0); m1=sm.Logit(y,X1).fit(disp=0)
from scipy.stats import chi2
lr=2*(m1.llf-m0.llf); p=chi2.sf(lr,1)
print(f"\n=== Likelihood-ratio test: MSI/EBV + NetImmune  vs  MSI/EBV alone ===")
print(f"LR chi2 = {lr:.3f}, df=1, p = {p:.4f}")
print(f"NetImmune coef (per SD) = {m1.params[2]:.3f}, OR = {np.exp(m1.params[2]):.2f}, p = {m1.pvalues[2]:.4f}")

r.to_csv(f"{O}/A3_cv_combo.csv",index=False)
with open(f"{O}/A3_lrt.txt","w") as f:
    f.write(f"LR chi2={lr:.4f} df=1 p={p:.5f}\nNetImmune per-SD coef={m1.params[2]:.4f} OR={np.exp(m1.params[2]):.3f} p={m1.pvalues[2]:.5f}\n")
print("\nsaved.")
