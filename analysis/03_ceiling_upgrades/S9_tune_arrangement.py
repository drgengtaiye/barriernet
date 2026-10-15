import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.linear_model import Ridge, LinearRegression
from scipy.stats import pearsonr
O="/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
M=pd.read_csv(f"{O}/S7_neigh_meta.csv"); X=np.load(f"{O}/S7_minibulk.npy"); genes=np.load(f"{O}/S7_genes.npy",allow_pickle=True)
secs=M.section.unique(); COMP=['m_imm','m_epi','m_bar','sd_imm','sd_epi']
Xz=np.zeros_like(X)
for s in secs:
    m=(M.section==s).values; sub=X[m]; sd=sub.std(0); sd[sd<1e-8]=1
    Xz[m]=(sub-sub.mean(0))/sd
ok=np.isfinite(Xz).all(0); Xz=Xz[:,ok]; genes=genes[ok]

def resid(col):
    y=M[col].values; f=np.isfinite(y)
    r=LinearRegression().fit(M.loc[f,COMP],y[f]); out=np.full(len(y),np.nan); out[f]=y[f]-r.predict(M.loc[f,COMP]); return out
M['excl_res']=resid('excl'); M['bar_res']=resid('bar_imm_excl')

def loso(y,K,alpha):
    f=np.isfinite(y); pr=np.full(len(y),np.nan)
    for s in secs:
        te=((M.section==s).values)&f; tr=(~(M.section==s).values)&f
        Xtr,ytr=Xz[tr],y[tr]
        c=(Xtr-Xtr.mean(0)).T@(ytr-ytr.mean())/(len(ytr)*Xtr.std(0)*ytr.std()+1e-12)  # marginal r, train only
        sel=np.argsort(-np.abs(c))[:K]
        pr[te]=Ridge(alpha=alpha).fit(Xtr[:,sel],ytr).predict(Xz[te][:,sel])
    per=np.array([pearsonr(y[((M.section==s).values)&f],pr[((M.section==s).values)&f])[0] for s in secs])
    return pearsonr(y[f],pr[f])[0], np.median(per), per

print("=== tuning: arrangement residual (composition removed) ===")
best=None
for K in [200,1000,4000]:
    for a in [30,300,3000]:
        r,med,_=loso(M.excl_res.values,K,a)
        print(f"  K={K:5d} alpha={a:6.0f}   LOSO r={r:+.3f}  per-section median={med:+.3f}")
        if best is None or med>best[0]: best=(med,K,a)
print(f"\nbest: K={best[1]} alpha={best[2]} (median r={best[0]:+.3f})")
r,med,per=loso(M.bar_res.values,best[1],best[2])
print(f"barrier-immune residual with same setting: LOSO r={r:+.3f} median={med:+.3f}")

K,a=best[1],best[2]
y=M.excl_res.values; f=np.isfinite(y)
c=(Xz[f]-Xz[f].mean(0)).T@(y[f]-y[f].mean())/(f.sum()*Xz[f].std(0)*y[f].std()+1e-12)
sel=np.argsort(-np.abs(c))[:K]
W=pd.Series(Ridge(alpha=a).fit(Xz[f][:,sel],y[f]).coef_,index=genes[sel]).sort_values()
W.to_csv(f"{O}/S9_arrangement_weights.csv")
print("\ninfiltrated-arrangement genes:",", ".join(W.index[:15]))
print("excluded-arrangement genes  :",", ".join(W.index[-15:][::-1]))
