import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.linear_model import Ridge, LinearRegression
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
O="/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
M=pd.read_csv(f"{O}/S7_neigh_meta.csv"); X=np.load(f"{O}/S7_minibulk.npy"); genes=np.load(f"{O}/S7_genes.npy",allow_pickle=True)
secs=M.section.unique(); print("sections:",len(secs)," neighborhoods:",len(M))

COMP=['m_imm','m_epi','m_bar','sd_imm','sd_epi']
# how much of the arrangement statistic is pure composition?
lm=LinearRegression().fit(M[COMP],M.excl)
print(f"\ncomposition-only R2 for exclusion statistic: {r2_score(M.excl,lm.predict(M[COMP])):.3f}")
M['excl_resid']=M.excl-lm.predict(M[COMP])
print(f"residual SD = {M.excl_resid.std():.3f} (raw SD {M.excl.std():.3f})")

# z-score genes WITHIN section (mirrors per-cohort standardisation in bulk)
Xz=np.zeros_like(X)
for s in secs:
    m=(M.section==s).values; sub=X[m]
    sd=sub.std(0); sd[sd<1e-8]=1
    Xz[m]=(sub-sub.mean(0))/sd
keep=np.isfinite(Xz).all(0) & (X.std(0)>1e-6)
Xz=Xz[:,keep]; genes=genes[keep]; print("genes used:",Xz.shape[1])

def loso(y,label,alpha=3000.0):
    pr=np.zeros(len(y)); 
    for s in secs:
        te=(M.section==s).values; tr=~te
        r=Ridge(alpha=alpha).fit(Xz[tr],y[tr]); pr[te]=r.predict(Xz[te])
    per=[pearsonr(y[(M.section==s).values],pr[(M.section==s).values])[0] for s in secs]
    print(f"  {label:34s} LOSO r = {pearsonr(y,pr)[0]:+.3f}   per-section median r = {np.median(per):+.3f}"
          f"   [{min(per):+.2f},{max(per):+.2f}]")
    return pr,np.array(per)

print("\n=== can aggregate expression predict spatial arrangement? (leave-one-section-out) ===")
p_raw,_ = loso(M.excl.values,      "raw exclusion (composition ok)")
p_res,per= loso(M.excl_resid.values,"exclusion RESIDUAL (arrangement only)")
p_bar,_ = loso(M.bar_imm_excl.values,"barrier-immune exclusion")

# sanity: shuffle target within section -> should collapse to ~0
rng=np.random.default_rng(0); ysh=M.excl_resid.values.copy()
for s in secs:
    m=(M.section==s).values; ysh[m]=rng.permutation(ysh[m])
loso(ysh,"NEGATIVE CONTROL (shuffled)")

np.save(f"{O}/S8_pred_resid.npy",p_res)
M['pred_resid']=p_res; M.to_csv(f"{O}/S8_neigh_pred.csv",index=False)

# final model on all sections -> gene weights for bulk scoring
r=Ridge(alpha=3000.0).fit(Xz,M.excl_resid.values)
W=pd.Series(r.coef_,index=genes).sort_values()
W.to_csv(f"{O}/S8_arrangement_weights.csv")
print("\ntop genes predicting INFILTRATED arrangement:", ", ".join(W.index[:12]))
print("top genes predicting EXCLUDED arrangement   :", ", ".join(W.index[-12:][::-1]))
