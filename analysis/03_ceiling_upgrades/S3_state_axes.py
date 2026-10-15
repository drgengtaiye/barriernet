import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
O="/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
lg=pd.read_parquet(f"{O}/S1_pseudobulk_logcpm.parquet")
meta=pd.read_csv(f"{O}/S1_pseudobulk_meta.csv").set_index('group').loc[lg.index]
lg=lg[meta.n_cells.values>=30]; meta=meta[meta.n_cells>=30]
CT=['Fibro','Macro','Tcell','NK','Bcell','Plasma','DC','Endo','Peri','Epi']
NPC=2
axes={}
print("=== state axes per cell type (PCA on across-patient variation, within cell type) ===")
for c in CT:
    try: sg=pd.read_csv(f"{O}/S2_state_{c}.csv").gene.tolist()
    except Exception: continue
    m=meta.celltype.values==c
    X=lg.loc[m,sg]
    if X.shape[0]<10 or X.shape[1]<15: continue
    Xc=X-X.mean(0)
    p=PCA(n_components=min(NPC,X.shape[0]-1,X.shape[1])).fit(Xc.values)
    # stability: leave-one-patient-out cosine similarity of PC1
    sims=[]
    for i in range(X.shape[0]):
        Xi=Xc.drop(Xc.index[i]); pi=PCA(n_components=1).fit(Xi.values-Xi.values.mean(0))
        sims.append(abs(np.dot(pi.components_[0],p.components_[0])))
    for k in range(p.n_components_):
        L=pd.Series(p.components_[k],index=sg)
        axes[f"{c}_S{k+1}"]=L
        top=L.sort_values(ascending=False)
        print(f"\n{c}_S{k+1}  var={p.explained_variance_ratio_[k]*100:.1f}%"
              + (f"  LOO-stability(PC1) r={np.mean(sims):.2f}" if k==0 else ""))
        print(f"   + {', '.join(top.index[:8])}")
        print(f"   - {', '.join(top.index[-8:][::-1])}")
A=pd.DataFrame(axes)
A.to_csv(f"{O}/S3_state_loadings.csv")
ident={c:pd.read_csv(f"{O}/S2_ident_{c}.csv").gene.tolist() for c in CT
       if __import__('os').path.exists(f"{O}/S2_ident_{c}.csv")}
pd.Series({k:",".join(v) for k,v in ident.items()}).to_csv(f"{O}/S3_identity_genes.csv")
print(f"\nsaved {A.shape[1]} state axes over {len(ident)} cell types")
