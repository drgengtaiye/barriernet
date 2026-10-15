import warnings; warnings.filterwarnings('ignore')
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp
R="/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
O="/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
f=f"{R}/out_scanpy_strict_GSE183904/GSE183904_TumorOnly_scanpy_strict_final_scored.h5ad"
a=ad.read_h5ad(f, backed='r')
obs=a.obs[['sample_id','anno_level1','tissue']].copy()
genes=np.array(a.var_names); n=a.shape[0]
key=(obs['sample_id'].astype(str)+"|"+obs['anno_level1'].astype(str)).values
groups=pd.unique(key); gidx={g:i for i,g in enumerate(groups)}
S=np.zeros((len(groups),len(genes)),dtype=np.float64); C=np.zeros(len(groups),dtype=np.int64)
CH=20000
for s in range(0,n,CH):
    e=min(s+CH,n)
    X=a[s:e].to_memory().layers['counts']
    X=X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    k=key[s:e]
    for g in pd.unique(k):
        m=(k==g); S[gidx[g]]+=np.asarray(X[m].sum(0)).ravel(); C[gidx[g]]+=m.sum()
    print(f"  {e}/{n}",flush=True)
meta=pd.DataFrame({'group':groups,'sample':[g.split('|')[0] for g in groups],
                   'celltype':[g.split('|')[1] for g in groups],'n_cells':C})
cpm=S/np.maximum(S.sum(1,keepdims=True),1)*1e6
lg=pd.DataFrame(np.log1p(cpm),index=groups,columns=genes)
meta.to_csv(f"{O}/S1_pseudobulk_meta.csv",index=False)
lg.to_parquet(f"{O}/S1_pseudobulk_logcpm.parquet")
print("\npseudobulk:",lg.shape)
print(meta.groupby('celltype')['n_cells'].agg(['count','sum']).to_string())
