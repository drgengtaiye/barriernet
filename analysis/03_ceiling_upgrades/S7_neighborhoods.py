import warnings; warnings.filterwarnings('ignore')
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp, glob, os
from collections import deque
R="/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/barriernet_spatial_demo"
O="/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
KHOP=3; STRIDE=2

# genes usable in every bulk cohort (score must be bulk-computable)
common=None  # gene intersection applied later, at scoring time
for t in ["TCGA","GSE84437","GSE62254"]:
    g=pd.read_csv(f"{O}/genes_{t}.txt",header=None)[0].tolist() if os.path.exists(f"{O}/genes_{t}.txt") else None
    if g is None: continue
    common=set(g) if common is None else common&set(g)
print("bulk-common genes preloaded:", len(common) if common else "none yet")

def neighbors(rows,cols):
    idx={(r,c):i for i,(r,c) in enumerate(zip(rows,cols))}
    adj=[[] for _ in range(len(rows))]
    for i,(r,c) in enumerate(zip(rows,cols)):
        for dr,dc in [(0,-2),(0,2),(-1,-1),(-1,1),(1,-1),(1,1)]:
            j=idx.get((r+dr,c+dc))
            if j is not None: adj[i].append(j)
    return adj

def khop(adj,s,k):
    seen={s}; q=deque([(s,0)])
    while q:
        v,d=q.popleft()
        if d==k: continue
        for w in adj[v]:
            if w not in seen: seen.add(w); q.append((w,d+1))
    return np.fromiter(seen,int)

FILES=sorted(glob.glob(f"{R}/*.with_barriernet.h5ad"))
gene_ref=None
for f in FILES:                                   # pass 1: common gene set across sections
    v=ad.read_h5ad(f, backed='r').var_names
    gene_ref=set(v) if gene_ref is None else gene_ref & set(v)
gene_ref=np.array(sorted(gene_ref)); print("common genes across sections:",len(gene_ref))

rows_out=[]; mb_list=[]
for f in FILES:
    a=ad.read_h5ad(f); sec=os.path.basename(f).split('.')[0]
    a=a[a.obs.in_tissue.astype(int)==1].copy()
    rr=a.obs.array_row.astype(int).values; cc=a.obs.array_col.astype(int).values
    adj=neighbors(rr,cc)
    imm=a.obs[['tg_Tcell','tg_NK','tg_Bcell','tg_Plasma']].sum(1).values
    epi=a.obs['tg_Epi'].values
    bar=a.obs[['tg_Fibro','tg_Macro']].sum(1).values
    a=a[:, gene_ref].copy()                       # align gene order across sections
    X=a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    centers=np.arange(a.n_obs)[::STRIDE]
    for s in centers:
        nb=khop(adj,s,KHOP)
        if len(nb)<12: continue
        i_,e_,b_=imm[nb],epi[nb],bar[nb]
        if i_.std()<1e-8 or e_.std()<1e-8: continue
        # --- spatial arrangement statistic: immune-epithelial mutual exclusion within the patch
        excl=-np.corrcoef(i_,e_)[0,1]
        # barrier interposition: does barrier sit where immune is low but epi is high?
        bi=-np.corrcoef(b_,i_)[0,1] if b_.std()>1e-8 else np.nan
        rows_out.append(dict(section=sec,center=int(s),n_spots=len(nb),
                             excl=excl, bar_imm_excl=bi,
                             m_imm=i_.mean(), m_epi=e_.mean(), m_bar=b_.mean(),
                             sd_imm=i_.std(), sd_epi=e_.std()))
        mb_list.append(np.asarray(X[nb].mean(0)).ravel())      # mini-bulk: spatial info destroyed
    print(f"  {sec}: {a.n_obs} spots -> {len(centers)} centers",flush=True)

M=pd.DataFrame(rows_out)
MB=np.vstack(mb_list).astype(np.float32)
np.save(f"{O}/S7_minibulk.npy",MB); np.save(f"{O}/S7_genes.npy",gene_ref)
M.to_csv(f"{O}/S7_neigh_meta.csv",index=False)
print("\nneighborhoods:",MB.shape," meta:",M.shape)
print(M[['excl','bar_imm_excl','m_imm','m_epi','m_bar']].describe().round(3).to_string())
