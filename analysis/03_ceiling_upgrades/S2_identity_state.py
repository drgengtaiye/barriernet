import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
O="/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
lg=pd.read_parquet(f"{O}/S1_pseudobulk_logcpm.parquet")
meta=pd.read_csv(f"{O}/S1_pseudobulk_meta.csv").set_index('group').loc[lg.index]
lg=lg[meta.n_cells.values>=30]; meta=meta[meta.n_cells>=30]     # need enough cells per profile
CT=sorted(meta.celltype.unique()); print("cell types:",CT)
print("profiles retained:",lg.shape[0])

ct_mean=lg.groupby(meta.celltype.values).mean()                  # celltype x gene
expressed=(ct_mean.max(0)>1.0)                                   # drop near-silent genes
lg=lg.loc[:,expressed]; ct_mean=ct_mean.loc[:,expressed]
print("genes after expression filter:",lg.shape[1])

rows=[]
for c in CT:
    m=meta.celltype.values==c
    if m.sum()<8: continue
    Xc=lg[m]                                                     # patients x genes for this cell type
    mu_c=ct_mean.loc[c]
    oth=ct_mean.drop(index=c)
    # robust specificity: tolerate 1-2 lineage-related cell types being similar
    spec=mu_c-oth.quantile(.78)
    spec_strict=mu_c-oth.max(0)
    sd_within=Xc.std(0)                                          # across-patient variability within c
    rows.append(pd.DataFrame({'celltype':c,'gene':lg.columns,'mean':mu_c.values,
                              'spec':spec.values,'spec_strict':spec_strict.values,'sd_within':sd_within.values}))
G=pd.concat(rows,ignore_index=True)
G.to_csv(f"{O}/S2_gene_stats.csv.gz",index=False,compression='gzip')

SPEC_MIN=0.8
print(f"\n=== genes specific to each cell type (spec>{SPEC_MIN} logCPM vs next best) ===")
summ=[]
for c in CT:
    g=G[(G.celltype==c)&(G.spec>SPEC_MIN)&(G['mean']>1.2)&(G.spec_strict>-1.0)]
    if len(g)<10: summ.append((c,len(g),0,0)); continue
    q=g.sd_within.quantile([.33,.67])
    ident=g[g.sd_within<=q.iloc[0]].nlargest(40,'spec')           # specific AND stable -> abundance
    state=g[g.sd_within>=q.iloc[1]].nlargest(80,'sd_within')      # specific AND variable -> state
    summ.append((c,len(g),len(ident),len(state)))
    ident.assign(role='identity').to_csv(f"{O}/S2_ident_{c}.csv",index=False)
    state.assign(role='state').to_csv(f"{O}/S2_state_{c}.csv",index=False)
s=pd.DataFrame(summ,columns=['celltype','n_specific','n_identity','n_state'])
print(s.to_string(index=False))
print("\nexample identity markers:")
for c in ['Fibro','Macro','Tcell','Epi']:
    try: print(f"  {c:7s}",list(pd.read_csv(f'{O}/S2_ident_{c}.csv').gene[:12]))
    except Exception: pass
print("\nexample state genes (high across-patient variance within cell type):")
for c in ['Fibro','Macro']:
    try: print(f"  {c:7s}",list(pd.read_csv(f'{O}/S2_state_{c}.csv').gene[:12]))
    except Exception: pass
