"""Sensitivity of the transport result to every arbitrary threshold."""
import warnings, functools; warnings.filterwarnings('ignore')
print = functools.partial(print, flush=True)
import anndata as ad, numpy as np, pandas as pd, glob, os, itertools
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra
from scipy.stats import rankdata, wilcoxon

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/barriernet_spatial_demo"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
NPERM = 50
rng = np.random.default_rng(1)
HEX = [(0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1)]

def nbrs(rr, cc):
    idx = {(r, c): i for i, (r, c) in enumerate(zip(rr, cc))}
    return [[idx[(r+dr, c+dc)] for dr, dc in HEX if (r+dr, c+dc) in idx]
            for r, c in zip(rr, cc)]

def gmat(adj, keep, cost=None):
    m = np.where(keep)[0]; pos = {v: k for k, v in enumerate(m)}
    I, J, W = [], [], []
    for i in m:
        for j in adj[i]:
            if keep[j]:
                I.append(pos[i]); J.append(pos[j])
                W.append(1.0 if cost is None else 0.5*(cost[i]+cost[j]))
    return sp.csr_matrix((W, (I, J)), shape=(len(m), len(m))), pos

def metric(adj, keep, barrier, src, tgt, hop, cmax):
    b = barrier[keep]
    c = np.zeros(len(barrier))
    c[keep] = 1.0 + (cmax-1.0)*(rankdata(b)-1)/max(len(b)-1, 1)
    G, pos = gmat(adj, keep, c)
    s = [pos[i] for i in np.where(src)[0] if i in pos]
    if not s: return np.nan
    d = dijkstra(G, indices=s, min_only=True)
    t = np.array([pos[i] for i in np.where(tgt)[0] if i in pos])
    if len(t) < 10: return np.nan
    ok = np.isfinite(d[t]) & np.isfinite(hop[t]) & (hop[t] > 0)
    return np.nan if ok.sum() < 10 else float(np.mean(d[t][ok]/hop[t][ok]))

def shift(rr, cc, v, dr, dc):
    r0, c0 = rr.min(), cc.min(); R_, C_ = rr.max()-r0+1, cc.max()-c0+1
    s = {((r-r0+dr) % R_, (c-c0+dc) % C_): v[i] for i, (r, c) in enumerate(zip(rr, cc))}
    o = np.array([s.get((r-r0, c-c0), np.nan) for r, c in zip(rr, cc)])
    o[np.isnan(o)] = np.nanmean(v); return o

SEC = []
for f in sorted(glob.glob(f"{R}/*.with_barriernet.h5ad")):
    a = ad.read_h5ad(f); a = a[a.obs.in_tissue.astype(int) == 1].copy()
    z = lambda v: (v-v.mean())/(v.std()+1e-9)
    SEC.append(dict(name=os.path.basename(f).split('.')[0],
                    rr=a.obs.array_row.astype(int).values, cc=a.obs.array_col.astype(int).values,
                    bar=z(a.obs[['tg_Fibro','tg_Macro']].sum(1).values),
                    imm=z(a.obs[['tg_Tcell','tg_NK','tg_Bcell','tg_Plasma']].sum(1).values),
                    epi=z(a.obs['tg_Epi'].values)))
print(f"loaded {len(SEC)} sections\n")

out = []
for qt, qi, cmax in itertools.product([0.60, 0.75, 0.90], [0.85, 0.90, 0.95], [3.0, 6.0, 12.0]):
    zs = []
    for S in SEC:
        rr, cc = S['rr'], S['cc']; adj = nbrs(rr, cc)
        tum = S['epi'] >= np.quantile(S['epi'], qt); D = ~tum
        bnd = np.array([D[i] and any(tum[j] for j in adj[i]) for i in range(len(rr))])
        if D.sum() < 50: continue
        src = D & (S['imm'] >= np.quantile(S['imm'][D], qi))
        if src.sum() < 5 or bnd.sum() < 15: continue
        Gu, pos = gmat(adj, D)
        s = [pos[i] for i in np.where(src)[0] if i in pos]
        if not s: continue
        hop = dijkstra(Gu, indices=s, min_only=True)
        obs = metric(adj, D, S['bar'], src, bnd, hop, cmax)
        nl = np.array([metric(adj, D, shift(rr, cc, S['bar'],
                        int(rng.integers(0, rr.max()-rr.min()+1)),
                        2*int(rng.integers(0, (cc.max()-cc.min())//2+1))), src, bnd, hop, cmax)
                       for _ in range(NPERM)])
        nl = nl[np.isfinite(nl)]
        if not np.isfinite(obs) or len(nl) < 20: continue
        zs.append((obs-nl.mean())/(nl.std()+1e-12))
    if len(zs) >= 8:
        zs = np.array(zs)
        out.append(dict(tumour_q=qt, immune_q=qi, CMAX=cmax, n_sec=len(zs),
                        median_z=float(np.median(zs)), frac_negative=float((zs < 0).mean()),
                        wilcox_p_less=float(wilcoxon(zs, alternative='less').pvalue)))
        print(f"  tumour_q={qt} immune_q={qi} CMAX={cmax:<5}  median z={np.median(zs):+6.2f}  "
              f"{(zs<0).mean()*100:3.0f}% of sections negative  p={out[-1]['wilcox_p_less']:.4f}")

D_ = pd.DataFrame(out); D_.to_csv(f"{O}/T4_sensitivity.csv", index=False)
print(f"\n=== {len(D_)} parameter settings ===")
print(f"median z across settings : {D_.median_z.median():+.2f}  "
      f"(range {D_.median_z.min():+.2f} to {D_.median_z.max():+.2f})")
print(f"settings where the majority of sections are negative: "
      f"{(D_.frac_negative > .5).sum()}/{len(D_)}")
print(f"settings with Wilcoxon p<0.05 (more permeable than null): "
      f"{(D_.wilcox_p_less < .05).sum()}/{len(D_)}")
