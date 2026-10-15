"""
Transport analysis re-run on MARKER-BASED scores instead of Tangram fractions.

Tangram fractions are compositional (sum to 1), so barrier and epithelium are
mechanically anti-correlated. Independent marker z-scores carry no such constraint,
so if the result survives here the compositional confound is excluded by
construction rather than by control.
"""
import warnings, functools; warnings.filterwarnings('ignore')
print = functools.partial(print, flush=True)
import anndata as ad, numpy as np, pandas as pd, glob, os
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra
from scipy.stats import rankdata, wilcoxon, pearsonr

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/barriernet_spatial_demo"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
NPERM, CMAX = 200, 6.0
rng = np.random.default_rng(0)
HEX = [(0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1)]

MK = dict(
    barrier=['COL1A1','COL1A2','COL3A1','LUM','DCN','POSTN','FN1','ACTA2','TAGLN','THY1','FAP',
             'SPP1','APOE','TREM2','C1QA','C1QB','C1QC','MARCO','CD68','CSF1R'],
    immune =['CD3D','CD3E','CD2','TRAC','CD8A','GZMB','NKG7','MS4A1','CD79A','CXCL13','LTB','MZB1','PTPRC'],
    epi    =['EPCAM','KRT8','KRT18','KRT19','CDH1','CLDN4','CLDN7','TFF1','TFF3','MUC1'])


def nbrs(rr, cc):
    idx = {(r, c): i for i, (r, c) in enumerate(zip(rr, cc))}
    return [[idx[(r+dr, c+dc)] for dr, dc in HEX if (r+dr, c+dc) in idx] for r, c in zip(rr, cc)]

def gmat(adj, keep, cost=None):
    m = np.where(keep)[0]; pos = {v: k for k, v in enumerate(m)}
    I, J, W = [], [], []
    for i in m:
        for j in adj[i]:
            if keep[j]:
                I.append(pos[i]); J.append(pos[j]); W.append(1.0 if cost is None else 0.5*(cost[i]+cost[j]))
    return sp.csr_matrix((W, (I, J)), shape=(len(m), len(m))), pos

def metric(adj, keep, barrier, src, tgt, hop):
    b = barrier[keep]
    c = np.zeros(len(barrier)); c[keep] = 1.0 + (CMAX-1.0)*(rankdata(b)-1)/max(len(b)-1, 1)
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


rows = []
for f in sorted(glob.glob(f"{R}/*.with_barriernet.h5ad")):
    a = ad.read_h5ad(f); sec = os.path.basename(f).split('.')[0]
    a = a[a.obs.in_tissue.astype(int) == 1].copy()
    rr = a.obs.array_row.astype(int).values; cc = a.obs.array_col.astype(int).values
    X = a.X.tocsc() if sp.issparse(a.X) else sp.csc_matrix(a.X)
    gi = {g: i for i, g in enumerate(a.var_names)}

    def score(genes):
        idx = [gi[g] for g in genes if g in gi]
        if len(idx) < 4: return None
        M = np.asarray(X[:, idx].todense())
        M = (M - M.mean(0)) / (M.std(0) + 1e-9)
        v = M.mean(1); return (v - v.mean()) / (v.std() + 1e-9)

    barrier, immune, epi = score(MK['barrier']), score(MK['immune']), score(MK['epi'])
    if barrier is None or immune is None or epi is None:
        print(f"  {sec}: skipped (markers missing)"); continue

    # the confound this analysis is designed to exclude
    r_be = float(np.corrcoef(epi, barrier)[0, 1])

    adj = nbrs(rr, cc)
    tumour = epi >= np.quantile(epi, .75); D = ~tumour
    bnd = np.array([D[i] and any(tumour[j] for j in adj[i]) for i in range(len(rr))])
    src = D & (immune >= np.quantile(immune[D], .90))
    if src.sum() < 5 or bnd.sum() < 15:
        print(f"  {sec}: skipped"); continue
    Gu, pos = gmat(adj, D)
    s = [pos[i] for i in np.where(src)[0] if i in pos]
    hop = dijkstra(Gu, indices=s, min_only=True)

    obs = metric(adj, D, barrier, src, bnd, hop)
    shell = np.array([3.0 if bnd[i] else -1.0 for i in range(len(rr))])
    pc = metric(adj, D, shell, src, bnd, hop)
    nl = np.array([metric(adj, D, shift(rr, cc, barrier,
                    int(rng.integers(0, rr.max()-rr.min()+1)),
                    2*int(rng.integers(0, (cc.max()-cc.min())//2+1))), src, bnd, hop)
                   for _ in range(NPERM)])
    nl = nl[np.isfinite(nl)]
    z = (obs - nl.mean())/(nl.std()+1e-12); p = (np.sum(nl >= obs)+1)/(len(nl)+1)
    rows.append(dict(section=sec, r_epi_barrier=r_be, detour=obs, pos_ctrl=pc,
                     null=float(nl.mean()), z=float(z), p=float(p)))
    print(f"  {sec}: r(Epi,Barrier)={r_be:+.2f}  obs={obs:.3f} null={nl.mean():.3f} "
          f"z={z:+5.2f} p={p:.3f}  [pos={pc:.3f}]")

D_ = pd.DataFrame(rows); D_.to_csv(f"{O}/T6_marker_based.csv", index=False)
print("\n=== Marker-based replication (no compositional constraint) ===")
print(D_.round(3).to_string(index=False))
print(f"\nr(Epi,Barrier) range: {D_.r_epi_barrier.min():+.2f} to {D_.r_epi_barrier.max():+.2f}"
      "   (Tangram version was -0.32 to -0.66 by construction)")
print(f"median z = {D_.z.median():+.2f}; {(D_.p<.05).sum()}/{len(D_)} sections p<0.05")
print(f"Wilcoxon(z>0) p={wilcoxon(D_.z, alternative='greater').pvalue:.4f}  "
      f"| Wilcoxon(z<0) p={wilcoxon(D_.z, alternative='less').pvalue:.4f}")
print(f"positive control median {D_.pos_ctrl.median():.3f} vs observed {D_.detour.median():.3f}")
