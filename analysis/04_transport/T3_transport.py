"""
Barrier as a transport property — compositional confound removed.

Tangram fractions sum to 1, so tumour spots are mechanically barrier-poor
(r(Epi,Barrier) = -0.3..-0.7). Routing *through* the tumour compartment therefore
makes the real arrangement look artefactually permeable. Fix: pose the transport
problem entirely inside the STROMAL compartment.

  domain D  : non-tumour spots (the tissue a lymphocyte must actually traverse)
  sources S : immune-rich spots within D
  targets B : tumour-boundary spots (in D, adjacent to the tumour compartment)
  cost      : rank-normalised barrier within D, bounded [1, CMAX]
  metric    : mean over B of (least-cost distance / hop distance)  -- dimensionless

  positive control : synthetic barrier shell around the tumour boundary
  null 1 (registration) : rigid translation of the barrier field
  null 2 (stratified)   : permute barrier within Epi deciles -> preserves the
                          compositional Epi-barrier relation, destroys geometry
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


def nbrs(rr, cc):
    idx = {(r, c): i for i, (r, c) in enumerate(zip(rr, cc))}
    adj = [[] for _ in range(len(rr))]
    for i, (r, c) in enumerate(zip(rr, cc)):
        for dr, dc in HEX:
            j = idx.get((r + dr, c + dc))
            if j is not None: adj[i].append(j)
    return adj


def gmat(adj, keep, cost=None):
    """graph restricted to `keep`; cost=None -> unweighted"""
    m = np.where(keep)[0]; pos = {v: k for k, v in enumerate(m)}
    I, J, W = [], [], []
    for i in m:
        for j in adj[i]:
            if keep[j]:
                w = 1.0 if cost is None else 0.5 * (cost[i] + cost[j])
                I.append(pos[i]); J.append(pos[j]); W.append(w)
    return sp.csr_matrix((W, (I, J)), shape=(len(m), len(m))), m, pos


def metric(adj, keep, barrier, src, tgt, hop):
    b = barrier[keep]
    c_local = 1.0 + (CMAX - 1.0) * (rankdata(b) - 1) / (max(len(b) - 1, 1))
    cost = np.zeros(len(barrier)); cost[keep] = c_local
    G, m, pos = gmat(adj, keep, cost)
    s = [pos[i] for i in np.where(src)[0] if i in pos]
    if not s: return np.nan
    d = dijkstra(G, indices=s, min_only=True)
    t = np.array([pos[i] for i in np.where(tgt)[0] if i in pos])
    if len(t) < 10: return np.nan
    wd, hd = d[t], hop[t]
    ok = np.isfinite(wd) & np.isfinite(hd) & (hd > 0)
    return np.nan if ok.sum() < 10 else float(np.mean(wd[ok] / hd[ok]))


def shift(rr, cc, v, dr, dc):
    r0, c0 = rr.min(), cc.min(); R_, C_ = rr.max() - r0 + 1, cc.max() - c0 + 1
    s = {((r - r0 + dr) % R_, (c - c0 + dc) % C_): v[i] for i, (r, c) in enumerate(zip(rr, cc))}
    o = np.array([s.get((r - r0, c - c0), np.nan) for r, c in zip(rr, cc)])
    o[np.isnan(o)] = np.nanmean(v); return o


def strat_perm(barrier, epi, k=10):
    o = barrier.copy()
    q = np.quantile(epi, np.linspace(0, 1, k + 1)); q[0] -= 1e-9
    for a_, b_ in zip(q[:-1], q[1:]):
        m = (epi > a_) & (epi <= b_)
        if m.sum() > 1: o[m] = rng.permutation(barrier[m])
    return o


rows = []
for f in sorted(glob.glob(f"{R}/*.with_barriernet.h5ad")):
    a = ad.read_h5ad(f); sec = os.path.basename(f).split('.')[0]
    a = a[a.obs.in_tissue.astype(int) == 1].copy()
    rr = a.obs.array_row.astype(int).values; cc = a.obs.array_col.astype(int).values
    z = lambda v: (v - v.mean()) / (v.std() + 1e-9)
    barrier = z(a.obs[['tg_Fibro', 'tg_Macro']].sum(1).values)
    immune  = z(a.obs[['tg_Tcell', 'tg_NK', 'tg_Bcell', 'tg_Plasma']].sum(1).values)
    epi     = z(a.obs['tg_Epi'].values)
    adj = nbrs(rr, cc)

    tumour = epi >= np.quantile(epi, .75)
    D = ~tumour                                          # travel domain = stroma
    bnd = np.array([D[i] and any(tumour[j] for j in adj[i]) for i in range(len(rr))])
    src = D & (immune >= np.quantile(immune[D], .90))
    if src.sum() < 5 or bnd.sum() < 15:
        print(f"  {sec}: skipped (src={src.sum()} bnd={bnd.sum()})"); continue

    Gu, m, pos = gmat(adj, D)
    s = [pos[i] for i in np.where(src)[0] if i in pos]
    hop = dijkstra(Gu, indices=s, min_only=True)         # unweighted, within D

    obs = metric(adj, D, barrier, src, bnd, hop)
    shell = np.array([3.0 if (D[i] and any(tumour[j] for j in adj[i])) else -1.0
                      for i in range(len(rr))])
    pos_ctrl = metric(adj, D, shell, src, bnd, hop)
    n1 = np.array([metric(adj, D, shift(rr, cc, barrier,
                    int(rng.integers(0, rr.max()-rr.min()+1)),
                    2*int(rng.integers(0, (cc.max()-cc.min())//2+1))), src, bnd, hop)
                   for _ in range(NPERM)])
    n2 = np.array([metric(adj, D, strat_perm(barrier, epi), src, bnd, hop)
                   for _ in range(NPERM)])
    n1 = n1[np.isfinite(n1)]; n2 = n2[np.isfinite(n2)]
    z1 = (obs - n1.mean()) / (n1.std() + 1e-12); p1 = (np.sum(n1 >= obs) + 1) / (len(n1) + 1)
    z2 = (obs - n2.mean()) / (n2.std() + 1e-12); p2 = (np.sum(n2 >= obs) + 1) / (len(n2) + 1)
    rows.append(dict(section=sec, n_bnd=int(bnd.sum()), n_src=int(src.sum()),
                     barrier_mean=float(a.obs[['tg_Fibro','tg_Macro']].sum(1).mean()),
                     detour=obs, pos_ctrl=pos_ctrl,
                     null_shift=float(n1.mean()), z_shift=float(z1), p_shift=float(p1),
                     null_strat=float(n2.mean()), z_strat=float(z2), p_strat=float(p2)))
    print(f"  {sec}: obs={obs:.3f}  shiftnull={n1.mean():.3f} z={z1:+5.2f} p={p1:.3f} | "
          f"stratnull={n2.mean():.3f} z={z2:+5.2f} p={p2:.3f}  [pos={pos_ctrl:.3f}]")

D_ = pd.DataFrame(rows); D_.to_csv(f"{O}/T3_transport.csv", index=False)
print("\n=== Topological excess exclusion within the stromal compartment ===")
print(D_[['section','barrier_mean','detour','null_shift','z_shift','p_shift',
          'null_strat','z_strat','p_strat','pos_ctrl']].round(3).to_string(index=False))
print(f"\npositive control median {D_.pos_ctrl.median():.3f} vs observed {D_.detour.median():.3f} "
      f"-> metric has power" if D_.pos_ctrl.median() > D_.detour.median() else "")
for lab, zc, pc in [("registration null", 'z_shift', 'p_shift'), ("stratified null", 'z_strat', 'p_strat')]:
    print(f"{lab:20s}: median z={D_[zc].median():+.2f}, {(D_[pc]<.05).sum()}/{len(D_)} sections p<0.05, "
          f"Wilcoxon(z>0) p={wilcoxon(D_[zc], alternative='greater').pvalue:.4f}")
r, pv = pearsonr(D_.barrier_mean, D_.z_shift)
print(f"barrier ABUNDANCE vs TOPOLOGY z: r={r:+.2f} p={pv:.3f}")
