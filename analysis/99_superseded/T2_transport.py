"""
Barrier as a transport property — numerically stable formulation.

metric  = detour ratio = mean over tumour spots of (least-cost distance / hop distance)
          with edge cost bounded in [1, CMAX] via rank normalisation of the barrier
          field. Dimensionless: geometry and overall barrier level are divided out,
          so what remains is whether the barrier lies ON the paths that matter.

controls
  positive : synthetic barrier shell placed around the tumour compartment
             -> the metric must rise sharply, otherwise it has no power
  null     : rigid translation of the real barrier field (preserves its abundance
             and its own autocorrelation, randomises only its registration)
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


def edges(rr, cc):
    idx = {(r, c): i for i, (r, c) in enumerate(zip(rr, cc))}
    E = []
    for i, (r, c) in enumerate(zip(rr, cc)):
        for dr, dc in HEX:
            j = idx.get((r + dr, c + dc))
            if j is not None and j > i:
                E.append((i, j))
    return np.array(E), idx


def graph(E, n, node_cost=None):
    i, j = E[:, 0], E[:, 1]
    w = np.ones(len(E)) if node_cost is None else 0.5 * (node_cost[i] + node_cost[j])
    return sp.csr_matrix((np.r_[w, w], (np.r_[i, j], np.r_[j, i])), shape=(n, n))


def detour(E, n, barrier, src, tum, hop_d):
    """mean (weighted distance / hop distance) over reachable tumour spots"""
    c = 1.0 + (CMAX - 1.0) * (rankdata(barrier) - 1) / (len(barrier) - 1)   # bounded [1, CMAX]
    d = dijkstra(graph(E, n, c), indices=np.where(src)[0], min_only=True)
    ok = tum & np.isfinite(d) & np.isfinite(hop_d) & (hop_d > 0)
    return np.nan if ok.sum() < 10 else float(np.mean(d[ok] / hop_d[ok]))


def shift(rr, cc, v, dr, dc):
    r0, c0 = rr.min(), cc.min(); R_, C_ = rr.max() - r0 + 1, cc.max() - c0 + 1
    src = {((r - r0 + dr) % R_, (c - c0 + dc) % C_): v[i] for i, (r, c) in enumerate(zip(rr, cc))}
    out = np.array([src.get((r - r0, c - c0), np.nan) for r, c in zip(rr, cc)])
    out[np.isnan(out)] = np.nanmean(v)
    return out


rows = []
for f in sorted(glob.glob(f"{R}/*.with_barriernet.h5ad")):
    a = ad.read_h5ad(f); sec = os.path.basename(f).split('.')[0]
    a = a[a.obs.in_tissue.astype(int) == 1].copy()
    rr = a.obs.array_row.astype(int).values; cc = a.obs.array_col.astype(int).values
    n = a.n_obs
    z = lambda v: (v - v.mean()) / (v.std() + 1e-9)
    barrier = z(a.obs[['tg_Fibro', 'tg_Macro']].sum(1).values)
    immune  = z(a.obs[['tg_Tcell', 'tg_NK', 'tg_Bcell', 'tg_Plasma']].sum(1).values)
    epi     = z(a.obs['tg_Epi'].values)
    src = immune >= np.quantile(immune, .90)
    tum = epi   >= np.quantile(epi, .75)
    if src.sum() < 5 or tum.sum() < 20:
        print(f"  {sec}: skipped"); continue
    E, _ = edges(rr, cc)
    hop = dijkstra(graph(E, n), indices=np.where(src)[0], min_only=True)   # unweighted

    obs = detour(E, n, barrier, src, tum, hop)

    # positive control: synthetic shell of barrier at hop distance 2-3 from the immune source
    synth = np.where((hop >= 2) & (hop <= 3), 3.0, -1.0).astype(float)
    pos = detour(E, n, synth, src, tum, hop)

    null = np.array([detour(E, n,
                            shift(rr, cc, barrier,
                                  int(rng.integers(0, rr.max() - rr.min() + 1)),
                                  2 * int(rng.integers(0, (cc.max() - cc.min()) // 2 + 1))),
                            src, tum, hop) for _ in range(NPERM)])
    null = null[np.isfinite(null)]
    zt = (obs - null.mean()) / (null.std() + 1e-12)
    p = (np.sum(null >= obs) + 1) / (len(null) + 1)
    rows.append(dict(section=sec, n_spots=n, n_tumour=int(tum.sum()),
                     barrier_mean=float(a.obs[['tg_Fibro', 'tg_Macro']].sum(1).mean()),
                     detour_obs=obs, detour_poscontrol=pos,
                     null_mean=float(null.mean()), null_sd=float(null.std()),
                     z_topology=float(zt), p_perm=float(p)))
    print(f"  {sec}: obs={obs:.3f}  null={null.mean():.3f}±{null.std():.3f}  "
          f"z={zt:+5.2f} p={p:.3f}   [pos-control={pos:.3f}]")

D = pd.DataFrame(rows); D.to_csv(f"{O}/T2_transport.csv", index=False)
print("\n=== Topological excess exclusion (z>0 = barrier positioned to obstruct) ===")
print(D[['section', 'barrier_mean', 'detour_obs', 'null_mean', 'z_topology', 'p_perm',
         'detour_poscontrol']].round(3).to_string(index=False))
print(f"\npositive control: median detour {D.detour_poscontrol.median():.3f} vs real "
      f"{D.detour_obs.median():.3f} vs null {D.null_mean.median():.3f}")
print(f"z across sections: median {D.z_topology.median():+.2f}; "
      f"{(D.p_perm < .05).sum()}/{len(D)} sections p<0.05")
print("one-sample Wilcoxon (z>0):", f"p={wilcoxon(D.z_topology, alternative='greater').pvalue:.4f}")
r, pv = pearsonr(D.barrier_mean, D.z_topology)
print(f"barrier ABUNDANCE vs TOPOLOGY z: r={r:+.2f} p={pv:.3f}  (low |r| = topology is a new axis)")
