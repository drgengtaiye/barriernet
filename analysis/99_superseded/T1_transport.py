"""
Barrier as a TRANSPORT property, not a concentration.

Edge cost on the Visium lattice rises with local barrier content; immune
accessibility of a tumour spot = least-cost path from the immune compartment.
The key control is a *registration* null: rigid translations of the barrier
field with wraparound preserve its abundance AND its own spatial autocorrelation,
and randomise only its positioning relative to tumour/immune. A significant
z-score therefore means the barrier is geometrically arranged to exclude, over
and above how much of it there is.
"""
import warnings, functools; warnings.filterwarnings('ignore')
print = functools.partial(print, flush=True)
import anndata as ad, numpy as np, pandas as pd, glob, os
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/barriernet_spatial_demo"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
NPERM, BETA = 200, 4.0
rng = np.random.default_rng(0)

HEX = [(0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1)]

def build_graph(rr, cc, cost_node):
    """Undirected lattice graph; edge cost = mean barrier cost of its two endpoints."""
    idx = {(r, c): i for i, (r, c) in enumerate(zip(rr, cc))}
    I, J, W = [], [], []
    for i, (r, c) in enumerate(zip(rr, cc)):
        for dr, dc in HEX:
            j = idx.get((r + dr, c + dc))
            if j is not None and j > i:
                w = 0.5 * (cost_node[i] + cost_node[j])
                I += [i, j]; J += [j, i]; W += [w, w]
    return sp.csr_matrix((W, (I, J)), shape=(len(rr), len(rr)))

def accessibility(rr, cc, barrier, immune_src, tumour):
    """Mean least-cost distance from the immune compartment to tumour spots."""
    cost = np.exp(BETA * barrier)                      # local impedance
    G = build_graph(rr, cc, cost)
    d = dijkstra(G, indices=np.where(immune_src)[0], min_only=True)
    d = d[tumour]
    d = d[np.isfinite(d)]
    return np.nan if len(d) == 0 else float(np.mean(d))

def shift_field(rr, cc, vals, dr, dc):
    """Rigid translation of a field on the lattice with wraparound (registration null)."""
    r0, r1 = rr.min(), rr.max(); c0, c1 = cc.min(), cc.max()
    R_, C_ = r1 - r0 + 1, c1 - c0 + 1
    src = {}
    for i, (r, c) in enumerate(zip(rr, cc)):
        src[((r - r0 + dr) % R_, (c - c0 + dc) % C_)] = vals[i]
    out = np.empty(len(rr))
    for i, (r, c) in enumerate(zip(rr, cc)):
        out[i] = src.get((r - r0, c - c0), np.nan)
    m = np.isnan(out); out[m] = np.nanmean(vals)       # spots with no source keep field mean
    return out

rows = []
for f in sorted(glob.glob(f"{R}/*.with_barriernet.h5ad")):
    a = ad.read_h5ad(f); sec = os.path.basename(f).split('.')[0]
    a = a[a.obs.in_tissue.astype(int) == 1].copy()
    rr = a.obs.array_row.astype(int).values; cc = a.obs.array_col.astype(int).values
    z = lambda v: (v - v.mean()) / (v.std() + 1e-9)
    barrier = z(a.obs[['tg_Fibro', 'tg_Macro']].sum(1).values)
    immune  = z(a.obs[['tg_Tcell', 'tg_NK', 'tg_Bcell', 'tg_Plasma']].sum(1).values)
    epi     = z(a.obs['tg_Epi'].values)
    src = immune >= np.quantile(immune, .90)
    tum = epi >= np.quantile(epi, .75)
    if src.sum() < 5 or tum.sum() < 20:
        print(f"  {sec}: skipped"); continue

    obs = accessibility(rr, cc, barrier, src, tum)
    null = []
    for _ in range(NPERM):
        dr = int(rng.integers(0, rr.max() - rr.min() + 1))
        dc = 2 * int(rng.integers(0, (cc.max() - cc.min()) // 2 + 1))
        b = shift_field(rr, cc, barrier, dr, dc)
        null.append(accessibility(rr, cc, b, src, tum))
    null = np.array([x for x in null if np.isfinite(x)])
    zsc = (obs - null.mean()) / (null.std() + 1e-9)
    p = (np.sum(null >= obs) + 1) / (len(null) + 1)
    rows.append(dict(section=sec, n_spots=a.n_obs, n_tumour=int(tum.sum()), n_src=int(src.sum()),
                     barrier_mean=float(a.obs[['tg_Fibro', 'tg_Macro']].sum(1).mean()),
                     access_obs=obs, null_mean=float(null.mean()), null_sd=float(null.std()),
                     z_topology=float(zsc), p_perm=float(p)))
    print(f"  {sec}: obs={obs:8.2f} null={null.mean():8.2f}±{null.std():6.2f}  z={zsc:+5.2f}  p={p:.3f}")

D = pd.DataFrame(rows)
D.to_csv(f"{O}/T1_transport.csv", index=False)
print("\n=== Topological excess exclusion (positive z = barrier positioned to exclude) ===")
print(D[['section', 'barrier_mean', 'access_obs', 'null_mean', 'z_topology', 'p_perm']].round(3).to_string(index=False))
from scipy.stats import wilcoxon, pearsonr
print(f"\nz across sections: median {D.z_topology.median():+.2f}, "
      f"{(D.p_perm < .05).sum()}/{len(D)} sections p<0.05")
print("one-sample Wilcoxon (z>0):", f"p={wilcoxon(D.z_topology, alternative='greater').pvalue:.4f}")
r, pv = pearsonr(D.barrier_mean, D.z_topology)
print(f"correlation between barrier ABUNDANCE and TOPOLOGY z: r={r:+.2f}, p={pv:.3f}"
      "   <- low |r| means topology is a genuinely new axis")
