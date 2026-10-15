"""Recompute TLS-centred ring zones with the correct Visium neighbourhood.

barriernet 0.2.0 (spatial.py) uses HEX_STEPS = (1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1): axial
offsets. Visium array coordinates are doubled-column (a spot's six neighbours are (0,+-2) and
(+-1,+-1)), so only the (1,-1)/(-1,1) steps ever match and the "hex graph" degenerates into
chains along one diagonal. This script (1) reproduces the stored ring labels with the old
offsets, as a check that the pipeline is understood, and (2) recomputes rings on the true
hexagonal graph, from the same TLS peaks, and the per-section Ring - Far contrasts and
BarrierIndex (z of CAF_ECM Ring1-Far + z of SPP1_TAM Ring1-Far, ddof 0, as stored).
Far: >= 3 hops (as in the stored files) and, as a sensitivity, >= 4 hops excluding tissue
edge spots (fewer than six in-tissue neighbours).
"""
import numpy as np, pandas as pd, anndata as ad, glob, os
from collections import deque
from scipy.stats import wilcoxon

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2"
OLD = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
HEX = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
PROGS = ["state_CAF_ECM", "state_SPP1_TAM", "state_IL10_STAT3", "state_CAF_TGFB", "state_STF3_FAP_INHBA"]


def dist_from_peaks(o, steps):
    rc = list(zip(o.array_row.astype(int), o.array_col.astype(int)))
    pos = {p: i for i, p in enumerate(rc)}
    nb = [[pos[(r + a, c + b)] for a, b in steps if (r + a, c + b) in pos] for r, c in rc]
    d = np.full(len(rc), np.inf); q = deque()
    for i in np.where(o.TLS_peaks.astype(int).values == 1)[0]:
        d[i] = 0; q.append(i)
    while q:
        u = q.popleft()
        for v in nb[u]:
            if d[v] > d[u] + 1:
                d[v] = d[u] + 1; q.append(v)
    deg = np.array([len(x) for x in nb])
    return d, deg


def zones(d, far_min=3, edge=None):
    z = np.full(len(d), "Other", dtype=object)
    z[d == 0] = "Hotspot"; z[d == 1] = "Ring1"; z[d == 2] = "Ring2"
    far = d >= far_min
    if edge is not None:
        far &= ~edge
    z[far] = "Far"
    return z


rows, check = [], []
for f in sorted(glob.glob(f"{R}/*.with_state_and_ring.h5ad")):
    s = os.path.basename(f).split(".")[0]
    o = ad.read_h5ad(f, backed="r").obs.copy()
    o = o[o.in_tissue.astype(int) == 1] if "in_tissue" in o else o
    d_old, _ = dist_from_peaks(o, OLD)
    check.append((s, (zones(d_old) == o.ring.astype(str).values).mean()))
    d, deg = dist_from_peaks(o, HEX)
    for tag, z in [("far3", zones(d, 3)), ("far4_noedge", zones(d, 4, edge=deg < 6)), ("stored", o.ring.astype(str).values)]:
        rec = {"slice": s, "def": tag, "n_Ring1": int((z == "Ring1").sum()), "n_Far": int((z == "Far").sum())}
        for p in PROGS:
            v = o[p].astype(float).values
            far = v[z == "Far"].mean()
            rec[f"{p}_Ring1_minus_Far"] = v[z == "Ring1"].mean() - far
            rec[f"{p}_Ring2_minus_Far"] = v[z == "Ring2"].mean() - far
        rows.append(rec)
    if s.startswith("GSM7990474") or s.startswith("GSM7990478") or s.startswith("GSM7990479"):
        pd.DataFrame({"ring_hex": zones(d, 3)}, index=o.index).to_csv(f"ring_hex_{s}.csv")

res = pd.DataFrame(rows)
z0 = lambda v: (v - v.mean()) / v.std(ddof=0)
res["BarrierIndex"] = np.nan
for tag, g in res.groupby("def"):
    res.loc[g.index, "BarrierIndex"] = z0(g.state_CAF_ECM_Ring1_minus_Far) + z0(g.state_SPP1_TAM_Ring1_minus_Far)
res.to_csv("ring_contrasts_recomputed.csv", index=False)
print("stored labels reproduced with the old offsets (fraction of spots):",
      ", ".join(f"{s[:10]} {x:.3f}" for s, x in check))
summ = []
for tag, g in res.groupby("def"):
    for p in ["state_SPP1_TAM", "state_CAF_ECM", "state_IL10_STAT3"]:
        for r in ["Ring1", "Ring2"]:
            v = g[f"{p}_{r}_minus_Far"].values
            summ.append({"def": tag, "prog": p, "contrast": f"{r}-Far", "n_pos": int((v > 0).sum()),
                         "median": np.median(v), "P": wilcoxon(v).pvalue})
S = pd.DataFrame(summ); S.to_csv("ring_contrasts_recomputed_summary.csv", index=False)
print(S.round(4).to_string(index=False))
print(res.pivot(index="slice", columns="def", values="BarrierIndex").round(2).to_string())
