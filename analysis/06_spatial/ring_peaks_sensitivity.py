"""TLS-centred ring contrasts with barriernet 0.2.1 under three TLS-peak definitions.

stored   the TLS_peaks column saved with the v13 analyses (local maxima against the 0.2.0
         broken graph, i.e. at most two diagonal neighbours) - reproduces figures/ring_recompute.py
knn9     local maxima among the 9 nearest spots (the definition stated in the Methods)
lattice  local maxima among the six lattice neighbours
Distances always use the corrected Visium lattice. Far: >= 3 hops (primary) and >= 4 hops
excluding tissue-edge spots (sensitivity). Outputs next to this script.
"""
import sys, glob, os
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import wilcoxon
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "barriernet"))
from barriernet.spatial import assign_rings, barrier_index

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2"
PROGS = ["state_CAF_ECM", "state_SPP1_TAM", "state_IL10_STAT3"]
rows = []
for f in sorted(glob.glob(f"{R}/*.with_state_and_ring.h5ad")):
    a = ad.read_h5ad(f, backed="r"); o = a.obs.copy(); xy = np.asarray(a.obsm["spatial"], float)
    sec = os.path.basename(f).split(".")[0]
    for peaks_def in ["stored", "knn9", "lattice"]:
        kw = dict(peaks=o.TLS_peaks.astype(int).values.astype(bool)) if peaks_def == "stored" else \
             dict(peak_method="knn", peak_k=9) if peaks_def == "knn9" else dict(peak_method="lattice")
        for far_def, fkw in [("far3", dict(far_min_dist=3)), ("far4_noedge", dict(far_min_dist=4, far_exclude_edge=True))]:
            rd = assign_rings(o, xy, "TLS_niche_score", **kw, **fkw)
            z = rd.ring.values
            rec = dict(section=sec, peaks=peaks_def, far=far_def, n_peaks=int(rd.tls_peak.sum()),
                       n_Ring1=int((z == "Ring1").sum()), n_Far=int((z == "Far").sum()))
            for p in PROGS:
                v = o[p].astype(float).values; far = v[z == "Far"].mean()
                rec[f"{p}_Ring1_minus_Far"] = v[z == "Ring1"].mean() - far
                rec[f"{p}_Ring2_minus_Far"] = v[z == "Ring2"].mean() - far
            rows.append(rec)
res = pd.DataFrame(rows)
res["BarrierIndex"] = np.nan
for key, g in res.groupby(["peaks", "far"]):
    res.loc[g.index, "BarrierIndex"] = barrier_index(g, "state_CAF_ECM_Ring1_minus_Far", "state_SPP1_TAM_Ring1_minus_Far")
res.to_csv(os.path.join(HERE, "ring_peaks_sensitivity.csv"), index=False)
summ = []
for (pk, fd), g in res.groupby(["peaks", "far"]):
    for p in PROGS:
        for r in ["Ring1", "Ring2"]:
            v = g[f"{p}_{r}_minus_Far"].values
            summ.append(dict(peaks=pk, far=fd, prog=p.replace("state_", ""), contrast=f"{r}-Far",
                             n_higher=int((v > 0).sum()), median=np.median(v), P=wilcoxon(v).pvalue))
S = pd.DataFrame(summ); S.to_csv(os.path.join(HERE, "ring_peaks_sensitivity_summary.csv"), index=False)
print(res.groupby("peaks").n_peaks.sum().to_string())
print(S[S.far == "far3"].pivot_table(index=["prog", "contrast"], columns="peaks", values="P").round(4).to_string())
print(S[S.far == "far3"].pivot_table(index=["prog", "contrast"], columns="peaks", values="n_higher").to_string())
print(res[res.far == "far3"].pivot(index="section", columns="peaks", values="BarrierIndex").round(2).to_string())
