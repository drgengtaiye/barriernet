"""Spot-level ring zones for all 10 sections under the primary definition (barriernet 0.2.1:
TLS peaks = local maxima among the 9 nearest spots above the 95th percentile; rings on the true
Visium lattice; Far >= 3 hops), for Supplementary Fig. S5."""
import sys, glob, os
import numpy as np, pandas as pd, anndata as ad
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/code/barriernet")
from barriernet.spatial import assign_rings
R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2"
out = []
for f in sorted(glob.glob(f"{R}/*.with_state_and_ring.h5ad")):
    a = ad.read_h5ad(f, backed="r"); o = a.obs.copy(); xy = np.asarray(a.obsm["spatial"], float)
    rd = assign_rings(o, xy, "TLS_niche_score", peak_method="knn", peak_k=9, far_min_dist=3)
    keep = o.in_tissue.astype(int).values == 1
    out.append(pd.DataFrame({"slice": os.path.basename(f).split(".")[0], "x": xy[keep, 0], "y": -xy[keep, 1],
                             "zone": rd.ring.values[keep]}))
pd.concat(out).to_parquet("/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/data/rings_knn9_far3.parquet")
print("sections", len(out))
