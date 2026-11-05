"""Re-derive the manuscript numbers that no analysis script prints directly (found while mapping the text to its
scripts), from the saved results; nothing is refitted. Writes results/text_numbers_check.csv.

1. Supplementary Table S26c 'predicted' column (E36_S4_power_vs_theory.csv): normal-theory power of the paired test,
   Phi(oracle_diff / median_se - 1.96), for each module as its own truth (computed in-session after E36, post hoc).
2. 'the share of signatures above random sets ranged from 3 of 52 in NeoTRIP to 28 of 52 in JAVELIN Renal 101':
   per comparison, signatures whose interaction exceeds the 95th percentile of all that comparison's random draws
   (as figures/cns/v34/make_fig3.py panel b, which first printed it).
3. 'three comparisons could detect only differences of 0.6 to 0.7 per SD': 2.80 x the median paired standard error of
   the signatures against GZMB, IDO1 and CXCL10 in the breast comparisons (E35_paired.csv).
"""
import numpy as np, pandas as pd
from scipy.stats import norm

O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
rows = []
# 1
pv = pd.read_csv(f"{O}/E36_S4_power_vs_theory.csv")
S4 = pd.read_csv(f"{O}/E36_S4.csv"); own = S4[S4.truth == S4.module].set_index("module")
pred = pd.Series(norm.cdf(own.oracle_diff / own.median_se - 1.959964), index=own.index)
d1 = float(np.max(np.abs(pv.set_index("module").predicted - pred.reindex(pv.module).values)))
for m, p in pred.items():
    rows.append(dict(item="S26c predicted power", key=m, value=round(float(p), 4), check=f"max |diff| vs file {d1:.1e}"))
assert d1 < 1e-3, d1
# 2
PER = pd.read_csv(f"{O}/E32_per_comparison.csv"); PER = PER[PER.variant == "primary"]
LIB = set(pd.read_csv(f"{O}/E32_primary.csv").query("kind == 'library'").signature)
Z = np.load(f"{O}/E32_reference_draws.npz")
for k in ["IMM150", "CM025", "JR101", "ISPY2_pembro", "ISPY2_durva", "NeoTRIP"]:
    ys = np.concatenate([Z[key][:, 0] for key in Z.files if key.startswith(k + "|")]); ys = ys[np.isfinite(ys)]
    hi95 = np.percentile(ys, 95); g = PER[(PER.comparison == k) & PER.signature.isin(LIB)]
    rows.append(dict(item="signatures above the 95th percentile of random sets", key=k, value=f"{int((g.y > hi95).sum())}/{len(g)}",
                     check="as figures/cns/v34/Fig3_numbers.md panel b"))
# 3
P = pd.read_csv(f"{O}/E35_paired.csv"); B = P[P.set == "breast three"]
for gname in ["GZMB", "IDO1", "CXCL10"]:
    rows.append(dict(item="breast detectable difference (80% power)", key=gname,
                     value=round(2.8016 * float(B[B.reference == gname].se.median()), 3), check="2.80 x median paired SE"))
out = pd.DataFrame(rows); out.to_csv(f"{O}/text_numbers_check.csv", index=False); print(out.to_string(index=False))
