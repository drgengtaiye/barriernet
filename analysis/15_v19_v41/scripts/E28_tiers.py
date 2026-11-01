"""E28 Q4 (E28_prereg.md): evidence tiers and multiplicity for the randomized-comparison results. No new
model: pooled results from E21_meta.csv (manuscript rule) and bootstrap-calibrated P from E28 Q3."""
import numpy as np, pandas as pd
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
M = pd.read_csv(f"{O}/E21_meta.csv"); Q3 = pd.read_csv(f"{O}/E28_Q3_calibration_pooled.csv").set_index("score")
TIER = {"tier 1: three pre-registered comparisons (CheckMate 025, JAVELIN Renal 101, NeoTRIP)": "3 new primary only",
        "first three comparisons (IMmotion150, I-SPY2 x 2), which motivated the extension": "reproduction: 3 published primary",
        "tier 2: all six": "primary: 6 comparisons"}
NINE = ["CAF_ECM", "TAM_SUP", "Barrier", "IFN", "Immune", "NetImmune", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1"]


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (m - r) * p[i]); adj[i] = min(1.0, run)
    return adj


def bh(p):
    p = np.asarray(p, float); o = np.argsort(p)[::-1]; m = len(p); adj = np.empty(m); run = 1.0
    for r, i in enumerate(o):
        run = min(run, p[i] * m / (m - r)); adj[i] = run
    return adj


rows = []
for lab, st in TIER.items():
    x = M[(M.set == st) & M.score.isin(NINE)].set_index("score").loc[NINE]
    for sc in NINE:
        r = x.loc[sc]
        rows.append(dict(tier=lab, score=sc, k=r.k, n_benefit=r.k_benefit, est=r.est, lo=r.lo, hi=r.hi, p=r.p,
                         stouffer_p=r.stouffer_p))
T = pd.DataFrame(rows)
six = T[T.tier.str.startswith("tier 2")].set_index("score")
ax = ["Barrier", "Immune", "NetImmune"]
six.loc[ax, "holm_3_axes"] = holm(six.loc[ax, "p"])
cal = [Q3.loc[s, "p_calibrated"] if s in Q3.index else six.loc[s, "p"] for s in ax]
six.loc[ax, "p_calibrated"] = cal; six.loc[ax, "holm_3_axes_calibrated"] = holm(cal)
six["bh_9"] = bh(six.p.values); six["bh_9_stouffer"] = bh(six.stouffer_p.values)
T = T.merge(six[["holm_3_axes", "p_calibrated", "holm_3_axes_calibrated", "bh_9", "bh_9_stouffer"]].reset_index(),
            on="score", how="left")
T.loc[~T.tier.str.startswith("tier 2"), ["holm_3_axes", "p_calibrated", "holm_3_axes_calibrated", "bh_9", "bh_9_stouffer"]] = np.nan
T.to_csv(f"{O}/E28_Q4_tiers.csv", index=False)
pd.set_option("display.width", 250)
print(T.round(4).to_string(index=False))
