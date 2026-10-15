"""Fig_E12 - checkpoint response across 21 cohorts: the 11 of E5 and 10 added in E12.
Per cohort: log-odds of response per SD (95% CI), five scores side by side. The added
(replication) cohorts sit on a tinted band in the score colour; the E5 cohorts are grey.
Pooled rows: DerSimonian-Laird estimate with its 95% CI (diamond) and the Hartung-Knapp CI
(thin bar beneath) for the E5 11, the 10 new cohorts alone, and all 21."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Polygon, Rectangle
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, GREY, INK, LIGHT, MM, W2, fmt_p

A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
per = pd.read_csv(f"{A}/E12_icb_per_cohort.csv")
pool = pd.read_csv(f"{A}/E12_icb_pooled.csv")

SCORES = [("CAF_ECM", "CAF_ECM", STROMA), ("TAM_SUP", "TAM_SUP", MACRO), ("IFN", "IFN", IMMUNE),
          ("Immune", "Immune", IMMUNE), ("NetImmune", "NetImmune", IMMUNE)]
NAME = {"ICB_Braun": "Braun", "ICB_Gide": "Gide", "ICB_Hugo": "Hugo", "ICB_Jung": "Jung", "ICB_Limagne1": "Limagne 1",
        "ICB_Liu": "Liu", "ICB_Mariathasan": "Mariathasan", "ICB_Miao1": "Miao", "ICB_Puch": "Puch",
        "ICB_Riaz": "Riaz", "ICB_Van_Allen": "Van Allen", "ICB_Nathanson": "Nathanson", "ICB_Snyder": "Snyder",
        "ICB_Limagne2": "Limagne 2", "ICB_Shiuan": "Shiuan", "GSE176307": "GSE176307", "GSE161537": "GSE161537",
        "GSE159067": "GSE159067", "GSE215868": "GSE215868", "GSE93157": "GSE93157",
        "brca_iatlas_anders_2022": "Anders (iAtlas)"}
TUM = {"kidney": "renal", "melanoma": "melanoma", "lung": "NSCLC", "bladder": "urothelial", "urothelial": "urothelial",
       "head and neck": "HNSCC", "mixed (lung, HNSCC, melanoma)": "mixed", "breast": "breast"}
PANEL = {"GSE161537", "GSE159067", "GSE215868", "GSE93157"}

old = per[per.set == "E5"]; new = per[per.set == "new"]
meta = per[per.set.isin(["E5", "new"])].drop_duplicates("cohort").set_index("cohort")
order_old = (old[old.score == "Immune"].assign(t=lambda d: d.tumour).sort_values(["t", "cohort"]).cohort.tolist())
order_new = (new[new.score == "Immune"].assign(t=lambda d: d.tumour).sort_values(["t", "cohort"]).cohort.tolist())
POOLS = [("E5 11 (reference)", "First 11, pooled", GREY), ("new only (replication)", "Replication 10, pooled", None),
         ("all", "All 21, pooled", INK)]

# y positions, top to bottom
rows = []
y = 0
rows.append(("head", "First 11 cohorts", y)); y -= 1
for c in order_old: rows.append(("old", c, y)); y -= 1
y -= 0.4; rows.append(("head", "Replication cohorts", y)); y -= 1
for c in order_new: rows.append(("new", c, y)); y -= 1
y -= 0.4
for key, lab, _ in POOLS: rows.append(("pool", key, y)); y -= 1.55
YMIN = y + 0.4
XL = (-1.6, 2.6)


def lab_of(kind, c):
    if kind == "head": return c
    if kind == "pool": return dict((k, l) for k, l, _ in POOLS)[c]
    m = meta.loc[c]
    star = "†" if c in PANEL else ""
    return f"{NAME[c]}{star}  {TUM.get(m.tumour, m.tumour)}  {int(m.n)}/{int(m.n_resp)}"


def draw(ax, k, col, first):
    for kind, c, yy in rows:
        if kind == "new":
            ax.add_patch(Rectangle((XL[0], yy - 0.5), XL[1] - XL[0], 1.0, color=col, alpha=0.07, lw=0, zorder=0))
    for kind, c, yy in rows:
        if kind in ("old", "new"):
            r = per[(per.cohort == c) & (per.score == k)]
            if r.empty:
                ax.text(0, yy, "n.c.", fontsize=5.0, color=GREY, ha="center", va="center"); continue
            r = r.iloc[0]; lo, hi = r.logOR - 1.96 * r.se, r.logOR + 1.96 * r.se
            cc = col if kind == "new" else GREY
            ax.plot([max(lo, XL[0]), min(hi, XL[1])], [yy, yy], color=cc, lw=0.6, solid_capstyle="butt")
            for edge, cond, mk in ((XL[0], lo < XL[0], "<"), (XL[1], hi > XL[1], ">")):
                if cond: ax.plot(edge, yy, mk, color=cc, ms=2.2, mew=0, clip_on=False)
            ms = float(np.clip(1.3 / r.se, 1.6, 4.2))
            ax.plot(r.logOR, yy, "s" if kind == "new" else "o", color=cc, ms=ms, mec="white", mew=0.3, zorder=3)
        elif kind == "pool":
            q = pool[(pool.set == c) & (pool.score == k)]
            if q.empty: continue
            q = q.iloc[0]; pc = dict((kk, cl) for kk, _, cl in POOLS)[c] or col
            h = 0.34
            ax.add_patch(Polygon([[q.lo_DL, yy], [q.logOR, yy + h], [q.hi_DL, yy], [q.logOR, yy - h]],
                                 closed=True, fc=pc, ec="none", zorder=3))
            ax.plot([q.lo_HK, q.hi_HK], [yy - 0.47] * 2, color=pc, lw=0.5)
            ax.text(XL[1], yy + 0.02, f"{q.logOR:+.2f}  {int(q.n_pos)}/{int(q.k)}", fontsize=5.0, ha="right",
                    va="bottom", color=INK, fontweight="bold" if max(q.p_DL, q.p_HK) < 0.05 else "normal")
            ax.text(XL[1], yy - 0.05, fmt_p(max(q.p_DL, q.p_HK)), fontsize=5.0, ha="right", va="top", color=GREY)
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.set_xlim(*XL); ax.set_ylim(YMIN, 0.8)
    ax.set_xticks([-1, 0, 1, 2]); ax.set_xlabel("log-odds of response per SD", fontsize=5.6)
    ax.spines["left"].set_visible(False); ax.set_yticks([])
    bs.title(ax, k, pad=2); ax.title.set_color(col)
    if first:
        tr = ax.get_yaxis_transform()
        for kind, c, yy in rows:
            ax.text(-0.04, yy, lab_of(kind, c), transform=tr, ha="right", va="center", clip_on=False,
                    fontsize=5.4 if kind == "head" else 5.0, fontweight="bold" if kind in ("head", "pool") else "normal",
                    color=INK if kind != "old" else "#4A545E")
        ax.text(-0.04, 0.9, "cohort  tumour  n/responders", transform=tr, ha="right", va="bottom", fontsize=5.0,
                color=GREY, clip_on=False)


fig = plt.figure(figsize=(W2 - 1.2 * MM, 126 * MM))
gs = GridSpec(1, 6, figure=fig, width_ratios=[1.32, 1, 1, 1, 1, 1], wspace=0.10, left=0.0, right=1.0, bottom=0.1)
for i, (k, _, col) in enumerate(SCORES):
    ax = fig.add_subplot(gs[0, i + 1]); draw(ax, k, col, i == 0)
fig.text(0.0, 0.0, "Squares, replication cohorts; circles, first 11 cohorts (larger marker, smaller SE). Diamond, DerSimonian\u2013Laird "
         "pooled log-OR with 95% CI; bar beneath, Hartung\u2013Knapp 95% CI;\nnumbers, pooled log-OR, cohorts > 0 of k, "
         "and the larger of the DL and HK P values. \u2020 targeted panel. n.c., module not computable (< 3 genes measured).",
         fontsize=5.0, color=GREY, ha="left", va="bottom", linespacing=1.3)
bs.save(fig, "Fig_E12_icb_expanded")
print("saved")
