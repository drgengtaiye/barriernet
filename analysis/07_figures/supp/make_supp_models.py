"""Supplementary Figs S9-S10: forest plots of the six random-effects meta-analyses, and the three
attempted upgrades (cell-state axes, abundance/state elastic net, spatially trained arrangement score)."""
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import pearsonr
import bnstyle as bs, figstyle as fs
from bnstyle import INK, GREY, STROMA, MACRO, IMMUNE, FITTED, NULL, title, fmt_p

F = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures"; S = f"{F}/supp"
R = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"; R2 = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
G = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2/BarrierIndex_by_slice_with_clinical.csv"
GC = dict(zip(pd.read_csv(G).slice, pd.read_csv(G).GC_raw))
REPORT = {}


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(f"{S}/{stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)


def minus(s): return s.replace("-", "−")


# ---------------- S9 ----------------
uni = pd.read_csv(f"{R}/B4_univariate.csv"); meta = pd.read_csv(f"{R}/B4_meta.csv")
SC = [("NetImmune", IMMUNE), ("Barrier", STROMA), ("Immune", IMMUNE), ("CD8A", INK), ("ESTIMATE(imm-str)", INK), ("NetImmune(ssGSEA)", IMMUNE)]
NAME = {"ESTIMATE(imm-str)": "ESTIMATE immune−stromal", "NetImmune(ssGSEA)": "NetImmune, ssGSEA scoring", "Immune": "Immune (no Barrier term)"}
COH = [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254")]
fig, axs = plt.subplots(2, 3, figsize=(bs.W2, 78 * bs.MM))
for k, (ax, (sc, col)) in enumerate(zip(axs.flat, SC)):
    for i, (c, cn) in enumerate(COH):
        r = uni[(uni.cohort == c) & (uni.score == sc)].iloc[0]
        y = 3 - i
        ax.plot([r.lo, r.hi], [y, y], color=col, lw=0.8)
        ax.plot(r.HR, y, "s", color=col, ms=2.4 + 1.8 * (1 / r.se) / 15)
    m = meta[meta.score == sc].iloc[0]
    ax.fill([m.lo, m.HR, m.hi, m.HR], [0, 0.22, 0, -0.22], color=col)
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_yticks([3, 2, 1, 0]); ax.set_yticklabels([cn for _, cn in COH] + ["Pooled"], fontsize=5.4)
    ax.tick_params(axis="y", length=0)
    ax.set_xscale("log"); ax.set_xticks([0.6, 0.8, 1, 1.25, 1.5]); ax.set_xticklabels(["0.6", "0.8", "1", "1.25", "1.5"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter()); ax.set_xlim(0.55, 1.6)
    ax.set_xlabel("HR per SD" if k >= 3 else "")
    # three-line title so it stays within its own column and clear of the next panel letter
    title(ax, NAME.get(sc, sc) + "\n" + minus(f"pooled {m.HR:.2f} ({m.lo:.2f}–{m.hi:.2f})") + "\n" + minus(f"{fmt_p(m.p)}, I² = {m.I2:.0f}%"))
    bs.panel(ax, "abcdef"[k], dx=-0.5, dy=1.4)
fig.subplots_adjust(hspace=0.95, wspace=0.6)
save(fig, "FigS9")

# ---------------- S10 ----------------
L = pd.read_csv(f"{R2}/S3_state_loadings.csv", index_col=0)
AX = [("Fibro_S2", "Fibroblast state axis 2"), ("Macro_S2", "Macrophage state axis 2"), ("Tcell_S1", "T-cell state axis 1"), ("NK_S1", "NK-cell state axis 1")]
NAMED = {"Fibro_S2": ["CFD", "TNXB", "C3", "PTGDS", "FAP", "TNC", "IL11"], "Macro_S2": ["CXCL5", "CCL18", "MMP9", "IL1A", "F13A1", "P2RY13"],
         "Tcell_S1": ["TNFRSF4", "TNFRSF18", "IL2RA", "IFNG", "CRTAM"], "NK_S1": ["PDCD1", "LAG3", "CXCR3"]}
coef = pd.read_csv(f"{R2}/S5_coefs.csv", index_col=0)
loco = pd.read_csv(f"{R2}/S5_loco_cindex.csv").set_index("model")
pred = pd.read_csv(f"{R2}/S8_neigh_pred.csv")
w = pd.read_csv(f"{R2}/S9_arrangement_weights.csv", index_col=0).iloc[:, 0]
sas = pd.read_csv(f"{R2}/S10_prognosis.csv")
fig = plt.figure(figsize=(bs.W2, 150 * bs.MM))
gsp = fig.add_gridspec(3, 4, height_ratios=[1.25, 1.05, 1], hspace=0.8, wspace=0.95)
for k, (a, lab) in enumerate(AX):
    ax = fig.add_subplot(gsp[0, k])
    v = L[a].dropna()
    keep = set(v.nlargest(5).index) | set(v.nsmallest(5).index) | set(NAMED[a])
    top = v[list(keep)].sort_values()
    ax.barh(range(len(top)), top.values, color=[GREY if x < 0 else INK for x in top.values], height=0.7)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=5.0, style="italic"); ax.tick_params(axis="y", length=0)
    for tl in ax.get_yticklabels():
        if tl.get_text() in NAMED[a]: tl.set_fontweight("bold")
    ax.axvline(0, color=INK, lw=0.4); ax.set_xlabel("loading", fontsize=5.4)
    title(ax, lab)
    if k == 0: bs.panel(ax, "a", dx=-0.6)
# b: selected elastic-net coefficients (abundance + state) by held-out cohort
ax = fig.add_subplot(gsp[1, :2])
cf = coef.loc[(coef != 0).any(axis=1)]
im = ax.imshow(cf.T.values, cmap=fs.DIV, vmin=-0.2, vmax=0.2, aspect="auto")
ax.set_xticks(range(len(cf))); ax.set_xticklabels([x.replace("ab_", "abund. ").replace("st_", "state ") for x in cf.index], rotation=60, ha="right", fontsize=5.0)
ax.set_yticks(range(3)); ax.set_yticklabels([f"held out: {c}" for c in cf.columns], fontsize=5.2); ax.tick_params(length=0)
for sp_ in ax.spines.values(): sp_.set_visible(False)
cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.0, length=1.5)
cb.set_label("coefficient", fontsize=5.2)
title(ax, "Elastic net on abundance + state (mean held-out C " + f"{loco.loc['EN: abundance + state (30)', 'mean']:.3f})")
bs.panel(ax, "b", dx=-0.12, dy=1.12)
# c: arrangement score, out-of-section r by section
ax = fig.add_subplot(gsp[1, 2])
rs = pred.groupby("section").apply(lambda g: pearsonr(g.excl_resid, g.pred_resid)[0], include_groups=False)
rs.index = rs.index.map(GC); rs = rs.sort_values()
ax.barh(range(len(rs)), rs.values, color=STROMA, height=0.7)
ax.set_yticks(range(len(rs))); ax.set_yticklabels(rs.index, fontsize=5.2); ax.tick_params(axis="y", length=0)
ax.axvline(0, color=INK, lw=0.4); ax.set_xlabel("r, predicted vs observed\n(held-out section)", fontsize=5.4)
title(ax, minus(f"Median r {rs.median():.2f}"))
bs.panel(ax, "c", dx=-0.45, dy=1.12)
REPORT["S10_arrangement_median_r"] = round(rs.median(), 3)
# d: top arrangement genes
ax = fig.add_subplot(gsp[1, 3])
top = pd.concat([w.nsmallest(6), w.nlargest(6)]).sort_values()
ax.barh(range(len(top)), top.values, color=[STROMA if x > 0 else GREY for x in top.values], height=0.7)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=5.0, style="italic"); ax.tick_params(axis="y", length=0)
ax.axvline(0, color=INK, lw=0.4); ax.set_xlabel("ridge weight", fontsize=5.4)
title(ax, "Arrangement genes")
bs.panel(ax, "d", dx=-0.55, dy=1.12)
# e: bulk performance of the arrangement score
ax = fig.add_subplot(gsp[2, :2])
coh = [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254")]
for i, (c, cn) in enumerate(coh):
    r = sas[sas.cohort == c].iloc[0]
    ax.plot(r.SAS_uni_HR, 2 - i + 0.12, "o", color=INK, ms=3.2)
    ax.plot(r.SAS_adj_HR, 2 - i - 0.12, "s", mfc="white", mec=INK, ms=3.0)
    ax.text(1.88, 2 - i, minus(f"P = {r.SAS_uni_p:.2f}; ΔC {r.dC:+.3f}*"), fontsize=5.0, va="center", ha="right")
ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
ax.set_yticks([2, 1, 0]); ax.set_yticklabels([cn for _, cn in coh]); ax.tick_params(axis="y", length=0)
ax.set_xlim(0.8, 1.9); ax.set_xticks([0.8, 1.0, 1.2, 1.4]); ax.set_xlabel("HR per SD of the arrangement score")
ax.plot([], [], "o", color=INK, ms=3, label="unadjusted"); ax.plot([], [], "s", mfc="white", mec=INK, ms=3, label="adjusted (stage, age)")
ax.legend(loc="upper center", bbox_to_anchor=(0.35, -0.3), ncol=2, fontsize=5.0)
ax.text(1.0, -0.47, "*ΔC over NetImmune, stage and age", transform=ax.transAxes, fontsize=5.0, ha="right", va="top", color=GREY)
title(ax, "Arrangement score in bulk cohorts")
bs.panel(ax, "e", dx=-0.12, dy=1.12)
# f: LOCO C of the upgrade models
ax = fig.add_subplot(gsp[2, 2:])
rows = [("EN: abundance + state (30)", "(A) abundance + state", FITTED), ("EN: abundance only (10)", "(A) abundance only", FITTED),
        ("NetImmune (fixed)", "NetImmune (fixed)", IMMUNE), ("EN: 6 barriernet modules", "(B) learned module weights", FITTED),
        ("EN: state only (20)", "(A) state only", FITTED)]
for i, (k, lab, col) in enumerate(rows):
    y = len(rows) - 1 - i
    for j, (c, mk) in enumerate([("TCGA", "o"), ("GSE84437", "^"), ("GSE62254", "v")]):
        ax.plot(loco.loc[k, c], y + 0.22, mk, color=GREY, ms=2.6, mew=0)
    ax.plot(loco.loc[k, "mean"], y - 0.12, "D", color=col, ms=3.2)
    ax.text(0.655, y, f"{loco.loc[k, 'mean']:.3f}", fontsize=5.0, va="center", ha="right", color=col)
ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([r[1] for r in rows], fontsize=5.2); ax.tick_params(axis="y", length=0)
ax.set_xlim(0.48, 0.66); ax.set_xlabel("held-out Harrell's C (grey: cohorts; diamond: mean)")
title(ax, "Leave-one-cohort-out C")
bs.panel(ax, "f", dx=-0.5, dy=1.12)
save(fig, "FigS10")
print(REPORT)
