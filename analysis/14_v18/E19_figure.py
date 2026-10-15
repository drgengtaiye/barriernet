"""E19 figure - robustness of the barriernet scores to gene loss and platform (Fig_E19_robustness).
(a) Spearman correlation of the perturbed with the full NetImmune score against the fraction of the
    48 genes dropped (1,000 random drops per fraction; per-cohort medians and the all-cohort median
    and 5th percentile).
(b) Pooled estimates under dropout (median and 2.5-97.5th percentile of the 1,000 draws), with the
    full-score estimate (dashed) and its 95% interval (band), the prespecified direction threshold
    (dotted), and the panel gene subsets.
(c) Leave-one-gene-out: each pooled estimate after removing one gene."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, GREY, INK, LIGHT, MM, W2

A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
RUN = pd.read_csv(f"{A}/E19_runs.csv")
SP = pd.read_csv(f"{A}/E19_spearman_summary.csv")
LG = pd.read_csv(f"{A}/E19_leave_one_gene_out.csv")
PAN = pd.read_csv(f"{A}/E19_panel_subset.csv").set_index("subset")
AUCCI = pd.read_csv(f"{A}/E7_auc_ci_log2.csv").set_index("Marker").loc["NetImmune"]
FS = 5.2
FULL = RUN[RUN.type == "full"].iloc[0]
FR = [("drop10", 10), ("drop20", 20), ("drop30", 30)]
MODCOL = {"CAF_ECM": STROMA, "TAM_SUP": MACRO, "APM": IMMUNE, "IFN": IMMUNE, "TLS_B": IMMUNE, "T_CELL": IMMUNE}

fig = plt.figure(figsize=(W2, 128 * MM))
gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1.32], hspace=0.5, left=0.07, right=0.99, top=0.95, bottom=0.1)
top = GridSpecFromSubplotSpec(1, 5, subplot_spec=gs[0], width_ratios=[1.35, 1, 1, 1, 1], wspace=0.55)

# ---------------- (a) Spearman vs dropout ----------------
ax = fig.add_subplot(top[0, 0])
coh = SP[(SP.perturbation == "drop10") & ~SP.cohort.str.contains("ALL|median over")].cohort.tolist()
xs = [0, 10, 20, 30]
for c in coh:
    v = [1.0] + [SP[(SP.perturbation == t) & (SP.cohort == c)].median_rho.iloc[0] for t, _ in FR]
    if np.isnan(v[1]):
        continue
    hi = c == "ICB_Limagne2"
    ax.plot(xs, v, color=GREY if hi else LIGHT, lw=0.8 if hi else 0.6, zorder=2 if hi else 1)
    if hi:
        ax.text(30.8, v[-1], "smallest cohort\n(lung, n = 22)", fontsize=FS - 0.6, va="center", ha="left", color=GREY)
allm = [1.0] + [SP[(SP.perturbation == t) & (SP.cohort == "ALL (cohort x draw)")].median_rho.iloc[0] for t, _ in FR]
allp = [1.0] + [SP[(SP.perturbation == t) & (SP.cohort == "ALL (cohort x draw)")].p5_rho.iloc[0] for t, _ in FR]
ax.plot(xs, allm, color=IMMUNE, lw=1.3, marker="o", ms=3.2, mec="white", mew=0.4, zorder=3)
ax.plot(xs, allp, color=IMMUNE, lw=0.8, ls=(0, (3, 1.5)), zorder=3)
ax.axhline(0.9, color=GREY, lw=0.5, ls=(0, (1, 1.5)))
ax.text(0.5, 0.902, "pre-specified 0.9", fontsize=FS - 0.4, color=GREY, va="bottom")
ax.text(30.8, allm[-1] + 0.004, f"median {allm[-1]:.3f}", fontsize=FS - 0.4, color=IMMUNE, va="bottom")
ax.text(30.8, allp[-1] - 0.002, f"5th pct {allp[-1]:.3f}", fontsize=FS - 0.4, color=IMMUNE, va="top")
ax.set_xlim(-1.5, 44); ax.set_ylim(0.86, 1.003); ax.spines["bottom"].set_bounds(0, 30)
ax.set_xticks(xs); ax.set_xticklabels(["0", "10\n(5)", "20\n(10)", "30\n(14)"])
ax.set_xlabel("genes dropped, % (n of 48)", fontsize=FS + 0.4)
ax.set_ylabel("Spearman ρ with full NetImmune", fontsize=FS + 0.4)
ax.set_title("NetImmune, 24 cohorts", fontsize=6, loc="left", pad=4)
ax.spines["right"].set_visible(False)
ax.text(-0.36, 1.08, "a", transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
ax.text(0.02, 0.03, "grey: per-cohort median\nblue: all cohorts, median and 5th pct",
        transform=ax.transAxes, fontsize=FS - 0.6, color=GREY, va="bottom")

# ---------------- (b) pooled estimates vs dropout ----------------
SPEC = [  # column, title, colour, transform, full CI, threshold, threshold label, ylim, log, panel column
    ("HR_NetImmune.est", "NetImmune HR\nbenchmark 3, OS", IMMUNE, np.exp,
     (np.exp(FULL["HR_NetImmune.lo"]), np.exp(FULL["HR_NetImmune.hi"])), 1.0, "HR = 1", (0.5, 1.12), True, "HR_NetImmune"),
    ("AUC_NetImmune_PRJEB", "NetImmune AUC\nPRJEB25780", IMMUNE, lambda x: x,
     (AUCCI.CI_low, AUCCI.CI_high), 0.5, "AUC = 0.5", (0.42, 1.0), False, "AUC_NetImmune_PRJEB"),
    ("OR_Immune.est", "Immune log-OR\n21 checkpoint cohorts", IMMUNE, lambda x: x,
     (FULL["OR_Immune.lo"], FULL["OR_Immune.hi"]), 0.0, "0", (-0.12, 0.72), False, "OR_Immune"),
    ("OR_CAF_ECM.est", "CAF_ECM log-OR\n21 checkpoint cohorts", STROMA, lambda x: x,
     (FULL["OR_CAF_ECM.lo"], FULL["OR_CAF_ECM.hi"]), 0.1, "0.1", (-0.32, 0.3), False, "OR_CAF_ECM"),
]
PX = [("panel_strict", "shared"), ("platform_HTG_OBP", "HTG"), ("platform_NanoString_PanCancerImmune", "PCI")]
for j, (col, ttl, c, tr, ci, thr, thl, yl, logy, pcol) in enumerate(SPEC):
    ax = fig.add_subplot(top[0, j + 1])
    full = tr(FULL[col])
    ax.axhspan(ci[0], ci[1], xmin=0, xmax=1, color=LIGHT, alpha=0.55, lw=0, zorder=0)
    ax.axhline(full, color=c, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
    ax.axhline(thr, color=GREY, lw=0.5, ls=(0, (1, 1.5)), zorder=1)
    ax.plot(0, full, "D", ms=3.4, color=c, mec="white", mew=0.4, zorder=3)
    for i, (t, _) in enumerate(FR, start=1):
        for sfx, dx, filled in (("", -0.13, True), ("_list", 0.17, False)):
            v = tr(RUN.loc[RUN.type == t + sfx, col].values)
            lo, md, hi = np.quantile(v, [0.025, 0.5, 0.975])
            ax.plot([i + dx, i + dx], [lo, hi], color=c if filled else GREY, lw=0.9, solid_capstyle="butt", zorder=2)
            ax.plot(i + dx, md, "o", ms=3.0, color=c if filled else "white", mec=c if filled else GREY,
                    mew=0.4 if filled else 0.7, zorder=3)
    for i, (p, lab) in enumerate(PX):
        x = 4.5 + i * 0.95
        v = PAN.loc[p, pcol]
        if pd.isna(v):
            ax.text(x, yl[0] + 0.06 * (yl[1] - yl[0]) if not logy else yl[0] * 1.04, "n.e.", fontsize=FS - 0.7,
                    color=GREY, ha="center", va="bottom", rotation=90)
        else:
            ax.plot(x, v, "s", ms=3.0, color=c, mec="white", mew=0.4, zorder=3)
    if logy:
        ax.set_yscale("log"); ax.set_yticks([0.5, 0.6, 0.7, 0.8, 1.0]); ax.set_yticklabels(["0.5", "0.6", "0.7", "0.8", "1"])
        ax.minorticks_off()
    ax.set_ylim(*yl); ax.set_xlim(-0.6, 4.5 + 2 * 0.95 + 0.5)
    ax.set_xticks([0, 1, 2, 3] + [4.5 + i * 0.95 for i in range(3)])
    ax.set_xticklabels(["full", "10", "20", "30"] + [l for _, l in PX], fontsize=FS - 0.5)
    for tl in ax.get_xticklabels()[4:]:
        tl.set_fontsize(FS - 0.8); tl.set_rotation(90)
    ax.axvline(3.75, color=LIGHT, lw=0.5, zorder=0)
    ax.text(3.95, thr, thl, fontsize=FS - 0.7, color=GREY, ha="right", va="bottom" if j != 3 else "top")
    ax.set_title(ttl, fontsize=6, loc="left", pad=4, color=c)
    ax.spines["right"].set_visible(False)
    if j == 0:
        ax.text(-0.5, 1.08, "b", transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
    ax.set_xlabel("% dropped      panel set", fontsize=FS - 0.4, labelpad=2)

# ---------------- (c) leave-one-gene-out ----------------
bot = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[1], hspace=0.22)
ROWS = [("HR_NetImmune", "NetImmune HR\n(benchmark 3)", FULL.HR_NetImmune, "HR_k", 3, "{:.3f}"),
        ("AUC_PRJEB", "NetImmune AUC\n(PRJEB25780)", FULL.AUC_NetImmune_PRJEB, None, None, "{:.3f}"),
        ("OR_Immune", "Immune log-OR\n(21 ICB)", FULL["OR_Immune.est"], "Immune_k", 21, "{:+.3f}"),
        ("OR_CAF_ECM", "CAF_ECM log-OR\n(21 ICB)", FULL["OR_CAF_ECM.est"], "CAF_k", 20, "{:+.3f}")]
x = np.arange(len(LG)); mods = LG.module.values
bounds = [i for i in range(1, len(LG)) if mods[i] != mods[i - 1]]
for r, (col, lab, full, kcol, kfull, fmt) in enumerate(ROWS):
    ax = fig.add_subplot(bot[r, 0])
    ax.axhline(full, color=INK, lw=0.5, ls=(0, (3, 1.5)), zorder=1)
    for b in bounds:
        ax.axvline(b - 0.5, color=LIGHT, lw=0.5, zorder=0)
    y = LG[col].values
    for i in range(len(LG)):
        reduced = kcol is not None and LG[kcol].iloc[i] < kfull
        cc = MODCOL[mods[i]]
        ax.plot([i, i], [full, y[i]], color=cc, lw=0.6, alpha=0.6, zorder=1)
        ax.plot(i, y[i], "o", ms=2.8, color="white" if reduced else cc, mec=cc, mew=0.7 if reduced else 0.3, zorder=2)
    span = max(np.abs(y - full).max(), 1e-3) * 1.35
    ax.set_ylim(full - span, full + 1.75 * span)
    ax.set_xlim(-0.8, len(LG) - 0.2)
    ax.set_ylabel(lab, fontsize=FS, rotation=0, ha="right", va="center", labelpad=4)
    ax.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax.tick_params(axis="y", labelsize=FS - 0.5)
    top1 = int(np.argmax(np.abs(y - full)))
    ax.text(1.0, 0.98, f"full {fmt.format(full)}; range {fmt.format(y.min())} to {fmt.format(y.max())}; largest change: {LG.gene.iloc[top1]}",
            transform=ax.transAxes, fontsize=FS - 0.6, color=GREY, ha="right", va="top")
    ax.spines["right"].set_visible(False)
    if r < 3:
        ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
    else:
        ax.set_xticks(x); ax.set_xticklabels(LG.gene, rotation=90, fontsize=FS - 0.9)
        ax.tick_params(axis="x", length=1.5, pad=1)
    if r == 0:
        starts = [0] + bounds; ends = bounds + [len(LG)]
        for s, e in zip(starts, ends):
            ax.text((s + e - 1) / 2, 1.12, mods[s], transform=ax.get_xaxis_transform(), fontsize=FS - 0.3,
                    color=MODCOL[mods[s]], ha="center", va="bottom", fontweight="bold")
        ax.text(-0.075, 1.25, "c", transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")
        ax.text(0.0, 1.47, "Leave-one-gene-out: pooled estimate after removing each gene (dashed, full score)",
                transform=ax.transAxes, fontsize=6, va="bottom", ha="left")

fig.text(0.0, 0.0,
         "b: filled, random drops with every scorable module keeping ≥ 3 genes in every cohort (main); open grey, drops constrained only "
         "on the 48-gene list (sensitivity); bars, 2.5–97.5th percentile of 1,000 draws;\nband, full-score 95% CI (wider of Hartung–Knapp "
         "and DerSimonian–Laird; AUC: bootstrap); squares, scores restricted to panel gene content (shared: the 22 genes measured in all four "
         "panel cohorts;\nCAF_ECM has 1 gene, so NetImmune and CAF_ECM are not estimable, n.e.; HTG, HTG EdgeSeq OBP; PCI, NanoString PanCancer "
         "Immune). c: open circles, removal leaves a module with < 3 genes in ≥ 1 cohort, which then leaves that estimate.",
         fontsize=4.6, color=GREY, ha="left", va="bottom")
bs.save(fig, "Fig_E19_robustness")
print("saved")
