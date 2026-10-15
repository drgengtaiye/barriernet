"""E15 figure - locked validation in the newly downloaded gastric cohorts (Fig_E15_gastric_validation).
(a) HR per SD for NetImmune and CAF_ECM: the three original cohorts, the three new cohorts, and
    pooled estimates (original 3, new 3, all 6) with the conservative (wider of HK/DL) interval.
(b) Structure-matched random-gene null for NetImmune in the new cohorts (per cohort and pooled).
(c) Paired difference in Harrell's C, NetImmune minus ESTIMATE immune-stromal and minus CD8A."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, GREY, INK, NULL, MM, W2, fmt_p

A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
PC = pd.read_csv(f"{A}/E15_per_cohort.csv"); PC = PC[PC.model == "unadjusted"]
PL = pd.read_csv(f"{A}/E15_pooled.csv"); PL = PL[PL.model == "unadjusted"]
OR = pd.read_csv(f"{A}/E8_univariate.csv")
NU = pd.read_csv(f"{A}/E15_null.csv")
import subprocess, io
DRAWS = pd.read_csv(io.StringIO(subprocess.run(
    ["Rscript", "-e", f"write.csv(readRDS('{A}/E15_null_draws.rds'), stdout(), row.names=FALSE)"],
    capture_output=True, text=True, check=True).stdout))
DC = pd.read_csv(f"{A}/E15_deltaC.csv")
FS = 5.2
ORIG = [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254")]
NEW = [("GSE15459", "GSE15459"), ("GSE34942", "GSE34942"), ("GSE26253", "GSE26253 (RFS)")]
POOLS = [("original 3 (reference)", "Pooled, benchmark 3"), ("new cohorts", "Pooled, locked 3"),
         ("new + original 3", "Pooled, all 6")]


def forest(ax, score, col, labels):
    rows = []
    for k, lab in ORIG:
        r = OR[(OR.cohort == k) & (OR.score == score)].iloc[0]
        rows.append((lab, np.exp(r.beta), np.exp(r.beta - 1.96 * r.se), np.exp(r.beta + 1.96 * r.se), None, "orig"))
    for k, lab in NEW:
        r = PC[(PC.cohort == k) & (PC.score == score)].iloc[0]
        rows.append((f"{lab}  {int(r.n)}/{int(r.events)}", r.HR, r.lo, r.hi, r.p, "new"))
    for k, lab in POOLS:
        r = PL[(PL.set == k) & (PL.score == score)].iloc[0]
        rows.append((lab, r.HR, r.lo, r.hi, r.p, "pool"))
    ys = [9.2, 8.2, 7.2, 5.6, 4.6, 3.6, 2.0, 1.0, 0.0]
    for (lab, hr, lo, hi, p, kind), y in zip(rows, ys):
        c = GREY if kind == "orig" else col
        ax.plot([lo, hi], [y, y], color=c, lw=0.9, solid_capstyle="butt")
        if kind == "pool":
            ax.plot(hr, y, "D", ms=4.2, color=INK if "all" in lab else c, mec="white", mew=0.4)
            ax.text(3.2, y, f"{hr:.2f} ({lo:.2f}–{hi:.2f})  {fmt_p(p)}", fontsize=FS, va="center", ha="left",
                    color=INK, fontweight="bold" if p < 0.05 else "normal", clip_on=False)
        else:
            ax.plot(hr, y, "o" if kind == "orig" else "s", ms=3.4, color=c, mec="white", mew=0.4)
            if kind == "new":
                ax.text(3.2, y, f"{hr:.2f} ({lo:.2f}–{hi:.2f})  {fmt_p(p)}", fontsize=FS, va="center", ha="left",
                        color=INK, clip_on=False)
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xscale("log"); ax.set_xlim(0.45, 2.6)
    ax.set_xticks([0.5, 0.7, 1, 1.4, 2]); ax.set_xticklabels(["0.5", "0.7", "1", "1.4", "2"])
    ax.minorticks_off()
    ax.set_ylim(-0.7, 10.2)
    ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows] if labels else [""] * len(rows), fontsize=FS)
    for t, (_, _, _, _, _, kind) in zip(ax.get_yticklabels(), rows):
        if kind == "pool": t.set_fontweight("bold")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("HR per SD (95% CI)", fontsize=FS + 0.4)
    ax.set_title(score, fontsize=6, loc="left", color=col, fontweight="bold", pad=4)
    if labels:
        ax.text(-0.02, 9.95, "benchmark cohorts", transform=ax.get_yaxis_transform(), ha="right", fontsize=FS, color=GREY)
        ax.text(-0.02, 6.35, "locked cohorts  n/events", transform=ax.get_yaxis_transform(), ha="right", fontsize=FS, color=GREY)


fig = plt.figure(figsize=(W2, 100 * MM))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.15, 1], width_ratios=[1, 1], hspace=0.55, wspace=0.9,
              left=0.17, right=0.86, top=0.94, bottom=0.12)
ax1 = fig.add_subplot(gs[0, 0]); forest(ax1, "NetImmune", IMMUNE, True)
ax2 = fig.add_subplot(gs[0, 1]); forest(ax2, "CAF_ECM", STROMA, False)
fig.text(0.005, 0.975, "a", fontsize=8, fontweight="bold", va="top")

# (b) null
sub = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[1, 0], hspace=0.25)
lev = [("GSE15459", "GSE15459"), ("GSE34942", "GSE34942"), ("GSE26253", "GSE26253"), ("pooled", "locked 3, pooled")]
for i, (k, lab) in enumerate(lev):
    ax = fig.add_subplot(sub[i, 0])
    v = DRAWS[k].values
    ax.hist(v, bins=np.linspace(0, 0.6, 49), color=NULL, lw=0)
    r = NU.iloc[i]
    ax.axvline(r.null_p95_abs, color=GREY, lw=0.6)
    ax.plot(abs(r.observed_logHR), ax.get_ylim()[1] * 0.35, "D", ms=3.6, color=IMMUNE, mec="white", mew=0.4)
    ax.text(0.99, 0.62, f"{lab}  {fmt_p(r.P)}", transform=ax.transAxes, ha="right", fontsize=FS,
            color=INK, fontweight="bold" if r.P < 0.05 else "normal")
    ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xlim(0, 0.6)
    if i < 3: ax.set_xticklabels([])
    if i == 0:
        ax.set_title("NetImmune vs 1,000 random axes (locked cohorts)", fontsize=6, loc="left", pad=4)
        fig.text(0.005, 0.44, "b", fontsize=8, fontweight="bold", va="top")
    if i == 3: ax.set_xlabel("absolute log HR per SD", fontsize=FS + 0.4)

# (c) delta C
ax = fig.add_subplot(gs[1, 1])
rows = []
for cmp, col in [("ESTIMATE", STROMA), ("CD8A", GREY)]:
    for k in ["GSE15459", "GSE34942", "GSE26253", "pooled (inverse variance)"]:
        r = DC[(DC.cohort == k) & (DC.comparator == cmp)].iloc[0]
        rows.append((f"vs {cmp}: {'pooled' if k.startswith('pooled') else k}", r.dC, r.lo, r.hi, r.p, k.startswith("pooled")))
for i, (lab, d, lo, hi, p, pooled) in enumerate(rows[::-1]):
    ax.plot([lo, hi], [i, i], color=INK if pooled else GREY, lw=0.9)
    ax.plot(d, i, "D" if pooled else "o", ms=4 if pooled else 3.2, color=IMMUNE if pooled else GREY, mec="white", mew=0.4)
    if pooled:
        ax.text(0.2, i, f"{d:+.3f} ({lo:+.3f} to {hi:+.3f})  {fmt_p(p)}", fontsize=FS, va="center", ha="left",
                clip_on=False, fontweight="bold" if p < 0.05 else "normal")
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=FS)
ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
ax.set_xlim(-0.2, 0.2); ax.set_xlabel("ΔC, NetImmune minus comparator", fontsize=FS + 0.4)
ax.set_title("Harrell's C, paired bootstrap", fontsize=6, loc="left", pad=4)
fig.text(0.47, 0.44, "c", fontsize=8, fontweight="bold", va="top")
fig.text(0.0, 0.0, "Squares, locked cohorts; circles, benchmark cohorts; diamonds, random-effects pooled estimates with the wider of the "
         "Hartung–Knapp and DerSimonian–Laird intervals. GSE26253 endpoint is recurrence-free survival.",
         fontsize=4.8, color=GREY, ha="left", va="bottom")
bs.save(fig, "Fig_E15_gastric_validation")
print("saved")
