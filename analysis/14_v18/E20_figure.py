"""E20 figure - extended locked replication in gastric cohorts (Fig_E20_extended_replication).
(a) HR per SD for NetImmune and CAF_ECM: benchmark cohorts, locked cohorts (first 3 and added 4), the
    sensitivity cohort GSE57303, and pooled estimates with the conservative (wider of HK/DL) interval.
(b) Structure-matched random-gene null for NetImmune in the added locked cohorts and pooled.
(c) Paired difference in Harrell's C, NetImmune minus ESTIMATE immune-stromal and minus CD8A.
Layout and style follow E15_figure.py."""
import sys, io, subprocess
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, GREY, INK, NULL, MM, W2, fmt_p

A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
RS = "/users/yc4229/.conda/envs/r45/bin/Rscript"
PRIM = "pre-registered exclusions (primary)"
P15 = pd.read_csv(f"{A}/E15_per_cohort.csv"); P15 = P15[P15.model == "unadjusted"]
PC = pd.read_csv(f"{A}/E20_per_cohort.csv")
PC = PC[(PC.model == "unadjusted") & (PC.endpoint == "OS") & (PC.duplicate_rule == PRIM)]
PL = pd.read_csv(f"{A}/E20_pooled.csv"); PL = PL[(PL.model == "unadjusted") & (PL.duplicate_rule == PRIM)]
OR = pd.read_csv(f"{A}/E8_univariate.csv")
NU = pd.read_csv(f"{A}/E20_null.csv").set_index("level")
DRAWS = pd.read_csv(io.StringIO(subprocess.run(
    [RS, "-e", f"write.csv(readRDS('{A}/E20_null_draws.rds'), stdout(), row.names=FALSE)"],
    capture_output=True, text=True, check=True).stdout))
EX = pd.read_csv(f"{A}/E20_duplicate_exclusions.csv")
N_EX_ADDED = int((EX.cohort != "GSE57303").sum()); N_EX_SENS = int((EX.cohort == "GSE57303").sum())
DC = pd.read_csv(f"{A}/E20_deltaC.csv")
FS = 5.2
BENCH = [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254")]
FIRST = [("GSE15459", "GSE15459"), ("GSE34942", "GSE34942"), ("GSE26253", "GSE26253 (RFS)")]
ADDED = [("GSE13861", "GSE13861"), ("GSE26901", "GSE26901"), ("GSE26899", "GSE26899"), ("GSE29272", "GSE29272")]
SENS = [("GSE57303", "GSE57303")]
POOLS = [("benchmark 3 (reference)", "Pooled, benchmark 3"),
         ("locked (first 3), E15 reference", "Pooled, locked (first 3)"),
         ("locked (added 4)", "Pooled, locked (added 4)"),
         ("locked 7 (primary)", "Pooled, locked 7 (primary)"),
         ("all 10 gastric", "Pooled, all 10 gastric")]
YS = [17.6, 16.6, 15.6,            # benchmark
      13.9, 12.9, 11.9,            # locked (first 3)
      10.2, 9.2, 8.2, 7.2,         # locked (added)
      5.5,                         # sensitivity
      3.8, 2.8, 1.8, 0.8, -0.2]    # pooled
XMIN, XMAX = 0.28, 3.6


def forest(ax, score, col, labels):
    rows = []
    for k, lab in BENCH:
        r = OR[(OR.cohort == k) & (OR.score == score)].iloc[0]
        rows.append((lab, np.exp(r.beta), np.exp(r.beta - 1.96 * r.se), np.exp(r.beta + 1.96 * r.se), None, "bench"))
    for k, lab in FIRST:
        r = P15[(P15.cohort == k) & (P15.score == score)].iloc[0]
        rows.append((f"{lab}  {int(r.n)}/{int(r.events)}", r.HR, r.lo, r.hi, r.p, "locked"))
    for k, lab in ADDED:
        r = PC[(PC.cohort == k) & (PC.score == score)].iloc[0]
        rows.append((f"{lab}  {int(r.n)}/{int(r.events)}", r.HR, r.lo, r.hi, r.p, "locked"))
    for k, lab in SENS:
        r = PC[(PC.cohort == k) & (PC.score == score)].iloc[0]
        rows.append((f"{lab}  {int(r.n)}/{int(r.events)}", r.HR, r.lo, r.hi, r.p, "sens"))
    for k, lab in POOLS:
        r = PL[(PL.set == k) & (PL.score == score)].iloc[0]
        rows.append((lab, r.HR, r.lo, r.hi, r.p, "pool"))
    for (lab, hr, lo, hi, p, kind), y in zip(rows, YS):
        c = GREY if kind in ("bench", "sens") else col
        ax.plot([lo, hi], [y, y], color=c, lw=0.9, solid_capstyle="butt")
        if kind == "pool":
            primary = "primary" in lab
            ax.plot(hr, y, "D", ms=4.6 if primary else 4.0, color=INK if primary else c, mec="white", mew=0.4)
        elif kind == "bench":
            ax.plot(hr, y, "o", ms=3.4, color=c, mec="white", mew=0.4)
        elif kind == "sens":
            ax.plot(hr, y, "s", ms=3.2, mfc="white", mec=GREY, mew=0.8)
        else:
            ax.plot(hr, y, "s", ms=3.4, color=c, mec="white", mew=0.4)
        if kind != "bench":
            ax.text(XMAX * 1.08, y, f"{hr:.2f} ({lo:.2f}–{hi:.2f})  {fmt_p(p)}", fontsize=FS, va="center", ha="left",
                    color=INK if kind != "sens" else GREY, clip_on=False,
                    fontweight="bold" if (kind == "pool" and p < 0.05) else "normal")
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xscale("log"); ax.set_xlim(XMIN, XMAX)
    ax.set_xticks([0.3, 0.5, 0.7, 1, 1.4, 2, 3]); ax.set_xticklabels(["0.3", "0.5", "0.7", "1", "1.4", "2", "3"])
    ax.minorticks_off()
    ax.set_ylim(-0.9, 18.6)
    ax.set_yticks(YS); ax.set_yticklabels([r[0] for r in rows] if labels else [""] * len(rows), fontsize=FS)
    for t, (_, _, _, _, _, kind) in zip(ax.get_yticklabels(), rows):
        if kind == "pool": t.set_fontweight("bold")
        if kind == "sens": t.set_color(GREY)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("HR per SD (95% CI)", fontsize=FS + 0.4)
    ax.set_title(score, fontsize=6, loc="left", color=col, fontweight="bold", pad=4)
    if labels:
        for y, t in [(18.35, "benchmark cohorts"), (14.65, "locked (first 3)  n/events"),
                     (10.95, "locked (added)  n/events"), (6.25, "sensitivity (secondary OS source)"),
                     (4.55, "pooled, random effects")]:
            ax.text(-0.02, y, t, transform=ax.get_yaxis_transform(), ha="right", fontsize=FS, color=GREY)


fig = plt.figure(figsize=(W2, 138 * MM))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.55, 1], width_ratios=[1, 1], hspace=0.26, wspace=1.05,
              left=0.19, right=0.84, top=0.955, bottom=0.09)
ax1 = fig.add_subplot(gs[0, 0]); forest(ax1, "NetImmune", IMMUNE, True)
ax2 = fig.add_subplot(gs[0, 1]); forest(ax2, "CAF_ECM", STROMA, False)
fig.text(0.005, 0.985, "a", fontsize=8, fontweight="bold", va="top")

# (b) null
lev = [("GSE13861", "GSE13861"), ("GSE26901", "GSE26901"), ("GSE26899", "GSE26899"), ("GSE29272", "GSE29272"),
       ("locked (added 4), pooled", "locked (added 4), pooled"), ("locked 7, pooled", "locked 7, pooled")]
sub = GridSpecFromSubplotSpec(len(lev), 1, subplot_spec=gs[1, 0], hspace=0.3)
for i, (k, lab) in enumerate(lev):
    ax = fig.add_subplot(sub[i, 0])
    v = DRAWS[k].values
    ax.hist(v, bins=np.linspace(0, 0.8, 65), color=NULL, lw=0)
    r = NU.loc[k]
    ax.axvline(r.null_p95_abs, color=GREY, lw=0.6)
    ax.plot(abs(r.observed_logHR), ax.get_ylim()[1] * 0.35, "D", ms=3.6, color=IMMUNE if "pooled" not in k else INK,
            mec="white", mew=0.4)
    ax.text(0.99, 0.97, f"{lab}  {fmt_p(r.P)}", transform=ax.transAxes, ha="right", va="top", fontsize=FS,
            color=INK, fontweight="bold" if r.P < 0.05 else "normal",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.9), zorder=5)
    ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xlim(0, 0.8)
    if i < len(lev) - 1: ax.set_xticklabels([])
    if i == 0:
        ax.set_title("NetImmune vs 1,000 random axes (locked cohorts)", fontsize=6, loc="left", pad=4)
    if i == len(lev) - 1: ax.set_xlabel("absolute log HR per SD", fontsize=FS + 0.4)
fig.text(0.005, 0.425, "b", fontsize=8, fontweight="bold", va="top")

# (c) delta C
ax = fig.add_subplot(gs[1, 1])
rows = []
DLEV = [("GSE13861", "GSE13861"), ("GSE26901", "GSE26901"), ("GSE26899", "GSE26899"), ("GSE29272", "GSE29272"),
        ("locked (added 4), pooled (inverse variance)", "pooled, locked (added 4)"),
        ("locked 7, pooled (inverse variance)", "pooled, locked 7")]
for cmp in ["ESTIMATE", "CD8A"]:
    for k, lab in DLEV:
        r = DC[(DC.cohort == k) & (DC.comparator == cmp)].iloc[0]
        rows.append((f"vs {cmp}: {lab}", r.dC, r.lo, r.hi, r.p, k.startswith("locked")))
ys = [13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]
for (lab, d, lo, hi, p, pooled), y in zip(rows, ys):
    ax.plot([lo, hi], [y, y], color=INK if pooled else GREY, lw=0.9)
    ax.plot(d, y, "D" if pooled else "o", ms=4 if pooled else 3.2, color=IMMUNE if pooled else GREY, mec="white", mew=0.4)
    if pooled:
        ax.text(0.33, y, f"{d:+.3f} ({lo:+.3f} to {hi:+.3f})  {fmt_p(p)}", fontsize=FS, va="center", ha="left",
                clip_on=False, fontweight="bold" if p < 0.05 else "normal")
ax.set_yticks(ys); ax.set_yticklabels([r[0] for r in rows], fontsize=FS)
for t, r in zip(ax.get_yticklabels(), rows):
    if r[5]: t.set_fontweight("bold")
ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
ax.set_xlim(-0.2, 0.3); ax.set_ylim(0.3, 13.7)
ax.set_xlabel("ΔC, NetImmune minus comparator", fontsize=FS + 0.4)
ax.set_title("Harrell's C, paired bootstrap (overall survival)", fontsize=6, loc="left", pad=4)
fig.text(0.47, 0.425, "c", fontsize=8, fontweight="bold", va="top")
fig.text(0.0, 0.0, "Squares, locked cohorts; circles, benchmark cohorts; open square, sensitivity cohort (survival from a third-party "
         "curation); diamonds, random-effects pooled estimates with the wider of the Hartung–Knapp and DerSimonian–Laird intervals. "
         "GSE26253 endpoint is recurrence-free survival;\nall others overall survival. Arrays flagged by the expression-correlation "
         "duplicate rule are excluded ("
         f"{N_EX_ADDED} in the added locked cohorts, {N_EX_SENS} in GSE57303). "
         "Locked 7 = locked (first 3) + locked (added 4); "
         "all 10 adds the benchmark cohorts.",
         fontsize=4.8, color=GREY, ha="left", va="bottom")
bs.save(fig, "Fig_E20_extended_replication")
print("saved")
