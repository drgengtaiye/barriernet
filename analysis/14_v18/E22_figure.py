"""E22 figure - gastric and gastro-oesophageal checkpoint cohorts (Fig_E22_gastric_icb).
(a) Log-odds of response (CR/PR vs SD/PD) per SD in PRJEB60680 (Samsung NCT04249739, chemo -> pembrolizumab,
    n = 46 baseline) and PRJEB25780 (Kim 2018, pembrolizumab, n = 44), and pooled (conservative rule: REML
    estimate, wider of Hartung-Knapp and DerSimonian-Laird intervals, larger P); open diamonds add PERFECT.
(b) ROC curves in PRJEB60680 with the structure-matched random-gene null for NetImmune.
(c) MSS/EBV-negative patients of PRJEB60680 (pre-registered sub-hypothesis), beside all patients.
(d) GSE183924 (adjuvant durvalumab, EAC/GEJ, n = 37): relapse-free survival HR per SD.
Style follows E21_fig_trials.py (bnstyle)."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from sklearn.metrics import roc_curve
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, GREY, INK, NULL, fmt_p, panel

import os
A = os.environ.get("E22_OUT", "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
R1 = pd.read_csv(f"{A}/E22_PRJEB60680_response.csv")
IN = pd.read_csv(f"{A}/E22_meta_inputs.csv")
P = pd.read_csv(f"{A}/E22_pooled.csv")
N1 = pd.read_csv(f"{A}/E22_PRJEB60680_null.csv").set_index("score")
ND = pd.read_csv(f"{A}/E22_PRJEB60680_null_draws.csv")
C2 = pd.read_csv(f"{A}/E22_GSE183924_cox.csv").set_index("score")
SC = pd.read_csv(f"{A}/E22_scores_PRJEB60680.csv")
FC = pd.read_csv(f"{A}/E22_firth_cps.csv")
PRIM = "PRJEB60680 + PRJEB25780 (primary)"; SENS = "PRJEB60680 + PRJEB25780 + PERFECT (sensitivity)"
ALL = "primary: all QC-pass baseline, z within library batch"; MSS = "MSS/EBV-negative (pre-registered subgroup)"

ROWS = [("CAF_ECM", "CAF_ECM", STROMA), ("Barrier", "Barrier", STROMA), ("TAM_SUP", "TAM_SUP", MACRO),
        ("IFN", "IFN", IMMUNE), ("APM", "APM", IMMUNE), ("TLS_B", "TLS_B", IMMUNE), ("T_CELL", "T_CELL", IMMUNE),
        ("Immune", "Immune", IMMUNE), ("NetImmune", "NetImmune", IMMUNE)]
NR = len(ROWS); FS = 5.2


def ypos(i):
    return NR - 1 - i


def rowaxis(ax, labels=True):
    ax.set_ylim(-0.6, NR - 0.4)
    if labels:
        ax.set_yticks([ypos(i) for i in range(NR)]); ax.set_yticklabels([r[1] for r in ROWS], fontsize=5.6)
        for lab, (_, _, c) in zip(ax.get_yticklabels(), ROWS):
            lab.set_color(c)
    else:
        ax.set_yticks([])
    ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", labelsize=5.2)
    for i in (3, 7):
        ax.axhline(ypos(i) + 0.5, color=NULL, lw=0.5)


def xaxis(ax, lim, ticks, ref=0.0, log=False):
    ax.set_xlim(*lim)
    if log:
        ax.set_xscale("log"); ax.set_xticks(ticks); ax.set_xticklabels([f"{t:g}" for t in ticks]); ax.minorticks_off()
    else:
        ax.set_xticks(ticks); ax.set_xticklabels([("−" if t < 0 else "") + f"{abs(t):g}" for t in ticks])
    ax.axvline(ref, color=GREY, lw=0.5, ls=(0, (2, 2)))


def point(ax, v, lo, hi, y, c, mk="o", ms=3.0, open_=False, lw=0.8, lim=None):
    if lim is not None:
        a, b = lim
        if lo < a:
            ax.plot(a, y, marker="<", ms=2.2, color=c, mew=0, clip_on=False); lo = a
        if hi > b:
            ax.plot(b, y, marker=">", ms=2.2, color=c, mew=0, clip_on=False); hi = b
    ax.plot([lo, hi], [y, y], color=c, lw=lw, solid_capstyle="butt")
    ax.plot(v, y, mk, ms=ms, color=c, mfc="white" if open_ else c, mec=c if open_ else "white",
            mew=0.7 if open_ else 0.4)


def head(ax, ttl, sub):
    ax.set_title(ttl, fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom")


def cohort_col(ax, cohort, ttl, sub, labels, lim=(-2.1, 3.1), ticks=(-2, -1, 0, 1, 2, 3)):
    T = IN[IN.cohort == cohort].set_index("score")
    for i, (k, _, c) in enumerate(ROWS):
        r = T.loc[k]
        point(ax, r.y, r.y - 1.96 * r.se, r.y + 1.96 * r.se, ypos(i), c, lim=lim)
        if r.p_cohort < 0.05:
            ax.plot(r.y, ypos(i), "o", ms=5.2, mfc="none", mec=INK, mew=0.5)
    xaxis(ax, lim, list(ticks)); head(ax, ttl, sub); rowaxis(ax, labels)
    ax.set_xlabel("log-OR of response per SD", fontsize=5.2, labelpad=1)


def pooled_col(ax):
    P1 = P[P.set == PRIM].set_index("score"); P2 = P[P.set == SENS].set_index("score")
    tr = ax.get_yaxis_transform(); lim = (-2.1, 3.1)
    for i, (k, _, c) in enumerate(ROWS):
        y = ypos(i); a, b = P1.loc[k], P2.loc[k]
        point(ax, a.est, a.lo, a.hi, y + 0.14, c, mk="D", ms=3.2, lim=lim)
        point(ax, b.est, b.lo, b.hi, y - 0.18, c, mk="D", ms=2.6, open_=True, lw=0.6, lim=lim)
        ax.text(1.2, y, fmt_p(a.p, lead=""), transform=tr, clip_on=False, fontsize=FS, ha="right", va="center",
                fontweight="bold" if a.p < 0.05 else "normal")
        ax.text(1.42, y, fmt_p(b.p, lead=""), transform=tr, clip_on=False, fontsize=FS, ha="right", va="center",
                color=GREY, fontweight="bold" if b.p < 0.05 else "normal")
    xaxis(ax, lim, [-2, -1, 0, 1, 2, 3]); head(ax, "Pooled", "random effects, conservative interval")
    ax.text(1.2, 1.012, "P", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    rowaxis(ax, False); ax.set_xlabel("log-OR of response per SD", fontsize=5.2, labelpad=1)
    ax.plot([], [], "D", color=INK, ms=3, label="PRJEB60680 + PRJEB25780")
    ax.plot([], [], "D", mfc="white", mec=INK, ms=2.6, mew=0.7, label="+ PERFECT (sensitivity; grey P)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1, fontsize=4.9, handletextpad=0.3,
              frameon=False, labelspacing=0.3)


def roc_panel(ax):
    y = SC.y.values
    for k, c, ls in [("NetImmune", IMMUNE, "-"), ("IFN", IMMUNE, (0, (3, 1.5))), ("CAF_ECM", STROMA, "-")]:
        f, t, _ = roc_curve(y, SC[k].values)
        a = R1[(R1.analysis == ALL) & (R1.score == k)].iloc[0]
        ax.plot(f, t, color=c, ls=ls, lw=0.9, label=f"{k}  {a.AUC:.2f} ({a.AUC_lo:.2f}–{a.AUC_hi:.2f})")
    ax.plot([0, 1], [0, 1], color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xlabel("1 − specificity", fontsize=5.4); ax.set_ylabel("Sensitivity", fontsize=5.4)
    ax.tick_params(labelsize=5.2)
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.2), fontsize=4.7, handlelength=1.8, title="AUC (95% CI, DeLong)",
              title_fontsize=4.7, frameon=False)
    n = N1.loc["NetImmune"]
    ax.text(-0.05, -0.62, f"NetImmune vs structure-matched random genes:\nnull 95th percentile {n.null_p95:.2f}; P = {n.p_emp:.2g}",
            transform=ax.transAxes, fontsize=4.7, va="top", color=INK)
    head(ax, "PRJEB60680 ROC", f"{int(y.sum())} CR/PR vs {int((1 - y).sum())} SD/PD")


def mss_panel(ax):
    T1 = R1[R1.analysis == MSS].set_index("score"); T0 = R1[R1.analysis == ALL].set_index("score")
    lim = (-1.6, 2.6); tr = ax.get_yaxis_transform()
    for i, (k, _, c) in enumerate(ROWS):
        y = ypos(i); a, b = T1.loc[k], T0.loc[k]
        point(ax, a.logOR, a.lo, a.hi, y + 0.14, c, lim=lim)
        point(ax, b.logOR, b.lo, b.hi, y - 0.18, c, ms=2.6, open_=True, lw=0.6, lim=lim)
        ax.text(1.22, y + 0.14, fmt_p(a.p, lead=""), transform=tr, clip_on=False, fontsize=FS, ha="right", va="center",
                fontweight="bold" if a.p < 0.05 else "normal")
    xaxis(ax, lim, [-1, 0, 1, 2]); rowaxis(ax, True)
    n = T1.loc["NetImmune"]
    head(ax, "MSS/EBV-negative", f"n = {int(n.n)}, {int(n.n_resp)} CR/PR")
    ax.text(1.22, 1.012, "P", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    ax.set_xlabel("log-OR of response per SD", fontsize=5.2, labelpad=1)
    f = FC[(FC.set == MSS) & (FC.term == "NetImmune per SD") & FC.analysis.str.contains("primary")].iloc[0]
    ax.plot([], [], "o", color=INK, ms=3, label="MSS/EBV-negative")
    ax.plot([], [], "o", mfc="white", mec=INK, ms=2.6, mew=0.7, label="all QC-pass patients")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=4.9, handletextpad=0.3,
              frameon=False, columnspacing=0.8)
    ax.text(0.0, -0.3, f"NetImmune beyond log(CPS + 1), Firth: {f.estimate:+.2f} per SD; penalised LRT P = {f.p:.2g}".replace("-", "−"),
            transform=ax.transAxes, fontsize=4.7, color=INK, va="top")


def rfs_panel(ax):
    lim = (0.25, 4.0); tr = ax.get_yaxis_transform()
    for i, (k, _, c) in enumerate(ROWS):
        r = C2.loc[k]; y = ypos(i)
        point(ax, r.HR, r.lo, r.hi, y, c, mk="s", ms=2.8, lim=lim)
        ax.text(1.22, y, fmt_p(r.p, lead=""), transform=tr, clip_on=False, fontsize=FS, ha="right", va="center",
                fontweight="bold" if r.p < 0.05 else "normal")
    xaxis(ax, lim, [0.25, 0.5, 1, 2, 4], ref=1.0, log=True); rowaxis(ax, True)
    r = C2.loc["NetImmune"]
    head(ax, "GSE183924 relapse-free survival", f"adjuvant durvalumab, n = {int(r.n)}, {int(r.events)} events")
    ax.text(1.22, 1.012, "P", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    ax.set_xlabel("HR per SD (< 1 = longer RFS)", fontsize=5.2, labelpad=1)


fig = plt.figure(figsize=(bs.W2, 5.3))
gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1.05], hspace=0.62)
top = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0], wspace=0.12)
r60 = IN[IN.cohort == "PRJEB60680"].iloc[0]; r25 = IN[IN.cohort == "PRJEB25780"].iloc[0]
a1 = fig.add_subplot(top[0]); cohort_col(a1, "PRJEB60680", "PRJEB60680 (Samsung)",
                                         f"5-FU/platinum, then + pembrolizumab; n = {r60.n}, {r60.n_resp} CR/PR", True)
a2 = fig.add_subplot(top[1]); cohort_col(a2, "PRJEB25780", "PRJEB25780 (Kim 2018)",
                                         f"pembrolizumab; n = {r25.n}, {r25.n_resp} CR/PR", False)
a3 = fig.add_subplot(top[2]); pooled_col(a3)
panel(a1, "a", dx=-0.32)
bot = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1], wspace=0.75, width_ratios=[1, 1.05, 1.05])
b = fig.add_subplot(bot[0]); roc_panel(b)
c = fig.add_subplot(bot[1]); mss_panel(c); panel(c, "c", dx=-0.42)
d = fig.add_subplot(bot[2]); rfs_panel(d); panel(d, "d", dx=-0.42)
fig.canvas.draw()                                  # align the b letter with c and d (ROC axes is square)
yb = c.transAxes.transform((0, 1.055))[1]; xb = b.transAxes.transform((-0.3, 0))[0]
fig.text(*fig.transFigure.inverted().transform((xb, yb)), "b", fontsize=8.2, fontweight="bold", va="top", ha="left", color=INK)
if "E22_OUT" in os.environ:                       # test run: keep the figure out of figures/
    for ext in ("pdf", "png"):
        fig.savefig(f"{A}/Fig_E22_gastric_icb.{ext}", dpi=400 if ext == "png" else None)
    print("written (test)", A)
else:
    bs.save(fig, "Fig_E22_gastric_icb")
    print("written figures/Fig_E22_gastric_icb.{png,pdf}")
