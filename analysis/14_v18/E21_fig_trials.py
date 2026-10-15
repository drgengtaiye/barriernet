"""E21 figure - score x checkpoint-arm interaction across six randomised comparisons (Fig_E21_trials).
(a) Benefit-oriented interaction per SD (-log interaction HR for PFS, +log interaction OR for pCR;
    right = more benefit from the checkpoint-containing arm) in each primary comparison, and the
    pooled estimate (REML; wider of the Hartung-Knapp and DerSimonian-Laird intervals, larger P),
    with PRINCE and DUTRENEO added as a sensitivity analysis.
(b) The two sensitivity-only comparisons.
(c) Robustness of the pooled NetImmune and CAF_ECM interaction: leave one trial out, tumour-type
    subsets, sensitivity sets and alternative endpoint / adjustment.
(d) Structure-matched random-gene nulls for the NetImmune and CAF_ECM interactions in the three
    newly added primary trials. Style follows E10_fig_trials.py; no internal labels."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, GREY, INK, NULL, fmt_p


_fmt_p_shared = fmt_p


def fmt_p(p, lead="P = "):
    """As the shared formatter, but keep trailing zeros between 0.001 and 0.1 (0.060, not 0.06), as in the text."""
    if p == p and p is not None and 0.001 <= p < 0.1:
        return f"{lead}{p:.3f}" if p >= 0.01 else f"{lead}{p:.4f}"
    return _fmt_p_shared(p, lead=lead)

A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
IN = pd.read_csv(f"{A}/E21_meta_inputs.csv")
M = pd.read_csv(f"{A}/E21_meta.csv")
ND = pd.read_csv(f"{A}/E21_null_draws.csv"); NS = pd.read_csv(f"{A}/E21_null_summary.csv")
PRIM = "primary: 6 comparisons"; SENS = "+ PRINCE + DUTRENEO (8)"

ROWS = [("CAF_ECM", "CAF_ECM", STROMA), ("Barrier", "Barrier", STROMA),
        ("TAM_SUP", "TAM_SUP", MACRO), ("C1Q_arm", "C1Q arm", MACRO), ("SPP1_arm", "SPP1 arm", MACRO),
        ("CXCL9_SPP1", "CXCL9:SPP1", MACRO),
        ("IFN", "IFN", IMMUNE), ("APM", "APM", IMMUNE), ("TLS_B", "TLS_B", IMMUNE), ("T_CELL", "T_CELL", IMMUNE),
        ("Immune", "Immune", IMMUNE), ("NetImmune", "NetImmune", IMMUNE)]
NR = len(ROWS)
FS = 5.2


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
    for i in range(NR):
        if i in (2, 6):
            ax.axhline(ypos(i) + 0.5, color=NULL, lw=0.5)


def xaxis(ax, lim, ticks):
    ax.set_xlim(*lim); ax.set_xticks(ticks)
    ax.set_xticklabels([("−" if t < 0 else "") + f"{abs(t):g}" for t in ticks])
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))


def point(ax, v, lo, hi, y, c, mk="o", ms=3.0, open_=False, lw=0.8, lim=None):
    if lim is not None:                       # clip, with an arrowhead where the interval runs off the axis
        a, b = lim
        if lo < a:
            ax.plot(a, y, marker="<", ms=2.2, color=c, mew=0); lo = a
        if hi > b:
            ax.plot(b, y, marker=">", ms=2.2, color=c, mew=0); hi = b
    ax.plot([lo, hi], [y, y], color=c, lw=lw, solid_capstyle="butt")
    ax.plot(v, y, mk, ms=ms, color=c, mfc="white" if open_ else c, mec=c if open_ else "white",
            mew=0.7 if open_ else 0.4)


def trial_col(ax, trial, lim, ticks, ttl, sub, color=None, labels=False):
    T = IN[(IN.set == SENS) & (IN.trial == trial)].set_index("score")
    for i, (k, _, c) in enumerate(ROWS):
        r = T.loc[k]; cc = color or c
        point(ax, r.y, r.y - 1.96 * r.se, r.y + 1.96 * r.se, ypos(i), cc, lim=lim)
        if r.p_trial < 0.05:
            ax.plot(r.y, ypos(i), "o", ms=5.2, mfc="none", mec=INK, mew=0.5)
    xaxis(ax, lim, ticks)
    ax.set_title(ttl, fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom")
    rowaxis(ax, labels)
    ax.set_xlabel("per SD", fontsize=5.2, labelpad=1)


def pooled(ax):
    P1 = M[M.set == PRIM].set_index("score"); P2 = M[M.set == SENS].set_index("score")
    tr = ax.get_yaxis_transform(); lim = (-0.65, 0.8)
    for i, (k, _, c) in enumerate(ROWS):
        y = ypos(i); a, b = P1.loc[k], P2.loc[k]
        point(ax, a.est, a.lo, a.hi, y + 0.13, c, mk="D", ms=3.2, lim=lim)
        point(ax, b.est, b.lo, b.hi, y - 0.17, c, mk="D", ms=2.6, open_=True, lw=0.6, lim=lim)
        ax.text(1.27, y, f"{int(a.k_benefit)}/6", transform=tr, clip_on=False, fontsize=FS, ha="right", va="center")
        ax.text(1.66, y, fmt_p(a.p, lead=""), transform=tr, clip_on=False, fontsize=FS, ha="right", va="center",
                fontweight="bold" if a.p < 0.05 else "normal")
    xaxis(ax, lim, [-0.5, 0, 0.5])
    ax.set_title("Pooled", fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.text(0, 1.012, "random effects", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom")
    ax.text(1.27, 1.012, "benefit", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    ax.text(1.66, 1.012, "P", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    rowaxis(ax, False)
    ax.set_xlabel("per SD", fontsize=5.2, labelpad=1)
    ax.plot([], [], "D", color=INK, ms=3, label="6 comparisons (primary)")
    ax.plot([], [], "D", mfc="white", mec=INK, ms=2.6, mew=0.7, label="+ PRINCE, DUTRENEO")
    ax.legend(loc="upper center", bbox_to_anchor=(0.75, -0.13), ncol=1, fontsize=4.9, handletextpad=0.3,
              frameon=False, labelspacing=0.3)


ROB = [(PRIM, "All six (primary)"),
       ("primary without IMmotion150", "− IMmotion150"),
       ("primary without CheckMate 025", "− CheckMate 025"),
       ("primary without JAVELIN Renal 101", "− JAVELIN Renal 101"),
       ("primary without I-SPY2 pembrolizumab", "− I-SPY2 pembro."),
       ("primary without I-SPY2 durvalumab-olaparib", "− I-SPY2 durva.–olap."),
       ("primary without NeoTRIP", "− NeoTRIP"),
       ("renal only (IMmotion150, CheckMate 025, JAVELIN)", "Renal only (3)"),
       ("breast only (I-SPY2 x 2, NeoTRIP)", "Breast only (3)"),
       (SENS, "+ PRINCE, DUTRENEO (8)"),
       ("primary, CheckMate 025 OS for PFS", "CheckMate 025 OS"),
       ("primary, NeoTRIP without batch", "NeoTRIP, no batch term")]


def robustness(ax):
    n = len(ROB); lim = (-1.25, 1.35)
    tr = ax.get_yaxis_transform()
    for i, (s, lab) in enumerate(ROB):
        y = n - 1 - i
        for k, c, dy in [("NetImmune", IMMUNE, 0.17), ("CAF_ECM", STROMA, -0.17)]:
            r = M[(M.set == s) & (M.score == k)].iloc[0]
            point(ax, r.est, r.lo, r.hi, y + dy, c, mk="D", ms=2.8, lim=lim)
            if r.p < 0.05:
                ax.plot(r.est, y + dy, "o", ms=4.8, mfc="none", mec=INK, mew=0.5)
        r = M[(M.set == s) & (M.score == "NetImmune")].iloc[0]
        ax.text(1.24, y + 0.17, fmt_p(r.p, lead=""), transform=tr, clip_on=False, fontsize=4.9, ha="right", va="center",
                color=IMMUNE, fontweight="bold" if r.p < 0.05 else "normal")
        r = M[(M.set == s) & (M.score == "CAF_ECM")].iloc[0]
        ax.text(1.24, y - 0.17, fmt_p(r.p, lead=""), transform=tr, clip_on=False, fontsize=4.9, ha="right", va="center",
                color=STROMA, fontweight="bold" if r.p < 0.05 else "normal")
        if i in (1, 7, 9):
            ax.axhline(y + 0.5, color=NULL, lw=0.5)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks([n - 1 - i for i in range(n)]); ax.set_yticklabels([r[1] for r in ROB], fontsize=5.2)
    ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", labelsize=5.2)
    xaxis(ax, lim, [-1, -0.5, 0, 0.5, 1])
    ax.set_xlabel("Pooled benefit-oriented log ratio per SD", fontsize=5.4)
    ax.text(1.24, 1.012, "P", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    ax.set_title("Robustness of the pooled interaction", fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.plot([], [], "D", color=IMMUNE, ms=2.8, label="NetImmune")
    ax.plot([], [], "D", color=STROMA, ms=2.8, label="CAF_ECM")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, fontsize=4.9, frameon=False,
              columnspacing=0.8, handletextpad=0.3)


def nulls(ax):
    items = [("NetImmune", "CM025", "CheckMate 025"), ("NetImmune", "JR101", "JAVELIN Renal 101"),
             ("NetImmune", "NeoTRIP", "NeoTRIP"), ("NetImmune", "combined", "three trials, pooled"),
             ("CAF_ECM", "CM025", "CheckMate 025"), ("CAF_ECM", "JR101", "JAVELIN Renal 101"),
             ("CAF_ECM", "NeoTRIP", "NeoTRIP"), ("CAF_ECM", "combined", "three trials, pooled")]
    hi = 0.9; bins = np.linspace(0, hi, 46)
    S = NS.set_index(["score", "contrast"])
    for i, (k, key, lab) in enumerate(items):
        yb = len(items) - 1 - i + (0.45 if i < 4 else 0)
        c = IMMUNE if k == "NetImmune" else STROMA
        h, e = np.histogram(np.clip(np.abs(ND[f"{k}_{key}"]), 0, hi - 1e-9), bins=bins); h = h / h.max() * 0.62
        ax.fill_between(e[:-1], yb, yb + h, step="post", color=NULL, lw=0)
        s = S.loc[(k, key)]
        ax.plot([s.null_p95_abs] * 2, [yb, yb + 0.7], color=GREY, lw=0.6)
        ax.plot(abs(s.obs_y), yb + 0.2, "D" if s.obs_y >= 0 else "d", color=c, ms=3.0, mec="white", mew=0.4)
        ax.text(hi * 0.985, yb + 0.35, fmt_p(s.p_emp), fontsize=FS, ha="right", va="center",
                fontweight="bold" if s.p_emp < 0.05 else "normal")
        ax.text(hi * 0.985, yb + 0.80, lab + ("" if s.obs_y >= 0 else " (against)"), fontsize=5.0, ha="right",
                va="center", color=INK)
    ax.text(0.01, 8.25, "NetImmune", fontsize=5.4, ha="left", va="center", color=IMMUNE, fontweight="bold")
    ax.text(0.01, 3.82, "CAF_ECM", fontsize=5.4, ha="left", va="center", color=STROMA, fontweight="bold")
    ax.set_ylim(-0.15, 8.5); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, hi); ax.tick_params(axis="x", labelsize=5.2)
    ax.set_xlabel("Absolute benefit-oriented interaction per SD", fontsize=5.4)
    ax.set_title("Structure-matched null, added trials", fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.text(0, 1.012, "direction credited; (against): observed favours the control arm", transform=ax.transAxes,
            fontsize=4.9, color=GREY, va="bottom")
    ax.plot([], [], "s", color=NULL, ms=4, label="2,000 random sets")
    ax.plot([], [], "|", color=GREY, ms=5, label="95th pct.")
    ax.plot([], [], "D", color=INK, ms=3, label="observed")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=3, fontsize=4.9, frameon=False,
              columnspacing=0.6, handletextpad=0.3)


fig = plt.figure(figsize=(bs.W2, 172 * bs.MM))
gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 0.98], hspace=0.44, left=0.085, right=0.955, top=0.925, bottom=0.075)
top = gs[0].subgridspec(1, 7, width_ratios=[1, 1, 1, 1, 1, 1, 1.05], wspace=0.13)
ax_a = [fig.add_subplot(top[0, j]) for j in range(7)]
REN = ((-0.9, 0.9), [-0.5, 0, 0.5]); BRE = ((-1.8, 1.8), [-1, 0, 1])
trial_col(ax_a[0], "IMmotion150", *REN, "IMmotion150", "PFS, −log HR", labels=True)
trial_col(ax_a[1], "CheckMate 025", *REN, "CheckMate 025", "PFS, −log HR")
trial_col(ax_a[2], "JAVELIN Renal 101", *REN, "JAVELIN Renal 101", "PFS, −log HR")
trial_col(ax_a[3], "I-SPY2 pembrolizumab", *BRE, "I-SPY2 pembro.", "pCR, log OR")
trial_col(ax_a[4], "I-SPY2 durvalumab-olaparib", *BRE, "I-SPY2 durva.–olap.", "pCR, log OR")
trial_col(ax_a[5], "NeoTRIP", *BRE, "NeoTRIP", "pCR, log OR")
pooled(ax_a[6])
for ax0, ax1, lab in [(ax_a[0], ax_a[2], "Renal cell carcinoma"), (ax_a[3], ax_a[5], "Breast cancer")]:
    b0, b1 = ax0.get_position(), ax1.get_position()
    fig.text((b0.x0 + b1.x1) / 2, b0.y1 + 0.032, lab, fontsize=5.6, ha="center", va="bottom", color=INK)
    fig.add_artist(plt.Line2D([b0.x0, b1.x1], [b0.y1 + 0.029] * 2, color=GREY, lw=0.5))
ax_a[0].text(0, -0.155, "Right of zero: the score predicts more benefit from the checkpoint-containing arm; "
             "ringed points, P < 0.05.", transform=ax_a[0].transAxes, fontsize=5.0, color=GREY, ha="left", va="top")

# bottom row placed in millimetres so the row labels of (c) never run into (b)
yb0 = gs[1].get_position(fig).y0; ybh = gs[1].get_position(fig).height
W = bs.W2 / bs.MM


def at(x0, w):
    return [x0 / W, yb0, w / W, ybh]


ax_b = [fig.add_axes(at(15.5, 20)), fig.add_axes(at(38.5, 20))]
trial_col(ax_b[0], "PRINCE", (-1.2, 1.2), [-1, 0, 1], "PRINCE", "PFS, −log HR", color=GREY, labels=True)
trial_col(ax_b[1], "DUTRENEO", (-2.4, 2.4), [-2, 0, 2], "DUTRENEO", "pCR, log OR", color=GREY)
ax_c = fig.add_axes(at(87, 32)); ax_d = fig.add_axes(at(133, 39))
robustness(ax_c); nulls(ax_d)
for ax, L, dx, dy in [(ax_a[0], "a", -0.75, 1.1), (ax_b[0], "b", -0.75, 1.07), (ax_c, "c", -0.84, 1.07),
                      (ax_d, "d", -0.1, 1.07)]:
    bs.panel(ax, L, dx=dx, dy=dy)
bs.save(fig, "Fig_E21_trials")
print("saved Fig_E21_trials")
