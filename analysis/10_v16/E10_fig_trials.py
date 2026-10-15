"""E10 figure - score x checkpoint-arm interaction across randomised trials (Fig_E10_trials).
(a) Interaction per SD in each randomised comparison on its own scale, every axis oriented so
    that benefit from the checkpoint arm lies to the right (HR axes reversed), and the pooled
    benefit-oriented log ratio (REML + Hartung-Knapp; PRINCE added as a sensitivity analysis).
(b, c) I-SPY2 within-arm OR for pCR per SD (HR-adjusted), checkpoint arm vs paclitaxel control.
(d) Structure-matched random-gene nulls for the NetImmune and CAF_ECM interactions in I-SPY2."""
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
E5 = pd.read_csv(f"{A}/E5_trials.csv")
R = pd.read_csv(f"{A}/E10_ispy2_interactions.csv")
W = pd.read_csv(f"{A}/E10_ispy2_within_arm.csv"); W = W[W.model == "HRadj"]
M = pd.read_csv(f"{A}/E10_meta.csv")
ND = pd.read_csv(f"{A}/E10_ispy2_null_draws.csv"); NS = pd.read_csv(f"{A}/E10_ispy2_null_summary.csv")

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


def logaxis(ax, lim, ticks, reverse=False):
    ax.set_xscale("log"); ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xticks(ticks); ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.set_xlim(*(lim[::-1] if reverse else lim))


def point(ax, v, lo, hi, y, c, mk="o", ms=3.0, open_=False, lw=0.8):
    ax.plot([lo, hi], [y, y], color=c, lw=lw, solid_capstyle="butt")
    ax.plot(v, y, mk, ms=ms, color=c, mfc="white" if open_ else c, mec=c if open_ else "white",
            mew=0.7 if open_ else 0.4)


def trial_col(ax, get, ref, lim, ticks, reverse, ttl, sub, color=None, labels=False):
    for i, (k, _, c) in enumerate(ROWS):
        v, lo, hi, p = get(k)
        cc = color or c
        point(ax, v, lo, hi, ypos(i), cc)
        if p < 0.05:
            ax.plot(v, ypos(i), "o", ms=5.2, mfc="none", mec=INK, mew=0.5)
    ax.axvline(ref, color=GREY, lw=0.5, ls=(0, (2, 2)))
    logaxis(ax, lim, ticks, reverse)
    ax.set_title(ttl, fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom")
    rowaxis(ax, labels)


def e5(pat):
    T = E5[E5.trial.str.startswith(pat)].set_index("score")
    return lambda k: (T.loc[k, "int_HR"], T.loc[k, "int_lo"], T.loc[k, "int_hi"], T.loc[k, "int_p"])


def e10(key):
    T = R[R.cohort == key].set_index("score")
    return lambda k: (T.loc[k, "HRadj_OR"], T.loc[k, "HRadj_lo"], T.loc[k, "HRadj_hi"], T.loc[k, "HRadj_p"])


def cons(T):
    """Conservative rule used throughout the manuscript: the wider of the Hartung-Knapp (REML) and
    fixed/DerSimonian-Laird intervals (identical when tau^2 = 0), with the larger P value."""
    T = T.copy()
    fe_wider = (T.fe_hi - T.fe_lo) > (T.hi - T.lo)
    T.loc[fe_wider, "lo"] = T.loc[fe_wider, "fe_lo"]; T.loc[fe_wider, "hi"] = T.loc[fe_wider, "fe_hi"]
    T["p"] = T[["p", "fe_p"]].max(axis=1)
    return T


def pooled(ax):
    P1 = cons(M[M.set.str.startswith("primary")].set_index("score")); P2 = cons(M[M.set == "+ PRINCE"].set_index("score"))
    tr = ax.get_yaxis_transform()
    for i, (k, _, c) in enumerate(ROWS):
        y = ypos(i); a, b = P1.loc[k], P2.loc[k]
        point(ax, a.est, a.lo, a.hi, y + 0.13, c, mk="D", ms=3.2)
        point(ax, b.est, b.lo, b.hi, y - 0.17, c, mk="D", ms=2.6, open_=True, lw=0.6)
        ax.text(1.20, y, f"{int(a.k_benefit)}/3", transform=tr, clip_on=False, fontsize=FS, ha="right", va="center")
        ax.text(1.53, y, fmt_p(a.p, lead=""), transform=tr, clip_on=False, fontsize=FS, ha="right", va="center",
                fontweight="bold" if a.p < 0.05 else "normal")
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(-1.25, 1.15); ax.set_xticks([-1, -0.5, 0, 0.5, 1]); ax.set_xticklabels(["−1", "−0.5", "0", "0.5", "1"])
    ax.set_title("Pooled", fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.text(0, 1.012, "benefit-oriented log ratio", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom")
    ax.text(1.20, 1.012, "benefit", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    ax.text(1.53, 1.012, "P", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    rowaxis(ax, False)
    ax.plot([], [], "D", color=INK, ms=3, label="3 comparisons (primary)")
    ax.plot([], [], "D", mfc="white", mec=INK, ms=2.6, mew=0.7, label="+ PRINCE")
    ax.legend(loc="upper center", bbox_to_anchor=(0.62, -0.13), ncol=2, fontsize=4.9, columnspacing=0.6,
              handletextpad=0.3, frameon=False)


def within(ax, key, ttl, lab_ici, labels):
    T = W[W.cohort == key]; I = R[R.cohort == key].set_index("score")
    tr = ax.get_yaxis_transform()
    for i, (k, _, c) in enumerate(ROWS):
        y = ypos(i)
        q = T[(T.score == k) & (T.arm == "control")].iloc[0]; e = T[(T.score == k) & (T.arm == "checkpoint")].iloc[0]
        ax.plot([q.lo, q.hi], [y - 0.17] * 2, color=GREY, lw=0.7)
        ax.plot(q.OR, y - 0.17, "s", mfc="white", mec=GREY, ms=2.7, mew=0.7)
        point(ax, e.OR, e.lo, e.hi, y + 0.17, c)
        r = I.loc[k]
        ax.text(1.20, y, f"{r.HRadj_OR:.2f}", transform=tr, clip_on=False, fontsize=FS, ha="right", va="center",
                fontweight="bold" if r.HRadj_p < 0.05 else "normal")
        ax.text(1.44, y, fmt_p(r.HRadj_p, lead=""), transform=tr, clip_on=False, fontsize=FS, ha="right", va="center")
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
    logaxis(ax, (0.15, 14), [0.2, 0.5, 1, 2, 5, 10])
    ax.set_xlabel("OR for pCR per SD, within arm", fontsize=5.4)
    ax.text(1.20, 1.012, "int. OR", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    ax.text(1.44, 1.012, "P", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom", ha="right")
    ax.set_title(ttl, fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    rowaxis(ax, labels)
    ax.plot([], [], "o", color=INK, ms=3, label=lab_ici)
    ax.plot([], [], "s", mfc="white", mec=GREY, ms=2.7, mew=0.7, label="paclitaxel control")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1, fontsize=4.9, frameon=False,
              handletextpad=0.3, labelspacing=0.3)


def nulls(ax):
    items = [("NetImmune", "ISPY2_pembro", "pembrolizumab"), ("NetImmune", "ISPY2_durva", "durvalumab–olaparib"),
             ("NetImmune", "combined", "both, pooled"),
             ("CAF_ECM", "ISPY2_pembro", "pembrolizumab"), ("CAF_ECM", "ISPY2_durva", "durvalumab–olaparib"),
             ("CAF_ECM", "combined", "both, pooled")]
    bins = np.linspace(0, 1.6, 49)
    S = NS.set_index(["score", "contrast"])
    for i, (k, key, lab) in enumerate(items):
        yb = len(items) - 1 - i + (0.45 if i < 3 else 0)
        c = IMMUNE if k == "NetImmune" else STROMA
        h, e = np.histogram(np.abs(ND[f"{k}_{key}"]), bins=bins); h = h / h.max() * 0.62
        ax.fill_between(e[:-1], yb, yb + h, step="post", color=NULL, lw=0)
        s = S.loc[(k, key)]
        ax.plot([s.null_p95_abs] * 2, [yb, yb + 0.7], color=GREY, lw=0.6)
        ax.plot(abs(s.obs_logOR), yb + 0.2, "D", color=c, ms=3.0, mec="white", mew=0.4)
        ax.text(1.58, yb + 0.35, fmt_p(s.p_emp), fontsize=FS, ha="right", va="center",
                fontweight="bold" if s.p_emp < 0.05 else "normal")
        ax.text(1.58, yb + 0.80, lab, fontsize=5.0, ha="right", va="center", color=INK)
    ax.text(0.02, 6.25, "NetImmune", fontsize=5.4, ha="left", va="center", color=IMMUNE,
            fontweight="bold")
    ax.text(0.02, 2.82, "CAF_ECM", fontsize=5.4, ha="left", va="center", color=STROMA, fontweight="bold")
    ax.set_ylim(-0.15, 6.5); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 1.6); ax.tick_params(axis="x", labelsize=5.2)
    ax.set_xlabel("Absolute interaction log-OR (hormone-receptor-adjusted)", fontsize=5.4)
    ax.set_title("Structure-matched null, I-SPY2", fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.plot([], [], "s", color=NULL, ms=4, label="2,000 random sets")
    ax.plot([], [], "|", color=GREY, ms=5, label="95th pct.")
    ax.plot([], [], "D", color=INK, ms=3, label="observed")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=3, fontsize=4.9, frameon=False,
              columnspacing=0.6, handletextpad=0.3)


fig = plt.figure(figsize=(bs.W2, 158 * bs.MM))
gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 0.98], hspace=0.42, left=0.1, right=0.955, top=0.93, bottom=0.075)
top = gs[0].subgridspec(1, 5, width_ratios=[1, 1, 1, 1, 1.15], wspace=0.14)
ax_a = [fig.add_subplot(top[0, j]) for j in range(5)]
trial_col(ax_a[0], e5("IMmotion150: atezolizumab (+/- bev)"), 1, (0.4, 1.9), [0.5, 0.7, 1, 1.4], True,
          "IMmotion150 (renal)", "PFS interaction HR, axis reversed", labels=True)
trial_col(ax_a[1], e10("ISPY2_pembro"), 1, (0.15, 6), [0.2, 0.5, 1, 2, 5], False,
          "I-SPY2 pembrolizumab", "pCR interaction OR")
trial_col(ax_a[2], e10("ISPY2_durva"), 1, (0.15, 6), [0.2, 0.5, 1, 2, 5], False,
          "I-SPY2 durvalumab–olaparib", "pCR interaction OR")
trial_col(ax_a[3], e5("PRINCE"), 1, (0.3, 3), [0.5, 1, 2], True,
          "PRINCE (sensitivity)", "PFS interaction HR, axis reversed", color=GREY)
pooled(ax_a[4])
for ax in ax_a[:4]:
    ax.set_xlabel("per SD", fontsize=5.2, labelpad=1)
ax_a[4].set_xlabel("per SD", fontsize=5.2, labelpad=1)
ax_a[0].text(0, -0.155, "Every axis is oriented so that more benefit from the checkpoint-containing arm lies to the "
             "right; ringed points, P < 0.05.", transform=ax_a[0].transAxes, fontsize=5.0, color=GREY, ha="left", va="top")
bot = gs[1].subgridspec(1, 3, width_ratios=[1, 1, 1.1], wspace=0.58)
ax_b = fig.add_subplot(bot[0, 0]); ax_c = fig.add_subplot(bot[0, 1]); ax_d = fig.add_subplot(bot[0, 2])
within(ax_b, "ISPY2_pembro", "I-SPY2 pembrolizumab: 69 vs 179", "pembrolizumab + paclitaxel", True)
within(ax_c, "ISPY2_durva", "I-SPY2 durvalumab–olaparib: 71 vs 29", "durvalumab + olaparib + paclitaxel", False)
nulls(ax_d)
for ax, L, dx, dy in [(ax_a[0], "a", -0.62, 1.1), (ax_b, "b", -0.62, 1.07), (ax_c, "c", -0.08, 1.07),
                      (ax_d, "d", -0.08, 1.07)]:
    bs.panel(ax, L, dx=dx, dy=dy)
bs.save(fig, "Fig_E10_trials")
print("saved Fig_E10_trials")
