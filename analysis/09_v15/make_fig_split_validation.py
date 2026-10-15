"""Fig_split_validation - the two halves of Barrier outside the discovery data.
(a) IMmotion150: within-arm HR per SD (atezolizumab +/- bevacizumab vs sunitinib) and the
    score x arm interaction. (b) 11 checkpoint cohorts in other cancers: log-odds of response
    per SD, per cohort and pooled. (c) Size-matched random-gene nulls for the pooled effects.
(d) The TAM_SUP arms and the CXCL9:SPP1 polarity score in every setting."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, GREY, INK, NULL, fmt_p

A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
T = pd.read_csv(f"{A}/E5_trials.csv"); T = T[T.trial.str.startswith("IMmotion150: atezolizumab (+/- bev)")].set_index("score")
per = pd.read_csv(f"{A}/E5_icb_per_cohort.csv")
pool = pd.read_csv(f"{A}/E5_icb_pooled.csv"); pool = pool[pool.set == "all 11"].set_index("score")
nd = pd.read_csv(f"{A}/E5_icb_null_draws.csv"); ns = pd.read_csv(f"{A}/E5_icb_null_summary.csv").set_index("score")
G = pd.read_csv(f"{A}/E5_gastric_arms_pooled.csv").set_index("score")
K = pd.read_csv(f"{A}/E5_prjeb_arms.csv"); K = K[K.scale == "linear"].set_index("score")

ROWS = [("CAF_ECM", "CAF_ECM", STROMA), ("TAM_SUP", "TAM_SUP", MACRO), ("C1Q_arm", "  C1Q arm", MACRO),
        ("SPP1_arm", "  SPP1 arm", MACRO), ("Barrier", "Barrier", STROMA), ("CXCL9_SPP1", "CXCL9:SPP1", INK),
        ("IFN", "IFN", IMMUNE), ("Immune", "Immune only", IMMUNE), ("NetImmune", "NetImmune", IMMUNE)]


def rowaxis(ax, rows=ROWS):
    ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([r[1] for r in rows], fontsize=5.6)
    ax.tick_params(axis="y", length=0)
    for lab, (k, _, c) in zip(ax.get_yticklabels(), rows):
        lab.set_color(c if c != INK else INK)


def trial(ax):
    for i, (k, _, c) in enumerate(ROWS):
        y = len(ROWS) - 1 - i; r = T.loc[k]
        ax.plot([r.ctrl_lo, r.ctrl_hi], [y - 0.16] * 2, color=GREY, lw=0.7)
        ax.plot(r.ctrl_HR, y - 0.16, "s", mfc="white", mec=GREY, ms=2.8, mew=0.7)
        ax.plot([r.ici_lo, r.ici_hi], [y + 0.16] * 2, color=c, lw=0.8)
        ax.plot(r.ici_HR, y + 0.16, "o", color=c, ms=3.0, mec="white", mew=0.4)
        tr = ax.get_yaxis_transform()
        ax.text(1.17, y, f"{r.int_HR:.2f}", transform=tr, clip_on=False, fontsize=5.1, ha="right", va="center",
                color=INK, fontweight="bold" if r.int_p < 0.05 else "normal")
        ax.text(1.33, y, fmt_p(r.int_p, lead=""), transform=tr, clip_on=False, fontsize=5.0, ha="right", va="center", color=INK)
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xscale("log"); ax.set_xlim(0.5, 1.9); ax.set_xticks([0.5, 0.7, 1, 1.4])
    ax.set_xticklabels(["0.5", "0.7", "1", "1.4"]); ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel("PFS HR per SD, within arm")
    tr = ax.get_yaxis_transform()
    ax.text(1.17, -0.95, "interaction\nHR", transform=tr, clip_on=False, fontsize=4.9, ha="right", va="top", color=GREY)
    ax.text(1.33, -0.95, "P", transform=tr, clip_on=False, fontsize=4.9, ha="right", va="top", color=GREY)
    ax.set_ylim(-0.6, len(ROWS) - 0.4); rowaxis(ax)
    ax.plot([], [], "o", color=INK, ms=3, label="atezolizumab ± bev (n = 174)")
    ax.plot([], [], "s", mfc="white", mec=GREY, ms=2.8, label="sunitinib (n = 89)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.42, -0.17), ncol=2, fontsize=5.0, columnspacing=0.8)
    bs.title(ax, "IMmotion150, randomised (renal)")


def icb(ax):
    for i, (k, _, c) in enumerate(ROWS):
        y = len(ROWS) - 1 - i; p = per[per.score == k]; q = pool.loc[k]
        ax.plot(p.logOR, y + 0.2 + np.linspace(-0.05, 0.05, len(p)), ".", color=GREY, ms=2.6, mew=0)
        ax.plot([q.lo, q.hi], [y - 0.12] * 2, color=c, lw=0.9)
        ax.plot(q.logOR, y - 0.12, "D", color=c, ms=3.0, mec="white", mew=0.4)
        ax.text(1.14, y, f"{int(q.n_pos)}/11", transform=ax.get_yaxis_transform(), clip_on=False,
                fontsize=5.0, ha="right", va="center", color=INK)
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(-1.3, 1.65); ax.set_xticks([-1, -0.5, 0, 0.5, 1])
    ax.set_xlabel("log-odds of response per SD")
    ax.text(1.14, -0.95, "cohorts\n> 0", transform=ax.get_yaxis_transform(), clip_on=False, fontsize=4.9,
            ha="right", va="top", color=GREY)
    ax.set_ylim(-0.6, len(ROWS) - 0.4); rowaxis(ax)
    ax.plot([], [], ".", color=GREY, ms=4, label="per cohort")
    ax.plot([], [], "D", color=INK, ms=3, label="pooled, 95% CI")
    ax.legend(loc="upper center", bbox_to_anchor=(0.42, -0.17), ncol=2, fontsize=5.0, columnspacing=0.8)
    bs.title(ax, "11 checkpoint cohorts, other cancers")


def nulls(ax):
    items = [("CAF_ECM", "size11", STROMA, "CAF_ECM (11 genes)"), ("TAM_SUP", "size11", MACRO, "TAM_SUP (11 genes)"),
             ("C1Q_arm", "size3", MACRO, "C1Q arm (3 genes)"), ("SPP1_arm", "size5", MACRO, "SPP1 arm (5 genes)")]
    bins = np.linspace(0, 0.4, 41)
    for i, (k, col, c, lab) in enumerate(items):
        yb = len(items) - 1 - i
        h, e = np.histogram(np.abs(nd[col]), bins=bins); h = h / h.max() * 0.5
        ax.fill_between(e[:-1], yb, yb + h, step="post", color=NULL, lw=0)
        ax.plot([ns.loc[k, "null_p95_abs"]] * 2, [yb, yb + 0.55], color=GREY, lw=0.6)
        ax.plot(abs(ns.loc[k, "obs"]), yb + 0.18, "D", color=c, ms=3.2, mec="white", mew=0.4)
        pv = ns.loc[k, "p_emp"]
        ax.text(0.005, yb + 0.84, lab, fontsize=5.5, fontweight="bold", va="center", color=c)
        ax.text(0.398, yb + 0.84, "P ≤ 0.001" if pv <= 0.001 else fmt_p(pv), fontsize=5.0, ha="right", va="center")
    ax.set_ylim(-0.1, len(items)); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 0.4); ax.set_xlabel("absolute pooled log-OR, random gene sets of the same size")
    bs.title(ax, "Size-matched random-gene null (1,000 sets)")


def arms(axs):
    rows = [("TAM_SUP", "TAM_SUP"), ("C1Q_arm", "C1Q arm"), ("SPP1_arm", "SPP1 arm"), ("CXCL9_SPP1", "CXCL9:SPP1")]
    cols = {"TAM_SUP": MACRO, "C1Q_arm": MACRO, "SPP1_arm": MACRO, "CXCL9_SPP1": INK}
    mk = {"TAM_SUP": "o", "C1Q_arm": "^", "SPP1_arm": "s", "CXCL9_SPP1": "D"}
    spec = [("Gastric prognosis, HR", lambda k: (G.loc[k, "HR"], G.loc[k, "lo"], G.loc[k, "hi"]), 1, (0.7, 1.4), True),
            ("PRJEB25780 response, AUC", lambda k: (K.loc[k, "AUC"], K.loc[k, "lo"], K.loc[k, "hi"]), 0.5, (0.3, 1.0), False),
            ("11 other cohorts, log-OR", lambda k: (pool.loc[k, "logOR"], pool.loc[k, "lo"], pool.loc[k, "hi"]), 0, (-0.3, 0.8), False),
            ("IMmotion150, interaction HR", lambda k: (T.loc[k, "int_HR"], T.loc[k, "int_lo"], T.loc[k, "int_hi"]), 1, (0.45, 1.6), True)]
    for j, (ax, (ttl, get, ref, xl, lg)) in enumerate(zip(axs, spec)):
        for i, (k, lab) in enumerate(rows):
            y = len(rows) - 1 - i; v, lo, hi = get(k)
            ax.plot([lo, hi], [y, y], color=cols[k], lw=0.8)
            ax.plot(v, y, mk[k], color=cols[k], ms=3.0, mec="white", mew=0.4)
        ax.axvline(ref, color=GREY, lw=0.5, ls=(0, (2, 2)))
        if lg:
            ax.set_xscale("log"); ax.xaxis.set_minor_formatter(plt.NullFormatter())
            ax.set_xticks([0.5, 0.7, 1, 1.4]); ax.set_xticklabels(["0.5", "0.7", "1", "1.4"])
        ax.set_xlim(*xl); ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_title(ttl, fontsize=5.4, loc="left", pad=3, color=INK)
        if j % 2 == 0:
            ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([r[1] for r in rows], fontsize=5.6)
        else:
            ax.set_yticks([])
        ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", labelsize=5.2)


fig = plt.figure(figsize=(bs.W2, 160 * bs.MM))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.25, 0.9], width_ratios=[1, 1], hspace=0.55, wspace=0.98)
ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[0, 1]); ax_c = fig.add_subplot(gs[1, 0])
sub = gs[1, 1].subgridspec(2, 2, wspace=0.12, hspace=0.75)
ax_d = [fig.add_subplot(sub[i // 2, i % 2]) for i in range(4)]
trial(ax_a); icb(ax_b); nulls(ax_c); arms(ax_d)
for ax, L, dx, dy in [(ax_a, "a", -0.3, 1.07), (ax_b, "b", -0.3, 1.07), (ax_c, "c", -0.06, 1.12), (ax_d[0], "d", -0.95, 1.3)]:
    bs.panel(ax, L, dx=dx, dy=dy)
bs.save(fig, "Fig_split_validation")
print("saved")
