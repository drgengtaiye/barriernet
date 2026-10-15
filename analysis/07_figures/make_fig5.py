"""Figure 5 - pembrolizumab response in PRJEB25780 (n = 44, 12 responders), published
linear-TPM scores: ROC curves, the structure-matched random-axis null, AUC with CI for all
markers, and the two Barrier components in responders and non-responders."""
import numpy as np, pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, GREY, INK, NULL, title, fmt_p

_orig_panel = bs.panel
bs.panel = lambda ax, l: _orig_panel(ax, l, dx=-0.30, dy=1.10)
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
R = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
b = pd.read_csv(f"{V}/benchmark_scores.tsv", sep="\t", index_col=0)
y = b.Responder.values.astype(int)


def roc(ax):
    curves = [("NetImmune", b.NetImmune, IMMUNE, "-", 1.2), ("IFN module", b.IFN, INK, "-", 0.8),
              ("GEP-18", b.GEP18, INK, (0, (3, 1.5)), 0.8), ("TIDE (inverted)", -b.TIDE, GREY, (0, (1, 1.2)), 0.9),
              ("IPS", b.IPS, GREY, (0, (4, 1.5, 1, 1.5)), 0.8)]
    for lab, s, col, ls, lw in curves:
        f, t, _ = roc_curve(y, s)
        ax.plot(f, t, color=col, ls=ls, lw=lw, label=f"{lab}  {roc_auc_score(y, s):.2f}", drawstyle="steps-post")
    ax.plot([0, 1], [0, 1], color=NULL, lw=0.5)
    ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.01); ax.set_aspect("equal")
    ax.set_xlabel("1 − specificity"); ax.set_ylabel("Sensitivity")
    ax.legend(loc="lower right", fontsize=4.9, title="AUC", title_fontsize=5.0, handlelength=1.8)
    title(ax, "ROC, 12 responders / 44")


def null(ax):
    d = pd.read_csv(f"{A}/E3_ici_null_draws_linear.csv")
    s = pd.read_csv(f"{A}/E3_ici_matched_null.csv"); s = s[s.scale == "linear"].set_index("score")
    items = [("NetImmune", "NetImmune (48 genes)"), ("Immune", "Immune only (26 genes)"), ("IFN", "IFN module (5 genes)")]
    bins = np.linspace(0.5, 0.95, 46)
    for i, (k, lab) in enumerate(items):
        yb = len(items) - 1 - i
        h, e = np.histogram(d[k], bins=bins); h = h / h.max() * 0.48
        ax.fill_between(e[:-1], yb, yb + h, step="post", color=NULL, lw=0)
        q95 = s.loc[k, "null_p95"]
        ax.plot([q95, q95], [yb, yb + 0.55], color=GREY, lw=0.6)
        ax.plot(s.loc[k, "AUC_obs"], yb + 0.18, "D", color=IMMUNE, ms=3.2, mec="white", mew=0.4)
        pv = s.loc[k, "p_emp"]
        ptxt = "P ≤ 0.0005" if pv <= 0.0005 else fmt_p(pv)
        ax.text(0.505, yb + 0.86, lab, fontsize=5.6, fontweight="bold", va="center")
        ax.text(0.945, yb + 0.66, f"AUC {s.loc[k, 'AUC_obs']:.3f}, {ptxt}", fontsize=5.0, ha="right", va="center")
    ax.set_ylim(-0.1, len(items)); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0.5, 0.95); ax.set_xlabel("AUC (random axes: better direction)")
    ax.plot([], [], "D", color=IMMUNE, ms=3, label="observed")
    ax.plot([], [], "s", color=NULL, ms=4, label="random axes")
    ax.plot([], [], "|", color=GREY, ms=5, label="95th pct.")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=3, fontsize=4.9, columnspacing=0.6, handletextpad=0.3)
    title(ax, "Structure-matched random-gene null")


def forest(ax):
    a = pd.read_csv(f"{R}/A2_auc_ci.csv"); dl = pd.read_csv(f"{R}/A2_delong.csv").set_index("Comparator")
    names = {"GEP18": "GEP-18", "TIDE(neg)": "TIDE (inverted)", "IFN": "IFN module", "APM": "APM module",
             "TAM_SUP": "TAM_SUP module*", "T_CELL": "T_CELL module", "TLS_B": "TLS_B module",
             "Immune": "Immune only", "Barrier": "Barrier (inverted)", "MSI/EBV": "MSI/EBV status",
             "MSI-H": "MSI-H status"}
    a = a.sort_values("AUC")
    for i, r in enumerate(a.itertuples()):
        col = {"NetImmune": IMMUNE, "Barrier": STROMA, "TAM_SUP": MACRO}.get(r.Marker, INK)
        ax.plot([r.CI_low, r.CI_high], [i, i], color=col, lw=0.8)
        ax.plot(r.AUC, i, "o", color=col, ms=3.0, mec="white", mew=0.4)
        if r.Marker in dl.index and dl.loc[r.Marker, "p_FDR"] < 0.05:
            ax.text(1.02, i, f"q = {dl.loc[r.Marker, 'p_FDR']:.3f}", fontsize=4.8, va="center", color=INK)
    ax.axvline(0.5, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_yticks(range(len(a))); ax.set_yticklabels([names.get(m, m) for m in a.Marker], fontsize=5.4)
    ax.set_xlim(0.25, 1.12); ax.set_xticks([0.3, 0.5, 0.7, 0.9]); ax.set_xlabel("AUC (95% CI)")
    ax.tick_params(axis="y", length=0)
    title(ax, "All markers vs NetImmune (DeLong)")
    ax.text(-0.075, 1.047, "c", transform=ax.transAxes, fontsize=8.2, fontweight="bold", va="top", color=INK)


def components(ax):
    rng = np.random.default_rng(1)
    groups = [("CAF_ECM", STROMA, 0), ("TAM_SUP", MACRO, 2.6)]
    for mod, col, x0 in groups:
        for j, (lab, g) in enumerate([("NR", 0), ("R", 1)]):
            v = b.loc[y == g, mod].values
            x = x0 + j
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.add_patch(__import__("matplotlib").patches.Rectangle((x - 0.28, q1), 0.56, q3 - q1, fc="none", ec=col, lw=0.7))
            ax.plot([x - 0.28, x + 0.28], [med, med], color=col, lw=1.0)
            ax.plot(x + rng.uniform(-0.18, 0.18, len(v)), v, "o", color=col, ms=2.0, alpha=0.75, mew=0)
        auc = roc_auc_score(y, b[mod])
        txt = f"AUC {1 - auc:.2f} for low scores" if auc < 0.5 else f"AUC {auc:.2f} for high scores"
        ax.text(x0 + 0.5, 2.55, f"{mod}\n{txt}", ha="center", va="bottom", fontsize=5.2, color=INK)
    ax.set_xticks([0, 1, 2.6, 3.6]); ax.set_xticklabels(["non-resp.", "resp.", "non-resp.", "resp."], fontsize=5.4)
    ax.set_xlim(-0.6, 4.2); ax.set_ylim(-2.3, 3.4); ax.set_ylabel("Module score (z)")
    ax.axhline(0, color=NULL, lw=0.5)
    title(ax, "The two halves of Barrier")


cells = [dict(pos=(0, 1, 0, 2), letter="a", draw=roc),
         dict(pos=(0, 1, 2, 4), letter="b", draw=null),
         dict(pos=(0, 2, 4, 6), letter="", draw=forest),
         dict(pos=(1, 2, 0, 3), letter="d", draw=components)]
print(bs.compose("Fig5", cells, grid=(2, 6), height_ratios=[1, 0.95], width=bs.W2,
                 height=120 * bs.MM, hspace=0.62, wspace=1.7))
