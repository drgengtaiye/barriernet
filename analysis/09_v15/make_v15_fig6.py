"""v15 Figure 6 - pembrolizumab response in PRJEB25780 (44 patients, 12 responders), all scores on
log2(TPM + 1): ROC curves, the structure-matched random-axis null, AUC with CI for every marker
(DeLong vs NetImmune, BH), the two Barrier components by MSI/EBV status, and cross-validated AUC
of fixed scores against a fitted full-transcriptome classifier through the same pipeline."""
import numpy as np, pandas as pd
import matplotlib.patches as mpatches
from sklearn.metrics import roc_curve, roc_auc_score
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, FITTED, GREY, INK, NULL, title, fmt_p

_orig_panel = bs.panel
bs.panel = lambda ax, l: _orig_panel(ax, l, dx=-0.30, dy=1.10)
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
S = pd.read_csv(f"{A}/E7_scores_log2.csv")
y = S.Responder.values.astype(int)


def roc(ax):
    curves = [("NetImmune", S.NetImmune, IMMUNE, "-", 1.2), ("IFN module", S.IFN, INK, "-", 0.8),
              ("GEP-18", S.GEP18, INK, (0, (3, 1.5)), 0.8), ("PD-L1 (CD274)", S.PDL1, GREY, (0, (1, 1.2)), 0.9),
              ("IPS", S.IPS_AZ, GREY, (0, (4, 1.5, 1, 1.5)), 0.8)]
    for lab, s, col, ls, lw in curves:
        f, t, _ = roc_curve(y, s)
        ax.plot(f, t, color=col, ls=ls, lw=lw, label=f"{lab}  {roc_auc_score(y, s):.2f}", drawstyle="steps-post")
    ax.plot([0, 1], [0, 1], color=NULL, lw=0.5)
    ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.01); ax.set_aspect("equal")
    ax.set_xlabel("1 − specificity"); ax.set_ylabel("Sensitivity")
    ax.legend(loc="lower right", fontsize=5.0, title="AUC", title_fontsize=5.0, handlelength=1.8)
    title(ax, "ROC, 12 responders / 44")


def null(ax):
    d = pd.read_csv(f"{A}/E3_ici_null_draws_log2.csv")
    s = pd.read_csv(f"{A}/E3_ici_matched_null.csv"); s = s[s.scale == "log2"].set_index("score")
    items = [("NetImmune", "NetImmune (48 genes)"), ("Immune", "Immune only (26 genes)"), ("IFN", "IFN module (5 genes)")]
    bins = np.linspace(0.5, 0.95, 46)
    for i, (k, lab) in enumerate(items):
        yb = len(items) - 1 - i
        h, e = np.histogram(d[k], bins=bins); h = h / h.max() * 0.48
        ax.fill_between(e[:-1], yb, yb + h, step="post", color=NULL, lw=0)
        q95 = s.loc[k, "null_p95"]
        ax.plot([q95, q95], [yb, yb + 0.52], color=GREY, lw=0.6)
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
    ax.legend(loc="upper center", bbox_to_anchor=(0.4, -0.26), ncol=3, fontsize=5.0, columnspacing=0.5, handletextpad=0.3)
    title(ax, "Structure-matched random-gene null")


def forest(ax):
    a = pd.read_csv(f"{A}/E7_auc_ci_log2.csv")
    dl = pd.read_csv(f"{A}/E7_delong_log2.csv").set_index("Comparator")
    names = {"GEP18": "GEP-18", "IFN": "IFN module", "APM": "APM module", "TAM_SUP": "TAM_SUP module*",
             "T_CELL": "T_CELL module", "TLS_B": "TLS_B module", "Immune": "Immune only",
             "Barrier(neg)": "Barrier (inverted)", "MSI/EBV": "MSI/EBV status", "MSI-H": "MSI-H status",
             "PD-L1": "PD-L1 (CD274)", "IPS": "IPS",
             "ExclusionProxy(neg)": "Exclusion proxy†"}
    a = a.sort_values("AUC")
    for i, r in enumerate(a.itertuples()):
        col = {"NetImmune": IMMUNE, "Barrier(neg)": STROMA, "TAM_SUP": MACRO}.get(r.Marker, INK)
        ax.plot([r.CI_low, r.CI_high], [i, i], color=col, lw=0.8)
        ax.plot(r.AUC, i, "o", color=col, ms=3.0, mec="white", mew=0.4)
        if r.Marker in dl.index and dl.loc[r.Marker, "q"] < 0.05:
            q = dl.loc[r.Marker, "q"]
            ax.text(1.02, i, "q < 0.001" if q < 0.001 else f"q = {q:.3f}", fontsize=5.0, va="center", color=INK)
    ax.axvline(0.5, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_yticks(range(len(a))); ax.set_yticklabels([names.get(m, m) for m in a.Marker], fontsize=5.4)
    ax.set_xlim(0.1, 1.16); ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_xlabel("AUC (95% CI)", labelpad=1)
    ax.tick_params(axis="y", length=0)
    title(ax, "All markers vs NetImmune (DeLong)")
    ax.text(0.0, -0.085, "* higher TAM_SUP predicts response\n† CAF/TAM exclusion proxy, 16 genes, inverted",
            transform=ax.transAxes, fontsize=5.0, color=GREY, va="top")
    ax.text(-0.075, 1.047, "c", transform=ax.transAxes, fontsize=8.2, fontweight="bold", va="top", color=INK)


def components(ax):
    rng = np.random.default_rng(1)
    hot = (S.MSI_H == 1) | (S.EBV == 1)
    for mod, col, x0 in [("CAF_ECM", STROMA, 0), ("TAM_SUP", MACRO, 2.6)]:
        for j, g in enumerate([0, 1]):
            sel = y == g
            v = S.loc[sel, mod].values
            x = x0 + j
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.add_patch(mpatches.Rectangle((x - 0.28, q1), 0.56, q3 - q1, fc="none", ec=col, lw=0.7))
            ax.plot([x - 0.28, x + 0.28], [med, med], color=col, lw=1.0)
            for hsel, mk in [(~hot, "o"), (hot, "^")]:
                vv = S.loc[sel & hsel, mod].values
                ax.plot(x + rng.uniform(-0.18, 0.18, len(vv)), vv, mk, color=col, ms=2.3, alpha=0.8, mew=0)
        auc = roc_auc_score(y, S[mod])
        txt = f"AUC {1 - auc:.2f} (low)" if auc < 0.5 else f"AUC {auc:.2f} (high)"
        ax.text(x0 + 0.5, 2.55, f"{mod}\n{txt}", ha="center", va="bottom", fontsize=5.0, color=INK)
    ax.set_xticks([0, 1, 2.6, 3.6]); ax.set_xticklabels(["non-resp.", "resp.", "non-resp.", "resp."], fontsize=5.4)
    ax.set_xlim(-0.6, 4.2); ax.set_ylim(-2.9, 3.5); ax.set_ylabel("Module score (z, log₂ scale)")
    ax.axhline(0, color=NULL, lw=0.5)
    ax.plot([], [], "o", color=GREY, ms=2.5, label="MSS, EBV−"); ax.plot([], [], "^", color=GREY, ms=2.5, label="MSI-H or EBV+")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=5.0, handletextpad=0.3, columnspacing=0.8)
    title(ax, "The two halves of Barrier")


def samecv(ax):
    cv = pd.read_csv(f"{A}/E7_samepipeline_cv.csv")
    best = cv.loc[cv.groupby("score").cv_auc.idxmax()].set_index("score")
    ce = pd.read_csv(f"{A}/E3_ici_ceiling.csv").sort_values("cv_auc").iloc[-1]
    rows = [("IFN module", best.loc["IFN"], False), ("Immune only", best.loc["Immune"], False),
            ("NetImmune", best.loc["NetImmune"], False), ("Full transcriptome\n(fitted)", ce, True)]
    for i, (lab, r, fitted) in enumerate(rows):
        yy = len(rows) - 1 - i
        col = FITTED if fitted else (IMMUNE if lab.startswith("NetImmune") else INK)
        ax.plot([r.lo, r.hi], [yy, yy], color=col, lw=0.9)
        ax.plot(r.cv_auc, yy, "D" if fitted else "o", color=col, ms=3.4, mec="white", mew=0.4)
        ax.text(0.905, yy, f"{r.cv_auc:.3f}", fontsize=5.2, va="center", ha="right", color=col)
    ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([r[0] for r in rows], fontsize=5.4)
    ax.set_xlim(0.62, 0.91); ax.set_xticks([0.65, 0.7, 0.75, 0.8, 0.85])
    ax.set_xlabel("Cross-validated AUC (40 × 5-fold)")
    ax.tick_params(axis="y", length=0)
    title(ax, "Same CV pipeline")


cells = [dict(pos=(0, 1, 0, 2), letter="a", draw=roc),
         dict(pos=(0, 1, 2, 4), letter="b", draw=null),
         dict(pos=(0, 2, 4, 6), letter="", draw=forest),
         dict(pos=(1, 2, 0, 2), letter="d", draw=components),
         dict(pos=(1, 2, 2, 4), letter="e", draw=samecv)]
print(bs.compose("v15_Fig6", cells, grid=(2, 6), height_ratios=[1, 0.95], width=bs.W2,
                 height=124 * bs.MM, hspace=0.62, wspace=1.8))
