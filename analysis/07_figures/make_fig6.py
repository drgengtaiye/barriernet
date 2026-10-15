"""Figure 6 - fitted models against the fixed score. (a) Leave-one-cohort-out Harrell's C for
fixed scores, the three upgrades and full-transcriptome elastic nets. (b) Pembrolizumab
response AUC in PRJEB25780, every score on the log2(TPM + 1) scale."""
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
import bnstyle as bs
from bnstyle import IMMUNE, FITTED, GREY, INK, NULL, title

_orig_panel = bs.panel
bs.panel = lambda ax, l: _orig_panel(ax, l, dx=-0.62, dy=1.07)
R2 = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
COHORTS = [("TCGA", "o"), ("GSE84437", "^"), ("GSE62254", "v")]


def loco(ax):
    s5 = pd.read_csv(f"{R2}/S5_loco_cindex.csv").set_index("model")
    s11 = pd.read_csv(f"{R2}/S11_ceiling_prognosis.csv").pivot(index="model", columns="test", values="C")
    rows = [("Full transcriptome, elastic net", s11.loc["full transcriptome EN(a=0.5)"], True),
            ("(A) abundance + state", s5.loc["EN: abundance + state (30)"], True),
            ("(A) abundance only", s5.loc["EN: abundance only (10)"], True),
            ("Full transcriptome, lasso", s11.loc["full transcriptome EN(a=1.0)"], True),
            ("NetImmune, 48 genes", s5.loc["NetImmune (fixed)"], False),
            ("ESTIMATE, 282 genes", s5.loc["ESTIMATE (fixed)"], False),
            ("(B) six modules, learned weights", s5.loc["EN: 6 barriernet modules"], True),
            ("(A) state only", s5.loc["EN: state only (20)"], True),
            ("CD8A", s5.loc["CD8A (fixed)"], False)]
    net = rows[4][1][["TCGA", "GSE84437", "GSE62254"]].mean()
    top = rows[0][1][["TCGA", "GSE84437", "GSE62254"]].mean()
    ax.axvspan(net, top, color=NULL, alpha=0.6, lw=0)
    for i, (lab, r, fitted) in enumerate(rows):
        y = len(rows) - 1 - i
        for j, (c, mk) in enumerate(COHORTS):
            ax.plot(r[c], y + 0.12 + 0.12 * j, mk, color=GREY, ms=2.6, mew=0)
        m = r[["TCGA", "GSE84437", "GSE62254"]].mean()
        col = IMMUNE if lab.startswith("NetImmune") else (FITTED if fitted else INK)
        ax.plot(m, y - 0.14, "D" if fitted else "o", color=col, ms=3.4, mec="white", mew=0.4)
        ax.text(0.66, y, f"{m:.3f}", fontsize=5.2, va="center", ha="right", color=col)
    ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([r[0] for r in rows], fontsize=5.6)
    ax.set_xlim(0.48, 0.665); ax.set_xticks([0.5, 0.55, 0.6, 0.65])
    ax.set_xlabel("Harrell's C in the held-out cohort")
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.75, len(rows) - 0.45)
    ax.text((net + top) / 2, -0.55, f"gap {top - net:.3f}", fontsize=5.0, ha="center", va="center", color=INK)
    for c, mk in COHORTS:
        ax.plot([], [], mk, color=GREY, ms=2.4, label=c.replace("TCGA", "TCGA-STAD"))
    ax.plot([], [], "D", color=FITTED, ms=3, label="fitted, mean")
    ax.plot([], [], "o", color=INK, ms=3, label="fixed, mean")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=5, fontsize=5.0, columnspacing=0.8, handletextpad=0.3)
    title(ax, "Prognosis, leave-one-cohort-out")


def ici(ax):
    nul = pd.read_csv(f"{A}/E3_ici_matched_null.csv"); nul = nul[nul.scale == "log2"].set_index("score")
    ceil = pd.read_csv(f"{A}/E3_ici_ceiling.csv").sort_values("cv_auc").iloc[-1]
    peri = pd.read_csv(f"{R2}/S6_ici_univariate.csv").set_index("feature").loc["ab_Peri", "AUC"]
    sas = pd.read_csv(f"{R2}/S10_ici_with_SAS.csv"); arr = roc_auc_score(sas.Responder, sas.SAS)
    rows = [("IFN module, 5 genes", nul.loc["IFN", "AUC_obs"], None, False),
            ("Immune only, 26 genes", nul.loc["Immune", "AUC_obs"], None, False),
            ("NetImmune, 48 genes", nul.loc["NetImmune", "AUC_obs"], None, False),
            ("(A) low pericyte abundance", peri, None, False),
            ("Full transcriptome, L2 logistic", ceil.cv_auc, (ceil.lo, ceil.hi), True),
            ("(C) spatial arrangement score", arr, None, True)]
    for i, (lab, v, rng, fitted) in enumerate(rows):
        y = len(rows) - 1 - i
        col = IMMUNE if lab.startswith("NetImmune") else (FITTED if fitted else INK)
        if rng:
            ax.plot(rng, [y, y], color=col, lw=0.8)
        ax.plot(v, y, "D" if fitted else "o", color=col, ms=3.4, mec="white", mew=0.4)
        ax.text(0.99, y, f"{v:.3f}", fontsize=5.2, va="center", ha="right", color=col)
    ax.axvline(0.5, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_yticks(range(len(rows))[::-1]); ax.set_yticklabels([r[0] for r in rows], fontsize=5.6)
    ax.set_xlim(0.45, 1.0); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9]); ax.set_xlabel("AUC, pembrolizumab response")
    ax.tick_params(axis="y", length=0)
    title(ax, "Checkpoint response, 12 responders / 44")


cells = [dict(pos=(0, 1, 0, 1), letter="a", draw=loco), dict(pos=(0, 1, 1, 2), letter="b", draw=ici)]
print(bs.compose("Fig6", cells, grid=(1, 2), width_ratios=[1.15, 1], width=bs.W2, height=62 * bs.MM, wspace=0.95))
