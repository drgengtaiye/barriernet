"""Figure 8 - perturbation-based drug prioritisation and GDSC2 assessment, redrawn in the
house style from the underlying tables (the v13 panels carried old figure numbers in their
titles). Panel d replaces the v13 'feature ablation' panel, whose data could not be found:
it compares the leave-one-section-out ranking model with a leave-one-section-out global-mean
baseline (rank_baseline.py)."""
import numpy as np, pandas as pd, textwrap
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch
from scipy import stats
import bnstyle as bs
import figstyle as fs

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
S = pd.read_csv(f"{R}/results/ai_recommendations_drug_top15_clean.csv", index_col=0)
FI = pd.read_csv(f"{R}/figures/Fig7_AI_rank/Fig7_rank_feature_importance.csv")
MECH = pd.read_csv(f"{R}/validation/drug/drug_mechanism_table.csv")
RB = pd.read_csv(f"{bs.OUT}/Fig8_rank_vs_baseline.csv")
G = pd.read_csv(f"{bs.OUT}/Fig8_gdsc_points.csv")
Q = pd.read_csv("/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results/C1_gdsc_fdr.csv").set_index("Drug")
short = lambda d: d.split(" (")[0]


def draw_workflow(ax):
    ax.set_xlim(-1.5, 101.5); ax.set_ylim(0, 18); ax.set_axis_off()
    steps = ["10 Visium sections:\nprogramme and ring metrics",
             "Barrier / NetImmune-\nassociated signatures",
             "Perturbation connectivity\n(pertpy), 21 compounds",
             "Drug score per section",
             "Leave-one-section-out\nranking model (17 features)"]
    w, gap, y0, h = 17.0, 3.75, 5, 11
    for i, t in enumerate(steps):
        x = i * (w + gap)
        ax.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0.2,rounding_size=1.2",
                                    fc=bs.PAPER, ec=fs.INK, lw=0.5))
        ax.text(x + w / 2, y0 + h / 2, t, ha="center", va="center", fontsize=5.2)
        if i < len(steps) - 1:
            ax.annotate("", xy=(x + w + gap - 0.3, y0 + h / 2), xytext=(x + w + 0.9, y0 + h / 2),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color=fs.INK, mutation_scale=6))
    ax.text(3 * (w + gap) + w / 2, 1.2, "tested against Barrier score in 28 GDSC2 gastric cell lines (f, g)",
            ha="center", va="center", fontsize=5.2, color=bs.GREY)
    bs.title(ax, "Workflow", pad=2)


def draw_ranking(ax):
    m = S.mean().sort_values(ascending=False).head(12)[::-1]; sd = S.std()[m.index]
    y = np.arange(len(m))
    ax.errorbar(m.values, y, xerr=sd.values, fmt="o", color=fs.INK, ms=3.2, elinewidth=0.6, capsize=0)
    ax.set_yticks(y); ax.set_yticklabels([short(d) for d in m.index])
    ax.set_xlabel("drug score, mean ± SD over 10 sections")
    ax.set_xlim(12, 42)
    bs.title(ax, "Compound ranking")


def feat_label(f):
    f = f.replace("state_", "").replace("_minus_", "−").replace("STF3_FAP_INHBA", "FAP⁺INHBA⁺ CAF")
    f = f.replace("tls_thr_q95", "TLS peak threshold").replace("n_peaks", "TLS peak count")
    f = f.replace("z_ecm", "ECM component (z)").replace("z_tam", "TAM component (z)")
    return f.replace("_", " ")


def feat_col(f):
    if "SPP1_TAM" in f or "tam" in f:
        return bs.MACRO
    if "CAF" in f or "STF3" in f or "ecm" in f:
        return bs.STROMA
    return bs.GREY


def draw_features(ax):
    d = FI.sort_values("mean_abs_weight").tail(11)
    y = np.arange(len(d))
    ax.barh(y, d.mean_abs_weight, color=[feat_col(f) for f in d.feature], height=0.66, lw=0)
    ax.set_yticks(y); ax.set_yticklabels([feat_label(f) for f in d.feature], fontsize=5.2)
    ax.set_xlabel("mean absolute weight in ranking model")
    bs.title(ax, "Features driving the ranking")


def draw_baseline(ax):
    d = RB.sort_values("top5_model")
    x = np.arange(len(d))
    for xi, (_, r) in zip(x, d.iterrows()):
        ax.plot([xi, xi], [r.top5_model, r.top5_baseline], color=bs.LIGHT, lw=1.2, zorder=1)
    jit = 0.12
    ax.scatter(x - jit, d.top5_model, s=14, color=bs.FITTED, zorder=3, label="ranking model", edgecolor="white", lw=0.4)
    ax.scatter(x + jit, d.top5_baseline, s=14, color=bs.GREY, zorder=3, marker="s", label="mean of other sections",
               edgecolor="white", lw=0.4)
    ax.set_xticks(x); ax.set_xticklabels(d.section, rotation=60, ha="right", fontsize=5.2)
    ax.set_ylabel("top-5 overlap with observed")
    ax.legend(loc="lower right", fontsize=4.9, handletextpad=0.2, borderaxespad=0.1, bbox_to_anchor=(1.0, 0.16))
    ax.text(0.02, 0.03, f"mean top-5 overlap {RB.top5_model.mean():.2f} vs {RB.top5_baseline.mean():.2f}\n"
                        f"mean NDCG@5 {RB.model.mean():.2f} vs {RB.baseline.mean():.2f}",
            transform=ax.transAxes, fontsize=4.9, va="bottom", color=fs.INK)
    ax.set_ylim(0.0, 1.08); ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    bs.title(ax, "Ranking against a simple baseline")


def draw_mech(ax):
    ax.set_axis_off()
    top = [c for c in S.mean().sort_values(ascending=False).index]
    M = MECH.set_index("Drug")
    rows = [d for d in top if d in M.index][:13]
    n = len(rows)
    ax.set_xlim(0, 1); ax.set_ylim(n - 0.4, -1.5)
    ax.text(0.03, -1.0, "Compound", fontsize=5.4, fontweight="bold", va="center")
    ax.text(0.27, -1.0, "Pathway", fontsize=5.4, fontweight="bold", va="center")
    ax.text(0.68, -1.0, "Proposed action", fontsize=5.4, fontweight="bold", va="center")
    for i, d in enumerate(rows):
        eff = M.loc[d, "Axis_effect"]
        col = bs.STROMA if "anti-Barrier" in eff else bs.IMMUNE
        ax.plot(0.012, i, "o", color=col, ms=3.0, clip_on=False)
        ax.text(0.03, i, short(d), fontsize=5.4, va="center")
        ax.text(0.27, i, M.loc[d, "Pathway"].replace("Immunomodulation / anti-angiogenesis", "IMiD, anti-angiogenic").replace("DNA damage / topoisomerase II", "DNA damage / TOP2").replace("VEGFR/PDGFR angiogenesis", "VEGFR/PDGFR"), fontsize=5.0, va="center")
        ax.text(0.68, i, eff.replace(" / immunomodulatory", " / immunomod."), fontsize=5.0, va="center", color=col)
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color=bs.PAPER, zorder=0, lw=0)
    bs.title(ax, "Mechanisms of the top-ranked compounds", pad=1)


def draw_gdsc(drug, ax):
    d = G[G.drug == drug]
    r, p = stats.spearmanr(d.barrier_z, d.ln_ic50)
    b = np.polyfit(d.barrier_z, d.ln_ic50, 1)
    xs = np.linspace(d.barrier_z.min(), d.barrier_z.max(), 10)
    ax.plot(xs, np.polyval(b, xs), color=bs.GREY, lw=0.8)
    ax.scatter(d.barrier_z, d.ln_ic50, s=10, color=fs.INK, edgecolor="white", lw=0.4, zorder=3)
    ax.set_xlabel("Barrier score (z)"); ax.set_ylabel("ln(IC$_{50}$)")
    q = Q.loc[drug, "Barrier_q"]
    ax.set_title(f"{drug}\nSpearman r = {r:+.2f}; P = {p:.3f}, q = {q:.2f}".replace("-", "−"),
                 fontsize=5.3, loc="left", pad=3, color=fs.INK)


fig = plt.figure(figsize=(bs.W2, 168 * bs.MM))
gs = GridSpec(2, 1, figure=fig, height_ratios=[1.25, 1.05], hspace=0.40)
top = gs[0].subgridspec(2, 1, height_ratios=[0.22, 1.0], hspace=0.38)
r2 = top[1].subgridspec(1, 3, width_ratios=[1, 1.05, 1.0], wspace=0.95)
r3 = gs[1].subgridspec(1, 3, width_ratios=[2.25, 1, 1], wspace=0.42)
ax_a = fig.add_subplot(top[0])
ax_b = fig.add_subplot(r2[0]); ax_c = fig.add_subplot(r2[1]); ax_d = fig.add_subplot(r2[2])
ax_e = fig.add_subplot(r3[0]); ax_f = fig.add_subplot(r3[1]); ax_g = fig.add_subplot(r3[2])
draw_workflow(ax_a); draw_ranking(ax_b); draw_features(ax_c); draw_baseline(ax_d); draw_mech(ax_e)
draw_gdsc("Dasatinib", ax_f); draw_gdsc("Sorafenib", ax_g)
for ax, L, dx, dy in [(ax_a, "a", -0.035, 1.24), (ax_b, "b", -0.5, 1.1), (ax_c, "c", -0.95, 1.1),
                      (ax_d, "d", -0.3, 1.1), (ax_e, "e", -0.045, 1.13), (ax_f, "f", -0.32, 1.17),
                      (ax_g, "g", -0.32, 1.17)]:
    bs.panel(ax, L, dx=dx, dy=dy)
bs.save(fig, "Fig8")
print("saved Fig8")
