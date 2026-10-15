"""E16 figure - PERFECT oesophageal checkpoint cohort and METABRIC bulk prognosis (Fig_E16_perfect_metabric).
(a) PERFECT (GSE165252; neoadjuvant atezolizumab + chemoradiotherapy): log-odds of response per SD,
    primary scale (log2 of the deposited 'vst' file) and log2(norm.cnt + 1) sensitivity.
(b) Pooled log-odds per SD across all checkpoint cohorts including chemo-ICB, as in E12 and with PERFECT
    added (wider of the Hartung-Knapp and DerSimonian-Laird intervals).
(c) METABRIC: HR per SD for disease-specific survival, adjusted for age, NPI and ER (all), and within ER+
    and ER- tumours (adjusted for age and NPI).
(d) METABRIC against the TCGA-BRCA stand-in used in E11 (adjusted log HR per SD, OS and DSS)."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, GREY, INK, NULL

A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
P = pd.read_csv(f"{A}/E16_perfect_per_score.csv")
Q = pd.read_csv(f"{A}/E16_perfect_pooled.csv")
C = pd.read_csv(f"{A}/E16_metabric_cox.csv")
V = pd.read_csv(f"{A}/E16_metabric_vs_tcga_brca.csv")

ROWS = [("CAF_ECM", "CAF_ECM", STROMA), ("Barrier", "Barrier", STROMA),
        ("TAM_SUP", "TAM_SUP", MACRO), ("C1Q_arm", "C1Q arm", MACRO), ("SPP1_arm", "SPP1 arm", MACRO),
        ("CXCL9_SPP1", "CXCL9:SPP1", MACRO),
        ("IFN", "IFN", IMMUNE), ("APM", "APM", IMMUNE), ("TLS_B", "TLS_B", IMMUNE), ("T_CELL", "T_CELL", IMMUNE),
        ("Immune", "Immune", IMMUNE), ("NetImmune", "NetImmune", IMMUNE)]
NR = len(ROWS); FS = 5.2
ypos = lambda i: NR - 1 - i


def rowaxis(ax, labels=True):
    ax.set_ylim(-0.7, NR - 0.3)
    if labels:
        ax.set_yticks([ypos(i) for i in range(NR)]); ax.set_yticklabels([r[1] for r in ROWS], fontsize=5.6)
        for lab, (_, _, c) in zip(ax.get_yticklabels(), ROWS):
            lab.set_color(c)
    else:
        ax.set_yticks([])
    ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", labelsize=5.2)
    for i in (2, 6):
        ax.axhline(ypos(i) + 0.5, color=NULL, lw=0.5)


def point(ax, v, lo, hi, y, c, mk="o", ms=3.0, open_=False, lw=0.8):
    ax.plot([lo, hi], [y, y], color=c, lw=lw, solid_capstyle="butt")
    ax.plot(v, y, mk, ms=ms, color=c, mfc="white" if open_ else c, mec=c if open_ else "white", mew=0.7 if open_ else 0.4)


def perfect(ax):
    for i, (k, _, c) in enumerate(ROWS):
        a = P[(P.scale == "primary_log2_vst") & (P.score == k)].iloc[0]
        b = P[(P.scale == "sens_log2_normcnt") & (P.score == k)].iloc[0]
        point(ax, a.logOR, a.lo, a.hi, ypos(i) + 0.14, c)
        point(ax, b.logOR, b.lo, b.hi, ypos(i) - 0.16, c, open_=True, ms=2.6, lw=0.6)
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(-1.6, 1.8); rowaxis(ax, True)
    ax.set_xlabel("log-odds of response per SD", fontsize=5.4)
    n = P.n.iloc[0]; r = P.n_resp.iloc[0]
    ax.set_title(f"PERFECT, oesophageal adenocarcinoma: {r}/{n} responders", fontsize=5.6, loc="left", pad=9,
                 color=INK, fontweight="bold")
    ax.text(0, 1.012, "atezolizumab + chemoradiotherapy (pre-treatment biopsies)", transform=ax.transAxes,
            fontsize=4.9, color=GREY, va="bottom")
    ax.plot([], [], "o", color=INK, ms=3, label=r"variance-stabilised log$_\mathrm{2}$ expression (primary)")
    ax.plot([], [], "o", mfc="white", mec=INK, ms=2.6, label=r"log$_\mathrm{2}$(normalised count + 1)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.13), ncol=2, fontsize=4.9, frameon=False,
              columnspacing=0.8, handletextpad=0.3)


def pooled(ax):
    before = Q[Q.set == "all + chemo-ICB (E12)"].set_index("score")
    after = Q[Q.set == "all + chemo-ICB + PERFECT"].set_index("score")
    tr = ax.get_yaxis_transform()
    for i, (k, _, c) in enumerate(ROWS):
        b, a = before.loc[k], after.loc[k]
        point(ax, b.logOR, b.lo, b.hi, ypos(i) + 0.14, GREY, open_=True, ms=2.6, lw=0.6)
        point(ax, a.logOR, a.lo, a.hi, ypos(i) - 0.16, c, mk="D", ms=2.8)
        ax.text(1.02, ypos(i), f"{a.logOR:+.2f}", transform=tr, fontsize=FS, ha="left", va="center", color=INK)
        ax.text(1.27, ypos(i), bs.fmt_p(a.p, lead=""), transform=tr, fontsize=FS, ha="left", va="center",
                color=INK, fontweight="bold" if a.p < 0.05 else "normal")
    ax.text(1.02, 1.012, "log-OR", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom")
    ax.text(1.27, 1.012, "P", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom")
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(-0.35, 0.75); rowaxis(ax, False)
    ax.set_xlabel("pooled log-odds of response per SD", fontsize=5.4)
    k_b = int(before.loc["Immune", "k"]); k_a = int(after.loc["Immune", "k"])
    ax.set_title("All checkpoint cohorts including chemotherapy–checkpoint cohorts", fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.text(0, 1.012, "wider of Hartung–Knapp and DerSimonian–Laird 95% CI; larger P", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom")
    ax.plot([], [], "o", mfc="white", mec=GREY, ms=2.6, label=f"before PERFECT ({k_b} cohorts)")
    ax.plot([], [], "D", color=INK, ms=2.8, label=f"+ PERFECT ({k_a})")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.13), ncol=2, fontsize=4.9, frameon=False,
              columnspacing=0.8, handletextpad=0.3)


def metabric(ax):
    d = C[C.endpoint == "DSS"]
    for i, (k, _, c) in enumerate(ROWS):
        a = d[(d.score == k) & (d.stratum == "all") & (d.model == "full")].iloc[0]
        p = d[(d.score == k) & (d.stratum == "ERpos")].iloc[0]
        n = d[(d.score == k) & (d.stratum == "ERneg")].iloc[0]
        point(ax, a.HR, a.lo, a.hi, ypos(i) + 0.22, c)
        point(ax, p.HR, p.lo, p.hi, ypos(i), c, mk="s", open_=True, ms=2.5, lw=0.55)
        point(ax, n.HR, n.lo, n.hi, ypos(i) - 0.22, c, mk="^", open_=True, ms=2.7, lw=0.55)
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xscale("log"); ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xticks([0.6, 0.8, 1, 1.25]); ax.set_xticklabels(["0.6", "0.8", "1", "1.25"]); ax.set_xlim(0.58, 1.38)
    rowaxis(ax, True)
    ax.set_xlabel("HR per SD, disease-specific survival", fontsize=5.4)
    na = int(d[(d.stratum == "all") & (d.model == "full")].n.iloc[0]); ea = int(d[(d.stratum == "all") & (d.model == "full")].events.iloc[0])
    ax.set_title(f"METABRIC: {na:,} patients, {ea} breast-cancer deaths", fontsize=5.6, loc="left", pad=9, color=INK,
                 fontweight="bold")
    ax.text(0, 1.012, "Illumina HT-12; Cox, adjusted", transform=ax.transAxes, fontsize=4.9, color=GREY, va="bottom")
    ax.plot([], [], "o", color=INK, ms=3, label="all (+ age, NPI, ER)")
    ax.plot([], [], "s", mfc="white", mec=INK, ms=2.5, label="ER+ (+ age, NPI)")
    ax.plot([], [], "^", mfc="white", mec=INK, ms=2.7, label="ER− (+ age, NPI)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.42, -0.13), ncol=3, fontsize=4.9, frameon=False,
              columnspacing=0.6, handletextpad=0.3)


def versus(ax):
    col = {k: c for k, _, c in ROWS}
    lim = (-0.4, 0.35)
    ax.plot(lim, lim, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.axhline(0, color=NULL, lw=0.5); ax.axvline(0, color=NULL, lw=0.5)
    for _, r in V.iterrows():
        ax.plot(np.log(r.HR_TCGA_BRCA), np.log(r.HR_METABRIC), "o" if r.endpoint == "OS" else "^",
                ms=3.2, color=col[r.score], mec="white", mew=0.4)
    for k, lab in [("NetImmune", "NetImmune"), ("CAF_ECM", "CAF_ECM")]:
        r = V[(V.score == k) & (V.endpoint == "DSS")].iloc[0]
        off = {"NetImmune": (-2, -11), "CAF_ECM": (4, 3)}[k]
        ax.annotate(lab, (np.log(r.HR_TCGA_BRCA), np.log(r.HR_METABRIC)), xytext=off, textcoords="offset points",
                    fontsize=5.0, color=col[k], ha="left")
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.tick_params(labelsize=5.2)
    ax.set_xlabel("TCGA-BRCA: adjusted log HR per SD", fontsize=5.4)
    ax.set_ylabel("METABRIC: adjusted log HR per SD", fontsize=5.4)
    agree = int((np.sign(np.log(V.HR_METABRIC)) == np.sign(np.log(V.HR_TCGA_BRCA))).sum())
    ax.set_title("METABRIC against TCGA-BRCA", fontsize=5.6, loc="left", pad=9, color=INK, fontweight="bold")
    ax.text(0.02, 0.97, f"same sign {agree}/{len(V)}", transform=ax.transAxes, fontsize=5.0, color=GREY, va="top")
    ax.plot([], [], "o", color=INK, ms=3, label="overall survival")
    ax.plot([], [], "^", color=INK, ms=3, label="disease-specific")
    ax.legend(loc="lower right", fontsize=4.9, frameon=False, handletextpad=0.3)


fig = plt.figure(figsize=(bs.W2, 150 * bs.MM))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 1], width_ratios=[1, 1], hspace=0.55, wspace=0.42,
              left=0.1, right=0.93, top=0.93, bottom=0.09)
ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0]); ax_d = fig.add_subplot(gs[1, 1])
perfect(ax_a); pooled(ax_b); metabric(ax_c); versus(ax_d)
for ax, L, dx in [(ax_a, "a", -0.3), (ax_b, "b", -0.08), (ax_c, "c", -0.3), (ax_d, "d", -0.2)]:
    bs.panel(ax, L, dx=dx, dy=1.09)
bs.save(fig, "Fig_E16_perfect_metabric")
print("saved Fig_E16_perfect_metabric")
