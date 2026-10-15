"""v17 Figure 1 - study design with the replication checkpoint cohorts, the I-SPY2 randomised
comparisons, TCGA pan-cancer prognosis and single-cell imaging (a), and the therapeutic hypothesis
that motivated the Barrier axis (b)."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, FITTED, GREY, INK, PAPER

_orig_panel = bs.panel
bs.panel = lambda ax, l: _orig_panel(ax, l, dx=0.0, dy=1.02)


def box(ax, x, y, w, h, head, body="", edge=INK, fill="white", hs=5.8, bs_=5.0, gap=0.075, top=0.018):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012",
                                fc=fill, ec=edge, lw=0.7))
    ax.text(x + w / 2, y + h - top, head, ha="center", va="top", fontsize=hs, fontweight="bold", color=edge if edge != INK else INK)
    if body:
        ax.text(x + w / 2, y + h - gap, body, ha="center", va="top", fontsize=bs_, color=INK, linespacing=1.35)


def arrow(ax, p, q, col=GREY, style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle=style, mutation_scale=6, color=col, lw=0.8, ls=ls,
                                 shrinkA=0, shrinkB=0))


def design(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
    # column 1: reference data
    ax.text(0.10, 1.0, "Define and characterise", ha="center", fontsize=5.8, color=GREY, va="top")
    box(ax, 0.005, 0.68, 0.19, 0.25, "Single-cell atlas", "GSE183904\n29 patients, 10 cell types\ncell-type attribution", gap=0.07)
    box(ax, 0.005, 0.37, 0.19, 0.25, "Visium sections", "GSE251950, 10 sections\nTLS rings, Moran's I\ntransport analysis", gap=0.07)
    box(ax, 0.005, 0.06, 0.19, 0.25, "Single-cell imaging", "IMC, breast and pancreas\n(not gastric); 1.2 M cells\nfibroblast-layer continuity", gap=0.07)
    # column 2: modules and axes
    ax.text(0.39, 1.0, "Six marker modules, two axes", ha="center", fontsize=5.8, color=GREY, va="top")
    mods = [("CAF_ECM", STROMA), ("TAM_SUP", MACRO), ("APM", IMMUNE), ("IFN", IMMUNE), ("TLS_B", IMMUNE), ("T_CELL", IMMUNE)]
    for i, (m, c) in enumerate(mods):
        xx, yy = 0.255 + (i % 2) * 0.135, 0.80 - (i // 2) * 0.095
        ax.add_patch(FancyBboxPatch((xx, yy), 0.125, 0.075, boxstyle="round,pad=0.003,rounding_size=0.01", fc="white", ec=c, lw=0.8))
        ax.text(xx + 0.0625, yy + 0.037, m, ha="center", va="center", fontsize=5.3, color=c, fontweight="bold")
    ax.text(0.39, 0.46, "Barrier = mean(CAF_ECM, TAM_SUP)", ha="center", fontsize=5.3, color=STROMA)
    ax.text(0.39, 0.39, "Immune = mean(APM, IFN, TLS_B, T_CELL)", ha="center", fontsize=5.3, color=IMMUNE)
    ax.text(0.39, 0.32, "NetImmune = Immune − Barrier", ha="center", fontsize=5.6, color=IMMUNE, fontweight="bold")
    ax.text(0.39, 0.24, "48 genes · z-mean · nothing fitted to outcome", ha="center", fontsize=5.0, color=GREY)
    arrow(ax, (0.198, 0.805), (0.248, 0.78)); arrow(ax, (0.198, 0.495), (0.248, 0.45))
    # column 3: scored unchanged
    ax.text(0.655, 1.0, "Scored unchanged in", ha="center", fontsize=5.8, color=GREY, va="top")
    H = 0.128
    ys = {"bulk": 0.815, "tcga": 0.657, "prot": 0.499, "pembro": 0.341, "ext": 0.183, "rct": 0.025}
    box(ax, 0.565, ys["bulk"], 0.18, H, "Gastric bulk RNA", "3 cohorts (1,097 patients)\n+ 3 locked (679 patients)", gap=0.052, top=0.012, hs=5.4, bs_=4.7)
    box(ax, 0.565, ys["tcga"], 0.18, H, "Other cancers, untreated", "TCGA 24 types (8,614);\nMETABRIC breast (1,979)", gap=0.052, top=0.012, hs=5.4, bs_=4.7)
    box(ax, 0.565, ys["prot"], 0.18, H, "Proteomics, gastric", "PDC000214, 130 pairs;\ndetection and co-variation only", gap=0.052, top=0.012, hs=5.4, bs_=4.7)
    box(ax, 0.565, ys["pembro"], 0.18, H, "Pembrolizumab, gastric", "PRJEB25780\n44 patients, 12 responders", gap=0.052, top=0.012, hs=5.4, bs_=4.7)
    box(ax, 0.565, ys["ext"], 0.18, H, "Checkpoint cohorts", "11 primary + 10 replication\n1,340 patients, 408 responders", gap=0.052, top=0.012, hs=5.4, bs_=4.7)
    box(ax, 0.565, ys["rct"], 0.18, H, "Randomised comparisons", "IMmotion150; I-SPY2 ×2\n611 patients", gap=0.052, top=0.012, hs=5.4, bs_=4.7)
    for k in ys:
        arrow(ax, (0.528, 0.40), (0.56, ys[k] + H / 2))
    # column 4: evaluation (dashed group); proteomics does not enter it
    ax.text(0.8975, 1.0, "Evaluated against", ha="center", fontsize=5.8, color=GREY, va="top")
    ax.add_patch(FancyBboxPatch((0.79, 0.015), 0.215, 0.945, boxstyle="round,pad=0.004,rounding_size=0.015",
                                fc="none", ec=GREY, lw=0.6, ls=(0, (3, 2))))
    H4 = 0.16
    y4 = [0.785, 0.600, 0.415, 0.230, 0.045]
    box(ax, 0.80, y4[0], 0.195, H4, "Random-gene nulls", "structure- or size-matched;\nall genes; immune/stromal", edge=INK, fill=PAPER, gap=0.06, top=0.014, hs=5.4, bs_=4.7)
    box(ax, 0.80, y4[1], 0.195, H4, "Comparators", "CD8A, ESTIMATE, TMEscore,\nGEP-18, PD-L1, IPS", edge=INK, fill=PAPER, gap=0.06, top=0.014, hs=5.4, bs_=4.7)
    box(ax, 0.80, y4[2], 0.195, H4, "Fitted models", "8 classes, leave-one-cohort-out;\nlearning curve", edge=FITTED, fill=PAPER, gap=0.06, top=0.014, hs=5.4, bs_=4.7)
    box(ax, 0.80, y4[3], 0.195, H4, "Replication", "10 new checkpoint cohorts;\nexpected directions in code", edge=INK, fill=PAPER, gap=0.06, top=0.014, hs=5.4, bs_=4.7)
    box(ax, 0.80, y4[4], 0.195, H4, "Treatment interaction", "score × arm in three\nrandomised comparisons", edge=INK, fill=PAPER, gap=0.06, top=0.014, hs=5.4, bs_=4.7)
    for k in ("bulk", "tcga", "pembro", "ext", "rct"):
        arrow(ax, (0.748, ys[k] + H / 2), (0.786, ys[k] + H / 2))


def hypothesis(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_axis_off()
    box(ax, 0.02, 0.36, 0.28, 0.50, "Barrier", "CAF/ECM remodelling\n+ suppressive TAMs", edge=STROMA, gap=0.25, top=0.05)
    box(ax, 0.38, 0.36, 0.28, 0.50, "Immune engagement", "antigen presentation, IFN,\nTLS/B cells, T cells", edge=IMMUNE, gap=0.25, top=0.05)
    box(ax, 0.74, 0.36, 0.24, 0.50, "Checkpoint response", "", edge=INK, top=0.2)
    ax.plot([0.305, 0.37], [0.61, 0.61], color=STROMA, lw=0.9)
    ax.plot([0.37, 0.37], [0.53, 0.69], color=STROMA, lw=1.2)
    arrow(ax, (0.665, 0.61), (0.735, 0.61), col=IMMUNE)
    ax.text(0.5, 0.20, "Hypothesis that motivated the axis. Tested here: whether it holds as composition;\n"
            "whether stroma obstructs by geometry: inconclusive on Visium; with single-cell imaging, CD8 exclusion\n"
            "tracked fibroblast number rather than layer arrangement (Supplementary Figs S8, S15).",
            ha="center", va="top", fontsize=5.0, color=GREY, linespacing=1.4)


cells = [dict(pos=(0, 1, 0, 1), letter="a", draw=design), dict(pos=(1, 2, 0, 1), letter="b", draw=hypothesis)]
print(bs.compose("v17_Fig1", cells, grid=(2, 1), height_ratios=[1, 0.36], width=bs.W2, height=112 * bs.MM, hspace=0.10))
