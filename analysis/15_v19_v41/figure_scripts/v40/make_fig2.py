"""v40 Figure 2 - the calibration framework and its operating characteristics (E36 simulation).
(a) the three references; (b) null and prognostic-only rejection rates; (c) power; (d) random-gene reference;
(e) paired test with CD8A as the truth; (f) paired-test power against the true difference, with normal theory.
Numbers from analysis/results/E36_*. Writes only into figures/cns/v40/.
"""
import sys, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.stats import norm
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, MACRO, MM, IMMUNE_LIGHT, PAPER

HERE = Path(__file__).resolve().parent
S.OUT = HERE
GUARD = {p: hashlib.md5(open(p, "rb").read()).hexdigest() for p in glob.glob(str(HERE.parent / "Fig*.p*"))}
R = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
FS, FSS = 6.0, 5.6
NUM = {}
rec = lambda k, v, s: NUM.__setitem__(k, (v, s))
LAB = {"CD8A": "CD8A", "T_CELL": "T_CELL", "IFN": "IFN", "APM": "APM", "TLS_B": "TLS_B", "CAF_ECM": "CAF_ECM",
       "NetImmune": "NetImmune"}
MK = ["CD8A", "T_CELL", "IFN", "NetImmune", "TLS_B", "APM", "CAF_ECM"]
COL = {m: (INK if m == "CD8A" else STROMA if m == "CAF_ECM" else IMMUNE) for m in MK}


def title(ax, s, x=0.0):
    ax.text(x, 1.05, s, transform=ax.transAxes, fontsize=6.8, fontweight="bold", color=INK, va="bottom", ha="left")


B = pd.read_csv(R / "E36_S0_S2.csv"); S3 = pd.read_csv(R / "E36_S3.csv"); S4 = pd.read_csv(R / "E36_S4.csv")
fig = plt.figure(figsize=(180 * MM, 184 * MM))

# ------------------------------------------------------------------ a  the framework
A = S.ax_mm(fig, 4, 6, 172, 40); S.clean(A); A.set_xlim(0, 172); A.set_ylim(40, 0); S.letter_at(fig, 0, 3.5, "a")
A.text(0, -1.5, "Three references for a predictive claim", fontsize=6.8, fontweight="bold", color=INK, va="bottom")


def box(x, y, w, h, head, lines, face=PAPER, edge=LIGHT):
    A.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6", fc=face, ec=edge, lw=0.6))
    A.text(x + 2.2, y + 3.4, head, fontsize=6.0, fontweight="bold", color=INK, va="center")
    for i, l in enumerate(lines):
        A.text(x + 2.2, y + 8.6 + 4.4 * i, l, fontsize=FSS, color=SLATE, va="center")


box(0, 3, 30, 34, "Input", ["Randomized comparisons", "with pre-treatment", "expression from both", "arms; a score fixed", "before it is tested"])
heads = [("1  Treatment interaction", ["score × arm per SD, benefit-", "oriented (−log HR, log OR);", "pooled over trials (wider of", "two random-effects intervals)", "Separates predictive from", "prognostic"]),
         ("2  Random genes", ["rank among 1,000 random", "gene sets of the same size", "and weights, drawn in each", "trial, pooled draw by draw", "Separates specific from", "shared signal"]),
         ("3  Single genes", ["paired difference from each", "gene of a T-cell and IFN-γ", "panel (within-arm bootstrap),", "with the difference detectable", "at 80% power", ""])]
for i, (h, ls) in enumerate(heads):
    x = 36 + 34 * i
    box(x, 3, 31, 34, h, ls, face=IMMUNE_LIGHT, edge="#B9D3EE")
    A.annotate("", xy=(x - 0.6, 20), xytext=(x - 5.2, 20), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.6, mutation_scale=5))
box(142, 3, 30, 34, "Output", ["Interaction (95% CI)", "Specific: yes / no", "Above or below each gene", "Detectable difference", "Operating characteristics", "by simulation (b–f)"])
A.annotate("", xy=(141.4, 20), xytext=(137.2, 20), arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.6, mutation_scale=5))

# ------------------------------------------------------------------ b  null and prognostic-only
axB = S.ax_mm(fig, 22, 60, 40, 30); S.letter_at(fig, 0, 55, "b"); title(axB, "False predictive claims", -0.45)
s0 = B[(B.scenario == "S0") & (B.test == "interaction")].set_index("marker")
s1 = B[(B.scenario == "S1")]
for k, m in enumerate(MK):
    v = 100 * s0.loc[m].reject_p05
    axB.barh(k, v, height=0.62, color=COL[m], alpha=0.85); axB.text(v + 0.3, k, f"{v:.1f}", fontsize=FSS, color=SLATE, va="center")
    rec(f"b null {m}", f"{v:.1f}%", "E36_S0_S2.csv")
axB.set_yticks(range(len(MK))); axB.set_yticklabels(MK, fontsize=FSS); axB.invert_yaxis()
axB.set_xlim(0, 10); axB.axvline(5, color=GREY, lw=0.5, ls=(0, (2, 2))); axB.set_xticks([0, 5, 10])
axB.text(1.0, 1.0, "null truth", transform=axB.transAxes, fontsize=FSS, color=SLATE, ha="right", va="bottom")
axB2 = S.ax_mm(fig, 22, 99, 40, 9)
caf_i = s1[(s1.marker == "CAF_ECM") & (s1.test == "interaction")].iloc[0].reject_p05
caf_w = s1[(s1.marker == "CAF_ECM") & (s1.test == "checkpoint arm alone")].iloc[0].reject_p05
axB2.barh(0, 100 * caf_i, height=0.62, color=STROMA); axB2.text(7, 0, f"{100 * caf_i:.1f}", fontsize=FSS, color=SLATE, va="center")
axB2.barh(1, 100 * caf_w, height=0.62, color=STROMA, hatch="//////", edgecolor="white", lw=0)
axB2.text(100 * caf_w - 2, 1, f"{100 * caf_w:.1f}", fontsize=FSS, color="white", va="center", ha="right", fontweight="bold")
rec("b prognostic CAF_ECM interaction", f"{100 * caf_i:.1f}%", "E36_S0_S2.csv"); rec("b prognostic CAF_ECM checkpoint arm alone", f"{100 * caf_w:.1f}%", "E36_S0_S2.csv")
axB2.set_yticks([0, 1]); axB2.set_yticklabels(["interaction", "checkpoint arm alone"], fontsize=FSS)
axB2.invert_yaxis(); axB2.set_xlim(0, 100); axB2.axvline(5, color=GREY, lw=0.5, ls=(0, (2, 2)))
axB2.set_xlabel("Rejected at P < 0.05 (%)")
axB2.text(1.0, 1.08, "prognostic-only CAF_ECM (adverse, 0.26/SD)", transform=axB2.transAxes, fontsize=FSS, color=SLATE, ha="right", va="bottom")

# ------------------------------------------------------------------ c  power
axC = S.ax_mm(fig, 82, 60, 36, 48); S.letter_at(fig, 66, 55, "c"); title(axC, "Power (truth: CD8A × arm)", -0.25)
bet = [0, 0.1, 0.2, 0.3]
ends = {}
for m in MK:
    ys = [100 * s0.loc[m].reject_p05] + [100 * B[(B.scenario == f"S2 beta {b}") & (B.marker == m)].iloc[0].reject_p05 for b in bet[1:]]
    lw = 1.1 if m == "CD8A" else 0.7
    axC.plot(bet, ys, color=COL[m], lw=lw, alpha=1 if m in ("CD8A", "CAF_ECM") else 0.75, marker="o", ms=1.8)
    ends[m] = ys[-1]; rec(f"c power {m}", ", ".join(f"{y:.1f}" for y in ys), "E36_S0_S2.csv")
order = sorted(MK, key=lambda m: -ends[m]); ylab = {}
last = 200
for m in order:
    y = min(ends[m], last - 5.6); ylab[m] = y; last = y
for m in MK:
    axC.text(0.307, ylab[m], f"{LAB[m]} {ends[m]:.0f}%", fontsize=5.2, color=INK if m == "CD8A" else SLATE, va="center", ha="left")
axC.set_xlim(-0.01, 0.3); axC.set_ylim(0, 100); axC.set_xticks(bet); axC.set_xticklabels(["0", "0.1", "0.2", "0.3"])
axC.set_xlabel("True interaction of CD8A (per SD)"); axC.set_ylabel("Pooled test, P < 0.05 (%)")
axC.axhline(5, color=GREY, lw=0.5, ls=(0, (2, 2)))

# ------------------------------------------------------------------ d  random-gene reference
axD = S.ax_mm(fig, 22, 130, 40, 42); S.letter_at(fig, 0, 124, "d"); title(axD, "Exceeds random genes", -0.45)
SC3 = [("S3 null", "Null", GREY), ("S3 CD8A truth", "CD8A truth", IMMUNE), ("S3d diffuse (PC1) truth", "Diffuse truth (PC1)", MACRO)]
for j, (sc, lab, c) in enumerate(SC3):
    t = S3[S3.scenario == sc].set_index("marker")
    for k, m in enumerate(MK):
        v = 100 * t.loc[m].exceed_rate
        axD.barh(k + (j - 1) * 0.27, v, height=0.25, color=c)
        rec(f"d {lab} {m}", f"{v:.1f}%", "E36_S3.csv")
axD.set_yticks(range(len(MK))); axD.set_yticklabels([LAB[m] + (" (−)" if m == "CAF_ECM" else "") for m in MK], fontsize=FSS)
axD.invert_yaxis(); axD.set_xlim(0, 100); axD.axvline(5, color=GREY, lw=0.5, ls=(0, (2, 2)))
axD.set_xlabel("Above the 95th percentile of 200 random sets (%)")
for j, (sc, lab, c) in enumerate(SC3):
    xx_, yy_ = [(42, 5.62), (68, 5.62), (42, 6.32)][j]
    axD.plot(xx_, yy_, "s", ms=2.4, color=c); axD.text(xx_ + 3, yy_, lab, fontsize=5.0, color=SLATE, va="center")

# ------------------------------------------------------------------ e  paired test, CD8A truth
axE = S.ax_mm(fig, 82, 130, 30, 42); S.letter_at(fig, 66, 124, "e"); title(axE, "Paired test, CD8A truth", -0.35)
t = S4[S4.truth == "CD8A"].set_index("module"); PM = ["T_CELL", "IFN", "APM", "TLS_B"]
for k, m in enumerate(PM):
    a_, b_ = 100 * t.loc[m].sig_above, 100 * t.loc[m].sig_below
    axE.barh(k, -b_, height=0.55, color=GREY); axE.barh(k, a_, height=0.55, color=IMMUNE)
    axE.text(-b_ - 1.5, k, f"{b_:.1f}", fontsize=FSS, color=SLATE, va="center", ha="right")
    axE.text(max(a_, 0) + 1.5, k, f"{a_:.1f}", fontsize=FSS, color=SLATE, va="center", ha="left")
    rec(f"e {m} above / below CD8A", f"{a_:.1f}% / {b_:.1f}%", "E36_S4.csv")
axE.axvline(0, color=INK, lw=0.6); axE.set_yticks(range(len(PM))); axE.set_yticklabels(PM, fontsize=FSS); axE.invert_yaxis()
axE.set_xlim(-62, 22); axE.set_xticks([-50, -25, 0, 20]); axE.set_xticklabels(["50", "25", "0", "20"])
axE.set_xlabel("Significantly below | above CD8A (%)")

# ------------------------------------------------------------------ f  power vs true difference
axF = S.ax_mm(fig, 128, 130, 46, 42); S.letter_at(fig, 116, 124, "f"); title(axF, "Paired test, module truth", -0.2)
xx = np.linspace(0, 0.4, 200)
for se_, lab, ls in [(0.029, "SE 0.029", "-"), (0.062, "SE 0.062", "--"), (0.076, "SE 0.076 (median, published signatures)", ":")]:
    axF.plot(xx, 100 * norm.cdf(xx / se_ - 1.96), color=GREY, lw=0.6, ls=ls)
own = S4[S4.truth == S4.module]
for r in own.itertuples():
    axF.plot(r.oracle_diff, 100 * r.sig_above, "o", ms=3.4, color=IMMUNE, mec="white", mew=0.4, zorder=4)
    axF.text(r.oracle_diff - 0.008, 100 * r.sig_above + 2.5, r.module, fontsize=5.2, color=SLATE, ha="right")
    rec(f"f {r.module} oracle / power / predicted", f"{r.oracle_diff:+.3f} / {100 * r.sig_above:.1f}% / "
        f"{100 * norm.cdf(r.oracle_diff / r.median_se - 1.96):.1f}%", "E36_S4.csv")
axF.axvline(0.213, color=MACRO, lw=0.6, ls=(0, (2, 2))); axF.axhline(80, color=LIGHT, lw=0.5)
axF.text(0.218, 4, "0.21", fontsize=5.2, color=MACRO)
axF.text(0.395, 30, "normal theory for\npaired SE 0.029 (solid),\n0.062 (dashed) and\n0.076 (dotted; median of\npublished signatures)", fontsize=4.9, color=SLATE, ha="right", va="top")
axF.set_xlim(0, 0.4); axF.set_ylim(0, 100); axF.set_xlabel("True difference from CD8A (per SD)"); axF.set_ylabel("Power (%)")

S.save(fig, "Fig2")
with open(HERE / "Fig2_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched")
