"""v37 Figure 3 - single genes as the reference (six panels; v35 Figure 2 panels a-d plus e, CD8A among the T-cell
genes in each comparison, and f, the breast comparisons against CD8A).

Previous docstring (v35 Figure 2):

(a) Pooled six-comparison interaction of every gene measured in all six comparisons (2,082), with the T-cell and
    IFN-gamma/checkpoint panels and the 52 signatures marked below.
(b) The 27 panel genes: pooled interaction (95% CI); number of the 52 signatures whose paired difference from the gene
    is > 0 at P < 0.05.
(c) The same genes in each comparison and pooled over the three breast comparisons.
(d) Number of signatures exceeding each gene against the gene's own interaction.
Numbers from analysis/results/E33_* (protocol analysis/E33_prereg.md). Writes only into figures/cns/v35/.
"""
import sys, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, MM, DIV

HERE = Path(__file__).resolve().parent
S.OUT = HERE
GUARD = {p: hashlib.md5(open(p, "rb").read()).hexdigest() for p in glob.glob(str(HERE.parent / "Fig*.p*"))}
R = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
TC, IC, SG = IMMUNE, "#7B5EA7", "#8E9AA6"
FS, FSS = 6.0, 5.6
NUM = {}
rec = lambda k, v, s: NUM.__setitem__(k, (v, s))


def title(ax, s, x=0.0):
    ax.text(x, 1.035, s, transform=ax.transAxes, fontsize=6.8, fontweight="bold", color=INK, va="bottom", ha="left")


G = pd.read_csv(R / "E33_pooled.csv")
six = G[(G.set == "six (primary)") & (G.kind == "gene")].set_index("name")
allg = six[six.all_six == True]
pan = six[six.panel.isin(["T cell", "IFN-γ/checkpoint"])].sort_values("est", ascending=False)
Q3 = pd.read_csv(R / "E33_Q3_by_reference.csv").set_index("reference")
L = pd.read_csv(R / "E32_primary.csv"); L = L[L.kind == "library"]
PG = pd.read_csv(R / "E33_per_comparison.csv"); PG = PG[(PG.kind == "gene") & (PG.variant == "primary")]
BR = G[(G.set == "breast three") & (G.kind == "gene")].set_index("name")
fig = plt.figure(figsize=(180 * MM, 206 * MM))
col_of = lambda g: TC if pan.loc[g, "panel"] == "T cell" else IC

# ------------------------------------------------------------------ a  all genes
axA = S.ax_mm(fig, 14, 9, 80, 38); S.letter_at(fig, 0, 3.0, "a"); title(axA, f"Every gene measured in all six comparisons ({len(allg):,})", -0.15)
bins = np.arange(-0.45, 0.46, 0.02)
axA.hist(allg.est.clip(-0.44, 0.44), bins=bins, color="#D5DBE1", lw=0)
p95, p99 = np.percentile(allg.est, [95, 99])
for v, lab, ha, dx in [(p95, "95th", "right", -0.006), (p99, "99th percentile", "left", 0.006)]:
    axA.axvline(v, color=SLATE, lw=0.5, ls=(0, (2, 2)))
    axA.text(v + dx, axA.get_ylim()[1] * 0.97, lab, fontsize=FSS, color=SLATE, ha=ha, va="top")
c8 = six.loc["CD8A", "est"]; pct = float(np.mean(allg.est < c8))
axA.annotate(f"CD8A\n{pct:.0%}", (c8, 3), xytext=(c8 + 0.07, 60), fontsize=FS, color=TC, ha="left", va="center", linespacing=0.95,
             arrowprops=dict(arrowstyle="-", lw=0.5, color=TC, shrinkA=0, shrinkB=0))
axA.set_xlim(-0.45, 0.45); axA.set_xticks([])
axA.set_ylabel("Genes"); axA.spines["bottom"].set_visible(False)
rec("a genes", f"{len(allg)}; median {allg.est.median():+.3f}; 95th {p95:+.3f}; 99th {p99:+.3f}", "E33_pooled.csv")
rec("a nominal", f"P<0.05 favoring {int(((allg.p < 0.05) & (allg.est > 0)).sum())}, against {int(((allg.p < 0.05) & (allg.est < 0)).sum())}", "E33_pooled.csv")
rec("a CD8A percentile", f"{pct:.3f}", "E33_pooled.csv")
axR = S.ax_mm(fig, 14, 48, 80, 9); S.clean(axR); axR.set_xlim(-0.45, 0.45); axR.set_ylim(3.2, -0.2)
for j, (lab, vals, col) in enumerate([("T-cell genes (13)", pan[pan.panel == "T cell"].est, TC),
                                      ("IFN-γ genes (14)", pan[pan.panel != "T cell"].est, IC),
                                      ("52 signatures", L.est, SG)]):
    for v in vals:
        axR.plot([v, v], [j + 0.1, j + 0.85], color=col, lw=0.7)
    axR.text(-0.45, j + 0.48, lab, fontsize=FSS, color=col, va="center", ha="left")
tm = np.median(pan[pan.panel == "T cell"].est); rec("a T-cell median percentile", f"{np.mean(allg.est < tm):.3f}", "E33_pooled.csv")
axX = S.ax_mm(fig, 14, 57.5, 80, 0.1); S.clean(axX); axX.set_xlim(-0.45, 0.45)
axX.spines["bottom"].set_visible(True); axX.set_xticks([-0.4, -0.2, 0, 0.2, 0.4]); axX.set_xticklabels(["−0.4", "−0.2", "0", "0.2", "0.4"])
axX.tick_params(axis="x", length=2); axX.set_xlabel("Pooled interaction per SD (positive favors the checkpoint arm)")

# ------------------------------------------------------------------ b  panel genes
TOP, H = 12, 104
n = len(pan)
axB = S.ax_mm(fig, 118, TOP, 26, H); S.letter_at(fig, 100, 3.0, "b"); title(axB, "T-cell and IFN-γ genes", -0.62)
axB.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
for i, (g, r) in enumerate(pan.iterrows()):
    c = col_of(g); sig = r.p < 0.05
    axB.plot([max(r.lo, -0.3), min(r.hi, 0.6)], [i, i], color=c, lw=0.7)
    axB.plot(r.est, i, "o", ms=3.0, mfc=c if sig else "white", mec=c, mew=0.7, zorder=3)
    axB.text(-0.34, i, g, fontsize=FSS, color=INK, va="center", ha="right", fontweight="bold" if g == "CD8A" else "normal")
    ex = int(Q3.loc[g, "exceeds"])
    axB.text(0.70, i, str(ex), fontsize=FSS, color=INK if ex <= 2 else "#C1401A", va="center", ha="center",
             fontweight="normal" if ex <= 2 else "bold")
    rec(f"b {g}", f"{r.est:+.3f} ({r.lo:.3f} to {r.hi:.3f}); P = {r.p:.3g}; k = {int(r.k)}; signatures above {ex}", "E33_pooled.csv; E33_Q3_by_reference.csv")
axB.set_ylim(n - 0.4, -0.7); axB.set_yticks([]); axB.spines["left"].set_visible(False)
axB.set_xlim(-0.3, 0.6); axB.set_xticks([0, 0.3, 0.6]); axB.set_xticklabels(["0", "0.3", "0.6"])
axB.set_xlabel("Pooled interaction")
axB.text(0.70, -1.0, "signa-\ntures\nabove", fontsize=FSS - 0.4, color=SLATE, ha="center", va="bottom", linespacing=0.9)

# ------------------------------------------------------------------ c  per comparison and breast pool
COMP = [("IMM150", "IMmotion150"), ("CM025", "CheckMate 025"), ("JR101", "JAVELIN 101"), ("ISPY2_pembro", "I-SPY2 pembro"),
        ("ISPY2_durva", "I-SPY2 durva"), ("NeoTRIP", "NeoTRIP"), ("breast", "Breast, pooled")]
axC = S.ax_mm(fig, 152, TOP, 25, H); S.letter_at(fig, 147.5, 3.0, "c"); title(axC, "Each comparison", 0.0)
norm = Normalize(-1.0, 1.0)
M = np.full((n, 7), np.nan); Pv = np.full((n, 7), np.nan)
for i, g in enumerate(pan.index):
    for j, (k, _) in enumerate(COMP):
        if k == "breast":
            if g in BR.index:
                M[i, j], Pv[i, j] = BR.loc[g, "est"], BR.loc[g, "p"]
        else:
            q = PG[(PG.name == g) & (PG.comparison == k)]
            if len(q):
                M[i, j], Pv[i, j] = q.y.iloc[0], q.p.iloc[0]
cmap = DIV.reversed().copy(); cmap.set_bad("#E3E7EB")
axC.imshow(np.ma.masked_invalid(np.clip(M, -1, 1)), cmap=cmap, norm=norm, aspect="auto", interpolation="nearest",
           extent=(-0.5, 6.5, n - 0.5, -0.5))
for j in range(8):
    axC.axvline(j - 0.5, color="white", lw=0.6 if j != 6 else 1.6)
yy, xx = np.where(Pv < 0.05); axC.scatter(xx, yy, s=2.0, c=INK, lw=0, zorder=3)
axC.set_xlim(-0.5, 6.5); axC.set_ylim(n - 0.4, -0.7); axC.set_yticks([])
axC.set_xticks(range(7)); axC.set_xticklabels([c[1] for c in COMP], rotation=90, fontsize=FSS)
axC.tick_params(axis="x", length=0, pad=1.5)
for s_ in axC.spines.values():
    s_.set_visible(False)
for g in ["CD8A", "GZMB", "IDO1", "CXCL10", "STAT1"]:
    if g in BR.index:
        rec(f"c breast {g}", f"{BR.loc[g, 'est']:+.3f} ({BR.loc[g, 'lo']:.3f} to {BR.loc[g, 'hi']:.3f}); P = {BR.loc[g, 'p']:.3g}", "E33_pooled.csv")
cb = S.ax_mm(fig, 152, TOP + H + 15.5, 25, 2.0)
cb.imshow(np.linspace(-1, 1, 200)[None, :], cmap=DIV.reversed(), norm=norm, aspect="auto", extent=(-1, 1, 0, 1))
cb.set_yticks([]); cb.set_xticks([-1, 0, 1]); cb.set_xticklabels(["≤ −1", "0", "≥ 1"], fontsize=FSS)
for s_ in cb.spines.values():
    s_.set_visible(False)
cb.tick_params(axis="x", length=1.5, pad=1)

# ------------------------------------------------------------------ d  signatures above each gene
axD = S.ax_mm(fig, 14, 76, 62, 40); S.letter_at(fig, 0, 69.5, "d"); title(axD, "Signatures significantly above each gene", -0.2)
for g, r in pan.iterrows():
    ex = int(Q3.loc[g, "exceeds"])
    axD.plot(r.est, ex, "o", ms=3.4, mfc=col_of(g), mec="white", mew=0.4, zorder=3)
    if ex > 2 or g == "CD8A":
        axD.annotate(g, (r.est, ex), xytext=(4, 3), textcoords="offset points", fontsize=FSS, color=INK)
axD.axhline(2.5, color=GREY, lw=0.5, ls=(0, (2, 2)))
axD.text(0.29, 3.2, "more than 2 of 52", fontsize=FSS, color=SLATE, ha="right", va="bottom")
axD.set_xlim(0, 0.30); axD.set_xticks([0, 0.1, 0.2, 0.3]); axD.set_ylim(-1, 19); axD.set_yticks([0, 5, 10, 15])
axD.set_xlabel("The gene's pooled interaction per SD"); axD.set_ylabel("Signatures above the gene\n(paired, P < 0.05; of 52)")
nt = int((Q3.loc[pan[pan.panel == "T cell"].index, "exceeds"] <= 2).sum()); na = int((Q3.loc[pan.index, "exceeds"] <= 2).sum())
rec("d genes exceeded by at most 2", f"T cell {nt}/13; all panel {na}/{n}", "E33_Q3_by_reference.csv")
K = S.ax_mm(fig, 80, 78, 18, 14); S.clean(K); K.set_xlim(0, 18); K.set_ylim(14, 0)
for j, (lab, c) in enumerate([("T cell", TC), ("IFN-γ /\ncheckpoint", IC)]):
    K.plot(1, 2 + 5 * j, "o", ms=3.4, mfc=c, mec="white"); K.text(2.6, 2 + 5 * j, lab, fontsize=FSS, color=SLATE, va="center", linespacing=0.95)
K.text(0, 12.5, "b, filled: P < 0.05\nc, dot: P < 0.05", fontsize=FSS, color=SLATE, va="center", linespacing=1.0)


# ------------------------------------------------------------------ e  CD8A among the T-cell genes in each comparison
T13 = list(pan[pan.panel == "T cell"].index)
ROWS = [("IMM150", "IMmotion150"), ("CM025", "CheckMate 025"), ("JR101", "JAVELIN Renal 101"), ("ISPY2_pembro", "I-SPY2 pembrolizumab"),
        ("ISPY2_durva", "I-SPY2 durvalumab–olaparib"), ("NeoTRIP", "NeoTRIP"), ("six", "Six pooled"), ("five", "Five, without IMmotion150"),
        ("breast", "Breast, pooled")]
axE = S.ax_mm(fig, 40, 152, 52, 44); S.letter_at(fig, 0, 145.5, "e"); title(axE, "CD8A among the 13 T-cell genes", -0.68)
G5 = G[(G.set == "five without IMmotion150") & (G.kind == "gene")].set_index("name")
for j, (k, lab) in enumerate(ROWS):
    if k == "six":
        v = six.loc[T13, "est"]
    elif k == "five":
        v = G5.loc[T13, "est"]
    elif k == "breast":
        v = BR.loc[T13, "est"]
    else:
        q = PG[PG.comparison == k].set_index("name"); v = q.loc[T13, "y"]
    rank = int((v > v["CD8A"]).sum()) + 1
    axE.scatter(v.drop("CD8A").clip(-0.9, 0.9), [j] * 12, s=7, c=GREY, edgecolors="none", zorder=2)
    axE.scatter([min(max(v["CD8A"], -0.9), 0.9)], [j], s=16, c=TC, edgecolors="white", linewidths=0.4, zorder=3)
    axE.text(-1.0, j, lab, fontsize=FSS, color=INK, va="center", ha="right", fontweight="bold" if k in ("six", "five", "breast") else "normal")
    axE.text(1.02, j, f"{rank}/13", fontsize=FSS, color=TC, va="center", ha="left")
    rec(f"e {lab}", f"CD8A {v['CD8A']:+.3f}, rank {rank}/13; top {v.idxmax()} {v.max():+.3f}", "E33_per_comparison.csv; E33_pooled.csv")
    if j == 5:
        axE.axhline(5.5, color=LIGHT, lw=0.6)
axE.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
axE.set_ylim(len(ROWS) - 0.4, -0.6); axE.set_yticks([]); axE.spines["left"].set_visible(False)
axE.set_xlim(-0.9, 0.9); axE.set_xticks([-0.5, 0, 0.5]); axE.set_xticklabels(["−0.5", "0", "0.5"])
axE.set_xlabel("Interaction per SD"); axE.text(1.02, -1.0, "CD8A\nrank", fontsize=FSS, color=SLATE, ha="left", va="bottom", linespacing=0.9)

# ------------------------------------------------------------------ f  breast comparisons: each signature minus CD8A
PR = pd.read_csv(R / "E33_paired.csv"); cb = PR[(PR.set == "breast three") & (PR.reference == "CD8A")].sort_values("diff", ascending=False).reset_index(drop=True)
FAMC = {"Inflamed/IFN": IMMUNE, "Other": "#8E9AA6", "Confounder": STROMA}
fam = dict(zip(L.signature, L.family))
axF = S.ax_mm(fig, 118, 152, 58, 44); S.letter_at(fig, 104, 145.5, "f"); title(axF, "Breast comparisons: each signature minus CD8A", -0.2)
axF.axhline(0, color=INK, lw=0.5)
NICE = {"IFN_Ayers": "IFN (Ayers)", "Teff-IFNG_Fehrenbacher": "Teff-IFNG (Fehrenbacher)", "VIGex_Hernando-Calvo": "VIGex",
        "Teff_McDermott": "Teff (McDermott)", "CTL_Jiang": "CTL (Jiang)", "Bcell_Budczies": "B cell (Budczies)"}
for i, r in cb.iterrows():
    c = FAMC[fam[r.signature]]; sig = (r["diff"] > 0) and (r.p < 0.05)
    axF.plot([i, i], [max(r.lo, -1.2), min(r.hi, 1.2)], color=c, lw=0.6)
    axF.plot(i, r["diff"], "o", ms=2.6, mfc=c if sig else "white", mec=c, mew=0.5, zorder=3)
    if sig:
        rec(f"f {r.signature}", f"{r['diff']:+.3f} ({r.lo:.3f} to {r.hi:.3f}); P = {r.p:.3g}", "E33_paired.csv")
labs = [r.signature for _, r in cb.iterrows() if r["diff"] > 0 and r.p < 0.05]
axF.text(0.98, 0.97, "above CD8A (filled):\n" + "\n".join(NICE.get(x, x) for x in labs), transform=axF.transAxes, fontsize=FSS - 0.2,
         color=INK, va="top", ha="right", linespacing=1.0)
axF.set_xlim(-1, len(cb)); axF.set_xticks([0, 9, 19, 29, 39, 49]); axF.set_xticklabels(["1", "10", "20", "30", "40", "50"])
axF.set_ylim(-1.2, 1.2); axF.set_yticks([-1, -0.5, 0, 0.5, 1]); axF.set_yticklabels(["−1", "−0.5", "0", "0.5", "1"])
axF.set_xlabel("Signatures ranked"); axF.set_ylabel("Signature − CD8A,\nthree breast comparisons")

S.save(fig, "Fig3")
with open(HERE / "Fig3_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched")
