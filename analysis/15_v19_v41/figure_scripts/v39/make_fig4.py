"""v39 Figure 4 - how far the result can be trusted (six panels; c adds NeoPembrOV, E34): (a) random gene sets reach significance; (b) every
single gene (volcano); (c) sensitivity analyses; (d) the breast comparisons; (e) adjustment for ESTIMATE immune; (f) I-SPY2
MammaPrint and processing checks (from v34 Fig. 3e).

Previous docstring (v35 Figure 3):

(a) Random gene sets of each signature's size reaching pooled P < 0.05 (1,000 draws each), by signature size.
(b) Signatures at P < 0.05 (favoring / against the checkpoint arm) and above CD8A under each sensitivity analysis.
(c) The three breast comparisons (chemotherapy with or without a checkpoint inhibitor) against all six.
(d) Pooled interaction before and after adjustment for the ESTIMATE immune score and its interaction with arm.
Numbers from analysis/results/E32_* and E33_*. The I-SPY2 checks and the per-comparison random bands are in
Supplementary Fig. S31 (v34 Fig. 3). Writes only into figures/cns/v35/.
"""
import sys, re, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, FAV, UNFAV, MM

HERE = Path(__file__).resolve().parent
S.OUT = HERE
GUARD = {p: hashlib.md5(open(p, "rb").read()).hexdigest() for p in glob.glob(str(HERE.parent / "Fig*.p*"))}
R = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
OTHER = "#8E9AA6"
FAMCOL = {"Inflamed/IFN": IMMUNE, "Other": OTHER, "Confounder": STROMA}
FS, FSS = 6.0, 5.6
NUM = {}
rec = lambda k, v, s: NUM.__setitem__(k, (v, s))


def nice(s):
    lab = {"TMEscore_z": "TMEscore", "Immune": "Immune axis"}
    if s in lab:
        return lab[s]
    m = re.match(r"^(.*)_([A-Z][A-Za-z\-]+)$", s)
    return f"{m.group(1).replace('_', ' ')} ({m.group(2)})" if m else s.replace("_", " ")


def title(ax, s, x=0.0):
    ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=6.8, fontweight="bold", color=INK, va="bottom", ha="left")


P = pd.read_csv(R / "E32_primary.csv"); L = P[P.kind == "library"].copy()
fig = plt.figure(figsize=(180 * MM, 214 * MM))

# ------------------------------------------------------------------ a  random genes reach significance
COV = pd.read_csv(R / "E32_coverage.csv")
SIX = ["IMM150", "ISPY2_pembro", "ISPY2_durva", "CM025", "JR101", "NeoTRIP"]
size = COV[COV.scored & COV.comparison.isin(SIX)].groupby("signature").present.mean()
L["size"] = L.signature.map(size)
axA = S.ax_mm(fig, 14, 10, 54, 46); S.letter_at(fig, 0, 3.5, "a"); title(axA, "Random genes reach significance", -0.2)
axA.axhline(5, color=GREY, lw=0.5, ls=(0, (2, 2)))
axA.scatter(L["size"], 100 * L.ref_frac_nominal_p05, s=11, c=[FAMCOL[f] for f in L.family], edgecolors="white", linewidths=0.3, zorder=3)
axA.set_xscale("log"); axA.set_xticks([3, 10, 30, 100, 300]); axA.set_xticklabels(["3", "10", "30", "100", "300"])
axA.set_xlabel("Genes in the signature"); axA.set_ylabel("Random sets of that size\nwith pooled P < 0.05 (%)")
axA.set_ylim(0, 23); axA.text(330, 5.6, "5% expected", fontsize=FSS, color=SLATE, ha="right")
rho = pd.Series(L["size"]).corr(L.ref_frac_nominal_p05, method="spearman")
axA.text(0.03, 0.97, f"Spearman ρ = {rho:.2f}", transform=axA.transAxes, fontsize=FSS, color=SLATE, va="top")
for lab, col in FAMCOL.items():
    axA.scatter([], [], s=11, c=col, label=lab)
axA.legend(loc="upper left", bbox_to_anchor=(0.0, 0.88), frameon=False, fontsize=FSS, handletextpad=0.2, borderaxespad=0.2)
rec("a random nominal P<0.05", f"median {100 * L.ref_frac_nominal_p05.median():.1f}%, range {100 * L.ref_frac_nominal_p05.min():.1f}-{100 * L.ref_frac_nominal_p05.max():.1f}%; rho {rho:.2f}", "E32_primary.csv; E32_coverage.csv")


# ------------------------------------------------------------------ b  every single gene (volcano)
Gp = pd.read_csv(R / "E33_pooled.csv"); g6 = Gp[(Gp.set == "six (primary)") & (Gp.kind == "gene") & (Gp.all_six == True)]
axV = S.ax_mm(fig, 92, 10, 56, 46); S.letter_at(fig, 78, 3.5, "b"); title(axV, "Every gene measured in all six comparisons", -0.2)
ly = -np.log10(g6.p.clip(lower=1e-6))
sig = g6.p < 0.05
axV.scatter(g6.est[~sig], ly[~sig], s=2.0, c="#C9D0D7", edgecolors="none", zorder=1)
axV.scatter(g6.est[sig], ly[sig], s=2.4, c=SLATE, edgecolors="none", zorder=2)
TCELL = ["CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "CD2", "CD247", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7"]
gt = g6[g6.name.isin(TCELL)]
axV.scatter(gt.est, -np.log10(gt.p), s=7, c=IMMUNE, edgecolors="white", linewidths=0.3, zorder=3)
c8 = g6[g6.name == "CD8A"].iloc[0]
axV.annotate("CD8A", (c8.est, -np.log10(c8.p)), xytext=(6, 6), textcoords="offset points", fontsize=FSS, color=IMMUNE,
             arrowprops=dict(arrowstyle="-", lw=0.4, color=IMMUNE, shrinkA=0, shrinkB=1.5))
axV.axhline(-np.log10(0.05), color=GREY, lw=0.5, ls=(0, (2, 2))); axV.axvline(0, color=GREY, lw=0.4)
nf, na = int((sig & (g6.est > 0)).sum()), int((sig & (g6.est < 0)).sum())
axV.text(0.03, 0.97, f"P < 0.05: {nf} favor, {na} against\n(of {len(g6):,}; blue, T-cell genes)", transform=axV.transAxes,
         fontsize=FSS, color=INK, va="top")
axV.set_xlim(-0.45, 0.45); axV.set_xticks([-0.4, -0.2, 0, 0.2, 0.4]); axV.set_xticklabels(["−0.4", "−0.2", "0", "0.2", "0.4"])
axV.set_xlabel("Pooled interaction per SD"); axV.set_ylabel("−log₁₀ P")
rec("b genes at P<0.05", f"favor {nf}, against {na} of {len(g6)}", "E33_pooled.csv")
# ------------------------------------------------------------------ b  sensitivity
A = pd.read_csv(R / "E32_pooled_all_sets.csv"); NV = pd.read_csv(R / "E32_native.csv")
G = pd.read_csv(R / "E33_pooled.csv"); PR = pd.read_csv(R / "E33_paired.csv")
above = lambda st: int(((PR.set == st) & (PR.reference == "CD8A") & (PR["diff"] > 0) & (PR.p < 0.05)).sum())
SETS = [("six (primary)", "Primary (six comparisons)", above("six (primary)")),
        ("six + PRINCE + DUTRENEO", "+ PRINCE and DUTRENEO", None),
        ("six, renal adjusted for angiogenesis", "Renal adjusted for angiogenesis", None),
        ("six, I-SPY2 hormone receptor only", "I-SPY2 without MammaPrint", None),
        ("native", "Published scoring engines", None),
        ("six, adjusted for ESTIMATE immune (Q4)", "Adjusted for ESTIMATE immune", None),
        ("five without IMmotion150", "Without IMmotion150", above("five without IMmotion150")),
        ("six, CheckMate 025 overall survival", "CheckMate 025 overall survival", None),
        ("breast three", "Breast comparisons only", above("breast three")),
        ("e34 seven", "+ NeoPembrOv (seven)", None), ("e34 four", "Chemotherapy backbone (four)", None)]
P34 = pd.read_csv(R / "E34_pooled.csv")
axB = S.ax_mm(fig, 50, 78, 30, 52); S.letter_at(fig, 0, 71.5, "c"); title(axB, "Sensitivity", -1.45)
for j, (key, lab, ab) in enumerate(SETS):
    if key.startswith("e34"):
        g = P34[P34.set.str.startswith("seven" if key == "e34 seven" else "four") & (P34.kind == "library")]
        ab = int(((g.diff_vs_CD8A > 0) & (g.cd8a_p < 0.05)).sum())
    elif key == "native":
        g = NV
    elif key == "breast three":
        g = G[(G.set == key) & (G.kind == "library")]
    else:
        g = A[(A.set == key) & (A.kind == "library")]
    fav = int(((g.p < 0.05) & (g.est > 0)).sum()); ag = int(((g.p < 0.05) & (g.est < 0)).sum())
    axB.barh(j, fav, height=0.62, color=FAV, lw=0); axB.barh(j, -ag, height=0.62, color=UNFAV, lw=0)
    axB.text(fav + 0.6, j, f"{fav} / {ag}", fontsize=FSS, color=INK, va="center")
    axB.text(-4.2, j, lab, fontsize=FSS, color=INK, va="center", ha="right")
    axB.text(30.5, j, "—" if ab is None else str(ab), fontsize=FS, color=INK if not ab else IMMUNE, va="center", ha="center",
             fontweight="bold" if ab else "normal")
    rec(f"c {lab}", f"favor {fav}, against {ag}; above CD8A {ab}", "E32_pooled_all_sets.csv / E32_native.csv / E33_pooled.csv / E33_paired.csv")
axB.set_ylim(len(SETS) - 0.4, -0.6); axB.set_yticks([]); axB.spines["left"].set_visible(False)
axB.set_xlim(-4, 27); axB.set_xticks([0, 10, 20]); axB.axvline(0, color=INK, lw=0.5)
axB.set_xlabel("Signatures at P < 0.05 (of 52)")
axB.text(25, -0.95, "favoring / against", fontsize=FSS, color=SLATE, ha="right", va="center")
axB.text(30.5, -1.05, "above\nCD8A", fontsize=FSS, color=SLATE, ha="center", va="bottom", linespacing=0.9)

# ------------------------------------------------------------------ c  breast against six
axC = S.ax_mm(fig, 118, 78, 56, 46); S.letter_at(fig, 104, 71.5, "d"); title(axC, "The three breast comparisons", -0.2)
br = G[(G.set == "breast three") & (G.kind == "library")].set_index("name")
BX = pd.read_csv(R / "E33_breast_extras.csv").set_index("name")
lim = (-0.85, 0.95)
axC.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2))); axC.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
c8s, c8b = P.set_index("signature").loc["CD8A", "est"], BX.loc["CD8A", "est"]
pairs = PR[(PR.set == "breast three") & (PR.reference == "CD8A")].set_index("signature")
for s, r in L.set_index("signature").iterrows():
    if s not in br.index:
        continue
    col = FAMCOL[r.family]; ab = s in pairs.index and pairs.loc[s, "diff"] > 0 and pairs.loc[s, "p"] < 0.05
    axC.plot(r.est, br.loc[s, "est"], "o", ms=3.4 if ab else 2.6, mfc=col if ab else "white", mec=col, mew=0.6, zorder=3)
    if ab:
        rec(f"d above CD8A in breast: {s}", f"breast {br.loc[s, 'est']:+.3f}; minus CD8A {pairs.loc[s, 'diff']:+.3f} (P = {pairs.loc[s, 'p']:.3g})", "E33_pooled.csv; E33_paired.csv")
for s, lab, dxy in [("CD8A", "CD8A", (5, 5)), ("NetImmune", "NetImmune", (-5, 5)), ("CAF_ECM", "CAF_ECM", (5, -2)),
                    ("ESTIMATE", "ESTIMATE", (-4, 8))]:
    x = P.set_index("signature").loc[s, "est"]; y = BX.loc[s, "est"]
    axC.plot(x, y, "D", ms=3.6, mfc=IMMUNE if s == "CD8A" else INK, mec="white", mew=0.4, zorder=4)
    axC.annotate(lab, (x, y), xytext=dxy, textcoords="offset points", fontsize=FSS, color=IMMUNE if s == "CD8A" else INK, va="center",
                 ha="right" if dxy[0] < 0 else "left")
    rec(f"d {s}", f"six {x:+.3f}; breast {y:+.3f} ({BX.loc[s, 'lo']:.3f} to {BX.loc[s, 'hi']:.3f}); P = {BX.loc[s, 'p']:.3g}", "E32_primary.csv; E33_breast_extras.csv")
axC.set_xlim(-0.3, 0.38); axC.set_ylim(*lim)
axC.set_xticks([-0.2, 0, 0.2]); axC.set_xticklabels(["−0.2", "0", "0.2"])
axC.set_yticks([-0.8, -0.4, 0, 0.4, 0.8]); axC.set_yticklabels(["−0.8", "−0.4", "0", "0.4", "0.8"])
axC.set_xlabel("Pooled interaction, six comparisons"); axC.set_ylabel("Pooled interaction,\nthree breast comparisons")
nb = int(((br.p < 0.05) & (br.est > 0)).sum()); nab = int(((pairs["diff"] > 0) & (pairs.p < 0.05)).sum())
axC.text(0.97, 0.04, f"{nb} of 52 signatures at P < 0.05 in the breast pool;\n{nab} above CD8A there (large filled markers)",
         transform=axC.transAxes, fontsize=FSS, color=INK, va="bottom", ha="right")

# ------------------------------------------------------------------ d  adjustment for general immune infiltration
axD = S.ax_mm(fig, 16, 146, 46, 46); S.letter_at(fig, 0, 139.5, "e"); title(axD, "Adjusted for ESTIMATE immune", -0.3)
lim = (-0.45, 0.75)
axD.plot(lim, lim, color=GREY, lw=0.5, ls=(0, (2, 2))); S.rule(axD, 0); axD.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
for _, r in L.iterrows():
    col = FAMCOL[r.family]; sig = r.p < 0.05
    axD.plot(r.est, r.est_estimate_adj, "o", ms=3.0 if sig else 2.2, mfc=col if r.p_estimate_adj < 0.05 else ("white" if sig else col),
             mec=col, mew=0.6, alpha=1 if sig else 0.45, zorder=3 if sig else 2)
s18 = L[L.p < 0.05]; ratio = (s18.est_estimate_adj / s18.est)
axD.set_xlim(*lim); axD.set_ylim(*lim)
axD.set_xticks([-0.4, 0, 0.4]); axD.set_yticks([-0.4, 0, 0.4]); axD.set_xticklabels(["−0.4", "0", "0.4"]); axD.set_yticklabels(["−0.4", "0", "0.4"])
axD.set_xlabel("Pooled interaction, unadjusted"); axD.set_ylabel("Adjusted for ESTIMATE immune × arm")
axD.text(0.03, 0.97, f"{len(s18)} signatures at P < 0.05 (large):\nestimate × {ratio.median():.2f} (median);\n"
         f"{int((s18.p_estimate_adj < 0.05).sum())} still at P < 0.05 (filled)", transform=axD.transAxes, fontsize=FSS, color=INK, va="top")
rec("e adjusted/unadjusted (18 significant)", f"median {ratio.median():.2f}, range {ratio.min():.2f}-{ratio.max():.2f}; kept P<0.05 {int((s18.p_estimate_adj < 0.05).sum())}/{len(s18)}", "E32_primary.csv")

# ------------------------------------------------------------------ e  I-SPY2 checks
Q = pd.read_csv(R / "E32_ispy2.csv")
axE = S.ax_mm(fig, 100, 146, 76, 46); S.letter_at(fig, 76, 139.5, "f"); title(axE, "I-SPY2: MammaPrint and processing", -0.3)
SC = ["NetImmune", "Immune", "CAF_ECM", "ESTIMATE", "CD8A"]
rowsE = []
for sc in SC:
    q5 = Q[(Q.question == "Q5") & (Q.score == sc)].set_index("model")
    q6 = Q[(Q.question == "Q6") & (Q.score == sc)].set_index("model")
    rowsE.append((sc, [("six-comparison pool, HR only", q5.loc["hormone receptor only (published)"], GREY, "o"),
                       ("six-comparison pool, HR + MammaPrint", q5.loc["HR + MammaPrint (E32 primary)"], INK, "o"),
                       ("I-SPY2 pembrolizumab, all 179 controls", q6.loc["all controls"], GREY, "s"),
                       ("I-SPY2 pembrolizumab, 40 Linear-adjusted controls", q6.loc["Linear-adjusted controls, z within subset"], INK, "s")]))
y = 0; YT = []
for sc, items in rowsE:
    YT.append((y + 1.5, sc))
    for k, (lab, r, col, mk) in enumerate(items):
        yy = y + k
        axE.plot([max(r.lo, -1.8), min(r.hi, 2.9)], [yy, yy], color=col, lw=0.6)
        axE.plot(r.est, yy, mk, ms=2.6, mfc=col if r.p < 0.05 else "white", mec=col, mew=0.6)
        rec(f"f {sc} {lab}", f"{r.est:+.3f} ({r.lo:.3f} to {r.hi:.3f}); P = {r.p:.3g}", "E32_ispy2.csv")
    if sc != SC[-1]:
        axE.axhline(y + 4.0, color=LIGHT, lw=0.4)
    y += 5
for yy, sc in YT:
    axE.text(-1.95, yy, {"Immune": "Immune axis", "CAF_ECM": "CAF_ECM"}.get(sc, sc), fontsize=FS, color=INK, va="center", ha="right",
             fontweight="bold")
axE.set_ylim(y - 1.4, -0.8); axE.set_yticks([]); axE.spines["left"].set_visible(False)
axE.set_xlim(-1.8, 2.9); axE.set_xticks([-1, 0, 1, 2]); axE.set_xticklabels(["−1", "0", "1", "2"]); S.rule(axE, 0)
axE.set_xlabel("Interaction per SD (positive favors checkpoint arm)")
K = S.ax_mm(fig, 76, 202, 100, 10); S.clean(K); K.set_xlim(0, 100); K.set_ylim(10, 0)
for j, (lab, col, mk) in enumerate([("six comparisons, I-SPY2 hormone receptor only", GREY, "o"),
                                     ("six comparisons, + MammaPrint (primary)", INK, "o"),
                                     ("I-SPY2 pembrolizumab, all controls", GREY, "s"),
                                     ("I-SPY2 pembrolizumab, Linear-adjusted controls only", INK, "s")]):
    xx, yy = (0, 1.5 + 3.6 * j) if j < 2 else (50, 1.5 + 3.6 * (j - 2))
    K.plot(xx + 1, yy, mk, ms=2.6, mfc=col, mec=col); K.text(xx + 2.6, yy, lab, fontsize=FSS, color=SLATE, va="center")
K.text(0, 8.9, "open, P ≥ 0.05", fontsize=FSS, color=SLATE, va="center")

S.save(fig, "Fig4")
with open(HERE / "Fig4_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched")
