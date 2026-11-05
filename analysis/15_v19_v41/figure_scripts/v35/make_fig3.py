"""v35 Figure 3 - how far the result can be trusted: random genes, sensitivity, the breast comparisons, immune infiltration.

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
fig = plt.figure(figsize=(180 * MM, 136 * MM))

# ------------------------------------------------------------------ a  random genes reach significance
COV = pd.read_csv(R / "E32_coverage.csv")
SIX = ["IMM150", "ISPY2_pembro", "ISPY2_durva", "CM025", "JR101", "NeoTRIP"]
size = COV[COV.scored & COV.comparison.isin(SIX)].groupby("signature").present.mean()
L["size"] = L.signature.map(size)
axA = S.ax_mm(fig, 14, 10, 50, 48); S.letter_at(fig, 0, 3.5, "a"); title(axA, "Random genes reach significance", -0.2)
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
        ("breast three", "Breast comparisons only", above("breast three"))]
axB = S.ax_mm(fig, 120, 10, 32, 48); S.letter_at(fig, 74, 3.5, "b"); title(axB, "Sensitivity", -1.38)
for j, (key, lab, ab) in enumerate(SETS):
    if key == "native":
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
    rec(f"b {lab}", f"favor {fav}, against {ag}; above CD8A {ab}", "E32_pooled_all_sets.csv / E32_native.csv / E33_pooled.csv / E33_paired.csv")
axB.set_ylim(len(SETS) - 0.4, -0.6); axB.set_yticks([]); axB.spines["left"].set_visible(False)
axB.set_xlim(-4, 27); axB.set_xticks([0, 10, 20]); axB.axvline(0, color=INK, lw=0.5)
axB.set_xlabel("Signatures at P < 0.05 (of 52)")
axB.text(25, -0.95, "favoring / against", fontsize=FSS, color=SLATE, ha="right", va="center")
axB.text(30.5, -1.05, "above\nCD8A", fontsize=FSS, color=SLATE, ha="center", va="bottom", linespacing=0.9)

# ------------------------------------------------------------------ c  breast against six
axC = S.ax_mm(fig, 14, 78, 62, 50); S.letter_at(fig, 0, 71.5, "c"); title(axC, "The three breast comparisons", -0.2)
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
        rec(f"c above CD8A in breast: {s}", f"breast {br.loc[s, 'est']:+.3f}; minus CD8A {pairs.loc[s, 'diff']:+.3f} (P = {pairs.loc[s, 'p']:.3g})", "E33_pooled.csv; E33_paired.csv")
for s, lab, dxy in [("CD8A", "CD8A", (5, 5)), ("NetImmune", "NetImmune", (5, 2)), ("CAF_ECM", "CAF_ECM", (5, -2)),
                    ("ESTIMATE", "ESTIMATE", (5, -3))]:
    x = P.set_index("signature").loc[s, "est"]; y = BX.loc[s, "est"]
    axC.plot(x, y, "D", ms=3.6, mfc=IMMUNE if s == "CD8A" else INK, mec="white", mew=0.4, zorder=4)
    axC.annotate(lab, (x, y), xytext=dxy, textcoords="offset points", fontsize=FSS, color=IMMUNE if s == "CD8A" else INK, va="center")
    rec(f"c {s}", f"six {x:+.3f}; breast {y:+.3f} ({BX.loc[s, 'lo']:.3f} to {BX.loc[s, 'hi']:.3f}); P = {BX.loc[s, 'p']:.3g}", "E32_primary.csv; E33_breast_extras.csv")
axC.set_xlim(-0.3, 0.38); axC.set_ylim(*lim)
axC.set_xticks([-0.2, 0, 0.2]); axC.set_xticklabels(["−0.2", "0", "0.2"])
axC.set_yticks([-0.8, -0.4, 0, 0.4, 0.8]); axC.set_yticklabels(["−0.8", "−0.4", "0", "0.4", "0.8"])
axC.set_xlabel("Pooled interaction, six comparisons"); axC.set_ylabel("Pooled interaction,\nthree breast comparisons")
nb = int(((br.p < 0.05) & (br.est > 0)).sum()); nab = int(((pairs["diff"] > 0) & (pairs.p < 0.05)).sum())
axC.text(0.03, 0.97, f"{nb} of 52 signatures at P < 0.05 in the breast pool;\n{nab} above CD8A there (large filled markers)",
         transform=axC.transAxes, fontsize=FSS, color=INK, va="top")

# ------------------------------------------------------------------ d  adjustment for general immune infiltration
axD = S.ax_mm(fig, 108, 78, 50, 50); S.letter_at(fig, 86, 71.5, "d"); title(axD, "Adjusted for ESTIMATE immune", -0.35)
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
rec("d adjusted/unadjusted (18 significant)", f"median {ratio.median():.2f}, range {ratio.min():.2f}-{ratio.max():.2f}; kept P<0.05 {int((s18.p_estimate_adj < 0.05).sum())}/{len(s18)}", "E32_primary.csv")

S.save(fig, "Fig3")
with open(HERE / "Fig3_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched")
