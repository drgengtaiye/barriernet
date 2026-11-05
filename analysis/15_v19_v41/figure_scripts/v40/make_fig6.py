"""v40 Figure 6 - barriernet: opposed components and a frozen predictive version (E32, E33, E37).
(a) NetImmune and CAF_ECM interaction per randomized comparison and pooled; (b) development of barriernet 0.3 in 20
single-arm cohorts (leave-one-cohort-out); (c) frozen module weights; (d) barriernet 0.3, NetImmune 0.2.1 and CD8A per
comparison and pooled; (e) barriernet 0.3 against its size-matched random gene sets; (f) paired difference from each
panel gene (six and breast comparisons; added after the test results). Writes only into figures/cns/v40/.
"""
import sys, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/repo_barriernet/barriernet")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, MACRO, MM, MODCOL
from barriernet.calibrate import pool

HERE = Path(__file__).resolve().parent
S.OUT = HERE
GUARD = {p: hashlib.md5(open(p, "rb").read()).hexdigest() for p in glob.glob(str(HERE.parent / "Fig*.p*"))}
R = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
FS, FSS = 6.0, 5.6
NUM = {}
rec = lambda k, v, s: NUM.__setitem__(k, (v, s))
LABEL = {"IMM150": "IMmotion150", "ISPY2_pembro": "I-SPY2 pembrolizumab", "ISPY2_durva": "I-SPY2 durvalumab–olaparib",
         "CM025": "CheckMate 025", "JR101": "JAVELIN Renal 101", "NeoTRIP": "NeoTRIP", "NeoPembrOV": "NeoPembrOv"}
SIX = ["IMM150", "CM025", "JR101", "ISPY2_pembro", "ISPY2_durva", "NeoTRIP"]
RENAL, BREAST = SIX[:3], SIX[3:]


def head(x_mm, y_mm, letter, s):
    """Panel letter and bold title on one line, in mm from the top-left corner."""
    S.letter_at(fig, x_mm, y_mm, letter)
    W, H = fig.get_size_inches()
    fig.text((x_mm + 4.2) * MM / W, 1 - (y_mm + 0.6) * MM / H, s, fontsize=6.8, fontweight="bold", color=INK, va="top", ha="left")


def fmtx(ax, lo, hi, ticks):
    ax.set_xlim(lo, hi); ax.set_xticks(ticks); ax.set_xticklabels([f"{t:g}".replace("-", "−") for t in ticks]); S.rule(ax, 0)


fig = plt.figure(figsize=(180 * MM, 206 * MM))

# ------------------------------------------------------------------ a  NetImmune and CAF_ECM per comparison
PC = pd.read_csv(R / "E32_per_comparison.csv"); PC = PC[PC.variant == "primary"]
E32P = pd.read_csv(R / "E32_primary.csv").set_index("signature"); BX = pd.read_csv(R / "E33_breast_extras.csv").set_index("name")
axA = S.ax_mm(fig, 30, 12, 56, 58); head(0, 3, "a", "Within barriernet: opposed components")
rows = [(LABEL[k], k) for k in SIX] + [("Renal three", "renal"), ("Breast three", "breast"), ("Six comparisons", "six")]
for i, (lab, k) in enumerate(rows):
    for j, (sc, col) in enumerate([("NetImmune", IMMUNE), ("CAF_ECM", STROMA)]):
        if k in SIX:
            r = PC[(PC.signature == sc) & (PC.comparison == k)].iloc[0]; est, lo, hi, p = r.y, r.y - 1.96 * r.se, r.y + 1.96 * r.se, r.p
        elif k == "six":
            r = E32P.loc[sc]; est, lo, hi, p = r.est, r.lo, r.hi, r.p
        elif k == "breast":
            r = BX.loc[sc]; est, lo, hi, p = r.est, r.lo, r.hi, r.p
        else:
            q = PC[(PC.signature == sc) & PC.comparison.isin(RENAL)]; q = pool(q.y.values, q.se.values); est, lo, hi, p = q["est"], q["lo"], q["hi"], q["p"]
        y = i + (j - 0.5) * 0.36; mk = "D" if k in ("renal", "breast", "six") else "o"
        axA.plot([max(lo, -2.25), min(hi, 2.4)], [y, y], color=col, lw=0.7)
        axA.plot(est, y, mk, ms=3.0 if mk == "o" else 3.2, mfc=col if p < 0.05 else "white", mec=col, mew=0.7)
        rec(f"a {sc} {lab}", f"{est:+.3f} ({lo:+.3f} to {hi:+.3f}); P = {p:.3g}", "E32/E33")
axA.axhline(len(SIX) - 0.5, color=LIGHT, lw=0.5)
axA.set_yticks(range(len(rows))); axA.set_yticklabels([r[0] for r in rows], fontsize=FSS); axA.invert_yaxis()
fmtx(axA, -2.25, 2.4, [-2, -1, 0, 1, 2]); axA.set_xlabel("Interaction per SD (positive favors checkpoint arm)")
for j, (sc, col) in enumerate([("NetImmune", IMMUNE), ("CAF_ECM", STROMA)]):
    axA.plot(1.05, 0.3 + 0.75 * j, "o", ms=2.8, color=col); axA.text(1.2, 0.3 + 0.75 * j, sc, fontsize=FSS, color=SLATE, va="center")
axA.text(2.4, 1.85, "filled, P < 0.05", fontsize=5.0, color=SLATE, ha="right", va="center")

# ------------------------------------------------------------------ b  development (single-arm cohorts only)
LO = pd.read_csv(R / "E37_loco_pooled.csv")
axB = S.ax_mm(fig, 128, 12, 48, 30); head(96, 3, "b", "barriernet 0.3: development (single-arm only)")
order = ["C0", "C1", "C2", "C4", "C3"]
lab = {"C0": "NetImmune 0.2.1", "C1": "Immune axis 0.2.1", "C2": "Module weights", "C4": "Gene weights, ridge",
       "C3": "Gene weights, elastic net"}
for i, c in enumerate(order):
    r = LO[LO.candidate == c].iloc[0]; col = IMMUNE if c == "C2" else SLATE
    axB.plot([r.lo, r.hi], [i, i], color=col, lw=0.7); axB.plot(r.est, i, "o", ms=3.2, mfc=col if c == "C2" else "white", mec=col, mew=0.7)
    axB.text(0.86, i, f"{r.est:+.2f}", fontsize=FSS, color=SLATE, va="center")
    rec(f"b {lab[c]}", f"{r.est:+.3f} ({r.lo:+.3f} to {r.hi:+.3f}); k = {int(r.k)}", "E37_loco_pooled.csv")
axB.set_yticks(range(5)); axB.set_yticklabels([lab[c] for c in order], fontsize=FSS); axB.invert_yaxis()
axB.set_xlim(0, 0.95); axB.set_xticks([0, 0.25, 0.5, 0.75]); axB.set_xticklabels(["0", "0.25", "0.5", "0.75"])
axB.set_xlabel("Held-out log-odds of response per SD\n(20 single-arm cohorts; NetImmune 19)")

# ------------------------------------------------------------------ c  frozen weights
FW = pd.read_csv(R / "E37_frozen_weights.tsv", sep="\t"); tot = FW.groupby("module").weight.sum()
axC = S.ax_mm(fig, 128, 60, 48, 20); head(96, 53, "c", "Frozen module weights")
mods = ["IFN", "APM", "TLS_B", "T_CELL", "CAF_ECM", "TAM_SUP"]
for i, m in enumerate(mods):
    axC.barh(i, tot[m], height=0.62, color=MODCOL[m]); axC.text(tot[m] + (0.012 if tot[m] > 0 else -0.012), i, f"{tot[m]:+.2f}",
                                                             fontsize=5.2, color=SLATE, va="center", ha="left" if tot[m] > 0 else "right")
    rec(f"c weight {m}", f"{tot[m]:+.4f}", "E37_frozen_weights.tsv")
axC.set_yticks(range(6)); axC.set_yticklabels(mods, fontsize=5.2); axC.invert_yaxis(); axC.axvline(0, color=INK, lw=0.5)
axC.set_xlim(-0.25, 0.45); axC.set_xticks([-0.2, 0, 0.2, 0.4]); axC.set_xticklabels(["−0.2", "0", "0.2", "0.4"])
axC.set_xlabel("Sum of gene weights (|w| sums to 1)")

# ------------------------------------------------------------------ d  test in the randomized comparisons
T7 = pd.read_csv(R / "E37_test_per_comparison.csv"); TP = pd.read_csv(R / "E37_test_pooled.csv")
axD = S.ax_mm(fig, 30, 94, 62, 64); head(0, 86, "d", "Frozen, then tested in randomized comparisons")
NETC = "#6FA6DC"
SCS = [("BN03", "barriernet 0.3", IMMUNE, "o", True), ("NetImmune", "NetImmune 0.2.1", NETC, "^", True), ("CD8A", "CD8A", INK, "s", True)]
PSET = [("Renal three", "renal three"), ("Breast three", "breast three"), ("Six comparisons", "six (primary)"),
        ("Seven (+ NeoPembrOv)", "seven (with NeoPembrOv)"), ("Chemotherapy backbone (four)", "four chemotherapy-backbone")]
rows = [(LABEL[k], k) for k in SIX + ["NeoPembrOV"]] + PSET
PLAB = {"BN03": "barriernet 0.3", "NetImmune": "NetImmune 0.2.1", "CD8A": "CD8A"}
for i, (lab_, k) in enumerate(rows):
    for j, (sc, nm, col, mk, filled) in enumerate(SCS):
        if isinstance(k, str) and k in LABEL:
            r = T7[(T7.comparison == k) & (T7.score == sc)].iloc[0]; est, lo, hi = r.y, r.y - 1.96 * r.se, r.y + 1.96 * r.se
            p = 2 * (1 - __import__("scipy.stats", fromlist=["norm"]).norm.cdf(abs(r.y / r.se)))
        else:
            r = TP[(TP.set == k) & (TP.score == PLAB[sc])].iloc[0]; est, lo, hi, p = r.est, r.lo, r.hi, r.p
        y = i + (j - 1) * 0.27; mk2 = "D" if not (isinstance(k, str) and k in LABEL) else mk
        axD.plot([max(lo, -1.6), min(hi, 2.6)], [y, y], color=col, lw=0.6)
        axD.plot(est, y, mk2, ms=2.7, mfc=col if p < 0.05 else "white", mec=col, mew=0.6)
        rec(f"d {nm} {lab_}", f"{est:+.3f} ({lo:+.3f} to {hi:+.3f}); P = {p:.3g}", "E37_test_*.csv")
axD.axhline(6.5, color=LIGHT, lw=0.5); axD.axhline(len(SIX) + 0.5, color=LIGHT, lw=0.5)
axD.set_yticks(range(len(rows))); axD.set_yticklabels([r[0] for r in rows], fontsize=FSS); axD.invert_yaxis()
fmtx(axD, -1.6, 2.6, [-1, 0, 1, 2]); axD.set_xlabel("Interaction per SD (positive favors checkpoint arm)")
for j, (sc, nm, col, mk, filled) in enumerate(SCS):
    axD.plot(1.35, -0.1 + 0.75 * j, mk, ms=2.7, mfc=col, mec=col, mew=0.6)
    axD.text(1.5, -0.1 + 0.75 * j, nm, fontsize=5.2, color=SLATE, va="center")
axD.text(1.35, 2.2, "filled, P < 0.05; diamonds, pooled", fontsize=5.0, color=SLATE, va="center")
six = TP[(TP.set == "six (primary)") & (TP.score == "barriernet 0.3")].iloc[0]

# ------------------------------------------------------------------ e  random reference
RD = pd.read_csv(R / "E37_reference_draws.csv")
axE = S.ax_mm(fig, 30, 178, 62, 18); head(0, 170, "e", "Specific to its genes (six comparisons)")
d6 = RD["six (primary)"].dropna().values
axE.hist(d6, bins=40, color=LIGHT, edgecolor="white", lw=0.3)
p95 = np.percentile(d6, 95); axE.axvline(p95, color=GREY, lw=0.6, ls=(0, (2, 2))); axE.axvline(six.est, color=IMMUNE, lw=1.1)
axE.text(six.est + 0.01, axE.get_ylim()[1] * 0.92, f"barriernet 0.3 {six.est:+.2f}\nempirical P = {six.emp_p:.3f}", fontsize=5.2, color=INK, va="top")
axE.text(p95 - 0.01, axE.get_ylim()[1] * 0.92, "95th\npercentile", fontsize=5.0, color=SLATE, va="top", ha="right")
axE.set_yticks([]); axE.spines["left"].set_visible(False)
axE.set_xlabel("Pooled interaction of 1,000 random gene sets (same size and weights)")
rec("e reference median / p95 / emp P", f"{np.median(d6):+.3f} / {p95:+.3f} / {six.emp_p:.4f}", "E37_reference_draws.csv")

# ------------------------------------------------------------------ f  single-gene panel
PX = pd.read_csv(R / "E37_panel_extra.csv")
g6 = PX[PX.set == "six (primary)"].sort_values("gene_est", ascending=False); genes = list(g6.gene)
axF1 = S.ax_mm(fig, 124, 100, 25, 92); axF2 = S.ax_mm(fig, 152, 100, 25, 92)
head(104, 86, "f", "Against each single gene (added after the test)")
for ax, setname, hd in [(axF1, "six (primary)", "Six"), (axF2, "breast three", "Breast three")]:
    t = PX[PX.set == setname].set_index("gene")
    for i, g in enumerate(genes):
        if g not in t.index:
            ax.text(0, i, "n.m.", fontsize=4.6, color=GREY, ha="center", va="center"); continue
        r = t.loc[g]; col = INK if g == "CD8A" else IMMUNE
        ax.plot([max(r.lo, -0.6), min(r.hi, 1.2)], [i, i], color=col, lw=0.55)
        ax.plot(r["diff"], i, "o", ms=2.3, mfc=col if r.p < 0.05 else "white", mec=col, mew=0.55)
        rec(f"f {hd} {g}", f"{r['diff']:+.3f} ({r.lo:+.3f} to {r.hi:+.3f}); P = {r.p:.3g}", "E37_panel_extra.csv")
    ax.set_ylim(len(genes) - 0.4, -0.6); fmtx(ax, -0.6, 1.2, [-0.5, 0, 0.5, 1])
    n_up = int(((t["diff"] > 0) & (t.p < 0.05)).sum()); n_dn = int(((t["diff"] < 0) & (t.p < 0.05)).sum())
    ax.text(0.5, 1.012, f"{hd}\nabove {n_up}, below {n_dn} of {len(t)}", transform=ax.transAxes, fontsize=5.2, color=SLATE, ha="center", va="bottom")
axF1.set_yticks(range(len(genes))); axF1.set_yticklabels(genes, fontsize=5.0); axF2.set_yticks([])
axF2.spines["left"].set_visible(False)
fig.text(150.5 * MM / (180 * MM), 1 - 203 * MM / (206 * MM), "barriernet 0.3 minus gene, per SD (filled, P < 0.05)", fontsize=FS, color=INK, ha="center", va="bottom")

S.save(fig, "Fig6")
with open(HERE / "Fig6_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched")
