"""v34 Figure 2 - published immune signatures against random genes and a single gene in six randomized comparisons.

(a) The six comparisons: patients per arm, cancer type, regimen, endpoint.
(b) Pooled score x treatment interaction (benefit-oriented, per SD; random effects, manuscript rule) for the 52 published
    signatures and six barriernet / comparator scores, ranked; grey bar, each score's size-matched random reference
    (median to 95th percentile in its declared direction); filled, beyond it; band, CD8A (95% CI).
(c) The same interaction in each comparison (dot, nominal P < 0.05).
(d) Paired difference from CD8A (300 within-arm bootstrap resamples per comparison, inverse-variance pooled).
(e) Counts over the 52 signatures.
Every number is read from analysis/results/E32_* (protocol analysis/E32_prereg.md). Writes only into figures/cns/v34/.
"""
import sys, re, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, PAPER, IMMUNE, STROMA, FAV, UNFAV, MM, DIV

HERE = Path(__file__).resolve().parent
S.OUT = HERE                                                    # never write into figures/cns
GUARD = {p: hashlib.md5(open(p, "rb").read()).hexdigest() for p in glob.glob(str(HERE.parent / "Fig*.p*"))}
R = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
FAMCOL = {"Inflamed/IFN": IMMUNE, "Other": "#8E9AA6", "Confounder": STROMA, "barriernet / comparator": INK}
REFC = "#C5CDD5"
NUM = {}
rec = lambda k, v, s: NUM.__setitem__(k, (v, s))


def nice(s):
    lab = {"TMEscore_z": "TMEscore", "Immune": "Immune axis", "CAF_ECM": "CAF_ECM (fibroblast)",
           "NetImmune": "NetImmune", "ESTIMATE": "ESTIMATE", "CD8A": "CD8A alone"}
    if s in lab:
        return lab[s]
    m = re.match(r"^(.*)_([A-Z][A-Za-z\-]+)$", s)
    return f"{m.group(1).replace('_', ' ')} ({m.group(2)})" if m else s.replace("_", " ")


P = pd.read_csv(R / "E32_primary.csv")
EXTRA = ["CD8A", "ESTIMATE", "TMEscore_z", "NetImmune", "Immune", "CAF_ECM"]
D = P[(P.kind == "library") | P.signature.isin(EXTRA)].copy().sort_values("est", ascending=False).reset_index(drop=True)
L = D[D.kind == "library"]
cd8 = P[P.signature == "CD8A"].iloc[0]
PER = pd.read_csv(R / "E32_per_comparison.csv"); PER = PER[PER.variant == "primary"]
COMP = [("IMM150", "IMmotion150"), ("CM025", "CheckMate 025"), ("JR101", "JAVELIN Renal 101"),
        ("ISPY2_pembro", "I-SPY2 pembrolizumab"), ("ISPY2_durva", "I-SPY2 durvalumab–olaparib"), ("NeoTRIP", "NeoTRIP")]
rec("CD8A pooled", f"{cd8.est:+.3f} ({cd8.lo:.3f} to {cd8.hi:.3f}); P = {cd8.p:.3g}", "E32_primary.csv")

fig = plt.figure(figsize=(180 * MM, 204 * MM))

# ------------------------------------------------------------------ a  the six comparisons
TRIALS = [("renal", "IMmotion150", "atezolizumab ± bevacizumab", "sunitinib", 174, 89),
          ("renal", "CheckMate 025", "nivolumab", "everolimus", 120, 130),
          ("renal", "JAVELIN Renal 101", "avelumab + axitinib", "sunitinib", 354, 372),
          ("breast", "I-SPY2", "pembrolizumab + paclitaxel", "paclitaxel", 69, 179),
          ("breast", "I-SPY2", "durvalumab + olaparib + paclitaxel", "paclitaxel", 71, 29),
          ("breast", "NeoTRIP", "atezolizumab + carbo/nab-paclitaxel", "carbo/nab-paclitaxel", 119, 122)]
axA = S.ax_mm(fig, 4, 9, 118, 36); S.clean(axA); axA.set_xlim(0, 118); axA.set_ylim(36, 0)
S.letter_at(fig, 0, 3.0, "a")
axA.text(0, -1.6, "Six randomized comparisons with expression in both arms", fontsize=6.5, fontweight="bold", color=INK, va="bottom")
XC, SC_ = 96.0, 19.0 / 372                                   # butterfly centre (mm) and mm per patient
for x, lab, ha in [(0, "Trial", "left"), (22, "Checkpoint arm", "left"), (57, "Control arm", "left"),
                   (XC - 1.2, "Control", "right"), (XC + 1.2, "Checkpoint", "left")]:
    axA.text(x, 1.2, lab, fontsize=5.0, color=SLATE, ha=ha, va="center")
axA.plot([0, 118], [2.6, 2.6], color=INK, lw=0.5)
yrow, last = 3.2, None
for ca, nm, rx1, rx0, n1, n0 in TRIALS:
    if ca != last:
        yrow += 1.6
        axA.text(0, yrow + 0.6, {"renal": "Renal cell carcinoma · progression-free survival",
                                 "breast": "Breast cancer · pathological complete response"}[ca],
                 fontsize=5.0, fontweight="bold", color=INK, va="center")
        yrow += 2.9; last = ca
    yc = yrow + 1.3
    axA.text(0, yc, nm, fontsize=5.2, color=INK, va="center")
    axA.text(22, yc, rx1, fontsize=4.8, color=INK, va="center")
    axA.text(57, yc, rx0, fontsize=4.8, color=SLATE, va="center")
    axA.add_patch(plt.Rectangle((XC - n0 * SC_, yc - 0.9), n0 * SC_, 1.8, fc=GREY, ec="none"))
    axA.add_patch(plt.Rectangle((XC, yc - 0.9), n1 * SC_, 1.8, fc=INK, ec="none"))
    axA.text(XC - n0 * SC_ - 0.8, yc, str(n0), fontsize=4.8, color=INK, va="center", ha="right")
    axA.text(XC + n1 * SC_ + 0.8, yc, str(n1), fontsize=4.8, color=INK, va="center", ha="left")
    yrow += 3.3
axA.plot([XC, XC], [3.0, yrow - 0.4], color=INK, lw=0.4)
axA.plot([0, 118], [yrow, yrow], color=INK, lw=0.5)
rec("a patients", f"{sum(t[4] + t[5] for t in TRIALS)}", "E21/E10/E5 trial tables (counts)")

# ------------------------------------------------------------------ e  counts (top right)
axE = S.ax_mm(fig, 146, 12, 30, 26)
S.letter_at(fig, 124, 3.0, "e")
axE.text(-0.86, 1.13, "Of 52 published signatures", transform=axE.transAxes, fontsize=6.5, fontweight="bold", color=INK, va="bottom")
nfav = int(((L.p < 0.05) & (L.est > 0)).sum()); nag = int(((L.p < 0.05) & (L.est < 0)).sum())
nbey = int((L.emp_p_declared < 0.05).sum())
nabove = int(((L.diff_vs_CD8A > 0) & (L.cd8a_p < 0.05)).sum()); nbelow = int(((L.diff_vs_CD8A < 0) & (L.cd8a_p < 0.05)).sum())
ROWS = [("favor checkpoint arm, P < 0.05", nfav, FAV), ("against it, P < 0.05", nag, UNFAV),
        ("beyond their random genes", nbey, INK), ("above CD8A", nabove, IMMUNE), ("below CD8A", nbelow, GREY)]
for j, (lab, v, col) in enumerate(ROWS):
    axE.barh(j, v, height=0.62, color=col, lw=0)
    axE.text(v + 0.8, j, str(v), fontsize=5.4, color=INK, va="center", fontweight="bold" if lab == "above CD8A" else "normal")
    axE.text(-1.2, j, lab, fontsize=5.0, color=INK, va="center", ha="right")
    rec(f"e {lab}", f"{v}/52", "E32_primary.csv")
axE.set_ylim(len(ROWS) - 0.4, -0.6); axE.set_yticks([]); axE.spines["left"].set_visible(False)
axE.set_xlim(0, 30); axE.set_xticks([0, 10, 20, 30]); axE.axvline(0, color=INK, lw=0.5)
axE.set_xlabel("Signatures")

# ------------------------------------------------------------------ b  ranked pooled interaction
TOP, H = 52.0, 136.0
n = len(D); XL = (-0.55, 0.62)
axB = S.ax_mm(fig, 38, TOP, 50, H)
S.letter_at(fig, 0, TOP - 7.0, "b")
axB.text(-0.70, 1.035, "Pooled over six comparisons", transform=axB.transAxes, fontsize=6.5, fontweight="bold", color=INK, va="bottom")
axB.axvspan(cd8.lo, cd8.hi, color="#EEF4FB", zorder=0, lw=0)
axB.axvline(cd8.est, color=IMMUNE, lw=0.6, zorder=1)
S.rule(axB, 0)
for i, r in D.iterrows():
    col = FAMCOL[r.family]
    if np.isfinite(r.ref_median) and np.isfinite(r.ref_p95_declared):
        a_, b_ = sorted([r.ref_median, r.ref_p95_declared])
        axB.plot([a_, b_], [i, i], color=REFC, lw=2.4, solid_capstyle="butt", zorder=1)
        axB.plot([r.ref_p95_declared] * 2, [i - 0.32, i + 0.32], color=GREY, lw=0.6, zorder=1)
    axB.plot([max(r.lo, XL[0]), min(r.hi, XL[1])], [i, i], color=col, lw=0.6, zorder=2)
    beyond = np.isfinite(r.emp_p_declared) and r.emp_p_declared < 0.05
    if r.signature == "CD8A":
        axB.plot(r.est, i, "D", ms=2.5, mfc=IMMUNE, mec=IMMUNE, zorder=3)
    else:
        mk = "D" if r.kind == "extra" else "o"
        axB.plot(r.est, i, mk, ms=2.5, mfc=col if beyond else "white", mec=col, mew=0.6, zorder=3)
    lab = nice(r.signature) + (" †" if r.declared_dir < 0 and r.kind == "library" else "")
    axB.text(-0.015, i, lab, transform=axB.get_yaxis_transform(), fontsize=4.5, ha="right", va="center",
             color=IMMUNE if r.signature == "CD8A" else INK, fontweight="bold" if r.kind == "extra" else "normal")
    rec(f"b {r.signature}", f"{r.est:+.3f} ({r.lo:.3f} to {r.hi:.3f}); P = {r.p:.3g}; random-reference P = {r.emp_p_declared:.3g}", "E32_primary.csv")
axB.set_ylim(n - 0.4, -0.7); axB.set_yticks([]); axB.spines["left"].set_visible(False)
axB.set_xlim(*XL); axB.set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6]); axB.set_xticklabels(["−0.4", "−0.2", "0", "0.2", "0.4", "0.6"])
axB.set_xlabel("Interaction per SD\n(positive favors checkpoint arm)")
axB.text(cd8.hi + 0.012, n - 1.2, "CD8A\n95% CI", fontsize=4.8, color=IMMUNE, ha="left", va="bottom")

# ------------------------------------------------------------------ c  each comparison (heat map, rows aligned with b)
axC = S.ax_mm(fig, 92, TOP, 30, H)
S.letter_at(fig, 89, TOP - 7.0, "c")
axC.text(0, 1.035, "Each comparison", transform=axC.transAxes, fontsize=6.5, fontweight="bold", color=INK, va="bottom")
norm = Normalize(-1.0, 1.0)
M = np.full((n, 6), np.nan); Pv = np.full((n, 6), np.nan)
for i, r in D.iterrows():
    for j, (k, _) in enumerate(COMP):
        q = PER[(PER.signature == r.signature) & (PER.comparison == k)]
        if len(q):
            M[i, j], Pv[i, j] = q.y.iloc[0], q.p.iloc[0]
cmap = DIV.reversed().copy(); cmap.set_bad("#E3E7EB")
axC.imshow(np.ma.masked_invalid(np.clip(M, -1, 1)), cmap=cmap, norm=norm, aspect="auto", interpolation="nearest",
           extent=(-0.5, 5.5, n - 0.5, -0.5))
for j in range(7):
    axC.axvline(j - 0.5, color="white", lw=0.6)
yy, xx = np.where(Pv < 0.05)
axC.scatter(xx, yy, s=1.6, c=INK, lw=0, zorder=3)
axC.set_xlim(-0.5, 5.5); axC.set_ylim(n - 0.4, -0.7); axC.set_yticks([])
axC.set_xticks(range(6)); axC.set_xticklabels(["IMmotion150", "CheckMate 025", "JAVELIN 101", "I-SPY2 pembro",
                                                 "I-SPY2 durva", "NeoTRIP"], rotation=90, fontsize=4.8)
axC.tick_params(axis="x", length=0, pad=1.5)
for s_ in axC.spines.values():
    s_.set_visible(False)
nanc = int(np.isnan(M).sum()); rec("c cells not scored", str(nanc), "E32_per_comparison.csv")

# ------------------------------------------------------------------ d  paired difference from CD8A
axD = S.ax_mm(fig, 128, TOP, 24, H)
S.letter_at(fig, 124, TOP - 7.0, "d")
axD.text(0, 1.035, "Against CD8A", transform=axD.transAxes, fontsize=6.5, fontweight="bold", color=INK, va="bottom")
S.rule(axD, 0)
for i, r in D.iterrows():
    if r.signature == "CD8A" or not np.isfinite(r.diff_vs_CD8A):
        continue
    col = FAMCOL[r.family]; sig = r.cd8a_p < 0.05
    axD.plot([max(r.cd8a_lo, -0.75), min(r.cd8a_hi, 0.35)], [i, i], color=col, lw=0.6, zorder=2)
    axD.plot(r.diff_vs_CD8A, i, "D" if r.kind == "extra" else "o", ms=2.4, mfc=col if sig else "white", mec=col, mew=0.6, zorder=3)
    rec(f"d {r.signature}", f"{r.diff_vs_CD8A:+.3f} ({r.cd8a_lo:.3f} to {r.cd8a_hi:.3f}); P = {r.cd8a_p:.3g}", "E32_primary.csv")
axD.set_ylim(n - 0.4, -0.7); axD.set_yticks([]); axD.spines["left"].set_visible(False)
axD.set_xlim(-0.75, 0.35); axD.set_xticks([-0.6, -0.3, 0, 0.3]); axD.set_xticklabels(["−0.6", "−0.3", "0", "0.3"])
axD.set_xlabel("Signature − CD8A\n(interaction per SD)")

# ------------------------------------------------------------------ keys (right column under e, and bottom)
K = S.ax_mm(fig, 156, TOP + 2, 22, 60); S.clean(K); K.set_xlim(0, 22); K.set_ylim(60, 0)
K.text(0, 1.0, "Family", fontsize=5.4, fontweight="bold", color=INK, va="center")
for j, (lab, col) in enumerate([("Inflamed/IFN", FAMCOL["Inflamed/IFN"]), ("Other", FAMCOL["Other"]),
                                ("Confounder", STROMA), ("barriernet /\ncomparator", INK)]):
    K.plot(1, 4.6 + 4.2 * j, "o" if j < 3 else "D", ms=2.6, mfc=col, mec=col); K.text(2.6, 4.6 + 4.2 * j, lab, fontsize=4.8, color=INK, va="center")
y0 = 23.5
K.text(0, y0, "b, d", fontsize=5.4, fontweight="bold", color=INK, va="center")
K.plot(1, y0 + 3.8, "o", ms=2.6, mfc=SLATE, mec=SLATE); K.text(2.6, y0 + 3.8, "beyond random (b)\nbelow CD8A (d)", fontsize=4.6, color=SLATE, va="center")
K.plot(1, y0 + 8.6, "o", ms=2.6, mfc="white", mec=SLATE, mew=0.6); K.text(2.6, y0 + 8.6, "not", fontsize=4.6, color=SLATE, va="center")
K.plot([0, 2.2], [y0 + 12.4] * 2, color=REFC, lw=2.4); K.text(2.6, y0 + 12.4, "random reference,\nmedian to 95th pct", fontsize=4.6, color=SLATE, va="center")
K.text(0, y0 + 17.2, "† declared resistance", fontsize=4.6, color=SLATE, va="center")
cb = S.ax_mm(fig, 157, TOP + 66, 18, 2.2)
cb.imshow(np.linspace(-1, 1, 200)[None, :], cmap=DIV.reversed(), norm=norm, aspect="auto", extent=(-1, 1, 0, 1))
cb.set_yticks([]); cb.set_xticks([-1, 0, 1]); cb.set_xticklabels(["≤ −1", "0", "≥ 1"], fontsize=4.6)
for s_ in cb.spines.values():
    s_.set_visible(False)
cb.tick_params(axis="x", length=1.5, pad=1)
cb.text(0, 1.9, "c: interaction per SD", transform=cb.transAxes, fontsize=4.8, color=INK, va="bottom")
cb.text(0, -2.6, "dot: P < 0.05", transform=cb.transAxes, fontsize=4.6, color=SLATE, va="top")
nb = S.ax_mm(fig, 157, TOP + 77, 3, 2.2); nb.add_patch(plt.Rectangle((0, 0), 1, 1, fc="#E3E7EB", ec="none")); S.clean(nb)
nb.set_xlim(0, 1); nb.set_ylim(0, 1)
nb.text(1.4, 0.5, "not scored\n(too few genes)", transform=nb.transAxes, fontsize=4.6, color=SLATE, va="center")

S.save(fig, "Fig2")
with open(HERE / "Fig2_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched")
