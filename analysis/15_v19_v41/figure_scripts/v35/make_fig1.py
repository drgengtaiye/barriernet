"""v35 Figure 1 - published immune signatures against random genes and CD8A in six randomized comparisons (print-size
version of v34 Fig. 2; the labelled 52-row version is Supplementary Fig. S30).

(a) Pooled score x treatment interaction (benefit-oriented, per SD) of the 52 signatures by family and of five composite
    scores; filled, beyond the size-matched random reference; ring, P < 0.05; band, CD8A (95% CI).
(b) Counts over the 52 signatures.
(c) Paired difference from CD8A for each signature, ranked, with 95% CI; band, the difference detectable with 80% power.
(d) Each comparison: the 52 signatures against single-comparison random sets (5th-95th percentile); line, CD8A.
Numbers from analysis/results/E32_* and E33_* (protocols E32_prereg.md, E33_prereg.md). Writes only into figures/cns/v35/.
"""
import sys, re, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, PAPER, IMMUNE, STROMA, FAV, UNFAV, MM

HERE = Path(__file__).resolve().parent
S.OUT = HERE
GUARD = {p: hashlib.md5(open(p, "rb").read()).hexdigest() for p in glob.glob(str(HERE.parent / "Fig*.p*"))}
R = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
OTHER = "#8E9AA6"
FAMCOL = {"Inflamed/IFN": IMMUNE, "Other": OTHER, "Confounder": STROMA, "barriernet / comparator": INK}
FS, FSS = 6.0, 5.6
NUM = {}
rec = lambda k, v, s: NUM.__setitem__(k, (v, s))


def title(ax, s, x=0.0):
    ax.text(x, 1.035, s, transform=ax.transAxes, fontsize=6.8, fontweight="bold", color=INK, va="bottom", ha="left")


def nice(s):
    lab = {"TMEscore_z": "TMEscore", "Immune": "Immune axis", "CAF_ECM": "CAF_ECM", "NetImmune": "NetImmune",
           "ESTIMATE": "ESTIMATE", "CD39-CD8Tcell_Chow": "CD39 CD8 T cell (Chow)", "Teff_McDermott": "Teff (McDermott)"}
    if s in lab:
        return lab[s]
    m = re.match(r"^(.*)_([A-Z][A-Za-z\-]+)$", s)
    return f"{m.group(1).replace('_', ' ')} ({m.group(2)})" if m else s.replace("_", " ")


P = pd.read_csv(R / "E32_primary.csv")
L = P[P.kind == "library"].copy()
cd8 = P[P.signature == "CD8A"].iloc[0]
COMPS = ["ESTIMATE", "TMEscore_z", "NetImmune", "Immune", "CAF_ECM"]
EX = P[P.signature.isin(COMPS)].set_index("signature").loc[COMPS].reset_index()
Q4 = pd.read_csv(R / "E33_expectations.csv").set_index("q").loc["Q4", "result"]
import json
Q4 = json.loads(Q4)
fig = plt.figure(figsize=(180 * MM, 150 * MM))

# ------------------------------------------------------------------ a  pooled interaction by family
axA = S.ax_mm(fig, 14, 9, 96, 58); S.letter_at(fig, 0, 3.0, "a"); title(axA, "Pooled over six randomized comparisons", -0.13)
axA.axhspan(cd8.lo, cd8.hi, color="#EAF2FA", lw=0, zorder=0)
axA.axhline(cd8.est, color=IMMUNE, lw=0.7, zorder=1)
axA.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=1)
GROUPS = [("Inflamed/IFN", "Inflamed/IFN"), ("Other", "Other"), ("Confounder", "Confounder"),
          ("barriernet / comparator", "Composite scores")]
LABEL = {"IMPRES_Auslander": (10, 0), "CYT_Davoli": (0.02, -0.12), "CD39-CD8Tcell_Chow": (0.02, -0.19)}
XMM, YMM, DM = 96 / 4.5, 58 / 1.0, 1.45                     # mm per x unit, per y unit; marker spacing (mm)


def swarm(ys):
    """Beeswarm x offsets (x units) for values ys, placed in order of y, closest free slot to the centre."""
    placed, off = [], np.zeros(len(ys))
    for i in np.argsort(ys):
        for k in range(40):
            dx = (k + 1) // 2 * DM * (1 if k % 2 else -1)
            if all((dx - px) ** 2 + ((ys[i] - py) * YMM) ** 2 >= DM ** 2 for px, py in placed):
                break
        placed.append((dx, ys[i])); off[i] = dx / XMM
    return off

for j, (fam, lab) in enumerate(GROUPS):
    g = EX if fam == "barriernet / comparator" else L[L.family == fam]
    g = g.sort_values("est").reset_index(drop=True)
    n = len(g); xs = j + swarm(g.est.values) if fam != "barriernet / comparator" else j + np.zeros(n)
    for x, (_, r) in zip(xs, g.iterrows()):
        col = FAMCOL[r.family] if fam != "barriernet / comparator" else INK
        beyond = np.isfinite(r.emp_p_declared) and r.emp_p_declared < 0.05
        mk = "D" if fam == "barriernet / comparator" else "o"
        axA.plot(x, r.est, mk, ms=3.6, mfc=col if beyond else "white", mec=col, mew=0.7, zorder=3)
        if r.p < 0.05:
            axA.plot(x, r.est, "o", ms=6.4, mfc="none", mec=INK, mew=0.5, zorder=4)
        if fam == "barriernet / comparator":
            ylab = {"NetImmune": 0.36, "TMEscore_z": 0.31, "ESTIMATE": 0.26}.get(r.signature, r.est)
            axA.annotate(nice(r.signature), (x, r.est), xytext=(x + 0.22, ylab), textcoords="data", fontsize=FSS, color=INK,
                         va="center", ha="left", arrowprops=dict(arrowstyle="-", lw=0.4, color=SLATE, shrinkA=0, shrinkB=2.5)
                         if abs(ylab - r.est) > 0.01 else None)
        elif r.signature in LABEL:
            dx, dy = LABEL[r.signature]
            data = r.signature != "IMPRES_Auslander"
            axA.annotate(nice(r.signature) + (" †" if r.declared_dir < 0 else ""), (x, r.est), xytext=(dx, dy),
                         textcoords="data" if data else "offset points", fontsize=FSS, color=INK,
                         ha="center" if data else "left", va="center",
                         arrowprops=dict(arrowstyle="-", lw=0.4, color=SLATE, shrinkA=1.5, shrinkB=2.5))
    axA.text(j, -0.47, f"{lab} ({n})", fontsize=FS, color=INK, ha="center", va="bottom")
axA.set_xlim(-0.55, 3.95); axA.set_ylim(-0.48, 0.52); axA.set_xticks([])
axA.spines["bottom"].set_visible(False)
axA.set_yticks([-0.4, -0.2, 0, 0.2, 0.4]); axA.set_yticklabels(["−0.4", "−0.2", "0", "0.2", "0.4"])
axA.set_ylabel("Interaction per SD\n(positive favors the checkpoint arm)")
axA.text(3.93, cd8.hi + 0.012, "CD8A alone, 95% CI", fontsize=FSS, color=IMMUNE, ha="right", va="bottom")
for _, r in L.iterrows():
    rec(f"a {r.signature}", f"{r.est:+.3f} ({r.lo:.3f} to {r.hi:.3f}); P = {r.p:.3g}; random-reference P = {r.emp_p_declared:.3g}", "E32_primary.csv")
rec("a CD8A", f"{cd8.est:+.3f} ({cd8.lo:.3f} to {cd8.hi:.3f}); P = {cd8.p:.3g}", "E32_primary.csv")
# key
K = S.ax_mm(fig, 14, 74, 96, 6); S.clean(K); K.set_xlim(0, 96); K.set_ylim(6, 0)
K.plot(1.5, 2, "o", ms=3.6, mfc=SLATE, mec=SLATE); K.text(3.2, 2, "beyond its size-matched random reference", fontsize=FSS, color=SLATE, va="center")
K.plot(1.5, 5.2, "o", ms=3.6, mfc="white", mec=SLATE, mew=0.7); K.text(3.2, 5.2, "not", fontsize=FSS, color=SLATE, va="center")
K.plot(52, 2, "o", ms=6.4, mfc="none", mec=INK, mew=0.5); K.text(54.2, 2, "P < 0.05", fontsize=FSS, color=SLATE, va="center")
K.text(52, 5.2, "† labelled a resistance signature in SignatureSets", fontsize=FSS, color=SLATE, va="center")

# ------------------------------------------------------------------ b  counts
axB = S.ax_mm(fig, 152, 12, 24, 40); S.letter_at(fig, 116, 3.0, "b"); title(axB, "Of 52 published signatures", -1.42)
nfav = int(((L.p < 0.05) & (L.est > 0)).sum()); nag = int(((L.p < 0.05) & (L.est < 0)).sum())
nbey = int((L.emp_p_declared < 0.05).sum())
nabove = int(((L.diff_vs_CD8A > 0) & (L.cd8a_p < 0.05)).sum()); nbelow = int(((L.diff_vs_CD8A < 0) & (L.cd8a_p < 0.05)).sum())
nbelow_res = int(((L.diff_vs_CD8A < 0) & (L.cd8a_p < 0.05) & (L.declared_dir < 0)).sum())
npt = int((L.diff_vs_CD8A > 0).sum())
ROWS = [("favor the checkpoint arm,\nP < 0.05", nfav, FAV), ("favor the control arm,\nP < 0.05", nag, UNFAV),
        ("beyond random genes\nof their size", nbey, INK), ("above CD8A, P < 0.05", nabove, IMMUNE),
        ("positive difference\nfrom CD8A, n.s.", npt, "#9CC3E6"), ("below CD8A, P < 0.05", nbelow, GREY)]
for j, (lab, v, col) in enumerate(ROWS):
    axB.barh(j, v, height=0.62, color=col, lw=0)
    axB.text(v + 0.8, j, str(v), fontsize=FS, color=INK, va="center", fontweight="bold" if lab.startswith("above CD8A, P") else "normal")
    axB.text(-1.5, j, lab, fontsize=FSS, color=INK, va="center", ha="right", linespacing=0.95)
    rec(f"b {lab}", f"{v}/52", "E32_primary.csv")
rec("b below CD8A, declared resistance", f"{nbelow_res}/{nbelow}", "E32_primary.csv")
axB.set_ylim(len(ROWS) - 0.4, -0.6); axB.set_yticks([]); axB.spines["left"].set_visible(False)
axB.set_xlim(0, 30); axB.set_xticks([0, 10, 20, 30]); axB.axvline(0, color=INK, lw=0.5)
axB.set_xlabel("Signatures")

# ------------------------------------------------------------------ c  paired difference from CD8A, ranked
axC = S.ax_mm(fig, 14, 94, 90, 44); S.letter_at(fig, 0, 87.5, "c"); title(axC, "Each signature minus CD8A (paired)", -0.13)
d = L.sort_values("diff_vs_CD8A", ascending=False).reset_index(drop=True)
mdd = Q4["mdd80"]
axC.axhspan(-mdd, mdd, color="#F1F3F5", lw=0, zorder=0)
axC.axhline(0, color=INK, lw=0.5, zorder=1)
for i, r in d.iterrows():
    col = FAMCOL[r.family]; sig = r.cd8a_p < 0.05
    axC.plot([i, i], [r.cd8a_lo, r.cd8a_hi], color=col, lw=0.7, zorder=2)
    axC.plot(i, r.diff_vs_CD8A, "o", ms=2.8, mfc=col if sig else "white", mec=col, mew=0.6, zorder=3)
    rec(f"c {r.signature}", f"{r.diff_vs_CD8A:+.3f} ({r.cd8a_lo:.3f} to {r.cd8a_hi:.3f}); P = {r.cd8a_p:.3g}", "E32_primary.csv")
axC.set_xlim(-1, len(d)); axC.set_xticks([0, 9, 19, 29, 39, 49]); axC.set_xticklabels(["1", "10", "20", "30", "40", "50"])
axC.set_xlabel("Signatures ranked by their difference from CD8A")
axC.set_ylim(-0.85, 0.42); axC.set_yticks([-0.8, -0.4, 0, 0.4]); axC.set_yticklabels(["−0.8", "−0.4", "0", "0.4"])
axC.set_ylabel("Signature − CD8A\n(interaction per SD)")
top = d.iloc[d.cd8a_hi.idxmax()]
axC.annotate(f"largest upper limit +{Q4['max_upper']:.2f}", (d.cd8a_hi.idxmax(), top.cd8a_hi), xytext=(14, 2), textcoords="offset points",
             fontsize=FSS, color=INK, va="center", arrowprops=dict(arrowstyle="-", lw=0.4, color=SLATE, shrinkA=0, shrinkB=1.5))
axC.text(len(d) - 0.5, mdd + 0.015, f"± {mdd:.2f}: difference detectable with 80% power", fontsize=FSS, color=SLATE, ha="right", va="bottom")
axC.text(len(d) - 0.5, -0.80, "filled, P < 0.05", fontsize=FSS, color=SLATE, ha="right", va="bottom")
rec("c detectable difference (80% power)", f"{mdd:.3f} (2.8 x median paired SE {Q4['median_se']:.3f})", "E33_expectations.csv")
rec("c largest upper limit", f"{Q4['max_upper']:.3f} ({Q4['max_upper_sig']}); median upper limit {Q4['median_upper']:.3f}", "E33_expectations.csv")

# ------------------------------------------------------------------ d  each comparison against random sets
PER = pd.read_csv(R / "E32_per_comparison.csv"); PER = PER[PER.variant == "primary"]
Z = np.load(R / "E32_reference_draws.npz")
ORDER = [("IMM150", "IMmotion150"), ("CM025", "CheckMate 025"), ("JR101", "JAVELIN Renal 101"),
         ("ISPY2_pembro", "I-SPY2 pembrolizumab"), ("ISPY2_durva", "I-SPY2 durvalumab–olaparib"), ("NeoTRIP", "NeoTRIP")]
axD = S.ax_mm(fig, 142, 94, 34, 44); S.letter_at(fig, 110, 87.5, "d"); title(axD, "Each comparison", -0.95)
rng = np.random.default_rng(1)
for j, (k, lab) in enumerate(ORDER):
    ys = np.concatenate([Z[key][:, 0] for key in Z.files if key.startswith(k + "|")]); ys = ys[np.isfinite(ys)]
    lo5, hi95 = np.percentile(ys, [5, 95])
    axD.add_patch(plt.Rectangle((lo5, j - 0.36), hi95 - lo5, 0.72, fc=LIGHT, ec="none", zorder=0))
    g = PER[(PER.comparison == k) & PER.signature.isin(L.signature)]
    fam = g.signature.map(dict(zip(L.signature, L.family)))
    axD.scatter(g.y.clip(-1.4, 1.4), j + rng.uniform(-0.24, 0.24, len(g)), s=3.2, c=[FAMCOL[f] for f in fam], edgecolors="none", zorder=2)
    c8 = PER[(PER.comparison == k) & (PER.signature == "CD8A")].y.iloc[0]
    axD.plot([c8, c8], [j - 0.42, j + 0.42], color=IMMUNE, lw=1.1, zorder=3)
    axD.text(-1.55, j, lab + (" (27 scored)" if k == "IMM150" else ""), fontsize=FSS, color=INK, va="center", ha="right")
    rec(f"d {lab}", f"random 5-95% {lo5:.3f} to {hi95:.3f}; CD8A {c8:+.3f}; scored {len(g)}", "E32_per_comparison.csv; E32_reference_draws.npz")
axD.set_ylim(len(ORDER) - 0.45, -0.6); axD.set_yticks([]); axD.spines["left"].set_visible(False)
axD.set_xlim(-1.4, 1.4); axD.set_xticks([-1, 0, 1]); axD.set_xticklabels(["−1", "0", "1"]); S.rule(axD, 0)
axD.set_xlabel("Interaction per SD")
K2 = S.ax_mm(fig, 112, 141, 64, 6); S.clean(K2); K2.set_xlim(0, 64); K2.set_ylim(6, 0)
for x0, (lab, col) in zip([0, 19, 31], [("Inflamed/IFN", IMMUNE), ("Other", OTHER), ("Confounder", STROMA)]):
    K2.plot(x0 + 1, 1.8, "o", ms=2.6, mfc=col, mec=col); K2.text(x0 + 2.6, 1.8, lab, fontsize=FSS, color=SLATE, va="center")
K2.add_patch(plt.Rectangle((0, 4.2), 2.2, 1.6, fc=LIGHT, ec="none")); K2.text(2.6, 5.0, "random sets, 5th–95th pct", fontsize=FSS, color=SLATE, va="center")
K2.plot([45, 45], [3.9, 6.0], color=IMMUNE, lw=1.1); K2.text(46.4, 5.0, "CD8A", fontsize=FSS, color=SLATE, va="center")

S.save(fig, "Fig1")
with open(HERE / "Fig1_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched")
