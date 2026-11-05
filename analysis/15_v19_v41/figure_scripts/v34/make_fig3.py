"""v34 Figure 3 - how far the signature result can be trusted: random genes, single trials, sensitivity, immune infiltration, I-SPY2.

(a) Random gene sets of each signature's size reaching pooled P < 0.05 (1,000 draws each), by signature size.
(b) Each signature in each comparison against single-comparison random sets (5th-95th percentile band); line, CD8A.
(c) Signatures at P < 0.05 under each pre-specified sensitivity analysis.
(d) Pooled interaction before and after adjustment for the ESTIMATE immune score and its interaction with arm.
(e) I-SPY2: MammaPrint adjustment (six-comparison pools) and the processing confound (pembrolizumab comparison, all
    controls against the 40 controls processed like the pembrolizumab arm).
Every number is read from analysis/results/E32_* (protocol analysis/E32_prereg.md). Writes only into figures/cns/v34/.
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
FAMCOL = {"Inflamed/IFN": IMMUNE, "Other": "#8E9AA6", "Confounder": STROMA, "barriernet / comparator": INK}
NUM = {}
rec = lambda k, v, s: NUM.__setitem__(k, (v, s))


def nice(s):
    m = re.match(r"^(.*)_([A-Z][A-Za-z\-]+)$", s)
    return f"{m.group(1).replace('_', ' ')} ({m.group(2)})" if m else s.replace("_", " ")


def title(ax, s, x=0.0):
    ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=6.5, fontweight="bold", color=INK, va="bottom", ha="left")


P = pd.read_csv(R / "E32_primary.csv"); L = P[P.kind == "library"].copy()
fig = plt.figure(figsize=(180 * MM, 150 * MM))

# ------------------------------------------------------------------ a  random genes reach significance
COV = pd.read_csv(R / "E32_coverage.csv")
SIX = ["IMM150", "ISPY2_pembro", "ISPY2_durva", "CM025", "JR101", "NeoTRIP"]
size = COV[COV.scored & COV.comparison.isin(SIX)].groupby("signature").present.mean()
L["size"] = L.signature.map(size)
axA = S.ax_mm(fig, 12, 10, 46, 46); S.letter_at(fig, 0, 3.5, "a"); title(axA, "Random genes reach significance", -0.2)
axA.axhline(5, color=GREY, lw=0.5, ls=(0, (2, 2)))
axA.scatter(L["size"], 100 * L.ref_frac_nominal_p05, s=9, c=[FAMCOL[f] for f in L.family], edgecolors="white", linewidths=0.3, zorder=3)
axA.set_xscale("log"); axA.set_xticks([3, 10, 30, 100, 300]); axA.set_xticklabels(["3", "10", "30", "100", "300"])
axA.set_xlabel("Genes in the signature"); axA.set_ylabel("Random sets of that size\nwith pooled P < 0.05 (%)")
axA.set_ylim(0, 23); axA.text(330, 5.6, "5% expected", fontsize=4.8, color=SLATE, ha="right")
top = L.sort_values("ref_frac_nominal_p05").iloc[-1]
axA.annotate(f"{nice(top.signature)}\n{100 * top.ref_frac_nominal_p05:.0f}%", (top["size"], 100 * top.ref_frac_nominal_p05),
             xytext=(-5, -1), textcoords="offset points", fontsize=4.8, color=INK, ha="right", va="center")
rho = pd.Series(L["size"]).corr(L.ref_frac_nominal_p05, method="spearman")
axA.text(0.03, 0.97, f"Spearman ρ = {rho:.2f}", transform=axA.transAxes, fontsize=4.8, color=SLATE, va="top")
for j, (lab, col) in enumerate([("Inflamed/IFN", FAMCOL["Inflamed/IFN"]), ("Other", FAMCOL["Other"]), ("Confounder", STROMA)]):
    axA.scatter([], [], s=9, c=col, label=lab)
axA.legend(loc="upper left", bbox_to_anchor=(0.0, 0.88), frameon=False, fontsize=4.8, handletextpad=0.2, borderaxespad=0.2)
rec("a random nominal P<0.05", f"median {100 * L.ref_frac_nominal_p05.median():.1f}%, range {100 * L.ref_frac_nominal_p05.min():.1f}-{100 * L.ref_frac_nominal_p05.max():.1f}%; rho {rho:.2f}", "E32_primary.csv; E32_coverage.csv")

# ------------------------------------------------------------------ b  each trial against its random sets
PER = pd.read_csv(R / "E32_per_comparison.csv"); PER = PER[PER.variant == "primary"]
Z = np.load(R / "E32_reference_draws.npz")
ORDER = [("IMM150", "IMmotion150"), ("CM025", "CheckMate 025"), ("JR101", "JAVELIN Renal 101"),
         ("ISPY2_pembro", "I-SPY2 pembrolizumab"), ("ISPY2_durva", "I-SPY2 durvalumab–olaparib"), ("NeoTRIP", "NeoTRIP")]
axB = S.ax_mm(fig, 92, 10, 34, 46); S.letter_at(fig, 64, 3.5, "b"); title(axB, "Each comparison", -0.82)
rng = np.random.default_rng(1)
for j, (k, lab) in enumerate(ORDER):
    ys = np.concatenate([Z[key][:, 0] for key in Z.files if key.startswith(k + "|")]); ys = ys[np.isfinite(ys)]
    lo5, hi95 = np.percentile(ys, [5, 95])
    axB.add_patch(plt.Rectangle((lo5, j - 0.36), hi95 - lo5, 0.72, fc=LIGHT, ec="none", zorder=0))
    g = PER[(PER.comparison == k) & PER.signature.isin(L.signature)]
    fam = g.signature.map(dict(zip(L.signature, L.family)))
    axB.scatter(g.y.clip(-1.6, 1.6), j + rng.uniform(-0.22, 0.22, len(g)), s=3.0, c=[FAMCOL[f] for f in fam], edgecolors="none", zorder=2)
    c8 = PER[(PER.comparison == k) & (PER.signature == "CD8A")].y.iloc[0]
    axB.plot([c8, c8], [j - 0.42, j + 0.42], color=IMMUNE, lw=1.0, zorder=3)
    above = int((g.y > hi95).sum())
    axB.text(-1.75, j, lab, fontsize=4.9, color=INK, va="center", ha="right")
    axB.text(1.95, j, f"{above}/{len(g)}", fontsize=4.8, color=INK, va="center", ha="center")
    rec(f"b {lab}", f"random 5-95% {lo5:.3f} to {hi95:.3f}; above 95th {above}/{len(g)}; CD8A {c8:+.3f}", "E32_per_comparison.csv; E32_reference_draws.npz")
axB.set_ylim(len(ORDER) - 0.4, -0.6); axB.set_yticks([]); axB.spines["left"].set_visible(False)
axB.set_xlim(-1.6, 2.3); axB.set_xticks([-1, 0, 1]); axB.set_xticklabels(["−1", "0", "1"]); S.rule(axB, 0)
axB.set_xlabel("Interaction per SD in that comparison")
axB.text(1.95, -0.95, "above\n95th", fontsize=4.4, color=SLATE, ha="center", va="bottom")
axB.text(0.02, -0.08, "band, random sets (5th–95th pct); blue line, CD8A", transform=axB.transAxes, fontsize=4.6, color=SLATE,
         va="top") if False else None

# ------------------------------------------------------------------ c  sensitivity
A = pd.read_csv(R / "E32_pooled_all_sets.csv"); NV = pd.read_csv(R / "E32_native.csv")
SETS = [("six (primary)", "Primary (six comparisons)"), ("six + PRINCE + DUTRENEO", "+ PRINCE and DUTRENEO"),
        ("six, renal adjusted for angiogenesis", "Renal adjusted for angiogenesis"),
        ("six, I-SPY2 hormone receptor only", "I-SPY2 without MammaPrint"), ("native", "Published scoring engines"),
        ("five without IMmotion150", "Without IMmotion150"), ("six, CheckMate 025 overall survival", "CheckMate 025 overall survival"),
        ("six, adjusted for ESTIMATE immune (Q4)", "Adjusted for ESTIMATE immune")]
axC = S.ax_mm(fig, 152, 10, 24, 46); S.letter_at(fig, 130, 3.5, "c"); title(axC, "Sensitivity", -0.92)
for j, (key, lab) in enumerate(SETS):
    g = NV if key == "native" else A[(A.set == key) & (A.kind == "library")]
    fav = int(((g.p < 0.05) & (g.est > 0)).sum()); ag = int(((g.p < 0.05) & (g.est < 0)).sum())
    axC.barh(j, fav, height=0.62, color=FAV, lw=0); axC.barh(j, -ag, height=0.62, color=UNFAV, lw=0)
    axC.text(fav + 0.6, j, f"{fav}" + (f" / {ag}" if ag else ""), fontsize=4.8, color=INK, va="center")
    axC.text(-3.2, j, lab, fontsize=4.9, color=INK, va="center", ha="right")
    rec(f"c {lab}", f"favor {fav}, against {ag}", "E32_pooled_all_sets.csv / E32_native.csv")
axC.set_ylim(len(SETS) - 0.4, -0.6); axC.set_yticks([]); axC.spines["left"].set_visible(False)
axC.set_xlim(-3, 27); axC.set_xticks([0, 10, 20]); axC.axvline(0, color=INK, lw=0.5)
axC.set_xlabel("Signatures at P < 0.05 (of 52)")
axC.text(27, -0.85, "favoring / against", fontsize=4.6, color=SLATE, ha="right", va="center")

# ------------------------------------------------------------------ d  adjustment for general immune infiltration
axD = S.ax_mm(fig, 12, 78, 46, 52); S.letter_at(fig, 0, 71.5, "d"); title(axD, "Adjusted for ESTIMATE immune", -0.2)
lim = (-0.45, 0.75)
axD.plot(lim, lim, color=GREY, lw=0.5, ls=(0, (2, 2))); S.rule(axD, 0); axD.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
for _, r in L.iterrows():
    col = FAMCOL[r.family]; sig = r.p < 0.05
    axD.plot(r.est, r.est_estimate_adj, "o", ms=2.8 if sig else 2.0, mfc=col if r.p_estimate_adj < 0.05 else ("white" if sig else col),
             mec=col, mew=0.6, alpha=1 if sig else 0.45, zorder=3 if sig else 2)
s18 = L[L.p < 0.05]; ratio = (s18.est_estimate_adj / s18.est)
axD.set_xlim(*lim); axD.set_ylim(*lim)
axD.set_xticks([-0.4, 0, 0.4]); axD.set_yticks([-0.4, 0, 0.4]); axD.set_xticklabels(["−0.4", "0", "0.4"]); axD.set_yticklabels(["−0.4", "0", "0.4"])
axD.set_xlabel("Pooled interaction, unadjusted"); axD.set_ylabel("Adjusted for ESTIMATE immune × arm")
axD.text(0.03, 0.97, f"{len(s18)} signatures at P < 0.05 (large markers):\nestimate × {ratio.median():.2f} (median) after adjustment;\n"
         f"{int((s18.p_estimate_adj < 0.05).sum())} still at P < 0.05 (filled)", transform=axD.transAxes, fontsize=4.8, color=INK, va="top")
axD.text(0.62, 0.70, "y = x", fontsize=4.6, color=SLATE, rotation=45, ha="center")
rec("d adjusted/unadjusted (18 significant)", f"median {ratio.median():.2f}, range {ratio.min():.2f}-{ratio.max():.2f}; kept P<0.05 {int((s18.p_estimate_adj < 0.05).sum())}/{len(s18)}", "E32_primary.csv")

# ------------------------------------------------------------------ e  I-SPY2 checks
Q = pd.read_csv(R / "E32_ispy2.csv")
axE = S.ax_mm(fig, 88, 78, 88, 52); S.letter_at(fig, 64, 71.5, "e"); title(axE, "I-SPY2: MammaPrint and processing", -0.27)
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
        rec(f"e {sc} {lab}", f"{r.est:+.3f} ({r.lo:.3f} to {r.hi:.3f}); P = {r.p:.3g}", "E32_ispy2.csv")
    if sc != SC[-1]:
        axE.axhline(y + 4.0, color=LIGHT, lw=0.4)
    y += 5
for yy, sc in YT:
    axE.text(-1.95, yy, {"Immune": "Immune axis", "CAF_ECM": "CAF_ECM"}.get(sc, sc), fontsize=5.2, color=INK, va="center", ha="right",
             fontweight="bold")
axE.set_ylim(y - 1.4, -0.8); axE.set_yticks([]); axE.spines["left"].set_visible(False)
axE.set_xlim(-1.8, 2.9); axE.set_xticks([-1, 0, 1, 2]); axE.set_xticklabels(["−1", "0", "1", "2"]); S.rule(axE, 0)
axE.set_xlabel("Interaction per SD (positive favors checkpoint arm)")
K = S.ax_mm(fig, 88, 136, 88, 10); S.clean(K); K.set_xlim(0, 88); K.set_ylim(10, 0)
for j, (lab, col, mk) in enumerate([("six comparisons, I-SPY2 hormone receptor only", GREY, "o"),
                                     ("six comparisons, + MammaPrint (primary)", INK, "o"),
                                     ("I-SPY2 pembrolizumab, all controls", GREY, "s"),
                                     ("I-SPY2 pembrolizumab, Linear-adjusted controls only", INK, "s")]):
    xx, yy = (0, 1.5 + 3.6 * j) if j < 2 else (44, 1.5 + 3.6 * (j - 2))
    K.plot(xx + 1, yy, mk, ms=2.6, mfc=col, mec=col); K.text(xx + 2.6, yy, lab, fontsize=4.6, color=SLATE, va="center")
K.text(0, 8.9, "open, P ≥ 0.05", fontsize=4.6, color=SLATE, va="center")

S.save(fig, "Fig3")
with open(HERE / "Fig3_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched")
