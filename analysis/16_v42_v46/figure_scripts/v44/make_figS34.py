"""v44 Supplementary Fig. S34: the I-SPY2 MammaPrint and processing checks, formerly Fig. 5f (v43/make_fig5.py, panel f,
unchanged except for its placement). Numbers from analysis/results/E32_ispy2.csv. Writes only into figures/cns/v44/.
"""
import sys, re, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import cns_style as S
S.apply_large()
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, FAV, UNFAV, MM

HERE = Path(__file__).resolve().parent
S.OUT = HERE
GUARD = {p: hashlib.md5(open(p, "rb").read()).hexdigest() for p in glob.glob(str(HERE.parent / "Fig*.p*"))}
R = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
OTHER = "#8E9AA6"
FAMCOL = {"Inflamed/IFN": IMMUNE, "Other": OTHER, "Confounder": STROMA}
FS, FSS = 7.0, 6.5
NUM = {}
rec = lambda k, v, s: NUM.__setitem__(k, (v, s))


def nice(s):
    lab = {"TMEscore_z": "TMEscore", "Immune": "Immune axis"}
    if s in lab:
        return lab[s]
    m = re.match(r"^(.*)_([A-Z][A-Za-z\-]+)$", s)
    return f"{m.group(1).replace('_', ' ')} ({m.group(2)})" if m else s.replace("_", " ")


def title(ax, s, x=0.0):
    ax.text(x, 1.04, s, transform=ax.transAxes, fontsize=7.8, fontweight="bold", color=INK, va="bottom", ha="left")


fig = plt.figure(figsize=(120 * MM, 82 * MM))

# ------------------------------------------------------------------ e  I-SPY2 checks
Q = pd.read_csv(R / "E32_ispy2.csv")
axE = S.ax_mm(fig, 26, 9, 90, 46); title(axE, "I-SPY2: MammaPrint and processing", -0.25)
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
        rec(f"S34 {sc} {lab}", f"{r.est:+.3f} ({r.lo:.3f} to {r.hi:.3f}); P = {r.p:.3g}", "E32_ispy2.csv")
    if sc != SC[-1]:
        axE.axhline(y + 4.0, color=LIGHT, lw=0.4)
    y += 5
for yy, sc in YT:
    axE.text(-1.95, yy, {"Immune": "Immune axis", "CAF_ECM": "CAF_ECM"}.get(sc, sc), fontsize=FS, color=INK, va="center", ha="right",
             fontweight="bold")
axE.set_ylim(y - 1.4, -0.8); axE.set_yticks([]); axE.spines["left"].set_visible(False)
axE.set_xlim(-1.8, 2.9); axE.set_xticks([-1, 0, 1, 2]); axE.set_xticklabels(["−1", "0", "1", "2"]); S.rule(axE, 0)
axE.set_xlabel("Interaction per SD (positive favors checkpoint arm)")
K = S.ax_mm(fig, 2, 64, 116, 16); S.clean(K); K.set_xlim(0, 100); K.set_ylim(16, 0)
for j, (lab, col, mk) in enumerate([("six comparisons, I-SPY2 hormone receptor only", GREY, "o"),
                                     ("six comparisons, + MammaPrint (primary)", INK, "o"),
                                     ("I-SPY2 pembrolizumab, all controls", GREY, "s"),
                                     ("I-SPY2 pembrolizumab, Linear-adjusted controls only", INK, "s")]):
    xx, yy = 0, 1.5 + 3.1 * j
    K.plot(xx + 1, yy, mk, ms=2.6, mfc=col, mec=col); K.text(xx + 2.6, yy, lab, fontsize=FSS, color=SLATE, va="center")
K.text(0, 1.5 + 3.1 * 4, "open, P ≥ 0.05", fontsize=FSS, color=SLATE, va="center")

S.save(fig, "FigS34")
with open(HERE / "FigS34_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched (v44)")
