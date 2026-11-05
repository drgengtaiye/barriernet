"""v45 Supplementary Fig. S35: barriernet 0.3 against each published signature (a; analysis/results/E42_paired.csv) and
ranking by absolute benefit against CD8A (b, c; E43b_contrasts.csv, null-centered, with the pre-specified
cross-validated values of E43_contrasts.csv). Writes only into figures/cns/v45/.
"""
import sys, re, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import cns_style as S
S.apply_large()
from cns_style import INK, SLATE, GREY, LIGHT, PAPER, IMMUNE, STROMA, MM

HERE = Path(__file__).resolve().parent
S.OUT = HERE
GUARD = {p: hashlib.md5(open(p, "rb").read()).hexdigest() for p in glob.glob(str(HERE.parent / "Fig*.p*"))}
R = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
OTHER = "#8E9AA6"
FAMCOL = {"Inflamed/IFN": IMMUNE, "Other": OTHER, "Confounder": STROMA}
FS, FSS, FN = 7.0, 6.5, 6.0
NUM = {}
rec = lambda k, v, s: NUM.__setitem__(k, (v, s))


def nice(s):
    m = re.match(r"^(.*)_([A-Z][A-Za-z\-]+)$", s)
    return f"{m.group(1).replace('_', ' ')} ({m.group(2)})" if m else s.replace("_", " ")


def title(ax, s, x=0.0):
    ax.text(x, 1.0, s, transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=INK, va="bottom", ha="left")


fig = plt.figure(figsize=(165 * MM, 170 * MM))
from matplotlib.transforms import blended_transform_factory as _bt

# ------------------------------------------------------------------ a  barriernet 0.3 and CD8A against each signature
P = pd.read_csv(R / "E42_paired.csv")
P = P[(P.set == "six (primary)") & (P.orientation == "conservative")].sort_values("sig_est_E32", ascending=False).reset_index(drop=True)
assert len(P) == 52
axA = S.ax_mm(fig, 46, 12, 53, 132); title(axA, "Each published signature, in its better direction", -0.78)
XL = (-0.32, 0.62)
axA.axvspan(-0.10, 0.10, color=PAPER, zorder=0, lw=0)
for i, r in P.iterrows():
    for dy, d, lo, hi, p, col in [(-0.17, r["diff"], r.lo, r.hi, r.p, INK), (0.17, r.cd8a_diff, r.cd8a_lo, r.cd8a_hi, r.cd8a_p, GREY)]:
        axA.plot([max(lo, XL[0]), min(hi, XL[1])], [i + dy, i + dy], color=col, lw=0.55, solid_capstyle="butt")
        axA.plot(d, i + dy, "o", ms=2.3, mfc=col if p < 0.05 else "white", mec=col, mew=0.55, zorder=3)
    axA.text(XL[0] - 0.012, i, nice(r.signature), fontsize=FN, color=INK, va="center", ha="right")
    axA.plot(3.5 / 165, i, "s", ms=2.2, color=FAMCOL[r.family], clip_on=False, transform=_bt(fig.transFigure, axA.transData))
    rec(f"S35a {r.signature}", f"barriernet 0.3 {r['diff']:+.3f} ({r.lo:+.3f} to {r.hi:+.3f}), P = {r.p:.3g}; CD8A "
        f"{r.cd8a_diff:+.3f} ({r.cd8a_lo:+.3f} to {r.cd8a_hi:+.3f}), P = {r.cd8a_p:.3g}", "E42_paired.csv")
axA.set_ylim(len(P) - 0.4, -0.6); axA.set_yticks([]); axA.spines["left"].set_visible(False)
axA.set_xlim(*XL); axA.set_xticks([-0.2, 0, 0.2, 0.4, 0.6]); axA.set_xticklabels(["−0.2", "0", "0.2", "0.4", "0.6"])
S.rule(axA, 0)
axA.set_xlabel("Paired difference in interaction per SD\n(score minus signature)")
na, nb = int(((P["diff"] > 0) & (P.p < 0.05)).sum()), int(((P["diff"] < 0) & (P.p < 0.05)).sum())
nc = int(((P.cd8a_diff > 0) & (P.cd8a_p < 0.05)).sum())
rec("S35a counts", f"barriernet 0.3 above {na}, below {nb}; CD8A above {nc}", "E42_paired.csv")
assert (na, nb, nc) == (6, 0, 4)
S.letter_at(fig, 2, 6, "a")

# ------------------------------------------------------------------ b, c  absolute benefit, difference from CD8A
C1 = pd.read_csv(R / "E43b_contrasts.csv"); C0 = pd.read_csv(R / "E43_contrasts.csv")
MOD = [("M_2AX", "Immune axis × arm,\n+ CAF_ECM (primary)"), ("M_BN03", "barriernet 0.3 × arm"),
       ("M_CD8A_CAF", "CD8A × arm,\n+ CAF_ECM"), ("M_2AXI", "Immune axis × arm,\nCAF_ECM × arm"),
       ("M_RISK", "Immune axis, CAF_ECM,\nno interaction")]
for j, (metric, ttl, xl, xt, y0) in enumerate([("C-for-benefit", "Ranking by absolute benefit: C-for-benefit", (-0.15, 0.15),
                                                 [-0.1, 0, 0.1], 12),
                                                ("gradient", "Benefit gradient (upper − lower half)", (-0.30, 0.30),
                                                 [-0.2, 0, 0.2], 92)]):
    ax = S.ax_mm(fig, 131, y0, 32, 44); title(ax, ttl, -1.02)
    for i, (m, lab) in enumerate(MOD):
        a = C1[(C1.set == "six") & (C1.contrast == f"{m} - M_CD8A") & (C1.metric == metric)].iloc[0]
        b = C0[(C0.set == "six (primary)") & (C0.contrast == f"{m} - M_CD8A") & (C0.metric == metric)].iloc[0]
        for dy, r, col, mk in [(-0.16, a, INK, "o"), (0.16, b, GREY, "s")]:
            ax.plot([max(r.lo, xl[0]), min(r.hi, xl[1])], [i + dy, i + dy], color=col, lw=0.6)
            ax.plot(r["diff"], i + dy, mk, ms=2.6, mfc=col if r.p < 0.05 else "white", mec=col, mew=0.6, zorder=3)
        ax.text(xl[0] - 0.03 * (xl[1] - xl[0]), i, lab, fontsize=FN, color=INK, va="center", ha="right", linespacing=1.0)
        rec(f"S35{'bc'[j]} {m} {metric}", f"null-centered {a['diff']:+.3f} ({a.lo:+.3f} to {a.hi:+.3f}), P = {a.p:.3g}; "
            f"cross-validated {b['diff']:+.3f} ({b.lo:+.3f} to {b.hi:+.3f})", "E43b_contrasts.csv; E43_contrasts.csv")
    ax.set_ylim(len(MOD) - 0.5, -0.5); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(*xl); ax.set_xticks(xt); ax.set_xticklabels([f"{v:g}".replace("-", "−") for v in xt]); S.rule(ax, 0)
    ax.set_xlabel("Difference from CD8A × arm")
    S.letter_at(fig, 101, y0 - 6, "bc"[j])

# ------------------------------------------------------------------ keys
K = S.ax_mm(fig, 0, 157, 100, 12); S.clean(K); K.set_xlim(0, 100); K.set_ylim(12, 0)
for x0, col, lab in [(2, INK, "barriernet 0.3 minus signature"), (50, GREY, "CD8A minus signature")]:
    K.plot([x0, x0 + 3], [2, 2], color=col, lw=0.6); K.plot(x0 + 1.5, 2, "o", ms=2.3, mfc=col, mec=col)
    K.text(x0 + 4.5, 2, lab, fontsize=FSS, color=SLATE, va="center")
for x0, (fam, col) in zip([2, 30, 46], FAMCOL.items()):
    K.plot(x0 + 1.5, 6.2, "s", ms=2.2, color=col); K.text(x0 + 3.5, 6.2, {"Inflamed/IFN": "inflamed / IFN", "Other": "other",
                                                                            "Confounder": "confounder"}[fam],
                                                     fontsize=FSS, color=SLATE, va="center")
K.add_patch(plt.Rectangle((66, 5.2), 3, 2, color=PAPER, lw=0)); K.text(70.5, 6.2, "±0.10 per SD", fontsize=FSS, color=SLATE, va="center")
K.text(2, 10.2, "filled, P < 0.05; open, P ≥ 0.05", fontsize=FSS, color=SLATE, va="center")
K2 = S.ax_mm(fig, 103, 150, 62, 14); S.clean(K2); K2.set_xlim(0, 62); K2.set_ylim(14, 0)
for yy, col, mk, lab in [(2, INK, "o", "full-data fit, centered on 200 permutations"), (6.2, GREY, "s", "cross-validated (pre-specified)")]:
    K2.plot(1.5, yy, mk, ms=2.6, mfc=col, mec=col); K2.text(3.5, yy, lab, fontsize=FSS, color=SLATE, va="center")
K2.text(0.5, 10.4, "six comparisons; fixed-effect pooling", fontsize=FSS, color=SLATE, va="center")

S.save(fig, "FigS35")
with open(HERE / "FigS35_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched (v45)")
