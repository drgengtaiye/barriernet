"""v47 Supplementary Fig. S36: eighteenth later protocol (analysis/E45_prereg.md). (a) Each gene's pooled single-arm log-odds against its
six-comparison interaction (E45_single_arm_gene_pooled.csv, E33_pooled.csv). (b) Each signature minus CD8A and minus CXCL10 (E45_A_paired.csv).
(c) Empirical P of the 52 signatures against the size-matched and the coexpression-matched reference (E45_B_coexpression.csv). (d) Within-arm
slopes of CD8A and of the Angio score in the renal comparisons (E45_C_within_arm.csv, E45_C_pooled.csv). Writes only into figures/cns/v47/."""
import sys, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import cns_style as S
S.apply_large()
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, MM
HERE = Path(__file__).resolve().parent; S.OUT = HERE
GUARD = {p: hashlib.md5(open(p, "rb").read()).hexdigest() for p in glob.glob(str(HERE.parent / "Fig*.p*"))}
R = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
OTHER = "#8E9AA6"; CHK = IMMUNE; CTL = "#B4532A"
FAMCOL = {"Inflamed/IFN": IMMUNE, "Other": OTHER, "Confounder": STROMA}
FS, FSS, FN = 7.0, 6.5, 6.0
NUM = {}; rec = lambda k, v, s: NUM.__setitem__(k, (v, s)); neg = lambda s: s.replace("-", "−")


def title(ax, s, x=0.0):
    ax.text(x, 1.0, s, transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=INK, va="bottom", ha="left")


def ticks(ax, xt, axis="x"):
    getattr(ax, f"set_{axis}ticks")(xt); getattr(ax, f"set_{axis}ticklabels")([neg(f"{v:g}") for v in xt])


fig = plt.figure(figsize=(180 * MM, 128 * MM))
# ---- a
SA = pd.read_csv(R / "E45_single_arm_gene_pooled.csv"); G = pd.read_csv(R / "E33_pooled.csv"); G = G[(G.set == "six (primary)") & (G.kind == "gene") & (G.all_six.astype(str) == "True")].set_index("name")
q = SA[SA.eligible].set_index("gene").join(G[["est"]], rsuffix="_rct"); assert len(q) == 1996
axA = S.ax_mm(fig, 16, 11, 62, 42); title(axA, "Genes: single-arm cohorts against randomized comparisons", -0.2)
axA.plot(q.est, q.est_rct, "o", ms=1.5, mfc="#C9D1D8", mec="none", zorder=1)
for g, c, dx, dy in [("CXCL10", CTL, 0.012, 0.035), ("CD8A", IMMUNE, -0.012, 0.045), ("STAT1", OTHER, 0.012, 0.045), ("CXCL13", OTHER, -0.01, 0.06), ("CXCL9", OTHER, -0.012, -0.05)]:
    if g in q.index:
        axA.plot(q.loc[g, "est"], q.loc[g, "est_rct"], "o", ms=3.0, mfc=c, mec="white", mew=0.3, zorder=3)
        axA.text(q.loc[g, "est"] + dx, q.loc[g, "est_rct"] + dy, g, fontsize=FN, color=INK, ha="left" if dx > 0 else "right", va="center")
        rec(f"a {g}", f"single-arm {q.loc[g, 'est']:+.3f}; randomized {q.loc[g, 'est_rct']:+.3f}", "E45_single_arm_gene_pooled.csv; E33_pooled.csv")
axA.set_xlim(-0.35, 0.55); axA.set_ylim(-0.35, 0.4); ticks(axA, [-0.2, 0, 0.2, 0.4]); ticks(axA, [-0.2, 0, 0.2], "y"); S.rule(axA, 0); S.hrule(axA, 0)
axA.set_xlabel("Pooled log-odds of response per SD, 20 single-arm cohorts"); axA.set_ylabel("Interaction per SD, six comparisons")
S.letter_at(fig, 2, 5, "a")
# ---- b
P = pd.read_csv(R / "E45_A_paired.csv"); P = P[(P.set == "six (primary)") & (P.family != "composite")]
c8 = P[P.reference == "CD8A"].set_index("score"); cx = P[P.reference == "CXCL10"].set_index("score"); j = c8[["diff", "p", "family"]].join(cx[["diff", "p"]], rsuffix="_x"); assert len(j) == 52
axB = S.ax_mm(fig, 112, 11, 62, 42); title(axB, "Each signature against CD8A and against CXCL10", -0.2)
for fam, col in FAMCOL.items():
    d = j[j.family == fam]
    axB.plot(d["diff"], d.diff_x, "o", ms=2.6, mfc="white", mec=col, mew=0.6, zorder=2)
    e = d[d.p_x < 0.05]; axB.plot(e["diff"], e.diff_x, "o", ms=2.6, mfc=col, mec=col, mew=0.6, zorder=3)
n_up, n_dn = int(((j.diff_x > 0) & (j.p_x < 0.05)).sum()), int(((j.diff_x < 0) & (j.p_x < 0.05)).sum()); assert (n_up, n_dn) == (8, 2)
assert int(((j["diff"] > 0) & (j.p < 0.05)).sum()) == 0
axB.text(0.97, 0.05, f"filled, P < 0.05 against CXCL10:\n{n_up} above, {n_dn} below", transform=axB.transAxes, fontsize=FSS, color=SLATE, va="bottom", ha="right", linespacing=1.1)
axB.set_xlim(-0.55, 0.15); axB.set_ylim(-0.45, 0.3); ticks(axB, [-0.4, -0.2, 0]); ticks(axB, [-0.4, -0.2, 0, 0.2], "y"); S.rule(axB, 0); S.hrule(axB, 0)
axB.set_xlabel("Signature minus CD8A, per SD"); axB.set_ylabel("Signature minus CXCL10, per SD"); rec("b counts", f"{n_up} above, {n_dn} below CXCL10; 0 above CD8A", "E45_A_paired.csv")
S.letter_at(fig, 96, 5, "b")
# ---- c
X = pd.read_csv(R / "E45_B_coexpression.csv")
axC = S.ax_mm(fig, 16, 72, 62, 42); title(axC, "Random reference: size only and size plus coexpression", -0.2)
for fam, col in FAMCOL.items():
    d = X[X.family == fam]; axC.plot(-np.log10(d.emp_p_size), -np.log10(d.emp_p_coex), "o", ms=2.6, mfc=col, mec="white", mew=0.3, zorder=3)
L5 = -np.log10(0.05); axC.axvline(L5, color=GREY, lw=0.5, ls=(0, (2, 2))); axC.axhline(L5, color=GREY, lw=0.5, ls=(0, (2, 2))); axC.plot([0, 3.1], [0, 3.1], color=LIGHT, lw=0.5, zorder=0)
nA, nB = int((X.emp_p_size < .05).sum()), int((X.emp_p_coex < .05).sum()); assert (nA, nB) == (26, 16)
axC.text(0.03, 0.95, f"beyond the reference (P < 0.05):\n{nA} of 52, size only; {nB}, with coexpression", transform=axC.transAxes, fontsize=FSS, color=SLATE, va="top", ha="left", linespacing=1.1)
axC.set_xlim(0, 3.1); axC.set_ylim(0, 3.1); axC.set_xticks([0, 1, 2, 3]); axC.set_yticks([0, 1, 2, 3]); axC.set_xlabel("−log10 empirical P, size-matched"); axC.set_ylabel("−log10 empirical P,\ncoexpression-matched")
rec("c counts", f"{nA}; {nB}", "E45_B_coexpression.csv"); S.letter_at(fig, 2, 66, "c")
# ---- d
W_ = pd.read_csv(R / "E45_C_within_arm.csv"); PL = pd.read_csv(R / "E45_C_pooled.csv"); PL = PL[PL.set == "renal three"].set_index("score")
ROWS = [("IMM150", "IMmotion150"), ("CM025", "CheckMate 025"), ("JR101", "JAVELIN Renal 101"), ("pool", "Three pooled")]
for j_, (sc, lab, x0) in enumerate([("CD8A", "CD8A", 124), ("Angio", "Angio score", 154)]):
    ax = S.ax_mm(fig, x0, 72, 24, 42); ax.text(0.5, 1.0, lab, transform=ax.transAxes, fontsize=FS, color=INK, ha="center", va="bottom", fontweight="bold")
    XL = (-1.3, 1.3)
    for i, (k, nm) in enumerate(ROWS):
        for dy, arm, col in [(-0.16, "icb", CHK), (0.16, "ctl", CTL)]:
            if k == "pool":
                r = PL.loc[sc]; e, lo, hi, p = r[f"{arm}_est"], r[f"{arm}_lo"], r[f"{arm}_hi"], r[f"{arm}_p"]; mk = "D"
            else:
                r = W_[(W_.score == sc) & (W_.comparison == k)].iloc[0]; e, s_ = r[f"{arm}_y"], r[f"{arm}_se"]; lo, hi = e - 1.96 * s_, e + 1.96 * s_; p = 2 * (1 - __import__("scipy.stats").stats.norm.cdf(abs(e / s_))); mk = "o"
            ax.plot([max(lo, XL[0]), min(hi, XL[1])], [i + dy, i + dy], color=col, lw=0.6); ax.plot(e, i + dy, mk, ms=2.4, mfc=col if p < 0.05 else "white", mec=col, mew=0.6, zorder=3)
            rec(f"d {sc} {k} {arm}", f"{e:+.3f} ({lo:+.3f} to {hi:+.3f}), P = {p:.3g}", "E45_C_within_arm.csv; E45_C_pooled.csv")
        if j_ == 0:
            ax.text(XL[0] - 0.08 * (XL[1] - XL[0]), i, nm, fontsize=FN, color=INK, ha="right", va="center")
    ax.axhline(2.5, color=LIGHT, lw=0.5); ax.set_ylim(3.5, -0.5); ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xlim(*XL); ticks(ax, [-1, 0, 1]); S.rule(ax, 0)
    ax.set_xlabel("Slope per SD")
S.letter_at(fig, 96, 66, "d")
K = S.ax_mm(fig, 0, 122, 180, 5); S.clean(K); K.patch.set_alpha(0); K.set_xlim(0, 180); K.set_ylim(5, 0); x0 = 16
for fam, col in FAMCOL.items():
    K.plot(x0, 2.5, "o", ms=2.6, mfc=col, mec="white", mew=0.3); lab = {"Inflamed/IFN": "inflamed / IFN signature", "Other": "other signature", "Confounder": "confounder"}[fam]
    K.text(x0 + 2, 2.5, lab, fontsize=FSS, color=SLATE, va="center"); x0 += 4 + 1.6 * len(lab)
for x0, col, lab in [(112, CHK, "checkpoint arm"), (142, CTL, "control arm")]:
    K.plot([x0, x0 + 3], [2.5, 2.5], color=col, lw=0.6); K.plot(x0 + 1.5, 2.5, "o", ms=2.3, mfc=col, mec=col); K.text(x0 + 4.5, 2.5, lab, fontsize=FSS, color=SLATE, va="center")
S.save(fig, "FigS36")
with open(HERE / "FigS36_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers")
