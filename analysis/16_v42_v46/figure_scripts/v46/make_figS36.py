"""v46 Supplementary Fig. S36: gastric and gastroesophageal cancer (seventeenth later protocol, analysis/E44_prereg.md).
(a) the two randomized gastric trials of the published-trials analysis (E44_A_gastric_trials.csv); (b) CD8A and
barriernet 0.3 per single-arm cohort and pooled (E44_B_per_cohort.csv, E44_B_pooled.csv, E44_B_adjusted.csv); (c) the
52 signatures in the three cohorts against their six-comparison interaction (E44_B_pooled.csv, E32_primary.csv); (d)
the two gastric cohorts before and after adjustment for MSI-H/EBV-positive status (E44_B_adjusted.csv). Writes only
into figures/cns/v46/.
"""
import sys, glob, hashlib
from pathlib import Path
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
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
neg = lambda s: s.replace("-", "−")


def title(ax, s, x=0.0):
    ax.text(x, 1.0, s, transform=ax.transAxes, fontsize=7.5, fontweight="bold", color=INK, va="bottom", ha="left")


def ticks(ax, xt, axis="x"):
    getattr(ax, f"set_{axis}ticks")(xt)
    getattr(ax, f"set_{axis}ticklabels")([neg(f"{v:g}") for v in xt])


def ci_line(ax, y, est, lo, hi, xl, col, mk="o", ms=2.6, filled=True):
    ax.plot([max(lo, xl[0]), min(hi, xl[1])], [y, y], color=col, lw=0.6, solid_capstyle="butt")
    for edge, lim, d in [(lo, xl[0], -1), (hi, xl[1], 1)]:            # clipped interval: open arrowhead at the axis edge
        if (d < 0 and edge < lim) or (d > 0 and edge > lim):
            ax.plot(lim, y, marker="<" if d < 0 else ">", ms=2.2, color=col, clip_on=False)
    ax.plot(est, y, mk, ms=ms, mfc=col if filled else "white", mec=col, mew=0.6, zorder=3)


fig = plt.figure(figsize=(165 * MM, 147 * MM))

# ------------------------------------------------------------------ a  randomized gastric trials (published subgroups)
A = pd.read_csv(R / "E44_A_gastric_trials.csv").set_index("row")
ROWS_A = [("CheckMate 649", "CheckMate 649 (n = 402)\nnivolumab + ipilimumab vs chemo"),
          ("KEYNOTE-061", "KEYNOTE-061 (n = 274)\npembrolizumab vs paclitaxel"),
          ("gastric two, fixed effect", "Both gastric trials"),
          ("other 16 trials, DerSimonian-Laird", "Other 16 trials")]
axA = S.ax_mm(fig, 40, 12, 32, 40); title(axA, "Published randomized gastric trials", -1.19)
XA = (-0.6, 0.6)
for i, (k, lab) in enumerate(ROWS_A):
    r = A.loc[k]
    ci_line(axA, i, r.per_sd, r.per_sd_lo, r.per_sd_hi, XA, INK if i < 3 else GREY, mk="o" if i < 2 else "D",
            filled=r.p < 0.05)
    axA.text(XA[0] - 0.04 * (XA[1] - XA[0]), i, lab, fontsize=FN, color=INK, va="center", ha="right", linespacing=1.0)
    axA.text(XA[1] + 0.05 * (XA[1] - XA[0]), i, f"{r.mdd80:.2f}", fontsize=FN, color=SLATE, va="center", ha="left")
    rec(f"S36a {k}", f"{r.per_sd:+.3f} ({r.per_sd_lo:+.3f} to {r.per_sd_hi:+.3f}), P = {r.p:.3g}; detectable {r.mdd80:.3f}",
        "E44_A_gastric_trials.csv")
axA.text(XA[1] + 0.05 * (XA[1] - XA[0]), -0.62, "detectable", fontsize=FN, color=SLATE, va="bottom", ha="left")
axA.axvline(0.2264, color=IMMUNE, lw=0.6, ls=(0, (1, 1.5)), zorder=0)
axA.text(0.19, 2.5, "CD8A, six\ncomparisons", fontsize=FN, color=IMMUNE, va="center", ha="right", linespacing=1.0)
axA.set_ylim(3.5, -0.5); axA.set_yticks([]); axA.spines["left"].set_visible(False)
axA.set_xlim(*XA); ticks(axA, [-0.4, 0, 0.4]); S.rule(axA, 0)
axA.set_xlabel("Interaction per SD (approximate)")
assert abs(A.loc["gastric two, fixed effect"].per_sd - 0.102) < 0.001
S.letter_at(fig, 2, 6, "a")

# ------------------------------------------------------------------ b  CD8A and barriernet 0.3 in the single-arm cohorts
PC = pd.read_csv(R / "E44_B_per_cohort.csv"); PO = pd.read_csv(R / "E44_B_pooled.csv"); AD = pd.read_csv(R / "E44_B_adjusted.csv")
THREE, TWO = "three gastroesophageal (primary)", "two gastric"
ROWS_B = [("PRJEB25780", "PRJEB25780, gastric\npembrolizumab (44; 12)"),
          ("PRJEB60680", "NCT04249739, gastric\nchemo, then + pembrolizumab (43; 29)"),
          ("PERFECT", "PERFECT, esophageal\natezolizumab + CRT (32; 12)"),
          ("fe3", "Three cohorts, fixed effect"), ("cons3", "Three cohorts, main-text rule"),
          ("fe2", "Two gastric"), ("adj2", "Two gastric, MSI-H/EBV-adjusted")]
axB = S.ax_mm(fig, 128, 12, 35, 58); title(axB, "Single-arm checkpoint cohorts", -1.05)
XB = (-1.0, 3.0)
for i, (k, lab) in enumerate(ROWS_B):
    for dy, nm, col in [(-0.17, "barriernet 0.3", INK), (0.17, "CD8A", GREY)]:
        if k in ("PRJEB25780", "PRJEB60680", "PERFECT"):
            r = PC[(PC.cohort == k) & (PC.name == nm)].iloc[0]
            e, lo, hi, p = r.logOR, r.logOR - 1.96 * r.se, r.logOR + 1.96 * r.se, r.p
        elif k == "adj2":
            r = AD[(AD.set == TWO) & (AD.name == nm)].iloc[0]; e, lo, hi, p = r.est_msiebv, r.lo_msiebv, r.hi_msiebv, r.p_msiebv
        else:
            r = PO[(PO.set == (THREE if k.endswith("3") else TWO)) & (PO.name == nm)].iloc[0]
            e, lo, hi, p = (r.cons_est, r.cons_lo, r.cons_hi, r.cons_p) if k == "cons3" else (r.est, r.lo, r.hi, r.p)
        ci_line(axB, i + dy, e, lo, hi, XB, col, mk="o" if i < 3 else "D", ms=2.3, filled=p < 0.05)
        rec(f"S36b {k} {nm}", f"{e:+.3f} ({lo:+.3f} to {hi:+.3f}), P = {p:.3g}", "E44_B_per_cohort / pooled / adjusted")
    axB.text(XB[0] - 0.04 * (XB[1] - XB[0]), i, lab, fontsize=FN, color=INK, va="center", ha="right", linespacing=1.0)
axB.axhline(2.5, color=LIGHT, lw=0.5)
axB.set_ylim(len(ROWS_B) - 0.5, -0.5); axB.set_yticks([]); axB.spines["left"].set_visible(False)
axB.set_xlim(*XB); ticks(axB, [-1, 0, 1, 2, 3]); S.rule(axB, 0)
axB.set_xlabel("Log-odds of response per SD")
S.letter_at(fig, 90, 6, "b")

# ------------------------------------------------------------------ c  agreement with the randomized interaction
E = pd.read_csv(R / "E32_primary.csv"); E = E[E.set == "six (primary)"].set_index("signature")
LIB = E[E.kind == "library"]
P3 = PO[(PO.set == THREE)].set_index("name")
q = LIB[["est", "family"]].join(P3[["est", "p"]], rsuffix="_g")
assert len(q) == 52 and q.est_g.notna().all()
rho = spearmanr(q.est, q.est_g); assert abs(rho.statistic - 0.760) < 0.001
B3 = pd.read_csv(R / "E37_test_pooled.csv"); bn6 = B3[(B3.score == "barriernet 0.3") & (B3.set == "six (primary)")].est.iloc[0]
axC = S.ax_mm(fig, 14, 86, 52, 44); title(axC, "Agreement with the randomized comparisons", -0.14)
for fam, col in FAMCOL.items():
    d = q[q.family == fam]
    axC.plot(d.est, d.est_g, "o", ms=2.6, mfc=col, mec="white", mew=0.3, zorder=3)
for nm, x, y, mk, dx, dy in [("CD8A", E.loc["CD8A", "est"], P3.loc["CD8A", "est"], "D", 0.012, -0.10),
                             ("barriernet 0.3", bn6, P3.loc["barriernet 0.3", "est"], "s", 0.014, -0.02),
                             ("TMEscore", E.loc["TMEscore_z", "est"], P3.loc["TMEscore_z", "est"], "^", 0.012, 0.0)]:
    axC.plot(x, y, mk, ms=3.0, mfc="white", mec=INK, mew=0.7, zorder=4)
    axC.text(x + dx, y + dy, nm, fontsize=FN, color=INK, va="center", ha="left")
    rec(f"S36c {nm}", f"randomized {x:+.3f}; gastroesophageal {y:+.3f}", "E32_primary.csv / E37_test_pooled.csv; E44_B_pooled.csv")
axC.text(0.97, 0.04, f"Spearman ρ = {rho.statistic:.2f}\n(52 signatures)", transform=axC.transAxes, fontsize=FSS,
         color=SLATE, va="bottom", ha="right", linespacing=1.1)
rec("S36c rho", f"{rho.statistic:.4f}, P = {rho.pvalue:.2g}", "E44_B_pooled.csv; E32_primary.csv")
axC.set_xlim(-0.35, 0.45); axC.set_ylim(-0.9, 1.1); ticks(axC, [-0.2, 0, 0.2, 0.4]); ticks(axC, [-0.5, 0, 0.5, 1], "y")
S.rule(axC, 0); S.hrule(axC, 0)
axC.set_xlabel("Interaction per SD, six randomized comparisons"); axC.set_ylabel("Log-odds of response per SD\n(three cohorts)")
S.letter_at(fig, 2, 80, "c")

# ------------------------------------------------------------------ d  MSI-H / EBV adjustment, two gastric cohorts
A2 = AD[AD.set == TWO].set_index("name")
d52 = A2.loc[LIB.index].join(LIB.family)
n_un, n_adj = int((d52.p < 0.05).sum()), int((d52.p_msiebv < 0.05).sum()); assert (n_un, n_adj) == (17, 1)
axD = S.ax_mm(fig, 108, 86, 52, 44); title(axD, "MSI-H/EBV adjustment, two gastric cohorts", -0.14)
LIMD = (-0.9, 1.3)
axD.plot(LIMD, LIMD, color=GREY, lw=0.5, zorder=0)
for fam, col in FAMCOL.items():
    d = d52[d52.family == fam]
    axD.plot(d.est, d.est_msiebv, "o", ms=2.6, mfc=col, mec="white", mew=0.3, zorder=3)
for nm, lab, mk, tx, ty in [("CD8A", "CD8A", "D", 1.0, 0.36), ("barriernet 0.3", "barriernet 0.3", "s", 0.86, 0.12),
                            ("TMEscore_z", "TMEscore", "^", 0.76, -0.08)]:
    x, y = A2.loc[nm, "est"], A2.loc[nm, "est_msiebv"]
    axD.plot(x, y, mk, ms=3.0, mfc="white", mec=INK, mew=0.7, zorder=4)
    axD.annotate(lab, (x, y), (tx, ty), fontsize=FN, color=INK, va="center", ha="left",
                 arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4, shrinkA=0, shrinkB=2.5))
    rec(f"S36d {lab}", f"unadjusted {x:+.3f} (P = {A2.loc[nm, 'p']:.3g}); adjusted {y:+.3f} (P = {A2.loc[nm, 'p_msiebv']:.3g})",
        "E44_B_adjusted.csv")
axD.text(0.97, 0.04, f"signatures at P < 0.05:\n{n_un} of 52 unadjusted, {n_adj} adjusted", transform=axD.transAxes,
         fontsize=FSS, color=SLATE, va="bottom", ha="right", linespacing=1.1)
rec("S36d counts", f"{n_un} unadjusted, {n_adj} adjusted ({', '.join(d52[d52.p_msiebv < 0.05].index)})", "E44_B_adjusted.csv")
axD.set_xlim(*LIMD); axD.set_ylim(*LIMD); ticks(axD, [-0.5, 0, 0.5, 1]); ticks(axD, [-0.5, 0, 0.5, 1], "y")
S.rule(axD, 0); S.hrule(axD, 0)
axD.set_xlabel("Log-odds per SD, unadjusted"); axD.set_ylabel("Adjusted for MSI-H/EBV")
S.letter_at(fig, 84, 80, "d")

# ------------------------------------------------------------------ keys
K = S.ax_mm(fig, 84, 1, 81, 5); S.clean(K); K.patch.set_alpha(0); K.set_xlim(0, 81); K.set_ylim(5, 0)
for x0, col, lab in [(44, INK, "barriernet 0.3"), (64, GREY, "CD8A")]:
    K.plot([x0, x0 + 3], [2.5, 2.5], color=col, lw=0.6); K.plot(x0 + 1.5, 2.5, "o", ms=2.3, mfc=col, mec=col)
    K.text(x0 + 4.5, 2.5, lab, fontsize=FSS, color=SLATE, va="center")
K2 = S.ax_mm(fig, 0, 140, 165, 6); S.clean(K2); K2.patch.set_alpha(0); K2.set_xlim(0, 165); K2.set_ylim(6, 0)
x0 = 14
for fam, col in FAMCOL.items():
    K2.plot(x0, 3, "o", ms=2.6, mfc=col, mec="white", mew=0.3)
    lab = {"Inflamed/IFN": "inflamed / IFN signature", "Other": "other signature", "Confounder": "confounder"}[fam]
    K2.text(x0 + 2, 3, lab, fontsize=FSS, color=SLATE, va="center"); x0 += 4 + 1.6 * len(lab)
for x0, mk, lab in [(108, "D", "CD8A"), (124, "s", "barriernet 0.3"), (148, "^", "TMEscore")]:
    K2.plot(x0, 3, mk, ms=3.0, mfc="white", mec=INK, mew=0.7); K2.text(x0 + 2.2, 3, lab, fontsize=FSS, color=SLATE, va="center")
K3 = S.ax_mm(fig, 2, 63, 80, 8); S.clean(K3); K3.patch.set_alpha(0); K3.set_xlim(0, 80); K3.set_ylim(8, 0)
K3.text(0, 2, "a, b: filled, P < 0.05; open, P ≥ 0.05; diamonds, pooled; 95% CI", fontsize=FSS, color=SLATE, va="center")
K3.text(0, 6, "b: (patients; responders); CRT, chemoradiotherapy", fontsize=FSS, color=SLATE, va="center")

S.save(fig, "FigS36")
with open(HERE / "FigS36_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers; figures/cns untouched (v46)")
