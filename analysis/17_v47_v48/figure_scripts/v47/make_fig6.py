"""v47 Figure 6: gastric and gastroesophageal cancer. (a) prognosis in the ten gastric survival cohorts without checkpoint inhibitors
(E30_dedup_runs/primary/E20_pooled.csv, as Fig. 7c of v46); (b) the two randomized gastric trials of the published-trials analysis
(E44_A_gastric_trials.csv); (c) CD8A and barriernet 0.3 per single-arm cohort and pooled (E44_B_*.csv); (d) response by MSI-H or
EBV-positive status in the two gastric cohorts (E44_labels.csv, counts only); (e) before and after adjustment for that status
(E44_B_adjusted.csv); (f) the 52 signatures in the three cohorts against their six-comparison interaction (E44_B_pooled.csv,
E32_primary.csv). Panels b, c, e and f were Supplementary Fig. S36a, b, d and c of v46. Writes only into figures/cns/v47/."""
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
    for edge, lim, d in [(lo, xl[0], -1), (hi, xl[1], 1)]:
        if (d < 0 and edge < lim) or (d > 0 and edge > lim):
            ax.plot(lim, y, marker="<" if d < 0 else ">", ms=2.2, color=col, clip_on=False)
    ax.plot(est, y, mk, ms=ms, mfc=col if filled else "white", mec=col, mew=0.6, zorder=3)


fig = plt.figure(figsize=(180 * MM, 176 * MM))

# ------------------------------------------------------------------ a  prognosis, ten gastric survival cohorts
E20 = pd.read_csv(R / "E30_dedup_runs/primary/E20_pooled.csv")
E20 = E20[(E20.set == "all 10 gastric") & (E20.model == "unadjusted") & (E20.duplicate_rule == "pre-registered exclusions (primary)")].set_index("score"); assert E20.index.is_unique and abs(E20.loc["CAF_ECM", "HR"] - 1.307) < 0.001
ROWS = [("NetImmune", "NetImmune", INK), ("ESTIMATE", "ESTIMATE", GREY), ("Immune", "Immune axis", IMMUNE), ("CD8A", "CD8A", IMMUNE), ("CAF_ECM", "CAF_ECM (fibroblast)", STROMA)]
axA = S.ax_mm(fig, 34, 12, 46, 30); title(axA, "Prognosis without checkpoint inhibitors", -0.72)
XP = (0.6, 1.6)
axA.set_xscale("log")
for i, (k, lab, col) in enumerate(ROWS):
    r = E20.loc[k]; ci_line(axA, i, r.HR, r.lo, r.hi, XP, col, mk="D", filled=r.p < 0.05)
    axA.text(XP[0] * 0.96, i, lab, fontsize=FN, color=INK, va="center", ha="right")
    axA.text(XP[1] * 1.03, i, f"{r.HR:.2f}", fontsize=FN, color=SLATE, va="center", ha="left")
    rec(f"a {k}", f"HR {r.HR:.3f} ({r.lo:.3f}-{r.hi:.3f}), P = {r.p:.3g}", "E30_dedup_runs/primary/E20_pooled.csv")
axA.axvline(1, color=GREY, lw=0.5, zorder=0); axA.set_ylim(len(ROWS) - 0.5, -0.5); axA.set_yticks([]); axA.spines["left"].set_visible(False)
axA.set_xlim(*XP); axA.set_xticks([0.6, 0.8, 1, 1.25, 1.6]); axA.set_xticklabels(["0.6", "0.8", "1", "1.25", "1.6"]); axA.minorticks_off()
axA.set_xlabel("HR per SD, ten gastric cohorts (2,144 patients)")
S.letter_at(fig, 2, 6, "a")

# ------------------------------------------------------------------ b  randomized gastric trials (published subgroups)
A = pd.read_csv(R / "E44_A_gastric_trials.csv").set_index("row")
ROWS_A = [("CheckMate 649", "CheckMate 649 (n = 402)\nnivolumab + ipilimumab vs chemo"), ("KEYNOTE-061", "KEYNOTE-061 (n = 274)\npembrolizumab vs paclitaxel"),
          ("gastric two, fixed effect", "Both gastric trials"), ("other 16 trials, DerSimonian-Laird", "Other 16 trials")]
axB = S.ax_mm(fig, 130, 12, 34, 30); title(axB, "Published randomized gastric trials", -1.13)
XA = (-0.6, 0.6)
for i, (k, lab) in enumerate(ROWS_A):
    r = A.loc[k]
    ci_line(axB, i, r.per_sd, r.per_sd_lo, r.per_sd_hi, XA, INK if i < 3 else GREY, mk="o" if i < 2 else "D", filled=r.p < 0.05)
    axB.text(XA[0] - 0.04 * (XA[1] - XA[0]), i, lab, fontsize=FN, color=INK, va="center", ha="right", linespacing=1.0)
    axB.text(XA[1] + 0.05 * (XA[1] - XA[0]), i, f"{r.mdd80:.2f}", fontsize=FN, color=SLATE, va="center", ha="left")
    rec(f"b {k}", f"{r.per_sd:+.3f} ({r.per_sd_lo:+.3f} to {r.per_sd_hi:+.3f}), P = {r.p:.3g}; detectable {r.mdd80:.3f}", "E44_A_gastric_trials.csv")
axB.text(XA[1] + 0.05 * (XA[1] - XA[0]), -0.62, "detectable", fontsize=FN, color=SLATE, va="bottom", ha="left")
axB.axvline(0.2264, color=IMMUNE, lw=0.6, ls=(0, (1, 1.5)), zorder=0)
axB.text(0.19, 2.5, "CD8A, six\ncomparisons", fontsize=FN, color=IMMUNE, va="center", ha="right", linespacing=1.0)
axB.set_ylim(3.5, -0.5); axB.set_yticks([]); axB.spines["left"].set_visible(False)
axB.set_xlim(*XA); ticks(axB, [-0.4, 0, 0.4]); S.rule(axB, 0)
axB.set_xlabel("Interaction per SD (approximate)")
assert abs(A.loc["gastric two, fixed effect"].per_sd - 0.102) < 0.001
S.letter_at(fig, 84, 6, "b")

# ------------------------------------------------------------------ c  CD8A and barriernet 0.3 in the single-arm cohorts
PC = pd.read_csv(R / "E44_B_per_cohort.csv"); PO = pd.read_csv(R / "E44_B_pooled.csv"); AD = pd.read_csv(R / "E44_B_adjusted.csv")
THREE, TWO = "three gastroesophageal (primary)", "two gastric"
ROWS_B = [("PRJEB25780", "PRJEB25780, gastric\npembrolizumab (44; 12)"), ("PRJEB60680", "NCT04249739, gastric\nchemo, then + pembrolizumab (43; 29)"),
          ("PERFECT", "PERFECT, esophageal\natezolizumab + CRT (32; 12)"), ("fe3", "Three cohorts, fixed effect"), ("cons3", "Three cohorts, main-text rule"),
          ("fe2", "Two gastric"), ("adj2", "Two gastric, MSI-H/EBV-adjusted")]
axC = S.ax_mm(fig, 48, 60, 58, 46); title(axC, "Single-arm checkpoint cohorts", -0.75)
XB = (-1.0, 3.0)
for i, (k, lab) in enumerate(ROWS_B):
    for dy, nm, col in [(-0.17, "barriernet 0.3", INK), (0.17, "CD8A", GREY)]:
        if k in ("PRJEB25780", "PRJEB60680", "PERFECT"):
            r = PC[(PC.cohort == k) & (PC.name == nm)].iloc[0]; e, lo, hi, p = r.logOR, r.logOR - 1.96 * r.se, r.logOR + 1.96 * r.se, r.p
        elif k == "adj2":
            r = AD[(AD.set == TWO) & (AD.name == nm)].iloc[0]; e, lo, hi, p = r.est_msiebv, r.lo_msiebv, r.hi_msiebv, r.p_msiebv
        else:
            r = PO[(PO.set == (THREE if k.endswith("3") else TWO)) & (PO.name == nm)].iloc[0]
            e, lo, hi, p = (r.cons_est, r.cons_lo, r.cons_hi, r.cons_p) if k == "cons3" else (r.est, r.lo, r.hi, r.p)
        ci_line(axC, i + dy, e, lo, hi, XB, col, mk="o" if i < 3 else "D", ms=2.3, filled=p < 0.05)
        rec(f"c {k} {nm}", f"{e:+.3f} ({lo:+.3f} to {hi:+.3f}), P = {p:.3g}", "E44_B_per_cohort / pooled / adjusted")
    axC.text(XB[0] - 0.03 * (XB[1] - XB[0]), i, lab, fontsize=FN, color=INK, va="center", ha="right", linespacing=1.0)
axC.axhline(2.5, color=LIGHT, lw=0.5); axC.set_ylim(len(ROWS_B) - 0.5, -0.5); axC.set_yticks([]); axC.spines["left"].set_visible(False)
axC.set_xlim(*XB); ticks(axC, [-1, 0, 1, 2, 3]); S.rule(axC, 0); axC.set_xlabel("Log-odds of response per SD")
for x0, col, lab in [(0.70, INK, "barriernet 0.3"), (0.70, GREY, "CD8A")]:
    y0 = 0.76 if lab != "CD8A" else 0.70
    axC.plot([x0, x0 + 0.05], [y0, y0], transform=axC.transAxes, color=col, lw=0.6); axC.plot(x0 + 0.025, y0, "o", transform=axC.transAxes, ms=2.3, mfc=col, mec=col)
    axC.text(x0 + 0.075, y0, lab, transform=axC.transAxes, fontsize=FN, color=SLATE, va="center")
S.letter_at(fig, 2, 54, "c")

# ------------------------------------------------------------------ d  response by MSI-H or EBV-positive status
LB = pd.read_csv(R / "E44_labels.csv"); G2 = LB[LB.cohort.isin(["PRJEB25780", "PRJEB60680"])]
c1 = G2[G2.msi_ebv == 1]; c0 = G2[G2.msi_ebv == 0]; n1, r1, n0, r0 = len(c1), int(c1.response.sum()), len(c0), int(c0.response.sum())
assert (r1, n1, r0, n0) == (13, 14, 28, 73)
axD = S.ax_mm(fig, 130, 60, 34, 46); title(axD, "Response by tumor status, two gastric cohorts", -0.42)
axD.bar([0, 1], [100 * r1 / n1, 100 * r0 / n0], width=0.55, color=[IMMUNE, OTHER])
for x, a_, b_ in [(0, r1, n1), (1, r0, n0)]:
    axD.text(x, 100 * a_ / b_ + 2, f"{a_} of {b_}", fontsize=FN, color=INK, ha="center", va="bottom")
axD.set_xticks([0, 1]); axD.set_xticklabels(["MSI-H or\nEBV-positive", "Other\ntumors"]); axD.set_ylim(0, 110); axD.set_yticks([0, 25, 50, 75, 100]); axD.set_ylabel("Responders (%)")
axD.set_xlim(-0.6, 1.6); rec("d counts", f"{r1} of {n1}; {r0} of {n0}", "E44_labels.csv")
S.letter_at(fig, 110, 54, "d")

# ------------------------------------------------------------------ e  MSI-H / EBV adjustment, two gastric cohorts
E = pd.read_csv(R / "E32_primary.csv"); E = E[E.set == "six (primary)"].set_index("signature"); LIB = E[E.kind == "library"]
A2 = AD[AD.set == TWO].set_index("name"); d52 = A2.loc[LIB.index].join(LIB.family)
n_un, n_adj = int((d52.p < 0.05).sum()), int((d52.p_msiebv < 0.05).sum()); assert (n_un, n_adj) == (17, 1)
axE = S.ax_mm(fig, 16, 124, 60, 40); title(axE, "MSI-H/EBV adjustment, two gastric cohorts", -0.14)
LIMD = (-0.9, 1.3)
axE.plot(LIMD, LIMD, color=GREY, lw=0.5, zorder=0)
for fam, col in FAMCOL.items():
    d = d52[d52.family == fam]; axE.plot(d.est, d.est_msiebv, "o", ms=2.6, mfc=col, mec="white", mew=0.3, zorder=3)
for nm, lab, mk, tx, ty in [("CD8A", "CD8A", "D", 1.0, 0.36), ("barriernet 0.3", "barriernet 0.3", "s", 0.86, 0.12), ("TMEscore_z", "TMEscore", "^", 0.76, -0.08)]:
    x, y = A2.loc[nm, "est"], A2.loc[nm, "est_msiebv"]; axE.plot(x, y, mk, ms=3.0, mfc="white", mec=INK, mew=0.7, zorder=4)
    axE.annotate(lab, (x, y), (tx, ty), fontsize=FN, color=INK, va="center", ha="left", arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4, shrinkA=0, shrinkB=2.5))
    rec(f"e {lab}", f"unadjusted {x:+.3f} (P = {A2.loc[nm, 'p']:.3g}); adjusted {y:+.3f} (P = {A2.loc[nm, 'p_msiebv']:.3g})", "E44_B_adjusted.csv")
axE.text(0.97, 0.04, f"signatures at P < 0.05:\n{n_un} of 52 unadjusted, {n_adj} adjusted", transform=axE.transAxes, fontsize=FSS, color=SLATE, va="bottom", ha="right", linespacing=1.1)
axE.set_xlim(*LIMD); axE.set_ylim(*LIMD); ticks(axE, [-0.5, 0, 0.5, 1]); ticks(axE, [-0.5, 0, 0.5, 1], "y"); S.rule(axE, 0); S.hrule(axE, 0)
axE.set_xlabel("Log-odds per SD, unadjusted"); axE.set_ylabel("Adjusted for MSI-H/EBV")
S.letter_at(fig, 2, 118, "e")

# ------------------------------------------------------------------ f  agreement with the randomized interaction
P3 = PO[(PO.set == THREE)].set_index("name"); q = LIB[["est", "family"]].join(P3[["est", "p"]], rsuffix="_g"); assert len(q) == 52 and q.est_g.notna().all()
rho = spearmanr(q.est, q.est_g); assert abs(rho.statistic - 0.760) < 0.001
B3 = pd.read_csv(R / "E37_test_pooled.csv"); bn6 = B3[(B3.score == "barriernet 0.3") & (B3.set == "six (primary)")].est.iloc[0]
axF = S.ax_mm(fig, 108, 124, 60, 40); title(axF, "Agreement with the randomized comparisons", -0.14)
for fam, col in FAMCOL.items():
    d = q[q.family == fam]; axF.plot(d.est, d.est_g, "o", ms=2.6, mfc=col, mec="white", mew=0.3, zorder=3)
for nm, x, y, mk, dx, dy in [("CD8A", E.loc["CD8A", "est"], P3.loc["CD8A", "est"], "D", 0.012, -0.10), ("barriernet 0.3", bn6, P3.loc["barriernet 0.3", "est"], "s", 0.014, -0.02),
                             ("TMEscore", E.loc["TMEscore_z", "est"], P3.loc["TMEscore_z", "est"], "^", 0.012, 0.0)]:
    axF.plot(x, y, mk, ms=3.0, mfc="white", mec=INK, mew=0.7, zorder=4); axF.text(x + dx, y + dy, nm, fontsize=FN, color=INK, va="center", ha="left")
    rec(f"f {nm}", f"randomized {x:+.3f}; gastroesophageal {y:+.3f}", "E32_primary.csv / E37_test_pooled.csv; E44_B_pooled.csv")
axF.text(0.97, 0.04, f"Spearman ρ = {rho.statistic:.2f}\n(52 signatures)", transform=axF.transAxes, fontsize=FSS, color=SLATE, va="bottom", ha="right", linespacing=1.1)
rec("f rho", f"{rho.statistic:.4f}, P = {rho.pvalue:.2g}", "E44_B_pooled.csv; E32_primary.csv")
axF.set_xlim(-0.35, 0.45); axF.set_ylim(-0.9, 1.1); ticks(axF, [-0.2, 0, 0.2, 0.4]); ticks(axF, [-0.5, 0, 0.5, 1], "y"); S.rule(axF, 0); S.hrule(axF, 0)
axF.set_xlabel("Interaction per SD, six randomized comparisons"); axF.set_ylabel("Log-odds of response per SD\n(three cohorts)")
S.letter_at(fig, 94, 118, "f")

# ------------------------------------------------------------------ keys
K2 = S.ax_mm(fig, 0, 170, 180, 5); S.clean(K2); K2.patch.set_alpha(0); K2.set_xlim(0, 180); K2.set_ylim(5, 0)
x0 = 16
for fam, col in FAMCOL.items():
    K2.plot(x0, 2.5, "o", ms=2.6, mfc=col, mec="white", mew=0.3)
    lab = {"Inflamed/IFN": "inflamed / IFN signature", "Other": "other signature", "Confounder": "confounder"}[fam]
    K2.text(x0 + 2, 2.5, lab, fontsize=FSS, color=SLATE, va="center"); x0 += 4 + 1.6 * len(lab)
for x0, mk, lab in [(116, "D", "CD8A"), (134, "s", "barriernet 0.3"), (160, "^", "TMEscore")]:
    K2.plot(x0, 2.5, mk, ms=3.0, mfc="white", mec=INK, mew=0.7); K2.text(x0 + 2.2, 2.5, lab, fontsize=FSS, color=SLATE, va="center")
K3 = S.ax_mm(fig, 2, 48.5, 120, 5); S.clean(K3); K3.patch.set_alpha(0); K3.set_xlim(0, 120); K3.set_ylim(5, 0)
K3.text(0, 2.5, "a–c: filled, P < 0.05; open, P ≥ 0.05; diamonds, pooled; 95% CI. c: (patients; responders); CRT, chemoradiotherapy", fontsize=FN, color=SLATE, va="center")

S.save(fig, "Fig6")
with open(HERE / "Fig6_numbers.md", "w") as fh:
    fh.write("| item | value | source |\n|---|---|---|\n" + "".join(f"| {k} | {v} | {s} |\n" for k, (v, s) in NUM.items()))
for p, h in GUARD.items():
    assert hashlib.md5(open(p, "rb").read()).hexdigest() == h, f"modified: {p}"
print("ok;", len(NUM), "numbers")
