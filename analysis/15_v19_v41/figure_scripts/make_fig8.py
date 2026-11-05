"""v19 Fig. 8: synthesis - each component against each question, on one favorable-to-the-right scale.

Every estimate is read from a stored result file (no new statistics):
  gastric survival   E20_pooled.csv      'all 10 gastric', unadjusted, pre-registered exclusions
  TCGA survival      E11_pooled.csv      '(i) SKCM KIRC LUAD LUSC BLCA', OS, adjusted
  response, 21       E12_icb_pooled.csv  'all'; the wider of the DL and HK intervals
  gastric response   E22_meta_inputs.csv per cohort, Wald 95% CI
  interaction        E21_meta.csv        'primary: 6 comparisons' (conservative interval)
Survival effects are shown as -log HR per SD so that, in every column, right means favorable:
longer survival, response, or more benefit from the checkpoint arm.
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
import cns_style as S

R = S.RES
ROWS = [("CAF_ECM", "CAF_ECM"), ("TAM_SUP", "TAM_SUP"), ("Barrier", "Barrier"),
        ("APM", "APM"), ("IFN", "IFN"), ("TLS_B", "TLS_B"), ("T_CELL", "T_CELL"),
        ("Immune", "Immune"), ("NetImmune", "NetImmune")]
GROUP_BREAKS = {3, 8}                 # thin rule before APM and before NetImmune


def gastric_os():
    p = pd.read_csv(R / "E20_pooled.csv")
    p = p[(p.set == "all 10 gastric") & (p.model == "unadjusted")
          & (p.duplicate_rule == "pre-registered exclusions (primary)")].set_index("score")
    return {s: (-np.log(p.at[s, "HR"]), -np.log(p.at[s, "hi"]), -np.log(p.at[s, "lo"]), p.at[s, "p"])
            for s, _ in ROWS}


def tcga_os():
    e = pd.read_csv(R / "E11_pooled.csv")
    e = e[(e.group == "(i) SKCM KIRC LUAD LUSC BLCA") & (e.endpoint == "OS")
          & (e.model == "adjusted")].set_index("score")
    return {s: (-np.log(e.at[s, "HR"]), -np.log(e.at[s, "hi"]), -np.log(e.at[s, "lo"]), e.at[s, "p"])
            for s, _ in ROWS}


def icb21():
    i = pd.read_csv(R / "E12_icb_pooled.csv")
    i = i[i.set == "all"].set_index("score")
    out = {}
    for s, _ in ROWS:
        r = i.loc[s]
        wide_hk = (r.hi_HK - r.lo_HK) >= (r.hi_DL - r.lo_DL)
        lo, hi = (r.lo_HK, r.hi_HK) if wide_hk else (r.lo_DL, r.hi_DL)
        out[s] = (r.logOR, lo, hi, max(r.p_DL, r.p_HK), int(r.n_pos), int(r.k))
    return out


def gastric_icb():
    g = pd.read_csv(R / "E22_meta_inputs.csv")
    out = {}
    for c in ("PRJEB25780", "PRJEB60680"):
        d = g[g.cohort == c].set_index("score")
        out[c] = {s: (d.at[s, "y"], d.at[s, "y"] - 1.96 * d.at[s, "se"], d.at[s, "y"] + 1.96 * d.at[s, "se"],
                      d.at[s, "p_cohort"]) for s, _ in ROWS}
    return out


def trials():
    m = pd.read_csv(R / "E21_meta.csv")
    m = m[m.set == "primary: 6 comparisons"].set_index("score")
    return {s: (m.at[s, "est"], m.at[s, "lo"], m.at[s, "hi"], m.at[s, "p"], int(m.at[s, "k_benefit"]),
                int(m.at[s, "k"])) for s, _ in ROWS}


COLS = [
    ("Survival,\ngastric cancer", "10 cohorts; 2,152 patients", "−log HR per SD", gastric_os(), (-0.45, 0.45)),
    ("Survival, five checkpoint\ncancer types (TCGA)", "untreated; adjusted", "−log HR per SD", tcga_os(),
     (-0.45, 0.45)),
    ("Response, 21 checkpoint\ncohorts", "1,340 patients", "log-odds per SD", icb21(), (-0.6, 0.9)),
    ("Response, two gastric\ncheckpoint cohorts", "44 and 43 patients", "log-odds per SD", gastric_icb(),
     (-1.9, 2.8)),
    ("Benefit, six randomized\ncomparisons", "score × arm; 1,828 patients", "interaction per SD",
     trials(), (-0.6, 0.9)),
]


def colour(lo, hi):
    if lo > 0:
        return S.FAV, True
    if hi < 0:
        return S.UNFAV, True
    return S.GREY, False


FIG_W, FIG_H = 180, 169
fig = plt.figure(figsize=(FIG_W * S.MM, FIG_H * S.MM))

# ---------------------------------------------------------------- a. forest matrix
LAB_W, X0, COL_W, GAP = 21, 23, 29.2, 2.2
TOP, ROW_H = 14.0, 6.4
H = ROW_H * len(ROWS)
S.letter_at(fig, 1, 1.5, "a")
fig.text(7 / FIG_W, 1 - 2.1 / FIG_H, "Effect per SD for each outcome", fontsize=6.5,
         fontweight="bold", va="top", ha="left", color=S.INK)


lab = S.ax_mm(fig, 0, TOP, LAB_W, H)
lab.set_xlim(0, LAB_W); lab.set_ylim(len(ROWS) - 0.5, -0.5); S.clean(lab)
for i, (s, name) in enumerate(ROWS):
    lab.add_patch(plt.Rectangle((LAB_W - 2.3, i - 0.22), 1.1, 0.44, fc=S.MODCOL[s], ec="none"))
    lab.text(LAB_W - 3.2, i, name, ha="right", va="center", fontsize=5.8,
             fontweight="bold" if s in ("CAF_ECM", "Immune", "NetImmune") else "normal", color=S.INK)

axes = []
for j, (head, sub, unit, data, xl) in enumerate(COLS):
    ax = S.ax_mm(fig, X0 + j * (COL_W + GAP), TOP, COL_W, H)
    axes.append(ax)
    ax.set_xlim(*xl); ax.set_ylim(len(ROWS) - 0.5, -0.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    S.rule(ax, 0)
    for b in GROUP_BREAKS:
        ax.axhline(b - 0.5, color=S.LIGHT, lw=0.5, zorder=0)
    for i in range(len(ROWS)):
        if i % 2 == 0:
            ax.axhspan(i - 0.5, i + 0.5, color=S.PAPER, zorder=-1, lw=0)
    ax.text(0, 1.075, head, transform=ax.transAxes, ha="left", va="bottom", fontsize=5.7,
            fontweight="bold", color=S.INK, linespacing=1.1)
    ax.set_xlabel(unit, fontsize=5.2, color=S.SLATE, labelpad=1.5)
    ax.tick_params(axis="x", labelsize=5.2)
    if xl[1] < 0.5:
        ax.set_xticks([-0.3, 0, 0.3])
    if j == 3:
        offs = {"PRJEB25780": -0.17, "PRJEB60680": 0.17}
        mk = {"PRJEB25780": "s", "PRJEB60680": "o"}
        for c, d in data.items():
            for i, (s, _) in enumerate(ROWS):
                y, lo, hi, p = d[s]
                col, sig = colour(lo, hi)
                yy = i + offs[c]
                ax.plot([max(lo, xl[0]), min(hi, xl[1])], [yy, yy], color=col, lw=0.7,
                        solid_capstyle="butt")
                ax.scatter([y], [yy], s=8, marker=mk[c], fc=col if sig else "white", ec=col,
                           lw=0.6, zorder=3)
    else:
        for i, (s, _) in enumerate(ROWS):
            v = data[s]
            y, lo, hi = v[0], v[1], v[2]
            col, sig = colour(lo, hi)
            ax.plot([lo, hi], [i, i], color=col, lw=0.9, solid_capstyle="butt")
            ax.scatter([y], [i], s=13, marker="o", fc=col if sig else "white", ec=col, lw=0.7, zorder=3)
            if len(v) == 6:                     # k positive / k
                ax.text(1.0, i, f"{v[4]}/{v[5]}", transform=ax.get_yaxis_transform(), ha="right",
                        va="center", fontsize=5.0, color=S.SLATE)

# a 'k favorable / k' note
axes[2].text(1.0, 1.0, "positive", transform=axes[2].transAxes, ha="right", va="bottom",
             fontsize=5.0, color=S.GREY)
axes[4].text(1.0, 1.0, "favoring", transform=axes[4].transAxes, ha="right", va="bottom",
             fontsize=5.0, color=S.GREY)

# key
key = S.ax_mm(fig, X0, TOP + H + 8.0, 157, 4)
key.set_xlim(0, 157); key.set_ylim(0, 4); S.clean(key)
xx = 0
for fc, ec, txt in [(S.FAV, S.FAV, "95% CI above zero (favorable)"),
                    (S.UNFAV, S.UNFAV, "95% CI below zero (unfavorable)"),
                    ("white", S.GREY, "95% CI includes zero")]:
    key.scatter([xx + 1], [2], s=13, fc=fc, ec=ec, lw=0.7)
    key.text(xx + 2.6, 2, txt, va="center", fontsize=5.0, color=S.SLATE)
    xx += 36
for mk, txt in [("s", "PRJEB25780"), ("o", "NCT04249739")]:
    key.scatter([xx + 1], [2], s=9, marker=mk, fc="white", ec=S.SLATE, lw=0.6)
    key.text(xx + 2.6, 2, txt, va="center", fontsize=5.0, color=S.SLATE)
    xx += 18

# ---------------------------------------------------------------- c. pre-registered expectations
def prereg_rows():
    """Every cell is read from the analyses' own expectation or result files, written when each
    analysis ran (expectations were fixed in the protocol before outcomes were related).
    held = direction / threshold as pre-specified; ci = 95% CI excludes the null (None = n.a.)."""
    fmt = lambda x: f"{x:.2f}".replace("-", "−")
    sg = lambda x: (f"{x:+.2f}").replace("-", "−")
    pct = lambda r: f"{r:.1%}".replace(".0%", "%")
    # E17 independent atlases (primary analysis)
    e17 = pd.read_csv(R / "E17_expectations.csv"); e17 = e17[e17.analysis == "primary"]
    top = e17[e17.expectation.str.contains("highest")]
    def rows17(key):
        return e17[e17.expectation.str.contains(key, regex=False)].set_index("dataset")
    fib, tnk, mac = rows17("fibroblasts contribute"), rows17("T + NK"), rows17("macrophages contribute")
    # E18 PD-L1 CPS
    e18 = pd.read_csv(R / "E18_cps.csv")
    rho = e18[e18.analysis == "Spearman with NetImmune"].iloc[0]
    firth = e18[e18.analysis.str.startswith("Firth: response ~ MSI/EBV + CPS + NetImmune (primary")
                & (e18.term == "NetImmune per SD")].iloc[0]
    # E19 gene loss (main rule, 30% dropout)
    e19 = pd.read_csv(R / "E19_expectations.csv"); e19 = e19[e19.analysis.str.startswith("main")]
    def r19(key):
        r = e19[e19.expectation.str.contains(key, regex=False)]
        return r.iloc[-1]
    rho19 = r19("median Spearman")
    hr19, imm19, caf19 = r19("NetImmune HR < 1"), r19("Immune log-OR > 0"), r19("CAF_ECM log-OR")
    # E20 locked extension (pre-registered exclusions)
    e20 = pd.read_csv(R / "E20_expectations.csv"); e20 = e20[e20.analysis.str.startswith("primary")]
    x1, x2 = e20.iloc[0], e20.iloc[1]
    p20 = pd.read_csv(R / "E20_pooled.csv")
    p20 = p20[(p20.set == "locked 7 (primary)") & (p20.score == "NetImmune") & (p20.model == "unadjusted")
              & (p20.duplicate_rule == "pre-registered exclusions (primary)")].iloc[0]
    n_sig = lambda r: re.search(r"\((\d)/4 with P", r.held).group(1)
    # E21 randomized extension
    m = pd.read_csv(R / "E21_expectations_pooled.csv").set_index("score")
    pt = pd.read_csv(R / "E21_expectations_per_trial.csv")
    prim = pt[pt.role == "primary"]
    cmos = pt[(pt.trial == "CheckMate 025") & (pt.endpoint == "OS")]
    dut = pt[(pt.trial == "DUTRENEO") & (pt.role == "sensitivity")]
    ci21 = lambda sc: f"{sg(m.loc[sc].est)} ({fmt(m.loc[sc].lo)} to {fmt(m.loc[sc].hi)})"
    # E22 gastro-esophageal
    e22 = pd.read_csv(R / "E22_expectations.csv").set_index("item")
    g = lambda k: e22.loc[k]
    grp = lambda prefix: e22[e22.index.str.startswith(prefix)]
    lit = grp("PRJEB60680 [literal]")
    lit = lit[~lit.index.str.contains("MSS")]
    pooled, perf, gse = grp("PRJEB60680 + PRJEB25780 (primary)"), grp("PRJEB60680 + PRJEB25780 + PERFECT"), grp("GSE183924")
    allheld = lambda d: bool(d.direction_held.all())
    cnt = lambda d: f"{int(d.direction_held.sum())}/{len(d)}"
    mss = g("PRJEB60680 MSS/EBV-negative NetImmune log-OR")
    return [
        ("Independent single-cell atlases", [
            ("Each module highest in its intended class", f"{int((top.held == 'held').sum())}/{len(top)} module × atlas",
             bool((top.held == "held").all()), None),
            ("Fibroblasts ≥ 50% of CAF_ECM signal, OMIX001073", fib.loc["OMIX001073"].result, fib.loc["OMIX001073"].held == "held", None),
            ("Fibroblasts ≥ 50% of CAF_ECM signal, GSE206785", fib.loc["GSE206785"].result, fib.loc["GSE206785"].held == "held", None),
            ("T and NK cells ≥ 50% of T_CELL signal", "; ".join(tnk.result), bool((tnk.held == "held").all()), None),
            ("Macrophages < 50% of TAM_SUP signal", "; ".join(mac.result), bool((mac.held == "held").all()), None)]),
        ("PD-L1 CPS, PRJEB25780", [
            ("NetImmune correlates positively with CPS", f"ρ {rho.estimate:.2f}", rho.estimate > 0, rho.p < 0.05),
            ("NetImmune > 0 beyond MSI/EBV and CPS", f"OR {firth.OR:.1f}, P = {firth.p:.2f}", firth.estimate > 0, firth.p < 0.05)]),
        ("Random gene loss (30% of genes)", [
            ("Correlation with full score ≥ 0.9", f"median ρ {rho19.value:.2f}", bool(rho19.held), None),
            ("Pooled NetImmune HR < 1 in ≥ 95% of draws", f"{pct(hr19.value)} of draws", bool(hr19.held), None),
            ("Pooled Immune log-odds > 0 in ≥ 95% of draws", f"{pct(imm19.value)} of draws", bool(imm19.held), None),
            ("Pooled CAF_ECM log-odds ≤ 0.1 in ≥ 95% of draws", f"{pct(caf19.value)} of draws", bool(caf19.held), None)]),
        ("Locked gastric extension", [
            ("NetImmune HR < 1 in each added cohort", f"{x1.held.split()[0]}; {n_sig(x1)} with P < 0.05", x1.held.startswith("4/4"), n_sig(x1) == "4"),
            ("CAF_ECM HR > 1 in each added cohort", f"{x2.held.split()[0]}; {n_sig(x2)} with P < 0.05", x2.held.startswith("4/4"), n_sig(x2) == "4"),
            ("Pooled NetImmune HR < 1, seven locked cohorts", f"{p20.HR:.2f} ({p20.lo:.2f}–{p20.hi:.2f})", p20.HR < 1, p20.hi < 1)]),
        ("Randomized extension (six comparisons)", [
            ("Pooled NetImmune interaction > 0", ci21("NetImmune"), bool(m.loc["NetImmune"].held), bool(m.loc["NetImmune"].supported_conservative)),
            ("Pooled Immune interaction > 0", ci21("Immune"), bool(m.loc["Immune"].held), bool(m.loc["Immune"].supported_conservative)),
            ("Pooled IFN interaction > 0", ci21("IFN"), bool(m.loc["IFN"].held), bool(m.loc["IFN"].supported_conservative)),
            ("Pooled CAF_ECM interaction ≤ 0", ci21("CAF_ECM"), bool(m.loc["CAF_ECM"].held), bool(m.loc["CAF_ECM"].supported_conservative)),
            ("Each added trial, four scores", f"{int(prim.held.sum())}/{len(prim)}; CI excludes 0 in {int(prim.supported.sum())}", bool(prim.held.all()), bool(prim.supported.all())),
            ("CheckMate 025 overall survival (secondary)", f"{int(cmos.held.sum())}/{len(cmos)}", bool(cmos.held.all()), bool(cmos.supported.all())),
            ("DUTRENEO (sensitivity)", f"{int(dut.held.sum())}/{len(dut)}", bool(dut.held.all()), bool(dut.supported.all()))]),
        ("Gastroesophageal checkpoint cohorts", [
            ("NCT04249739: Immune and IFN log-odds > 0",
             f"{sg(g('PRJEB60680 Immune log-OR').estimate)}; {sg(g('PRJEB60680 IFN log-OR').estimate)}",
             bool(g("PRJEB60680 Immune log-OR").direction_held and g("PRJEB60680 IFN log-OR").direction_held),
             bool(g("PRJEB60680 Immune log-OR").ci_excludes_null and g("PRJEB60680 IFN log-OR").ci_excludes_null)),
            ("NCT04249739: NetImmune log-odds > 0", sg(g("PRJEB60680 NetImmune log-OR").estimate),
             bool(g("PRJEB60680 NetImmune log-OR").direction_held), bool(g("PRJEB60680 NetImmune log-OR").ci_excludes_null)),
            ("NCT04249739: CAF_ECM log-odds ≤ 0", sg(g("PRJEB60680 CAF_ECM log-OR").estimate),
             bool(g("PRJEB60680 CAF_ECM log-OR").direction_held), bool(g("PRJEB60680 CAF_ECM log-OR").ci_excludes_null)),
            ("NCT04249739, MSS/EBV-negative: NetImmune > 0", sg(mss.estimate), bool(mss.direction_held), bool(mss.ci_excludes_null)),
            ("NCT04249739, pre-specified processing", f"{cnt(lit)} scores", allheld(lit), bool(lit.ci_excludes_null.all())),
            ("Two-cohort pooled log-odds, four scores", f"{cnt(pooled)}; CIs include 0" if not pooled.ci_excludes_null.any() else cnt(pooled),
             allheld(pooled), bool(pooled.ci_excludes_null.all())),
            ("With PERFECT (sensitivity), four scores", cnt(perf), allheld(perf), bool(perf.ci_excludes_null.all())),
            ("GSE183924: relapse HR as expected, four scores", cnt(gse), allheld(gse), bool(gse.ci_excludes_null.all()))]),
    ]


ROWS_B = prereg_rows()
SC_Y = TOP + H + 17                     # panel b (scorecard)
S.letter_at(fig, 1, SC_Y - 1.0, "b")
fig.text(7 / FIG_W, 1 - (SC_Y - 0.4) / FIG_H, "Pre-registered expectations",
         fontsize=6.5, fontweight="bold", va="top", ha="left", color=S.INK)
LH = 3.55
cols = [ROWS_B[:4], ROWS_B[4:]]
heights = [sum(1 + len(r) for _, r in c) * LH for c in cols]
HB = max(heights) + 6
for ci_, col in enumerate(cols):
    x0 = 0 if ci_ == 0 else 91
    ax = S.ax_mm(fig, x0, SC_Y + 3.0, 89, HB)
    ax.set_xlim(0, 89); ax.set_ylim(HB, 0); S.clean(ax)
    ax.text(77.5, 1.2, "held", ha="center", va="center", fontsize=5.0, color=S.SLATE)
    ax.text(85.5, 1.2, "CI", ha="center", va="center", fontsize=5.0, color=S.SLATE)
    y = 4.2
    for head, rows in col:
        ax.text(1.5, y, head, fontsize=5.6, fontweight="bold", color=S.INK, va="center")
        y += LH
        for k, (stmt, res, d, c) in enumerate(rows):
            if k % 2 == 0:
                ax.add_patch(plt.Rectangle((0.5, y - LH / 2), 88, LH, fc=S.PAPER, ec="none", zorder=0))
            ax.text(2.5, y, stmt, fontsize=5.0, color=S.INK, va="center")
            ax.text(73.0, y, res, fontsize=5.0, color=S.SLATE, va="center", ha="right")
            ax.scatter([77.5], [y], s=12, marker="o", fc=S.INK if d else "white", ec=S.INK if d else S.UNFAV,
                       lw=0.7, zorder=3)
            if c is None:
                ax.text(85.5, y, "–", fontsize=5.0, color=S.GREY, ha="center", va="center")
            else:
                ax.scatter([85.5], [y], s=10, marker="s", fc=S.INK if c else "white", ec=S.INK, lw=0.6, zorder=3)
            y += LH
leg = S.ax_mm(fig, 0, SC_Y + 3.0 + HB + 0.5, 180, 3.5)
leg.set_xlim(0, 180); leg.set_ylim(0, 3.5); S.clean(leg)
xx = 2
for mk, fc, ec, txt in [("o", S.INK, S.INK, "held as pre-specified"), ("o", "white", S.UNFAV, "not held"),
                        ("s", S.INK, S.INK, "95% CI excludes the null"), ("s", "white", S.INK, "95% CI includes the null")]:
    leg.scatter([xx + 1], [1.75], s=11, marker=mk, fc=fc, ec=ec, lw=0.7)
    leg.text(xx + 2.6, 1.75, txt, va="center", fontsize=5.0, color=S.SLATE)
    xx += 40

print(S.save(fig, "Fig8"))
