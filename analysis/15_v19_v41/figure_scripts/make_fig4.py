"""v20 Figure 4 - gastric checkpoint cohorts: immune modules lean towards response, as in other cancers.

(a,b) ROC curves in PRJEB25780 (pembrolizumab; 44 patients, 12 responders) and PRJEB60680 / NCT04249739
      (5-FU/platinum, then + pembrolizumab; 43 patients, 29 CR/PR), same four markers in both.
(c)   AUC of clinical markers and scores in the two cohorts.
(d)   Log-odds of response per SD in each cohort against the pooled estimate of the 21 checkpoint cohorts of
      other cancers (grey band; E23, post hoc), with the NCT04249739 consistency P and power, and GSE183924
      relapse-free survival (adjuvant durvalumab), on shared rows.
(e)   Score distributions by response, marker shape for MSI-H/EBV+.
All estimates are read from analysis/results (E7, E18, E22); ROC curves are descriptive re-plots of the score
tables."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from sklearn.metrics import roc_curve, roc_auc_score
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, MACRO, fmt_p

A = S.RES
FS = 5.2
NUM = {}
FIGH = 221.0


LET = []           # (letter, title artist, x in mm): letters are placed after layout, top level with the title


def place_letters(fig):
    """Letter ascender line level with the cap height of the panel's bold title (glyph metrics, not text boxes)."""
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    fig.canvas.draw()
    r = fig.canvas.get_renderer(); k = fig.dpi / 72.0
    lp = FontProperties(family=plt.rcParams["font.family"], size=8, weight="bold")
    asc = TextPath((0, 0), "l", prop=lp).get_extents().y1
    for s, art, x in LET:
        prop = art.get_fontproperties()
        d = max(r.get_text_width_height_descent("lp", prop, ismath=False)[2],
                r.get_text_width_height_descent(art.get_text(), prop, ismath=False)[2])
        base = art.get_window_extent(r).y0 + d
        cap = TextPath((0, 0), "H", prop=prop).get_extents().y1
        y = (base + (cap - asc) * k) / fig.bbox.height
        fig.text(x / 180.0, y, s, fontsize=8, fontweight="bold", va="baseline", ha="left", color=INK)
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import Affine2D


def _mm(fig):
    return Affine2D().scale(S.MM) + fig.dpi_scale_trans          # millimetres from the bottom-left corner


def _ymm(fig, ax, y):
    bb = ax.get_position(); y0, y1 = ax.get_ylim()
    return (bb.y0 + (y - y0) / (y1 - y0) * bb.height) * FIGH


def chip(fig, x_mm, y_mm, col):
    """Family-colour chip (1.1 x 1.6 mm rounded bar), left edge x_mm, centered on y_mm (mm from bottom)."""
    fig.add_artist(FancyBboxPatch((x_mm, y_mm - 0.8), 1.1, 1.6, boxstyle="round,pad=0,rounding_size=0.3",
                                  transform=_mm(fig), fc=col, ec="none"))


def row_label(fig, ax, y, text, col, gap=1.2, fs=FS):
    """Ink label right-aligned left of the axes, with the family chip between label and axes."""
    xr = ax.get_position().x0 * 180.0 - gap
    ym = _ymm(fig, ax, y)
    if col:
        chip(fig, xr - 1.1, ym, col)
    fig.text((xr - 2.0) / 180.0, ym / FIGH, text, ha="right", va="center", fontsize=fs, color=INK)


def chip_title(fig, ax, text, col, dy_mm=3.6, fs=5.6):
    """Chip + bold ink name, left-aligned above the axes."""
    bb = ax.get_position(); x0 = bb.x0 * 180.0; ym = bb.y1 * FIGH + dy_mm
    chip(fig, x0, ym, col)
    return fig.text((x0 + 2.0) / 180.0, ym / FIGH, text, ha="left", va="center", fontsize=fs, fontweight="bold",
                    color=INK)


def rec(k, v, s):
    NUM[k] = (v, s)


def neg(x, fmt="{:+.2f}"):
    return fmt.format(x).replace("-", "−")


K = pd.read_csv(A / "E7_scores_log2.csv")                       # PRJEB25780
K["y"] = K.Responder.astype(int); K["hot"] = K.MSI_EBV.astype(int)
SM = pd.read_csv(A / "E22_scores_PRJEB60680.csv")                # PRJEB60680
SM["hot"] = SM.msi_ebv.astype(int)
R60 = pd.read_csv(A / "E22_PRJEB60680_response.csv")
R60 = R60[R60.analysis == "primary: all QC-pass baseline, z within library batch"].set_index("score")
IN = pd.read_csv(A / "E22_meta_inputs.csv")
P21 = pd.read_csv(A / "E23_pooled_with_gastric.csv").set_index("score")          # post hoc (E23)
PW = pd.read_csv(A / "E23_power_consistency.csv"); PW = PW[PW.cohort == "PRJEB60680"].set_index("score")
AU = pd.read_csv(A / "E7_auc_ci_log2.csv").set_index("Marker")
E18 = pd.read_csv(A / "E18_cps.csv")
FC = pd.read_csv(A / "E22_firth_cps.csv")
FC = FC[FC.set == "all QC-pass baseline (n = 43)"]
C2 = pd.read_csv(A / "E22_GSE183924_cox.csv").set_index("score")

COH = [("PRJEB25780", "PRJEB25780", "pembrolizumab", K, "Responder"),
       ("PRJEB60680", "NCT04249739 (PRJEB60680)", "chemotherapy, then + pembrolizumab", SM, "y")]
fig = plt.figure(figsize=(S.W2, FIGH * S.MM))

# ------------------------------------------------------------------ a, b  ROC
CURVES = [("NetImmune", "NetImmune", IMMUNE, "-", 1.1), ("IFN", "IFN", IMMUNE, (0, (2.5, 1.2)), 0.8),
          ("CAF_ECM", "CAF_ECM", STROMA, "-", 0.8), ("hot", "MSI-H/EBV+", GREY, "-", 0.8)]
for i, (key, name, tx, D, ycol) in enumerate(COH):
    ax = S.ax_mm(fig, 13 + i * 60, 10, 33, 33)
    yv = D[ycol].astype(int).values
    for k, lab, col, ls, lw in CURVES:
        f, t, _ = roc_curve(yv, D[k].values)
        auc = roc_auc_score(yv, D[k].values)
        ax.plot(f, t, color=col, ls=ls, lw=lw, drawstyle="steps-post" if k != "hot" else "default", zorder=3)
        ax.plot([], [], color=col, ls=ls, lw=lw, label=f"{lab}  {auc:.2f}")
        rec(f"{'ab'[i]} {key} AUC {lab}", f"{auc:.3f}", f"recomputed from {'E7_scores_log2' if i == 0 else 'E22_scores_PRJEB60680'}.csv; "
            f"matches {'E7_auc_ci_log2 / E22_meta_inputs' if i == 0 else 'E22_PRJEB60680_response'}")
    ax.plot([0, 1], [0, 1], color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.set_xlabel("1 − specificity", labelpad=1)
    if i == 0:
        ax.set_ylabel("Sensitivity", labelpad=1)
    lg = ax.legend(loc="lower left", bbox_to_anchor=(1.01, -0.02), fontsize=5.0, handlelength=1.5,
                   handletextpad=0.4, borderaxespad=0.0, labelspacing=0.35, title="AUC", title_fontsize=5.0)
    lg._legend_box.align = "left"
    ax.set_title(name, loc="left", fontsize=6.2, fontweight="bold", pad=8)

    LET.append(("ab"[i], ax._left_title, 1.0 if i == 0 else 65.0))

# ------------------------------------------------------------------ c  AUC against clinical markers
axC = S.ax_mm(fig, 141, 10, 35, 33)
cps25 = E18[(E18.analysis == "AUC (n = 40)") & (E18.term == "log(CPS+1)")].iloc[0]
cps60 = FC[(FC.analysis == "AUC (n = 43)") & (FC.term == "log(CPS+1)")].iloc[0]
msi60 = FC[(FC.analysis == "AUC (n = 43)") & (FC.term == "MSI/EBV")].iloc[0]
CR = [("PD-L1 CPS*", (cps25.estimate, cps25.lo, cps25.hi), (cps60.estimate, cps60.lo, cps60.hi), INK),
      ("MSI-H/EBV+", tuple(AU.loc["MSI/EBV", ["AUC", "CI_low", "CI_high"]]), (msi60.estimate, msi60.lo, msi60.hi), INK),
      ("IFN", tuple(AU.loc["IFN", ["AUC", "CI_low", "CI_high"]]), tuple(R60.loc["IFN", ["AUC", "AUC_lo", "AUC_hi"]]),
       IMMUNE),
      ("Immune", tuple(AU.loc["Immune", ["AUC", "CI_low", "CI_high"]]),
       tuple(R60.loc["Immune", ["AUC", "AUC_lo", "AUC_hi"]]), IMMUNE),
      ("NetImmune", tuple(AU.loc["NetImmune", ["AUC", "CI_low", "CI_high"]]),
       tuple(R60.loc["NetImmune", ["AUC", "AUC_lo", "AUC_hi"]]), IMMUNE)]
for j, (lab, a, b, col) in enumerate(CR):
    y = j
    axC.plot([a[0], b[0]], [y, y], color=LIGHT, lw=2.2, solid_capstyle="round", zorder=1)
    axC.plot([a[1], a[2]], [y - 0.17, y - 0.17], color=col, lw=0.6, zorder=2)
    axC.plot([b[1], b[2]], [y + 0.17, y + 0.17], color=col, lw=0.6, zorder=2)
    axC.plot(a[0], y - 0.17, "o", ms=3.2, color=col, mec="white", mew=0.4, zorder=3)
    axC.plot(b[0], y + 0.17, "o", ms=3.2, mfc="white", mec=col, mew=0.8, zorder=3)
    rec(f"c {lab} PRJEB25780", f"{a[0]:.3f} ({a[1]:.3f}–{a[2]:.3f})", "E18_cps (CPS, n = 40) / E7_auc_ci_log2")
    rec(f"c {lab} PRJEB60680", f"{b[0]:.3f} ({b[1]:.3f}–{b[2]:.3f})", "E22_firth_cps (CPS, MSI/EBV) / E22_PRJEB60680_response")
S.rule(axC, 0.5)
axC.set_yticks([])
axC.tick_params(axis="y", length=0); axC.spines["left"].set_visible(False)
axC.set_ylim(len(CR) - 0.45, -0.75)
for j, c in enumerate(CR):
    row_label(fig, axC, j, c[0], c[3] if c[3] != INK else None)
axC.set_xlim(0.2, 1.02); axC.set_xticks([0.25, 0.5, 0.75, 1]); axC.set_xticklabels(["0.25", "0.5", "0.75", "1"])
axC.set_xlabel("AUC (95% CI)", labelpad=1)
tC = axC.text(-0.33, 1.13, "Markers in both cohorts", transform=axC.transAxes, fontsize=6.2, fontweight="bold", va="bottom")
axC.plot(-0.30, 1.045, "o", ms=3.0, color=INK, mec="white", mew=0.3, transform=axC.transAxes, clip_on=False)
axC.text(-0.26, 1.015, "PRJEB25780", transform=axC.transAxes, fontsize=5.0, color=SLATE, va="bottom")
axC.plot(0.36, 1.045, "o", ms=3.0, mfc="white", mec=INK, mew=0.7, transform=axC.transAxes, clip_on=False)
axC.text(0.40, 1.015, "NCT04249739", transform=axC.transAxes, fontsize=5.0, color=SLATE, va="bottom")
LET.append(("c", tC, 122.0))

# ------------------------------------------------------------------ d  log-odds per SD, pooled, RFS
ROWS = [("CAF_ECM", STROMA), ("TAM_SUP", MACRO), ("APM", IMMUNE), ("IFN", IMMUNE), ("TLS_B", IMMUNE),
        ("T_CELL", IMMUNE), ("Immune", IMMUNE), ("NetImmune", IMMUNE)]
NR = len(ROWS)
TOPD, HD = 63.0, 44.0
LIM = (-2.2, 3.4)


def point(ax, v, lo, hi, y, c, mk="o", ms=3.0, open_=False, lim=LIM, log=False):
    a, b = lim
    if lo < a:
        ax.plot(a, y, "<", ms=2.2, color=c, mew=0, clip_on=False); lo = a
    if hi > b:
        ax.plot(b, y, ">", ms=2.2, color=c, mew=0, clip_on=False); hi = b
    ax.plot([lo, hi], [y, y], color=c, lw=0.8, solid_capstyle="butt", zorder=2)
    ax.plot(v, y, mk, ms=ms, color=c, mfc="white" if open_ else c, mec=c if open_ else "white",
            mew=0.8 if open_ else 0.4, zorder=3)


def rows_axis(ax, labels):
    ax.set_ylim(NR - 0.45, -0.6)
    ax.set_yticks([])
    if labels:
        for i, (k, c) in enumerate(ROWS):
            row_label(fig, ax, i, k, c)
    ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
    for yy in (0.5, 1.5):
        ax.axhline(yy, color=LIGHT, lw=0.4, zorder=0)


def head(ax, t, sub):
    ax.set_title(t, loc="left", fontsize=6.0, fontweight="bold", pad=8)



XS = [(22, "PRJEB25780"), (65, "PRJEB60680")]
AXD = []
for j, (x0, coh) in enumerate(XS):
    ax = S.ax_mm(fig, x0, TOPD, 37, HD)
    AXD.append(ax)
    T = IN[IN.cohort == coh].set_index("score")
    for i, (k, c) in enumerate(ROWS):
        q = P21.loc[k]                                   # 21 checkpoint cohorts of other cancers, pooled
        ax.add_patch(mpatches.Rectangle((q.pan21_lo, i - 0.30), q.pan21_hi - q.pan21_lo, 0.60, facecolor=LIGHT,
                                        edgecolor="none", zorder=0.5))
        ax.plot([q.pan21_est] * 2, [i - 0.30, i + 0.30], color=GREY, lw=0.9, solid_capstyle="butt", zorder=1)
        if j == 0:
            rec(f"d 21-cohort pooled {k}", f"{q.pan21_est:+.2f} ({q.pan21_lo:+.2f} to {q.pan21_hi:+.2f})",
                "E23_pooled_with_gastric.csv (pan21)")
        r = T.loc[k]; lo, hi = r.y - 1.96 * r.se, r.y + 1.96 * r.se
        point(ax, r.y, lo, hi, i, c)
        if r.p_cohort < 0.05:
            ax.plot(r.y, i, "o", ms=5.6, mfc="none", mec=INK, mew=0.5, zorder=4)
        rec(f"d {coh} {k} log-OR", f"{r.y:+.2f} ({lo:+.2f} to {hi:+.2f}), P {fmt_p(r.p_cohort, '')}",
            "E22_meta_inputs.csv")
    S.rule(ax, 0)
    ax.set_xlim(*LIM); ax.set_xticks([-2, -1, 0, 1, 2, 3]); ax.set_xticklabels(["−2", "−1", "0", "1", "2", "3"])
    rows_axis(ax, j == 0)
    n = T.iloc[0]
    head(ax, "PRJEB25780" if j == 0 else "NCT04249739", f"{int(n.n)} patients, {int(n.n_resp)} CR/PR")
    ax.set_xlabel("log-odds of response per SD", labelpad=1)
    if j == 0:
        LET.append(("d", ax._left_title, 1.0))
# NCT04249739 against the 21 other cohorts (post hoc, E23): P for a difference and power to detect their effect
trP = ax.get_yaxis_transform()
for i, (k, c) in enumerate(ROWS):
    r = PW.loc[k]
    ax.text(1.17, i, fmt_p(r.p_differs_from_pan21, ""), transform=trP, fontsize=5.0, va="center", ha="right",
            color=INK, clip_on=False)
    ax.text(1.40, i, f"{100 * r.power_at_pan21:.0f}%", transform=trP, fontsize=5.0, va="center", ha="right",
            color=INK, clip_on=False)
    rec(f"d NCT04249739 vs 21 {k}", f"P {fmt_p(r.p_differs_from_pan21, '')}; power {100 * r.power_at_pan21:.0f}%; "
        f"n for 80% {int(r.n_for_80pct)}", "E23_power_consistency.csv")
ax.text(1.17, 1.015, "P vs\n21", transform=ax.transAxes, fontsize=5.0, color=SLATE, va="bottom", ha="right",
        linespacing=1.1)
ax.text(1.40, 1.015, "power", transform=ax.transAxes, fontsize=5.0, color=SLATE, va="bottom", ha="right")


# GSE183924 relapse-free survival
axR = S.ax_mm(fig, 134, TOPD, 36, HD)
RL = (0.3, 3.0)
for i, (k, c) in enumerate(ROWS):
    r = C2.loc[k]
    point(axR, r.HR, r.lo, r.hi, i, c, mk="s", ms=2.8, lim=RL)
    axR.text(1.04, i, fmt_p(r.p, ""), transform=axR.get_yaxis_transform(), fontsize=5.0, va="center", ha="left",
             color=INK, clip_on=False)
    rec(f"d GSE183924 {k} HR", f"{r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f}), P {fmt_p(r.p, '')}", "E22_GSE183924_cox.csv")
S.rule(axR, 1)
axR.set_xscale("log"); axR.set_xlim(*RL); axR.set_xticks([0.5, 1, 2]); axR.set_xticklabels(["0.5", "1", "2"])
axR.xaxis.set_minor_formatter(NullFormatter())
rows_axis(axR, False)
r0 = C2.loc["NetImmune"]
head(axR, "GSE183924, relapse-free", f"adjuvant durvalumab, {int(r0.n)} patients")
axR.text(1.04, 1.015, "P", transform=axR.transAxes, fontsize=5.0, color=SLATE, va="bottom", ha="left")
axR.set_xlabel("HR per SD", labelpad=1)
LET.append(("e", axR._left_title, 126.0))


# ------------------------------------------------------------------ row 3: f CPS, g null, h MSS/EBV-negative, i PFS
TOP3, H3 = 127.0, 30.0
# f  NetImmune against PD-L1 CPS, PRJEB25780
CP = pd.read_csv(A / "E18_cps_patients.csv")
rho = E18[(E18.analysis == "Spearman with NetImmune") & (E18.term == "log(CPS+1)")].iloc[0]
rho60 = FC[(FC.analysis == "Spearman with NetImmune") & (FC.term == "log(CPS+1)")].iloc[0]
axF3 = S.ax_mm(fig, 14, TOP3, 30, H3)
CP = CP.assign(lx=np.log1p(CP.cps))
ties = CP.lx.map(CP.lx.value_counts()) > 1                       # spread tied CPS values (mostly CPS = 0)
CP["jx"] = CP.lx + np.where(ties, np.random.default_rng(11).uniform(-0.22, 0.22, len(CP)), 0.0)
for resp, hot in [(0, 0), (0, 1), (1, 0), (1, 1)]:
    d = CP[(CP.Responder == resp) & (CP.MSI_EBV == hot)]
    mk = "^" if hot else "o"
    if resp:
        axF3.plot(d.jx, d.NetImmune, mk, ms=2.9 if hot else 2.5, color=INK, mec="white", mew=0.3, zorder=3)
    else:
        axF3.plot(d.jx, d.NetImmune, mk, ms=2.9 if hot else 2.5, mfc="white", mec=GREY if not hot else INK,
                  mew=0.6, zorder=4 if hot else 2)
axF3.set_xlabel("log(CPS + 1)", labelpad=1); axF3.set_ylabel("NetImmune (z)", labelpad=1)
axF3.set_xlim(-0.45, 5.0); axF3.set_xticks([0, 2, 4]); axF3.set_ylim(-2.6, 2.6); axF3.set_yticks([-2, 0, 2])
axF3.set_yticklabels(["−2", "0", "2"])
axF3.text(0.03, 0.97, f"ρ = {rho.estimate:.2f}, {fmt_p(rho.p)}", transform=axF3.transAxes, fontsize=5.0, va="top",
          color=INK)
axF3.set_title("NetImmune and PD-L1 CPS", loc="left", fontsize=6.0, fontweight="bold", pad=8)

rec("f Spearman NetImmune vs log(CPS+1), PRJEB25780", f"rho {rho.estimate:.2f}, P {fmt_p(rho.p, '')}, n {int(rho.n)}",
    "E18_cps.csv")
rec("f (legend) Spearman, NCT04249739", f"rho {rho60.estimate:+.2f}, P {fmt_p(rho60.p, '')}, n {int(rho60.n)}",
    "E22_firth_cps.csv")
for x0, mk, kw, t in [(0.0, "o", dict(color=INK, mec="white", mew=0.3), "CR/PR"),
                      (0.33, "o", dict(mfc="white", mec=GREY, mew=0.6), "SD/PD")]:
    axF3.plot(x0 + 0.03, -0.30, mk, ms=2.5, transform=axF3.transAxes, clip_on=False, **kw)
    axF3.text(x0 + 0.08, -0.30, t, transform=axF3.transAxes, fontsize=5.0, color=SLATE, va="center")
axF3.plot(0.66, -0.30, "^", ms=2.9, color=INK, mec="white", mew=0.3, transform=axF3.transAxes, clip_on=False)
axF3.plot(0.71, -0.30, "^", ms=2.9, mfc="white", mec=INK, mew=0.6, transform=axF3.transAxes, clip_on=False)
axF3.text(0.76, -0.30, "MSI/EBV+", transform=axF3.transAxes, fontsize=5.0, color=SLATE, va="center")
LET.append(("f", axF3._left_title, 1.0))

# g  structure-matched random-gene null for NetImmune, both cohorts
N25 = pd.read_csv(A / "E3_ici_matched_null.csv"); N25 = N25[(N25.scale == "log2") & (N25.score == "NetImmune")].iloc[0]
D25 = pd.read_csv(A / "E3_ici_null_draws_log2.csv").NetImmune.values
N60 = pd.read_csv(A / "E22_PRJEB60680_null.csv").set_index("score").loc["NetImmune"]
D60 = pd.read_csv(A / "E22_PRJEB60680_null_draws.csv").NetImmune.values
axG = S.ax_mm(fig, 57, TOP3, 33, H3)
bins = np.linspace(0.45, 0.95, 41)
for i, (lab, dr, ob, p95, pe) in enumerate([("PRJEB25780", D25, N25.AUC_obs, N25.null_p95, N25.p_emp),
                                            ("NCT04249739", D60, N60.AUC_obs, N60.null_p95, N60.p_emp)]):
    yb = 1 - i
    h, e = np.histogram(dr, bins=bins); h = h / h.max() * 0.55
    axG.fill_between(e[:-1], yb, yb + h, step="post", color=LIGHT, lw=0)
    axG.plot([p95, p95], [yb, yb + 0.6], color=GREY, lw=0.6)
    axG.plot(ob, yb + 0.12, "D", ms=3.2, color=IMMUNE, mec="white", mew=0.4, zorder=4)
    axG.text(0.45, yb + 0.86, lab, fontsize=5.2, fontweight="bold", va="center", ha="left")
    axG.text(0.95, yb + 0.86, f"AUC {ob:.2f}, {fmt_p(pe)}", fontsize=5.0, va="center", ha="right", color=INK)
    rec(f"g {lab} NetImmune null", f"AUC {ob:.3f}; null 95th pct {p95:.3f}; P {fmt_p(pe, '')}; {len(dr)} draws",
        "E3_ici_matched_null / E22_PRJEB60680_null")
axG.set_ylim(-0.05, 2.0); axG.set_yticks([]); axG.spines["left"].set_visible(False)
axG.set_xlim(0.45, 0.95); axG.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9]); axG.set_xticklabels(["0.5", "", "0.7", "", "0.9"])
axG.set_xlabel("AUC, NetImmune ◇ vs random axes", labelpad=1)
axG.set_xlabel("AUC", labelpad=1)
axG.set_title("Random-gene null", loc="left", fontsize=6.0, fontweight="bold", pad=8)

LET.append(("g", axG._left_title, 50.0))

# h  all patients against MSS/EBV-negative patients
SG25 = pd.read_csv(A / "E7_subgroup_auc.csv").set_index("score")
SG60 = R60_MSS = pd.read_csv(A / "E22_PRJEB60680_response.csv")
SG60 = SG60[SG60.analysis == "MSS/EBV-negative (pre-registered subgroup)"].set_index("score")
axH = S.ax_mm(fig, 116, TOP3, 25, H3)
HR_ = [("PRJEB25780", None), ("NetImmune", "25"), ("Immune", "25"), ("IFN", "25"),
       ("NCT04249739", None), ("NetImmune", "60"), ("Immune", "60"), ("IFN", "60")]
for j, (k, c) in enumerate(HR_):
    if c is None:
        continue
    if c == "25":
        a = tuple(AU.loc[k, ["AUC", "CI_low", "CI_high"]]); b = tuple(SG25.loc[k, ["AUC", "lo", "hi"]])
        nb = f"{int(SG25.loc[k, 'n'])} patients, {int(SG25.loc[k, 'responders'])} CR/PR"
    else:
        a = tuple(R60.loc[k, ["AUC", "AUC_lo", "AUC_hi"]]); b = tuple(SG60.loc[k, ["AUC", "AUC_lo", "AUC_hi"]])
        nb = f"{int(SG60.loc[k, 'n'])} patients, {int(SG60.loc[k, 'n_resp'])} CR/PR"
    axH.plot([a[1], a[2]], [j - 0.17] * 2, color=IMMUNE, lw=0.6)
    axH.plot(a[0], j - 0.17, "o", ms=2.8, color=IMMUNE, mec="white", mew=0.3, zorder=3)
    if c == "25":                                   # 2 responders: bootstrap CI degenerate, drawn dotted
        axH.plot([b[1], b[2]], [j + 0.17] * 2, color=IMMUNE, lw=0.6, ls=(0, (1, 1.2)))
    else:
        axH.plot([b[1], b[2]], [j + 0.17] * 2, color=IMMUNE, lw=0.6)
    axH.plot(b[0], j + 0.17, "o", ms=2.8, mfc="white", mec=IMMUNE, mew=0.7, zorder=3)
    rec(f"h {'PRJEB25780' if c == '25' else 'NCT04249739'} {k}",
        f"all {a[0]:.3f} ({a[1]:.3f}–{a[2]:.3f}); MSS/EBV-negative {b[0]:.3f} ({b[1]:.3f}–{b[2]:.3f}), {nb}",
        "E7_auc_ci_log2 + E7_subgroup_auc / E22_PRJEB60680_response")
axH.set_ylim(len(HR_) - 0.5, -0.6); axH.set_yticks([]); axH.spines["left"].set_visible(False)
S.rule(axH, 0.5)
axH.set_xlim(0.2, 1.02); axH.set_xticks([0.25, 0.5, 0.75, 1]); axH.set_xticklabels(["0.25", "0.5", "0.75", "1"])
axH.set_xlabel("AUC (95% CI)", labelpad=1)
for j, (k, c) in enumerate(HR_):
    if c is None:
        ym = _ymm(fig, axH, j)
        fig.text((116 - 15.5) / 180.0, ym / FIGH, k, fontsize=5.0, fontweight="bold", ha="left", va="center",
                 color=INK)
    else:
        row_label(fig, axH, j, k, IMMUNE, fs=5.0)
axH.text(1.0, 0, "MSS/EBV−: 2 responders", fontsize=5.0, color=SLATE, ha="right", va="center",
         transform=axH.get_yaxis_transform(), bbox=dict(fc="white", ec="none", pad=0.4), zorder=5)
axH.set_title("MSS/EBV-negative tumors", loc="left", fontsize=6.0, fontweight="bold", pad=8, x=-0.62)
axH.plot(-0.58, 1.045, "o", ms=2.6, color=IMMUNE, mec="white", mew=0.3, transform=axH.transAxes, clip_on=False)
axH.text(-0.52, 1.015, "all", transform=axH.transAxes, fontsize=5.0, color=SLATE, va="bottom")
axH.plot(-0.25, 1.045, "o", ms=2.6, mfc="white", mec=IMMUNE, mew=0.7, transform=axH.transAxes, clip_on=False)
axH.text(-0.19, 1.015, "MSS/EBV−", transform=axH.transAxes, fontsize=5.0, color=SLATE, va="bottom")
LET.append(("h", axH._left_title, 96.0))

# i  progression-free survival by the Immune axis, NCT04249739 (not pre-specified)
from lifelines import KaplanMeierFitter
CX = pd.read_csv(A / "E22_PRJEB60680_cox.csv")
cx = CX[(CX.endpoint == "PFS") & (CX.score == "Immune")].iloc[0]
axI = S.ax_mm(fig, 150, TOP3, 26, H3)
m = SM.pfs_d / 30.44
hi = SM.Immune > SM.Immune.median()
for grp, col, lab in [(hi, IMMUNE, "high"), (~hi, "#7A8794", "low")]:
    k = KaplanMeierFitter().fit(m[grp], SM.pfs_ev[grp])
    sf = k.survival_function_
    axI.step(sf.index.values, sf.iloc[:, 0].values, where="post", color=col, lw=0.9)
    cens = m[grp & (SM.pfs_ev == 0)]
    axI.plot(cens, k.survival_function_at_times(cens), "|", color=col, ms=1.8, mew=0.4)
    axI.text(sf.index.values[-1] + 0.4, sf.iloc[-1, 0], lab, color=col, fontsize=5.0, va="center", ha="left")
axI.set_xlim(0, 27); axI.set_xticks([0, 12, 24]); axI.set_ylim(0, 1.02); axI.set_yticks([0, 0.5, 1])
axI.set_yticklabels(["0", "0.5", "1"])
axI.set_xlabel("Months", labelpad=1); axI.set_ylabel("Progression-free", labelpad=1)
axI.text(0.98, 0.99, f"Immune axis, per SD\nHR {cx.HR:.2f} ({cx.lo:.2f}–{cx.hi:.2f})\n{fmt_p(cx.p)}",
         transform=axI.transAxes, fontsize=5.0, ha="right", va="top", color=INK, linespacing=1.15)
axI.set_title("PFS, NCT04249739", loc="left", fontsize=6.0, fontweight="bold", pad=8)

rec("i Immune axis PFS HR per SD", f"{cx.HR:.2f} ({cx.lo:.2f}–{cx.hi:.2f}), P {fmt_p(cx.p, '')}; {int(cx.n)} patients, "
    f"{int(cx.events)} events", "E22_PRJEB60680_cox.csv (PFS, Immune)")
rec("i KM groups", f"median split of Immune: {int(hi.sum())} high / {int((~hi).sum())} low",
    "descriptive, E22_scores_PRJEB60680.csv pfs_d/pfs_ev")
LET.append(("i", axI._left_title, 143.0))

# ------------------------------------------------------------------ j  distributions by response
SC3 = [("NetImmune", IMMUNE), ("IFN", IMMUNE), ("CAF_ECM", STROMA)]
TOPF, HF = 183.0, 26.0
rng = np.random.default_rng(7)
k = 0
for j, (key, name, tx, D, ycol) in enumerate(COH):
    for m, (sc, col) in enumerate(SC3):
        x0 = 13 + j * 85 + m * 27.5
        ax = S.ax_mm(fig, x0, TOPF, 21, HF)
        yv = D[ycol].astype(int).values
        for g in (0, 1):
            sel = yv == g
            v = D.loc[sel, sc].values
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.add_patch(mpatches.Rectangle((g - 0.27, q1), 0.54, q3 - q1, fc="none", ec=col, lw=0.6, zorder=2))
            ax.plot([g - 0.27, g + 0.27], [med, med], color=col, lw=1.0, zorder=3)
            for hot, mk in [(0, "o"), (1, "^")]:
                vv = D.loc[sel & (D.hot == hot), sc].values
                ax.plot(g + rng.uniform(-0.17, 0.17, len(vv)), vv, mk, ms=2.0 if mk == "o" else 2.4,
                        color=col, alpha=0.75, mew=0, zorder=4)
        auc = roc_auc_score(yv, D[sc].values)
        ax.text(0.0, 1.02, f"AUC {auc:.2f}", transform=ax.transAxes, fontsize=5.0, ha="left", va="bottom",
                color=SLATE)
        ax.set_xlim(-0.6, 1.6); ax.set_xticks([0, 1]); ax.set_xticklabels(["SD/PD", "CR/PR"], fontsize=5.0)
        ax.set_ylim(-3.2, 3.2); ax.set_yticks([-2, 0, 2]); ax.set_yticklabels(["−2", "0", "2"])
        S.hrule(ax, 0)
        if m == 0:
            ax.set_ylabel("Score (z)", labelpad=1)
        else:
            ax.set_yticklabels([])
        chip_title(fig, ax, sc, col, dy_mm=4.0)
        if m == 0:
            tF = fig.text(x0 * S.MM / S.W2, 1 - (TOPF - 9.5) * S.MM / (FIGH * S.MM),
                     "PRJEB25780" if j == 0 else "NCT04249739", fontsize=6.0, fontweight="bold", va="bottom",
                     ha="left", color=INK)
            fig.add_artist(plt.Line2D([x0 * S.MM / S.W2, (x0 + 76) * S.MM / S.W2],
                                      [1 - (TOPF - 9.0) * S.MM / (FIGH * S.MM)] * 2, color=LIGHT, lw=0.5))
        if j == 0 and m == 0:
            LET.append(("j", tF, 1.0))
# marker-shape key
kx = S.ax_mm(fig, 13, TOPF + HF + 6.5, 80, 3); S.clean(kx)
kx.set_xlim(0, 1); kx.set_ylim(0, 1)
kx.plot(0.0, 0.5, "o", ms=2.2, color=GREY, mew=0, clip_on=False); kx.text(0.025, 0.5, "MSS, EBV−", fontsize=5.0, va="center",
                                                              color=SLATE)
kx.plot(0.2, 0.5, "^", ms=2.6, color=GREY, mew=0, clip_on=False); kx.text(0.225, 0.5, "MSI-H or EBV+", fontsize=5.0, va="center",
                                                             color=SLATE)


place_letters(fig)
S.save(fig, "Fig4")
with open(S.OUT / "Fig4_numbers.md", "w") as fh:
    fh.write("# Fig. 4 - numbers printed on the figure and their sources\n\n| item | value | source |\n|---|---|---|\n")
    for kk, (v, s) in NUM.items():
        fh.write(f"| {kk} | {v} | {s} |\n")
print("ok")
