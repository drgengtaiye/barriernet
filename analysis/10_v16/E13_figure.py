"""Fig_E13 - single-cell tests of the barrier idea in public IMC (PDAC, METABRIC breast).
(a,b) PDAC ROI crops with the tumour region, fibroblasts, the contiguous rind layer and CD8 T cells.
(c) Positive control: continuity coefficient when random walkers are blocked by the rind layer.
(d) Observed coefficients vs a simulated complete barrier, raw and decomposed into composition-
    expected and excess (arrangement) parts.
(e) METABRIC partial-residual plot for raw continuity.
(f) Macrophage and comparator densities next to B/T aggregates vs > 150 um away, per patient.
(g) METABRIC disease-specific survival, HR per SD adjusted for age, NPI and ER."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, GREY, INK, NULL, LIGHT, MM, W2, fmt_p
import statsmodels.api as sm

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
MO = pd.read_csv(f"{R}/E13_models.csv")
PO = pd.read_csv(f"{R}/E13_poscontrol.csv")
MA = pd.read_csv(f"{R}/E13_macro_aggregates.csv")
MAP = pd.read_csv(f"{R}/E13_macro_aggregates_patients.csv")
SV = pd.read_csv(f"{R}/E13_survival.csv")
PM = pd.read_csv(f"{R}/E13_patient_metrics.csv")
ROI = pd.read_csv(f"{R}/E13_pdac_roi_model_data.csv").set_index("image")
CELLS = pd.read_csv(f"{R}/E13_pdac_cells_annotated.csv.gz")
TUM = "#E39A80"; FIBL = "#CFE9DD"; OTH = "#E4E8EC"
EX = [("P1-R1", "high continuity"), ("P1-R3", "low continuity")]


def roi_map(ax, key, lab):
    g = CELLS[CELLS.image == key]
    tumour, fib, cd8 = g.cls == 0, g.cls == 1, g.cls == 4      # class codes of E13_singlecell_barrier
    other = ~(tumour | fib | cd8)
    ax.scatter(g.x[other], g.y[other], s=0.35, c=OTH, lw=0, rasterized=True)
    ax.scatter(g.x[fib & ~g.goodrind], g.y[fib & ~g.goodrind], s=0.5, c=FIBL, lw=0, rasterized=True)
    ax.scatter(g.x[tumour], g.y[tumour], s=0.9, c=TUM, lw=0, rasterized=True)
    ax.scatter(g.x[fib & g.goodrind], g.y[fib & g.goodrind], s=1.6, c=STROMA, lw=0, rasterized=True)
    ax.scatter(g.x[cd8], g.y[cd8], s=2.2, c=IMMUNE, lw=0, rasterized=True)
    x0, x1, y0, y1 = g.x.min(), g.x.max(), g.y.min(), g.y.max()
    ax.set_xlim(x0 - 5, x1 + 5); ax.set_ylim(y1 + 5, y0 - 5); ax.set_aspect("equal"); ax.set_anchor("N")
    bs.fs.clean(ax)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(LIGHT); sp.set_linewidth(0.5)
    ax.plot([x1 - 215, x1 - 15], [y1 + 30] * 2, color=INK, lw=1.0, clip_on=False)
    ax.text(x1 - 115, y1 + 40, "200 µm", ha="center", va="top", fontsize=4.8)
    r = ROI.loc[key]
    bs.title(ax, f"PDAC {key}: {lab}", pad=2)
    ax.text(0.0, -0.035, f"continuity {r.cont:.2f} (null {r.cont_null:.2f})\n"
            f"CD8 in/out density, log₂ {r.lir:+.2f}", transform=ax.transAxes,
            va="top", ha="left", fontsize=5.0, color=bs.fs.SLATE)


def map_legend(fig, x, y):
    from matplotlib.lines import Line2D
    items = [(TUM, "tumour/ductal cell"), (FIBL, "fibroblast"),
             (STROMA, "fibroblast in contiguous boundary layer (run ≥ 5)"), (IMMUNE, "CD8 T cell"),
             (OTH, "other")]
    h = [Line2D([], [], ls="", marker="o", ms=3, mfc=c, mec="none") for c, _ in items]
    fig.legend(h, [l for _, l in items], loc="upper left", bbox_to_anchor=(x, y), ncol=3,
               fontsize=5.0, columnspacing=1.0, handletextpad=0.2, frameon=False)


def rowlabels(ax, yt, yl, heads):
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=5.3); ax.tick_params(axis="y", length=0)
    for t, l in zip(ax.get_yticklabels(), yl):
        if l in heads: t.set_fontweight("bold"); t.set_fontsize(5.4)


def posctrl(ax):
    yt, yl, heads, y = [], [], [], 0
    for ds, hd in (("PDAC", "PDAC, 33 ROIs"), ("METABRIC", "METABRIC, 361 patients")):
        yt.append(-y); yl.append(hd); heads.append(hd); y += 0.9
        q = PO[(PO.dataset == ds) & (PO.spec == "continuity") & (PO.term == "cont")]
        for p in (0.0, 0.5, 1.0):
            z = q[q.p_block == p][["coef", "lo", "hi"]].mean()
            col = STROMA if p == 1 else GREY
            ax.plot([z.lo, z.hi], [-y] * 2, color=col, lw=0.8)
            ax.plot(z.coef, -y, "o", ms=2.8, mfc="white", mec=col, mew=0.8)
            yt.append(-y); yl.append(f"simulated, block p = {p:g}"); y += 1
        o = MO[(MO.dataset == ds) & (MO.spec == "continuity") & (MO.term == "cont") &
               MO.level.str.startswith("ROI" if ds == "PDAC" else "patient")].iloc[0]
        ax.plot([o.lo, o.hi], [-y] * 2, color=INK, lw=0.9)
        ax.plot(o.coef, -y, "D", ms=2.8, color=INK)
        yt.append(-y); yl.append("observed"); y += 1.3
    rowlabels(ax, yt, yl, heads)
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(-0.8, 0.3); ax.set_ylim(-y + 0.9, 0.5)
    ax.set_xlabel("continuity coefficient\n(log₂ infiltration per SD)")
    bs.title(ax, "Positive control", pad=2)


def decomp(ax):
    terms = [("continuity", "cont", "raw continuity"),
             ("decomposed", "cont_null", "expected from composition"),
             ("decomposed", "cont_exc", "excess over null (arrangement)")]
    yt, yl, heads, y = [], [], [], 0
    tr = ax.get_yaxis_transform()
    for ds, hd in (("PDAC", "PDAC (ROIs, patient random intercept)"),
                   ("METABRIC", "METABRIC (patients, OLS)")):
        yt.append(-y); yl.append(hd); heads.append(hd); y += 0.9
        lev = "ROI" if ds == "PDAC" else "patient"
        for sp, t, lab in terms:
            o = MO[(MO.dataset == ds) & (MO.spec == sp) & (MO.term == t) & MO.level.str.startswith(lev)].iloc[0]
            z = PO[(PO.dataset == ds) & (PO.spec == sp) & (PO.term == t) & (PO.p_block == 1.0)][["coef", "lo", "hi"]].mean()
            ax.plot([z.lo, z.hi], [-y + 0.17] * 2, color=STROMA, lw=0.7)
            ax.plot(z.coef, -y + 0.17, "o", ms=2.7, mfc="white", mec=STROMA, mew=0.8)
            ax.plot([o.lo, o.hi], [-y - 0.17] * 2, color=INK, lw=0.9)
            ax.plot(o.coef, -y - 0.17, "D", ms=2.7, color=INK)
            ax.text(1.16, -y - 0.17, fmt_p(o.p, lead=""), transform=tr, fontsize=4.9, va="center",
                    ha="right", color=INK, clip_on=False)
            yt.append(-y); yl.append(lab); y += 1
        y += 0.3
    rowlabels(ax, yt, yl, heads)
    ax.text(1.16, 0, "P (obs.)", transform=tr, fontsize=4.8, ha="right", va="center", color=GREY)
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(-0.95, 0.5); ax.set_ylim(-y + 0.9, 0.5)
    ax.set_xlabel("coefficient (log₂ infiltration per SD), adjusted for\nfibroblast and CD8 density per mm²")
    ax.plot([], [], "D", color=INK, ms=2.8, label="observed")
    ax.plot([], [], "o", mfc="white", mec=STROMA, ms=2.8, label="simulated complete barrier")
    ax.legend(loc="center right", fontsize=4.9, handletextpad=0.3, bbox_to_anchor=(1.0, (4.2 - y + 0.9) / (y - 0.4) * -1))
    bs.title(ax, "Composition vs arrangement", pad=2)


def partial(ax):
    d = PM[PM.dataset == "METABRIC"].copy()
    for c in ("cont", "logS", "logCD8"):
        d["z_" + c] = (d[c] - d[c].mean()) / d[c].std()
    X = sm.add_constant(d[["z_logS", "z_logCD8"]])
    ry = sm.OLS(d.lir, X).fit().resid; rx = sm.OLS(d.z_cont, X).fit().resid
    ax.scatter(rx, ry, s=3, c=GREY, alpha=0.55, lw=0, rasterized=True)
    b = np.polyfit(rx, ry, 1); xx = np.linspace(rx.min(), rx.max(), 50)
    ax.plot(xx, np.polyval(b, xx), color=STROMA, lw=1.1)
    o = MO[(MO.dataset == "METABRIC") & (MO.spec == "continuity") & (MO.term == "cont")].iloc[0]
    ax.text(0.98, 0.97, f"{o.coef:+.2f} per SD ({o.lo:+.2f} to {o.hi:+.2f})\n{fmt_p(o.p)}; n = {int(o.n)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.1)
    ax.set_xlabel("raw continuity (SD), residual on\nfibroblast and CD8 density")
    ax.set_ylabel("log₂ CD8 in/out density, residual")
    ax.set_ylim(-5.5, 5)
    bs.title(ax, "METABRIC breast, per patient", pad=2)


def macro(ax):
    groups = [("METABRIC", "METABRIC (86 patients)",
               [("M", "macrophages", MACRO), ("M163", "CD163-high macrophages", MACRO),
                ("M163_share", "CD163-high share of macrophages", MACRO), ("S", "fibroblasts", STROMA),
                ("scaf", "all non-immune cells", GREY), ("CD8", "CD8 T cells", IMMUNE)]),
              ("PDAC", "PDAC (7 patients)",
               [("M", "macrophages", MACRO), ("S", "fibroblasts", STROMA),
                ("scaf", "all non-immune cells", GREY), ("CD8", "CD8 T cells", IMMUNE)])]
    rng = np.random.default_rng(0)
    tr = ax.get_yaxis_transform()
    yt, yl, heads, y = [], [], [], 0
    for ds, hd, rows in groups:
        yt.append(-y); yl.append(hd); heads.append(hd); y += 0.9
        P = MAP[(MAP.dataset == ds) & (MAP.region == "all")]
        for c, lab, col in rows:
            v = P[f"lr_{c}"].dropna().values
            ax.scatter(np.clip(v, -4.8, 5.8), -y + rng.uniform(-0.22, 0.22, len(v)),
                       s=2.2 if ds == "METABRIC" else 5, c=col,
                       alpha=0.35 if ds == "METABRIC" else 0.8, lw=0, rasterized=True)
            z = MA[(MA.dataset == ds) & (MA.region == "all") & (MA.contrast == c)].iloc[0]
            ax.plot([z.lo, z.hi], [-y] * 2, color=INK, lw=1.0)
            ax.plot(z.median_log2, -y, "|", color=INK, ms=6, mew=1.2)
            ax.text(1.13, -y, f"{int(z.n_lower_near)}/{int(z.n_pat)}", transform=tr, fontsize=4.9,
                    ha="right", va="center", clip_on=False)
            ax.text(1.38, -y, fmt_p(z.p_wilcoxon, lead=""), transform=tr, fontsize=4.9, ha="right",
                    va="center", clip_on=False)
            yt.append(-y); yl.append(lab); y += 1
        y += 0.3
    rowlabels(ax, yt, yl, heads)
    ax.text(1.13, 0, "lower near", transform=tr, fontsize=4.7, ha="right", va="center", color=GREY)
    ax.text(1.38, 0, "P", transform=tr, fontsize=4.7, ha="right", va="center", color=GREY)
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(-5, 6); ax.set_xticks([-4, -2, 0, 2, 4]); ax.set_ylim(-y + 0.9, 0.5)
    ax.set_xlabel("log₂ density ≤ 50 µm vs > 150 µm\nfrom a B/T-cell aggregate")
    bs.title(ax, "Macrophages next to B/T-cell aggregates", pad=2)


def surv(ax):
    d = SV[SV.endpoint == "dss"]
    rows = [("cont (+age, NPI, ER)", "cont", "raw continuity", STROMA),
            ("cont_exc (+age, NPI, ER)", "cont_exc", "excess continuity", STROMA),
            ("logS (+age, NPI, ER)", "logS", "fibroblast density", STROMA),
            ("lir (+age, NPI, ER)", "lir", "CD8 infiltration ratio", IMMUNE),
            ("logCD8 (+age, NPI, ER)", "logCD8", "CD8 density", IMMUNE),
            ("joint: cont + logS + logCD8 (+age, NPI, ER)", "cont", "raw continuity, joint model", STROMA)]
    tr = ax.get_yaxis_transform()
    for i, (m, t, lab, col) in enumerate(rows):
        r = d[(d.model == m) & (d.term == t)].iloc[0]
        ax.plot([r.lo, r.hi], [-i] * 2, color=col, lw=0.8)
        ax.plot(r.HR, -i, "o", color=col, ms=3.0, mec="white", mew=0.4)
        ax.text(1.17, -i, f"{r.HR:.2f}", transform=tr, fontsize=4.9, ha="right", va="center", clip_on=False)
        ax.text(1.40, -i, fmt_p(r.p, lead=""), transform=tr, fontsize=4.9, ha="right", va="center", clip_on=False)
    rowlabels(ax, -np.arange(len(rows)), [r[2] for r in rows], [])
    ax.text(1.17, 0.9, "HR", transform=tr, fontsize=4.7, ha="right", color=GREY)
    ax.text(1.40, 0.9, "P", transform=tr, fontsize=4.7, ha="right", color=GREY)
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xscale("log"); ax.set_xlim(0.72, 1.7); ax.set_xticks([0.8, 1, 1.25, 1.6])
    ax.set_xticklabels(["0.8", "1", "1.25", "1.6"]); ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_ylim(-len(rows) + 0.4, 1.3)
    r0 = d.iloc[0]
    ax.set_xlabel(f"disease-specific HR per SD, + age, NPI, ER\n(n = {int(r0.n)}, {int(r0.events)} deaths)")
    bs.title(ax, "METABRIC survival", pad=2)


FW, FH = 180.0, 200.0
fig = plt.figure(figsize=(FW * MM, FH * MM))


def axmm(x, y, w, h):
    return fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH])


def letter(x, y, L):
    fig.text(x / FW, 1 - y / FH, L, fontsize=8.2, fontweight="bold", va="top", ha="left", color=INK)


axa = axmm(4, 8, 52, 52); roi_map(axa, *EX[0])
axb = axmm(60, 8, 52, 52); roi_map(axb, *EX[1])
map_legend(fig, 4 / FW, 1 - 71 / FH)
axc = axmm(144, 8, 34, 60); posctrl(axc)
axd = axmm(52, 92, 48, 50); decomp(axd)
axe = axmm(128, 92, 48, 50); partial(axe)
axf = axmm(40, 158, 44, 34); macro(axf)
axg = axmm(136, 158, 30, 34); surv(axg)
for x, y, L in [(0, 2, "a"), (57, 2, "b"), (114, 2, "c"), (0, 86, "d"), (114, 86, "e"),
                (0, 152, "f"), (106, 152, "g")]:
    letter(x, y, L)
bs.save(fig, "Fig_E13_singlecell_barrier")
print("saved")
