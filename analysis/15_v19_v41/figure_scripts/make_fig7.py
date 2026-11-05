"""v19 Fig. 7 - the Barrier axis is compositional, not a demonstrated physical barrier (10 panels).

a-c  Visium, representative primary tumor (the primary whose SPP1_arm Ring1-Far contrast is closest to
     the median of the 10 published sections): rings around TLS peaks called on the TLS_B marker module,
     SPP1+ macrophage (SPP1_arm) and CAF_ECM module scores (analysis/E27_visium_marker_tls.py).
d    Profile of each score from TLS peak to far tissue (zone mean minus far mean) in the 10 published
     sections and the independent atlas GSE270678 (donor means); counts and two-sided Wilcoxon P from
     E27_rings_summary.csv. T cells are the positive control for the ring centers.
e    Section-level BarrierIndex of the published sections (E27_rings_sections.csv) with EBV, MSI, Lauren
     and pM from results/ring_barrier_state_v2/BarrierIndex_by_slice_with_clinical.csv.
f    Least-cost transport: detour z against the rigid-translation null with Tangram fractions
     (revision2/results/T3_transport.csv, z_shift) and with independent marker scores (T6_marker_based.csv, z).
g    Imaging mass cytometry, the two pancreatic regions of patient P1 with the highest and lowest
     fibroblast-layer continuity (E13_pdac_cells_annotated.csv.gz, E13_pdac_roi_model_data.csv).
h    Positive control: continuity coefficient when simulated CD8 walkers are blocked with probability
     0, 0.5, 1 (E13_poscontrol.csv, mean over replicates) against the observed coefficient (E13_models.csv).
i    Continuity split into the part expected from fibroblast number and the excess (arrangement):
     observed coefficients (E13_models.csv) against a simulated complete barrier (E13_poscontrol.csv,
     p_block = 1, mean over replicates).
j    Cell densities within 50 um vs beyond 150 um of B/T-cell aggregates, per patient
     (E13_macro_aggregates_patients.csv) with median and 95% CI and Wilcoxon P (E13_macro_aggregates.csv).
Zone means and contrasts are read from E27 and checked against each other on load.
The micro sign is written as U+03BC: Helvetica's U+00B5 glyph is slanted.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.stats import wilcoxon
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, MACRO, MM, fmt_p

ROOT = S.ROOT
DATA = S.OUT / "data"
RES = S.RES
VIS = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2")
SENS = ROOT / "code/analysis/06_spatial"
PITCH_UM = 100.0                               # Visium center-to-center spot distance
STATES = [("state_SPP1_TAM", "SPP1⁺ TAM", MACRO), ("state_CAF_ECM", "CAF_ECM", STROMA),
          ("state_IL10_STAT3", "IL-10/STAT3", SLATE)]
ZONE_COL = {"Hotspot": "#0A3D6E", "Ring1": IMMUNE, "Ring2": "#9CC3E6", "Far": "#E6EAEE"}
ZONES = ["Hotspot", "Ring1", "Ring2", "Far"]
minus = lambda s: s.replace("-", "−")


# ------------------------------------------------------------------ inputs (analysis/E27_visium_marker_tls.py)
# TLS peaks on the TLS_B marker module (the published Tangram TLS_niche_score did not track lymphoid
# markers; E26b), in the published sections (GSE251950) and the independent atlas (GSE270678, donor means).
SEC = pd.read_csv(RES / "E27_rings_sections.csv")
ZM_ = pd.read_csv(RES / "E27_zone_means.csv")
TS = pd.read_csv(RES / "E27_rings_summary.csv")
SP = pd.read_parquet(RES / "E27_spots_GSE251950.parquet")
STATES = [("SPP1_arm", "SPP1 subset", MACRO), ("CAF_ECM", "CAF_ECM", STROMA), ("T_CELL", "T cells (control)", IMMUNE)]
met = (SEC[(SEC.atlas == "GSE251950") & (SEC.peaks == "knn9") & (SEC.far == "far3")]
       .rename(columns={"section": "slice"}).reset_index(drop=True))
ind = SEC[SEC.atlas == "GSE270678"]
N_DON = ind.donor.nunique()
v_ = met.set_index("slice").SPP1_arm_Ring1_minus_Far
rep = (v_[met[met.is_metastasis == 0].slice] - v_.median()).abs().idxmin()   # primary closest to the median
REP_GC = dict(zip(met.slice, met.GC_raw))[rep]
spots = SP[SP.slice == rep].reset_index(drop=True)
zm = ZM_.copy()
far_ = zm[zm.zone == "Far"].set_index(["atlas", "unit", "score"])["mean"]
zm["rel"] = zm["mean"].values - far_.reindex(pd.MultiIndex.from_frame(zm[["atlas", "unit", "score"]])).values
summ = TS[(TS.peaks == "knn9") & (TS.far == "far3") & (TS.contrast == "Ring1-Far")]
MO = pd.read_csv(RES / "E13_models.csv")
PO = pd.read_csv(RES / "E13_poscontrol.csv")
MA = pd.read_csv(RES / "E13_macro_aggregates.csv")
MAP = pd.read_csv(RES / "E13_macro_aggregates_patients.csv")
ROI = pd.read_csv(RES / "E13_pdac_roi_model_data.csv").set_index("image")
CELLS = pd.read_csv(RES / "E13_pdac_cells_annotated.csv.gz")
CLIN = pd.read_csv(VIS / "BarrierIndex_by_slice_with_clinical.csv").drop(columns=["BarrierIndex", "z_ecm", "z_tam", "z_bar"])
R2 = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results")
T3 = pd.read_csv(R2 / "T3_transport.csv"); T6 = pd.read_csv(R2 / "T6_marker_based.csv")
NUM = []                                         # (panel, quantity, value, source) for Fig7_numbers.md

# consistency: zone means reproduce the E27 Ring1 - Far contrasts and the printed P values
for st, _, _ in STATES:
    for atl, ref in [("GSE251950", met.set_index("GC_raw")[f"{st}_Ring1_minus_Far"]),
                     ("GSE270678", ind.groupby("donor")[f"{st}_Ring1_minus_Far"].mean())]:
        r1 = zm[(zm.atlas == atl) & (zm.score == st) & (zm.zone == "Ring1")].set_index("unit").rel
        assert np.allclose(r1.sort_index(), ref.sort_index(), atol=1e-9), (st, atl)
        k_ = summ[(summ.atlas == atl) & (summ.score == st)].iloc[0]
        assert abs(wilcoxon(ref.values).pvalue - k_.P) < 1e-12, (st, atl)

# ------------------------------------------------------------------ figure
FW, FH = 180.0, 230.0
fig = plt.figure(figsize=(FW * MM, FH * MM))
ax_mm = lambda x, y, w, h: S.ax_mm(fig, x, y, w, h)
HEAD = {}                                        # panel letter -> artist carrying its bold title


def scalebar(ax, x1, y0, length, label, dy=0.0):
    ax.plot([x1 - length, x1], [y0, y0], color=INK, lw=1.0, solid_capstyle="butt", clip_on=False)
    ax.text(x1 - length / 2, y0 - dy, label, ha="center", va="top", fontsize=5.0, color=INK, clip_on=False)


def btitle(ax, s, pad=3, sub=None):
    """Bold panel title (optionally a muted subtitle directly above the axes)."""
    t = ax.set_title(s, loc="left", fontsize=6.5, fontweight="bold", color=INK, pad=pad)
    return t                                     # subtitles are carried by the legend


def rowlabels(ax, yt, yl, heads):
    ax.set_yticks(yt); ax.set_yticklabels(yl); ax.tick_params(axis="y", length=0)
    for tl in ax.get_yticklabels():
        if tl.get_text() in heads:
            tl.set_fontweight("bold")
    ax.spines["left"].set_visible(False)


# ================================================================ row 1: a-c Visium maps
MAPW, MAPY, MAPH = 50.0, 9.0, 44.0
pad = 0.035 * np.ptp(spots.x)
xlim = (spots.x.min() - pad, spots.x.max() + pad)
ylim = (spots.y.min() - pad - 0.10 * np.ptp(spots.y), spots.y.max() + pad)     # room for scale bar
from scipy.spatial import cKDTree
pitch = np.median(cKDTree(spots[["x", "y"]].values).query(spots[["x", "y"]].values, k=2)[0][:, 1])
px_per_um = pitch / PITCH_UM


def visium_axes(x):
    ax = ax_mm(x, MAPY, MAPW, MAPH)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.set_aspect("equal", adjustable="box"); ax.set_anchor("N")
    S.clean(ax)
    return ax


def spot_area(ax, fill=0.93):
    """Marker area (pt^2) so neighbouring spots almost touch at this panel's scale."""
    bb = ax.get_position(); W, H = fig.get_size_inches()
    units_per_in = max(np.diff(xlim)[0] / (bb.width * W), np.diff(ylim)[0] / (bb.height * H))
    return (pitch / units_per_in * 72 * fill) ** 2


seq = lambda c: LinearSegmentedColormap.from_list("s", ["#F3F5F7", c])
xs = [6.0, 66.0, 126.0]
axa = visium_axes(xs[0]); sa = spot_area(axa)
for z in ["Far", "Ring2", "Ring1", "Hotspot"]:
    k = spots.ring == z
    axa.scatter(spots.x[k], spots.y[k], s=sa, c=ZONE_COL[z], lw=0, rasterized=True)
cnt = spots.ring.value_counts()
NUM += [("a", f"{REP_GC} spots {z}", int(cnt[z]), "E27_spots_GSE251950.parquet (TLS_B peaks, knn9, far3)") for z in ZONES]
HEAD["a"] = btitle(axa, "Rings around lymphoid aggregates", pad=9, sub=f"{REP_GC}, primary tumor (Visium)")

maps = []
for L, ax_x, (st, name, col) in zip("bc", xs[1:], STATES[:2]):
    ax = visium_axes(ax_x)
    v = spots[st].astype(float)
    lo, hi = np.nanpercentile(v, [2, 98])
    sc = ax.scatter(spots.x, spots.y, s=sa, c=v, cmap=seq(col), vmin=lo, vmax=hi, lw=0, rasterized=True)
    k = spots.ring == "Hotspot"
    ax.scatter(spots.x[k], spots.y[k], s=sa * 0.55, facecolors="none", edgecolors=INK, linewidths=0.3)
    HEAD[L] = btitle(ax, name, pad=9, sub="module score; outlined, TLS peaks")
    maps.append((ax, sc, lo, hi, name))
    NUM += [(L, f"{name} colour range (2nd-98th pct)", f"{lo:.2f} to {hi:.2f}", "E27_spots_GSE251950.parquet")]

for ax in [axa] + [m[0] for m in maps]:
    xr = xlim[1] - 0.03 * np.ptp(spots.x)
    scalebar(ax, xr, ylim[0] + 0.05 * np.ptp(spots.y), 1000 * px_per_um, "1 mm", dy=0.025 * np.ptp(spots.y))

LY = MAPY + MAPH + 1.8
lax = ax_mm(xs[0], LY, MAPW, 3.0); S.clean(lax); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
for i, (z, lab) in enumerate(zip(ZONES, ["TLS peak", "ring 1", "ring 2", "far (≥ 3)"])):
    x0 = 0.02 + i * 0.245
    lax.scatter([x0], [0.5], s=9, c=ZONE_COL[z], marker="o", lw=0)
    lax.text(x0 + 0.035, 0.5, lab, va="center", fontsize=5.0, color=SLATE)
for ax, sc, lo, hi, name in maps:
    x0 = ax.get_position().x0 * FW
    cax = ax_mm(x0 + 2, LY + 0.6, 20, 1.6)
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.3)
    cb.set_ticks([lo, hi]); cb.set_ticklabels([minus(f"{lo:.2f}"), minus(f"{hi:.2f}")])
    cb.ax.tick_params(labelsize=5.0, length=1.2, width=0.3, pad=0.8)
    fig.text((x0 + 24) / FW, 1 - (LY + 1.4) / FH, "module score", fontsize=5.0, color=SLATE, va="center")

# ================================================================ row 2: d profiles, e BarrierIndex
DY, DH, DW = 78.0, 25.0, 22.0
dxs = [14.0, 44.0, 74.0]
ismet = dict(zip(met.GC_raw, met.is_metastasis.astype(bool)))
dax = []
for x, (st, name, col) in zip(dxs, STATES):
    ax = ax_mm(x, DY, DW, DH); dax.append(ax)
    w = zm[(zm.atlas == "GSE251950") & (zm.score == st)].pivot(index="unit", columns="zone", values="rel")[ZONES]
    wi = zm[(zm.atlas == "GSE270678") & (zm.score == st)].pivot(index="unit", columns="zone", values="rel")[ZONES]
    for _, row in wi.iterrows():
        ax.plot(range(4), row.values, color=GREY, lw=0.35, alpha=0.55, zorder=1)
    for gc, row in w.iterrows():
        if ismet[gc]:
            ax.plot(range(4), row.values, color=col, lw=0.6, ls=(0, (2, 1.2)), alpha=0.9, zorder=2)
            ax.scatter(range(4), row.values, s=5, facecolors="white", edgecolors=col, linewidths=0.5, zorder=3)
        else:
            ax.plot(range(4), row.values, color=col, lw=0.5, alpha=0.55, zorder=2)
    med, medi = w.median(), wi.median()
    ax.plot(range(4), med.values, color=INK, lw=1.1, zorder=4, solid_capstyle="round")
    ax.scatter(range(4), med.values, s=7, color=INK, zorder=5, lw=0)
    ax.plot(range(4), medi.values, color=INK, lw=0.9, ls=(0, (2.4, 1.2)), zorder=4)
    S.hrule(ax, 0)
    allv = np.r_[w.values.ravel(), wi.values.ravel()]; lo_, hi_ = np.nanmin(allv), np.nanmax(allv)
    ax.set_xlim(-0.3, 3.3); ax.set_ylim(lo_ - 0.06 * (hi_ - lo_), hi_ + 0.06 * (hi_ - lo_))
    ax.set_xticks(range(4)); ax.set_xticklabels(["peak", "1", "2", "far"])
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: minus(f"{v:g}")))
    fig.add_artist(Line2D([(x + 0.3) / FW, (x + 2.6) / FW], [1 - (DY - 2.3) / FH] * 2, color=col, lw=2.2,
                          solid_capstyle="butt"))
    fig.text((x + 3.6) / FW, 1 - (DY - 2.3) / FH, name, fontsize=6.0, color=INK, va="center", ha="left")
    k1 = summ[(summ.atlas == "GSE251950") & (summ.score == st)].iloc[0]
    k2 = summ[(summ.atlas == "GSE270678") & (summ.score == st)].iloc[0]
    txt = f"sections {fmt_p(k1.P)}\ndonors {fmt_p(k2.P)}"
    fig.text(x / FW, 1 - (DY + DH + 3.6) / FH, txt, fontsize=5.0, color=INK, va="top", ha="left", linespacing=1.2)
    NUM += [("d", f"{name} sections ring 1 > far; Wilcoxon P", f"{int(k1.n_above)}/{int(k1.n)}; {fmt_p(k1.P, lead='')}",
             "E27_rings_summary.csv GSE251950 knn9/far3"),
            ("d", f"{name} donors ring 1 > far; Wilcoxon P", f"{int(k2.n_above)}/{int(k2.n)}; {fmt_p(k2.P, lead='')}",
             "E27_rings_summary.csv GSE270678 (donor means)")]
dax[0].set_ylabel("score minus far mean")
HEAD["d"] = fig.text(dxs[0] / FW, 1 - (DY - 13.6) / FH, "From TLS peak to far tissue, two gastric atlases",
                     fontsize=6.5, fontweight="bold", color=INK, va="baseline")
h = [Line2D([], [], color=GREY, lw=0.6, label="published sections (10)"),
     Line2D([], [], color=GREY, lw=0.6, ls=(0, (2, 1.2)), marker="o", ms=2.2, mfc="white", mew=0.5,
            label="peritoneal metastasis"),
     Line2D([], [], color=INK, lw=1.1, marker="o", ms=2.2, label="median, sections"),
     Line2D([], [], color=GREY, lw=0.35, alpha=0.8, label=f"independent atlas, donors ({N_DON})"),
     Line2D([], [], color=INK, lw=0.9, ls=(0, (2.4, 1.2)), label="median, donors")]
fig.legend(handles=h, loc="upper left", bbox_to_anchor=(dxs[0] / FW, 1 - (DY - 12.2) / FH), ncol=3,
           fontsize=5.0, frameon=False, handlelength=2.0, columnspacing=1.0, borderaxespad=0, labelspacing=0.3)

# e: BarrierIndex by section with clinical strip (published sections)
EXX, EYY, EWW, EHH = 120.0, 71.0, 57.0, 22.0
ax = ax_mm(EXX, EYY, EWW, EHH)
bi = met[["slice", "GC_raw", "is_metastasis", "BarrierIndex"]].merge(
    CLIN[["slice", "EBV", "MSI", "Type", "pM"]], on="slice").sort_values("BarrierIndex", ascending=False).reset_index(drop=True)
xb = np.arange(len(bi))
ax.bar(xb, bi.BarrierIndex, width=0.68, color=[STROMA if m == 0 else "white" for m in bi.is_metastasis],
       edgecolor=STROMA, linewidth=0.7)
ax.axhline(0, color=INK, lw=0.5)
ax.set_xlim(-0.6, len(bi) - 0.4)
ax.set_ylim(min(-1.0, bi.BarrierIndex.min() - 0.4), max(1.0, bi.BarrierIndex.max() + 0.6))
ax.set_xticks(xb); ax.set_xticklabels([g_.replace("-PM", "\nPM") for g_ in bi.GC_raw], fontsize=5.0, linespacing=0.95, va="top")
ax.tick_params(axis="x", length=0, pad=1.2)
ax.set_ylabel("BarrierIndex (sum of z)")
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: minus(f"{v:g}")))
HEAD["e"] = btitle(ax, "Section-level BarrierIndex", pad=9, sub="CAF_ECM + SPP1 subset, ring 1 − far; open, metastasis")
STRIP = [("EBV", lambda r: ("+", True) if r.EBV == "positive" else ("−", False)),
         ("MSI", lambda r: ("H", True) if r.MSI == "MSI" else ("S", False)),
         ("Lauren", lambda r: ("I", True) if r.Type == "Intestinal" else ("D", False)),
         ("pM", lambda r: (str(int(r.pM)), int(r.pM) == 1))]
sax = ax_mm(EXX, EYY + EHH + 5.2, EWW, 11.0)
sax.set_xlim(-0.6, len(bi) - 0.4); sax.set_ylim(len(STRIP) - 0.5, -0.5); S.clean(sax)
for k_, (name, f) in enumerate(STRIP):
    sax.text(-0.85, k_, name, ha="right", va="center", fontsize=5.0, color=INK)
    for i, r in bi.iterrows():
        txt, on = f(r)
        sax.add_patch(plt.Rectangle((i - 0.42, k_ - 0.42), 0.84, 0.84, color=SLATE if on else S.PAPER, lw=0))
        sax.text(i, k_, txt, ha="center", va="center", fontsize=5.0, color="white" if on else INK)

for _, r in bi.iterrows():
    NUM += [("e", f"BarrierIndex {r.GC_raw}", f"{r.BarrierIndex:.2f}", "E27_rings_sections.csv knn9/far3")]


# ================================================================ row 3: f transport, g IMC maps, h positive control
R3 = 126.0
# f: least-cost transport
gcmap = dict(zip(met.slice, met.GC_raw))
tj = T3.assign(GC=T3.section.map(gcmap)).merge(T6.assign(GC=T6.section.map(gcmap)), on="GC", suffixes=("_tg", "_mk"))
ax = ax_mm(14.0, R3, 22.0, 32.0)
for _, r in tj.iterrows():
    ax.plot([0, 1], [r.z_shift, r.z], color=GREY, lw=0.5, alpha=0.7, zorder=1)
pm = tj.GC.eq("GC6-PM").values
for x_, v_, c_ in [(0, tj.z_shift.values, STROMA), (1, tj.z.values, INK)]:
    ax.scatter(np.full((~pm).sum(), x_), v_[~pm], s=9, color=c_, zorder=3, lw=0)
    ax.scatter(np.full(pm.sum(), x_), v_[pm], s=9, facecolors="white", edgecolors=c_, linewidths=0.7, zorder=3)
for x_, v_ in [(0, tj.z_shift.median()), (1, tj.z.median())]:
    ax.plot([x_ - 0.2, x_ + 0.2], [v_] * 2, color=INK, lw=1.1, zorder=4)
S.hrule(ax, 0)
ax.set_xlim(-0.45, 1.45); ax.set_ylim(-6.2, 2.6)
ax.set_xticks([0, 1]); ax.set_xticklabels(["Tangram\nfractions", "marker\nscores"], linespacing=1.0)
ax.tick_params(axis="x", length=0)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: minus(f"{v:g}")))
ax.set_ylabel("detour z vs null (+ = obstruction)")
n_tg, n_mk = int((tj.z_shift < 0).sum()), int((tj.z > 0).sum())

HEAD["f"] = btitle(ax, "Least-cost transport", pad=9, sub="10 sections; bars, medians")
NUM += [("f", "Tangram z: sections < 0; median", f"{n_tg}/10; {tj.z_shift.median():.2f}", "revision2/results/T3_transport.csv z_shift"),
        ("f", "marker-score z: sections > 0; median", f"{n_mk}/10; {tj.z.median():.2f}", "revision2/results/T6_marker_based.csv z")]

# g: IMC maps (pancreas)
TUM, FIBL, OTH = "#EDB5A1", "#BFE0D0", "#E4E8EC"
EX = [("P1-R1", "High continuity"), ("P1-R3", "Low continuity")]
EW = 34.0; exs = [54.0, 91.0]
eax = []
for x, (key, lab) in zip(exs, EX):
    ax = ax_mm(x, R3, EW, EW); eax.append(ax)
    g = CELLS[CELLS.image == key]
    tumor, fib, cd8 = g.cls == 0, g.cls == 1, g.cls == 4
    other = ~(tumor | fib | cd8)
    ax.scatter(g.x[other], g.y[other], s=0.10, c=OTH, lw=0, rasterized=True)
    ax.scatter(g.x[fib & ~g.goodrind], g.y[fib & ~g.goodrind], s=0.22, c=FIBL, lw=0, rasterized=True)
    ax.scatter(g.x[tumor], g.y[tumor], s=0.30, c=TUM, lw=0, rasterized=True)
    ax.scatter(g.x[fib & g.goodrind], g.y[fib & g.goodrind], s=0.8, c=STROMA, lw=0, rasterized=True)
    ax.scatter(g.x[cd8], g.y[cd8], s=1.1, c=IMMUNE, lw=0, rasterized=True)
    ax.set_xlim(0, 1000); ax.set_ylim(1000, 0); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(LIGHT); sp.set_linewidth(0.5)
    r = ROI.loc[key]
    t = btitle(ax, lab, pad=9, sub=f"pancreatic cancer, region {key}")
    if key == "P1-R1":
        HEAD["g"] = t
    ax.text(0.0, -0.03, f"continuity {r.cont:.2f} (expected {r.cont_null:.2f})\n"
            f"CD8 in/out of tumor, log₂ {minus(f'{r.lir:+.2f}')}", transform=ax.transAxes,
            va="top", ha="left", fontsize=5.0, color=SLATE, linespacing=1.2)
    ax.plot([1000 - 200, 1000], [1040] * 2, color=INK, lw=1.0, clip_on=False, solid_capstyle="butt")
    ax.text(900, 1055, "200 μm", ha="center", va="top", fontsize=5.0, clip_on=False)
    NUM += [("g", f"{key} continuity / expected / log2 CD8 in-out", f"{r.cont:.2f} / {r.cont_null:.2f} / {r.lir:+.2f}",
             "E13_pdac_roi_model_data.csv cont, cont_null, lir")]

leg = ax_mm(exs[0], R3 + EW + 12.6, exs[1] + EW - exs[0], 5.0); S.clean(leg); leg.set_xlim(0, 1); leg.set_ylim(0, 1)
for (c, l), (cx, cy) in zip([(TUM, "tumor cell"), (FIBL, "other fibroblast"), (STROMA, "boundary-layer fibroblast"),
                             (OTH, "other cell"), (IMMUNE, "CD8 T cell")],
                            [(0, 0.8), (0, 0.2), (0.31, 0.8), (0.31, 0.2), (0.74, 0.8)]):
    leg.scatter([cx + 0.012], [cy], s=7, c=c, lw=0)
    leg.text(cx + 0.04, cy, l, fontsize=5.0, va="center", color=SLATE)

# h: positive control, dose-response of the continuity coefficient
ax = ax_mm(141.0, R3, 36.0, 32.0)
DSN = [("METABRIC", "patient", "breast", INK, "o"), ("PDAC", "ROI", "pancreas", GREY, "s")]
for j_, (ds, lev, nm, col, mk) in enumerate(DSN):
    off = (j_ - 0.5) * 0.07
    q = PO[(PO.dataset == ds) & (PO.spec == "continuity") & (PO.term == "cont")].groupby("p_block")[["coef", "lo", "hi"]].mean()
    xx = q.index.values + off
    ax.plot(xx, q.coef, color=col, lw=0.8, zorder=2)
    ax.vlines(xx, q.lo, q.hi, color=col, lw=0.7, zorder=2)
    ax.scatter(xx, q.coef, s=10, marker=mk, facecolors="white", edgecolors=col, linewidths=0.8, zorder=3)
    o = MO[(MO.dataset == ds) & (MO.spec == "continuity") & (MO.term == "cont") & MO.level.str.startswith(lev)].iloc[0]
    ax.vlines(1.42 + off, o.lo, o.hi, color=col, lw=0.9)
    ax.scatter(1.42 + off, o.coef, s=12, marker="D", color=col, zorder=3, lw=0)
    for pb, rr in q.iterrows():
        NUM += [("h", f"{ds} simulated block p={pb:g}: coef (CI)", f"{rr.coef:.2f} ({rr.lo:.2f} to {rr.hi:.2f})",
                 "E13_poscontrol.csv spec=continuity term=cont, mean over reps")]
    NUM += [("h", f"{ds} observed coef (CI)", f"{o.coef:.2f} ({o.lo:.2f} to {o.hi:.2f})", "E13_models.csv continuity cont")]
ax.axvline(1.24, color=LIGHT, lw=0.6)
S.hrule(ax, 0)
ax.set_xlim(-0.15, 1.6); ax.set_ylim(-0.8, 0.25)
ax.set_xticks([0, 0.5, 1, 1.42]); ax.set_xticklabels(["0", "0.5", "1", "obs."])
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: minus(f"{v:g}")))
ax.set_xlabel("simulated block probability")
ax.set_ylabel("continuity coefficient")
HEAD["h"] = btitle(ax, "Positive control", pad=9, sub="simulated CD8 walkers (open) vs observed")
ax.legend(handles=[Line2D([], [], color=c_, lw=0.8, marker=m_, ms=2.8, mfc="white", mec=c_, mew=0.8, label=n_)
                   for _, _, n_, c_, m_ in DSN], loc="lower left", fontsize=5.0, frameon=False, handlelength=1.6,
          borderaxespad=0.3)

# ================================================================ row 4: i number vs arrangement, j aggregates
R4 = 190.0
ax = ax_mm(38.0, R4, 44.0, 33.0)
TERMS = [("continuity", "cont", "raw continuity"),
         ("decomposed", "cont_null", "expected from fibroblast number"),
         ("decomposed", "cont_exc", "excess (arrangement)")]
DS = [("METABRIC", "patient", "Breast, 361 patients"), ("PDAC", "ROI", "Pancreas, 33 regions (9 patients)")]
tr = ax.get_yaxis_transform()
yt, yl, heads, y = [], [], [], 0.0
for ds, lev, hd in DS:
    yt.append(-y); yl.append(hd); heads.append(hd); y += 1.0
    for sp, t, lab in TERMS:
        o = MO[(MO.dataset == ds) & (MO.spec == sp) & (MO.term == t) & MO.level.str.startswith(lev)].iloc[0]
        z = PO[(PO.dataset == ds) & (PO.spec == sp) & (PO.term == t) & (PO.p_block == 1.0)][["coef", "lo", "hi"]].mean()
        ax.plot([z.lo, z.hi], [-y + 0.2] * 2, color=GREY, lw=0.7)
        ax.scatter(z.coef, -y + 0.2, s=9, facecolors="white", edgecolors=GREY, linewidths=0.7, zorder=3)
        ax.plot([o.lo, o.hi], [-y - 0.15] * 2, color=INK, lw=1.0)
        ax.scatter(o.coef, -y - 0.15, s=12, marker="D", color=INK, zorder=4, lw=0)
        ax.text(1.03, -y - 0.15, fmt_p(o.p, lead=""), transform=tr, fontsize=5.0, va="center", ha="left",
                color=INK, clip_on=False)
        yt.append(-y); yl.append(lab); y += 1.0
        NUM += [("i", f"{ds} {t} observed", f"{o.coef:.2f} ({o.lo:.2f} to {o.hi:.2f}); P {fmt_p(o.p, lead='')}",
                 f"E13_models.csv level {lev}*"),
                ("i", f"{ds} {t} simulated complete barrier", f"{z.coef:.2f} ({z.lo:.2f} to {z.hi:.2f})",
                 "E13_poscontrol.csv p_block=1, mean over reps")]
    y += 0.35
rowlabels(ax, yt, yl, heads)
ax.text(1.03, 0.25, "P", transform=tr, fontsize=5.0, ha="left", va="center", color=SLATE)
S.rule(ax, 0)
ax.set_xlim(-0.85, 0.25); ax.set_ylim(-y + 0.55, 0.55)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: minus(f"{v:g}")))
ax.set_xlabel("log₂ CD8 in/out ratio per SD of continuity")
HEAD["i"] = btitle(ax, "Fibroblast number and arrangement", pad=11)
hh = [Line2D([], [], color=INK, lw=1.0, marker="D", ms=2.6, label="observed"),
      Line2D([], [], color=GREY, lw=0.7, marker="o", ms=2.8, mfc="white", mew=0.7, label="simulated complete barrier")]
ax.legend(handles=hh, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=5.0, frameon=False,
          handlelength=1.8, columnspacing=1.0, borderaxespad=0.15)

# j: densities near B/T aggregates
GX, GW = 127.0, 31.0
ax = ax_mm(GX, R4, GW, 33.0)
ROWS = [("M", "macrophages", MACRO), ("S", "fibroblasts", STROMA), ("CD8", "CD8 T cells", IMMUNE),
        ("scaf", "non-immune cells", GREY)]
rng = np.random.default_rng(0)
tr = ax.get_yaxis_transform()
C1 = (166.5 - GX) / GW; C2 = (168.5 - GX) / GW
yt, yl, heads, y = [], [], [], 0.0
for ds, hd in (("METABRIC", "Breast, 86 patients"), ("PDAC", "Pancreas, 7 patients")):
    yt.append(-y); yl.append(hd); heads.append(hd); y += 1.0
    P = MAP[(MAP.dataset == ds) & (MAP.region == "all")]
    for c, lab, col in ROWS:
        v = P[f"lr_{c}"].dropna().values
        ax.scatter(np.clip(v, -4.3, 5.3), -y + rng.uniform(-0.2, 0.2, len(v)), s=1.8 if ds == "METABRIC" else 4.0,
                   c=col, alpha=0.35 if ds == "METABRIC" else 0.75, lw=0, rasterized=True)
        z = MA[(MA.dataset == ds) & (MA.region == "all") & (MA.contrast == c)].iloc[0]
        ax.plot([z.lo, z.hi], [-y] * 2, color=INK, lw=1.0, zorder=4)
        ax.plot([z.median_log2] * 2, [-y - 0.28, -y + 0.28], color=INK, lw=1.2, zorder=4)
        ax.text(C1, -y, f"{int(z.n_lower_near)}/{int(z.n_pat)}", transform=tr, fontsize=5.0, ha="right",
                va="center", clip_on=False)
        ax.text(C2, -y, fmt_p(z.p_wilcoxon, lead=""), transform=tr, fontsize=5.0, ha="left", va="center",
                clip_on=False)
        yt.append(-y); yl.append(lab); y += 1.0
        NUM += [("j", f"{ds} {c} median log2 near/far (95% CI); lower near n/N; P",
                 f"{z.median_log2:.2f} ({z.lo:.2f} to {z.hi:.2f}); {int(z.n_lower_near)}/{int(z.n_pat)}; "
                 f"{fmt_p(z.p_wilcoxon, lead='')}", "E13_macro_aggregates.csv region=all")]
    y += 0.35
rowlabels(ax, yt, yl, heads)
ax.text(C1, 0.25, "lower\nnear", transform=tr, fontsize=5.0, ha="right", va="center", color=SLATE, linespacing=1.0)
ax.text(C2, 0.25, "P", transform=tr, fontsize=5.0, ha="left", va="center", color=SLATE)
S.rule(ax, 0)
ax.set_xlim(-4.5, 5.5); ax.set_xticks([-4, -2, 0, 2, 4]); ax.set_ylim(-y + 0.55, 0.55)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: minus(f"{v:g}")))
ax.set_xlabel("log₂ density, ≤ 50 μm vs > 150 μm")
HEAD["j"] = btitle(ax, "Near B/T-cell aggregates", pad=11, sub="per patient; median and 95% CI")

# ================================================================ panel letters
CAP = 0.717                                       # Helvetica cap height, em
fig.canvas.draw()
rend = fig.canvas.get_renderer()
inv = fig.transFigure.inverted()


def letter_on(x_mm, art, L, size=8.0):
    """Letter whose cap height is level with the cap height of the panel's bold title."""
    bb = art.get_window_extent(rend)
    _, _, desc = rend.get_text_width_height_descent(art.get_text().split("\n")[0], art.get_fontproperties(), ismath=False)
    y = bb.y0 + desc + (art.get_fontsize() - size) * CAP * fig.dpi / 72
    fig.text(x_mm / FW, inv.transform((0, y))[1], L, fontsize=size, fontweight="bold", va="baseline", ha="left",
             color=INK)


for L, x in [("a", 1.0), ("b", 61.0), ("c", 121.0), ("d", 1.0), ("e", 106.0), ("f", 1.0), ("g", 46.0),
             ("h", 127.0), ("i", 1.0), ("j", 98.0)]:
    letter_on(x, HEAD[L], L)

S.save(fig, "Fig7")

with open(S.OUT / "Fig7_numbers.md", "w") as f:
    f.write("# Fig. 7 - numbers printed on the figure\n\n| Panel | Quantity | Value | Source |\n|---|---|---|---|\n")
    for p_, q_, v_, s_ in NUM:
        f.write(f"| {p_} | {q_} | {v_} | {s_} |\n")
print("saved", S.OUT / "Fig7.png")
