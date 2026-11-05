"""v19 supplementary spatial figures in the cns house style.

  FigS_visium       a,b barrier state in the GC6 primary tumor and its peritoneal metastasis; c BarrierIndex
                    components (GC8 maps, ring contrasts and BarrierIndex by section are Fig. 7a-e).
  FigS_visium_sens  a ring zones of the nine other sections; b ring contrasts under two alternative
                    definitions (lattice peaks; Far >= 4 hops excluding tissue edges).
  FigS_barrierindex BarrierIndex against age, pT and pN in the nine primary tumors (Spearman).
  FigS_moran        a Moran's I per state and section; b feature comparison with spatial tools.
  FigS_transport    a marker-based detour ratios; b 27 parameter settings (Tangram vs marker z is Fig. 7f).
  FigS_imc          METABRIC imaging mass cytometry: a per-patient continuity vs CD8 infiltration (partial
                    residuals); b disease-specific survival (the other imaging panels are Fig. 7g-j).

Sources (no new tests; the Spearman correlations of S_barrierindex and the one-sided Wilcoxon P of
S_transport a are the ones the previous supplementary script figures/supp/make_supp_spatial.py
computed and printed on the same panels; the transport P is no longer shown):
  ring zones          analysis/results/E27_spots_GSE251950.parquet (TLS_B marker peaks, knn9, far >= 3)
  ring contrasts      analysis/results/E27_rings_sections.csv, E27_rings_summary.csv
  clinical            GC_TME_Reprogrammer/results/ring_barrier_state_v2/BarrierIndex_by_slice_with_clinical.csv
  Moran's I           GC_TME_Reprogrammer/validation/bench/spatial_benchmark_results.csv
  transport           revision2/results/T3_transport.csv, T4_sensitivity.csv, T6_marker_based.csv
  imaging             analysis/results/E13_*.csv
The micro sign is written as U+03BC: Helvetica's U+00B5 glyph is slanted.
Usage: python make_supp_spatial.py [key ...]   (default: all six)
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, wilcoxon
import statsmodels.api as sm
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, PAPER, IMMUNE, STROMA, MACRO, MM
from matplotlib.transforms import ScaledTranslation


def fmt_p(p, lead="P = "):
    """House P format (cns_style.fmt_p) with trailing zeros kept to two significant figures (0.020, not 0.02)."""
    if p is not None and p == p and 0.001 <= p < 0.1:
        return f"{lead}{p:#.2g}"
    return S.fmt_p(p, lead=lead)

ROOT = S.ROOT
DATA = S.OUT / "data"
RES = S.RES
VIS = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2")
BENCH = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation/bench")
R2 = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results")
SENS = ROOT / "code/analysis/06_spatial"
PITCH_UM = 100.0
ZONE_COL = {"Hotspot": "#0A3D6E", "Ring1": IMMUNE, "Ring2": "#9CC3E6", "Far": "#E6EAEE"}
ZONES = ["Hotspot", "Ring1", "Ring2", "Far"]
ZLAB = ["TLS peak", "ring 1", "ring 2", "far (≥ 3)"]
STATES = [("SPP1_arm", "SPP1 subset", MACRO), ("CAF_ECM", "CAF_ECM", STROMA), ("state_IL10_STAT3", "IL-10/STAT3", SLATE)]
CAP = 0.717                                       # Helvetica cap height, em
minus = lambda s: s.replace("-", "−")
fmtm = plt.FuncFormatter(lambda v, _: minus(f"{v:g}"))
seq = lambda c: LinearSegmentedColormap.from_list("s", ["#F3F5F7", c])

# ------------------------------------------------------------------ shared data
clin = pd.read_csv(VIS / "BarrierIndex_by_slice_with_clinical.csv").drop(columns=["BarrierIndex", "z_ecm", "z_tam", "z_bar"])
# rings around TLS peaks on the TLS_B marker module (analysis/E27_visium_marker_tls.py; the published
# Tangram TLS_niche_score did not track lymphoid markers)
E27 = pd.read_csv(RES / "E27_rings_sections.csv")
rs = E27[E27.atlas == "GSE251950"].rename(columns={"section": "slice"})
sm_ = pd.read_csv(RES / "E27_rings_summary.csv").query("atlas == 'GSE251950'")
met = rs.query("peaks == 'knn9' and far == 'far3'").reset_index(drop=True)
z0 = lambda v: (v - v.mean()) / v.std(ddof=0)
assert np.allclose(z0(met.CAF_ECM_Ring1_minus_Far) + z0(met.SPP1_arm_Ring1_minus_Far), met.BarrierIndex)
GC = dict(zip(met.slice, met.GC_raw))
ORDER = sorted(GC, key=lambda s: (GC[s].endswith("PM"), GC[s]))           # GC1..GC9, then GC6-PM
rings = pd.read_parquet(RES / "E27_spots_GSE251950.parquet")


# ------------------------------------------------------------------ helpers
class Fig:
    def __init__(self, w, h):
        self.W, self.H = w, h
        self.fig = plt.figure(figsize=(w * MM, h * MM))
        self.head = {}
        self.num = []

    def ax(self, x, y, w, h):
        return S.ax_mm(self.fig, x, y, w, h)

    def text(self, x, y, s, **kw):
        kw.setdefault("fontsize", 5.0); kw.setdefault("color", SLATE); kw.setdefault("va", "top")
        return self.fig.text(x / self.W, 1 - y / self.H, s, **kw)

    def letters(self, spec):
        """Letters (x mm, panel title artist) with cap height level with the bold title."""
        f = self.fig
        f.canvas.draw()
        rend = f.canvas.get_renderer(); inv = f.transFigure.inverted()
        for L, x in spec:
            art = self.head[L]
            bb = art.get_window_extent(rend)
            _, _, d = rend.get_text_width_height_descent(art.get_text().split("\n")[0], art.get_fontproperties(), ismath=False)
            y = bb.y0 + d + (art.get_fontsize() - 8.0) * CAP * f.dpi / 72
            f.text(x / self.W, inv.transform((0, y))[1], L, fontsize=8.0, fontweight="bold", va="baseline",
                   ha="left", color=INK)

    def save(self, key):
        out = HERE / f"FigS_{key}"
        self.fig.savefig(f"{out}.pdf"); self.fig.savefig(f"{out}.png", dpi=600)
        plt.close(self.fig)
        with open(HERE / f"FigS_{key}_numbers.md", "w") as fh:
            fh.write(f"# Supplementary Fig. ({key}) - numbers printed on the figure\n\n"
                     "| Panel | Quantity | Value | Source |\n|---|---|---|---|\n")
            for p, q, v, s in self.num:
                fh.write(f"| {p} | {q} | {v} | {s} |\n")
        print("saved", out)


SUBGAP = 1.6                                      # pt between the top of the axes and a subtitle


def btitle(ax, s, pad=3, sub=None, x=0.0):
    """Bold title starting at axes-x `x`; subtitles go to the legend."""
    sub = None
    t = ax.set_title(s, loc="left", fontsize=6.5, fontweight="bold", color=INK, pad=pad + (SUBGAP if sub else 0), x=x)
    if sub:
        ax.text(x, 1.0, sub, fontsize=5.0, color=SLATE, va="bottom", ha="left",
                transform=ax.transAxes + ScaledTranslation(0, SUBGAP / 72, ax.figure.dpi_scale_trans))
    return t


def rowlabels(ax, yt, yl, heads):
    ax.set_yticks(yt); ax.set_yticklabels(yl); ax.tick_params(axis="y", length=0)
    for tl in ax.get_yticklabels():
        if tl.get_text() in heads:
            tl.set_fontweight("bold")
    ax.spines["left"].set_visible(False)


def chip(F, x, y, col, label, size=6.0):
    """Family-colour chip followed by an ink label (x, y in mm; y is the text center)."""
    F.fig.add_artist(Line2D([(x) / F.W, (x + 2.3) / F.W], [1 - y / F.H] * 2, color=col, lw=2.2,
                            solid_capstyle="butt"))
    return F.text(x + 3.3, y, label, fontsize=size, color=INK, va="center", ha="left")


class Section:
    """Spot table of one Visium section with its own frame, spot pitch and marker size."""

    def __init__(self, df):
        self.d = df.reset_index(drop=True)
        P = self.d[["x", "y"]].values
        self.pitch = np.median(cKDTree(P).query(P, k=2)[0][:, 1])
        pad = 0.035 * max(np.ptp(self.d.x), np.ptp(self.d.y))
        self.xlim = (self.d.x.min() - pad, self.d.x.max() + pad)
        self.ylim = (self.d.y.min() - pad - 0.10 * np.ptp(self.d.y), self.d.y.max() + pad)

    def axes(self, F, x, y, w, h, anchor="N"):
        ax = F.ax(x, y, w, h)
        ax.set_xlim(*self.xlim); ax.set_ylim(*self.ylim); ax.set_aspect("equal", adjustable="box")
        ax.set_anchor(anchor); S.clean(ax)
        F.fig.canvas.draw()                       # resolve the aspect before sizing the spots
        bb = ax.get_position(); W, H = F.fig.get_size_inches()
        upi = max(np.diff(self.xlim)[0] / (bb.width * W), np.diff(self.ylim)[0] / (bb.height * H))
        self.s = (self.pitch / upi * 72 * 0.93) ** 2
        return ax

    def zones(self, ax):
        for z in ["Far", "Ring2", "Ring1", "Hotspot"]:
            k = self.d.ring == z
            ax.scatter(self.d.x[k], self.d.y[k], s=self.s, c=ZONE_COL[z], lw=0, rasterized=True)

    def values(self, ax, v, cmap, lo, hi, peaks=True):
        sc = ax.scatter(self.d.x, self.d.y, s=self.s, c=v, cmap=cmap, vmin=lo, vmax=hi, lw=0, rasterized=True)
        if peaks:
            k = self.d.ring == "Hotspot"
            ax.scatter(self.d.x[k], self.d.y[k], s=self.s * 0.55, facecolors="none", edgecolors=INK, linewidths=0.3)
        return sc

    def scalebar(self, ax):
        x1 = self.xlim[1] - 0.03 * np.ptp(self.d.x)
        y0 = self.ylim[0] + 0.05 * np.ptp(self.d.y)
        L = 1000 * self.pitch / PITCH_UM
        ax.plot([x1 - L, x1], [y0, y0], color=INK, lw=1.0, solid_capstyle="butt", clip_on=False)
        ax.text(x1 - L / 2, y0 - 0.025 * np.ptp(self.d.y), "1 mm", ha="center", va="top", fontsize=5.0,
                color=INK, clip_on=False)


def zone_legend(F, x, y, w):
    lax = F.ax(x, y, w, 3.0); S.clean(lax); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
    for i, (z, lab) in enumerate(zip(ZONES, ZLAB)):
        x0 = 0.02 + i * 0.245
        lax.scatter([x0], [0.5], s=9, c=ZONE_COL[z], marker="o", lw=0)
        lax.text(x0 + 0.035, 0.5, lab, va="center", fontsize=5.0, color=SLATE)


def colorbar(F, sc, x, y, lo, hi, label):
    cax = F.ax(x, y, 20, 1.6)
    cb = F.fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.outline.set_linewidth(0.3)
    cb.set_ticks([lo, hi]); cb.set_ticklabels([minus(f"{lo:.2f}"), minus(f"{hi:.2f}")])
    cb.ax.tick_params(labelsize=5.0, length=1.2, width=0.3, pad=0.8)
    F.text(x + 22, y + 0.8, label, va="center")


def bi_strip(F, ax, x, y, w, bi):
    """Clinical strip under a BarrierIndex bar chart, with its key."""
    STRIP = [("EBV", lambda r: ("+", True) if r.EBV == "positive" else ("−", False)),
             ("MSI", lambda r: ("H", True) if r.MSI == "MSI" else ("S", False)),
             ("Lauren", lambda r: ("I", True) if r.Type == "Intestinal" else ("D", False)),
             ("pM", lambda r: (str(int(r.pM)), int(r.pM) == 1))]
    sax = F.ax(x, y, w, 11.0)
    sax.set_xlim(-0.6, len(bi) - 0.4); sax.set_ylim(len(STRIP) - 0.5, -0.5); S.clean(sax)
    for k, (name, f) in enumerate(STRIP):
        sax.text(-0.85, k, name, ha="right", va="center", fontsize=5.0, color=INK)
        for i, r in bi.iterrows():
            txt, on = f(r)
            sax.add_patch(plt.Rectangle((i - 0.42, k - 0.42), 0.84, 0.84, color=SLATE if on else PAPER, lw=0))
            sax.text(i, k, txt, ha="center", va="center", fontsize=5.0, color="white" if on else INK)
    F.text(x - 10.0, y + 12.0, "EBV +/−; MSI: S stable, H high; Lauren: I intestinal, D diffuse;\n"
           "pM 0/1. Dark tile: +, H, I or 1.", linespacing=1.2)


# ================================================================== FigS_visium
def fig_visium():
    """One patient's primary tumor and peritoneal metastasis (the other Visium panels are Fig. 7a-e)."""
    F = Fig(180.0, 62.0)
    zs = lambda v: (v - v.mean()) / v.std()
    secs = []
    for tag in ("gc6", "gc6pm"):
        d = rings[rings.GC == ("GC6" if tag == "gc6" else "GC6-PM")].copy()
        d["barrier"] = (zs(d.CAF_ECM.astype(float)) + zs(d.SPP1_arm.astype(float))) / 2
        secs.append(Section(d))
    lo, hi = np.nanpercentile(np.r_[secs[0].d.barrier, secs[1].d.barrier], [2, 98])
    for L, x, sec, ttl in [("a", 6.0, secs[0], "GC6, primary tumor"), ("b", 62.0, secs[1], "GC6, peritoneal metastasis")]:
        ax = sec.axes(F, x, 9.5, 50.0, 40.0)
        sc = sec.values(ax, sec.d.barrier, seq(STROMA), lo, hi); sec.scalebar(ax)
        F.head[L] = btitle(ax, ttl, pad=9, sub="barrier state; outlined, TLS peaks")
    colorbar(F, sc, 8.0, 54.0, lo, hi, "barrier state: mean z of CAF_ECM and the SPP1 subset, within section")
    F.num += [("a-b", "barrier-state colour range (2nd-98th pct, both sections)", f"{lo:.2f} to {hi:.2f}",
               "E27_spots_GSE251950.parquet (GC6, GC6-PM)")]
    ax = F.ax(128.0, 9.5, 48.0, 40.0)
    p = met.set_index("GC_raw")
    comps = [("z_ecm", "CAF_ECM\ncomponent"), ("z_tam", "SPP1 subset\ncomponent"), ("BarrierIndex", "BarrierIndex")]
    for k, (c, lab) in enumerate(comps):
        for j, (g, fill) in enumerate([("GC6", True), ("GC6-PM", False)]):
            v = p.loc[g, c]
            ax.bar(k + (j - 0.5) * 0.36, v, width=0.32, color=STROMA if fill else "white", edgecolor=STROMA, lw=0.7)
            F.num += [("c", f"{g} {c}", f"{v:.2f}", "E27_rings_sections.csv knn9/far3, z across 10 sections")]
    ax.axhline(0, color=INK, lw=0.5)
    ax.set_xticks(range(3)); ax.set_xticklabels([c[1] for c in comps], linespacing=1.0); ax.tick_params(axis="x", length=0)
    ax.set_ylabel("z across 10 sections"); ax.yaxis.set_major_formatter(fmtm)
    vv = p.loc[["GC6", "GC6-PM"], [c for c, _ in comps]].values
    ax.set_ylim(min(-1.0, vv.min() - 0.5), max(1.0, vv.max() + 1.4))
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=STROMA, ec=STROMA, lw=0.7, label="primary tumor"),
                       plt.Rectangle((0, 0), 1, 1, fc="white", ec=STROMA, lw=0.7, label="peritoneal metastasis")],
              loc="upper left", fontsize=5.0, frameon=False, handlelength=1.2, borderaxespad=0.2)
    F.head["c"] = btitle(ax, "One patient, two sites", pad=9, sub="BarrierIndex and its components")
    F.letters([("a", 1.0), ("b", 57.0), ("c", 113.0)])
    F.save("visium")


# ================================================================== FigS_visium_sens
def fig_visium_sens():
    F = Fig(180.0, 150.0)
    others = [s for s in ORDER if GC[s] != "GC8"]
    c = met.set_index("slice")
    for k, s in enumerate(others):
        sec = Section(rings[rings.slice == s].rename(columns={"zone": "ring"}))
        x = 6.0 + (k % 3) * 36.0; y = 12.0 + (k // 3) * 42.0
        ax = sec.axes(F, x, y, 33.0, 30.0)
        sec.zones(ax); sec.scalebar(ax)
        btitle(ax, GC[s], pad=9, sub=f"ring 1 {int(c.loc[s, 'n_Ring1']):,} · far {int(c.loc[s, 'n_Far']):,} spots")
        if k == 0:
            F.head["a"] = F.text(6.0, 4.0, "Ring zones, nine other sections", fontsize=6.5, fontweight="bold",
                                 color=INK, va="baseline")
        F.num += [("a", f"{GC[s]} ring 1 / far spots", f"{int(c.loc[s, 'n_Ring1'])} / {int(c.loc[s, 'n_Far'])}",
                   "E27_rings_sections.csv knn9/far3 n_Ring1, n_Far")]
    zone_legend(F, 6.0, 140.0, 60.0)
    # b: three small multiples
    DEFS = [("lattice", "far3", "lattice peaks", "o"), ("knn9", "far4_noedge", "far ≥ 4, no edge", "s")]
    F.head["b"] = F.text(126.0, 4.0, "Alternative definitions", fontsize=6.5, fontweight="bold", color=INK,
                         va="baseline")
    rng = np.random.default_rng(3)
    for i, (st, name, col) in enumerate(STATES):
        y0 = 22.0 + i * 38.0
        ax = F.ax(130.0, y0, 47.0, 26.0)
        xt, xl = [], []
        for j, ring in enumerate(["Ring1", "Ring2"]):
            for k, (pk, fr, dl, mk) in enumerate(DEFS):
                x = j * 2.4 + k
                g = rs[(rs.peaks == pk) & (rs.far == fr)]
                v = g[f"{st}_{ring}_minus_Far"].values
                ax.scatter(x + rng.uniform(-0.15, 0.15, len(v)), v, s=6, marker=mk, color=col, lw=0,
                           alpha=1.0 if k == 0 else 0.75, zorder=3)
                ax.plot([x - 0.25, x + 0.25], [np.median(v)] * 2, color=INK, lw=1.0, zorder=4)
                P = sm_[(sm_.peaks == pk) & (sm_.far == fr) & (sm_.score == st) & (sm_.contrast == f"{ring}-Far")].P.iloc[0]
                ax.text(x, 1.02, fmt_p(P, lead=""), transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                        fontsize=5.0, color=INK)
                xt.append(x); xl.append("lattice" if k == 0 else "far ≥ 4")
                F.num += [("b", f"{name} {ring}-Far, {dl}: median; P", f"{np.median(v):.3f}; {fmt_p(P, lead='')}",
                           "E27_rings_summary.csv")]
            if i == 0:
                ax.text(j * 2.4 + 0.5, 1.17, f"ring {j + 1} − far", transform=ax.get_xaxis_transform(), ha="center",
                        va="bottom", fontsize=5.6, color=INK)
        S.hrule(ax, 0)
        ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=5.0); ax.tick_params(axis="x", length=0)
        ax.set_xlim(-0.5, 3.9)
        ax.yaxis.set_major_formatter(fmtm); ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        chip(F, 121.0, y0 - (9.5 if i == 0 else 5.5), col, name)
        if i == 1:
            ax.set_ylabel("ring minus far (per section)")
        ax.text(-0.02, 1.02, "P", transform=ax.transAxes, ha="right", va="bottom", fontsize=5.0, color=SLATE)
    F.letters([("a", 1.0), ("b", 121.0)])
    F.save("visium_sens")


# ================================================================== FigS_barrierindex
def fig_barrierindex():
    F = Fig(180.0, 57.0)
    g = met.merge(clin[["slice", "Age", "pT", "pN"]], on="slice")
    g = g[g.is_metastasis == 0].copy()
    g["pT_o"] = g.pT.str.extract(r"T(\d)").astype(float); g["pN_o"] = g.pN.str.extract(r"N(\d)").astype(float)
    rng = np.random.default_rng(2)
    for k, (v, lab, ttl) in enumerate([("Age", "age (years)", "Age"), ("pT_o", "pT", "Tumor stage (pT)"),
                                       ("pN_o", "pN", "Nodal stage (pN)")]):
        x0 = 14.0 + k * 58.0
        ax = F.ax(x0, 8.5, 42.0, 38.0)
        r, p = spearmanr(g.BarrierIndex, g[v])
        xx = g[v] + (rng.uniform(-0.08, 0.08, len(g)) if v != "Age" else 0)
        ax.scatter(xx, g.BarrierIndex, s=12, color=STROMA, lw=0, zorder=3)
        for xv, yv, nm in zip(xx, g.BarrierIndex, g.GC_raw):
            ax.text(xv + (0.12 if v != "Age" else 0.6), yv, nm, fontsize=5.0, va="center", color=SLATE)
        S.hrule(ax, 0)
        if v != "Age":
            t = sorted(g[v].unique()); ax.set_xticks(t); ax.set_xticklabels([f"{lab[1]}{int(q)}" for q in t])
            ax.set_xlim(min(t) - 0.4, max(t) + 0.7)
        else:
            ax.set_xlim(46, 74)
        ax.set_xlabel(lab); ax.yaxis.set_major_formatter(fmtm)
        if k == 0:
            ax.set_ylabel("BarrierIndex (sum of z)")
        F.head["abc"[k]] = btitle(ax, ttl, pad=9, sub=minus(f"Spearman ρ = {r:.2f}, {fmt_p(p)}"))
        F.num += [("abc"[k], f"Spearman rho, P (BarrierIndex vs {v}, 9 primaries)", f"{r:.2f}, {fmt_p(p, lead='')}",
                   "E27_rings_sections.csv BarrierIndex, Spearman")]
    F.letters([("a", 1.0), ("b", 59.0), ("c", 117.0)])
    F.save("barrierindex")


# ================================================================== FigS_moran
def fig_moran():
    F = Fig(180.0, 87.0)
    mb = pd.read_csv(BENCH / "spatial_benchmark_results.csv")
    st = [("moran_state_CAF_ECM", "CAF_ECM", STROMA), ("moran_state_CAF_TGFB", "TGF-β CAF", STROMA),
          ("moran_state_STF3_FAP_INHBA", "FAP⁺INHBA⁺ CAF", STROMA), ("moran_state_SPP1_TAM", "SPP1⁺ TAM", MACRO),
          ("moran_state_IL10_STAT3", "IL-10/STAT3", SLATE)]
    M = mb.set_index(mb.slice.map(GC))[[c for c, _, _ in st]]
    M = M.loc[[GC[s] for s in ORDER]]
    ax = F.ax(16.0, 9.0, 50.0, 62.0)
    im = ax.imshow(M.values, cmap=S.SEQ, vmin=0, vmax=0.7, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M.values[i, j]:.2f}", ha="center", va="center", fontsize=5.0,
                    color="white" if M.values[i, j] > 0.42 else INK)
    ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index); ax.set_xticks([])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for j, (_, lab, col) in enumerate(st):                 # column labels with family chips, rotated
        xc = 16.0 + (j + 0.5) * 10.0
        F.fig.add_artist(Line2D([(xc - 1.2) / F.W, (xc + 1.2) / F.W], [1 - 72.6 / F.H] * 2, color=col, lw=2.2,
                                solid_capstyle="butt"))
        F.text(xc + 0.8, 74.6, lab, color=INK, rotation=40, ha="right", va="top", rotation_mode="anchor")
    cax = F.ax(68.5, 25.0, 1.8, 30.0)
    cb = F.fig.colorbar(im, cax=cax); cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=5.0, length=1.2, width=0.3, pad=0.8); cb.set_label("Moran's I", fontsize=5.0)
    nsig, ntest = int(mb.n_moran_significant.sum()), int(mb.n_states_tested.sum())
    F.head["a"] = btitle(ax, "Spatial autocorrelation", pad=9, sub=f"Moran's I; {nsig} of {ntest} P < 0.05")
    F.num += [("a", "Moran's I significant / tested", f"{nsig} / {ntest}", "spatial_benchmark_results.csv")]
    F.num += [("a", f"Moran's I {g} {lab}", f"{M.loc[g, c]:.2f}", "spatial_benchmark_results.csv")
              for g in M.index for c, lab, _ in st]
    # b: feature comparison
    ax = F.ax(94.0, 9.0, 84.0, 66.0); ax.set_xlim(0, 1); ax.set_ylim(0, 1); S.clean(ax)
    tools = ["barriernet", "Squidpy", "BANKSY", "Spatial-LDA", "CellCharter"]
    rows = [("Primary output", ["marker axes", "statistics,\ngraphs", "domains", "topics", "niches"]),
            ("Needs spatial coordinates", ["optional", "yes", "yes", "yes", "yes"]),
            ("Scores bulk RNA or protein", ["yes", "no", "no", "no", "no"]),
            ("Scores an IHC/IF panel", ["yes", "no", "no", "no", "no"]),
            ("Marker-level interpretability", ["yes", "partial", "no", "partial", "partial"]),
            ("Applies across cohorts unchanged", ["yes", "no", "no", "no", "no"]),
            ("Spatial domain detection", ["via rings", "no", "yes", "yes", "yes"]),
            ("Spatial autocorrelation", ["via rings", "yes", "built in", "implicit", "implicit"]),
            ("Neighborhood enrichment", ["via rings", "yes", "indirect", "indirect", "yes"])]
    cw = [0.36] + [0.128] * 5
    xs = np.cumsum([0] + cw[:-1]); nr = len(rows) + 1; h = 1.0 / nr
    for j, cname in enumerate(["Feature"] + tools):
        ax.text(xs[j] + (0 if j == 0 else cw[j] / 2), 1 - h / 2, cname, fontsize=5.0, fontweight="bold",
                va="center", ha="left" if j == 0 else "center", color=INK)
    ax.plot([0, 1], [1 - h, 1 - h], color=INK, lw=0.5)
    for i, (feat, vals) in enumerate(rows):
        yy = 1 - (i + 1.5) * h
        if i % 2 == 0:
            ax.add_patch(plt.Rectangle((0, yy - h / 2), 1, h, color=PAPER, lw=0))
        ax.text(xs[0] + 0.01, yy, feat, fontsize=5.0, va="center", color=INK)
        for j, v in enumerate(vals, start=1):
            ax.text(xs[j] + cw[j] / 2, yy, v, fontsize=5.0, va="center", ha="center", linespacing=1.0,
                    color=GREY if v == "no" else INK, fontweight="bold" if v == "yes" else "normal")
    F.head["b"] = btitle(ax, "Feature comparison with spatial analysis tools", pad=9,
                         sub="yes in bold; no in gray")
    F.letters([("a", 1.0), ("b", 87.0)])
    F.save("moran")


# ================================================================== FigS_transport
def fig_transport():
    """Marker-based detour ratios and the parameter sweep (Tangram vs marker z is Fig. 7f)."""
    F = Fig(120.0, 61.0)
    t6 = pd.read_csv(R2 / "T6_marker_based.csv"); t4 = pd.read_csv(R2 / "T4_sensitivity.csv")
    t6["GC"] = t6.section.map(GC)
    j = t6.set_index("GC").loc[[GC[s] for s in ORDER]]
    ax = F.ax(17.0, 11.5, 50.0, 40.0)
    y = np.arange(len(j))[::-1]
    for yy, a_, b_ in zip(y, j.null, j.pos_ctrl):
        ax.plot([a_, b_], [yy, yy], color=LIGHT, lw=0.5, zorder=1)
    ax.scatter(j.null, y, s=14, marker="|", color=GREY, lw=0.9, zorder=2)
    ax.scatter(j.detour, y, s=10, color=INK, lw=0, zorder=3)
    ax.scatter(j.pos_ctrl, y, s=10, marker="D", color=STROMA, lw=0, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(j.index); ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("detour ratio (marker scores)")
    ax.legend(handles=[Line2D([], [], ls="", marker="|", color=GREY, ms=4, mew=0.9, label="null (translated)"),
                       Line2D([], [], ls="", marker="o", color=INK, ms=2.8, label="observed"),
                       Line2D([], [], ls="", marker="D", color=STROMA, ms=2.6, label="synthetic shell")],
              loc="lower left", bbox_to_anchor=(-10.0 / 50.0, 1.0), ncol=3, fontsize=5.0, frameon=False, handletextpad=0.2,
              columnspacing=0.8, borderaxespad=0.1)
    F.head["a"] = btitle(ax, "Detour ratios by section", pad=11, x=-10.0 / 50.0)
    F.num += [("a", f"{g} null / observed / shell", f"{r.null:.2f} / {r.detour:.2f} / {r.pos_ctrl:.2f}", "T6_marker_based.csv")
              for g, r in j.iterrows()]
    ax = F.ax(91.0, 11.5, 26.0, 40.0)
    ax.scatter(np.random.default_rng(4).uniform(-0.15, 0.15, len(t4)), t4.median_z, s=8, color=STROMA, lw=0)
    ax.plot([-0.3, 0.3], [t4.median_z.median()] * 2, color=INK, lw=1.1)
    S.hrule(ax, 0)
    ax.set_xticks([0]); ax.set_xticklabels(["27 settings"]); ax.set_xlim(-0.6, 0.6); ax.set_ylim(-4.8, 0.4)
    ax.tick_params(axis="x", length=0); ax.yaxis.set_major_formatter(fmtm)
    ax.set_ylabel("median z across sections (Tangram)")
    F.head["b"] = btitle(ax, "Parameter sweep", pad=9.4, x=-9.0 / 26.0, sub=f"{int((t4.median_z < 0).sum())} of {len(t4)} below 0")
    F.num += [("b", "settings with median z < 0; range; median", f"{int((t4.median_z < 0).sum())}/{len(t4)}; "
               f"{t4.median_z.min():.2f} to {t4.median_z.max():.2f}; {t4.median_z.median():.2f}", "T4_sensitivity.csv median_z")]
    F.letters([("a", 1.0), ("b", 79.0)])
    F.save("transport")


# ================================================================== FigS_imc
def fig_imc():
    """METABRIC per-patient association and survival (maps, positive control, decomposition and
    aggregate densities are Fig. 7g-j)."""
    MO = pd.read_csv(RES / "E13_models.csv"); SV = pd.read_csv(RES / "E13_survival.csv")
    PM = pd.read_csv(RES / "E13_patient_metrics.csv")
    F = Fig(180.0, 62.0)
    ax = F.ax(16.0, 11.5, 50.0, 40.0)
    d = PM[PM.dataset == "METABRIC"].copy()
    for c_ in ("cont", "logS", "logCD8"):
        d["z_" + c_] = (d[c_] - d[c_].mean()) / d[c_].std()
    X = sm.add_constant(d[["z_logS", "z_logCD8"]])
    ry = sm.OLS(d.lir, X).fit().resid; rx = sm.OLS(d.z_cont, X).fit().resid
    ax.scatter(rx, ry, s=3, c=GREY, alpha=0.6, lw=0, rasterized=True)
    b = np.polyfit(rx, ry, 1); xx = np.linspace(rx.min(), rx.max(), 50)
    ax.plot(xx, np.polyval(b, xx), color=STROMA, lw=1.1)
    o = MO[(MO.dataset == "METABRIC") & (MO.spec == "continuity") & (MO.term == "cont")].iloc[0]
    assert len(d) == int(o.n)
    ax.text(0.98, 0.97, minus(f"{o.coef:+.2f} per SD ({o.lo:+.2f} to {o.hi:+.2f})") + f"\n{fmt_p(o.p)}; n = {int(o.n)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.0, color=INK)
    ax.set_xlabel("raw continuity (SD), residual on fibroblast and CD8 density")
    ax.set_ylabel("log₂ CD8 in/out, residual"); ax.set_ylim(min(ry.min(), -5) - 0.3, max(ry.max(), 4) + 1.6)
    ax.xaxis.set_major_formatter(fmtm); ax.yaxis.set_major_formatter(fmtm)
    F.head["a"] = btitle(ax, "Continuity and CD8 infiltration, per patient", pad=9, sub="METABRIC imaging mass cytometry, 361 patients")
    F.num += [("a", "METABRIC continuity coefficient (CI); P; n", f"{o.coef:+.2f} ({o.lo:+.2f} to {o.hi:+.2f}); {fmt_p(o.p, lead='')}; {int(o.n)}",
               "E13_models.csv")]
    ax = F.ax(126.0, 11.5, 32.0, 40.0)
    dd = SV[SV.endpoint == "dss"]
    rows = [("cont (+age, NPI, ER)", "cont", "raw continuity", STROMA),
            ("cont_exc (+age, NPI, ER)", "cont_exc", "excess continuity", STROMA),
            ("logS (+age, NPI, ER)", "logS", "fibroblast density", STROMA),
            ("lir (+age, NPI, ER)", "lir", "CD8 in/out ratio", IMMUNE),
            ("logCD8 (+age, NPI, ER)", "logCD8", "CD8 density", IMMUNE),
            ("joint: cont + logS + logCD8 (+age, NPI, ER)", "cont", "raw continuity, joint model", STROMA)]
    tr = ax.get_yaxis_transform()
    for i, (m, t, lab, col) in enumerate(rows):
        r = dd[(dd.model == m) & (dd.term == t)].iloc[0]
        ax.plot([r.lo, r.hi], [-i] * 2, color=INK, lw=0.9)
        ax.scatter(r.HR, -i, s=11, color=col, lw=0, zorder=3)
        ax.text(1.06, -i, f"{r.HR:.2f}", transform=tr, fontsize=5.0, ha="left", va="center", clip_on=False)
        ax.text(1.30, -i, fmt_p(r.p, lead=""), transform=tr, fontsize=5.0, ha="left", va="center", clip_on=False)
        F.num += [("b", f"DSS HR per SD, {lab}", f"{r.HR:.2f} ({r.lo:.2f} to {r.hi:.2f}); P {fmt_p(r.p, lead='')}", "E13_survival.csv endpoint=dss")]
    rowlabels(ax, -np.arange(len(rows)), [r_[2] for r_ in rows], [])
    ax.text(1.06, 0.8, "HR", transform=tr, fontsize=5.0, ha="left", color=SLATE)
    ax.text(1.30, 0.8, "P", transform=tr, fontsize=5.0, ha="left", color=SLATE)
    S.rule(ax, 1)
    ax.set_xscale("log"); ax.set_xlim(0.72, 1.7); ax.set_xticks([0.8, 1, 1.25, 1.6])
    ax.set_xticklabels(["0.8", "1", "1.25", "1.6"]); ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_ylim(-len(rows) + 0.4, 1.5)
    r0 = dd.iloc[0]
    ax.set_xlabel("HR per SD, adjusted for age, NPI and ER")
    F.head["b"] = btitle(ax, "Disease-specific survival", pad=9, x=-25.0 / 32.0,
                         sub=f"METABRIC; {int(r0.n)} patients, {int(r0.events)} deaths")
    F.letters([("a", 1.0), ("b", 94.0)])
    F.save("imc")


FIGS = {"visium": fig_visium, "visium_sens": fig_visium_sens, "barrierindex": fig_barrierindex,
        "moran": fig_moran, "transport": fig_transport, "imc": fig_imc}
if __name__ == "__main__":
    for k in (sys.argv[1:] or FIGS):
        FIGS[k]()
