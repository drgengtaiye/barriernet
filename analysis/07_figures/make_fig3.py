"""v14 Figure 3 - spatial organisation of the modules.

Programme scores and TLS peaks come from the ring_barrier_state_v2 outputs, but ring zones,
ring contrasts and BarrierIndex are RECOMPUTED on the true Visium hexagonal neighbourhood
(ring_recompute.py): the stored zones were built with axial hex offsets that degenerate into
one-diagonal chains on Visium array coordinates.
  a-c  one representative primary tumour (the primary whose corrected SPP1_TAM Ring1-Far is
       closest to the median of the 10 sections)
  d    per-section ring-minus-far contrasts across all 10 sections, paired Wilcoxon vs 0
  e    BarrierIndex by section with clinical annotation
  f-h  the one patient with a primary tumour (GC6) and matched peritoneal metastasis (GC6-PM)
Ring zones: Hotspot = TLS peak, Ring1/Ring2 = 1/2 hops, Far = >= 3 hops (as in the stored files;
ring_recompute.py also reports >= 4 hops excluding tissue edges, with the same conclusions).
"""
import sys
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import wilcoxon
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import bnstyle as bs
from bnstyle import INK, GREY, LIGHT, PAPER, IMMUNE, STROMA, MACRO, fmt_p, title
from rowpack import pack

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2"
PRI, MET = "GSM7990478_21_01251_LI_SING", "GSM7990479_21_01252_LI_SING"
ZONE_COL = {"Hotspot": "#0A3D6E", "Ring1": IMMUNE, "Ring2": "#9CC3E6", "Far": "#E6EAEE"}
seq = lambda c: LinearSegmentedColormap.from_list("s", ["#F1F3F5", c])

clin = pd.read_csv(f"{R}/BarrierIndex_by_slice_with_clinical.csv").drop(columns=["z_ecm", "z_tam", "z_bar", "BarrierIndex"])
# Primary definition (Methods): TLS peaks = spots above the 95th percentile of the TLS score that
# are not smaller than any of their 9 nearest spots; rings on the true Visium lattice; Far >= 3 hops.
# Section-level contrasts from code/analysis/06_spatial/ring_peaks_sensitivity.py (barriernet 0.2.1).
SENS = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/code/analysis/06_spatial/ring_peaks_sensitivity.csv"
met = pd.read_csv(SENS).query("peaks == 'knn9' and far == 'far3'").rename(columns={"section": "slice"}).reset_index(drop=True)
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/code/barriernet")
from barriernet.spatial import assign_rings
z0 = lambda v: (v - v.mean()) / v.std(ddof=0)
met["z_ecm"] = z0(met.state_CAF_ECM_Ring1_minus_Far); met["z_tam"] = z0(met.state_SPP1_TAM_Ring1_minus_Far)
bi = met[["slice", "z_ecm", "z_tam", "BarrierIndex"]].merge(clin, on="slice")
lab = dict(zip(bi.slice, bi.GC_raw))
prim = bi[bi.is_metastasis == 0].slice
v = met.set_index("slice").state_SPP1_TAM_Ring1_minus_Far
REP = (v[prim] - v.median()).abs().idxmin()


def load(s):
    a = ad.read_h5ad(f"{R}/{s}.with_state_and_ring.h5ad", backed="r")
    o = a.obs.copy()
    xy = np.asarray(a.obsm["spatial"], float)
    rd = assign_rings(o, xy, "TLS_niche_score", peak_method="knn", peak_k=9, far_min_dist=3)
    o["ring"] = rd.ring.values
    o["x"], o["y"] = xy[:, 0], -xy[:, 1]
    o = o[o.in_tissue.astype(int) == 1].copy()
    return o


SEC = {s: load(s) for s in (REP, PRI, MET)}
print("representative primary:", lab[REP], REP)


def spot_size(ax, o, fill=0.95):
    """Marker area so neighbouring Visium spots nearly touch at the panel's scale."""
    from scipy.spatial import cKDTree
    d = np.median(cKDTree(o[["x", "y"]].values).query(o[["x", "y"]].values, k=2)[0][:, 1])
    ax.figure.canvas.draw()
    bb = ax.get_window_extent()
    span = max(np.ptp(o.x) / bb.width, np.ptp(o.y) / bb.height)   # data units per pixel
    px = d / span * fill
    return (px * 72 / ax.figure.dpi) ** 2


def spatial(ax, o, col=None, cmap=None, zones=False, ttl="", cbar=None, ring_outline=True):
    ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_xlim(o.x.min() - 60, o.x.max() + 60); ax.set_ylim(o.y.min() - 60, o.y.max() + 60)
    s = spot_size(ax, o)
    if zones:
        for z in ["Other", "Far", "Ring2", "Ring1", "Hotspot"]:
            if z not in ZONE_COL:
                continue
            k = o.ring == z
            ax.scatter(o.x[k], o.y[k], s=s, c=ZONE_COL[z], linewidths=0, rasterized=True)
        h = [Patch(color=ZONE_COL[z], label=z) for z in ["Hotspot", "Ring1", "Ring2", "Far"]]
        ax.legend(handles=h, loc="lower left", bbox_to_anchor=(0.0, -0.13), ncol=4, fontsize=5.2,
                  handlelength=0.8, columnspacing=0.8, handletextpad=0.3, borderaxespad=0)
    else:
        v = o[col].astype(float)
        lo, hi = np.nanpercentile(v, [2, 98])
        sc = ax.scatter(o.x, o.y, s=s, c=v, cmap=cmap, vmin=lo, vmax=hi, linewidths=0, rasterized=True)
        if ring_outline:
            k = o.ring == "Hotspot"
            ax.scatter(o.x[k], o.y[k], s=s * 0.55, facecolors="none", edgecolors=INK, linewidths=0.25)
        cax = ax.inset_axes([0.05, -0.08, 0.42, 0.035])
        cb = ax.figure.colorbar(sc, cax=cax, orientation="horizontal")
        cb.outline.set_linewidth(0.3); cb.ax.tick_params(labelsize=5.0, length=1.5, width=0.3, pad=1)
        cb.set_ticks([lo, (lo + hi) / 2, hi]); cb.set_ticklabels([f"{t:.2f}".replace("-", "−") for t in (lo, (lo + hi) / 2, hi)])
        if cbar:
            ax.text(0.5, -0.07, cbar, transform=ax.transAxes, fontsize=5.0, va="center", ha="left",
                    color=INK)
    title(ax, ttl, pad=2)


def zscore(v):
    v = v.astype(float); return (v - v.mean()) / v.std()


def barrier_map(ax, o, ttl):
    o = o.copy()
    o["barrier_state"] = (zscore(o.state_CAF_ECM) + zscore(o.state_SPP1_TAM)) / 2
    spatial(ax, o, "barrier_state", seq(STROMA), ttl=ttl, cbar="barrier state (mean z)")


# ---------------------------------------------------------------- d: ring contrasts
PROG = [("state_SPP1_TAM", "SPP1⁺ TAM", MACRO), ("state_CAF_ECM", "CAF_ECM", STROMA),
        ("state_IL10_STAT3", "IL-10/STAT3", GREY)]


def contrasts(ax):
    rng = np.random.default_rng(3)
    ismet = met.slice.map(lab).eq("GC6-PM").values
    xt, xl = [], []
    for i, (col, name, c) in enumerate(PROG):
        for j, ring in enumerate(["Ring1", "Ring2"]):
            x = i * 2.6 + j
            v = met[f"{col}_{ring}_minus_Far"].values
            jit = rng.uniform(-0.18, 0.18, len(v))
            ax.scatter(x + jit[~ismet], v[~ismet], s=7, color=c, edgecolors="white", linewidths=0.3, zorder=3)
            ax.scatter(x + jit[ismet], v[ismet], s=7, facecolors="white", edgecolors=c, linewidths=0.6, zorder=3)
            ax.plot([x - 0.28, x + 0.28], [np.median(v)] * 2, color=INK, lw=0.9, zorder=4)
            p = wilcoxon(v).pvalue
            ax.text(x, 0.225, fmt_p(p, lead="P = "), ha="center", va="bottom", fontsize=5.0, color=INK)
            xt.append(x); xl.append(f"R{j + 1}−Far")
        ax.text(i * 2.6 + 0.5, -0.13, name, ha="center", va="top", fontsize=6.0, color=INK,
                transform=ax.get_xaxis_transform())
    ax.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=5.3)
    lo = min(met[[f"{c}_{r}_minus_Far" for c, _, _ in PROG for r in ("Ring1", "Ring2")]].min().min(), -0.3)
    ax.set_xlim(-0.6, 2 * 2.6 + 1.6); ax.set_ylim(lo - 0.02, 0.26)
    ax.set_ylabel("Score near TLS − far (per section)")
    h = [Line2D([], [], marker="o", ls="", color=GREY, markersize=3, label="primary tumour"),
         Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=GREY, markersize=3,
                label="peritoneal metastasis"),
         Line2D([], [], color=INK, lw=0.9, label="median")]
    ax.legend(handles=h, loc="lower right", fontsize=5.0, ncol=1, borderaxespad=0.2)
    title(ax, "Near TLS minus far, 10 sections (paired Wilcoxon vs 0)")


# ---------------------------------------------------------------- e: BarrierIndex
def barrier_index(ax):
    d = bi.sort_values("BarrierIndex", ascending=False).reset_index(drop=True)
    x = np.arange(len(d))
    cols = [STROMA if m == 0 else "white" for m in d.is_metastasis]
    ax.bar(x, d.BarrierIndex, width=0.72, color=cols, edgecolor=STROMA, linewidth=0.6)
    ax.axhline(0, color=INK, lw=0.5)
    ax.set_xticks(x); ax.set_xticklabels(d.GC_raw, fontsize=5.4)
    ax.set_ylabel("BarrierIndex (sum of z)")
    ax.set_xlim(-0.6, len(d) - 0.4); ax.set_ylim(-3.6, 3.7)
    rows = [("EBV", lambda r: "+" if r.EBV == "positive" else "−"),
            ("MSI", lambda r: "H" if r.MSI == "MSI" else "S"),
            ("Lauren", lambda r: "I" if r.Type == "Intestinal" else "D"),
            ("pM", lambda r: str(int(r.pM)))]
    for k, (name, f) in enumerate(rows):
        yy = -0.17 - 0.075 * k
        ax.text(-0.9, yy, name, transform=ax.get_xaxis_transform(), ha="right", va="center", fontsize=5.0)
        for i, r in d.iterrows():
            ax.text(i, yy, f(r), transform=ax.get_xaxis_transform(), ha="center", va="center", fontsize=5.0)
    title(ax, "BarrierIndex by section (open bar: metastasis)")


# ---------------------------------------------------------------- h: case bars
def case(ax):
    d = bi.set_index("GC_raw").loc[["GC6", "GC6-PM"]]
    comps = [("z_ecm", "CAF_ECM\ncomponent"), ("z_tam", "SPP1⁺ TAM\ncomponent"), ("BarrierIndex", "BarrierIndex")]
    for i, (c, _) in enumerate(comps):
        ax.bar(i - 0.19, d.loc["GC6", c], width=0.36, color=STROMA, edgecolor=STROMA, lw=0.6)
        ax.bar(i + 0.19, d.loc["GC6-PM", c], width=0.36, color="white", edgecolor=STROMA, lw=0.6)
    ax.axhline(0, color=INK, lw=0.5)
    ax.set_xticks(range(3)); ax.set_xticklabels([n for _, n in comps], fontsize=5.4)
    ax.set_ylabel("z (across 10 sections)")
    ax.legend(handles=[Patch(facecolor=STROMA, edgecolor=STROMA, label="GC6 primary"),
                       Patch(facecolor="white", edgecolor=STROMA, label="GC6 peritoneal metastasis")],
              loc="upper left", fontsize=5.0, borderaxespad=0.1)
    ax.set_ylim(-2.4, 5.0)
    title(ax, "One patient, primary vs metastasis")


rows = [
    [dict(letter="a", aspect=1.0, m=(1, 6, 1, 4), draw=lambda ax: spatial(ax, SEC[REP], zones=True,
          ttl=f"{lab[REP]}: TLS-centred ring zones")),
     dict(letter="b", aspect=1.0, m=(1, 6, 1, 4), draw=lambda ax: spatial(ax, SEC[REP], "state_SPP1_TAM", seq(MACRO),
          ttl=f"{lab[REP]}: SPP1⁺ TAM state", cbar="state score; rings = TLS peaks")),
     dict(letter="c", aspect=1.0, m=(1, 6, 1, 4), draw=lambda ax: spatial(ax, SEC[REP], "state_CAF_ECM", seq(STROMA),
          ttl=f"{lab[REP]}: CAF_ECM state", cbar="state score; rings = TLS peaks"))],
    [dict(letter="d", aspect=1.55, m=(10, 11, 1, 4), draw=contrasts),
     dict(letter="e", aspect=1.45, m=(12, 15, 1, 4), draw=barrier_index)],
    [dict(letter="f", aspect=1.0, m=(1, 6, 1, 4), draw=lambda ax: barrier_map(ax, SEC[PRI], "GC6 primary tumour")),
     dict(letter="g", aspect=1.0, m=(1, 6, 1, 4), draw=lambda ax: barrier_map(ax, SEC[MET], "GC6 peritoneal metastasis")),
     dict(letter="h", aspect=1.0, m=(10, 9, 1, 4), draw=case)],
]
if __name__ == "__main__":
    print("height mm", round(pack("Fig3", rows), 1))
