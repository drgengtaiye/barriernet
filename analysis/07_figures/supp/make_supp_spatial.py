"""Supplementary Figs S5-S8 (spatial): ring zones on the corrected lattice for the other nine
sections and the stricter Far definition; corrected BarrierIndex vs clinical variables; Moran's I
and feature comparison; least-cost transport with Tangram fractions and independent markers."""
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import spearmanr, wilcoxon
from matplotlib.patches import Patch
import bnstyle as bs, figstyle as fs
from bnstyle import INK, GREY, STROMA, MACRO, NULL, title, fmt_p

F = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures"; S = f"{F}/supp"
G = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
R2 = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
clin = pd.read_csv(f"{G}/results/ring_barrier_state_v2/BarrierIndex_by_slice_with_clinical.csv")
GC = dict(zip(clin.slice, clin.GC_raw))
ZC = {"Hotspot": "#0A3D6E", "Ring1": "#0F6FC5", "Ring2": "#9CC3E6", "Far": "#E3E8EC", "Other": "#F6F7F8"}
PROG = [("state_SPP1_TAM", "SPP1⁺ TAM", MACRO), ("state_CAF_ECM", "CAF_ECM", STROMA), ("state_IL10_STAT3", "IL-10/STAT3", GREY)]
REPORT = {}


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(f"{S}/{stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)


def minus(s): return s.replace("-", "−")


# ---------------- S5 ----------------
rings = pd.read_parquet(f"{F}/data/rings_knn9_far3.parquet")
SENS = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/code/analysis/06_spatial"
rs = pd.read_csv(f"{SENS}/ring_peaks_sensitivity.csv").rename(columns={"section": "slice"})
sm = pd.read_csv(f"{SENS}/ring_peaks_sensitivity_summary.csv")
rc = rs.query("peaks == 'knn9' and far == 'far3'").copy()
secs = [s for s in sorted(rings.slice.unique(), key=lambda s: (GC[s].replace("GC", "").replace("-PM", ".5").zfill(4))) if GC[s] != "GC8"]
fig = plt.figure(figsize=(bs.W2, 150 * bs.MM))
gsp = fig.add_gridspec(3, 5, width_ratios=[1, 1, 1, 0.32, 1.6], hspace=0.3, wspace=0.1)
for k, s in enumerate(secs):
    ax = fig.add_subplot(gsp[k // 3, k % 3])
    r = rings[rings.slice == s]
    ax.scatter(r.x, r.y, s=0.9, c=[ZC.get(z, ZC["Other"]) for z in r.zone], marker="h", lw=0, rasterized=True)
    ax.set_aspect("equal"); ax.set_axis_off()
    n = rc[rc.slice == s].iloc[0]
    title(ax, GC[s])
    ax.text(0.0, -0.04, f"Ring1 {int(n.n_Ring1):,} · Far {int(n.n_Far):,} spots", transform=ax.transAxes, fontsize=5.0,
            color=GREY, va="top")
    if k == 0: bs.panel(ax, "a", dx=-0.16, dy=1.1)
fig.legend(handles=[Patch(color=ZC[z], label=z) for z in ["Hotspot", "Ring1", "Ring2", "Far"]], loc="lower left",
           bbox_to_anchor=(0.12, 0.0), ncol=4, fontsize=5.4)
ax = fig.add_subplot(gsp[:, 4])
DEFS = [("lattice", "far3", "lattice peaks (upper P)", "o"), ("knn9", "far4_noedge", "Far ≥ 4 hops, no edge (lower P)", "s")]
x = 0; ticks, labs = [], []
for p, lab, col in PROG:
    short = p.replace("state_", "")
    for rr in ["Ring1", "Ring2"]:
        for k, (pk, fr, dl, mk) in enumerate(DEFS):
            g = rs[(rs.peaks == pk) & (rs.far == fr)]
            v = g[f"{p}_{rr}_minus_Far"].values
            xx = x + (k - 0.5) * 0.42
            jit = np.random.default_rng(3 + k).uniform(-0.08, 0.08, len(v))
            ax.scatter(xx + jit, v, s=6, color=col, marker=mk, lw=0, alpha=0.9 if k == 0 else 0.6)
            ax.plot([xx - 0.14, xx + 0.14], [np.median(v)] * 2, color=INK, lw=0.9)
            P = sm[(sm.peaks == pk) & (sm.far == fr) & (sm.prog == short) & (sm.contrast == f"{rr}-Far")].P.iloc[0]
            ax.text(xx, 0.236 if k == 0 else 0.216, f"{P:.2g}", ha="center", fontsize=5.0)
        ticks.append(x); labs.append(f"{rr[0]}{rr[-1]}−Far")
        x += 1.35
    ax.text(x - 2.0, -0.405, lab, ha="center", fontsize=5.4, color=col if col != GREY else INK)
    x += 0.45
ax.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
ax.set_xticks(ticks); ax.set_xticklabels(labs, fontsize=5.0, rotation=45)
ax.set_ylim(-0.34, 0.25); ax.set_xlim(left=-0.95); ax.set_ylabel("score next to TLS minus far (per section)")
ax.text(-1.5, 0.226, "P", fontsize=5.0, fontstyle="italic", ha="center")
ax.legend(handles=[plt.Line2D([], [], marker=m, ls="", color=GREY, ms=3, label=d) for _, _, d, m in DEFS],
          loc="lower right", fontsize=5.0, borderaxespad=0.2)
title(ax, "Sensitivity definitions")
bs.panel(ax, "b", dx=-0.3, dy=1.05)
save(fig, "FigS5")

# ---------------- S6 ----------------
g3 = rc.merge(clin[["slice", "GC_raw", "is_metastasis", "Age", "pT", "pN"]], on="slice")
g3 = g3[g3.is_metastasis == 0].copy()
g3["pT_o"] = g3.pT.str.extract(r"T(\d)").astype(float); g3["pN_o"] = g3.pN.str.extract(r"N(\d)").astype(float)
fig, axs = plt.subplots(1, 3, figsize=(bs.W2 * 0.85, 55 * bs.MM))
for ax, (v, lab, L) in zip(axs, [("Age", "age (years)", "a"), ("pT_o", "pT", "b"), ("pN_o", "pN", "c")]):
    r, p = spearmanr(g3.BarrierIndex, g3[v])
    ax.scatter(g3[v] + (np.random.default_rng(2).uniform(-0.08, 0.08, len(g3)) if v != "Age" else 0), g3.BarrierIndex,
               s=14, color=STROMA, edgecolor="white", lw=0.4)
    for _, q in g3.iterrows():
        ax.text(q[v] + (0.12 if v != "Age" else 0.6), q.BarrierIndex, q.GC_raw, fontsize=5.0, va="center", color=GREY)
    ax.axhline(0, color=GREY, lw=0.4, ls=(0, (2, 2)))
    if v != "Age": ax.set_xticks(sorted(g3[v].unique())); ax.set_xticklabels([f"{lab[1]}{int(t)}" for t in sorted(g3[v].unique())])
    ax.set_xlabel(lab); ax.set_ylabel("BarrierIndex" if L == "a" else "")
    title(ax, minus(f"Spearman ρ = {r:.2f}, {fmt_p(p)}")); bs.panel(ax, L, dx=-0.25)
    REPORT[f"S6_{v}"] = (round(r, 2), round(p, 2))
fig.subplots_adjust(wspace=0.4)
save(fig, "FigS6")

# ---------------- S7 ----------------
mb = pd.read_csv(f"{G}/validation/bench/spatial_benchmark_results.csv")
st = [("moran_state_CAF_ECM", "CAF_ECM"), ("moran_state_CAF_TGFB", "TGF-β CAF"), ("moran_state_STF3_FAP_INHBA", "FAP⁺INHBA⁺ CAF"),
      ("moran_state_SPP1_TAM", "SPP1⁺ TAM"), ("moran_state_IL10_STAT3", "IL-10/STAT3")]
M = mb.set_index(mb.slice.map(GC))[[c for c, _ in st]]
M = M.loc[sorted(M.index, key=lambda s: s.replace("GC", "").replace("-PM", ".5").zfill(4))]
fig = plt.figure(figsize=(bs.W2, 100 * bs.MM))
gsp = fig.add_gridspec(1, 2, width_ratios=[0.8, 1.6], wspace=0.35)
ax = fig.add_subplot(gsp[0])
im = ax.imshow(M.values, cmap=fs.SEQ, vmin=0, vmax=0.7, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f"{M.values[i, j]:.2f}", ha="center", va="center", fontsize=5.0, color="white" if M.values[i, j] > 0.45 else INK)
ax.set_xticks(range(len(st))); ax.set_xticklabels([l for _, l in st], rotation=45, ha="right", fontsize=5.2)
ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=5.2); ax.tick_params(length=0)
for sp_ in ax.spines.values(): sp_.set_visible(False)
cb = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.03); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.0, length=1.5)
cb.set_label("Moran's I", fontsize=5.2)
title(ax, f"Moran's I; {int(mb.n_moran_significant.sum())} of {int(mb.n_states_tested.sum())} P < 0.05")
bs.panel(ax, "a", dx=-0.35, dy=1.07)
ax = fig.add_subplot(gsp[1]); ax.set_axis_off()
tools = ["barriernet", "Squidpy", "BANKSY", "SpatialLDA", "CellCharter"]
rows = [("Primary output", ["marker axes", "statistics,\ngraphs", "domains", "topics", "niches"]),
        ("Needs spatial coordinates", ["optional", "yes", "yes", "yes", "yes"]),
        ("Scores bulk RNA or protein", ["yes", "no", "no", "no", "no"]),
        ("Scores an IHC/IF panel", ["yes", "no", "no", "no", "no"]),
        ("Marker-level interpretability", ["yes", "partial", "no", "partial", "partial"]),
        ("Applies across cohorts unchanged", ["yes", "no", "no", "no", "no"]),
        ("Spatial domain detection", ["via rings", "no", "yes", "yes", "yes"]),
        ("Spatial autocorrelation", ["via rings", "yes", "built in", "implicit", "implicit"]),
        ("Neighbourhood enrichment", ["via rings", "yes", "indirect", "indirect", "yes"])]
x0, cw = 0.0, [0.34] + [0.132] * 5
xs = np.cumsum([x0] + cw[:-1]); nr = len(rows) + 1; h = 1.0 / nr
for jj, c in enumerate(["Feature"] + tools):
    ax.text(xs[jj] + (0 if jj == 0 else cw[jj] / 2), 1 - h / 2, c, fontsize=5.3, fontweight="bold", va="center",
            ha="left" if jj == 0 else "center", transform=ax.transAxes)
for ii, (feat, vals) in enumerate(rows):
    yy = 1 - (ii + 1.5) * h
    if ii % 2 == 0:
        ax.add_patch(plt.Rectangle((0, yy - h / 2), sum(cw), h, color=bs.PAPER, transform=ax.transAxes, lw=0))
    ax.text(xs[0], yy, feat, fontsize=5.1, va="center", transform=ax.transAxes)
    for jj, v in enumerate(vals, start=1):
        ax.text(xs[jj] + cw[jj] / 2, yy, v, fontsize=5.0, va="center", ha="center", transform=ax.transAxes,
                color=STROMA if v == "yes" else (GREY if v == "no" else INK), fontweight="bold" if v == "yes" else "normal")
title(ax, "Feature comparison with spatial analysis tools")
bs.panel(ax, "b", dx=-0.02, dy=1.07)
save(fig, "FigS7")
REPORT["S7_moran_median"] = {l: round(M[c].median(), 3) for c, l in st}

# ---------------- S8 ----------------
t3 = pd.read_csv(f"{R2}/T3_transport.csv"); t6 = pd.read_csv(f"{R2}/T6_marker_based.csv"); t4 = pd.read_csv(f"{R2}/T4_sensitivity.csv")
t3["GC"] = t3.section.map(GC); t6["GC"] = t6.section.map(GC)
j = t3.merge(t6, on="GC", suffixes=("_tg", "_mk"))
j = j.set_index("GC").loc[sorted(j.GC, key=lambda s: s.replace("GC", "").replace("-PM", ".5").zfill(4))]
fig, axs = plt.subplots(1, 3, figsize=(bs.W2, 62 * bs.MM), gridspec_kw=dict(width_ratios=[1.1, 1.3, 0.8]))
ax = axs[0]
for i, (g, r) in enumerate(j.iterrows()):
    ax.plot([0, 1], [r.z_shift, r.z], color=NULL, lw=0.7)
ax.scatter(np.zeros(len(j)), j.z_shift, s=12, color=STROMA, zorder=3, label="Tangram fractions")
ax.scatter(np.ones(len(j)), j.z, s=12, color=INK, zorder=3, label="marker scores")
ax.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
ax.set_xticks([0, 1]); ax.set_xticklabels(["Tangram\nfractions", "independent\nmarker scores"]); ax.set_xlim(-0.4, 1.4)
ax.set_ylabel("z, observed vs null (positive = obstruction)")
pw = wilcoxon(j.z, alternative="greater").pvalue
title(ax, minus(f"median z {j.z_shift.median():.2f} vs +{j.z.median():.2f}") + f"\nmarker scores, one-sided {fmt_p(pw)}")
bs.panel(ax, "a", dx=-0.32)
ax = axs[1]
y = np.arange(len(j))[::-1]
ax.scatter(j.null, y, s=12, marker="|", color=GREY, label="null (repositioned)")
ax.scatter(j.detour_mk, y, s=12, color=INK, label="observed")
ax.scatter(j.pos_ctrl_mk, y, s=12, marker="D", color=STROMA, label="positive control (shell)")
ax.set_yticks(y); ax.set_yticklabels(j.index, fontsize=5.2); ax.tick_params(axis="y", length=0)
ax.set_xlabel("detour ratio (marker scores)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3, fontsize=5.0, columnspacing=0.6, handletextpad=0.2)
title(ax, "Marker-based detour ratio by section"); bs.panel(ax, "b", dx=-0.22)
ax = axs[2]
ax.scatter(np.random.default_rng(4).uniform(-0.15, 0.15, len(t4)), t4.median_z, s=10, color=STROMA, lw=0)
ax.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
ax.set_xticks([0]); ax.set_xticklabels(["27 settings"]); ax.set_xlim(-0.6, 0.6)
ax.set_ylabel("median z across sections (Tangram)")
title(ax, "Parameter sweep"); bs.panel(ax, "c", dx=-0.55)
fig.subplots_adjust(wspace=0.55)
save(fig, "FigS8")
REPORT["S8"] = dict(tangram_median=round(j.z_shift.median(), 2), marker_median=round(j.z.median(), 2),
                    marker_pos=int((j.z > 0).sum()), wilcoxon_greater=round(pw, 3), sweep_neg=int((t4.median_z < 0).sum()))
print(REPORT)
