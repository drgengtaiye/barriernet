"""v19 Fig. 2 - the marker modules are expressed by their intended cell types.

a  UMAP of the discovery atlas (GSE183904) by cell class, with the shared cell-class key.
b  Per-cell module scores on the same UMAP (six modules).
c  Module genes x cell class dot plot (discovery atlas).
d  Fibroblast states: UMAP, and program effect sizes, ECM-remodeling CAFs vs iCAFs.
e  Macrophage states: UMAP, and program effect sizes, SPP1+ TAMs vs inflammatory macrophages.
f  Cell-class composition of each tumor sample, by TCGA subtype.
g  Share of module signal by cell class in the discovery and two independent atlases.
h  Module enrichment (z across classes) in the three atlases, rows aligned with g.
i  Check of the B-cell annotation (E17b): reassigned cells, their markers, and the shares.

Additional sources for d-f and i: figures/data/sc/fibro_obs.parquet, macro_obs.parquet,
cells_expr.parquet, sample_clinical_fibro.csv; figures/data/sc/Fig2_effect_sizes.csv (the
effect sizes are recomputed here with the same code as figures/make_fig2.py and asserted equal);
analysis/results/E17b_bcell_reassignment.csv (shares; the reassignment rule is re-applied here
only to locate the cells and count marker detection, and the count is asserted).

Sources: figures/data/sc/cells_obs.parquet (UMAP, classes); analysis/cache_E17/GSE183904_expr.npz
(42 module genes, ln(1+CP10k), written by E17_singlecell_attribution.py --stage extract);
analysis/results/E17_attribution_by_class.csv and E17_gene_expression_by_class.csv.
Per-cell scores are recomputed exactly as E17 cell_scores() does (mean of the module genes
present, >= 3) and cached to figures/cns/data/fig2_cell_module_scores.parquet; the script
checks that their class means reproduce E17_attribution_by_class.csv before drawing.
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import cns_style as S
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

ROOT = S.ROOT
DATA = S.OUT / "data"
DATA.mkdir(parents=True, exist_ok=True)

MODULES = {
    "CAF_ECM": ["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
    "TAM_SUP": ["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
    "APM": ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    "IFN": ["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
    "TLS_B": ["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
    "T_CELL": ["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"],
}
MODS = list(MODULES)

# E17 class names -> display label and colour key in S.CT
CLASSES = ["T cell", "NK cell", "B cell", "plasma cell", "macrophage/monocyte", "dendritic cell",
           "mast cell", "fibroblast", "endothelial", "pericyte/smooth muscle", "epithelial/tumour"]
LAB = {"T cell": "T cell", "NK cell": "NK cell", "B cell": "B cell", "plasma cell": "Plasma cell",
       "macrophage/monocyte": "Macrophage/monocyte", "dendritic cell": "Dendritic cell",
       "mast cell": "Mast cell", "fibroblast": "Fibroblast", "endothelial": "Endothelial",
       "pericyte/smooth muscle": "Pericyte/SMC", "epithelial/tumour": "Epithelial/tumor"}
SHORT = {"T cell": "T", "NK cell": "NK", "B cell": "B", "plasma cell": "Plasma",
         "macrophage/monocyte": "Mac/mono", "dendritic cell": "DC", "mast cell": "Mast",
         "fibroblast": "Fibro", "endothelial": "Endo", "pericyte/smooth muscle": "Peri/SMC",
         "epithelial/tumour": "Epi"}
KEY = {"T cell": "Tcell", "NK cell": "NK", "B cell": "Bcell", "plasma cell": "Plasma",
       "macrophage/monocyte": "Macro", "dendritic cell": "DC", "fibroblast": "Fibro",
       "endothelial": "Endo", "pericyte/smooth muscle": "Peri", "epithelial/tumour": "Epi"}
CCOL = {c: S.CT[k] for c, k in KEY.items()}
CCOL["mast cell"] = "#9A5B13"                    # myeloid family, darker than macrophage / DC
ANNO2CLASS = {v: k for k, v in KEY.items()}
INTENDED = {"CAF_ECM": ["fibroblast"], "TAM_SUP": ["macrophage/monocyte"],
            "TLS_B": ["B cell", "plasma cell"], "T_CELL": ["T cell", "NK cell"]}
DATASETS = ["GSE183904", "OMIX001073", "GSE206785"]
DSLAB = {"GSE183904": "Discovery", "OMIX001073": "OMIX001073", "GSE206785": "GSE206785"}

# ------------------------------------------------------------------ data
obs = pd.read_parquet(ROOT / "figures/data/sc/cells_obs.parquet")
obs["class"] = obs.anno_level1.map(ANNO2CLASS)
cache = DATA / "fig2_cell_module_scores.parquet"
if cache.exists():
    sc = pd.read_parquet(cache)
else:
    z = np.load(ROOT / "analysis/cache_E17/GSE183904_expr.npz", allow_pickle=False)
    o = pd.read_csv(ROOT / "analysis/cache_E17/GSE183904_obs.csv.gz", usecols=["cell", "class"])
    assert (o.cell.values == obs.index.values).all()
    genes, X = list(z["genes"]), z["harm"]
    sc = pd.DataFrame({m: X[:, [genes.index(g) for g in MODULES[m] if g in genes]].mean(1, dtype=np.float64)
                       for m in MODS}, index=obs.index)
    sc.to_parquet(cache)
attr = pd.read_csv(S.RES / "E17_attribution_by_class.csv")
attr = attr[attr.analysis == "primary"]
chk = sc.groupby(obs["class"].values, observed=True).mean()
ref = attr[attr.dataset == "GSE183904"].pivot(index="class", columns="module", values="mean_score")
assert np.allclose(chk.loc[ref.index, MODS].values, ref[MODS].values, atol=1e-6), "per-cell scores do not reproduce E17"

gex = pd.read_csv(S.RES / "E17_gene_expression_by_class.csv")
gex = gex[gex.dataset == "GSE183904"]
present = [g for m in MODS for g in MODULES[m] if g in set(gex.gene)]
missing = [g for m in MODS for g in MODULES[m] if g not in set(gex.gene)]
assert len(present) == 42, len(present)

# ------------------------------------------------------------------ extra data
fib = pd.read_parquet(ROOT / "figures/data/sc/fibro_obs.parquet")
mac = pd.read_parquet(ROOT / "figures/data/sc/macro_obs.parquet")
expr = pd.read_parquet(ROOT / "figures/data/sc/cells_expr.parquet")
clin = pd.read_csv(ROOT / "figures/data/sc/sample_clinical_fibro.csv")
eff_ref = pd.read_csv(ROOT / "figures/data/sc/Fig2_effect_sizes.csv", index_col=0).iloc[:, 0]
e17b = pd.read_csv(S.RES / "E17b_bcell_reassignment.csv")


def cohen(x, y):
    x, y = np.asarray(x), np.asarray(y); nx, ny = len(x), len(y)
    sp = np.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2))
    return (x.mean() - y.mean()) / sp


def effects(o, lab, focal, comp, progs, key):
    """Cohen's d between patient pseudobulks (patients with >= 10 cells of the state) and between cells,
    exactly as figures/make_fig2.py (v14-v18 Fig. 2e,g); asserted against Fig2_effect_sizes.csv."""
    pb = o.groupby(["Patient_ID", lab], observed=True)[progs].mean()
    n = o.groupby(["Patient_ID", lab], observed=True).size(); pb = pb[n >= 10]
    A, B = pb.xs(focal, level=1), pb.xs(comp, level=1)
    dp = [cohen(A[p], B[p]) for p in progs]
    dc = [cohen(o.loc[o[lab] == focal, p], o.loc[o[lab] == comp, p]) for p in progs]
    ref = eval(eff_ref[key].replace("np.float64", "float"))
    assert all(abs(round(dp[i], 2) - ref[p][0]) < 1e-9 and abs(round(dc[i], 2) - ref[p][1]) < 1e-9
               for i, p in enumerate(progs)), key
    assert (len(A), len(B)) == eval(eff_ref[key + "_n"]), key
    return np.array(dp), np.array(dc), (len(A), len(B))


FPROG = ["prog_CAF_ECM", "prog_iCAF", "prog_STF3_FAP_INHBA", "prog_myCAF", "prog_CAF_TGFB"]
FNAME = ["ECM remodeling", "Inflammatory", "FAP⁺INHBA⁺", "Myofibroblastic", "TGF-β"]
MPROG = ["prog_SPP1_TAM", "prog_Inflam", "prog_HLAII_APC", "prog_IL10_STAT3", "prog_IFN"]
MNAME = ["SPP1⁺ TAM", "Inflammatory", "HLA-II APC", "IL-10/STAT3", "IFN response"]
fd = effects(fib, "fibro_subtype", "CAF_ECM", "iCAF", FPROG, "fibro")
md = effects(mac, "macro_state", "SPP1_TAM", "Inflam", MPROG, "macro")

# E17b: locate the reassigned cells (same rule as analysis/E17b_bcell_reassignment.py)
lab0 = obs.anno_level1.astype(str)
myeloid = (expr[["CD68", "LYZ", "C1QA", "FCER1G", "LST1"]] > 0).any(axis=1)
bmark = (expr[["MS4A1", "CD79A"]] > 0).any(axis=1)
move = ((lab0 == "Bcell") & myeloid & ~bmark).values
n_move = int(move.sum())
assert n_move == 12958 - int(e17b[(e17b.labels == "reassigned") & (e17b["class"] == "Bcell")].n_cells.iloc[0])
MARK = ["MS4A1", "CD79A", "CD68", "LYZ", "C1QA"]
mark_pct = {g: 100 * float((expr.loc[move, g] > 0).mean()) for g in MARK}

# ------------------------------------------------------------------ figure
FIG_W, FIG_H = 180, 230
fig = plt.figure(figsize=(FIG_W * S.MM, FIG_H * S.MM))
FW, FH = fig.get_size_inches()
x_um, y_um = obs.umap1.values, obs.umap2.values
xr = (x_um.min() - 0.6, x_um.max() + 0.6)
yr = (y_um.min() - 0.6, y_um.max() + 0.6)
aspect = (yr[1] - yr[0]) / (xr[1] - xr[0])
halo = [pe.withStroke(linewidth=1.6, foreground="white")]
COMP = "#6A4C9C"           # comparator state (iCAF, inflammatory macrophages); passes vs STROMA and MACRO
OTHER = "#D3D9DF"


def ftext(x_mm, y_mm, s, **kw):
    kw.setdefault("color", S.INK)
    return fig.text(x_mm * S.MM / FW, 1 - y_mm * S.MM / FH, s, **kw)


def head(xl, xt, yb, letter, title):
    """Panel letter and bold title on a common baseline (yb, mm from the top)."""
    ftext(xl, yb, letter, fontsize=8, fontweight="bold", va="baseline", ha="left")
    if title:
        ftext(xt, yb, title, fontsize=6.5, fontweight="bold", va="baseline", ha="left")


def chip(x_mm, y_mm, col, w=1.8, h=1.8):
    a_ = S.ax_mm(fig, x_mm, y_mm, w, h)
    a_.add_patch(Rectangle((0, 0), 1, 1, color=col, lw=0)); S.clean(a_)
    a_.set_xlim(0, 1); a_.set_ylim(0, 1)


# row baselines (mm from top)
Y1, Y2, Y3, Y4 = 5.0, 77.5, 125.5, 173.5

# a ---------------------------------------------------------------- UMAP by class + class key
head(1, 7, Y1, "a", "Discovery atlas (GSE183904)")
A_X, A_Y, A_W = 6, 8.0, 60
A_H = A_W * aspect
ax = S.ax_mm(fig, A_X, A_Y, A_W, A_H)
rng = np.random.default_rng(1)
order = rng.permutation(len(obs))
ax.scatter(x_um[order], y_um[order], s=0.3, lw=0, c=[CCOL[c] for c in obs["class"].values[order]],
           rasterized=True)
ax.set_xlim(*xr); ax.set_ylim(*yr); S.clean(ax)
LABPOS = {"plasma cell": ((12.7, 12.2), None), "epithelial/tumour": ((0.6, 8.6), None),
          "dendritic cell": ((12.7, 8.3), None), "B cell": ((17.6, 3.5), None),
          "NK cell": ((15.6, -3.0), None), "T cell": ((9.6, 0.9), None),
          "fibroblast": ((-4.2, 6.9), None), "pericyte/smooth muscle": ((4.6, 2.1), (-0.8, 2.6)),
          "endothelial": ((-2.0, -2.3), None), "macrophage/monocyte": ((3.6, -4.9), None)}
for c, ((tx, ty), tip) in LABPOS.items():
    if tip:
        ax.plot([tip[0], tx - 2.6], [tip[1], ty], color=S.INK, lw=0.4, zorder=5)
    ax.text(tx, ty, LAB[c], fontsize=5.4, ha="center", va="center", color=S.INK, path_effects=halo, zorder=6)
ax.text(0.0, 1.0, f"{len(obs):,} cells", transform=ax.transAxes,
        fontsize=5.2, color=S.SLATE, ha="left", va="top")
# cell-class key (shared by a, f, g and h)
KY = A_Y + A_H + 2.2
for k, c in enumerate(CLASSES):
    r, q = divmod(k, 4)
    x0, y0 = A_X + q * 16.2, KY + r * 3.0
    chip(x0, y0 + 0.2, CCOL[c], 1.8, 1.8)
    ftext(x0 + 2.5, y0 + 1.1, LAB[c].replace("Macrophage/monocyte", "Macro/monocyte"), fontsize=5.0,
          color=S.SLATE, va="center", ha="left")

# b ---------------------------------------------------------------- module scores on the UMAP
B_X0 = 78
head(72, B_X0, Y1, "b", "Per-cell module score")
bw = 31.5
bh = bw * aspect
gx, gy = 2.0, 4.0
fam = {"CAF_ECM": S.STROMA, "TAM_SUP": S.MACRO}
cms = {}
for i, m in enumerate(MODS):
    r, cidx = divmod(i, 3)
    x0, y0 = B_X0 + cidx * (bw + gx), A_Y + r * (bh + gy)
    axb = S.ax_mm(fig, x0, y0, bw, bh)
    v = sc[m].values
    lo, hi = np.nanpercentile(v, [1, 99])
    cm = LinearSegmentedColormap.from_list(m, ["#E4E7EB", fam.get(m, S.IMMUNE)])
    cms[m] = cm
    o_ = np.argsort(v, kind="stable")
    axb.scatter(x_um[o_], y_um[o_], s=0.1, lw=0, c=np.clip(v[o_], lo, hi), cmap=cm, vmin=lo, vmax=hi,
                rasterized=True)
    axb.set_xlim(*xr); axb.set_ylim(*yr); S.clean(axb)
    chip(x0, y0 + 0.5, fam.get(m, S.IMMUNE), 1.6, 1.6)
    ftext(x0 + 2.3, y0 + 1.3, m, fontsize=5.8, va="center", ha="left")
kx = B_X0 + 3 * bw + 2 * gx - 23.5
for k, m in enumerate(["CAF_ECM", "TAM_SUP", "APM"]):
    kax = S.ax_mm(fig, kx + 6 + k * 5.0, Y1 - 1.8, 4.7, 1.6)
    kax.imshow(np.linspace(0, 1, 64)[None, :], aspect="auto", cmap=cms[m]); S.clean(kax)
ftext(kx + 5.4, Y1 - 1.0, "low", fontsize=5.0, color=S.SLATE, ha="right", va="center")
ftext(kx + 6 + 15.2, Y1 - 1.0, "high", fontsize=5.0, color=S.SLATE, ha="left", va="center")
B_BOTTOM = A_Y + 2 * bh + gy

# c ---------------------------------------------------------------- dot plot
head(1, 7, Y2, "c", "Module genes by cell class")
C_X, C_W = 30, 133
C_Y = Y2 + 7.0
cls_disc = [c for c in CLASSES if c in set(gex["class"])]
nrow, ncol = len(cls_disc), len(present)
C_H = nrow * 2.65
axc = S.ax_mm(fig, C_X, C_Y, C_W, C_H)
G = gex.set_index(["gene", "class"])
mean = np.array([[G.loc[(g, c), "mean_expr"] for g in present] for c in cls_disc])
pct = np.array([[G.loc[(g, c), "pct_expressing"] for g in present] for c in cls_disc])
scaled = mean / mean.max(0, keepdims=True)
xx, yy = np.meshgrid(np.arange(ncol), np.arange(nrow))
smax = 11.0
axc.scatter(xx.ravel(), yy.ravel(), s=(pct.ravel() / 100) * smax, c=scaled.ravel(), cmap=S.SEQ, vmin=0, vmax=1,
            lw=0.25, edgecolors="#9AA5B0", zorder=3)
axc.set_xlim(-0.6, ncol - 0.4); axc.set_ylim(nrow - 0.5, -0.5)
axc.set_yticks(range(nrow)); axc.set_yticklabels([LAB[c] for c in cls_disc], fontsize=5.4)
axc.set_xticks(range(ncol)); axc.set_xticklabels(present, rotation=90, fontstyle="italic", fontsize=5.1)
axc.tick_params(length=0, pad=2.0)
for s_ in axc.spines.values():
    s_.set_visible(False)
start = 0
for m in MODS:
    n = sum(g in present for g in MODULES[m])
    axc.add_patch(Rectangle((start - 0.45, -1.15), n - 0.1, 0.36, color=S.MODCOL[m], clip_on=False, lw=0))
    axc.text(start + (n - 1) / 2, -1.33, m, ha="center", va="bottom", fontsize=5.5, fontweight="bold", color=S.INK)
    if start:
        axc.axvline(start - 0.5, color=S.LIGHT, lw=0.5, zorder=0)
    for c in INTENDED.get(m, []):
        if c in cls_disc:
            j = cls_disc.index(c)
            axc.add_patch(Rectangle((start - 0.5, j - 0.5), n, 1, color=S.PAPER, lw=0, zorder=0))
    start += n
lx = C_X + C_W + 4
lax = S.ax_mm(fig, lx, C_Y - 1.0, 13, 13)
S.clean(lax); lax.set_xlim(0, 1); lax.set_ylim(0, 1)
lax.text(0, 1.0, "Cells\nexpressing", fontsize=5.0, color=S.SLATE, va="top", linespacing=1.1)
for k, p_ in enumerate([25, 50, 100]):
    yv = 0.58 - k * 0.22
    lax.scatter([0.1], [yv], s=(p_ / 100) * smax, color="white", edgecolors="#9AA5B0", lw=0.25)
    lax.text(0.3, yv, f"{p_}%", fontsize=5.0, color=S.SLATE, va="center")
cbx = S.ax_mm(fig, lx + 0.4, C_Y + 17.5, 1.8, 7.5)
cbx.imshow(np.linspace(1, 0, 64)[:, None], aspect="auto", cmap=S.SEQ)
S.clean(cbx)
cbx.text(2.0, 1.0, "max", transform=cbx.transAxes, fontsize=5.0, color=S.SLATE, va="top")
cbx.text(2.0, 0.0, "0", transform=cbx.transAxes, fontsize=5.0, color=S.SLATE, va="bottom")
cbx.text(-0.2, 1.1, "Mean,\nscaled", transform=cbx.transAxes, fontsize=5.0, color=S.SLATE, va="bottom",
         linespacing=1.1)


# d, e -------------------------------------------------------------- fibroblast and macrophage states
def states(x0, letter, title, o, lab, focal, comp, fcol, names, eff, flab, clab):
    head(x0 - 5, x0, Y3, letter, title)
    u1, u2 = o.umap1.values, o.umap2.values
    xr_ = (u1.min() - 0.4, u1.max() + 0.4); yr_ = (u2.min() - 0.4, u2.max() + 0.4)
    asp = (yr_[1] - yr_[0]) / (xr_[1] - xr_[0])
    uw = 26.0
    uh = min(uw * asp, 27.0); uw = uh / asp
    ua = S.ax_mm(fig, x0, Y3 + 3.5, uw, uh)
    st = o[lab].astype(str).values
    col = np.where(st == focal, fcol, np.where(st == comp, COMP, OTHER))
    rk = np.where(st == focal, 2, np.where(st == comp, 1, 0))
    ordr = np.lexsort((np.random.default_rng(2).random(len(o)), rk))
    ua.scatter(u1[ordr], u2[ordr], s=0.5, lw=0, c=col[ordr], rasterized=True)
    ua.set_xlim(*xr_); ua.set_ylim(*yr_); S.clean(ua)
    ua.text(0.0, 1.0, f"{len(o):,} cells", transform=ua.transAxes, fontsize=5.0, color=S.SLATE, va="top")
    # key under the UMAP: chips + ink labels
    ky = Y3 + 3.5 + uh + 1.2
    for k, (cc, tt) in enumerate([(fcol, f"{flab} ({(st == focal).sum():,})"),
                                   (COMP, f"{clab} ({(st == comp).sum():,})"), (OTHER, "Other states")]):
        chip(x0, ky + k * 2.6, cc, 1.6, 1.6)
        ftext(x0 + 2.3, ky + k * 2.6 + 0.8, tt, fontsize=5.0, color=S.SLATE, va="center")
    # effect sizes
    dp, dc, (na, nb) = eff
    ex = x0 + uw + 17.0
    ew = 14.5
    ea = S.ax_mm(fig, ex, Y3 + 4.5, ew, 24.0)
    yy_ = np.arange(len(names))
    ea.barh(yy_, dp, height=0.62, color=[fcol if v > 0 else COMP for v in dp], lw=0)
    ea.plot(dc, yy_, "o", mfc="white", mec=S.INK, ms=2.4, mew=0.5, zorder=4)
    ea.axvline(0, color=S.INK, lw=0.5)
    ea.set_ylim(len(names) - 0.5, -0.5)
    ea.set_yticks(yy_); ea.set_yticklabels(names, fontsize=5.1)
    ea.tick_params(axis="y", length=0, pad=1.5)
    ea.spines["left"].set_visible(False)
    lim = np.ceil(np.abs(np.r_[dp, dc]).max())
    ea.set_xlim(-lim - 0.3, lim + 0.3)
    ea.set_xticks([-lim, 0, lim]); ea.set_xticklabels([f"{-lim:.0f}".replace("-", "−"), "0", f"{lim:.0f}"],
                                                      fontsize=5.0)
    ea.set_xlabel("Cohen's d", fontsize=5.4, labelpad=1.5)
    ea.set_title(f"{flab} vs {clab}", loc="right", fontsize=5.0, color=S.SLATE, fontweight="normal", pad=2.0)



states(7, "d", "Fibroblast states", fib, "fibro_subtype", "CAF_ECM", "iCAF", S.STROMA, FNAME, fd,
       "ECM-remodeling CAF", "iCAF")
states(66, "e", "Macrophage states", mac, "macro_state", "SPP1_TAM", "Inflam", S.MACRO, MNAME, md,
       "SPP1⁺ TAM", "Inflammatory")

# f ---------------------------------------------------------------- per-sample composition
F_X, F_W = 132, 47
head(F_X - 5, F_X, Y3, "f", "Composition per sample")
fr = obs.groupby(["sample_id", "class"], observed=True).size().unstack(fill_value=0)
fr = fr.div(fr.sum(1), axis=0)
cl = clin.set_index("sample_id")
sub = cl.TCGA.fillna("n.a.").replace({"N.A": "n.a."}).reindex(fr.index)
tissue = obs.groupby("sample_id", observed=True).tissue.first().reindex(fr.index)
GROUPS = ["CIN", "GS", "MSI", "n.a."]
ordr = []
for g in GROUPS:
    ids = fr.index[sub.values == g]
    ordr += list(fr.loc[ids, "fibroblast"].sort_values(ascending=False).index)
fa = S.ax_mm(fig, F_X + 6, Y3 + 4.5, F_W - 6, 24.0)
xs, x = [], 0.0
gpos = {}
for g in GROUPS:
    ids = [i for i in ordr if sub[i] == g]
    gpos[g] = (x, x + len(ids) - 1)
    for i in ids:
        xs.append(x); x += 1
    x += 0.8
for xi, sid in zip(xs, ordr):
    bottom = 0.0
    for c in [c for c in CLASSES if c in fr.columns]:
        v = fr.loc[sid, c]
        fa.bar(xi, v, bottom=bottom, width=0.86, color=CCOL[c], lw=0)
        bottom += v
    if tissue[sid] == "Peritoneum":
        fa.plot(xi, 1.045, marker="v", ms=2.2, color=S.INK, mec="none", clip_on=False)
fa.set_xlim(-0.7, xs[-1] + 0.7); fa.set_ylim(0, 1)
fa.set_xticks([])
fa.set_yticks([0, 0.5, 1]); fa.set_yticklabels(["0", "0.5", "1"])
fa.set_ylabel("Fraction of cells", labelpad=1.5)
fa.spines["bottom"].set_visible(False)
for g in GROUPS:
    a_, b_ = gpos[g]
    fa.plot([a_ - 0.4, b_ + 0.4], [-0.04, -0.04], color=S.INK, lw=0.5, clip_on=False)
    fa.text((a_ + b_) / 2, -0.07, g, fontsize=5.0, ha="center", va="top", color=S.INK)
ftext(F_X + 6, Y3 + 33.4, "TCGA subtype", fontsize=5.0, color=S.SLATE, va="top")
pk = S.ax_mm(fig, F_X + 22.0, Y3 + 33.5, 1.8, 1.8); S.clean(pk); pk.set_xlim(0, 1); pk.set_ylim(0, 1)
pk.plot(0.5, 0.5, marker="v", ms=2.2, color=S.INK, mec="none")
ftext(F_X + 24.3, Y3 + 33.4, "peritoneal sample", fontsize=5.0, color=S.SLATE, va="top")

# g / h ------------------------------------------------------------- 18 aligned rows: module x atlas
head(1, 7, Y4, "g", "Signal by cell class in three atlases")
R_Y = Y4 + 9.5
barh, gap_in, gap_out = 1.8, 0.4, 1.75
rows = []
y = 0.0
for m in MODS:
    for ds in DATASETS:
        rows.append((m, ds, y)); y += barh + gap_in
    y += gap_out - gap_in
R_H = y - gap_out + gap_in
D_X, D_W = 34, 46
axd = S.ax_mm(fig, D_X, R_Y, D_W, R_H)
for m, ds, y0 in rows:
    d = attr[(attr.dataset == ds) & (attr.module == m)].set_index("class").share_pct
    first = [c for c in INTENDED.get(m, []) if c in d.index]
    rest = [c for c in CLASSES if c in d.index and c not in first]
    left = 0.0
    for c in first + rest:
        axd.add_patch(Rectangle((left, y0), d[c], barh, facecolor=CCOL[c], edgecolor="white", lw=0.25))
        left += d[c]
    if first:
        s_int = d[first].sum()
        axd.add_patch(Rectangle((0, y0), s_int, barh, facecolor="none", edgecolor=S.INK, lw=0.55, zorder=4))
        axd.text(105.5, y0 + barh / 2, f"{s_int:.0f}%", fontsize=5.0, va="center", ha="left", color=S.INK)
    axd.text(-1.5, y0 + barh / 2, DSLAB[ds], fontsize=5.0, va="center", ha="right", color=S.SLATE)
for m in MODS:
    ys = [y0 for mm, ds, y0 in rows if mm == m]
    yc = R_Y + (ys[0] + ys[-1] + barh) / 2
    chip(2.0, R_Y + ys[0] + 0.2, S.MODCOL[m], 1.3, ys[-1] + barh - ys[0] - 0.4)
    ftext(4.3, yc, m, fontsize=5.5, fontweight="bold", va="center", ha="left")
axd.axvline(50, color=S.INK, lw=0.5, ls=(0, (2, 2)), zorder=6)
axd.set_xlim(0, 100); axd.set_ylim(R_H, 0)
axd.set_yticks([]); axd.spines["left"].set_visible(False)
axd.spines["bottom"].set_visible(False); axd.spines["top"].set_visible(True)
axd.xaxis.tick_top(); axd.xaxis.set_label_position("top")
axd.set_xticks([0, 50, 100]); axd.set_xticklabels(["0", "50", "100"], fontsize=5.0)
axd.set_xlabel("Share of module signal (%)", fontsize=5.4, labelpad=1.5)
axd.text(105.5, -1.2, "Intended\nclass", fontsize=5.0, color=S.SLATE, ha="left", va="bottom", linespacing=1.1)

E_X = 97
cw = 2.55
head(E_X - 6, E_X, Y4, "h", "Enrichment (z)")
axe = S.ax_mm(fig, E_X, R_Y, len(CLASSES) * cw, R_H)
norm = Normalize(-1, 3)
for m, ds, y0 in rows:
    Z = attr[(attr.dataset == ds) & (attr.module == m)].set_index("class").z
    for j, c in enumerate(CLASSES):
        if c not in Z.index:
            axe.add_patch(Rectangle((j, y0), 1, barh, facecolor="white", edgecolor=S.LIGHT, lw=0.3, hatch="//////"))
        else:
            axe.add_patch(Rectangle((j, y0), 1, barh, facecolor=S.SEQ(norm(Z[c])), edgecolor="white", lw=0.3))
    for c in INTENDED.get(m, []):
        j = CLASSES.index(c)
        axe.add_patch(Rectangle((j + 0.06, y0 + 0.08), 0.88, barh - 0.16, facecolor="none", edgecolor=S.INK,
                                lw=0.55, zorder=4))
    jt = CLASSES.index(Z.idxmax())
    axe.plot(jt + 0.5, y0 + barh / 2, marker="o", ms=1.2, color="white", mec="none", zorder=5)
axe.set_xlim(0, len(CLASSES)); axe.set_ylim(R_H, 0)
axe.set_yticks([])
axe.xaxis.tick_top()
axe.set_xticks(np.arange(len(CLASSES)) + 0.5)
HSHORT = {"plasma cell": "Pl", "macrophage/monocyte": "Mac", "fibroblast": "Fib", "pericyte/smooth muscle": "Peri"}
axe.set_xticklabels([HSHORT.get(c, SHORT[c]) for c in CLASSES], rotation=90, fontsize=5.0)
axe.tick_params(length=0, pad=3.6)
for s_ in axe.spines.values():
    s_.set_visible(False)
for j, c in enumerate(CLASSES):
    axe.add_patch(Rectangle((j + 0.15, -0.95), 0.7, 0.7, color=CCOL[c], clip_on=False, lw=0))
cbe = S.ax_mm(fig, E_X + len(CLASSES) * cw + 2.0, R_Y + 1, 1.6, 12)
cbe.imshow(np.linspace(3, -1, 64)[:, None], aspect="auto", cmap=S.SEQ, norm=norm)
S.clean(cbe)
for v_, t_ in [(3, "≥3"), (0, "0"), (-1, "−1")]:
    cbe.text(2.0, (v_ + 1) / 4, t_, transform=cbe.transAxes, fontsize=5.0, color=S.SLATE, va="center")
cbe.text(0.0, 1.06, "z", transform=cbe.transAxes, fontsize=5.0, color=S.SLATE, va="bottom", fontstyle="italic")

# i ---------------------------------------------------------------- B-cell annotation check (E17b)
I_X = 137
head(I_X - 5, I_X, Y4, "i", "B-cell label check")
iu = S.ax_mm(fig, I_X, Y4 + 2.5, 22, 22 * aspect)
bl = (lab0 == "Bcell").values & ~move
iu.scatter(x_um, y_um, s=0.06, lw=0, c="#E4E7EB", rasterized=True)
iu.scatter(x_um[bl], y_um[bl], s=0.12, lw=0, c=CCOL["B cell"], rasterized=True)
iu.scatter(x_um[move], y_um[move], s=0.12, lw=0, c=S.INK, rasterized=True)
iu.set_xlim(*xr); iu.set_ylim(*yr); S.clean(iu)
ky = Y4 + 2.5 + 22 * aspect + 1.0
for k, (cc, tt) in enumerate([(S.INK, f"Reassigned ({n_move:,})"), (CCOL["B cell"], f"B cell, kept ({int(bl.sum()):,})")]):
    chip(I_X, ky + k * 2.6, cc, 1.6, 1.6)
    ftext(I_X + 2.3, ky + k * 2.6 + 0.8, tt, fontsize=5.0, color=S.SLATE, va="center")
im_ = S.ax_mm(fig, I_X + 30.5, Y4 + 3.0, 10.5, 14.0)
yy_ = np.arange(len(MARK))
im_.barh(yy_, [mark_pct[g] for g in MARK], height=0.62,
         color=[CCOL["B cell"] if g in ("MS4A1", "CD79A") else CCOL["macrophage/monocyte"] for g in MARK], lw=0)
for k, g in enumerate(MARK):
    im_.text(mark_pct[g] + 3, k, f"{mark_pct[g]:.0f}", fontsize=5.0, va="center", color=S.INK)
im_.set_ylim(len(MARK) - 0.5, -0.5); im_.set_xlim(0, 100)
im_.set_yticks(yy_); im_.set_yticklabels(MARK, fontstyle="italic", fontsize=5.0)
im_.tick_params(axis="y", length=0, pad=1.2)
im_.spines["left"].set_visible(False)
im_.set_xticks([0, 100]); im_.set_xticklabels(["0", "100"], fontsize=5.0)
im_.set_xlabel("% detected", fontsize=5.0, labelpad=1.0)
# share of the intended class, as annotated -> reassigned
db = S.ax_mm(fig, I_X + 14, Y4 + 30.0, 18.5, 13.5)
DMOD = [("CAF_ECM", ["Fibro"]), ("TAM_SUP", ["Macro"]), ("TLS_B", ["Bcell", "Plasma"]), ("T_CELL", ["Tcell", "NK"])]
shares = {}
for k, (m, cls) in enumerate(DMOD):
    d = e17b[(e17b.module == m) & e17b["class"].isin(cls)].groupby("labels").share_pct.sum()
    a0, a1 = d["as annotated"], d["reassigned"]
    shares[m] = (a0, a1)
    db.plot([a0, a1], [k, k], color=S.GREY, lw=0.8, zorder=1)
    db.plot(a1, k, "o", color=S.MODCOL[m], ms=2.4, mec="none", zorder=3)          # reassigned, filled
    db.plot(a0, k, "o", mfc="none", mec=S.INK, ms=3.4, mew=0.55, zorder=4)         # as annotated, ring on top
    db.text(104, k, f"{a0:.0f} / {a1:.0f}", fontsize=5.0, va="center", ha="left", color=S.INK)
db.set_ylim(len(DMOD) - 0.5, -0.5); db.set_xlim(0, 100)
db.set_yticks(range(len(DMOD))); db.set_yticklabels([m for m, _ in DMOD], fontsize=5.0)
db.tick_params(axis="y", length=0, pad=1.5)
db.spines["left"].set_visible(False)
db.set_xticks([0, 50, 100]); db.set_xticklabels(["0", "50", "100"], fontsize=5.0)
db.set_xlabel("Intended-class share (%)", fontsize=5.0, labelpad=1.0)
db.axvline(50, color=S.GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
for k, (mfc, tt) in enumerate([("none", "as annotated"), (S.GREY, "filled, reassigned")]):
    mk = S.ax_mm(fig, I_X + 6 + k * 15.5, Y4 + 50.2, 1.8, 1.8); S.clean(mk); mk.set_xlim(0, 1); mk.set_ylim(0, 1)
    if mfc == "none":
        mk.plot(0.5, 0.5, "o", mfc="none", mec=S.INK, ms=3.4, mew=0.55)
    else:
        mk.plot(0.5, 0.5, "o", color=mfc, mec="none", ms=2.4)
    ftext(I_X + 8.2 + k * 15.5, Y4 + 51.1, tt, fontsize=5.0, color=S.SLATE, va="center")

print("bottom of g/h rows (mm):", R_Y + R_H)
print("reassigned:", n_move, {g: round(v, 1) for g, v in mark_pct.items()}, {m: tuple(round(x, 1) for x in v) for m, v in shares.items()})
print("effect sizes fibro:", np.round(fd[0], 2), np.round(fd[1], 2), fd[2], "macro:", np.round(md[0], 2), np.round(md[1], 2), md[2])
print("missing genes in discovery atlas:", missing)
print(S.save(fig, "Fig2"))
