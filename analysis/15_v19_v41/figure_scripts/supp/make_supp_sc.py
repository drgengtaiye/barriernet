"""v19 supplementary single-cell figures in the cns house style.

    FigS_sc_atlas         (a-b) canonical markers, program enrichment (annotation of the discovery atlas)
    FigS_sc_programmes    (a-g) program scores by cell type (box plots)
    FigS_sc_composition   (a-c) mean composition by AJCC stage, Lauren type, TCGA subtype
    FigS_sc_fractions     (a-b) fibroblast vs macrophage state fractions across patients
    FigS_sc_indep         (a-b) per-patient intended-class shares in three atlases; OMIX001073 module genes

Panel content and letters follow the legends in manuscript/v19_supp_src.md. Every plotted value comes
from an existing table or is a descriptive summary of figures/data/sc; the two statistics printed in
the images (Kruskal-Wallis P in sc_composition, Spearman rho/P in sc_fractions) re-run the original
code of figures/supp/make_supp_sc.py on the same inputs, and the script asserts that they reproduce
the values in the previous figures. Usage: python make_supp_sc.py [key ...]
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
from scipy.stats import kruskal, spearmanr

ROOT = S.ROOT
D = ROOT / "figures/data/sc"
GC = ROOT.parent / "GC_TME_Reprogrammer"
OUTDIR = S.OUT / "supp"
OUTDIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------ shared encodings (as main Fig. 2)
MODULES = {
    "CAF_ECM": ["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
    "TAM_SUP": ["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
    "APM": ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    "IFN": ["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
    "TLS_B": ["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
    "T_CELL": ["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"],
}
MODS = list(MODULES)
GENES48 = [g for m in MODS for g in MODULES[m]]
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
CCOL["mast cell"] = "#9A5B13"
ANNO2CLASS = {v: k for k, v in KEY.items()}
DISC = [c for c in CLASSES if c != "mast cell"]          # the discovery atlas has no mast class
INTENDED = {"CAF_ECM": ["fibroblast"], "TAM_SUP": ["macrophage/monocyte"],
            "TLS_B": ["B cell", "plasma cell"], "T_CELL": ["T cell", "NK cell"]}
DATASETS = ["GSE183904", "OMIX001073", "GSE206785"]
DSLAB = {"GSE183904": "GSE183904 (discovery)", "OMIX001073": "OMIX001073", "GSE206785": "GSE206785"}
FAM = {m: S.MODCOL[m] for m in MODS}
COMP = "#6A4C9C"          # comparator state (iCAF, inflammatory macrophages)
OTHER = "#D3D9DF"
HALO = [pe.withStroke(linewidth=1.6, foreground="white")]
ZNORM = Normalize(-1, 3)

attr = pd.read_csv(S.RES / "E17_attribution_by_class.csv")
attr = attr[attr.analysis == "primary"]
gexall = pd.read_csv(S.RES / "E17_gene_expression_by_class.csv")
REPORT = {}


# ------------------------------------------------------------------ helpers (mm layout, top-left origin)
class Page:
    def __init__(self, w_mm, h_mm):
        self.fig = plt.figure(figsize=(w_mm * S.MM, h_mm * S.MM))
        self.W, self.H = self.fig.get_size_inches()

    def ax(self, x, y, w, h):
        return S.ax_mm(self.fig, x, y, w, h)

    def text(self, x, y, s, **kw):
        kw.setdefault("color", S.INK)
        return self.fig.text(x * S.MM / self.W, 1 - y * S.MM / self.H, s, **kw)

    def head(self, xl, xt, yb, letter, title=None, sub=None):
        """Panel letter and bold title on one baseline (yb mm from the top); subtitles go to the legend."""
        sub = None
        if letter:
            self.text(xl, yb, letter, fontsize=8, fontweight="bold", va="baseline", ha="left")
        if title:
            self.text(xt, yb, title, fontsize=6.5, fontweight="bold", va="baseline", ha="left")
        if sub:
            self.text(xt, yb + 3.0, sub, fontsize=5.2, color=S.SLATE, va="baseline", ha="left")

    def chip(self, x, y, col, w=1.8, h=1.8):
        a = self.ax(x, y, w, h)
        a.add_patch(Rectangle((0, 0), 1, 1, color=col, lw=0)); S.clean(a)
        a.set_xlim(0, 1); a.set_ylim(0, 1)

    def key(self, x, y, items, ncol=1, dx=16.0, dy=2.9, size=5.0):
        for k, (col, lab) in enumerate(items):
            r, q = divmod(k, ncol)
            self.chip(x + q * dx, y + r * dy, col, 1.8, 1.8)
            self.text(x + q * dx + 2.5, y + r * dy + 0.9, lab, fontsize=size, color=S.SLATE, va="center")

    def save(self, stem):
        self.fig.savefig(OUTDIR / f"{stem}.pdf")
        self.fig.savefig(OUTDIR / f"{stem}.png", dpi=600)
        plt.close(self.fig)
        print("wrote", OUTDIR / f"{stem}.png")


def heat(ax, Z, rows, cols, rowlab, collab, boxes=(), numbers=True, rot=45, fs=5.0):
    """Cells coloured on the shared z scale (S.SEQ, -1..3), values printed, NaN hatched."""
    for i, r in enumerate(rows):
        for j, c in enumerate(cols):
            v = Z.loc[r, c] if (r in Z.index and c in Z.columns) else np.nan
            if np.isnan(v):
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor="white", edgecolor=S.LIGHT, lw=0.3, hatch="//////"))
                continue
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=S.SEQ(ZNORM(v)), edgecolor="white", lw=0.4))
            if numbers:
                ax.text(j + 0.5, i + 0.52, f"{v:.1f}".replace("-", "−").replace("−0.0", "0.0"), ha="center",
                        va="center", fontsize=fs, color="white" if v > 1.6 else S.INK)
    for (i, j) in boxes:
        ax.add_patch(Rectangle((j + 0.04, i + 0.05), 0.92, 0.9, facecolor="none", edgecolor=S.INK, lw=0.7, zorder=4))
    ax.set_xlim(0, len(cols)); ax.set_ylim(len(rows), 0)
    ax.set_yticks(np.arange(len(rows)) + 0.5); ax.set_yticklabels(rowlab, fontsize=5.2)
    ax.set_xticks(np.arange(len(cols)) + 0.5)
    ax.set_xticklabels(collab, rotation=rot, ha="right" if rot not in (0, 90) else "center",
                       rotation_mode="anchor" if rot not in (0, 90) else "default", fontsize=5.2)
    ax.tick_params(length=0, pad=1.5)
    for s_ in ax.spines.values():
        s_.set_visible(False)


def zbar(pg, x, y, h=12, label="z across classes"):
    cb = pg.ax(x, y, 1.6, h)
    cb.imshow(np.linspace(3, -1, 64)[:, None], aspect="auto", cmap=S.SEQ, norm=ZNORM)
    S.clean(cb)
    for v_, t_ in [(3, "≥3"), (1, "1"), (-1, "≤−1")]:
        cb.text(2.0, (v_ + 1) / 4, t_, transform=cb.transAxes, fontsize=5.0, color=S.SLATE, va="center")
    pg.text(x, y - 1.2, label, fontsize=5.0, color=S.SLATE, va="baseline")


def dotplot(ax, mean, pct, rows, cols, smax=11.0, edge="#9AA5B0"):
    scaled = mean / np.nanmax(mean, axis=0, keepdims=True)
    xx, yy = np.meshgrid(np.arange(len(cols)), np.arange(len(rows)))
    ok = ~np.isnan(pct.ravel())
    ax.scatter(xx.ravel()[ok], yy.ravel()[ok], s=(pct.ravel()[ok] / 100) * smax, c=scaled.ravel()[ok], cmap=S.SEQ,
               vmin=0, vmax=1, lw=0.25, edgecolors=edge, zorder=3)
    ax.set_xlim(-0.6, len(cols) - 0.4); ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=90, fontstyle="italic", fontsize=5.0)
    ax.tick_params(length=0, pad=2.0)
    for s_ in ax.spines.values():
        s_.set_visible(False)


def dotkeys(pg, x, y, smax=11.0, sizes=(25, 50, 100)):
    pg.text(x, y, "Cells\nexpressing", fontsize=5.0, color=S.SLATE, va="top", linespacing=1.1)
    for k, p_ in enumerate(sizes):
        a = pg.ax(x, y + 5.3 + k * 2.6, 2.5, 2.5); S.clean(a); a.set_xlim(0, 1); a.set_ylim(0, 1)
        a.scatter([0.5], [0.5], s=(p_ / 100) * smax, color="white", edgecolors="#9AA5B0", lw=0.25, clip_on=False)
        pg.text(x + 3.2, y + 6.6 + k * 2.6, f"{p_}%", fontsize=5.0, color=S.SLATE, va="center")
    cb = pg.ax(x + 0.4, y + 18.5, 1.6, 7.5)
    cb.imshow(np.linspace(1, 0, 64)[:, None], aspect="auto", cmap=S.SEQ); S.clean(cb)
    cb.text(2.0, 1.0, "max", transform=cb.transAxes, fontsize=5.0, color=S.SLATE, va="top")
    cb.text(2.0, 0.0, "0", transform=cb.transAxes, fontsize=5.0, color=S.SLATE, va="bottom")
    pg.text(x, y + 17.4, "Mean, scaled", fontsize=5.0, color=S.SLATE, va="baseline")


def module_bands(ax, genes, y=-1.15, h=0.36, label_y=-1.33):
    start = 0
    for m in MODS:
        n = sum(g in genes for g in MODULES[m])
        if n == 0:
            continue
        ax.add_patch(Rectangle((start - 0.45, y), n - 0.1, h, color=FAM[m], clip_on=False, lw=0))
        ax.text(start + (n - 1) / 2, label_y, m, ha="center", va="bottom", fontsize=5.5, fontweight="bold",
                color=S.INK)
        if start:
            ax.axvline(start - 0.5, color=S.LIGHT, lw=0.5, zorder=0)
        start += n


def cohen(x, y):
    x, y = np.asarray(x), np.asarray(y); nx, ny = len(x), len(y)
    sp = np.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2))
    return (x.mean() - y.mean()) / sp


def effects(o, lab, focal, comp, progs, key):
    """Cohen's d, patient pseudobulks (>= 10 cells of the state) and cells; asserted against
    figures/data/sc/Fig2_effect_sizes.csv (written by the original figures/make_fig2.py)."""
    ref_all = pd.read_csv(D / "Fig2_effect_sizes.csv", index_col=0).iloc[:, 0]
    pb = o.groupby(["Patient_ID", lab], observed=True)[progs].mean()
    n = o.groupby(["Patient_ID", lab], observed=True).size(); pb = pb[n >= 10]
    A, B = pb.xs(focal, level=1), pb.xs(comp, level=1)
    dp = np.array([cohen(A[p], B[p]) for p in progs])
    dc = np.array([cohen(o.loc[o[lab] == focal, p], o.loc[o[lab] == comp, p]) for p in progs])
    ref = eval(ref_all[key].replace("np.float64", "float"))
    assert all(abs(round(dp[i], 2) - ref[p][0]) < 1e-9 and abs(round(dc[i], 2) - ref[p][1]) < 1e-9
               for i, p in enumerate(progs)), key
    assert (len(A), len(B)) == eval(ref_all[key + "_n"])
    return dp, dc, (len(A), len(B))


def state_umap(pg, ax, o, lab, focal, comp, fcol, flab, clab, x, ykey):
    u1, u2 = o.umap1.values, o.umap2.values
    st = o[lab].astype(str).values
    col = np.where(st == focal, fcol, np.where(st == comp, COMP, OTHER))
    rk = np.where(st == focal, 2, np.where(st == comp, 1, 0))
    ordr = np.lexsort((np.random.default_rng(2).random(len(o)), rk))
    ax.scatter(u1[ordr], u2[ordr], s=0.6, lw=0, c=col[ordr], rasterized=True)
    ax.set_xlim(u1.min() - 0.4, u1.max() + 0.4); ax.set_ylim(u2.min() - 0.4, u2.max() + 0.4); S.clean(ax)
    pg.key(x, ykey, [(fcol, f"{flab} ({(st == focal).sum():,})"), (COMP, f"{clab} ({(st == comp).sum():,})"),
                     (OTHER, "Other states")], ncol=1, dy=2.7)


def effect_bars(pg, ea, names, eff, fcol, focal_lab, comp_lab, xnote, ynote):
    dp, dc, (na, nb) = eff
    yy = np.arange(len(names))
    ea.barh(yy, dp, height=0.62, color=[fcol if v > 0 else COMP for v in dp], lw=0)
    ea.plot(dc, yy, "o", mfc="white", mec=S.INK, ms=2.4, mew=0.5, zorder=4)
    ea.axvline(0, color=S.INK, lw=0.5)
    ea.set_ylim(len(names) - 0.5, -0.5)
    ea.set_yticks(yy); ea.set_yticklabels(names, fontsize=5.2)
    ea.tick_params(axis="y", length=0, pad=1.5)
    ea.spines["left"].set_visible(False)
    lim = np.ceil(np.abs(np.r_[dp, dc]).max())
    ea.set_xlim(-lim - 0.3, lim + 0.3)
    step = 1 if lim <= 3 else (2 if lim <= 6 else 3)
    k_ = np.floor(lim / step)
    ticks = np.arange(-k_ * step, k_ * step + 0.1, step)
    ea.set_xticks(ticks); ea.set_xticklabels([f"{t:.0f}".replace("-", "−") for t in ticks], fontsize=5.0)
    ea.set_xlabel(f"Cohen's d, {focal_lab} vs {comp_lab}", fontsize=5.4, labelpad=1.5)
    pg.text(xnote, ynote, f"Bars, patient pseudobulks ({na} vs {nb})", fontsize=5.0, color=S.SLATE, va="baseline")
    pg.text(xnote, ynote + 2.5, "Circles, cells", fontsize=5.0, color=S.SLATE, va="baseline")


# ================================================================== sc_atlas
def fig_sc_atlas():
    """Annotation of the discovery atlas: (a) canonical markers, (b) program enrichment.
    The atlas UMAP, fibroblast/macrophage states, module attribution and module genes are in main Fig. 2."""
    obs = pd.read_parquet(D / "cells_obs.parquet")
    obs["class"] = obs.anno_level1.map(ANNO2CLASS)
    E = pd.read_parquet(D / "cells_expr.parquet")
    pg = Page(180, 62)
    Y1 = 5.0
    # a: canonical markers
    CANON = [("T cell", ["CD3D", "CD3E"]), ("NK cell", ["NKG7", "GNLY"]), ("B cell", ["MS4A1", "CD79A"]),
             ("plasma cell", ["MZB1", "JCHAIN"]), ("macrophage/monocyte", ["CD68", "C1QA", "LYZ"]),
             ("dendritic cell", ["LST1", "CLEC10A"]), ("fibroblast", ["COL1A1"]), ("endothelial", ["PECAM1", "VWF"]),
             ("pericyte/smooth muscle", ["RGS5", "MCAM"]), ("epithelial/tumour", ["EPCAM", "KRT19", "KRT8"])]
    genes = [g for _, gs in CANON for g in gs]
    cl = obs["class"].values
    mean = np.array([[E.loc[cl == c, g].mean() for g in genes] for c in DISC])
    pct = np.array([[(E.loc[cl == c, g] > 0).mean() * 100 for g in genes] for c in DISC])
    pg.head(1, 7, Y1, "a", "Canonical markers", f"GSE183904, {len(obs):,} cells, 10 types")
    bx = pg.ax(27, Y1 + 12.0, 70, len(DISC) * 3.1)
    dotplot(bx, mean, pct, DISC, genes, smax=15.0)
    bx.set_yticks(range(len(DISC))); bx.set_yticklabels([LAB[c] for c in DISC], fontsize=5.2)
    start = 0
    for c, gs in CANON:
        n = len(gs)
        bx.add_patch(Rectangle((start - 0.42, -1.05), n - 0.16, 0.34, color=CCOL[c], clip_on=False, lw=0))
        bx.text(start + (n - 1) / 2, -1.25, {"Mac/mono": "Mac", "Fibro": "Fib", "Peri/SMC": "Peri"}.get(SHORT[c], SHORT[c]),
                ha="center", va="bottom", fontsize=5.0, color=S.INK)
        start += n
    dotkeys(pg, 100, Y1 + 11.0, smax=15.0)

    # b: program enrichment (z of class means across classes)
    P = [("score_TLS", "TLS"), ("score_CYTO", "Cytotoxic"), ("score_EXHAUST", "Exhaustion"), ("score_IFN", "IFN"),
         ("score_CAF_ECM", "CAF/ECM"), ("score_TGFb", "TGF-β"), ("score_TAM_SUP", "Supp. TAM")]
    mm = obs.groupby("class", observed=True)[[p for p, _ in P]].mean().reindex(DISC)
    Zc = (mm - mm.mean()) / mm.std()
    REPORT["sc_atlas_b_z"] = Zc.round(2).to_dict()
    pg.head(116, 122, Y1, "b", "Program enrichment", "z across cell types")
    cx = pg.ax(143, Y1 + 9.0, 7 * 3.85, len(DISC) * 3.1)
    heat(cx, Zc, DISC, [p for p, _ in P], [LAB[c] for c in DISC], [l for _, l in P])
    zbar(pg, 172.6, Y1 + 11.5, h=11, label="z")
    pg.save("FigS_sc_atlas")


# ================================================================== sc_programmes
def fig_sc_programmes():
    obs = pd.read_parquet(D / "cells_obs.parquet")
    obs["class"] = obs.anno_level1.map(ANNO2CLASS)
    P = [("score_TLS", "TLS"), ("score_CYTO", "Cytotoxic"), ("score_EXHAUST", "Exhaustion"), ("score_IFN", "IFN"),
         ("score_CAF_ECM", "CAF/ECM"), ("score_TGFb", "TGF-β"), ("score_TAM_SUP", "Suppressive TAM")]
    pg = Page(180, 102)
    cw, ch, gx = 35.0, 30.0, 8.5
    for k, (col, lab) in enumerate(P):
        r, q = divmod(k, 4)
        x0 = 10 + q * (cw + gx); yb = 5.0 + r * 50.0
        pg.head(x0 - 9, x0 - 3, yb, "abcdefg"[k], lab)
        ax = pg.ax(x0, yb + 4.0, cw, ch)
        data = [obs.loc[obs["class"] == c, col].values for c in DISC]
        bp = ax.boxplot(data, showfliers=False, widths=0.64, patch_artist=True,
                        medianprops=dict(color=S.INK, lw=0.8), whiskerprops=dict(lw=0.5, color=S.GREY),
                        capprops=dict(lw=0.5, color=S.GREY), boxprops=dict(lw=0.0))
        for b, c in zip(bp["boxes"], DISC):
            b.set_facecolor(CCOL[c])
        ax.axhline(0, color=S.GREY, lw=0.4, ls=(0, (2, 2)), zorder=0)
        ax.set_xticks(range(1, len(DISC) + 1))
        ax.set_xticklabels([SHORT[c] for c in DISC], rotation=90, fontsize=5.0)
        ax.tick_params(axis="x", length=0, pad=1.8)
        if q == 0:
            ax.set_ylabel("Program score", fontsize=5.6)
    pg.save("FigS_sc_programmes")


# ================================================================== sc_composition
def fig_sc_composition():
    comp = pd.read_csv(GC / "out_scanpy_strict_GSE183904/composition_by_sample_anno_level1.csv").set_index("sample_id")
    cl = pd.concat([pd.read_csv(D / "sample_clinical_fibro.csv"), pd.read_csv(D / "sample_clinical_macro.csv")]) \
        .drop_duplicates("sample_id").set_index("sample_id")
    cl = cl[cl.sample_type == "Primary_Tumor"]
    m = comp.join(cl, how="inner")
    m["Lauren"] = m.Laurens.replace({"Mixed intestinal and diffuse": "Mixed", "N.A": np.nan})
    m["TCGA_"] = m.TCGA.replace({"N.A": np.nan})
    m["Stage_"] = m.Stage.replace({"DM": np.nan})
    ORDER = [KEY[c] for c in DISC]
    pg = Page(180, 66)
    specs = [("Stage_", "AJCC stage", ["I", "II", "III", "IV"]), ("Lauren", "Lauren type", ["Intestinal", "Diffuse", "Mixed"]),
             ("TCGA_", "TCGA subtype", ["CIN", "GS", "MSI"])]
    xs = [12, 66, 110]; ws = [44, 34, 34]
    minp = {}
    for k, ((var, lab, order), x0, w) in enumerate(zip(specs, xs, ws)):
        order = [o for o in order if (m[var] == o).sum() > 0]
        mean = m.groupby(var)[ORDER].mean().loc[order]
        ps = []
        for t in ORDER:
            groups = [m.loc[m[var] == o, t].values for o in order if (m[var] == o).sum() >= 2]
            if len(groups) >= 2:
                ps.append(kruskal(*groups).pvalue)
        minp[lab] = min(ps)
        pg.head(x0 - 10, x0 - 4, 5.0, "abc"[k], lab, f"Kruskal–Wallis, smallest {S.fmt_p(min(ps))}")
        ax = pg.ax(x0, 12.0, w, 43)
        bottom = np.zeros(len(order))
        for t, c in zip(ORDER, DISC):
            ax.bar(range(len(order)), mean[t], bottom=bottom, color=CCOL[c], width=0.74, edgecolor="white", lw=0.35)
            bottom += mean[t].values
        ns = m[var].value_counts()
        ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{o}\n(n = {ns.get(o, 0)})" for o in order], fontsize=5.2)
        ax.tick_params(axis="x", length=0, pad=2.0)
        ax.set_ylim(0, 1); ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(["0", "0.5", "1"])
        if k == 0:
            ax.set_ylabel("Mean fraction of cells", fontsize=5.6)
        ax.spines["bottom"].set_visible(False)
    pg.key(150, 14.0, [(CCOL[c], LAB[c]) for c in DISC[::-1]], ncol=1, dy=2.8)
    ref = {"AJCC stage": 0.096, "Lauren type": 0.041, "TCGA subtype": 0.051}
    assert all(abs(round(minp[k], 3) - v) < 1e-9 for k, v in ref.items()), minp
    REPORT["sc_composition"] = {"n_primary": len(m), "min_P": {k: round(v, 4) for k, v in minp.items()},
                                "n_by_group": {v: m[v].value_counts().to_dict() for v in ("Stage_", "Lauren", "TCGA_")}}
    pg.save("FigS_sc_composition")


# ================================================================== sc_fractions
def fig_sc_fractions():
    co = pd.read_csv(GC / "results/Fig5_FibroMacro/FibroMacro_patient_coenrichment.csv")
    ff = pd.read_csv(GC / "results/Fig5_Fibro/Fibro_patient_subtype_fractions.csv")
    mf = pd.read_csv(GC / "results/Fig5_Macro/Macro_patient_state_fractions.csv")
    j = ff.merge(mf, on="Patient_ID")
    pg = Page(120, 57)
    specs = [(co.BarrierCAF_frac, co.SupTAM_frac, "Barrier-type CAF fraction of fibroblasts",
              "Suppressive-TAM fraction of macrophages", "a", "Barrier-type CAFs and suppressive TAMs", S.STROMA),
             (j.CAF_ECM, j.SPP1_TAM, "ECM-remodeling CAF fraction of fibroblasts", "SPP1⁺ TAM fraction of macrophages",
              "b", "ECM-remodeling CAFs and SPP1⁺ TAMs", S.MACRO)]
    ref = {"a": (27, 0.15, 0.46), "b": (27, -0.11, 0.60)}
    for k, (x, y, xl, yl, L, tt, col) in enumerate(specs):
        x0 = 12 + k * 60
        r, p = spearmanr(x, y)
        assert (len(x), round(r, 2), round(p, 2)) == ref[L], (L, len(x), r, p)
        REPORT[f"sc_fractions_{L}"] = (len(x), round(r, 3), round(p, 3))
        pg.head(x0 - 11, x0 - 5, 5.0, L, tt,
                f"{len(x)} patients; Spearman ρ = {r:.2f}, {S.fmt_p(p)}".replace("-", "−"))
        ax = pg.ax(x0, 11.5, 44, 38)
        ax.scatter(x, y, s=9, color=S.INK, edgecolor="white", lw=0.35, zorder=3)
        ax.set_xlabel(xl, fontsize=5.4); ax.set_ylabel(yl, fontsize=5.4)
        ax.set_xlim(left=min(0, x.min()) - 0.02 * x.max()); ax.set_ylim(bottom=min(0, y.min()) - 0.02 * y.max())
    pg.save("FigS_sc_fractions")


# ================================================================== sc_indep
def fig_sc_indep():
    PS = pd.read_csv(S.RES / "E17_share_by_patient_summary.csv")
    PS = PS[PS.analysis == "primary"]
    pg = Page(180, 132)
    Y2, Y3 = 5.0, 74.5
    # a: per-patient (discovery: per-sample) intended-class share, each patient a dot, median and IQR
    PP = pd.read_csv(S.RES / "E17_share_by_patient.csv")
    PP = PP[PP.analysis == "primary"]
    GRP = {"CAF_ECM": (["fibroblast"], "fibroblast", "Fibroblasts"),
           "TAM_SUP": (["macrophage/monocyte"], "macrophage/monocyte", "Macrophages/monocytes"),
           "TLS_B": (["B cell", "plasma cell"], "B+plasma", "B + plasma cells"),
           "T_CELL": (["T cell", "NK cell"], "T+NK", "T + NK cells")}
    pg.head(1, 7, Y2, "a", "Intended-class share per patient",
            "Each dot one patient (discovery atlas: one tumor sample); bar, median; line, IQR")
    rowh, gap_out = 3.4, 2.6
    rows, y = [], 0.0
    for m in GRP:
        for ds in DATASETS:
            rows.append((m, ds, y)); y += rowh
        y += gap_out
    RH = y - gap_out
    RY = Y2 + 10.0
    ax = pg.ax(52, RY, 72, RH)
    rng = np.random.default_rng(3)
    stats = {}
    for m, ds, y0 in rows:
        cs, grp, _ = GRP[m]
        d = PP[(PP.dataset == ds) & (PP.module == m)]
        v = d[d["class"].isin(cs)].groupby("patient").share_pct.sum().reindex(d.patient.unique()).fillna(0)
        ref = PS[(PS.dataset == ds) & (PS.module == m) & (PS.class_group == grp)].iloc[0]
        assert abs(v.median() - ref["median"]) < 1e-9 and abs(v.quantile(.25) - ref.q1) < 1e-9 and (v >= 50).sum() == ref.n_ge50
        yc = y0 + rowh / 2
        ax.scatter(v.values, yc + rng.uniform(-0.75, 0.75, len(v)), s=4.5, color=FAM[m], alpha=0.55, lw=0, zorder=2)
        ax.plot([ref.q1, ref.q3], [yc, yc], color=S.INK, lw=1.1, solid_capstyle="butt", zorder=3)
        ax.plot([ref["median"]] * 2, [yc - 1.1, yc + 1.1], color=S.INK, lw=1.3, zorder=4)
        ax.text(103, yc, f"{ref['median']:.0f} ({ref.q1:.0f}–{ref.q3:.0f})", fontsize=5.0, va="center", ha="left")
        ax.text(127, yc, f"{int(ref.n_ge50)}/{int(ref.n_patients)}", fontsize=5.0, va="center", ha="left", color=S.SLATE)
        lab = {"GSE183904": "GSE183904 (samples)", "OMIX001073": "OMIX001073", "GSE206785": "GSE206785"}[ds]
        ax.text(-2, yc, lab, fontsize=5.0, va="center", ha="right", color=S.SLATE)
        stats[(m, ds)] = (round(ref["median"], 1), round(ref.q1, 1), round(ref.q3, 1), int(ref.n_ge50), int(ref.n_patients))
    for m in GRP:
        ys = [y0 for mm, ds, y0 in rows if mm == m]
        pg.chip(2.0, RY + ys[0] + 0.3, FAM[m], 1.3, ys[-1] + rowh - ys[0] - 0.6)
        pg.text(4.3, RY + (ys[0] + ys[-1] + rowh) / 2 - 1.3, m, fontsize=5.5, fontweight="bold", va="center")
        pg.text(4.3, RY + (ys[0] + ys[-1] + rowh) / 2 + 1.4, GRP[m][2], fontsize=5.0, color=S.SLATE, va="center")
        if m != "T_CELL":
            ax.axhline(ys[-1] + rowh + gap_out / 2, color=S.LIGHT, lw=0.4, zorder=0)
    ax.axvline(50, color=S.INK, lw=0.5, ls=(0, (2, 2)), zorder=1)
    ax.set_xlim(0, 100); ax.set_ylim(RH, 0)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Share of module signal in the intended class (%)", fontsize=5.4, labelpad=1.5)
    ax.text(103, -0.8, "Median (IQR)", fontsize=5.0, color=S.SLATE, va="bottom")
    ax.text(127, -0.8, "≥ 50%", fontsize=5.0, color=S.SLATE, va="bottom")
    REPORT["sc_indep_a"] = {f"{m}|{ds}": v for (m, ds), v in stats.items()}

    # c: OMIX001073 gene dot plot, intended classes boxed
    e = gexall[gexall.dataset == "OMIX001073"]
    mean = e.pivot(index="class", columns="gene", values="mean_expr").reindex(index=CLASSES, columns=GENES48).values
    pct = e.pivot(index="class", columns="gene", values="pct_expressing").reindex(index=CLASSES, columns=GENES48).values
    pg.head(1, 7, Y3, "b", "Module genes by class, OMIX001073", "Boxes, pre-specified intended class")
    cx = pg.ax(30, Y3 + 11.5, 133, len(CLASSES) * 3.1)
    dotplot(cx, mean, pct, CLASSES, GENES48, smax=12.0)
    cx.set_yticks(range(len(CLASSES))); cx.set_yticklabels([LAB[c] for c in CLASSES], fontsize=5.2)
    module_bands(cx, GENES48)
    for m in MODS:
        js = [GENES48.index(g) for g in MODULES[m]]
        for c in INTENDED.get(m, []):
            i = CLASSES.index(c)
            cx.add_patch(Rectangle((min(js) - 0.5, i - 0.5), len(js), 1, fill=False, edgecolor=S.INK, lw=0.6, zorder=4))
    dotkeys(pg, 167, Y3 + 10.5, smax=12.0)
    pg.save("FigS_sc_indep")


FIGS = {"sc_atlas": fig_sc_atlas, "sc_programmes": fig_sc_programmes, "sc_composition": fig_sc_composition,
        "sc_fractions": fig_sc_fractions, "sc_indep": fig_sc_indep}

if __name__ == "__main__":
    for k in (sys.argv[1:] or FIGS):
        FIGS[k]()
    import json
    print(json.dumps(REPORT, indent=1, default=str)[:4000])
