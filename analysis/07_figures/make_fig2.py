"""Figure 2 - single-cell definition and cellular attribution of the marker modules (GSE183904,
126,559 tumour cells), drawn from data (sc_extract.py outputs) in the house style.

Effect sizes (e, g) are Cohen's d between patient pseudobulks (mean programme score per patient
and state, patients with >= 10 cells of the state); open circles give the cell-level d."""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import bnstyle as bs
import figstyle as fs
from bnstyle import INK, GREY, STROMA, MACRO, title

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/data/sc"
COMP = "#6A4C9C"          # comparator state (iCAF, inflammatory macrophages); checked vs STROMA/MACRO
OTHER = ["#C9D1D8", "#B4BEC7", "#DDE2E7", "#A7B2BC"]
CT = fs.CT; ORDER = fs.CTORDER
CTNAME = {"Tcell": "T cells", "NK": "NK cells", "Bcell": "B cells", "Plasma": "Plasma cells", "Macro": "Macrophages",
          "DC": "Dendritic cells", "Fibro": "Fibroblasts", "Endo": "Endothelial", "Peri": "Pericytes", "Epi": "Epithelial"}
MOD = dict(CAF_ECM=["COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"],
           TAM_SUP=["SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"],
           APM=["B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"], IFN=["STAT1","IRF1","CXCL9","CXCL10","IFNG"],
           TLS_B=["CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"], T_CELL=["CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"])
MODCOL = {"CAF_ECM": STROMA, "TAM_SUP": MACRO, "APM": bs.IMMUNE, "IFN": bs.IMMUNE, "TLS_B": bs.IMMUNE, "T_CELL": bs.IMMUNE}
CANON = [("Tcell", ["CD3D", "CD3E"]), ("NK", ["NKG7", "GNLY"]), ("Bcell", ["MS4A1", "CD79A"]), ("Plasma", ["MZB1", "JCHAIN"]),
         ("Macro", ["CD68", "C1QA", "LYZ"]), ("DC", ["LST1", "CLEC10A"]), ("Fibro", ["COL1A1", "DCN"]),
         ("Endo", ["PECAM1", "VWF"]), ("Peri", ["RGS5", "MCAM"]), ("Epi", ["EPCAM", "KRT19", "KRT8"])]

obs = pd.read_parquet(f"{D}/cells_obs.parquet")
E = pd.read_parquet(f"{D}/cells_expr.parquet")
ct = obs.anno_level1.astype(str)
fib = pd.read_parquet(f"{D}/fibro_obs.parquet"); mac = pd.read_parquet(f"{D}/macro_obs.parquet")


def umap_all(ax):
    rng = np.random.default_rng(0); o = obs.iloc[rng.permutation(len(obs))]
    ax.scatter(o.umap1, o.umap2, s=0.03, c=[CT[t] for t in o.anno_level1.astype(str)], lw=0, rasterized=True)
    for t in ORDER:
        ax.plot([], [], "o", color=CT[t], ms=3.2, label=CTNAME[t])
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=5.0, handletextpad=0.15,
              columnspacing=0.6, labelspacing=0.25)
    ax.set_axis_off(); ax.set_aspect("equal", adjustable="datalim")
    title(ax, f"{len(obs):,} cells, 10 types")


def dotplot(ax, groups, sub_E, colmap=None, gene_label_rot=90, bracket=None):
    genes = [g for _, gs in groups for g in gs if g in sub_E.columns]
    mean = sub_E[genes].groupby(ct).mean().loc[ORDER]
    frac = (sub_E[genes] > 0).groupby(ct).mean().loc[ORDER]
    scaled = (mean - mean.min()) / (mean.max() - mean.min()).replace(0, 1)
    for j, g in enumerate(genes):
        ax.scatter([j] * len(ORDER), range(len(ORDER)), s=frac[g].values * 26 + 0.3, c=scaled[g].values,
                   cmap=fs.SEQ, vmin=0, vmax=1, lw=0.2, edgecolors=GREY)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=gene_label_rot, fontsize=5.0, style="italic")
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([CTNAME[t] for t in ORDER], fontsize=5.4)
    ax.set_xlim(-0.7, len(genes) - 0.3); ax.set_ylim(len(ORDER) - 0.5, -0.5)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    # group brackets above
    x = 0
    for lab, gs in groups:
        n = len([g for g in gs if g in sub_E.columns])
        if n == 0: continue
        col = (colmap or {}).get(lab, CT.get(lab, INK))
        ax.plot([x - 0.35, x + n - 0.65], [-0.95, -0.95], color=col, lw=1.2, clip_on=False)
        ax.text(x + (n - 1) / 2, -1.25, bracket(lab) if bracket else CTNAME.get(lab, lab), fontsize=5.0, ha="center",
                va="bottom", color=col, clip_on=False)
        x += n
    return genes


def canon(ax):
    dotplot(ax, CANON, E, bracket=lambda l: {"Tcell": "T", "NK": "NK", "Bcell": "B", "Plasma": "Plasma", "Macro": "Macro",
                                              "DC": "DC", "Fibro": "Fibro", "Endo": "Endo", "Peri": "Peri", "Epi": "Epi"}[l])
    title(ax, "Canonical markers", pad=16)
    # size / colour keys
    for i, fr in enumerate([0.25, 0.5, 0.75, 1.0]):
        ax.scatter(len(sum([g for _, g in CANON], [])) + 1.2, 1.2 + i * 1.3, s=fr * 26 + 0.3, c=[GREY], clip_on=False)
        ax.text(len(sum([g for _, g in CANON], [])) + 2.0, 1.2 + i * 1.3, f"{int(fr*100)}%", fontsize=5.0, va="center", clip_on=False)
    ax.text(len(sum([g for _, g in CANON], [])) + 0.8, 0.1, "cells\nexpressing", fontsize=5.0, va="center", clip_on=False)


def prog_heat(ax):
    P = {"score_TLS": "TLS", "score_CYTO": "Cytotoxic", "score_EXHAUST": "Exhaustion", "score_IFN": "IFN",
         "score_CAF_ECM": "CAF/ECM", "score_TGFb": "TGF-β", "score_TAM_SUP": "Supp. TAM"}
    m = obs[list(P)].groupby(ct).mean().loc[ORDER]
    z = (m - m.mean()) / m.std()
    im = ax.imshow(z.values, cmap=fs.DIV, vmin=-2.5, vmax=2.5, aspect="auto")
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            v = z.values[i, j]
            ax.text(j, i, f"{v:.1f}".replace("-", "−").replace("−0.0", "0.0"), ha="center", va="center", fontsize=5.0, color="white" if abs(v) > 1.6 else INK)
    ax.set_xticks(range(len(P))); ax.set_xticklabels(P.values(), rotation=45, ha="right", fontsize=5.2)
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([CTNAME[t] for t in ORDER], fontsize=5.2)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03); cb.ax.tick_params(labelsize=5.0, length=1.5); cb.outline.set_visible(False)
    cb.set_label("z across cell types", fontsize=5.0)
    title(ax, "Programme enrichment")


def umap_sub(o, lab, focal, comp, fcol, name):
    def draw(ax):
        others = [s for s in o[lab].unique() if s not in (focal, comp)]
        cmap = {s: OTHER[i % len(OTHER)] for i, s in enumerate(sorted(others))}
        cmap[focal] = fcol; cmap[comp] = COMP
        rng = np.random.default_rng(1); oo = o.iloc[rng.permutation(len(o))]
        oo = pd.concat([oo[~oo[lab].isin([focal, comp])], oo[oo[lab].isin([focal, comp])]])
        ax.scatter(oo.umap1, oo.umap2, s=0.6, c=[cmap[s] for s in oo[lab]], lw=0, rasterized=True)
        nice = {"CAF_ECM": "ECM-remodelling CAFs", "iCAF": "iCAFs", "SPP1_TAM": "SPP1⁺ TAMs", "Inflam": "inflammatory"}
        for s_, c_ in [(focal, fcol), (comp, COMP)]:
            ax.plot([], [], "o", color=c_, ms=3.2, label=f"{nice[s_]} ({(o[lab] == s_).sum():,})")
        ax.plot([], [], "o", color=OTHER[0], ms=3.2, label="other states")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=5.0, handletextpad=0.15, columnspacing=0.6)
        ax.set_axis_off(); ax.set_aspect("equal", adjustable="datalim")
        title(ax, f"{name}, {len(o):,} cells")
    return draw


def cohen(x, y):
    x, y = np.asarray(x), np.asarray(y); nx, ny = len(x), len(y)
    sp = np.sqrt(((nx - 1) * x.var(ddof=1) + (ny - 1) * y.var(ddof=1)) / (nx + ny - 2)); return (x.mean() - y.mean()) / sp


EFFECTS = {}


def effects(o, lab, focal, comp, progs, fcol, names, xlab, key):
    def draw(ax):
        pb = o.groupby(["Patient_ID", lab], observed=True)[progs].mean()
        n = o.groupby(["Patient_ID", lab], observed=True).size(); pb = pb[n >= 10]
        A, B = pb.xs(focal, level=1), pb.xs(comp, level=1)
        dp = [cohen(A[p], B[p]) for p in progs]
        dc = [cohen(o.loc[o[lab] == focal, p], o.loc[o[lab] == comp, p]) for p in progs]
        EFFECTS[key] = dict(zip(progs, zip(np.round(dp, 2), np.round(dc, 2)))); EFFECTS[key + "_n"] = (len(A), len(B))
        y = np.arange(len(progs))[::-1]
        ax.barh(y, dp, color=[fcol if v > 0 else COMP for v in dp], height=0.62)
        ax.plot(dc, y, "o", mfc="white", mec=INK, ms=2.6, mew=0.6)
        ax.axvline(0, color=INK, lw=0.5)
        ax.set_yticks(y); ax.set_yticklabels(names, fontsize=5.3); ax.tick_params(axis="y", length=0)
        ax.set_xlabel(xlab, fontsize=5.6)
        ax.plot([], [], "s", color=GREY, ms=4, label=f"patients ({len(A)} vs {len(B)})")
        ax.plot([], [], "o", mfc="white", mec=INK, ms=2.6, label="cells")
        ax.legend(loc="lower left", fontsize=5.0, handletextpad=0.3, borderaxespad=0.3)
        title(ax, "Programme scores, effect sizes")
    return draw


def attribution(ax):
    sc = pd.DataFrame({m: E[[g for g in G if g in E.columns]].mean(axis=1) for m, G in MOD.items()})
    m = sc.groupby(ct).mean().loc[ORDER]
    z = (m - m.mean()) / m.std()
    z.to_csv("data/sc/attribution_z_recomputed.csv")
    im = ax.imshow(z.values, cmap=fs.DIV, vmin=-2.8, vmax=2.8, aspect="auto")
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            v = z.values[i, j]
            ax.text(j, i, f"{v:.1f}".replace("-", "−").replace("−0.0", "0.0"), ha="center", va="center", fontsize=5.0, color="white" if abs(v) > 1.8 else INK)
    ax.set_xticks(range(6)); ax.set_xticklabels(list(MOD), rotation=45, ha="right", fontsize=5.2)
    for t, c in zip(ax.get_xticklabels(), MODCOL.values()): t.set_color(c)
    ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([CTNAME[t] for t in ORDER], fontsize=5.2)
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.plot([-0.4, 1.4], [-0.85, -0.85], color=STROMA, lw=1.1, clip_on=False)
    ax.text(0.5, -1.05, "Barrier", fontsize=5.2, ha="center", va="bottom", color=STROMA, clip_on=False)
    ax.plot([1.6, 5.4], [-0.85, -0.85], color=bs.IMMUNE, lw=1.1, clip_on=False)
    ax.text(3.5, -1.05, "Immune", fontsize=5.2, ha="center", va="bottom", color=bs.IMMUNE, clip_on=False)
    cb = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.03); cb.ax.tick_params(labelsize=5.0, length=1.5); cb.outline.set_visible(False)
    cb.set_label("z across cell types", fontsize=5.0)
    title(ax, "Module attribution", pad=12)


def modgenes(ax):
    genes = dotplot(ax, [(m, G) for m, G in MOD.items()], E, colmap=MODCOL, bracket=lambda l: l)
    title(ax, "The barriernet module genes (42 of 48 present in the atlas)", pad=14)


fig = plt.figure(figsize=(bs.W2, 200 * bs.MM))
gs = GridSpec(4, 1, figure=fig, height_ratios=[1.0, 0.95, 0.95, 0.72], hspace=0.62)
r1 = gs[0].subgridspec(1, 2, width_ratios=[0.85, 1.6], wspace=0.28)
r2 = gs[1].subgridspec(1, 3, width_ratios=[1.0, 0.9, 0.95], wspace=0.55)
r3 = gs[2].subgridspec(1, 3, width_ratios=[0.9, 0.95, 1.0], wspace=0.6)
axes = dict(a=fig.add_subplot(r1[0]), b=fig.add_subplot(r1[1]), c=fig.add_subplot(r2[0]), d=fig.add_subplot(r2[1]),
            e=fig.add_subplot(r2[2]), f=fig.add_subplot(r3[0]), g=fig.add_subplot(r3[1]), h=fig.add_subplot(r3[2]),
            i=fig.add_subplot(gs[3]))
umap_all(axes["a"]); canon(axes["b"]); prog_heat(axes["c"])
umap_sub(fib, "fibro_subtype", "CAF_ECM", "iCAF", STROMA, "Fibroblasts")(axes["d"])
effects(fib, "fibro_subtype", "CAF_ECM", "iCAF", ["prog_CAF_ECM", "prog_iCAF", "prog_STF3_FAP_INHBA", "prog_myCAF", "prog_CAF_TGFB"],
        STROMA, ["ECM-remodelling", "iCAF", "FAP⁺INHBA⁺ CAF", "myCAF", "TGF-β CAF"], "Cohen's d, ECM-remodelling CAFs vs iCAFs", "fibro")(axes["e"])
umap_sub(mac, "macro_state", "SPP1_TAM", "Inflam", MACRO, "Macrophages")(axes["f"])
effects(mac, "macro_state", "SPP1_TAM", "Inflam", ["prog_SPP1_TAM", "prog_Inflam", "prog_HLAII_APC", "prog_IL10_STAT3", "prog_IFN"],
        MACRO, ["SPP1⁺ TAM", "Inflammatory", "HLA-II APC", "IL-10/STAT3", "IFN"], "Cohen's d, SPP1⁺ TAMs vs inflammatory", "macro")(axes["g"])
attribution(axes["h"]); modgenes(axes["i"])
for L, ax in axes.items():
    dx = {"a": -0.1, "b": -0.2, "c": -0.42, "d": -0.12, "e": -0.55, "f": -0.12, "g": -0.55, "h": -0.52, "i": -0.075}[L]
    dy = {"b": 1.3, "i": 1.42, "h": 1.24}.get(L, 1.1)
    bs.panel(ax, L, dx=dx, dy=dy)
bs.save(fig, "Fig2")
pd.Series({k: str(v) for k, v in EFFECTS.items()}).to_csv("data/sc/Fig2_effect_sizes.csv")
print(EFFECTS)
