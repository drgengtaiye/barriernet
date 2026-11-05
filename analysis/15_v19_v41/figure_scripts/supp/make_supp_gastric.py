"""v19 Supplementary Figures, gastric prognosis and gastro-esophageal response, in the cns house style.

meta_forest   random-effects meta-analyses of six scores in the benchmark cohorts
fitted_models fitted transcriptome models against the fixed score
upgrades      attempted upgrades (cell-state axes, abundance/state elastic net, arrangement score)
robust        robustness to gene loss (E19)
locked        locked replication in seven gastric cohorts (E15, E20)
prjeb_detail  pembrolizumab response in PRJEB25780 in detail (E3, E7)
gastric_icb   gastro-esophageal checkpoint cohorts (E22)

Run:  python make_supp_gastric.py [key ...]   (no key = all). Every estimate is read from existing result files
(analysis/results, revision/results, revision2/results); ROC curves are descriptive re-plots of score tables."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
import io, subprocess
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import Affine2D
from matplotlib.ticker import NullFormatter
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, MACRO, FITTED, fmt_p

A = S.RES
REV = S.ROOT.parent / "revision/results"
REV2 = S.ROOT.parent / "revision2/results"
OUT = S.OUT / "supp"
OUT.mkdir(parents=True, exist_ok=True)
FS = 5.2
LOW = "#7A8794"
NEUTRAL = "#8E9AA6"                  # chip for comparators that belong to no module family


# ------------------------------------------------------------------ shared helpers
class Fig:
    """A figure laid out in millimetres, with letters placed after layout and a numbers log."""

    def __init__(self, key, h_mm, w_mm=180.0):
        self.key, self.H, self.W = key, h_mm, w_mm
        self.fig = plt.figure(figsize=(w_mm * S.MM, h_mm * S.MM))
        self.LET, self.NUM = [], {}

    def ax(self, x, y, w, h):
        return S.ax_mm(self.fig, x, y, w, h)

    def rec(self, k, v, src):
        self.NUM[k] = (v, src)

    def _mm(self):
        return Affine2D().scale(S.MM) + self.fig.dpi_scale_trans

    def chip(self, x_mm, y_mm, col):
        """1.1 x 1.6 mm rounded bar; x from left, y from the bottom (mm)."""
        self.fig.add_artist(FancyBboxPatch((x_mm, y_mm - 0.8), 1.1, 1.6, boxstyle="round,pad=0,rounding_size=0.3",
                                           transform=self._mm(), fc=col, ec="none"))

    def ymm(self, ax, y):
        bb = ax.get_position(); y0, y1 = ax.get_ylim()
        if ax.get_yscale() == "log":
            y, y0, y1 = np.log(y), np.log(y0), np.log(y1)
        return (bb.y0 + (y - y0) / (y1 - y0) * bb.height) * self.H

    def row_label(self, ax, y, text, col=None, gap=1.2, fs=FS, bold=False, color=INK):
        """Ink label right-aligned left of the axes; family chip between label and axes."""
        xr = ax.get_position().x0 * self.W - gap
        ym = self.ymm(ax, y)
        if col:
            self.chip(xr - 1.1, ym, col)
        return self.fig.text((xr - (2.0 if col else 0.0)) / self.W, ym / self.H, text, ha="right", va="center",
                             fontsize=fs, color=color, fontweight="bold" if bold else "normal", linespacing=1.1)

    def title(self, ax, text, col=None, sub=None, dy=None, fs=6.0, x_mm=None):
        """Bold ink title above the axes (chip first when a family colour is given); subtitles go to the legend."""
        sub = None
        bb = ax.get_position()
        x0 = bb.x0 * self.W if x_mm is None else x_mm
        top = bb.y1 * self.H
        yt = top + (dy if dy is not None else (4.8 if sub else 2.6))
        if col:
            self.chip(x0, yt + 0.9, col)
        t = self.fig.text((x0 + (2.0 if col else 0.0)) / self.W, yt / self.H, text, ha="left", va="baseline",
                          fontsize=fs, fontweight="bold", color=INK)
        if sub:
            self.fig.text(x0 / self.W, (top + 1.4) / self.H, sub, ha="left", va="baseline", fontsize=5.0,
                          color=SLATE)
        return t

    def text_at(self, x_mm, y_mm_top, s, **kw):
        kw.setdefault("fontsize", 5.0); kw.setdefault("color", INK)
        return self.fig.text(x_mm / self.W, 1 - y_mm_top / self.H, s, **kw)

    def letter(self, s, art, x_mm):
        self.LET.append((s, art, x_mm))

    def save(self):
        fig = self.fig
        fig.canvas.draw()
        r = fig.canvas.get_renderer(); k = fig.dpi / 72.0
        asc = TextPath((0, 0), "l", prop=FontProperties(family=plt.rcParams["font.family"], size=8,
                                                         weight="bold")).get_extents().y1
        for s, art, x in self.LET:
            prop = art.get_fontproperties()
            d = max(r.get_text_width_height_descent("lp", prop, ismath=False)[2],
                    r.get_text_width_height_descent(art.get_text(), prop, ismath=False)[2])
            base = art.get_window_extent(r).y0 + d
            cap = TextPath((0, 0), "H", prop=prop).get_extents().y1
            fig.text(x / self.W, (base + (cap - asc) * k) / fig.bbox.height, s, fontsize=8, fontweight="bold",
                     va="baseline", ha="left", color=INK)
        stem = f"FigS_{self.key}"
        fig.savefig(OUT / f"{stem}.pdf"); fig.savefig(OUT / f"{stem}.png", dpi=600)
        plt.close(fig)
        with open(OUT / f"{stem}_numbers.md", "w") as fh:
            fh.write(f"# Supplementary figure '{self.key}' - numbers printed on the figure and their sources\n\n"
                     "| item | value | source |\n|---|---|---|\n")
            for kk, (v, src) in self.NUM.items():
                fh.write(f"| {kk} | {v} | {src} |\n")
        print("saved", stem)


def rule(ax, x):
    S.rule(ax, x)


def logx(ax, lim, ticks):
    ax.set_xscale("log"); ax.set_xlim(*lim); ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks]); ax.xaxis.set_minor_formatter(NullFormatter())


def point(ax, v, lo, hi, y, c, mk="o", ms=3.0, open_=False, lim=None, lw=0.8, ls="-"):
    if lim is not None:
        a, b = lim
        if lo < a:
            ax.plot(a, y, "<", ms=2.2, color=c, mew=0, clip_on=False); lo = a
        if hi > b:
            ax.plot(b, y, ">", ms=2.2, color=c, mew=0, clip_on=False); hi = b
    ax.plot([lo, hi], [y, y], color=c, lw=lw, ls=ls, solid_capstyle="butt", zorder=2)
    ax.plot(v, y, mk, ms=ms, color=c, mfc="white" if open_ else c, mec=c if open_ else "white",
            mew=0.8 if open_ else 0.4, zorder=3)


def diamond(ax, est, lo, hi, y, c, h=0.30, filled=True, lw=0.7):
    ax.fill([lo, est, hi, est], [y, y - h, y, y + h], color=c if filled else "white", ec=c, lw=lw, zorder=3)


def noy(ax):
    ax.set_yticks([]); ax.spines["left"].set_visible(False)


def neg(s):
    return s.replace("-", "−")


def rds_csv(path):
    out = subprocess.run(["/users/yc4229/.conda/envs/r45/bin/Rscript", "-e",
                          f"write.csv(as.data.frame(readRDS('{path}')), stdout(), row.names=FALSE)"],
                         capture_output=True, text=True, check=True).stdout
    return pd.read_csv(io.StringIO(out))


# ================================================================== meta_forest
def meta_forest():
    F = Fig("meta_forest", 88.0)
    uni = pd.read_csv(REV / "B4_univariate.csv"); meta = pd.read_csv(REV / "B4_meta.csv")
    MK = pd.read_csv(A / "E8_meta_DL_HK.csv").set_index("score")
    KEY = {"ESTIMATE(imm-str)": "ESTIMATE", "NetImmune(ssGSEA)": "NetImmune_ssGSEA"}

    def conservative(sc):
        r = MK.loc[KEY.get(sc, sc)]
        dl, hk = np.log(r.DL_hi / r.DL_lo), np.log(r.HK_hi / r.HK_lo)
        lo, hi, rule_ = (r.HK_lo, r.HK_hi, "HK") if hk >= dl else (r.DL_lo, r.DL_hi, "DL")
        return r.HR, lo, hi, max(r.DL_p, r.HK_p), rule_

    SC = [("Barrier", "Barrier", STROMA),
          ("Immune", "Immune (no Barrier term)", IMMUNE), ("CD8A", "CD8A", NEUTRAL),
          ("ESTIMATE(imm-str)", "ESTIMATE immune−stromal", NEUTRAL),
          ("NetImmune(ssGSEA)", "NetImmune, ssGSEA scoring", IMMUNE)]
    COH = [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254")]
    wmax = max(1 / uni.se.values ** 2)
    for k, (sc, name, col) in enumerate(SC):
        r_, c_ = divmod(k, 3)
        ax = F.ax(22 + c_ * 60, 14 + r_ * 42, 34, 22)
        hrp, plo, phi, pp, prule = conservative(sc)
        mc = INK if col == NEUTRAL else col
        for i, (c, cn) in enumerate(COH):
            u = uni[(uni.cohort == c) & (uni.score == sc)].iloc[0]
            ax.plot([u.lo, u.hi], [i, i], color=mc, lw=0.8, zorder=2)
            ax.plot(u.HR, i, "s", ms=1.8 + 2.4 * np.sqrt((1 / u.se ** 2) / wmax), color=mc, mec="white", mew=0.3,
                    zorder=3)
            F.rec(f"{name} {cn}", f"{u.HR:.2f} ({u.lo:.2f}–{u.hi:.2f})", "revision/results/B4_univariate.csv")
        m = meta[meta.score == sc].iloc[0]
        diamond(ax, hrp, plo, phi, 3.15, mc, h=0.32, filled=pp < 0.05)
        rule(ax, 1)
        logx(ax, (0.55, 1.6), [0.6, 0.8, 1, 1.25, 1.5])
        ax.set_ylim(3.7, -0.6); noy(ax)
        for i, (c, cn) in enumerate(COH):
            F.row_label(ax, i, cn)
        F.row_label(ax, 3.15, "Pooled", bold=True)
        if r_ == 1 or k == 2:
            ax.set_xlabel("HR per SD", labelpad=1)
        t = F.title(ax, name, col=col, sub=neg(f"pooled {hrp:.2f} ({plo:.2f}–{phi:.2f}); {fmt_p(pp)}; "
                                                f"I² = {m.I2:.0f}%"), x_mm=22 + c_ * 60 - 17)
        F.rec(f"{name} pooled (conservative, {prule} interval)", f"{hrp:.2f} ({plo:.2f}–{phi:.2f}); {fmt_p(pp)}; "
              f"I2 {m.I2:.0f}%", "E8_meta_DL_HK.csv (wider interval, larger P); I2 from revision/results/B4_meta.csv")
        F.letter("abcdef"[k], t, 1.0 if c_ == 0 else 22 + c_ * 60 - 21.5)
    F.save()


# ================================================================== fitted_models
def fitted_models():
    F = Fig("fitted_models", 78.0)
    COH = [("TCGA", "o", "TCGA-STAD"), ("GSE84437", "^", "GSE84437"), ("GSE62254", "v", "GSE62254")]
    # a  model classes
    d = pd.read_csv(A / "E6_model_classes.csv")
    w = d.pivot(index="model", columns="test", values="C")
    w["mean"] = w[[c for c, _, _ in COH]].mean(1)
    w = w.sort_values("mean", ascending=False)
    ax = F.ax(30, 12, 52, 48)
    ax.set_ylim(len(w) - 0.5, -0.7)
    net = w.loc["NetImmune (fixed)", "mean"]
    ax.axvline(net, color=IMMUNE, lw=0.6, ls=(0, (2, 2)), zorder=0)
    names = {"Elastic net (alpha 0.5)": "Elastic net (α = 0.5)", "Boosted Cox (xgboost)": "Boosted Cox",
             "PC-Cox": "Cox on 20 PCs"}
    for y, (m, r) in enumerate(w.iterrows()):
        fixed = "(fixed)" in m
        for j, (c, mk, _) in enumerate(COH):
            ax.plot(r[c], y - 0.22 + 0.11 * j, mk, color=GREY, ms=2.2, mew=0, zorder=2)
        col = IMMUNE if m.startswith("NetImmune") else (INK if fixed else FITTED)
        ax.plot(r["mean"], y + 0.18, "o" if fixed else "D", color=col, ms=3.2, mec="white", mew=0.4, zorder=3)
        ax.text(1.02, y, f"{r['mean']:.3f}", transform=ax.get_yaxis_transform(), fontsize=5.0, va="center", ha="left",
                color=INK, fontweight="bold" if m.startswith("NetImmune") else "normal")
        F.rec(f"a {m} mean held-out C", f"{r['mean']:.3f}", "E6_model_classes.csv")
    noy(ax)
    for y, m in enumerate(w.index):
        F.row_label(ax, y, names.get(m, m), col=IMMUNE if m.startswith("NetImmune") else None, fs=5.2)
    ax.set_xlim(0.49, 0.67); ax.set_xticks([0.5, 0.55, 0.6, 0.65])
    ax.set_xlabel("Harrell's C in the held-out cohort", labelpad=1)
    t = F.title(ax, "Model classes", sub="leave one cohort out", x_mm=4.0)
    F.letter("a", t, 1.0)
    ax.text(1.02, 1.0, "mean", transform=ax.transAxes, fontsize=5.0, color=SLATE, va="bottom")
    # key
    kx = 6.0
    for i, (c, mk, lab) in enumerate(COH):
        ax.plot(0.02 + i * 0.26, -0.20, mk, ms=2.2, color=GREY, mew=0, transform=ax.transAxes, clip_on=False)
        ax.text(0.045 + i * 0.26, -0.20, lab, fontsize=5.0, color=SLATE, va="center", transform=ax.transAxes)
    ax.plot(0.02, -0.265, "D", ms=2.8, color=FITTED, transform=ax.transAxes, clip_on=False)
    ax.text(0.045, -0.265, "fitted, mean of three", fontsize=5.0, color=SLATE, va="center", transform=ax.transAxes)
    ax.plot(0.54, -0.265, "o", ms=2.8, color=INK, transform=ax.transAxes, clip_on=False)
    ax.text(0.565, -0.265, "fixed, mean", fontsize=5.0, color=SLATE, va="center", transform=ax.transAxes)
    # b  gain over the fixed score
    dc = pd.read_csv(A / "E6_review_dC.csv"); wc = pd.read_csv(A / "E6_within_cohort_cv.csv")
    rows = [("Elastic net", "Elastic net", None), ("Best alternative: pcacox", "Cox on 20 PCs", None),
            ("NetImmune offset + elastic net", "NetImmune + elastic net\n(offset)", IMMUNE),
            ("Elastic net, inverse-normal genes", "Elastic net,\ninverse-normal genes", None)]
    ax = F.ax(120, 12, 42, 40)
    ax.set_ylim(len(rows) + 0.6, -0.6)
    for i, (key, lab, chipc) in enumerate(rows):
        s = dc[dc.model == key]
        for j, (c, mk, _) in enumerate(COH):
            r = s[s.test == c].iloc[0]
            ax.plot(r.dC, i - 0.22 + 0.11 * j, mk, color=GREY, ms=2.2, mew=0, zorder=2)
        p = s[s.test == "pooled"].iloc[0]
        col = IMMUNE if "offset" in key else FITTED
        ax.plot([p.lo, p.hi], [i + 0.18] * 2, color=col, lw=0.9, zorder=2)
        ax.plot(p.dC, i + 0.18, "D", color=col, ms=3.0, mec="white", mew=0.4, zorder=3)
        ax.text(1.03, i, neg(f"{p.dC:+.3f}") + f"\n{fmt_p(p.p)}", transform=ax.get_yaxis_transform(), fontsize=5.0,
                ha="left", va="center", color=INK, linespacing=1.15)
        F.rec(f"b {lab.replace(chr(10), ' ')}", f"{p.dC:+.3f} ({p.lo:+.3f} to {p.hi:+.3f}), {fmt_p(p.p)}",
              "E6_review_dC.csv pooled")
    y = len(rows)
    for j, (c, mk, _) in enumerate(COH):
        r = wc[wc.test == c].iloc[0]
        ax.plot(r.C_cv_mean - r.C_net, y - 0.22 + 0.11 * j, mk, mfc="white", mec=GREY, ms=2.4, mew=0.6, zorder=2)
    m = (wc.C_cv_mean - wc.C_net).mean()
    ax.plot(m, y + 0.18, "D", mfc="white", mec=FITTED, ms=3.0, mew=0.8, zorder=3)
    ax.text(1.03, y, neg(f"{m:+.3f}") + "\nno interval", transform=ax.get_yaxis_transform(), fontsize=5.0, ha="left",
            va="center", color=INK, linespacing=1.15)
    F.rec("b within-cohort CV elastic net, mean dC", f"{m:+.3f}", "E6_within_cohort_cv.csv")
    rule(ax, 0)
    noy(ax)
    for i, (key, lab, chipc) in enumerate(rows):
        F.row_label(ax, i, lab, col=chipc, fs=5.2)
    F.row_label(ax, len(rows), "Elastic net, within-\ncohort 10-fold CV", fs=5.2)
    ax.set_xlim(-0.06, 0.09); ax.set_xticks([-0.05, 0, 0.05]); ax.set_xticklabels(["−0.05", "0", "0.05"])
    ax.set_xlabel("ΔC versus NetImmune (95% CI)", labelpad=1)
    t = F.title(ax, "Gain over the fixed score", sub="paired bootstrap, pooled over held-out cohorts", x_mm=95.0)
    F.letter("b", t, 90.5)
    F.save()



# ================================================================== upgrades
def upgrades():
    from scipy.stats import pearsonr
    F = Fig("upgrades", 146.0)
    G = S.ROOT.parent / "GC_TME_Reprogrammer/results/ring_barrier_state_v2/BarrierIndex_by_slice_with_clinical.csv"
    GC = dict(zip(pd.read_csv(G).slice, pd.read_csv(G).GC_raw))
    L = pd.read_csv(REV2 / "S3_state_loadings.csv", index_col=0)
    coef = pd.read_csv(REV2 / "S5_coefs.csv", index_col=0)
    loco = pd.read_csv(REV2 / "S5_loco_cindex.csv").set_index("model")
    pred = pd.read_csv(REV2 / "S8_neigh_pred.csv")
    w = pd.read_csv(REV2 / "S9_arrangement_weights.csv", index_col=0).iloc[:, 0]
    sas = pd.read_csv(REV2 / "S10_prognosis.csv")
    AX = [("Fibro_S2", "Fibroblast axis 2", STROMA), ("Macro_S2", "Macrophage axis 2", MACRO),
          ("Tcell_S1", "T-cell axis 1", IMMUNE), ("NK_S1", "NK-cell axis 1", IMMUNE)]
    NAMED = {"Fibro_S2": ["CFD", "TNXB", "C3", "PTGDS", "FAP", "TNC", "IL11"],
             "Macro_S2": ["CXCL5", "CCL18", "MMP9", "IL1A", "F13A1", "P2RY13"],
             "Tcell_S1": ["TNFRSF4", "TNFRSF18", "IL2RA", "IFNG", "CRTAM"], "NK_S1": ["PDCD1", "LAG3", "CXCR3"]}
    # a  cell-state axis loadings
    for k, (a, lab, col) in enumerate(AX):
        ax = F.ax(17 + k * 42, 11, 24, 36)
        v = L[a].dropna()
        keep = set(v.nlargest(5).index) | set(v.nsmallest(5).index) | set(NAMED[a])
        top = v[list(keep)].sort_values(ascending=False)
        ax.barh(range(len(top)), top.values, color=[GREY if x < 0 else INK for x in top.values], height=0.7)
        ax.set_ylim(len(top) - 0.5, -0.6)
        ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=5.0, style="italic")
        for tl in ax.get_yticklabels():
            if tl.get_text() in NAMED[a]:
                tl.set_fontweight("bold")
        ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False)
        ax.axvline(0, color=INK, lw=0.4)
        ax.set_xlabel("Loading", labelpad=1)
        t = F.title(ax, lab, col=col, x_mm=17 + k * 42 - 11.5)
        if k == 0:
            F.letter("a", t, 1.0)
        F.rec(f"a {lab} genes shown", ", ".join(top.index), "revision2/results/S3_state_loadings.csv")
    # b  elastic-net coefficients on abundance and state
    ax = F.ax(30, 62, 78, 17)
    cf = coef.loc[(coef != 0).any(axis=1)]
    im = ax.imshow(cf.T.values, cmap=S.DIV, vmin=-0.2, vmax=0.2, aspect="auto")
    ax.set_xticks(range(len(cf)))
    ax.set_xticklabels([x.replace("ab_", "abundance ").replace("st_", "state ").replace("_S", " S")
                        for x in cf.index], rotation=60, ha="right", fontsize=5.0, rotation_mode="anchor")
    ax.set_yticks(range(3)); ax.set_yticklabels([{"TCGA": "TCGA-STAD"}.get(c, c) for c in cf.columns], fontsize=5.0)
    ax.tick_params(length=0)
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    kx = F.ax(110, 62, 1.8, 17)
    cb = plt.colorbar(im, cax=kx); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.0, length=1.5)
    cb.set_ticks([-0.2, 0, 0.2]); cb.set_ticklabels(["−0.2", "0", "0.2"])
    cb.set_label("Coefficient", fontsize=5.0, labelpad=1)
    mC = loco.loc["EN: abundance + state (30)", "mean"]
    t = F.title(ax, "Elastic net on abundance and state", sub=f"rows, held-out cohort; mean held-out C {mC:.3f}",
                x_mm=6.0)
    F.letter("b", t, 1.0)
    F.rec("b mean held-out C, abundance + state", f"{mC:.3f}", "revision2/results/S5_loco_cindex.csv")
    # c  arrangement score, out-of-section r
    ax = F.ax(134, 62, 14, 30)
    rs = pred.groupby("section").apply(lambda g: pearsonr(g.excl_resid, g.pred_resid)[0], include_groups=False)
    rs.index = rs.index.map(GC); rs = rs.sort_values(ascending=False)
    ax.barh(range(len(rs)), rs.values, color=STROMA, height=0.7)
    ax.set_ylim(len(rs) - 0.5, -0.6)
    ax.set_yticks(range(len(rs))); ax.set_yticklabels(rs.index, fontsize=5.0); ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False); ax.axvline(0, color=INK, lw=0.4)
    ax.set_xlim(-0.2, 0.3); ax.set_xticks([-0.2, 0, 0.2]); ax.set_xticklabels(["−0.2", "0", "0.2"])
    ax.set_xlabel("r, held-out section", labelpad=1)
    t = F.title(ax, "Arrangement score", sub=neg(f"median r {rs.median():.2f}"), x_mm=124.0)
    F.letter("c", t, 119.5)
    F.rec("c median r across sections", f"{rs.median():.3f}", "revision2/results/S8_neigh_pred.csv (Pearson per section)")
    # d  arrangement genes
    ax = F.ax(164, 62, 12, 30)
    top = pd.concat([w.nsmallest(6), w.nlargest(6)]).sort_values(ascending=False)
    ax.barh(range(len(top)), top.values, color=[STROMA if x > 0 else GREY for x in top.values], height=0.7)
    ax.set_ylim(len(top) - 0.5, -0.6)
    ax.set_yticks(range(len(top))); ax.set_yticklabels(top.index, fontsize=5.0, style="italic")
    ax.tick_params(axis="y", length=0); ax.spines["left"].set_visible(False); ax.axvline(0, color=INK, lw=0.4)
    ax.set_xticks([-0.02, 0]); ax.set_xticklabels(["−0.02", "0"])
    ax.set_xlabel("Ridge weight", labelpad=1)
    t = F.title(ax, "Genes", sub="ridge weights", x_mm=152.0)
    F.letter("d", t, 147.5)
    # e  arrangement score in the bulk cohorts
    ax = F.ax(25, 113, 42, 20)
    coh = [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254")]
    ax.set_ylim(2.6, -0.6)
    for i, (c, cn) in enumerate(coh):
        r = sas[sas.cohort == c].iloc[0]
        ax.plot(r.SAS_uni_HR, i - 0.13, "o", color=INK, ms=3.0)
        ax.plot(r.SAS_adj_HR, i + 0.13, "s", mfc="white", mec=INK, ms=2.8, mew=0.7)
        ax.text(1.03, i, neg(f"P = {r.SAS_uni_p:.2f}; ΔC {r.dC:+.3f}*"), transform=ax.get_yaxis_transform(),
                fontsize=5.0, va="center", ha="left")
        F.rec(f"e {cn}", f"HR {r.SAS_uni_HR:.2f} (adj {r.SAS_adj_HR:.2f}); P {r.SAS_uni_p:.2f}; dC {r.dC:+.3f}",
              "revision2/results/S10_prognosis.csv")
    rule(ax, 1); noy(ax)
    for i, (c, cn) in enumerate(coh):
        F.row_label(ax, i, cn)
    ax.set_xlim(0.8, 1.3); ax.set_xticks([0.8, 1.0, 1.2])
    ax.set_xlabel("HR per SD of the arrangement score", labelpad=1)
    t = F.title(ax, "Arrangement score in bulk cohorts", sub="HR per SD, unadjusted and adjusted", x_mm=6.0)
    F.letter("e", t, 1.0)
    ax.plot(0.0, -0.40, "o", color=INK, ms=2.8, transform=ax.transAxes, clip_on=False)
    ax.text(0.04, -0.40, "unadjusted", fontsize=5.0, color=SLATE, va="center", transform=ax.transAxes)
    ax.plot(0.32, -0.40, "s", mfc="white", mec=INK, ms=2.6, mew=0.7, transform=ax.transAxes, clip_on=False)
    ax.text(0.36, -0.40, "adjusted for stage and age", fontsize=5.0, color=SLATE, va="center", transform=ax.transAxes)
    ax.text(0.0, -0.51, "*ΔC over NetImmune, stage and age", fontsize=5.0, color=SLATE, va="center",
            transform=ax.transAxes)
    # f  leave-one-cohort-out C of the upgrade models
    rows = [("EN: abundance + state (30)", "(A) abundance + state", FITTED, None),
            ("EN: abundance only (10)", "(A) abundance only", FITTED, None),
            ("NetImmune (fixed)", "NetImmune (fixed)", IMMUNE, IMMUNE),
            ("EN: 6 barriernet modules", "(B) learned module weights", FITTED, None),
            ("EN: state only (20)", "(A) state only", FITTED, None)]
    ax = F.ax(134, 113, 34, 20)
    ax.set_ylim(len(rows) - 0.5, -0.6)
    for i, (k, lab, col, chipc) in enumerate(rows):
        for j, (c, mk) in enumerate([("TCGA", "o"), ("GSE84437", "^"), ("GSE62254", "v")]):
            ax.plot(loco.loc[k, c], i - 0.2 + 0.1 * j, mk, color=GREY, ms=2.2, mew=0)
        ax.plot(loco.loc[k, "mean"], i + 0.18, "D", color=col, ms=3.0, mec="white", mew=0.4)
        ax.text(1.03, i, f"{loco.loc[k, 'mean']:.3f}", transform=ax.get_yaxis_transform(), fontsize=5.0,
                va="center", ha="left")
        F.rec(f"f {lab}", f"{loco.loc[k, 'mean']:.3f}", "revision2/results/S5_loco_cindex.csv")
    noy(ax)
    for i, (k, lab, col, chipc) in enumerate(rows):
        F.row_label(ax, i, lab, col=chipc)
    ax.set_xlim(0.49, 0.65); ax.set_xticks([0.5, 0.55, 0.6, 0.65])
    ax.set_xlabel("Held-out Harrell's C", labelpad=1)
    t = F.title(ax, "Leave-one-cohort-out C", sub="gray, cohorts; diamond, mean", x_mm=104.0)
    F.letter("f", t, 99.5)
    F.save()


# ================================================================== robust
def robust():
    F = Fig("robust", 156.0)
    RUN = pd.read_csv(A / "E19_runs.csv"); SP = pd.read_csv(A / "E19_spearman_summary.csv")
    LG = pd.read_csv(A / "E19_leave_one_gene_out.csv"); PAN = pd.read_csv(A / "E19_panel_subset.csv").set_index("subset")
    AUCCI = pd.read_csv(A / "E7_auc_ci_log2.csv").set_index("Marker").loc["NetImmune"]
    FULL = RUN[RUN.type == "full"].iloc[0]
    FR = [("drop10", 10), ("drop20", 20), ("drop30", 30)]
    # a  pooled estimates against dropout
    SPEC = [
            ("AUC_NetImmune_PRJEB", "NetImmune AUC", "PRJEB25780", IMMUNE, lambda x: x,
             (AUCCI.CI_low, AUCCI.CI_high), 0.5, "AUC = 0.5", (0.42, 1.0), False, "AUC_NetImmune_PRJEB"),
            ("OR_Immune.est", "Immune log-OR", "21 checkpoint cohorts", IMMUNE, lambda x: x,
             (FULL["OR_Immune.lo"], FULL["OR_Immune.hi"]), 0.0, "0", (-0.12, 0.72), False, "OR_Immune"),
            ("OR_CAF_ECM.est", "CAF_ECM log-OR", "21 checkpoint cohorts", STROMA, lambda x: x,
             (FULL["OR_CAF_ECM.lo"], FULL["OR_CAF_ECM.hi"]), 0.1, "0.1", (-0.32, 0.3), False, "OR_CAF_ECM")]
    PX = [("panel_strict", "shared"), ("platform_HTG_OBP", "HTG"), ("platform_NanoString_PanCancerImmune", "PCI")]
    for j, (col, ttl, sub, c, tr, ci, thr, thl, yl, logy, pcol) in enumerate(SPEC):
        x0 = 22 + j * 58
        ax = F.ax(x0, 13, 34, 38)
        full = tr(FULL[col])
        ax.axhspan(ci[0], ci[1], color=LIGHT, alpha=0.55, lw=0, zorder=0)
        ax.axhline(full, color=c, lw=0.6, ls=(0, (3, 1.5)), zorder=1)
        ax.axhline(thr, color=GREY, lw=0.5, ls=(0, (1, 1.5)), zorder=1)
        ax.plot(0, full, "D", ms=3.0, color=c, mec="white", mew=0.4, zorder=3)
        for i, (t_, _) in enumerate(FR, start=1):
            for sfx, dx, filled in (("", -0.15, True), ("_list", 0.18, False)):
                v = tr(RUN.loc[RUN.type == t_ + sfx, col].values)
                lo, md, hi = np.quantile(v, [0.025, 0.5, 0.975])
                ax.plot([i + dx] * 2, [lo, hi], color=c if filled else GREY, lw=0.8, solid_capstyle="butt", zorder=2)
                ax.plot(i + dx, md, "o", ms=2.6, color=c if filled else "white", mec=c if filled else GREY,
                        mew=0.4 if filled else 0.7, zorder=3)
                if t_ == "drop30":
                    F.rec(f"a {ttl} drop30{sfx or ' (main)'}", f"median {md:.3f} ({lo:.3f}–{hi:.3f})", "E19_runs.csv")
        for i, (pp, lab) in enumerate(PX):
            x = 4.5 + i * 0.95
            v = PAN.loc[pp, pcol]
            if pd.isna(v):
                ax.text(x, yl[0] + 0.05 * (yl[1] - yl[0]) if not logy else yl[0] * 1.03, "n.e.", fontsize=5.0,
                        color=SLATE, ha="center", va="bottom", rotation=90)
            else:
                ax.plot(x, v, "s", ms=2.6, color=c, mec="white", mew=0.4, zorder=3)
                F.rec(f"a {ttl} panel {lab}", f"{v:.3f}", "E19_panel_subset.csv")
        if logy:
            ax.set_yscale("log"); ax.set_yticks([0.5, 0.6, 0.7, 0.8, 1.0])
            ax.set_yticklabels(["0.5", "0.6", "0.7", "0.8", "1"]); ax.minorticks_off()
        ax.set_ylim(*yl); ax.set_xlim(-0.6, 4.5 + 2 * 0.95 + 0.5)
        ax.set_xticks([0, 1, 2, 3] + [4.5 + i * 0.95 for i in range(3)])
        ax.set_xticklabels(["full", "10", "20", "30"] + [l for _, l in PX], fontsize=5.0)
        for tl in ax.get_xticklabels()[4:]:
            tl.set_rotation(90)
        ax.axvline(3.75, color=LIGHT, lw=0.5, zorder=0)
        ax.text(3.6, thr, thl, fontsize=5.0, color=SLATE, ha="right", va="bottom" if j != 3 else "top")
        ax.set_xlabel("Genes dropped (%)                    Panel gene set", labelpad=1, fontsize=5.0)
        t = F.title(ax, ttl, col=c, sub=sub, x_mm=x0 - 6)
        if j == 0:
            F.letter("a", t, 1.0)
    # b  leave-one-gene-out
    MODC = {"CAF_ECM": STROMA, "TAM_SUP": MACRO, "APM": IMMUNE, "IFN": IMMUNE, "TLS_B": IMMUNE, "T_CELL": IMMUNE}
    ROWS = [("HR_NetImmune", "NetImmune HR\n(benchmark 3)", FULL.HR_NetImmune if "HR_NetImmune" in FULL else
             np.exp(FULL["HR_NetImmune.est"]), "HR_k", 3, "{:.3f}"),
            ("AUC_PRJEB", "NetImmune AUC\n(PRJEB25780)", FULL.AUC_NetImmune_PRJEB, None, None, "{:.3f}"),
            ("OR_Immune", "Immune log-OR\n(21 cohorts)", FULL["OR_Immune.est"], "Immune_k", 21, "{:+.3f}"),
            ("OR_CAF_ECM", "CAF_ECM log-OR\n(21 cohorts)", FULL["OR_CAF_ECM.est"], "CAF_k", 20, "{:+.3f}")]
    n = len(LG); mods = LG.module.values
    bounds = [i for i in range(1, n) if mods[i] != mods[i - 1]]
    for r_, (col, lab, full, kcol, kfull, fmt) in enumerate(ROWS):
        ax = F.ax(26, 76 + r_ * 18, 150, 14)
        ax.axhline(full, color=INK, lw=0.5, ls=(0, (3, 1.5)), zorder=1)
        for b in bounds:
            ax.axvline(b - 0.5, ymax=0.74, color=LIGHT, lw=0.5, zorder=0)
        y = LG[col].values
        for i in range(n):
            reduced = kcol is not None and LG[kcol].iloc[i] < kfull
            cc = MODC[mods[i]]
            ax.plot([i, i], [full, y[i]], color=cc, lw=0.6, alpha=0.6, zorder=1)
            ax.plot(i, y[i], "o", ms=2.5, color="white" if reduced else cc, mec=cc, mew=0.7 if reduced else 0.3,
                    zorder=2)
        span = max(np.abs(y - full).max(), 1e-3) * 1.35
        ax.set_ylim(full - span, full + 1.9 * span); ax.set_xlim(-0.8, n - 0.2)
        ax.yaxis.set_major_locator(plt.MaxNLocator(3))
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: neg(f"{v:g}")))
        top1 = int(np.argmax(np.abs(y - full)))
        ax.text(1.0, 0.98, neg(f"full {fmt.format(full)}; range {fmt.format(y.min())} to {fmt.format(y.max())}; "
                               f"largest change {LG.gene.iloc[top1]}"),
                transform=ax.transAxes, fontsize=5.0, color=SLATE, ha="right", va="top",
                bbox=dict(fc="white", ec="none", pad=0.3))
        F.rec(f"b {lab.replace(chr(10), ' ')}", neg(f"full {fmt.format(full)}; range {fmt.format(y.min())} to "
                                                    f"{fmt.format(y.max())}; largest {LG.gene.iloc[top1]}"),
              "E19_leave_one_gene_out.csv")
        ym = F.ymm(ax, full)
        F.fig.text(4.0 / F.W, ym / F.H, lab, fontsize=5.0, ha="left", va="center", color=INK, linespacing=1.1)
        if r_ < 3:
            ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
        else:
            ax.set_xticks(range(n)); ax.set_xticklabels(LG.gene, rotation=90, fontsize=5.0, style="italic")
            ax.tick_params(axis="x", length=1.5, pad=1)
        if r_ == 0:
            starts = [0] + bounds; ends = bounds + [n]
            bb = ax.get_position()
            for s_, e_ in zip(starts, ends):
                xm = (bb.x0 + (((s_ + e_ - 1) / 2 + 0.8) / (n - 0.2 + 0.8)) * bb.width) * F.W
                ytop = bb.y1 * F.H + 2.2
                tt = F.fig.text((xm + 1.0) / F.W, ytop / F.H, mods[s_], fontsize=5.2, ha="center", va="center",
                                color=INK, fontweight="bold")
                F.chip(xm - 0.5 - 0.55 * len(mods[s_]) * 1.0 - 1.1, ytop, MODC[mods[s_]])
            t = F.text_at(4.0, 76 - 6.6, "Leave-one-gene-out", fontsize=6.0, fontweight="bold", va="baseline")
            F.letter("b", t, 1.0)
    F.save()


# ================================================================== locked
def locked():
    F = Fig("locked", 72.0)
    PRIM = "pre-registered exclusions (primary)"
    U = pd.read_csv(A / "E8_univariate.csv")
    P15 = pd.read_csv(A / "E15_per_cohort.csv"); P15 = P15[P15.model == "unadjusted"]
    PC = pd.read_csv(A / "E20_per_cohort.csv")
    PC = PC[(PC.model == "unadjusted") & (PC.endpoint == "OS") & (PC.duplicate_rule == PRIM)]
    PL = pd.read_csv(A / "E20_pooled.csv"); PL = PL[(PL.model == "unadjusted") & (PL.duplicate_rule == PRIM)]
    NU = pd.read_csv(A / "E20_null.csv").set_index("level")
    DR = pd.read_csv(S.OUT / "data/E20_null_draws.csv")
    DC = pd.read_csv(A / "E20_deltaC.csv")
    # a  structure-matched null per added cohort and pooled
    LEV = [("GSE13861", "GSE13861"), ("GSE26899", "GSE26899"), ("GSE26901", "GSE26901"), ("GSE29272", "GSE29272"),
           ("locked (added 4), pooled", "Added four, pooled")]
    for i, (k, lab) in enumerate(LEV):
        ax = F.ax(10, 13 + i * 7.6, 58, 6.2)
        v = DR[k].values
        ax.hist(v, bins=np.linspace(0, 0.8, 65), color=LIGHT, lw=0)
        r = NU.loc[k]
        ax.axvline(r.null_p95_abs, color=GREY, lw=0.6)
        ax.plot(abs(r.observed_logHR), ax.get_ylim()[1] * 0.35, "D", ms=3.0, color=IMMUNE, mec="white", mew=0.4)
        ptxt = "P ≤ 0.001" if int(r["rank"]) == 1 else fmt_p(r.P)
        ax.text(0.99, 0.95, f"{lab}  {ptxt}", transform=ax.transAxes, ha="right", va="top", fontsize=5.0,
                color=INK, fontweight="bold" if r.P < 0.05 else "normal")
        ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xlim(0, 0.8)
        if i < len(LEV) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Absolute log HR per SD", labelpad=1)
        if i == 0:
            t = F.title(ax, "NetImmune against 1,000 random axes", sub="structure-matched, added locked cohorts",
                        x_mm=6.0)
            F.letter("a", t, 1.0)
        F.rec(f"a {lab}", f"|log HR| {abs(r.observed_logHR):.3f}; 95th pct {r.null_p95_abs:.3f}; rank {int(r['rank'])}; {ptxt}",
              "E20_null.csv; draws E20_null_draws.rds")
    # b  paired difference in C
    DLEV = [("GSE13861", "GSE13861"), ("GSE26899", "GSE26899"), ("GSE26901", "GSE26901"), ("GSE29272", "GSE29272"),
            ("locked (added 4), pooled (inverse variance)", "Added four, pooled"),
            ("locked 7, pooled (inverse variance)", "Locked 7, pooled")]
    ax = F.ax(106, 13, 30, 46)
    yy, yl = 0.0, []
    for cmp in ["ESTIMATE", "CD8A"]:
        for k, lab in DLEV:
            r = DC[(DC.cohort == k) & (DC.comparator == cmp)].iloc[0]
            pooled = k.startswith("locked")
            ax.plot([r.lo, r.hi], [yy, yy], color=IMMUNE if pooled else GREY, lw=0.8)
            ax.plot(r.dC, yy, "D" if pooled else "o", ms=3.0 if pooled else 2.6, color=IMMUNE if pooled else GREY,
                    mec="white", mew=0.4)
            if pooled:
                ax.text(1.04, yy, neg(f"{r.dC:+.3f} ({r.lo:+.3f} to {r.hi:+.3f})"), transform=ax.get_yaxis_transform(),
                        fontsize=5.0, va="center", ha="left", fontweight="bold" if r.p < 0.05 else "normal")
                ax.text(2.4, yy, fmt_p(r.p, ""), transform=ax.get_yaxis_transform(), fontsize=5.0, va="center",
                        ha="right", fontweight="bold" if r.p < 0.05 else "normal")
            F.rec(f"b vs {cmp} {lab}", f"{r.dC:+.3f} ({r.lo:+.3f} to {r.hi:+.3f}), P {fmt_p(r.p, '')}", "E20_deltaC.csv")
            yl.append((yy, lab, pooled, cmp)); yy += 1.0
        yy += 1.5
    ax.set_ylim(yy - 0.5, -1.4)
    rule(ax, 0); noy(ax)
    for v, lab, pooled, cmp in yl:
        F.row_label(ax, v, lab, bold=pooled)
    for cmp, v in [("vs ESTIMATE", -0.95), ("vs CD8A", 6.55)]:
        F.row_label(ax, v, cmp, bold=True, color=SLATE)
    ax.set_xlim(-0.12, 0.3); ax.set_xticks([-0.1, 0, 0.1, 0.2, 0.3])
    ax.set_xticklabels(["−0.1", "0", "0.1", "0.2", "0.3"])
    ax.set_xlabel("ΔC, NetImmune − comparator", labelpad=1)
    ax.text(2.4, 1.0, "P", transform=ax.transAxes, fontsize=5.0, color=SLATE, va="bottom", ha="right")
    t = F.title(ax, "Harrell's C, head to head", sub="paired bootstrap; GSE26253 recurrence-free", x_mm=84.0)
    F.letter("b", t, 79.5)
    F.save()


# ================================================================== prjeb_detail
def prjeb_detail():
    from sklearn.metrics import roc_curve, roc_auc_score
    F = Fig("prjeb_detail", 128.0)
    Sx = pd.read_csv(A / "E7_scores_log2.csv"); y = Sx.Responder.values.astype(int)
    # a  ROC
    ax = F.ax(14, 13, 34, 34)
    curves = [("NetImmune", Sx.NetImmune, IMMUNE, "-", 1.1), ("IFN module", Sx.IFN, IMMUNE, (0, (2.5, 1.2)), 0.8),
              ("GEP-18", Sx.GEP18, INK, (0, (3, 1.5)), 0.8), ("PD-L1 (CD274)", Sx.PDL1, GREY, (0, (1, 1.2)), 0.9),
              ("IPS", Sx.IPS_AZ, GREY, (0, (4, 1.5, 1, 1.5)), 0.8)]
    for lab, sc, col, ls, lw in curves:
        f, t_, _ = roc_curve(y, sc)
        auc = roc_auc_score(y, sc)
        ax.plot(f, t_, color=col, ls=ls, lw=lw, drawstyle="steps-post", label=f"{lab}  {auc:.2f}")
        F.rec(f"a AUC {lab}", f"{auc:.3f}", "recomputed from E7_scores_log2.csv; matches E7_auc_ci_log2")
    ax.plot([0, 1], [0, 1], color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.set_xlabel("1 − specificity", labelpad=1); ax.set_ylabel("Sensitivity", labelpad=1)
    lg = ax.legend(loc="lower left", bbox_to_anchor=(1.02, -0.02), fontsize=5.0, handlelength=1.6, title="AUC",
                   title_fontsize=5.0, labelspacing=0.35, borderaxespad=0)
    lg._legend_box.align = "left"
    t = F.title(ax, "ROC", sub=f"{int(y.sum())} responders of {len(y)}", x_mm=6.0)
    F.letter("a", t, 1.0)
    # b  structure-matched null
    d = pd.read_csv(A / "E3_ici_null_draws_log2.csv")
    sn = pd.read_csv(A / "E3_ici_matched_null.csv"); sn = sn[sn.scale == "log2"].set_index("score")
    ax = F.ax(80, 13, 38, 30)
    items = [("Immune", "Immune only (26 genes)"), ("IFN", "IFN module (5 genes)")]
    bins = np.linspace(0.5, 0.95, 46)
    for i, (k, lab) in enumerate(items):
        yb = len(items) - 1 - i
        h, e = np.histogram(d[k], bins=bins); h = h / h.max() * 0.5
        ax.fill_between(e[:-1], yb, yb + h, step="post", color=LIGHT, lw=0)
        q95 = sn.loc[k, "null_p95"]
        ax.plot([q95, q95], [yb, yb + 0.55], color=GREY, lw=0.6)
        ax.plot(sn.loc[k, "AUC_obs"], yb + 0.15, "D", color=IMMUNE, ms=3.0, mec="white", mew=0.4)
        pv = sn.loc[k, "p_emp"]
        ptxt = "P ≤ 0.0005" if pv <= 0.0005 else fmt_p(pv)
        ax.text(0.505, yb + 0.85, lab, fontsize=5.2, fontweight="bold", va="center")
        ax.text(0.945, yb + 0.66, f"AUC {sn.loc[k, 'AUC_obs']:.3f}, {ptxt}", fontsize=5.0, ha="right", va="center")
        F.rec(f"b {lab}", f"AUC {sn.loc[k, 'AUC_obs']:.3f}; 95th pct {q95:.3f}; {ptxt}", "E3_ici_matched_null.csv (log2)")
    ax.set_ylim(-0.1, len(items)); noy(ax)
    ax.set_xlim(0.5, 0.95); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_xlabel("AUC (random axes: better direction)", labelpad=1)
    t = F.title(ax, "Structure-matched random-gene null", sub="2,000 random axes", x_mm=72.0)
    F.letter("b", t, 67.5)
    # c  every marker against NetImmune
    a = pd.read_csv(A / "E7_auc_ci_log2.csv").sort_values("AUC", ascending=False)
    dl = pd.read_csv(A / "E7_delong_log2.csv").set_index("Comparator")
    names = {"GEP18": ("GEP-18", None), "IFN": ("IFN module", IMMUNE), "APM": ("APM module", IMMUNE),
             "TAM_SUP": ("TAM_SUP module*", MACRO), "T_CELL": ("T_CELL module", IMMUNE),
             "TLS_B": ("TLS_B module", IMMUNE), "Immune": ("Immune only", IMMUNE), "NetImmune": ("NetImmune", IMMUNE),
             "Barrier(neg)": ("Barrier (inverted)", STROMA), "MSI/EBV": ("MSI/EBV status", None),
             "MSI-H": ("MSI-H status", None), "PD-L1": ("PD-L1 (CD274)", None), "IPS": ("IPS", None),
             "ExclusionProxy(neg)": ("Exclusion proxy†", None)}
    ax = F.ax(146, 13, 22, 100)
    ax.set_ylim(len(a) - 0.5, -0.6)
    for i, r in enumerate(a.itertuples()):
        lab, col = names.get(r.Marker, (r.Marker, None))
        mc = col or INK
        ax.plot([r.CI_low, r.CI_high], [i, i], color=mc, lw=0.8)
        ax.plot(r.AUC, i, "o", color=mc, ms=2.8, mec="white", mew=0.4)
        if r.Marker in dl.index and dl.loc[r.Marker, "q"] < 0.05:
            q = dl.loc[r.Marker, "q"]
            ax.text(1.03, i, f"q = {q:.3f}", transform=ax.get_yaxis_transform(), fontsize=5.0, va="center")
            F.rec(f"c DeLong q {lab}", f"{q:.3f}", "E7_delong_log2.csv")
        F.rec(f"c AUC {lab}", f"{r.AUC:.3f} ({r.CI_low:.3f}–{r.CI_high:.3f})", "E7_auc_ci_log2.csv")
    rule(ax, 0.5); noy(ax)
    for i, r in enumerate(a.itertuples()):
        lab, col = names.get(r.Marker, (r.Marker, None))
        F.row_label(ax, i, lab, col=col, fs=5.0)
    ax.set_xlim(0.2, 1.02); ax.set_xticks([0.25, 0.5, 0.75, 1]); ax.set_xticklabels(["0.25", "0.5", "0.75", "1"])
    ax.set_xlabel("AUC (95% CI)", labelpad=1)
    t = F.title(ax, "Every marker", sub="DeLong q against NetImmune", x_mm=124.0)
    F.letter("c", t, 119.5)
    ax.text(-1.05, -0.06, "*higher TAM_SUP predicts response\n†CAF/TAM exclusion proxy, 16 genes, inverted",
            transform=ax.transAxes, fontsize=5.0, color=SLATE, va="top", linespacing=1.2)
    # d  same cross-validation pipeline
    cv = pd.read_csv(A / "E7_samepipeline_cv.csv")
    best = cv.loc[cv.groupby("score").cv_auc.idxmax()].set_index("score")
    ce = pd.read_csv(A / "E3_ici_ceiling.csv").sort_values("cv_auc").iloc[-1]
    rows = [("IFN module", best.loc["IFN"], IMMUNE, False), ("Immune only", best.loc["Immune"], IMMUNE, False),
            ("NetImmune", best.loc["NetImmune"], IMMUNE, False), ("Full transcriptome\n(fitted)", ce, None, True)]
    ax = F.ax(32, 74, 58, 34)
    ax.set_ylim(len(rows) - 0.5, -0.6)
    for i, (lab, r, chipc, fitted) in enumerate(rows):
        col = FITTED if fitted else IMMUNE
        ax.plot([r.lo, r.hi], [i, i], color=col, lw=0.9)
        ax.plot(r.cv_auc, i, "D" if fitted else "o", color=col, ms=3.0, mec="white", mew=0.4)
        ax.text(1.03, i, f"{r.cv_auc:.3f}", transform=ax.get_yaxis_transform(), fontsize=5.0, va="center")
        F.rec(f"d {lab.replace(chr(10), ' ')}", f"{r.cv_auc:.3f} ({r.lo:.3f}–{r.hi:.3f})",
              "E7_samepipeline_cv.csv / E3_ici_ceiling.csv")
    noy(ax)
    for i, (lab, r, chipc, fitted) in enumerate(rows):
        F.row_label(ax, i, lab, col=chipc, fs=5.0)
    ax.set_xlim(0.64, 0.88); ax.set_xticks([0.65, 0.75, 0.85])
    ax.set_xlabel("Cross-validated AUC", labelpad=1)
    t = F.text_at(6.0, 74 - 10.5, "Same CV pipeline", fontsize=6.0, fontweight="bold", va="baseline")
    F.letter("d", t, 1.0)
    F.save()


# ================================================================== gastric_icb
def gastric_icb():
    F = Fig("gastric_icb", 84.0)
    R1 = pd.read_csv(A / "E22_PRJEB60680_response.csv"); P = pd.read_csv(A / "E22_pooled.csv")
    FC = pd.read_csv(A / "E22_firth_cps.csv")
    PRIM = "PRJEB60680 + PRJEB25780 (primary)"; SENS = "PRJEB60680 + PRJEB25780 + PERFECT (sensitivity)"
    ALL = "primary: all QC-pass baseline, z within library batch"; MSS = "MSS/EBV-negative (pre-registered subgroup)"
    ROWS = [("CAF_ECM", STROMA), ("Barrier", STROMA), ("TAM_SUP", MACRO), ("IFN", IMMUNE), ("APM", IMMUNE),
            ("TLS_B", IMMUNE), ("T_CELL", IMMUNE), ("Immune", IMMUNE), ("NetImmune", IMMUNE)]
    NR = len(ROWS)

    def frame(ax, lim, ticks):
        ax.set_ylim(NR - 0.45, -0.6)
        ax.set_xlim(*lim); ax.set_xticks(ticks); ax.set_xticklabels([neg(f"{t:g}") for t in ticks])
        rule(ax, 0.0); noy(ax)
        for i in (2.5, 6.5):
            ax.axhline(i, color=LIGHT, lw=0.4, zorder=0)
        for i, (k, c) in enumerate(ROWS):
            F.row_label(ax, i, k, col=c)

    # a  pooled over the gastro-esophageal cohorts, without and with PERFECT
    LIM = (-2.1, 3.1)
    ax = F.ax(24, 14, 40, 48)
    P1 = P[P.set == PRIM].set_index("score"); P2 = P[P.set == SENS].set_index("score")
    frame(ax, LIM, [-2, -1, 0, 1, 2, 3])
    tr = ax.get_yaxis_transform()
    for i, (k, c) in enumerate(ROWS):
        a_, b_ = P1.loc[k], P2.loc[k]
        point(ax, a_.est, a_.lo, a_.hi, i - 0.15, c, mk="D", ms=2.8, lim=LIM)
        point(ax, b_.est, b_.lo, b_.hi, i + 0.17, c, mk="D", ms=2.4, open_=True, lw=0.6, lim=LIM)
        ax.text(1.15, i, fmt_p(a_.p, ""), transform=tr, fontsize=5.0, ha="right", va="center", clip_on=False)
        ax.text(1.47, i, fmt_p(b_.p, ""), transform=tr, fontsize=5.0, ha="right", va="center", color=SLATE,
                clip_on=False)
        F.rec(f"a pooled {k}", f"{a_.est:+.2f} ({a_.lo:+.2f} to {a_.hi:+.2f}), P {fmt_p(a_.p, '')}; + PERFECT "
              f"{b_.est:+.2f} ({b_.lo:+.2f} to {b_.hi:+.2f}), P {fmt_p(b_.p, '')}", "E22_pooled.csv")
    ax.text(1.15, 1.0, "P", transform=ax.transAxes, fontsize=5.0, color=SLATE, va="bottom", ha="right")
    ax.text(1.47, 1.0, "+PERFECT", transform=ax.transAxes, fontsize=5.0, color=SLATE, va="bottom", ha="right")
    ax.set_xlabel("log-odds of response per SD", labelpad=1)
    t = F.title(ax, "Pooled gastroesophageal response", sub="random effects, conservative interval", x_mm=6.0)
    F.letter("a", t, 1.0)
    ax.plot(0.0, -0.21, "D", ms=2.8, color=INK, transform=ax.transAxes, clip_on=False)
    ax.text(0.05, -0.21, "NCT04249739 + PRJEB25780", fontsize=5.0, color=SLATE, va="center", transform=ax.transAxes)
    ax.plot(0.0, -0.275, "D", ms=2.4, mfc="white", mec=INK, mew=0.8, transform=ax.transAxes, clip_on=False)
    ax.text(0.05, -0.275, "+ PERFECT (sensitivity)", fontsize=5.0, color=SLATE, va="center", transform=ax.transAxes)
    # b  MSS/EBV-negative patients of NCT04249739
    ax = F.ax(124, 14, 34, 48)
    T1 = R1[R1.analysis == MSS].set_index("score"); T0 = R1[R1.analysis == ALL].set_index("score")
    lim = (-1.6, 2.6)
    frame(ax, lim, [-1, 0, 1, 2])
    for i, (k, c) in enumerate(ROWS):
        a_, b_ = T1.loc[k], T0.loc[k]
        point(ax, a_.logOR, a_.lo, a_.hi, i - 0.15, c, lim=lim)
        point(ax, b_.logOR, b_.lo, b_.hi, i + 0.17, c, ms=2.4, open_=True, lw=0.6, lim=lim)
        ax.text(1.04, i, fmt_p(a_.p, ""), transform=ax.get_yaxis_transform(), fontsize=5.0, ha="left", va="center")
        F.rec(f"b MSS/EBV-negative {k}", f"{a_.logOR:+.2f} ({a_.lo:+.2f} to {a_.hi:+.2f}), P {fmt_p(a_.p, '')}",
              "E22_PRJEB60680_response.csv")
    ax.text(1.04, 1.0, "P", transform=ax.transAxes, fontsize=5.0, color=SLATE, va="bottom")
    ax.set_xlabel("log-odds of response per SD", labelpad=1)
    nn = T1.loc["NetImmune"]
    t = F.title(ax, "NCT04249739, MSS/EBV-negative", sub=f"{int(nn.n)} patients, {int(nn.n_resp)} CR/PR", x_mm=104.0)
    F.letter("b", t, 99.5)
    f_ = FC[(FC.set == MSS) & (FC.term == "NetImmune per SD") & FC.analysis.str.contains("primary")].iloc[0]
    ax.plot(-0.55, -0.21, "o", ms=2.6, color=INK, transform=ax.transAxes, clip_on=False)
    ax.text(-0.50, -0.21, "MSS/EBV-negative", fontsize=5.0, color=SLATE, va="center", transform=ax.transAxes)
    ax.plot(0.32, -0.21, "o", ms=2.4, mfc="white", mec=INK, mew=0.7, transform=ax.transAxes, clip_on=False)
    ax.text(0.37, -0.21, "all patients", fontsize=5.0, color=SLATE, va="center", transform=ax.transAxes)
    F.rec("b Firth NetImmune beyond log(CPS+1), MSS/EBV-negative", f"{f_.estimate:+.2f}; P {fmt_p(f_.p, '')}",
          "E22_firth_cps.csv")
    F.save()

FIGS = {"meta_forest": meta_forest, "fitted_models": fitted_models, "upgrades": upgrades, "robust": robust,
        "locked": locked, "prjeb_detail": prjeb_detail, "gastric_icb": gastric_icb}

if __name__ == "__main__":
    keys = sys.argv[1:] or list(FIGS)
    for k in keys:
        FIGS[k]()
