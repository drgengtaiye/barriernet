"""v19 supplementary figures, pan-cancer and trial block, redrawn in the cns house style.

    python make_supp_pancancer.py <key> [<key> ...]      (no key: all)

keys: icb21_forest, icb11, tcga, perfect_metabric, trials_detail, trials_forest, proteomics, drugs.
Output: figures/cns/supp/FigS_<key>.pdf/.png. Every estimate is read from analysis/results (or the
tables the earlier figure scripts read); nothing here fits a model or runs a test, except the
descriptive Spearman correlations printed on panels that already carried them (tcga b, drugs d-e).
No supplementary-figure numbers are drawn inside the images.
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, FancyBboxPatch
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.colors import TwoSlopeNorm, to_rgb
from scipy import stats
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, PAPER, IMMUNE, STROMA, MACRO, fmt_p

R = S.RES
OUTD = S.OUT / "supp"
PT = 25.4 / 72
MINUS = "−"
NUMBERS = {}                        # key -> list of (panel, text) for FigS_<key>_numbers.md


# ------------------------------------------------------------------------------------------
# shared helpers
# ------------------------------------------------------------------------------------------
def fmtp(p, lead=""):
    """Two significant figures, keeping a trailing zero (0.060); a x 10^-b below 0.001."""
    if p is None or not np.isfinite(p):
        return lead + "n.a."
    if 0.001 <= p < 0.1:
        return lead + f"{p:#.2g}"
    return fmt_p(p, lead=lead)


def sg(v, nd=2):
    return f"{v:+.{nd}f}".replace("-", MINUS)


def mn(s):
    return str(s).replace("-", MINUS)


DISP = {"CXCL9_SPP1": "CXCL9:SPP1", "C1Q_arm": "C1Q subset", "SPP1_arm": "SPP1 subset", "Immune": "Immune"}


def disp(k):
    return DISP.get(k, k)


def colof(k):
    return S.MODCOL.get(disp(k), MACRO if k == "CXCL9_SPP1" else IMMUNE)


def mrk(k):
    return "s" if k in ("Barrier", "Immune", "NetImmune") else "o"


class Fig:
    def __init__(self, h, w=180):
        self.W, self.H = w, h
        self.fig = plt.figure(figsize=(w * S.MM, h * S.MM))
        self.cv = self.fig.add_axes([0, 0, 1, 1]); self.cv.set_xlim(0, w); self.cv.set_ylim(h, 0); self.cv.axis("off")

    def txt(self, x, y, s, size=5.0, color=INK, ha="left", va="center", **kw):
        return self.cv.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, **kw)

    @staticmethod
    def _ink(s, size):
        return TextPath((0, 0), s, size=size, prop=FontProperties(family=S.fs.SANS, weight="bold")).get_extents().y1 * PT

    def heading(self, xl, xt, ytop, letter, title, sub=None):
        """Letter (8 pt) and bold title (6.5 pt) with inked tops level at ytop; subtitles go to the legend."""
        sub = None
        if letter:
            self.fig.text(xl / self.W, 1 - (ytop + self._ink(letter, 8)) / self.H, letter, fontsize=8, fontweight="bold",
                          va="baseline", ha="left", color=INK)
        if title:
            self.cv.text(xt, ytop + self._ink(title, 6.5), title, fontsize=6.5, fontweight="bold", va="baseline",
                         ha="left", color=INK)
        if sub:
            self.txt(xt, ytop + 4.2, sub, size=5.2, color=SLATE)

    def chip(self, x, y, col, s=1.5):
        self.cv.add_patch(Rectangle((x, y - s / 2), s, s, color=col, lw=0))

    def label(self, x, y, k, text=None, size=5.0, ha="right", bold=False):
        """Ink label with a family-colour chip on its left (ha='right': chip left of the text block)."""
        t = text if text is not None else disp(k)
        if ha == "right":
            w = TextPath((0, 0), t, size=size, prop=FontProperties(family=S.fs.SANS,
                         weight="bold" if bold else "normal")).get_extents().width * PT
            self.chip(x - w - 2.4, y, colof(k))
        else:
            self.chip(x, y, colof(k)); x = x + 2.4
        self.txt(x, y, t, size=size, ha=ha, fontweight="bold" if bold else "normal")

    def ax(self, x, y, w, h, rows=None):
        a = S.ax_mm(self.fig, x, y, w, h)
        a.patch.set_alpha(0)                      # let canvas text and washes show through
        if rows is not None:
            rows_axis(a, rows)
        return a

    def coltitle(self, ax, title, sub=None, size=5.4):
        """Bold column title a fixed distance above an axes; subtitles go to the legend."""
        sub = None
        top = self.H * (1 - ax.get_position().y1)
        x = self.W * ax.get_position().x0
        self.txt(x, top - (4.4 if sub else 1.2), title, size=size, fontweight="bold", va="bottom", linespacing=1.1)
        if sub:
            self.txt(x, top - 2.2, sub, size=5.0, color=SLATE)

    def ymm(self, ax, yv):
        """Figure-mm y of a data y in axes ax."""
        bb = ax.get_position(); y0, y1 = ax.get_ylim()
        f = (yv - y0) / (y1 - y0)
        return self.H * (1 - (bb.y0 + f * bb.height))

    def save(self, key):
        OUTD.mkdir(parents=True, exist_ok=True)
        self.fig.savefig(OUTD / f"FigS_{key}.pdf")
        self.fig.savefig(OUTD / f"FigS_{key}.png", dpi=600)
        plt.close(self.fig)
        print("saved", OUTD / f"FigS_{key}.png")


def rows_axis(ax, n, pad=0.55):
    ax.set_ylim(-(n - 1) - pad, pad); ax.set_yticks([]); ax.spines["left"].set_visible(False)


def xlog(ax, ticks, lim):
    ax.set_xscale("log"); ax.set_xlim(*lim); ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks]); ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.xaxis.set_minor_locator(plt.NullLocator())


def xlin(ax, ticks, lim, fmt="{:g}"):
    ax.set_xlim(*lim); ax.set_xticks(ticks); ax.set_xticklabels([mn(fmt.format(t)) for t in ticks])


def ci(ax, lo, hi, y, col, lw=0.7, lim=None):
    """Horizontal CI, clipped to lim with an arrowhead where it runs off."""
    a, b = lo, hi
    if lim is not None:
        a, b = max(lo, lim[0]), min(hi, lim[1])
    ax.plot([a, b], [y, y], color=col, lw=lw, solid_capstyle="butt", zorder=2)
    if lim is not None:
        if lo < lim[0]:
            ax.plot(lim[0], y, "<", color=col, ms=2.0, mew=0, clip_on=False, zorder=3)
        if hi > lim[1]:
            ax.plot(lim[1], y, ">", color=col, ms=2.0, mew=0, clip_on=False, zorder=3)


def diamond(ax, est, lo, hi, y, col, h=0.30, filled=True, lim=None):
    if lim is not None:
        lo, hi = max(lo, lim[0]), min(hi, lim[1])
    ax.add_patch(Polygon([[lo, y], [est, y + h], [hi, y], [est, y - h]], closed=True,
                         fc=col if filled else "white", ec=col, lw=0.6, zorder=3))


def null_hist(ax, draws, xmax, bins=36):
    a = np.abs(np.asarray(draws, float))
    ax.hist(a, bins=np.linspace(0, xmax, bins), color=LIGHT, lw=0)
    top = ax.get_ylim()[1]
    p95 = np.quantile(a, 0.95)
    ax.plot([p95, p95], [0, top * 0.92], color=GREY, lw=0.6)
    ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xlim(0, xmax)
    return top, p95


R5 = pd.read_csv(R / "E12_icb_per_cohort.csv")
NAME = {"ICB_Braun": "Braun", "ICB_Gide": "Gide", "ICB_Hugo": "Hugo", "ICB_Jung": "Jung", "ICB_Limagne1": "Limagne 1",
        "ICB_Liu": "Liu", "ICB_Mariathasan": "Mariathasan", "ICB_Miao1": "Miao", "ICB_Puch": "Puch",
        "ICB_Riaz": "Riaz", "ICB_Van_Allen": "Van Allen", "ICB_Nathanson": "Nathanson", "ICB_Snyder": "Snyder",
        "ICB_Limagne2": "Limagne 2", "ICB_Shiuan": "Shiuan", "GSE176307": "GSE176307", "GSE161537": "GSE161537",
        "GSE159067": "GSE159067", "GSE215868": "GSE215868", "GSE93157": "GSE93157",
        "brca_iatlas_anders_2022": "Anders"}
TUM = {"kidney": "renal", "melanoma": "melanoma", "lung": "NSCLC", "bladder": "urothelial", "urothelial": "urothelial",
       "head and neck": "HNSCC", "mixed (lung, HNSCC, melanoma)": "mixed", "breast": "breast"}
PANEL = {"GSE161537", "GSE159067", "GSE215868", "GSE93157"}
TORD = ["melanoma", "NSCLC", "urothelial", "renal", "HNSCC", "breast", "mixed"]


# ==========================================================================================
# icb21_forest
# ==========================================================================================
def icb21_forest():
    key = "icb21_forest"
    per = R5[R5.set.isin(["E5", "new"])]
    pool = pd.read_csv(R / "E12_icb_pooled.csv")
    meta = per.drop_duplicates("cohort").set_index("cohort"); meta["tum"] = meta.tumour.map(TUM)

    def order(s):
        m = meta[meta.set == s].copy(); m["t"] = m.tum.map(TORD.index)
        return m.sort_values(["t", "n"], ascending=[True, False]).index.tolist()

    SC = ["CAF_ECM", "TAM_SUP", "IFN", "Immune", "NetImmune"]
    POOLS = [("E5 11 (reference)", "Pooled, first 11"), ("new only (replication)", "Pooled, replication 10"),
             ("all", "Pooled, all 21")]
    rows = [("head", "First 11 cohorts")] + [("E5", c) for c in order("E5")] + [("gap", None)] + \
           [("head", "Replication cohorts")] + [("new", c) for c in order("new")] + [("gap", None)] + \
           [("pool", p) for p in POOLS]
    H_ROW = {"head": 3.4, "E5": 3.1, "new": 3.1, "gap": 1.2, "pool": 10.0}
    Y1_ = 16.0 + sum(H_ROW[k_] for k_, _ in rows)
    F = Fig(h=Y1_ + 18.0)
    YT = 1.7
    F.heading(1.0, 5.0, YT, "", "Checkpoint response in 21 cohorts, per cohort",
              "log-odds of response per SD; circles, first 11 cohorts; squares, replication cohorts")
    Y0 = 16.0                                     # top of rows
    ys, y = [], Y0
    for kind, v in rows:
        ys.append(y + H_ROW[kind] / 2); y += H_ROW[kind]
    Y1 = y
    XL0, CWID, GAP = 49.0, 23.0, 2.6
    LIM = (-1.6, 2.6)
    # column headers
    F.txt(4.2, Y0 - 1.6, "cohort", size=5.0, color=GREY); F.txt(22.0, Y0 - 1.6, "tumor", size=5.0, color=GREY)
    F.txt(45.5, Y0 - 1.6, "n / resp.", size=5.0, color=GREY, ha="right")
    nums = []
    for j, k in enumerate(SC):
        x0 = XL0 + j * (CWID + GAP)
        ax = F.ax(x0, Y0, CWID, Y1 - Y0)
        ax.set_ylim(Y1, Y0)                       # data y in figure mm
        ax.set_yticks([]); ax.spines["left"].set_visible(False)
        col = colof(k)
        for (kind, v), yy in zip(rows, ys):
            if kind == "new":
                ax.add_patch(Rectangle((LIM[0], yy - 1.55), LIM[1] - LIM[0], 3.1, color=col, alpha=0.07, lw=0, zorder=0))
            if kind in ("E5", "new"):
                r = per[(per.cohort == v) & (per.score == k)]
                if r.empty or not np.isfinite(r.iloc[0].logOR):
                    ax.text(-0.35, yy, "n.c.", fontsize=5.0, color=GREY, ha="right", va="center"); continue
                r = r.iloc[0]; lo, hi = r.logOR - 1.96 * r.se, r.logOR + 1.96 * r.se
                ci(ax, lo, hi, yy, col, lw=0.6, lim=LIM)
                ms = float(np.clip(1.25 / r.se, 1.8, 4.2))
                ax.plot(r.logOR, yy, "s" if kind == "new" else "o", color=col, ms=ms, mec="white", mew=0.3, zorder=3)
            elif kind == "pool":
                q = pool[(pool.set == v[0]) & (pool.score == k)].iloc[0]
                pmax = max(q.p_DL, q.p_HK)
                diamond(ax, q.logOR, q.lo_DL, q.hi_DL, yy - 3.0, col, h=1.0, lim=LIM)
                ax.plot([q.lo_HK, q.hi_HK], [yy - 1.3, yy - 1.3], color=col, lw=0.5, zorder=2)
                ax.text(0.5, yy + 1.1, f"{sg(q.logOR)}   {int(q.n_pos)}/{int(q.k)}", fontsize=5.0, ha="center",
                        va="center", fontweight="bold" if pmax < 0.05 else "normal", transform=ax.get_yaxis_transform(),
                        bbox=dict(fc="white", ec="none", pad=0.4), zorder=5)
                ax.text(0.5, yy + 3.4, fmtp(pmax, "P = "), fontsize=5.0, ha="center", va="center", color=SLATE,
                        transform=ax.get_yaxis_transform(), bbox=dict(fc="white", ec="none", pad=0.4), zorder=5)
                nums.append(f"{k} {v[1]}: {sg(q.logOR)} (DL {mn(f'{q.lo_DL:.2f}')} to {mn(f'{q.hi_DL:.2f}')}; "
                            f"HK {mn(f'{q.lo_HK:.2f}')} to {mn(f'{q.hi_HK:.2f}')}), {int(q.n_pos)}/{int(q.k)}, P = {fmtp(pmax)}")
        S.rule(ax, 0)
        xlin(ax, [-1, 0, 1, 2], LIM)
        F.label(x0 + 0.5, Y0 - 1.6, k, ha="left", size=5.6, bold=True)
    for (kind, v), yy in zip(rows, ys):
        if kind == "head":
            F.txt(2.0, yy, v, size=5.4, fontweight="bold")
        elif kind in ("E5", "new"):
            m = meta.loc[v]
            F.txt(4.2, yy, NAME[v] + ("†" if v in PANEL else ""), size=5.0)
            F.txt(22.0, yy, m.tum, size=5.0, color=SLATE)
            F.txt(45.5, yy, f"{int(m.n)}/{int(m.n_resp)}", size=5.0, color=SLATE, ha="right")
        elif kind == "pool":
            F.txt(4.2, yy - 2.2, v[1], size=5.2, fontweight="bold")
    F.txt(XL0 + (5 * CWID + 4 * GAP) / 2, Y1 + 8.0, "log-odds of response per SD", size=6.0, ha="center")
    F.save(key)
    NUMBERS[key] = [("all", "Per-cohort log-OR ± 1.96 SE: `E12_icb_per_cohort.csv` (set E5, new). Pooled: `E12_icb_pooled.csv` "
                     "sets `E5 11 (reference)`, `new only (replication)`, `all`; diamond = DL (`lo_DL`, `hi_DL`), bar = HK "
                     "(`lo_HK`, `hi_HK`), P = max(`p_DL`, `p_HK`).")] + [("pooled", s) for s in nums]


# ==========================================================================================
# icb11
# ==========================================================================================
def icb11():
    key = "icb11"
    T = pd.read_csv(R / "E5_trials.csv"); T = T[T.trial.str.startswith("IMmotion150: atezolizumab (+/- bev)")].set_index("score")
    pool = pd.read_csv(R / "E5_icb_pooled.csv"); pool = pool[pool.set == "all 11"].set_index("score")
    G = pd.read_csv(R / "E5_gastric_arms_pooled.csv").set_index("score")
    K = pd.read_csv(R / "E5_prjeb_arms.csv"); K = K[K.scale == "log2"].set_index("score")
    ROWS = ["CAF_ECM", "TAM_SUP", "Barrier", "IFN", "Immune", "NetImmune"]
    F = Fig(h=88)
    nums = []
    YA = 1.7; YAX = YA + 13.0; HAX = 54.0
    # ---- a: IMmotion150 within arm -------------------------------------------------------------
    F.heading(1.0, 5.0, YA, "a", "IMmotion150, within arm",
              "renal; PFS HR per SD in each arm and the interaction HR")
    axa = F.ax(21.5, YAX, 46.0, HAX, rows=len(ROWS))
    LC = (0.5, 1.9)
    for i, k in enumerate(ROWS):
        y = -i; r = T.loc[k]; col = colof(k)
        ci(axa, r.ici_lo, r.ici_hi, y + 0.16, col, lw=0.8, lim=LC)
        axa.scatter([r.ici_HR], [y + 0.16], s=12, color=col, ec="white", lw=0.4, zorder=3)
        ci(axa, r.ctrl_lo, r.ctrl_hi, y - 0.16, GREY, lw=0.7, lim=LC)
        axa.scatter([r.ctrl_HR], [y - 0.16], s=11, marker="s", color="white", ec=GREY, lw=0.7, zorder=3)
        F.label(20.0, F.ymm(axa, y), k)
        sig = r.int_p < 0.05
        F.txt(77.0, F.ymm(axa, y), f"{r.int_HR:.2f}", ha="right", fontweight="bold" if sig else "normal")
        F.txt(86.0, F.ymm(axa, y), fmtp(r.int_p), ha="right")
        nums.append(("a", f"{k}: atezolizumab ± bev HR {r.ici_HR:.2f}, sunitinib {r.ctrl_HR:.2f}, interaction HR {r.int_HR:.2f}, P = {fmtp(r.int_p)}"))
    axa.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
    xlog(axa, [0.5, 0.7, 1, 1.4], LC)
    F.txt(77.0, YAX - 2.0, "int. HR", color=GREY, ha="right"); F.txt(86.0, YAX - 2.0, "P", color=GREY, ha="right")
    F.txt(44.5, YAX + HAX + 7.0, "PFS HR per SD, within arm", size=6.0, ha="center")
    ky = YAX + HAX + 12.0
    F.cv.scatter([4.0], [ky], s=11, color=INK, ec="white", lw=0.4); F.txt(5.6, ky, "atezolizumab ± bevacizumab (174)")
    F.cv.scatter([52.0], [ky], s=11, marker="s", color="white", ec=GREY, lw=0.7); F.txt(53.6, ky, "sunitinib (89)")
    # ---- b: macrophage arms and polarity in four settings --------------------------------------
    X0 = 92.0
    F.heading(X0 - 3.0, X0 + 1.0, YA, "b", "Macrophage arms and polarity",
              "TAM_SUP, its C1Q and SPP1 subsets, CXCL9:SPP1; 95% CIs")
    AR = ["TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1"]
    MK = {"TAM_SUP": "o", "C1Q_arm": "^", "SPP1_arm": "v", "CXCL9_SPP1": "D"}
    spec = [("Gastric prognosis", "HR, 3 cohorts", lambda k: (G.loc[k, "HR"], G.loc[k, "lo"], G.loc[k, "hi"]), 1, (0.7, 1.4), True, [0.7, 1, 1.4]),
            ("Gastric pembrolizumab", "AUC, PRJEB25780", lambda k: (K.loc[k, "AUC"], K.loc[k, "lo"], K.loc[k, "hi"]), 0.5, (0.3, 1.0), False, [0.5, 0.75, 1.0]),
            ("First 11 cohorts", "log-odds of response", lambda k: (pool.loc[k, "logOR"], pool.loc[k, "lo"], pool.loc[k, "hi"]), 0, (-0.3, 0.8), False, [0, 0.5]),
            ("IMmotion150", "interaction HR", lambda k: (T.loc[k, "int_HR"], T.loc[k, "int_lo"], T.loc[k, "int_hi"]), 1, (0.45, 1.6), True, [0.5, 0.7, 1, 1.4])]
    wd, hd, gx, gy = 30.0, 20.0, 6.0, 15.0
    for j, (ttl, sub, get, ref, lim, lg, tk) in enumerate(spec):
        cx, cy = j % 2, j // 2
        ax = F.ax(X0 + 20.5 + cx * (wd + gx), YAX + 4.0 + cy * (hd + gy), wd, hd, rows=len(AR))
        for i, k in enumerate(AR):
            v, lo, hi = get(k); y = -i
            ci(ax, lo, hi, y, colof(k), lw=0.8, lim=lim)
            ax.scatter([v], [y], s=12, marker=MK[k], color=colof(k), ec="white", lw=0.4, zorder=3)
            if cx == 0:
                F.label(X0 + 19.0, F.ymm(ax, y), k)
            nums.append(("b", f"{ttl} ({sub}) — {disp(k)}: {v:.2f} ({lo:.2f}–{hi:.2f})"))
        ax.axvline(ref, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
        (xlog if lg else xlin)(ax, tk, lim)
        F.coltitle(ax, ttl, sub)
    F.save(key)
    NUMBERS[key] = [("a", "`E5_trials.csv`, trial `IMmotion150: atezolizumab (+/- bev) vs sunitinib` (`ici_*`, `ctrl_*`, `int_HR`, `int_p`).")] + \
        [("b", "`E5_gastric_arms_pooled.csv` (DL), `E5_prjeb_arms.csv` scale log2, `E5_icb_pooled.csv` all 11, `E5_trials.csv`.")] + nums


# ==========================================================================================
# tcga
# ==========================================================================================
def tcga():
    key = "tcga"
    per = pd.read_csv(R / "E11_per_type.csv"); per = per[per.endpoint == "OS"]
    pool = pd.read_csv(R / "E11_pooled.csv"); pool = pool[pool.endpoint == "OS"]
    dis = pd.read_csv(R / "E11_dissociation_scores.csv").set_index("score")
    G1 = ["SKCM", "KIRC", "LUAD", "LUSC", "BLCA"]; G2 = ["ESCA_AC", "ESCA_SCC", "STAD"]
    G3 = sorted(["ACC", "BRCA", "CESC", "COAD", "GBM", "HNSC", "KIRP", "LGG", "LIHC", "MESO", "OV", "PAAD", "PRAD",
                 "READ", "SARC", "THCA", "UCEC"])
    LAB = {"SKCM": "SKCM*", "ESCA_AC": "ESCA adeno.", "ESCA_SCC": "ESCA squam."}
    F = Fig(h=168)
    nums = []
    YA = 1.7; Y0 = YA + 13.0
    F.heading(1.0, 5.0, YA, "a", "Untreated TCGA, 24 cancer types",
              "overall-survival HR per SD, adjusted for age and stage; pooled, wider of Hartung–Knapp and DerSimonian–Laird CIs")
    LA = (0.42, 2.6)

    def forest(x0, groups, unit_h, labx):
        rows = []
        for head, units, gk, gl in groups:
            rows.append(("head", head)); rows += [("unit", u) for u in units]; rows.append(("pool", (gk, gl))); rows.append(("gap", None))
        rows = rows[:-1]
        HR_ = {"head": 3.6, "unit": unit_h, "pool": 6.4, "gap": 1.6}
        ys, y = [], Y0
        for kind, v in rows:
            ys.append(y + HR_[kind] / 2); y += HR_[kind]
        Y1 = y
        ax = F.ax(x0 + 21.5, Y0, 38.0, Y1 - Y0); ax.set_ylim(Y1, Y0); ax.set_yticks([]); ax.spines["left"].set_visible(False)
        for (kind, v), yy in zip(rows, ys):
            if kind == "head":
                F.txt(x0 + 1.0, yy, v, size=5.4, fontweight="bold")
            elif kind == "unit":
                F.txt(x0 + 20.0, yy, LAB.get(v, v), ha="right")
                for k, dy in [("CAF_ECM", -0.9), ("Immune", 0.9)]:
                    r = per[(per.unit == v) & (per.score == k) & (per.model == "adjusted")].iloc[0]
                    ci(ax, r.lo, r.hi, yy + dy, colof(k), lw=0.7, lim=LA)
                    ax.scatter([r.HR], [yy + dy], s=8, marker=mrk(k), color=colof(k), ec="white", lw=0.3, zorder=3)
            elif kind == "pool":
                gk, gl = v
                ax.add_patch(Rectangle((LA[0], yy - 3.0), LA[1] - LA[0], 6.0, color=PAPER, lw=0, zorder=0))
                F.txt(x0 + 20.0, yy, gl, ha="right", fontweight="bold")
                for k, dy in [("CAF_ECM", -1.3), ("Immune", 1.3)]:
                    for model, filled, ddy in [("adjusted", True, -0.45), ("unadjusted", False, 0.45)]:
                        q = pool[(pool.group == gk) & (pool.score == k) & (pool.model == model)].iloc[0]
                        diamond(ax, q.HR, q.lo, q.hi, yy + dy + ddy, colof(k), h=0.5, filled=filled, lim=LA)
                    qa = pool[(pool.group == gk) & (pool.score == k) & (pool.model == "adjusted")].iloc[0]
                    qu = pool[(pool.group == gk) & (pool.score == k) & (pool.model == "unadjusted")].iloc[0]
                    F.chip(x0 + 61.5, yy + dy, colof(k), s=1.3)
                    F.txt(x0 + 63.8, yy + dy, f"{fmtp(qu.p)} / {fmtp(qa.p)}", size=5.0)
                    nums.append(("a", f"{gl} {k}: adjusted HR {qa.HR:.2f} ({qa.lo:.2f}–{qa.hi:.2f}), unadjusted {qu.HR:.2f}; P unadj/adj {fmtp(qu.p)}/{fmtp(qa.p)}"))
        ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
        xlog(ax, [0.5, 0.7, 1, 1.4, 2], LA)
        F.txt(x0 + 61.5, Y0 - 2.0, "pooled P, unadj. / adj.", size=5.0, color=GREY)
        F.txt(x0 + 40.5, Y1 + 7.0, "OS HR per SD (adjusted)", size=6.0, ha="center")
        return Y1

    YL1 = forest(0.0, [("Checkpoint cancer types (i)", G1, "(i) SKCM KIRC LUAD LUSC BLCA", "Pooled (i)"),
                       ("Esophagogastric (ii)", G2, "(ii) ESCA-AC ESCA-SCC STAD", "Pooled (ii)")], 5.2, 20.0)
    YR1 = forest(91.0, [("Other cancer types", G3, "(iii) all 24", "Pooled, all 24 (iii)")], 6.4, 111.0)
    kx, ky = 93.0, YR1 + 13.0
    for j, (lab, k) in enumerate([("CAF_ECM", "CAF_ECM"), ("Immune axis", "Immune")]):
        F.cv.scatter([kx + j * 18], [ky], s=9, marker=mrk(k), color=colof(k), ec="white", lw=0.3); F.txt(kx + 1.6 + j * 18, ky, lab)
    F.cv.add_patch(Polygon([[kx + 38, ky], [kx + 39.2, ky + 0.8], [kx + 40.4, ky], [kx + 39.2, ky - 0.8]], fc=INK, lw=0))
    F.txt(kx + 41.6, ky, "pooled, adjusted")
    F.cv.add_patch(Polygon([[kx + 62, ky], [kx + 63.2, ky + 0.8], [kx + 64.4, ky], [kx + 63.2, ky - 0.8]], fc="white", ec=INK, lw=0.6))
    F.txt(kx + 65.6, ky, "unadjusted")
    # ---- b: group (i) prognosis vs first-11 response --------------------------------------------
    YB = YL1 + 15.0
    F.heading(1.0, 5.0, YB, "b", "Same cancer types, two outcomes",
              "TCGA group (i) prognosis against response, first 11 cohorts")
    axb = F.ax(17.0, YB + 11.0, 68.0, F.H - (YB + 11.0) - 12.0)
    ks = list(dis.index)
    for k in ks:
        r = dis.loc[k]; col = colof(k); big = k in ("CAF_ECM", "Immune")
        if big:
            axb.plot([r.prog_lo, r.prog_hi], [r.resp_logOR] * 2, color=col, lw=0.7, alpha=0.6)
            axb.plot([r.prog_logHR] * 2, [r.resp_lo, r.resp_hi], color=col, lw=0.7, alpha=0.6)
        axb.scatter([r.prog_logHR], [r.resp_logOR], s=22 if big else 11, marker=mrk(k), color=col, ec="white", lw=0.4,
                    zorder=3 if big else 4)                     # small markers above the large Immune square
    rho = stats.spearmanr(dis.prog_logHR, dis.resp_logOR)[0]
    POS = {"CAF_ECM": (0.13, -0.12), "Barrier": (0.13, 0.06), "TAM_SUP": (0.10, 0.20), "SPP1_arm": (0.10, 0.13),
           "C1Q_arm": (0.10, 0.27), "CXCL9_SPP1": (-0.42, 0.165), "NetImmune": (-0.42, 0.255), "T_CELL": (-0.04, 0.31),
           "APM": (-0.04, 0.42), "IFN": (-0.06, 0.55), "TLS_B": (-0.33, 0.345), "Immune": (-0.43, 0.48)}
    for k in ks:
        r = dis.loc[k]; lx, ly = POS[k]
        axb.annotate(disp(k), (r.prog_logHR, r.resp_logOR), (lx, ly), fontsize=5.0, va="center",
                     fontweight="bold" if k in ("CAF_ECM", "Immune") else "normal",
                     arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4, shrinkA=0.5, shrinkB=2.0))
    axb.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0); axb.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
    xlin(axb, [-0.4, -0.2, 0, 0.2], (-0.45, 0.25), "{:.1f}"); axb.set_ylim(-0.25, 0.72)
    axb.set_yticks([-0.2, 0, 0.2, 0.4, 0.6]); axb.set_yticklabels([mn(f"{t:.1f}") for t in [-0.2, 0, 0.2, 0.4, 0.6]])
    axb.text(0.02, 0.97, f"Spearman ρ = {mn(f'{rho:.2f}')} (12 scores)", transform=axb.transAxes, ha="left", va="top",
             fontsize=5.0, color=SLATE)
    axb.set_xlabel("TCGA group (i), OS log HR per SD (adjusted)"); axb.set_ylabel("Response, log-odds per SD")
    nums.append(("b", f"Spearman ρ over 12 scores = {rho:.2f} (`E11_dissociation_scores.csv`; as in the earlier figure)"))
    F.save(key)
    NUMBERS[key] = [("a", "`E11_per_type.csv` (OS, adjusted), `E11_pooled.csv` (OS; `HR`, `lo`, `hi`, `p` = conservative).")] + \
        [("b", "`E11_dissociation_scores.csv` (prog = group (i) adjusted log HR; resp = first-11 pooled log-OR, DL).")] + nums


# ==========================================================================================
# perfect_metabric
# ==========================================================================================
def perfect_metabric():
    key = "perfect_metabric"
    P = pd.read_csv(R / "E16_perfect_per_score.csv")
    Q = pd.read_csv(R / "E16_perfect_pooled.csv")
    C = pd.read_csv(R / "E16_metabric_cox.csv")
    V = pd.read_csv(R / "E16_metabric_vs_tcga_brca.csv")
    SC = ["CAF_ECM", "Barrier", "TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1", "IFN", "APM", "TLS_B", "T_CELL",
          "Immune", "NetImmune"]
    F = Fig(h=190)
    nums = []
    YA = 1.7; YR = YA + 12.0; HR_ = 64.0

    def rowlabels(ax, x):
        for i, k in enumerate(SC):
            F.label(x, F.ymm(ax, -i), k)

    # a PERFECT per score
    F.heading(1.0, 5.0, YA, "a", "PERFECT, esophageal adenocarcinoma",
              "atezolizumab + chemoradiotherapy; 32 biopsies, 12 responders")
    axa = F.ax(24.0, YR, 58.0, HR_, rows=len(SC))
    LA = (-1.6, 1.8)
    for i, k in enumerate(SC):
        for sc_, dy, filled in [("primary_log2_vst", 0.17, True), ("sens_log2_normcnt", -0.17, False)]:
            r = P[(P.scale == sc_) & (P.score == k)].iloc[0]
            ci(axa, r.lo, r.hi, -i + dy, colof(k), lw=0.6, lim=LA)
            axa.scatter([r.logOR], [-i + dy], s=10, marker=mrk(k), color=colof(k) if filled else "white", ec=colof(k), lw=0.6, zorder=3)
            if filled:
                nums.append(("a", f"{k}: {sg(r.logOR)} ({mn(f'{r.lo:.2f}')} to {mn(f'{r.hi:.2f}')}), P = {fmtp(r.p)}"))
    rows_axis(axa, len(SC)); S.rule(axa, 0); xlin(axa, [-1.5, -1, -0.5, 0, 0.5, 1, 1.5], LA)
    axa.set_xticklabels([mn(t) for t in ["-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5"]])
    for yy in (-1.5, -5.5):
        axa.axhline(yy, color=LIGHT, lw=0.4)
    rowlabels(axa, 22.5)
    F.txt(24.0 + 29.0, YR + HR_ + 7.0, "log-odds of response per SD", size=6.0, ha="center")
    ky = YR + HR_ + 11.5
    F.cv.scatter([25.0], [ky], s=10, color=INK, lw=0); F.txt(26.4, ky, "variance-stabilised log₂ (primary)")
    F.cv.scatter([62.0], [ky], s=10, color="white", ec=INK, lw=0.6); F.txt(63.4, ky, "log₂(count + 1)")
    # b pooled with/without PERFECT
    XB = 104.0
    F.heading(XB - 4.0, XB, YA, "b", "All checkpoint cohorts with PERFECT",
              "including chemotherapy–checkpoint cohorts; wider of HK and DL CIs")
    axb = F.ax(XB + 20.0, YR, 27.0, HR_, rows=len(SC))
    LB = (-0.3, 0.7)
    for i, k in enumerate(SC):
        for st, dy, filled in [("all + chemo-ICB (E12)", 0.17, False), ("all + chemo-ICB + PERFECT", -0.17, True)]:
            q = Q[(Q.set == st) & (Q.score == k)].iloc[0]
            ci(axb, q.lo, q.hi, -i + dy, colof(k) if filled else GREY, lw=0.6, lim=LB)
            axb.scatter([q.logOR], [-i + dy], s=10, marker="D", color=colof(k) if filled else "white",
                        ec=colof(k) if filled else GREY, lw=0.6, zorder=3)
        q = Q[(Q.set == "all + chemo-ICB + PERFECT") & (Q.score == k)].iloc[0]
        q0 = Q[(Q.set == "all + chemo-ICB (E12)") & (Q.score == k)].iloc[0]
        yy = F.ymm(axb, -i)
        F.txt(159.0, yy, f"{int(q0.k)} / {int(q.k)}", ha="right", color=SLATE)
        F.txt(167.5, yy, sg(q.logOR), ha="right")
        F.txt(179.5, yy, fmtp(q.p), ha="right", fontweight="bold" if q.p < 0.05 else "normal")
        nums.append(("b", f"{k}: cohorts {int(q0.k)} before / {int(q.k)} with PERFECT; before {sg(q0.logOR)} (P = {fmtp(q0.p)}), + PERFECT {sg(q.logOR)} (P = {fmtp(q.p)})"))
    rows_axis(axb, len(SC)); S.rule(axb, 0); xlin(axb, [-0.2, 0, 0.2, 0.4, 0.6], LB, "{:.1f}")
    axb.set_xticklabels([mn(t) for t in ["-0.2", "0", "0.2", "0.4", "0.6"]])
    for yy in (-1.5, -5.5):
        axb.axhline(yy, color=LIGHT, lw=0.4)
    rowlabels(axb, XB + 18.5)
    F.txt(159.0, YR - 4.6, "cohorts,", color=GREY, ha="right"); F.txt(159.0, YR - 2.0, "before / with", color=GREY, ha="right")
    F.txt(167.5, YR - 2.0, "log-OR", color=GREY, ha="right"); F.txt(179.5, YR - 2.0, "P", color=GREY, ha="right")
    F.txt(XB + 33.5, YR + HR_ + 7.0, "pooled log-odds per SD", size=6.0, ha="center")
    F.cv.scatter([XB + 10.0], [ky], s=10, marker="D", color="white", ec=GREY, lw=0.6); F.txt(XB + 11.4, ky, "before PERFECT")
    F.cv.scatter([XB + 36.0], [ky], s=10, marker="D", color=INK, lw=0); F.txt(XB + 37.4, ky, "+ PERFECT")
    # c METABRIC
    YC = YA + 94.0; YRC = YC + 12.0; HC = 64.0
    r0 = C[(C.stratum == "all") & (C.endpoint == "DSS") & (C.model == "full")].iloc[0]
    F.heading(1.0, 5.0, YC, "c", "METABRIC breast cancer",
              f"{int(r0.n):,} patients, {int(r0.events)} breast-cancer deaths; disease-specific survival, Cox")
    axc = F.ax(24.0, YRC, 58.0, HC, rows=len(SC))
    LC = (0.6, 1.4)
    for i, k in enumerate(SC):
        for st, mdl, dy, m_ in [("all", "full", 0.24, "o"), ("ERpos", "age_npi", 0.0, "s"), ("ERneg", "age_npi", -0.24, "^")]:
            r = C[(C.stratum == st) & (C.score == k) & (C.endpoint == "DSS") & (C.model == mdl)].iloc[0]
            ci(axc, r.lo, r.hi, -i + dy, colof(k), lw=0.55, lim=LC)
            axc.scatter([r.HR], [-i + dy], s=8, marker=m_, color=colof(k) if st == "all" else "white", ec=colof(k), lw=0.55, zorder=3)
            if st == "all":
                nums.append(("c", f"{k} all (age, NPI, ER): DSS HR {r.HR:.2f} ({r.lo:.2f}–{r.hi:.2f}), P = {fmtp(r.p)}; n {int(r.n)}, events {int(r.events)}"))
    rows_axis(axc, len(SC)); axc.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
    xlog(axc, [0.6, 0.8, 1, 1.25], LC)
    for yy in (-1.5, -5.5):
        axc.axhline(yy, color=LIGHT, lw=0.4)
    rowlabels(axc, 22.5)
    F.txt(24.0 + 29.0, YRC + HC + 7.0, "HR per SD, disease-specific survival", size=6.0, ha="center")
    ky2 = YRC + HC + 11.5
    for j, (m_, lab, filled) in enumerate([("o", "all (+ age, NPI, ER)", True), ("s", "ER+ (+ age, NPI)", False),
                                           ("^", "ER− (+ age, NPI)", False)]):
        F.cv.scatter([25.0 + j * 22.0], [ky2], s=9, marker=m_, color=INK if filled else "white", ec=INK, lw=0.55)
        F.txt(26.4 + j * 22.0, ky2, lab)
    # d METABRIC vs TCGA-BRCA
    F.heading(XB - 4.0, XB, YC, "d", "METABRIC against TCGA-BRCA", "adjusted log HR per SD, 12 scores, two endpoints")
    axd = F.ax(XB + 10.0, YRC, 62.0, HC)
    same = 0
    for _, r in V.iterrows():
        x, y = np.log(r.HR_TCGA_BRCA), np.log(r.HR_METABRIC)
        same += int(np.sign(x) == np.sign(y))
        axd.scatter([x], [y], s=12, marker="o" if r.endpoint == "OS" else "^", color=colof(r.score), ec="white", lw=0.4, zorder=3)
    axd.plot([-0.4, 0.35], [-0.4, 0.35], color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
    axd.axvline(0, color=LIGHT, lw=0.5, zorder=0); axd.axhline(0, color=LIGHT, lw=0.5, zorder=0)
    xlin(axd, [-0.4, -0.2, 0, 0.2], (-0.4, 0.35), "{:.1f}"); axd.set_ylim(-0.4, 0.35)
    axd.set_yticks([-0.4, -0.2, 0, 0.2]); axd.set_yticklabels([mn(f"{t:.1f}") for t in [-0.4, -0.2, 0, 0.2]])
    axd.set_xlabel("TCGA-BRCA, adjusted log HR per SD"); axd.set_ylabel("METABRIC, adjusted log HR per SD")
    axd.text(0.03, 0.97, f"same sign {same}/{len(V)}", transform=axd.transAxes, fontsize=5.0, color=SLATE, va="top")
    for k, (lx, ly) in {"CAF_ECM": (0.2, 0.17), "NetImmune": (-0.33, -0.25)}.items():
        r = V[(V.score == k) & (V.endpoint == "DSS")].iloc[0]
        axd.annotate(k, (np.log(r.HR_TCGA_BRCA), np.log(r.HR_METABRIC)), (lx, ly), fontsize=5.0, fontweight="bold",
                     arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4, shrinkB=2))
    F.cv.scatter([XB + 12.0], [ky2], s=9, marker="o", color=GREY, lw=0); F.txt(XB + 13.4, ky2, "overall survival")
    F.cv.scatter([XB + 38.0], [ky2], s=9, marker="^", color=GREY, lw=0); F.txt(XB + 39.4, ky2, "disease-specific")
    nums.append(("d", f"same sign {same}/{len(V)} (`E16_metabric_vs_tcga_brca.csv`)"))
    F.save(key)
    NUMBERS[key] = [("a", "`E16_perfect_per_score.csv` (scales primary_log2_vst, sens_log2_normcnt).")] + \
        [("b", "`E16_perfect_pooled.csv` sets `all + chemo-ICB (E12)` and `all + chemo-ICB + PERFECT` (`lo`, `hi`, `p` conservative).")] + \
        [("c", "`E16_metabric_cox.csv` (DSS; all: model full; ER strata: age_npi).")] + [("d", "`E16_metabric_vs_tcga_brca.csv`.")] + nums


# ==========================================================================================
# trials (shared)
# ==========================================================================================
TRSC = ["CAF_ECM", "Barrier", "TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1", "IFN", "APM", "TLS_B", "T_CELL",
        "Immune", "NetImmune"]


def trial_column(F, x, y, w, h, vals, title, sub, lim, ticks, labels_x=None, grey=False, ring=True, tsize=5.4):
    """vals: dict score -> (y, se, p). Benefit-oriented estimates per SD with 95% CI."""
    ax = F.ax(x, y, w, h, rows=len(TRSC))
    for i, k in enumerate(TRSC):
        if k not in vals:
            continue
        v, se, p = vals[k]
        col = GREY if grey else colof(k)
        ci(ax, v - 1.96 * se, v + 1.96 * se, -i, col, lw=0.7, lim=lim)
        ax.scatter([v], [-i], s=10, color=col, ec="white", lw=0.3, zorder=3)
        if ring and p < 0.05:
            ax.scatter([v], [-i], s=40, facecolor="none", ec=INK, lw=0.5, zorder=4)
    rows_axis(ax, len(TRSC)); S.rule(ax, 0); xlin(ax, ticks, lim)
    for yy in (-1.5, -5.5):
        ax.axhline(yy, color=LIGHT, lw=0.4)
    F.coltitle(ax, title, sub, size=tsize)
    if labels_x is not None:
        for i, k in enumerate(TRSC):
            F.label(labels_x, F.ymm(ax, -i), k)
    return ax


def pooled_column(F, x, y, w, h, P1, P2, lim, ticks, txt_x, kname):
    ax = F.ax(x, y, w, h, rows=len(TRSC))
    for i, k in enumerate(TRSC):
        r1, r2 = P1.loc[k], P2.loc[k]
        col = colof(k)
        ci(ax, r1.lo, r1.hi, -i + 0.17, col, lw=0.7, lim=lim)
        ax.scatter([r1.est], [-i + 0.17], s=12, marker="D", color=col, ec="white", lw=0.3, zorder=3)
        ci(ax, r2.lo, r2.hi, -i - 0.17, col, lw=0.6, lim=lim)
        ax.scatter([r2.est], [-i - 0.17], s=11, marker="D", color="white", ec=col, lw=0.6, zorder=3)
        yy = F.ymm(ax, -i)
        F.txt(txt_x[0], yy, f"{int(r1[kname])}/{int(r1.k)}", ha="right")
        F.txt(txt_x[1], yy, fmtp(r1.p), ha="right", fontweight="bold" if r1.p < 0.05 else "normal")
    rows_axis(ax, len(TRSC)); S.rule(ax, 0); xlin(ax, ticks, lim)
    for yy in (-1.5, -5.5):
        ax.axhline(yy, color=LIGHT, lw=0.4)
    return ax


# ==========================================================================================
# trials_detail (the first three randomized comparisons, E10)
# ==========================================================================================
def trials_detail():
    key = "trials_detail"
    RI = pd.read_csv(R / "E10_ispy2_interactions.csv")
    W = pd.read_csv(R / "E10_ispy2_within_arm.csv"); W = W[W.model == "HRadj"]
    ND = pd.read_csv(R / "E10_ispy2_null_draws.csv"); NS = pd.read_csv(R / "E10_ispy2_null_summary.csv")
    F = Fig(h=100)
    nums = []
    YB = 1.7; YRB = YB + 12.0; HB = 64.0
    for j, (coh, letter, ttl, n1, n0) in enumerate([("ISPY2_pembro", "a", "I-SPY2 pembrolizumab: 69 vs 179", 69, 179),
                                                    ("ISPY2_durva", "b", "I-SPY2 durvalumab–olaparib: 71 vs 29", 71, 29)]):
        x0 = 1.0 + j * 62.0
        F.heading(x0, x0 + 4.0, YB, letter, ttl, "within-arm OR for pCR per SD, hormone-receptor adjusted")
        ax = F.ax(x0 + 22.5, YRB, 26.0, HB, rows=len(TRSC))
        LO = (0.15, 12.0)
        for i, k in enumerate(TRSC):
            for arm, dy, filled in [("checkpoint", 0.17, True), ("control", -0.17, False)]:
                r = W[(W.cohort == coh) & (W.score == k) & (W.arm == arm)].iloc[0]
                col = colof(k) if filled else GREY
                ci(ax, r.lo, r.hi, -i + dy, col, lw=0.6, lim=LO)
                ax.scatter([r.OR], [-i + dy], s=9, marker="o" if filled else "s", color=col if filled else "white", ec=col,
                           lw=0.6, zorder=3)
            ri = RI[(RI.cohort == coh) & (RI.score == k)].iloc[0]
            yy = F.ymm(ax, -i)
            if j == 0:
                F.label(x0 + 21.5, yy, k)
            F.txt(x0 + 54.5, yy, f"{ri.HRadj_OR:.2f}", ha="right", fontweight="bold" if ri.HRadj_p < 0.05 else "normal")
            F.txt(x0 + 61.5, yy, fmtp(ri.HRadj_p), ha="right")
            nums.append((letter, f"{k}: within-arm OR checkpoint {W[(W.cohort == coh) & (W.score == k) & (W.arm == 'checkpoint')].OR.iloc[0]:.2f}, "
                                 f"control {W[(W.cohort == coh) & (W.score == k) & (W.arm == 'control')].OR.iloc[0]:.2f}; "
                                 f"interaction OR {ri.HRadj_OR:.2f}, P = {fmtp(ri.HRadj_p)}"))
        rows_axis(ax, len(TRSC)); ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
        xlog(ax, [0.2, 0.5, 1, 2, 5, 10], LO)
        for yy in (-1.5, -5.5):
            ax.axhline(yy, color=LIGHT, lw=0.4)
        F.txt(x0 + 54.5, YRB - 2.0, "int. OR", color=GREY, ha="right"); F.txt(x0 + 61.5, YRB - 2.0, "P", color=GREY, ha="right")
        F.txt(x0 + 35.5, YRB + HB + 7.0, "OR for pCR per SD, within arm", size=6.0, ha="center")
    ky2 = YRB + HB + 12.0
    F.cv.scatter([25.0], [ky2], s=9, color=INK, lw=0); F.txt(26.4, ky2, "checkpoint-containing arm")
    F.cv.scatter([60.0], [ky2], s=9, marker="s", color="white", ec=GREY, lw=0.6); F.txt(61.4, ky2, "paclitaxel control")
    # c structure-matched null, I-SPY2
    XD = 128.0
    F.heading(XD - 4.0, XD, YB, "c", "Structure-matched reference", "2,000 random gene sets, I-SPY2")
    items = [("NetImmune", "ISPY2_pembro", "pembrolizumab"), ("NetImmune", "ISPY2_durva", "durvalumab–olaparib"),
             ("NetImmune", "combined", "both, pooled"), ("CAF_ECM", "ISPY2_pembro", "pembrolizumab"),
             ("CAF_ECM", "ISPY2_durva", "durvalumab–olaparib"), ("CAF_ECM", "combined", "both, pooled")]
    hd = 7.8
    for j, (k, c, lab) in enumerate(items):
        yb = YRB + j * (hd + 2.4) + (2.0 if j >= 3 else 0)
        ax = F.ax(XD, yb, 51.0, hd)
        top, p95 = null_hist(ax, ND[f"{k}_{c}"], 1.6, bins=49)
        r = NS[(NS.score == k) & (NS.contrast == c)].iloc[0]
        ax.set_ylim(0, top * 1.95)
        ax.scatter([abs(r.obs_logOR)], [top * 1.15], s=12, marker="D", color=colof(k), ec="white", lw=0.4, zorder=3)
        ax.text(0.01, top * 1.65, f"{k}, {lab}   {fmtp(r.p_emp, 'P = ')}", fontsize=5.0, ha="left", va="center",
                fontweight="bold" if r.p_emp < 0.05 else "normal")
        xlin(ax, [0, 0.5, 1, 1.5], (0, 1.6))
        if j < 5:
            ax.set_xticklabels([])
        nums.append(("c", f"{k} {lab}: |obs| {abs(r.obs_logOR):.2f}, null 95th {p95:.2f}, P = {fmtp(r.p_emp)}"))
    F.txt(XD + 25.5, YRB + HB + 7.0, "Absolute interaction log-OR", size=6.0, ha="center")
    F.save(key)
    NUMBERS[key] = [("a, b", "`E10_ispy2_within_arm.csv` (model HRadj), `E10_ispy2_interactions.csv` (`HRadj_OR`, `HRadj_p`).")] + \
        [("c", "`E10_ispy2_null_draws.csv`, `E10_ispy2_null_summary.csv`.")] + nums


# ==========================================================================================
# trials_forest (six primary comparisons, E21)
# ==========================================================================================
def trials_forest():
    key = "trials_forest"
    IN = pd.read_csv(R / "E21_meta_inputs.csv"); M = pd.read_csv(R / "E21_meta.csv")
    PRIM, SENS = "primary: 6 comparisons", "+ PRINCE + DUTRENEO (8)"
    F = Fig(h=104)
    nums = []
    YA = 1.7; YR = YA + 21.0; HA = 64.0
    F.heading(1.0, 5.0, YA, "a", "Score × treatment interaction, per comparison",
              "benefit-oriented per SD (−log HR for PFS, log OR for pCR); right of zero favors the checkpoint arm")
    cols = [("IMmotion150", "IMmotion150", "renal; PFS", False), ("CheckMate 025", "CheckMate\n025", "renal; PFS", False),
            ("JAVELIN Renal 101", "JAVELIN\nRenal 101", "renal; PFS", False),
            ("I-SPY2 pembrolizumab", "I-SPY2\npembro.", "breast; pCR", False),
            ("I-SPY2 durvalumab-olaparib", "I-SPY2\ndurva.–olap.", "breast; pCR", False),
            ("NeoTRIP", "NeoTRIP", "breast; pCR", False),
            ("PRINCE", "PRINCE", "pancreas; PFS", True), ("DUTRENEO", "DUTRENEO", "bladder; pCR", True)]
    LIM = (-1.6, 1.6)
    CWID, CGAP = 11.6, 2.6
    for j, (tk, tl, sub, grey) in enumerate(cols):
        d = IN[(IN.set == SENS) & (IN.trial == tk)].set_index("score")
        vals = {k: (d.loc[k, "y"], d.loc[k, "se"], d.loc[k, "p_trial"]) for k in TRSC}
        x = 23.5 + j * (CWID + CGAP) + (2.4 if grey else 0)
        trial_column(F, x, YR, CWID, HA, vals, tl, sub, LIM, [-1, 0, 1], labels_x=22.0 if j == 0 else None,
                     grey=grey, tsize=5.0)
        nums.append(("a", f"{tl.replace(chr(10), ' ')}: NetImmune {sg(vals['NetImmune'][0])} (P = {fmtp(vals['NetImmune'][2])}); "
                          f"CAF_ECM {sg(vals['CAF_ECM'][0])} (P = {fmtp(vals['CAF_ECM'][2])})"))
    xs = 23.5 + 6 * (CWID + CGAP) + 2.4
    F.txt(xs, YA + 10.6, "sensitivity", size=5.0, color=GREY)
    F.cv.plot([xs - 1.6, xs - 1.6], [YA + 9.5, YR + HA], color=LIGHT, lw=0.4)
    F.txt(23.5 + (6 * (CWID + CGAP)) / 2, YR + HA + 7.0, "interaction per SD", size=6.0, ha="center")
    ky = YR + HA + 12.0
    F.cv.scatter([25.0], [ky], s=40, facecolor="none", ec=INK, lw=0.5); F.txt(26.8, ky, "P < 0.05")
    # ---- b: pooled ------------------------------------------------------------------------------
    XB = 145.5
    F.heading(XB - 3.5, XB + 0.5, YR - 4.4 - 0.21 * 5.4 * PT - F._ink("Pooled", 5.4), "b", "", None)   # va=bottom sits on the descent
    P1 = M[M.set == PRIM].set_index("score"); P2 = M[M.set == SENS].set_index("score")
    axp = pooled_column(F, XB + 0.5, YR, 16.0, HA, P1, P2, (-0.8, 0.8), [-0.5, 0, 0.5], (169.5, 179.5), "k_benefit")
    F.coltitle(axp, "Pooled", "random effects")
    F.txt(169.5, YR - 2.0, "fav.", color=GREY, ha="right"); F.txt(179.5, YR - 2.0, "P", color=GREY, ha="right")
    for k in TRSC:
        r = P1.loc[k]; r2 = P2.loc[k]
        nums.append(("b", f"pooled 6 {k}: {sg(r.est)} ({mn(f'{r.lo:.2f}')} to {mn(f'{r.hi:.2f}')}), {int(r.k_benefit)}/6, P = {fmtp(r.p)}; "
                          f"with PRINCE, DUTRENEO {sg(r2.est)}, P = {fmtp(r2.p)}"))
    F.txt(XB + 8.5, YR + HA + 7.0, "per SD", size=6.0, ha="center")
    F.cv.scatter([XB - 6.0], [ky], s=11, marker="D", color=INK, lw=0); F.txt(XB - 4.6, ky, "6 primary")
    F.cv.scatter([XB + 8.5], [ky], s=11, marker="D", color="white", ec=INK, lw=0.6); F.txt(XB + 9.9, ky, "+ PRINCE, DUTRENEO")
    F.save(key)
    NUMBERS[key] = [("a", "`E21_meta_inputs.csv` set `+ PRINCE + DUTRENEO (8)` (`y`, `se`, `p_trial`; the six primary rows equal set `primary: 6 comparisons`).")] + \
        [("b", "`E21_meta.csv` sets `primary: 6 comparisons` and `+ PRINCE + DUTRENEO (8)` (`est`, `lo`, `hi`, `p` = conservative; `k_benefit`).")] + nums


# ==========================================================================================
# proteomics
# ==========================================================================================
def proteomics():
    key = "proteomics"
    RR = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
    MOD = dict(CAF_ECM=["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
               TAM_SUP=["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
               APM=["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"], IFN=["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
               TLS_B=["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"], T_CELL=["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"])
    E = pd.read_csv(f"{RR}/expr.tsv", sep="\t", index_col=0)
    E = E.drop(columns=[c for c in ["Mean", "Median", "StdDev"] if c in E.columns])
    det = {g: int(E[g].notna().sum()) if g in E.columns else 0 for G in MOD.values() for g in G}
    Z = (E - E.mean()) / E.std()
    used = {m: [g for g in G if det[g] >= 30] for m, G in MOD.items()}
    score = pd.DataFrame({m: Z[G].mean(1, skipna=True) for m, G in used.items()})
    coh = pd.read_csv(S.ROOT / "figures/Fig7_module_coherence.csv")      # written by figures/prot_stats.py (seed 2026)
    mm = score.corr(method="spearman")
    N = len(Z)
    F = Fig(h=196)
    nums = [("source", f"`GC_TME_Reprogrammer/expr.tsv` ({N} pairs); coherence from `figures/Fig7_module_coherence.csv` "
                       "(prot_stats.py, 500 random protein sets, seed 2026); module correlation recomputed here as in prot_stats.py.")]
    # a detection
    YA = 1.7
    F.heading(1.0, 5.0, YA, "a", "Detection of the 48 marker proteins",
              f"samples quantified of {N} tumor/adjacent pairs; proteins in at least 30 samples enter module scores")
    axa = F.ax(14.0, YA + 16.0, 150.0, 36.0)
    x, pos, genes, cent = 0.0, [], [], []
    for m, G in MOD.items():
        s0 = x
        for g in G:
            v = det[g]
            axa.bar(x, v, width=0.78, color=colof(m) if v >= 30 else LIGHT, lw=0)
            if v == 0:
                axa.text(x, 3, "nd", ha="center", va="bottom", fontsize=5.0, color=GREY, rotation=90)
            pos.append(x); genes.append(g); x += 1
        cent.append((m, (s0 + x - 1) / 2)); x += 0.9
    axa.set_xticks(pos); axa.set_xticklabels(genes, rotation=90, fontsize=5.0)
    axa.axhline(30, color=GREY, lw=0.5, ls=(0, (2, 2)))
    axa.text(pos[-1] + 0.9, 30, "30", fontsize=5.0, color=GREY, va="center")
    axa.set_xlim(-0.8, pos[-1] + 0.8); axa.set_ylim(0, N * 1.05); axa.set_yticks([0, 30, 65, N])
    axa.set_ylabel("samples quantified")
    for m, c in cent:
        xx = 14.0 + 150.0 * (c + 0.8) / (pos[-1] + 1.6)
        F.label(xx - 4.0, YA + 13.5, m, ha="left", size=5.4)
    nums.append(("a", "proteins used: " + ", ".join(f"{m} {len(used[m])}/{len(MOD[m])}" for m in MOD)))
    # b coherence
    YB = YA + 70.0
    F.heading(1.0, 5.0, YB, "b", "Within-module co-variation",
              "median Spearman ρ among a module's proteins; gray, 500 random sets")
    axb = F.ax(22.0, YB + 12.0, 50.0, 40.0, rows=len(coh))
    for i, (_, r) in enumerate(coh.iterrows()):
        axb.plot([r.null_median, r.null_p95], [-i, -i], color=LIGHT, lw=3.0, solid_capstyle="butt", zorder=1)
        axb.plot([r.null_median] * 2, [-i - 0.2, -i + 0.2], color=GREY, lw=0.7)
        axb.scatter([r.rho_obs], [-i], s=14, color=colof(r.module), ec="white", lw=0.4, zorder=3)
        yy = F.ymm(axb, -i)
        F.label(20.5, yy, r.module)
        F.txt(82.0, yy, f"{int(r.n_used)}/{int(r.n_defined)}", ha="right", color=SLATE)
        nums.append(("b", f"{r.module}: ρ {r.rho_obs:.2f}, null median {r.null_median:.2f}, 95th {r.null_p95:.2f}, P {r.p_emp:.3f}"))
    rows_axis(axb, len(coh)); axb.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0)
    xlin(axb, [-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0], (-0.25, 1.0), "{:.1f}")
    F.txt(82.0, YB + 10.0, "proteins", color=GREY, ha="right")
    F.txt(47.0, YB + 12.0 + 40.0 + 7.0, "median Spearman ρ", size=6.0, ha="center")
    # c module correlation
    XC = 100.0
    F.heading(XC - 4.0, XC, YB, "c", "Correlation between module scores", "Spearman ρ across tumor/adjacent pairs")
    mods = list(mm.columns); n = len(mods); cs = 6.4
    X0, Y0c = XC + 18.0, YB + 12.0
    norm = TwoSlopeNorm(vcenter=0, vmin=-1, vmax=1)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            v = mm.values[i, j]; fc = S.DIV(norm(v))
            F.cv.add_patch(Rectangle((X0 + j * cs, Y0c + i * cs), cs - 0.3, cs - 0.3, fc=fc, lw=0))
            r_, g_, b_ = to_rgb(fc)
            F.txt(X0 + j * cs + cs / 2 - 0.15, Y0c + i * cs + cs / 2 - 0.15, mn(f"{v:.2f}"), ha="center",
                  color="white" if 0.2126 * r_ + 0.7152 * g_ + 0.0722 * b_ < 0.45 else INK)
    for i, m in enumerate(mods):
        F.label(X0 - 1.0, Y0c + i * cs + cs / 2, m)
        F.txt(X0 + i * cs + cs / 2, Y0c + n * cs + 1.0, m, rotation=90, ha="center", va="top")
    for t in range(60):
        v = -1 + 2 * (t + 0.5) / 60
        F.cv.add_patch(Rectangle((X0 + n * cs + 3.0, Y0c + n * cs * (1 - (t + 1) / 60)), 2.0, n * cs / 60 + 0.02,
                                 color=S.DIV(norm(v)), lw=0))
    for v, s_ in [(-1, "−1"), (0, "0"), (1, "1")]:
        F.txt(X0 + n * cs + 5.8, Y0c + n * cs * (1 - (v + 1) / 2), s_)
    nums.append(("c", "CAF_ECM–TAM_SUP ρ " + f"{mm.loc['CAF_ECM', 'TAM_SUP']:.2f}; APM–IFN {mm.loc['APM', 'IFN']:.2f}; "
                 f"CAF_ECM–APM {mm.loc['CAF_ECM', 'APM']:.2f}; TLS_B–T_CELL {mm.loc['TLS_B', 'T_CELL']:.2f}"))
    # d heat map
    YD = YB + 67.0
    F.heading(1.0, 5.0, YD, "d", "Marker proteins across tumors",
              f"protein z-scores; {N} tumor/adjacent pairs ordered by protein Barrier score, high to low; missing, white")
    order = score[["CAF_ECM", "TAM_SUP"]].mean(1).sort_values(ascending=False).index
    G_ = [g for m in used for g in used[m]]
    Mx = Z.loc[order, G_].T.clip(-2.5, 2.5)
    axd = F.ax(26.0, YD + 12.0, 140.0, 38.0)
    axd.imshow(Mx.values, aspect="auto", cmap=S.DIV, vmin=-2.5, vmax=2.5, interpolation="nearest")
    start = 0
    for m in used:
        k_ = len(used[m])
        if start > 0:
            axd.axhline(start - 0.5, color="white", lw=1.0)
        y_top = YD + 12.0 + 38.0 * start / len(G_); y_bot = YD + 12.0 + 38.0 * (start + k_) / len(G_)
        F.cv.add_patch(Rectangle((24.3, y_top + 0.2), 1.0, y_bot - y_top - 0.4, color=colof(m), lw=0))
        F.txt(23.3, (y_top + y_bot) / 2, f"{m} ({k_})", ha="right")
        start += k_
    axd.set_xticks([]); axd.set_yticks([])
    for s_ in axd.spines.values():
        s_.set_visible(False)
    for t in range(60):
        v = -2.5 + 5 * (t + 0.5) / 60
        F.cv.add_patch(Rectangle((168.0, YD + 12.0 + 38.0 * (1 - (t + 1) / 60)), 2.0, 38.0 / 60 + 0.02,
                                 color=S.DIV((v + 2.5) / 5), lw=0))
    for v, s_ in [(-2, "−2"), (0, "0"), (2, "2")]:
        F.txt(171.0, YD + 12.0 + 38.0 * (1 - (v + 2.5) / 5), s_)
    F.txt(169.0, YD + 9.5, "z", ha="center", color=SLATE)
    F.save(key)
    NUMBERS[key] = nums


# ==========================================================================================
# drugs
# ==========================================================================================
def drugs():
    key = "drugs"
    RR = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
    SD = pd.read_csv(f"{RR}/results/ai_recommendations_drug_top15_clean.csv", index_col=0)
    MECH = pd.read_csv(f"{RR}/validation/drug/drug_mechanism_table.csv")
    GD = pd.read_csv(S.ROOT / "figures/Fig8_gdsc_points.csv")
    Q = pd.read_csv("/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results/C1_gdsc_fdr.csv").set_index("Drug")
    short = lambda d: d.split(" (")[0]
    F = Fig(h=150)
    nums = []
    YA = 1.7
    F.heading(1.0, 5.0, YA, "a", "Workflow", "")
    steps = ["10 Visium sections:\nprogram and ring metrics", "Barrier- and NetImmune-\nassociated signatures",
             "Perturbation connectivity\n(pertpy), 21 compounds", "Drug score per section\n(10 sections)"]
    w, gap, y0, h = 36.0, 7.5, YA + 6.0, 9.0
    for i, t in enumerate(steps):
        x = 5.0 + i * (w + gap)
        F.cv.add_patch(FancyBboxPatch((x, y0), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=PAPER, ec=INK, lw=0.5))
        F.txt(x + w / 2, y0 + h / 2, t, ha="center", linespacing=1.25)
        if i < len(steps) - 1:
            F.cv.annotate("", xy=(x + w + gap - 0.6, y0 + h / 2), xytext=(x + w + 0.6, y0 + h / 2),
                          arrowprops=dict(arrowstyle="-|>", lw=0.6, color=INK, mutation_scale=6))
    # b ranking
    YB = YA + 26.0
    F.heading(1.0, 5.0, YB, "b", "Compound ranking", "drug score, mean ± SD over 10 sections")
    m = SD.mean().sort_values(ascending=False).head(12); sd = SD.std()[m.index]
    axb = F.ax(26.0, YB + 11.0, 38.0, 100.0, rows=len(m))
    for i, d in enumerate(m.index):
        axb.plot([m[d] - sd[d], m[d] + sd[d]], [-i, -i], color=INK, lw=0.7)
        axb.scatter([m[d]], [-i], s=12, color=INK, ec="white", lw=0.4, zorder=3)
        F.txt(24.5, F.ymm(axb, -i), short(d), ha="right")
        nums.append(("b", f"{short(d)}: {m[d]:.1f} ± {sd[d]:.1f}"))
    rows_axis(axb, len(m)); xlin(axb, [15, 20, 25, 30, 35, 40], (12, 42))
    F.txt(45.0, YB + 11.0 + 100.0 + 7.0, "drug score", size=6.0, ha="center")
    # c mechanism table
    XC = 74.0
    F.heading(XC - 4.0, XC, YB, "c", "Mechanisms of the compounds in b", "")
    rows = list(m.index)
    MT = MECH.set_index("Drug").reindex(rows)
    MT.loc["Daphnetin", ["Pathway", "Axis_effect"]] = ["Protein kinases (EGFR, PKA, PKC)", "not annotated"]
    RT = 4.2; Y0 = YB + 11.0
    F.txt(XC + 2.5, Y0, "Compound", fontweight="bold"); F.txt(XC + 27.0, Y0, "Pathway", fontweight="bold")
    F.txt(XC + 72.0, Y0, "Proposed action", fontweight="bold")
    for i, d in enumerate(rows):
        yy = Y0 + 4.0 + i * RT
        if i % 2 == 0:
            F.cv.add_patch(Rectangle((XC, yy - RT / 2), 104.0, RT, color=PAPER, lw=0, zorder=0))
        eff = MT.loc[d, "Axis_effect"]
        col = STROMA if "anti-Barrier" in eff else (GREY if eff == "not annotated" else IMMUNE)
        F.chip(XC + 0.4, yy, col, s=1.4)
        F.txt(XC + 2.5, yy, short(d))
        pw = (MT.loc[d, "Pathway"].replace("Immunomodulation / anti-angiogenesis", "IMiD, anti-angiogenic")
              .replace("DNA damage / topoisomerase II", "DNA damage / TOP2").replace("VEGFR/PDGFR angiogenesis", "VEGFR/PDGFR")
              .replace("signaling", "signalling"))
        F.txt(XC + 27.0, yy, pw)
        F.txt(XC + 72.0, yy, eff.replace(" / immunomodulatory", " / immunomod."), color=INK)
    # d, e GDSC
    YD = YB + 70.0
    for j, (drug, letter) in enumerate([("Dasatinib", "d"), ("Sorafenib", "e")]):
        x0 = XC - 4.0 + j * 54.0
        d = GD[GD.drug == drug]
        r, p = stats.spearmanr(d.barrier_z, d.ln_ic50)
        q = Q.loc[drug, "Barrier_q"]
        F.heading(x0, x0 + 4.0, YD, letter, drug, f"Spearman r = {mn(f'{r:+.2f}')}; P = {fmtp(p)}; q = {q:.2f}")
        ax = F.ax(x0 + 12.0, YD + 11.0, 38.0, 30.0)
        ax.patch.set_alpha(1)
        b = np.polyfit(d.barrier_z, d.ln_ic50, 1); xs = np.linspace(d.barrier_z.min(), d.barrier_z.max(), 10)
        ax.plot(xs, np.polyval(b, xs), color=GREY, lw=0.8)
        ax.scatter(d.barrier_z, d.ln_ic50, s=9, color=INK, ec="white", lw=0.3, zorder=3)
        ax.set_xlabel("Barrier score (z)"); ax.set_ylabel("ln(IC₅₀)")
        ax.set_yticklabels([mn(t.get_text()) for t in ax.get_yticklabels()]) if False else None
        nums.append((letter, f"{drug}: n = {len(d)}, Spearman r = {r:+.2f}, P = {p:.3f}, q = {q:.2f}"))
    F.save(key)
    NUMBERS[key] = [("source", "`GC_TME_Reprogrammer/results/ai_recommendations_drug_top15_clean.csv` (b), "
                               "`validation/drug/drug_mechanism_table.csv` (c), `figures/Fig8_gdsc_points.csv` and "
                               "`revision/results/C1_gdsc_fdr.csv` (d, e).")] + nums


FUN = dict(icb21_forest=icb21_forest, icb11=icb11, tcga=tcga, perfect_metabric=perfect_metabric,
           trials_detail=trials_detail, trials_forest=trials_forest, proteomics=proteomics, drugs=drugs)

if __name__ == "__main__":
    keys = sys.argv[1:] or list(FUN)
    for k in keys:
        FUN[k]()
        with open(OUTD / f"FigS_{k}_numbers.md", "w") as fh:
            fh.write(f"# Supplementary figure {k} — numbers and sources\n\nScript: `figures/cns/supp/make_supp_pancancer.py {k}`.\n\n")
            for panel, s in NUMBERS.get(k, []):
                fh.write(f"- **{panel}**: {s}\n")
