"""v19 Figure 3 - NetImmune is prognostic in ten gastric cohorts.

(a-c) Kaplan-Meier curves, NetImmune above vs below the cohort median, benchmark cohorts.
(d)   HR per SD for NetImmune and CAF_ECM in 3 benchmark, 3 + 4 locked and 1 sensitivity cohort, with pooled
      random-effects estimates (wider of Hartung-Knapp and DerSimonian-Laird intervals).
(e)   Random-gene nulls for the pooled log HR.
(f)   Paired difference in Harrell's C, NetImmune minus comparator.
(g)   Robustness to random gene loss (E19).
(h)   Fitted transcriptome models: learning curves against the fixed score (E6).
Every estimate is read from analysis/results; KM curves and log-rank P are descriptive re-plots of the
benchmark survival data (figures/data/km_*.csv), as in the v18 Fig. 2."""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-com-gastric-2/727be15e-8d67-407e-a66a-a54d26cc1ebb/scratchpad/v33fig/code")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, IMMUNE, STROMA, FITTED, fmt_p

A = S.RES
D = S.ROOT / "figures/data"
CD = S.OUT / "data"
FS = 5.2
LOW = "#7A8794"
NUM = {}           # every number printed, for Fig3_numbers.md


LET = []           # (letter, title artist, x in mm): letters are placed after layout, top level with the title


def place_letters(fig):
    """Letter ascender line level with the cap height of the panel's bold title (glyph metrics, not text boxes)."""
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    fig.canvas.draw()
    r = fig.canvas.get_renderer(); k = fig.dpi / 72.0
    lp = FontProperties(family=plt.rcParams["font.family"], size=8, weight="bold")
    asc = TextPath((0, 0), "l", prop=lp).get_extents().y1
    for s, art, x in LET:
        prop = art.get_fontproperties()
        d = max(r.get_text_width_height_descent("lp", prop, ismath=False)[2],
                r.get_text_width_height_descent(art.get_text(), prop, ismath=False)[2])
        base = art.get_window_extent(r).y0 + d
        cap = TextPath((0, 0), "H", prop=prop).get_extents().y1
        y = (base + (cap - asc) * k) / fig.bbox.height
        fig.text(x / 180.0, y, s, fontsize=8, fontweight="bold", va="baseline", ha="left", color=INK)


def rec(key, val, src):
    NUM[key] = (val, src)


def neg(s):
    return s.replace("-", "−")


FIGH = 228.0
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import Affine2D


def chip_title(fig, ax, text, col, dy_mm=2.4, fs=6.2):
    """Family-colour chip (1.1 x 1.6 mm rounded bar) + bold ink name, left-aligned above the axes."""
    bb = ax.get_position(); x0 = bb.x0 * 180.0; ym = bb.y1 * FIGH + dy_mm
    fig.add_artist(FancyBboxPatch((x0, ym - 0.8), 1.1, 1.6, boxstyle="round,pad=0,rounding_size=0.3",
                                  transform=Affine2D().scale(S.MM) + fig.dpi_scale_trans, fc=col, ec="none"))
    return fig.text((x0 + 2.0) / 180.0, ym / FIGH, text, ha="left", va="center", fontsize=fs, fontweight="bold",
                    color=INK)


fig = plt.figure(figsize=(S.W2, FIGH * S.MM))

# ------------------------------------------------------------------ a-c  Kaplan-Meier
KM = [("TCGA", "TCGA-STAD", 1 / 30.44), ("GSE84437", "GSE84437", 1.0), ("GSE62254", "GSE62254 (ACRG)", 1.0)]
U = pd.read_csv(A / "E8_univariate.csv")
for i, (tag, name, scale) in enumerate(KM):
    ax = S.ax_mm(fig, 13 + i * 60, 7, 44, 26)
    d = pd.read_csv(D / f"km_{tag}.csv").dropna(); d = d[d.time > 0].copy(); d["m"] = d.time * scale
    hi = d.NetImmune > d.NetImmune.median()
    ticks = [0, 24, 48, 72, 96]
    for grp, col, lab in [(hi, IMMUNE, "high"), (~hi, LOW, "low")]:
        k = KaplanMeierFitter().fit(d.m[grp], d.event[grp])
        tcut = min(96, np.sort(d.m[grp].values)[::-1][9])          # stop where < 10 remain at risk
        sf = k.survival_function_; sf = sf[sf.index <= tcut]
        xs = np.r_[sf.index.values, tcut]; ys = np.r_[sf.iloc[:, 0].values, sf.iloc[-1, 0]]
        ax.step(xs, ys, where="post", color=col, lw=0.9)
        cens = d[grp & (d.event == 0) & (d.m <= tcut)].m
        ax.plot(cens, k.survival_function_at_times(cens), "|", color=col, ms=1.8, mew=0.4)
        ax.text(tcut + 1.5, ys[-1], lab, color=col, fontsize=FS, va="center", ha="left")
        for t in ticks:
            yy = -0.40 if lab == "high" else -0.52
            ax.text(t, yy, str(int((d.m[grp] >= t).sum())), transform=ax.get_xaxis_transform(), ha="center",
                    fontsize=5.0, color=col)
    p = logrank_test(d.m[hi], d.m[~hi], d.event[hi], d.event[~hi]).p_value
    u = U[(U.cohort == tag) & (U.score == "NetImmune")].iloc[0]
    hr, lo, hi_ = np.exp(u.beta), np.exp(u.beta - 1.96 * u.se), np.exp(u.beta + 1.96 * u.se)
    ax.text(0.03, 0.05, f"HR per SD {hr:.2f} ({lo:.2f}–{hi_:.2f})\nlog-rank {fmt_p(p)}",
            transform=ax.transAxes, fontsize=FS, va="bottom", color=INK, linespacing=1.25)
    rec(f"{name}: n/events", f"{len(d)}/{int(d.event.sum())}", f"figures/data/km_{tag}.csv")
    rec(f"{name}: HR per SD", f"{hr:.2f} ({lo:.2f}–{hi_:.2f})", "E8_univariate.csv beta, se (NetImmune)")
    rec(f"{name}: log-rank P", fmt_p(p, ""), f"recomputed from figures/data/km_{tag}.csv (median split)")
    ax.set_xlim(0, 96); ax.set_ylim(0, 1.02); ax.set_xticks(ticks); ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.set_xlabel("Months", labelpad=0.5)
    if i == 0:
        ax.set_ylabel("Overall survival")
    ax.text(-0.075, -0.40, "high", transform=ax.transAxes, ha="right", fontsize=5.0, color=IMMUNE)
    ax.text(-0.075, -0.52, "low", transform=ax.transAxes, ha="right", fontsize=5.0, color=LOW)
    ax.text(-0.075, -0.29, "At risk", transform=ax.transAxes, ha="right", fontsize=5.0, color=SLATE)
    ax.set_title(name, loc="left", fontsize=6.2, fontweight="bold", pad=2)
    ax.text(1.0, 1.015, f"{len(d)} patients, {int(d.event.sum())} deaths", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=5.0, color=SLATE)
    LET.append(("abc"[i], ax._left_title, 1.0 if i == 0 else 13 + i * 60 - 7))

# ------------------------------------------------------------------ d  forest over ten cohorts
E15 = pd.read_csv(A / "E15_per_cohort.csv"); E15 = E15[E15.model == "unadjusted"]
PRIM = "pre-registered exclusions (primary)"
E20 = pd.read_csv(A / "E20_per_cohort.csv")
E20 = E20[(E20.model == "unadjusted") & (E20.duplicate_rule == PRIM)]
E20 = E20[~((E20.cohort == "GSE57303") & (E20.endpoint != "OS"))]
PL = pd.read_csv(A / "E20_pooled.csv"); PL = PL[(PL.model == "unadjusted") & (PL.duplicate_rule == PRIM)]
# v33: the first three locked cohorts as first analyzed (with the eight GSE26253 arrays later found to be shared with
# GSE62254) - the value on which the extension was decided - drawn beside the duplicate-free estimate
PL_OFF = pd.read_csv(S.ROOT / "analysis/results/E30_dedup_runs/off/E20_pooled.csv")
PL_OFF = PL_OFF[(PL_OFF.model == "unadjusted") & (PL_OFF.duplicate_rule == PRIM)]
NE = {k: pd.read_csv(D / f"km_{k}.csv").dropna().query("time > 0") for k, _, _ in KM}

rows = []          # (kind, label, sub, {score: (hr, lo, hi, p, se)})


def est(src, k, sc):
    if src == "bench":
        u = U[(U.cohort == k) & (U.score == sc)].iloc[0]
        return np.exp(u.beta), np.exp(u.beta - 1.96 * u.se), np.exp(u.beta + 1.96 * u.se), None, u.se
    T = E15 if src == "e15" else E20
    r = T[(T.cohort == k) & (T.score == sc) & (T.endpoint.isin(["OS", "RFS"]))].iloc[0]
    return r.HR, r.lo, r.hi, r.p, r.se


def ne(src, k):
    if src == "bench":
        d = NE[k]; return f"{len(d)}/{int(d.event.sum())}"
    T = E15 if src == "e15" else E20
    r = T[(T.cohort == k) & (T.score == "NetImmune")].iloc[0]
    return f"{int(r.n)}/{int(r.events)}"


def ci2(v):
    """Two decimals, or three where two would print a limit of 1.00 on the other side of 1."""
    return f"{v:.3f}" if round(v, 2) == 1.0 and v != 1.0 else f"{v:.2f}"


def pooled(setname, mode=None):
    out = {}
    T = PL_OFF if mode == "off" else PL
    for sc in ("NetImmune", "CAF_ECM"):
        r = T[(T.set == setname) & (T.score == sc)].iloc[0]
        out[sc] = (r.HR, r.lo, r.hi, r.p, None)
    return out


GROUPS = [("Benchmark", "bench", [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254")],
           ("benchmark 3 (reference)", "Pooled, benchmark")),
          ("Locked, first three", "e15", [("GSE15459", "GSE15459"), ("GSE34942", "GSE34942"),
                                          ("GSE26253", "GSE26253 (RFS)")],
           [("locked (first 3), E15 reference", "Pooled, as first analyzed", "off"),
            ("locked (first 3), E15 reference", "Pooled, without shared patients", None)]),
          ("Locked, pre-registered extension", "e20", [("GSE13861", "GSE13861"), ("GSE26899", "GSE26899"),
                                                       ("GSE26901", "GSE26901"), ("GSE29272", "GSE29272")],
           ("locked 7 (primary)", "Pooled, locked 7")),
          ("Sensitivity (third-party survival)", "e20s", [("GSE57303", "GSE57303")], None),
          ("", None, [], ("all 10 gastric", "Pooled, all 10 gastric"))]
y = 0.0; YS = []
for head, src, coh, pl in GROUPS:
    if head:
        rows.append(("head", head, "", None, y)); y += 1.0
    else:
        rows.append(("gap", "", "", None, y - 0.6))
    for k, lab in coh:
        s = "e20" if src == "e20s" else src
        rows.append(("sens" if src == "e20s" else "coh", lab, ne(s, k),
                     {sc: est(s, k, sc) for sc in ("NetImmune", "CAF_ECM")}, y)); y += 1.0
    if pl:
        pls = pl if isinstance(pl, list) else [(pl[0], pl[1], None)]
        for i, (setname, lab, mode) in enumerate(pls):
            rows.append(("pool", lab, "", pooled(setname, mode), y + 0.15))
            y += 1.75
    else:
        y += 0.75
NROW = y

FX = dict(lab=2, w1=27.5, x1=31, w=27, gap=20, x2=0)
X1 = 30.0; WF = 27.0; TX = 1.2; X2 = X1 + WF + 19.5
TOP, H = 55.0, 66.0
axN = S.ax_mm(fig, X1, TOP, WF, H)
axC = S.ax_mm(fig, X2, TOP, WF, H)
XL = (0.30, 3.2)
wmax = max(1 / r[3][sc][4] ** 2 for r in rows if r[0] in ("coh", "sens") for sc in ("NetImmune", "CAF_ECM"))
for ax, sc, col in [(axN, "NetImmune", IMMUNE), (axC, "CAF_ECM", STROMA)]:
    for kind, lab, sub, e, yy in rows:
        if kind in ("head", "gap"):
            continue
        hr, lo, hi_, p, se = e[sc]
        if kind == "pool":
            primary = "locked 7" in lab
            dh = 0.36 if primary else 0.30
            ax.fill([lo, hr, hi_, hr], [yy, yy - dh, yy, yy + dh], color=col if p < 0.05 else "white",
                    ec=col, lw=0.7, zorder=3)
            tr = ax.get_yaxis_transform()
            ax.text(1.04, yy - 0.05, f"{hr:.2f} ({ci2(lo)}–{ci2(hi_)})", transform=tr, fontsize=FS, va="bottom",
                    ha="left", color=INK, fontweight="bold" if primary else "normal", clip_on=False)
            ax.text(1.04, yy + 0.05, fmt_p(p), transform=tr, fontsize=5.0, va="top", ha="left", color=SLATE,
                    clip_on=False)
            rec(f"d {sc} {lab}", f"{hr:.2f} ({ci2(lo)}–{ci2(hi_)}), {fmt_p(p)}",
                "E30_dedup_runs/off/E20_pooled.csv (as first analyzed)" if "first analyzed" in lab else "E20_pooled.csv (unadjusted, primary)")
        else:
            ms = 1.6 + 2.6 * np.sqrt((1 / se ** 2) / wmax)
            a, b = max(lo, XL[0]), min(hi_, XL[1])
            ax.plot([a, b], [yy, yy], color=col, lw=0.7, solid_capstyle="butt", zorder=2)
            if lo < XL[0]:
                ax.plot(XL[0], yy, "<", ms=2, color=col, mew=0, clip_on=False)
            if kind == "sens":
                ax.plot(hr, yy, "s", ms=ms, mfc="white", mec=col, mew=0.7, zorder=3)
            else:
                ax.plot(hr, yy, "s", ms=ms, color=col, mec="white", mew=0.3, zorder=3)
            rec(f"d {sc} {lab}", f"{hr:.2f} ({lo:.2f}–{hi_:.2f})", "E8_univariate / E15_per_cohort / E20_per_cohort")
    S.rule(ax, 1.0)
    ax.set_xscale("log"); ax.set_xlim(*XL)
    ax.set_xticks([0.5, 1, 2]); ax.set_xticklabels(["0.5", "1", "2"]); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(NROW - 0.4, -0.7); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlabel("HR per SD (95% CI)", labelpad=1)
    t = chip_title(fig, ax, sc, col)
    if sc == "NetImmune":
        LET.append(("d", t, 1.0))
    ax.text(1.04, 1.0, "pooled HR, P", transform=ax.transAxes, fontsize=5.0, color=SLATE, va="bottom", ha="left")
    for kind, lab, sub, e, yy in rows:                        # light bands behind group headers
        if kind == "head":
            ax.axhline(yy + 0.5, color=LIGHT, lw=0.4, zorder=0)
        if kind == "gap":
            ax.axhline(yy, color=LIGHT, lw=0.4, zorder=0)
# labels column, shared
tr = axN.get_yaxis_transform()
for kind, lab, sub, e, yy in rows:
    if kind == "gap":
        continue
    if kind == "head":
        axN.text(-0.995, yy + 0.1, lab, transform=tr, fontsize=FS, fontweight="bold", ha="left", va="center",
                 color=INK, clip_on=False)
    elif kind == "pool":
        axN.text(-0.995, yy, lab, transform=tr, fontsize=FS, ha="left", va="center", color=INK, clip_on=False,
                 fontweight="bold" if "locked 7" in lab else "normal")
    else:
        axN.text(-0.94, yy, lab, transform=tr, fontsize=FS, ha="left", va="center",
                 color=SLATE if kind == "sens" else INK, clip_on=False)
        axN.text(-0.03, yy, sub, transform=tr, fontsize=5.0, ha="right", va="center", color=SLATE, clip_on=False)
axN.text(-0.03, 1.0, "n/events", transform=axN.transAxes, fontsize=5.0, color=SLATE, ha="right", va="bottom")
# the decision to extend followed the first-three result
yp3 = [yy for kind, lab, _, _, yy in rows if lab == "Pooled, without shared patients"][0]
yl7 = [yy for kind, lab, _, _, yy in rows if lab == "Pooled, locked 7"][0]
# v33: no arrow; the row "Pooled, as first analyzed" is the result on which the extension was decided


# ------------------------------------------------------------------ e  random-gene nulls
axE = S.ax_mm(fig, 136, TOP, 40, H)
wa = pd.read_csv(D / "null_draws.csv"); wb = pd.read_csv(D / "bio_null_draws.csv")
ja = pd.read_csv(A / "E2_null_joint.csv").iloc[0]
pb = pd.read_csv(A / "E8_biology_null.csv").set_index("cohort")
n20 = pd.read_csv(A / "E20_null.csv").set_index("level")
dr20 = pd.read_csv(CD / "E20_null_draws.csv")
items = [("Benchmark, pooled", "1,000 axes from all genes", np.abs(wa.lhr_pooled.values), abs(np.log(ja.obs_pooled_HR)),
          int(ja.rank_pooled), ja.p_emp_pooled),
         ("Benchmark, pooled", "1,000 axes from immune and\nstromal genes", np.abs(wb.lhr_pooled.values),
          abs(np.log(ja.obs_pooled_HR)), int(pb.loc["pooled", "rank"]), pb.loc["pooled", "p_emp"]),
         ("Locked 7, pooled", "1,000 structure-matched axes", np.abs(dr20["locked 7, pooled"].values),
          abs(n20.loc["locked 7, pooled", "observed_logHR"]), int(n20.loc["locked 7, pooled", "rank"]),
          n20.loc["locked 7, pooled", "P"])]
bins = np.linspace(0, 0.42, 43)
for i, (t1, t2, draws, ob, rk, p) in enumerate(items):
    yb = (len(items) - 1 - i) * 1.0
    h, e = np.histogram(draws, bins=bins); h = h / h.max() * 0.50
    axE.fill_between(e[:-1], yb, yb + h, step="post", color=LIGHT, lw=0)
    q = np.quantile(draws, 0.95)
    axE.plot([q, q], [yb, yb + 0.55], color=GREY, lw=0.6)
    axE.plot(ob, yb + 0.12, "D", ms=3.4, color=IMMUNE, mec="white", mew=0.4, zorder=4)
    short = {0: "All genes", 1: "Immune and stromal genes", 2: "Locked 7, structure-matched"}[i]
    axE.text(0.0, yb + 0.86, short, fontsize=FS, va="center", ha="left", color=INK)
    ptxt = "P ≤ 0.001" if rk == 1 else fmt_p(p)
    axE.text(0.42, yb + 0.30, f"rank {rk}\n{ptxt}", fontsize=5.0, ha="right", va="bottom", color=INK,
             linespacing=1.15)
    rec(f"e {t1} / {t2.replace(chr(10), ' ')}", f"rank {rk} of 1,001; {ptxt}; observed |log HR| {ob:.3f}",
        "E2_null_joint / E8_biology_null / E20_null")
axE.set_ylim(-0.08, len(items) - 0.02); axE.set_xlim(0, 0.42); axE.set_yticks([])
axE.spines["left"].set_visible(False)
axE.set_xticks([0, 0.1, 0.2, 0.3, 0.4]); axE.set_xticklabels(["0", "0.1", "0.2", "0.3", "0.4"])
axE.set_xlabel("Absolute log HR per SD", labelpad=1)
axE.set_title("Random-gene nulls", loc="left", fontsize=6.2, fontweight="bold", pad=2)
axE.plot(0.70, 1.022, "D", ms=2.8, color=IMMUNE, mec="white", mew=0.3, transform=axE.transAxes, clip_on=False)
axE.text(0.735, 1.0, "NetImmune", transform=axE.transAxes, fontsize=5.0, color=SLATE, ha="left", va="bottom")
LET.append(("e", axE._left_title, 132.0))

# ------------------------------------------------------------------ f  head-to-head ΔC
TOP3, H3 = 134.0, 33.0
axF = S.ax_mm(fig, 25, TOP3, 30, H3)
dc = pd.read_csv(A / "E1_paired_dC.csv"); dc = dc[dc.a == "NetImmune"]
tm = pd.read_csv(A / "E8_h2h_TMEscore.csv"); tm = tm[tm.comparator == "TMEscore_PC"].assign(b="TMEscore")
dc = pd.concat([dc[["cohort", "b", "dC", "lo", "hi", "p"]], tm[["cohort", "b", "dC", "lo", "hi", "p"]]])
d15 = pd.read_csv(A / "E15_deltaC.csv"); d20 = pd.read_csv(A / "E20_deltaC.csv")
FROWS = [("bench", "CD8A", "CD8A"), ("bench", "Immune", "Immune only"), ("bench", "Barrier", "Barrier"),
         ("bench", "CAF_ECM", "CAF_ECM"), ("bench", "NetImmune_ssGSEA", "ssGSEA scoring"),
         ("bench", "ESTIMATE", "ESTIMATE"), ("bench", "TMEscore", "TMEscore"),
         ("lock", "CD8A", "CD8A"), ("lock", "ESTIMATE", "ESTIMATE")]
yy = 0.0; ylab = []
for j, (grp, b, lab) in enumerate(FROWS):
    if j == 7:
        yy += 0.9
    if grp == "bench":
        s = dc[dc.b == b]; per = s[s.cohort != "pooled"].dC.values; pr = s[s.cohort == "pooled"].iloc[0]
    else:
        per = np.r_[d15[(d15.comparator == b) & ~d15.cohort.str.startswith("pooled")].dC.values,
                    d20[(d20.comparator == b) & (d20.role == "locked (added)")].dC.values]
        pr = d20[(d20.comparator == b) & (d20.cohort == "locked 7, pooled (inverse variance)")].iloc[0]
    axF.plot(per, np.full(len(per), yy), "o", ms=1.6, color=GREY, mew=0, alpha=0.9, zorder=2)
    axF.plot([pr.lo, pr.hi], [yy, yy], color=IMMUNE, lw=0.8, zorder=3)
    axF.plot(pr.dC, yy, "o", ms=3.0, color=IMMUNE, mec="white", mew=0.4, zorder=4)
    axF.text(1.03, yy, fmt_p(pr.p, ""), transform=axF.get_yaxis_transform(), fontsize=5.0, va="center",
             ha="left", color=INK, fontweight="bold" if pr.p < 0.05 else "normal", clip_on=False)
    rec(f"f {grp} vs {lab}", f"{pr.dC:+.3f} ({pr.lo:+.3f} to {pr.hi:+.3f}), P {fmt_p(pr.p, '')}",
        "E1_paired_dC / E8_h2h_TMEscore / E20_deltaC")
    ylab.append((yy, lab)); yy += 1.0
axF.set_yticks([v for v, _ in ylab]); axF.set_yticklabels([l for _, l in ylab], fontsize=FS)
axF.tick_params(axis="y", length=0); axF.spines["left"].set_visible(False)
axF.set_ylim(yy - 0.4, -1.1)
S.rule(axF, 0)
axF.set_xlim(-0.06, 0.17); axF.set_xticks([0, 0.05, 0.1, 0.15]); axF.set_xticklabels(["0", "0.05", "0.1", "0.15"])
axF.set_xlabel("ΔC, NetImmune − comparator", labelpad=1)
tr = axF.get_yaxis_transform()
axF.text(-0.03, -0.75, "benchmark 3", transform=tr, fontsize=5.0, color=SLATE, ha="right", va="center",
         fontweight="bold", clip_on=False)
axF.text(-0.03, 7.25, "locked 7", transform=tr, fontsize=5.0, color=SLATE, ha="right", va="center",
         fontweight="bold", clip_on=False)
axF.axhline(6.95, color=LIGHT, lw=0.4)
axF.text(1.03, 1.0, "P", transform=axF.transAxes, fontsize=5.0, color=SLATE, ha="left", va="bottom")
axF.set_title("Head-to-head", loc="left", fontsize=6.2, fontweight="bold", pad=7)
for x0, mk, ms, c, t in [(0.03, "o", 3.0, IMMUNE, "pooled (95% CI)"), (0.52, "o", 1.6, GREY, "cohorts")]:
    axF.plot(x0, 1.03, mk, ms=ms, color=c, mec="white", mew=0.3, transform=axF.transAxes, clip_on=False)
    axF.text(x0 + 0.04, 1.0, t, transform=axF.transAxes, fontsize=5.0, color=SLATE, ha="left", va="bottom")
LET.append(("f", axF._left_title, 1.0))

# ------------------------------------------------------------------ g  gene-loss robustness
sp = pd.read_csv(A / "E19_spearman_summary.csv")
su = pd.read_csv(A / "E19_summary.csv")
axG1 = S.ax_mm(fig, 82, TOP3, 36, 14.0)
axG2 = S.ax_mm(fig, 82, TOP3 + 19.0, 36, 14.0)
X = [0, 10, 20, 30]
coh = sp[(sp.perturbation == "drop10") & ~sp.cohort.str.contains("ALL|median")].cohort.unique()
for c in coh:
    v = [1.0] + [sp[(sp.perturbation == f"drop{k}") & (sp.cohort == c)].median_rho.iloc[0] for k in (10, 20, 30)]
    if np.isnan(v[1]):
        continue
    axG1.plot(X, v, color=LIGHT, lw=0.5, zorder=1)
allr = sp[sp.cohort == "ALL (cohort x draw)"].set_index("perturbation")
med = [1.0] + [allr.loc[f"drop{k}", "median_rho"] for k in (10, 20, 30)]
p5 = [1.0] + [allr.loc[f"drop{k}", "p5_rho"] for k in (10, 20, 30)]
axG1.fill_between(X, p5, med, color=IMMUNE, alpha=0.15, lw=0)
axG1.plot(X, med, "-o", color=IMMUNE, ms=2.4, lw=0.9, mec="white", mew=0.3, zorder=3)
axG1.text(30, 1.0005, f"median ρ {med[-1]:.2f}", fontsize=5.0, ha="right", va="center", color=INK)
rec("g median rho at 30% removed", f"{med[-1]:.3f}", "E19_spearman_summary.csv ALL (cohort x draw), drop30")
rec("g 5th percentile rho at 30% removed", f"{p5[-1]:.3f}", "E19_spearman_summary.csv ALL, drop30")
axG1.set_ylim(0.9, 1.012); axG1.set_yticks([0.9, 0.95, 1]); axG1.set_yticklabels(["0.9", "0.95", "1"])
axG1.set_xlim(-2, 32); axG1.set_xticks(X); axG1.set_xticklabels([])
axG1.set_ylabel("ρ with full", labelpad=1)
axG1.set_title("Random gene loss", loc="left", fontsize=6.2, fontweight="bold", pad=2)

hr = su[su.estimate == "NetImmune HR per SD (benchmark 3, OS)"].set_index("perturbation")
for x, k in zip(X, ["full", "drop10", "drop20", "drop30"]):
    r = hr.loc[k]
    if k != "full":
        axG2.plot([x, x], [r["p2.5"], r["p97.5"]], color=IMMUNE, lw=0.8)
        rec(f"g pooled HR {k}", f"median {r['median']:.3f} ({r['p2.5']:.3f}–{r['p97.5']:.3f}); direction kept "
            f"{r.frac_direction_kept:.0%}", "E19_summary.csv")
    axG2.plot(x, r["median"], "o", ms=2.6, color=IMMUNE, mec="white", mew=0.3)
S.hrule(axG2, 1.0)
axG2.set_ylim(0.62, 1.04); axG2.set_yticks([0.7, 0.85, 1]); axG2.set_yticklabels(["0.7", "0.85", "1"])
axG2.set_xlim(-2, 32); axG2.set_xticks(X); axG2.set_xticklabels(["0", "10", "20", "30"])
axG2.set_xlabel("Genes removed (%)", labelpad=1); axG2.set_ylabel("pooled HR", labelpad=1)

LET.append(("g", axG1._left_title, 72.0))

# ------------------------------------------------------------------ h  learning curves
axH = S.ax_mm(fig, 136, TOP3, 40, H3 - 3)
NET = pd.read_csv(A / "E6_model_classes.csv")
netC = NET[NET.model == "NetImmune (fixed)"].C.mean()
rec("h NetImmune mean held-out C", f"{netC:.3f}", "E6_model_classes.csv NetImmune (fixed), mean of 3 cohorts")
for key, col, mk, lab, (lx, ly) in [("pcacox", FITTED, "o", "Cox on 20 PCs", (560, 0.611)),
                                    ("en", "#8E9AA6", "s", "Elastic net", (560, 0.576))]:
    d = pd.read_csv(A / f"E6_learning_{key}.csv")
    d["size"] = np.where(d.full, -1, d.n)
    pts = []
    for s_, g in d.groupby("size"):
        per = g.groupby("rep").C.mean()
        x = g.n.mean() if s_ < 0 else s_
        pts.append((x, per.mean(), per.quantile(0.1), per.quantile(0.9)))
    pts = sorted(pts); xs = [p[0] for p in pts]
    axH.fill_between(xs, [p[2] for p in pts], [p[3] for p in pts], color=col, alpha=0.13, lw=0)
    axH.plot(xs, [p[1] for p in pts], "-", marker=mk, color=col, lw=0.9, ms=2.4, mec="white", mew=0.3)
    axH.text(lx, ly, lab, fontsize=5.0, color=INK, va="center", ha="center")
    rec(f"h {lab} mean C at full size", f"{pts[-1][1]:.3f} (n = {xs[-1]:.0f})", f"E6_learning_{key}.csv")
axH.axhline(netC, color=IMMUNE, lw=0.8, ls=(0, (2.5, 1.5)))
axH.text(100, netC + 0.003, f"NetImmune, fixed ({netC:.3f})", fontsize=5.0, color=INK, va="bottom")
axH.set_xlim(80, 760); axH.set_xticks([100, 200, 400, 700]); axH.set_ylim(0.52, 0.64)
axH.set_yticks([0.52, 0.56, 0.6, 0.64])
axH.set_xlabel("Training patients", labelpad=1); axH.set_ylabel("Held-out C", labelpad=1)
axH.set_title("Fitted models", loc="left", fontsize=6.2, fontweight="bold", pad=2)

LET.append(("h", axH._left_title, 127.0))

# ------------------------------------------------------------------ i  module-level prognosis, ten cohorts
EP = pd.read_csv(A / "E20_estimates_used_for_pooling.csv")
EP = EP[(EP.model == "unadjusted") & (EP.duplicate_rule == PRIM) &
        EP.group.isin(["benchmark", "locked (first 3)", "locked (added)"])]
P10 = PL[PL.set == "all 10 gastric"].set_index("score")
HCOLS = [("CAF_ECM", STROMA), ("TAM_SUP", S.MACRO), ("APM", IMMUNE), ("IFN", IMMUNE), ("TLS_B", IMMUNE),
         ("T_CELL", IMMUNE), ("Immune", IMMUNE), ("NetImmune", IMMUNE)]
HROWS = [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254"), None,
         ("GSE15459", "GSE15459"), ("GSE34942", "GSE34942"), ("GSE26253", "GSE26253 (RFS)"), None,
         ("GSE13861", "GSE13861"), ("GSE26899", "GSE26899"), ("GSE26901", "GSE26901"), ("GSE29272", "GSE29272"), None,
         ("pooled", "Pooled, 10 cohorts"), ("P", "P")]
TOP4, CW, RH = 193.0, 8.4, 2.25
XH = 89.0
gaps = {2: 0.0}
ys, y = [], 0.0
for r in HROWS:
    if r is None:
        y += 0.45; continue
    ys.append((r, y)); y += 1.0
NHR = y
axI = S.ax_mm(fig, XH, TOP4, CW * len(HCOLS), RH * NHR)
axI.set_xlim(0, len(HCOLS)); axI.set_ylim(NHR, 0); S.clean(axI)
VM = 0.5
for (k, lab), yy in ys:
    for j, (sc, col) in enumerate(HCOLS):
        if k == "P":
            pv = P10.loc[sc, "p"]
            axI.text(j + 0.5, yy + 0.5, fmt_p(pv, ""), ha="center", va="center", fontsize=5.0, color=SLATE)
            continue
        if k == "pooled":
            hr = P10.loc[sc, "HR"]
            rec(f"j pooled {sc}", f"{hr:.2f} ({P10.loc[sc, 'lo']:.2f}–{P10.loc[sc, 'hi']:.2f}), P {fmt_p(P10.loc[sc, 'p'], '')}",
                "E20_pooled.csv all 10 gastric")
        else:
            e = EP[(EP.cohort == k) & (EP.score == sc)].iloc[0]
            hr = np.exp(e.beta)
            rec(f"j {lab} {sc}", f"{hr:.2f}", "E20_estimates_used_for_pooling.csv (unadjusted, primary)")
        v = np.log(hr)
        axI.add_patch(plt.Rectangle((j + 0.04, yy + 0.05), 0.92, 0.90, fc=S.DIV(0.5 + v / (2 * VM)), ec="none"))
        axI.text(j + 0.5, yy + 0.52, f"{hr:.2f}", ha="center", va="center", fontsize=5.0,
                 color="white" if abs(v) > 0.33 else INK, fontweight="bold" if k == "pooled" else "normal")
    fig.text((XH - 1.2) / 180.0, 1 - (TOP4 + (yy + 0.5) * RH) / FIGH, lab, ha="right", va="center",
             fontsize=5.0 if k != "pooled" else 5.2, color=SLATE if k == "P" else INK,
             fontweight="bold" if k == "pooled" else "normal")
for j, (sc, col) in enumerate(HCOLS):                      # chip strip + rotated ink header
    x0 = XH + j * CW
    fig.add_artist(FancyBboxPatch((x0 + CW * 0.08, FIGH - TOP4 + 0.5), CW * 0.84, 1.0,
                                  boxstyle="round,pad=0,rounding_size=0.3",
                                  transform=Affine2D().scale(S.MM) + fig.dpi_scale_trans, fc=col, ec="none"))
    fig.text((x0 + CW * 0.30) / 180.0, (FIGH - TOP4 + 2.1) / FIGH, sc, rotation=35, ha="left", va="bottom",
             fontsize=5.0, color=INK, rotation_mode="anchor")
tI = fig.text(1.0 / 180.0 + 0.02, 1 - (TOP4 - 11.5) / FIGH, "Each module and the Immune and NetImmune axes", fontsize=6.2,
              fontweight="bold", color=INK, ha="left", va="baseline")
tI.set_x(68.0 / 180.0)
# colour key
axK = S.ax_mm(fig, XH + CW * len(HCOLS) + 5.0, TOP4 + 5, 2.0, 16)
grad = np.linspace(-VM, VM, 101).reshape(-1, 1)
axK.imshow(grad, cmap=S.DIV, vmin=-VM, vmax=VM, aspect="auto", origin="lower", extent=(0, 1, -VM, VM))
axK.set_xticks([]); axK.yaxis.tick_right()
axK.set_yticks(np.log([0.67, 1, 1.5])); axK.set_yticklabels(["0.67", "1", "1.5"], fontsize=5.0)
for sp in axK.spines.values():
    sp.set_visible(False)
axK.text(0.0, 1.06, "HR per SD", transform=axK.transAxes, fontsize=5.0, ha="left", va="bottom", color=SLATE)
LET.append(("j", tI, 63.5))

# ------------------------------------------------------------------ j  beyond age and stage
LR = pd.read_csv(A / "E8_lrt.csv")
CI = pd.read_csv(S.ROOT.parent / "revision/results/B4_cindex.csv")
axJ = S.ax_mm(fig, 20, TOP4 - 2.5, 32, 30)
JR = [("TCGA", "TCGA-STAD"), ("GSE84437", "GSE84437"), ("GSE62254", "GSE62254")]
for i, (k, lab) in enumerate(JR):
    c = CI[(CI.cohort == k) & (CI.score == "NetImmune")].iloc[0]
    l = LR[(LR.cohort == k) & (LR.base == "age+stage") & (LR.added == "NetImmune")].iloc[0]
    axJ.plot([c.C_clin, c.C_full], [i, i], color=GREY, lw=0.9, solid_capstyle="butt", zorder=2)
    axJ.plot(c.C_clin, i, "o", ms=3.0, mfc="white", mec=INK, mew=0.7, zorder=3)
    axJ.plot(c.C_full, i, "o", ms=3.2, color=IMMUNE, mec="white", mew=0.4, zorder=3)
    axJ.text(1.03, i, f"χ² {l.chisq:.1f}\n{fmt_p(l.p)}", transform=axJ.get_yaxis_transform(), fontsize=5.0,
             va="center", ha="left", color=INK, linespacing=1.1, clip_on=False)
    rec(f"i {lab}", f"C {c.C_clin:.3f} -> {c.C_full:.3f}; LR chi2 {l.chisq:.1f}, P {fmt_p(l.p, '')} (n = {int(l.n)})",
        "revision/results/B4_cindex.csv; E8_lrt.csv")
axJ.set_yticks(range(len(JR))); axJ.set_yticklabels([l for _, l in JR], fontsize=5.2)
axJ.tick_params(axis="y", length=0); axJ.spines["left"].set_visible(False)
axJ.set_ylim(len(JR) - 0.5, -0.6)
axJ.set_xlim(0.65, 0.755); axJ.set_xticks([0.66, 0.70, 0.74]); axJ.set_xticklabels(["0.66", "0.70", "0.74"])
axJ.set_xlabel("Apparent Harrell's C", labelpad=1)
tJ = fig.text(5.5 / 180.0, 1 - (TOP4 - 11.5) / FIGH, "Beyond age and stage", fontsize=6.2, fontweight="bold",
               color=INK, ha="left", va="baseline")
axJ.plot(0.0, 1.13, "o", ms=2.8, mfc="white", mec=INK, mew=0.7, transform=axJ.transAxes, clip_on=False)
axJ.text(0.04, 1.13, "age + stage", transform=axJ.transAxes, fontsize=5.0, color=SLATE, va="center")
axJ.plot(0.52, 1.13, "o", ms=3.0, color=IMMUNE, mec="white", mew=0.3, transform=axJ.transAxes, clip_on=False)
axJ.text(0.56, 1.13, "+ NetImmune", transform=axJ.transAxes, fontsize=5.0, color=SLATE, va="center")
axJ.text(1.03, 1.13, "LR test", transform=axJ.transAxes, fontsize=5.0, color=SLATE, va="center")
LET.append(("i", tJ, 1.0))

place_letters(fig)
S.save(fig, "Fig3")
with open(S.OUT / "Fig3_numbers.md", "w") as fh:
    fh.write("# Fig. 3 - numbers printed on the figure and their sources\n\n| item | value | source |\n|---|---|---|\n")
    for k, (v, s) in NUM.items():
        fh.write(f"| {k} | {v} | {s} |\n")
print("ok")
