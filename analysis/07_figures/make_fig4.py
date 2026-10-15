"""Figure 4 - NetImmune and overall survival: KM curves, adjusted HRs, the structure-matched
random-gene null (same gene sets in all three cohorts) and head-to-head differences in C."""
import numpy as np, pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import bnstyle as bs
import figstyle as _fs
_orig_panel = bs.panel
bs.panel = lambda ax, l: _orig_panel(ax, l, dx=-0.30, dy=1.10)
from bnstyle import IMMUNE, STROMA, GREY, INK, NULL, title, fmt_p

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
D = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/data"
COH = [("TCGA", "TCGA-STAD", 1 / 30.44), ("GSE84437", "GSE84437", 1.0), ("GSE62254", "GSE62254", 1.0)]
LOW = "#6B7784"


def km(tag, name, scale):
    def draw(ax):
        d = pd.read_csv(f"{D}/km_{tag}.csv").dropna()
        d = d[d.time > 0]
        d["m"] = d.time * scale
        hi = d.NetImmune > d.NetImmune.median()
        ticks = [0, 24, 48, 72, 96]
        for grp, col, lab in [(hi, IMMUNE, "NetImmune high"), (~hi, LOW, "NetImmune low")]:
            k = KaplanMeierFitter().fit(d.m[grp], d.event[grp])
            # truncate where fewer than 10 patients remain at risk
            tcut = min(96, np.sort(d.m[grp].values)[::-1][9])
            sf = k.survival_function_; sf = sf[sf.index <= tcut]
            xs = np.r_[sf.index.values, tcut]; ys = np.r_[sf.iloc[:, 0].values, sf.iloc[-1, 0]]
            ax.step(xs, ys, where="post", color=col, lw=1.0)
            cens = d[grp & (d.event == 0) & (d.m <= tcut)].m
            ax.plot(cens, k.survival_function_at_times(cens), "|", color=col, ms=2.2, mew=0.5)
            ax.text(tcut + 1.5, ys[-1], lab.split()[1], color=col, fontsize=5.6, va="center")
        p = logrank_test(d.m[hi], d.m[~hi], d.event[hi], d.event[~hi]).p_value
        ax.text(0.03, 0.05, f"log-rank {fmt_p(p)}\n{len(d)} patients, {int(d.event.sum())} deaths",
                transform=ax.transAxes, fontsize=5.6, color=INK, va="bottom")
        ax.set_xlim(0, 96); ax.set_ylim(0, 1.02); ax.set_xticks(ticks)
        ax.set_xlabel("Months"); ax.set_ylabel("Overall survival")
        title(ax, name)
        # numbers at risk under the axis
        for grp, col, yy in [(hi, IMMUNE, -0.36), (~hi, LOW, -0.46)]:
            for t in ticks:
                ax.text(t, yy, str(int((d.m[grp] >= t).sum())), transform=ax.get_xaxis_transform(),
                        ha="center", fontsize=5.0, color=col)
        ax.text(-0.02, -0.26, "At risk", transform=ax.transAxes, ha="right", fontsize=5.0, color=INK)
    return draw


def forest(ax):
    uni = pd.read_csv(f"{R}/B4_univariate.csv"); mv = pd.read_csv(f"{R}/B4_multivariate.csv")
    meta = pd.read_csv(f"{R}/B4_meta.csv")
    rows = []
    for tag, name, _ in COH:
        for sc, col in [("NetImmune", IMMUNE), ("Barrier", STROMA)]:
            u = uni[(uni.cohort == tag) & (uni.score == sc)].iloc[0]
            m = mv[(mv.cohort == tag) & (mv.score == sc)].iloc[0]
            rows.append((f"{name}  {sc}", col, (u.HR, u.lo, u.hi), (m.HR, m.lo, m.hi)))
    for sc, col in [("NetImmune", IMMUNE), ("Barrier", STROMA)]:
        m = meta[meta.score == sc].iloc[0]
        rows.append((f"Pooled  {sc}", col, (m.HR, m.lo, m.hi), None))
    y = np.arange(len(rows))[::-1].astype(float)
    y[-2:] -= 0.0
    for yi, (lab, col, u, a) in zip(y, rows):
        ax.plot([u[1], u[2]], [yi + 0.14] * 2, color=col, lw=0.8)
        ax.plot(u[0], yi + 0.14, "o", color=col, ms=3.2, mec="white", mew=0.4)
        if a:
            ax.plot([a[1], a[2]], [yi - 0.14] * 2, color=col, lw=0.8, alpha=0.9)
            ax.plot(a[0], yi - 0.14, "s", mfc="white", mec=col, ms=2.8, mew=0.8)
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.axhline(1.5, color=NULL, lw=0.5)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=5.6)
    ax.set_xscale("log"); ax.set_xticks([0.5, 0.7, 1, 1.4]); ax.set_xticklabels(["0.5", "0.7", "1", "1.4"])
    ax.xaxis.set_minor_formatter(__import__("matplotlib").ticker.NullFormatter())
    ax.set_xlim(0.5, 1.6); ax.set_xlabel("HR per SD")
    ax.plot([], [], "o", color=INK, ms=3, label="unadjusted"); ax.plot([], [], "s", mfc="white", mec=INK, ms=2.8, label="adjusted for stage, age")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=5.0, columnspacing=0.8)
    ax.tick_params(axis="y", length=0)
    title(ax, "Hazard ratios")


def null(ax):
    w = pd.read_csv(f"{D}/null_draws.csv"); o = pd.read_csv(f"{D}/null_obs.csv").set_index("cohort").lhr
    wts = None
    per = pd.read_csv(f"{A}/E2_null_per_cohort.csv").set_index("cohort")
    items = [("TCGA-STAD", w.lhr_TCGA, o["TCGA"], per.loc["TCGA"]),
             ("GSE84437", w.lhr_GSE84437, o["GSE84437"], per.loc["GSE84437"]),
             ("GSE62254", w.lhr_GSE62254, o["GSE62254"], per.loc["GSE62254"])]
    joint = pd.read_csv(f"{A}/E2_null_joint.csv").iloc[0]
    items.append(("Pooled", w.lhr_pooled, np.log(joint.obs_pooled_HR), None))
    bins = np.linspace(0, 0.45, 46)
    for i, (lab, r, ob, p) in enumerate(items):
        yb = len(items) - 1 - i
        h, e = np.histogram(np.abs(r), bins=bins)
        h = h / h.max() * 0.52
        ax.fill_between(e[:-1], yb, yb + h, step="post", color=NULL, lw=0)
        q95 = np.quantile(np.abs(r), 0.95)
        ax.plot([q95, q95], [yb, yb + 0.55], color=GREY, lw=0.6)
        ax.plot(abs(ob), yb + 0.18, "D", color=IMMUNE, ms=3.2, mec="white", mew=0.4)
        if p is not None:
            t1, t2 = f"rank {int(p['rank'])} of 1,001", f"{int(round(p.frac_random_logrank_sig * 100))}% of random axes P < 0.05"
        else:
            t1, t2 = "rank 1 of 1,001", f"{joint.random_sig_all3_same_direction * 100:.1f}% significant in all three"
        ax.text(0.452, yb + 0.84, t1, fontsize=5.0, color=INK, ha="right", va="center")
        ax.text(0.452, yb + 0.64, t2, fontsize=5.0, color=GREY, ha="right", va="center")
        ax.text(0.004, yb + 0.84, lab, fontsize=5.6, ha="left", va="center", fontweight="bold")
    ax.set_ylim(-0.1, len(items)); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 0.46); ax.set_xlabel("Absolute log HR per SD")
    ax.plot([], [], "D", color=IMMUNE, ms=3, label="NetImmune")
    ax.plot([], [], "s", color=NULL, ms=4, label="1,000 random axes")
    ax.plot([], [], "|", color=GREY, ms=5, label="95th percentile")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=5.0, columnspacing=0.8)
    title(ax, "Structure-matched random-gene null")


def h2h(ax):
    d = pd.read_csv(f"{A}/E1_paired_dC.csv")
    comps = [("CD8A", "CD8A"), ("Immune", "Immune only"), ("Barrier", "Barrier"),
             ("CAF_ECM", "CAF_ECM"), ("NetImmune_ssGSEA", "ssGSEA"), ("ESTIMATE", "ESTIMATE")]
    d = d[d.a == "NetImmune"]
    for i, (b, lab) in enumerate(comps):
        y = len(comps) - 1 - i
        s = d[d.b == b]
        for j, (c, mk) in enumerate([("TCGA", "o"), ("GSE84437", "^"), ("GSE62254", "v")]):
            r = s[s.cohort == c].iloc[0]
            ax.plot(r.dC, y + 0.14 + 0.13 * j, mk, color=GREY, ms=2.6, mew=0)
        p = s[s.cohort == "pooled"].iloc[0]
        ax.plot([p.lo, p.hi], [y - 0.14] * 2, color=IMMUNE, lw=0.9)
        ax.plot(p.dC, y - 0.14, "o", color=IMMUNE, ms=3.2, mec="white", mew=0.4)
        ax.text(0.198, y, fmt_p(p.p), fontsize=5.0, ha="right", va="center", color=INK)
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_yticks(range(len(comps))[::-1]); ax.set_yticklabels([c[1] for c in comps], fontsize=5.6)
    ax.set_xlim(-0.03, 0.2); ax.set_xticks([0, 0.05, 0.1]); ax.set_xlabel("ΔC, NetImmune minus comparator")
    ax.tick_params(axis="y", length=0)
    ax.plot([], [], "o", color=IMMUNE, ms=3, label="pooled (95% CI)")
    ax.plot([], [], "o", color=GREY, ms=2.6, label="TCGA-STAD"); ax.plot([], [], "^", color=GREY, ms=2.6, label="GSE84437")
    ax.plot([], [], "v", color=GREY, ms=2.6, label="GSE62254")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=2, fontsize=5.0, columnspacing=0.8)
    title(ax, "Head-to-head, Harrell's C")


cells = [dict(pos=(0, 1, 0, 2), letter="a", draw=km(*COH[0])),
         dict(pos=(0, 1, 2, 4), letter="b", draw=km(*COH[1])),
         dict(pos=(0, 1, 4, 6), letter="c", draw=km(*COH[2])),
         dict(pos=(2, 3, 0, 2), letter="d", draw=forest),
         dict(pos=(2, 3, 2, 4), letter="e", draw=null),
         dict(pos=(2, 3, 4, 6), letter="f", draw=h2h)]
print(bs.compose("Fig4", cells, grid=(3, 6), height_ratios=[1, 0.46, 1.25],
                 width=bs.W2, height=132 * bs.MM, hspace=0.25, wspace=1.6))
