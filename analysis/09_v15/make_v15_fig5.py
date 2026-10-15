"""v15 Figure 5 - fitted models against the fixed score (prognosis only, leave-one-cohort-out).
(a) held-out Harrell's C by model class and cohort; (b) learning curves; (c) paired differences
in C versus NetImmune (bootstrap over held-out patients, pooled by inverse variance)."""
import numpy as np, pandas as pd
import bnstyle as bs
from bnstyle import IMMUNE, FITTED, GREY, INK, NULL, title, fmt_p

_orig_panel = bs.panel
bs.panel = lambda ax, l: _orig_panel(ax, l, dx=-0.62, dy=1.07)
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
BEST = "pcacox"
LABEL = {"en": "Elastic net", "ridge": "Ridge", "pcacox": "PC-Cox", "screenridge": "Screened ridge Cox",
         "superpc": "Supervised PC", "rsf": "Random survival forest", "xgb": "Boosted Cox"}
COH = [("TCGA", "o", "TCGA-STAD"), ("GSE84437", "^", "GSE84437"), ("GSE62254", "v", "GSE62254")]
NET = {"TCGA": 0.5631, "GSE84437": 0.5871, "GSE62254": 0.6176}


def classes(ax):
    d = pd.read_csv(f"{A}/E6_model_classes.csv")
    w = d.pivot(index="model", columns="test", values="C")
    w["mean"] = w[[c for c, _, _ in COH]].mean(1)
    w = w.sort_values("mean")
    net = w.loc["NetImmune (fixed)", "mean"]
    ax.axvline(net, color=IMMUNE, lw=0.6, ls=(0, (2, 2)))
    for y, (m, r) in enumerate(w.iterrows()):
        fixed = "(fixed)" in m
        for j, (c, mk, _) in enumerate(COH):
            ax.plot(r[c], y + 0.2 + 0.1 * j, mk, color=GREY, ms=2.4, mew=0)
        col = IMMUNE if m.startswith("NetImmune") else (INK if fixed else FITTED)
        ax.plot(r["mean"], y - 0.12, "o" if fixed else "D", color=col, ms=3.4, mec="white", mew=0.4)
        ax.text(0.699, y, f"{r['mean']:.3f}", fontsize=5.2, va="center", ha="right", color=col)
    ax.set_yticks(range(len(w))); ax.set_yticklabels([i.replace("(alpha 0.5)", "(α = 0.5)").replace(" (xgboost)", "") for i in w.index], fontsize=5.6)
    ax.set_xlim(0.49, 0.70); ax.set_xticks([0.5, 0.55, 0.6, 0.65]); ax.set_xlabel("Harrell's C in the held-out cohort")
    ax.tick_params(axis="y", length=0)
    for c, mk, lab in COH:
        ax.plot([], [], mk, color=GREY, ms=2.4, label=lab)
    ax.plot([], [], "D", color=FITTED, ms=3, label="fitted, mean"); ax.plot([], [], "o", color=INK, ms=3, label="fixed, mean")
    ax.legend(loc="upper center", bbox_to_anchor=(0.42, -0.095), ncol=5, fontsize=5.0, columnspacing=0.7, handletextpad=0.3)
    title(ax, "Model classes, leave-one-cohort-out")


def curves(ax):
    for key, col, mk in [("en", FITTED, "D"), (BEST, INK, "s")]:
        try:
            d = pd.read_csv(f"{A}/E6_learning_{key}.csv")
        except FileNotFoundError:
            continue
        d["size"] = np.where(d.full, -1, d.n)
        # mean over the three held-out cohorts, per repeat; x = mean training size of that step
        rows = []
        for s, g in d.groupby("size"):
            per = g.groupby("rep").C.mean() if s > 0 else g.groupby("rep").C.mean()
            x = g.n.mean() if s < 0 else s
            rows.append((x, per.mean(), per.quantile(0.1), per.quantile(0.9)))
        rows = sorted(rows)
        xs = [r[0] for r in rows]
        ax.fill_between(xs, [r[2] for r in rows], [r[3] for r in rows], color=col, alpha=0.12, lw=0)
        ax.plot(xs, [r[1] for r in rows], "-", color=col, lw=0.9, marker=mk, ms=3, mec="white", mew=0.4, label=LABEL[key])
    ax.axhline(np.mean(list(NET.values())), color=IMMUNE, lw=0.8, ls=(0, (2, 2)))
    ax.text(110, np.mean(list(NET.values())) + 0.002, "NetImmune, fixed 48 genes", color=IMMUNE, fontsize=5.2, va="bottom")
    ax.set_xlabel("Training patients (two cohorts pooled)"); ax.set_ylabel("Mean held-out C")
    ax.set_xlim(80, 780); ax.set_xticks([100, 200, 400, 700])
    ax.legend(loc="lower right", fontsize=5.2)
    title(ax, "Learning curve")



def deltas(ax):
    d = pd.read_csv(f"{A}/E6_review_dC.csv")
    rows = [("Elastic net", "Elastic net"), ("Best alternative: pcacox", "PC-Cox (20 PCs)"),
            ("NetImmune offset + elastic net", "NetImmune + elastic net\n(offset)"),
            ("Elastic net, inverse-normal genes", "Elastic net,\ninverse-normal genes")]
    wc = pd.read_csv(f"{A}/E6_within_cohort_cv.csv")
    n = len(rows) + 1
    for i, (key, lab) in enumerate(rows):
        y = n - 1 - i
        s = d[d.model == key]
        for j, (c, mk, _) in enumerate(COH):
            r = s[s.test == c].iloc[0]
            ax.plot(r.dC, y + 0.16 + 0.12 * j, mk, color=GREY, ms=2.4, mew=0)
        p = s[s.test == "pooled"].iloc[0]
        col = IMMUNE if "offset" in key else FITTED
        ax.plot([p.lo, p.hi], [y - 0.14] * 2, color=col, lw=0.9)
        ax.plot(p.dC, y - 0.14, "D", color=col, ms=3.2, mec="white", mew=0.4)
        ax.text(0.099, y, f"{p.dC:+.3f}".replace("-", "−") + f"\n{fmt_p(p.p)}", fontsize=5.0, ha="right", va="center", color=INK, linespacing=1.2)
    # within-cohort 10-fold CV elastic net (no bootstrap interval): per-cohort points only
    y = 0
    for j, (c, mk, _) in enumerate(COH):
        r = wc[wc.test == c].iloc[0]
        ax.plot(r.C_cv_mean - r.C_net, y + 0.16 + 0.12 * j, mk, mfc="white", mec=GREY, ms=2.6, mew=0.6)
    m = (wc.C_cv_mean - wc.C_net).mean()
    ax.plot(m, y - 0.14, "D", mfc="white", mec=FITTED, ms=3.2, mew=0.8)
    ax.text(0.099, y, f"{m:+.3f}".replace("-", "−") + "\nno interval", fontsize=5.0, ha="right", va="center", color=INK, linespacing=1.2)
    ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_yticks(range(n)[::-1])
    ax.set_yticklabels([r[1] for r in rows] + ["Elastic net, within-\ncohort 10-fold CV"], fontsize=5.4)
    ax.set_xlim(-0.06, 0.1); ax.set_xticks([-0.05, 0, 0.05]); ax.set_xlabel("ΔC versus NetImmune (95% CI)")
    ax.tick_params(axis="y", length=0)
    title(ax, "Gain over the fixed score")


cells = [dict(pos=(0, 2, 0, 1), letter="a", draw=classes), dict(pos=(0, 1, 1, 2), letter="b", draw=curves),
         dict(pos=(1, 2, 1, 2), letter="c", draw=deltas)]
print(bs.compose("v15_Fig5", cells, grid=(2, 2), width_ratios=[1.1, 1], height_ratios=[1, 1.15], width=bs.W2,
                 height=112 * bs.MM, wspace=0.85, hspace=0.55))
