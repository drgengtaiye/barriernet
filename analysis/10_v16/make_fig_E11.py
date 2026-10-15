"""Fig E11 - prognosis versus checkpoint response in the same cancer types.
(a) Untreated TCGA, OS HR per SD adjusted for age and stage, for CAF_ECM and the Immune axis in every
    cancer type, with pooled estimates for the pre-specified groups (Hartung-Knapp or DerSimonian-Laird,
    whichever interval is wider). (b) Every score: pooled prognostic log HR in the checkpoint cancer types
    (group i) against the pooled log-odds of response in the 11 checkpoint cohorts (E5). (c) Random-gene
    nulls for the group (i) pooled effects. (d) The same comparison within each tumour type."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import bnstyle as bs
from bnstyle import IMMUNE, STROMA, MACRO, GREY, INK, NULL, fmt_p

A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
per = pd.read_csv(f"{A}/E11_per_type.csv")
per = per[(per.endpoint == "OS") & (per.model == "adjusted")].set_index(["unit", "score"])
pool = pd.read_csv(f"{A}/E11_pooled.csv")
poolu = pool[(pool.endpoint == "OS") & (pool.model == "unadjusted")].set_index(["group", "score"])
pool = pool[(pool.endpoint == "OS") & (pool.model == "adjusted")].set_index(["group", "score"])
dis = pd.read_csv(f"{A}/E11_dissociation_scores.csv").set_index("score")
tum = pd.read_csv(f"{A}/E11_dissociation_by_tumour.csv")
nul = pd.read_csv(f"{A}/E11_null.csv"); nd = pd.read_csv(f"{A}/E11_null_draws.csv")

GI, GII, GIII = "(i) SKCM KIRC LUAD LUSC BLCA", "(ii) ESCA-AC ESCA-SCC STAD", "(iii) all 24"
G1 = ["SKCM", "KIRC", "LUAD", "LUSC", "BLCA"]; G2 = ["ESCA_AC", "ESCA_SCC", "STAD"]
OTHER = ["ACC", "BRCA", "CESC", "COAD", "GBM", "HNSC", "KIRP", "LGG", "LIHC", "MESO", "OV", "PAAD",
         "PRAD", "READ", "SARC", "THCA", "UCEC"]
LAB = {"ESCA_AC": "ESCA adeno.", "ESCA_SCC": "ESCA squam.", "SKCM": "SKCM*"}
SCORES = [("CAF_ECM", STROMA, "o", -0.17), ("Immune", IMMUNE, "o", 0.17)]
MINEV = 10


def forest(ax):
    rows = []                                   # (label, kind, key)
    rows += [("Checkpoint cancer types (i)", "head", None)] + [(LAB.get(u, u), "type", u) for u in G1]
    rows += [("Pooled (i)", "pool", GI), ("", "gap", None)]
    rows += [("Oesophagogastric (ii)", "head", None)] + [(LAB.get(u, u), "type", u) for u in G2]
    rows += [("Pooled (ii)", "pool", GII), ("", "gap", None)]
    rows += [("Other types", "head", None)] + [(u, "type", u) for u in OTHER]
    rows += [("Pooled, all 24 (iii)", "pool", GIII)]
    n = len(rows); lo_x, hi_x = 0.45, 2.6
    for i, (lab, kind, key) in enumerate(rows):
        y = n - 1 - i
        if kind == "head":
            ax.text(lo_x * 1.02, y, lab, fontsize=5.4, fontweight="bold", va="center", ha="left", color=INK)
            continue
        if kind == "gap":
            continue
        for k, c, m, dy in SCORES:
            if kind == "type":
                r = per.loc[(key, k)]
                if r.events < MINEV or not np.isfinite(r.se):
                    continue
                ax.plot([max(r.lo, lo_x), min(r.hi, hi_x)], [y + dy] * 2, color=c, lw=0.6, alpha=0.9)
                ax.plot(r.HR, y + dy, m, color=c, ms=2.4, mec="white", mew=0.3)
            else:
                u = poolu.loc[(key, k)]; r = pool.loc[(key, k)]
                ax.plot([u.lo, u.hi], [y + 1.9 * dy] * 2, color=c, lw=0.6)
                ax.plot(u.HR, y + 1.9 * dy, "D", mfc="white", mec=c, ms=2.8, mew=0.7)
                ax.plot([r.lo, r.hi], [y + 0.6 * dy] * 2, color=c, lw=1.0)
                ax.plot(r.HR, y + 0.6 * dy, "D", color=c, ms=3.2, mec="white", mew=0.4)
        if kind == "pool":
            ax.axhspan(y - 0.5, y + 0.5, color=bs.PAPER, zorder=-1, lw=0)
            for j, (k, c, _, _) in enumerate(SCORES):
                u = poolu.loc[(key, k)]; r = pool.loc[(key, k)]
                ax.text(1.03, y - 0.28 + 0.56 * j, f"{fmt_p(u.p, lead='')} / {fmt_p(r.p, lead='')}",
                        transform=ax.get_yaxis_transform(), clip_on=False, fontsize=4.8, ha="left", va="center", color=c)
    ax.set_yticks([n - 1 - i for i, r in enumerate(rows) if r[1] in ("type", "pool")])
    ax.set_yticklabels([r[0] for r in rows if r[1] in ("type", "pool")], fontsize=5.2)
    for t, r in zip(ax.get_yticklabels(), [r for r in rows if r[1] in ("type", "pool")]):
        if r[1] == "pool":
            t.set_fontweight("bold")
    ax.tick_params(axis="y", length=0)
    ax.axvline(1, color=GREY, lw=0.5, ls=(0, (2, 2)))
    ax.set_xscale("log"); ax.set_xlim(lo_x, hi_x); ax.set_ylim(-0.7, n - 0.3)
    ax.set_xticks([0.5, 0.7, 1, 1.4, 2]); ax.set_xticklabels(["0.5", "0.7", "1", "1.4", "2"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_xlabel("OS HR per SD (age + stage adjusted)")
    tr = ax.get_yaxis_transform()
    ax.text(1.03, n - 0.4, "pooled P\nunadj. / adj.", transform=tr, clip_on=False, fontsize=4.8, ha="left", va="center", color=GREY)
    ax.plot([], [], "o", color=STROMA, ms=3, label="CAF_ECM"); ax.plot([], [], "o", color=IMMUNE, ms=3, label="Immune axis")
    ax.plot([], [], "D", color=GREY, ms=3, label="pooled, adjusted")
    ax.plot([], [], "D", mfc="white", mec=GREY, ms=2.8, label="pooled, unadjusted")
    ax.legend(loc="upper center", bbox_to_anchor=(0.45, -0.075), ncol=4, fontsize=5.0, columnspacing=0.8, handletextpad=0.3)
    bs.title(ax, "Untreated TCGA, 24 cancer types")


COL = {"CAF_ECM": STROMA, "Barrier": STROMA, "TAM_SUP": MACRO, "C1Q_arm": MACRO, "SPP1_arm": MACRO,
       "CXCL9_SPP1": INK, "APM": IMMUNE, "IFN": IMMUNE, "TLS_B": IMMUNE, "T_CELL": IMMUNE,
       "Immune": IMMUNE, "NetImmune": IMMUNE}
NAME = {"C1Q_arm": "C1Q arm", "SPP1_arm": "SPP1 arm", "CXCL9_SPP1": "CXCL9:SPP1", "Immune": "Immune",
        "NetImmune": "NetImmune"}
LBL = {  # label position (x, y) in data units; a thin leader joins label and point
    "CAF_ECM": (0.085, -0.075), "Barrier": (0.10, 0.075), "SPP1_arm": (0.10, 0.135), "TAM_SUP": (0.10, 0.195),
    "C1Q_arm": (0.10, 0.255), "Immune": (-0.44, 0.53), "IFN": (-0.08, 0.53), "APM": (-0.02, 0.44), "TLS_B": (-0.23, 0.47),
    "T_CELL": (-0.04, 0.33), "NetImmune": (-0.44, 0.33), "CXCL9_SPP1": (-0.44, 0.17)}
LEFT = {"Immune", "NetImmune", "CXCL9_SPP1"}


def scatter(ax):
    ax.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2))); ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    for k, r in dis.iterrows():
        c = COL[k]; big = k in ("CAF_ECM", "Immune")
        if big:
            ax.plot([r.prog_lo, r.prog_hi], [r.resp_logOR] * 2, color=c, lw=0.8)
            ax.plot([r.prog_logHR] * 2, [r.resp_lo, r.resp_hi], color=c, lw=0.8)
        ax.plot(r.prog_logHR, r.resp_logOR, "o", color=c, ms=3.8 if big else 2.8, mec="white", mew=0.4, zorder=3)
        lx, ly = LBL[k]
        ax.annotate(NAME.get(k, k), (r.prog_logHR, r.resp_logOR), xytext=(lx, ly), textcoords="data", fontsize=5.0,
                    color=c, ha="left" if (lx > r.prog_logHR or k in LEFT) else "right", va="center", fontweight="bold" if big else "normal",
                    arrowprops=dict(arrowstyle="-", color=c, lw=0.35, alpha=0.7, shrinkA=0.5, shrinkB=2))
    rho = spearman(dis.prog_logHR, dis.resp_logOR)
    ax.text(0.99, 0.97, f"Spearman ρ = {rho:.2f}\n(12 overlapping scores)", transform=ax.transAxes, ha="right", va="top",
            fontsize=4.9, color=GREY)
    ax.set_xlim(-0.46, 0.27); ax.set_ylim(-0.25, 0.75)
    ax.set_xlabel("Prognosis: OS log HR per SD\n(group i, age + stage adjusted)")
    ax.set_ylabel("Response: log-OR per SD\n(11 checkpoint cohorts)")
    bs.title(ax, "Same cancer types, two outcomes")


def spearman(a, b):
    return np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1]


def nulls(ax):
    items = [("CAF_ECM", "unadjusted", STROMA, "CAF_ECM, unadjusted"), ("CAF_ECM", "adjusted", STROMA, "CAF_ECM, adjusted"),
             ("Immune", "unadjusted", IMMUNE, "Immune, unadjusted"), ("Immune", "adjusted", IMMUNE, "Immune, adjusted")]
    o = nul.set_index(["score", "model"])
    xmax = max(np.abs(nd).quantile(0.995).max(), np.abs(o.b).max()) * 1.12
    bins = np.linspace(0, xmax, 41)
    for i, (k, md, c, lab) in enumerate(items):
        yb = len(items) - 1 - i
        h, e = np.histogram(np.abs(nd[f"{k}_{md}"]), bins=bins); h = h / h.max() * 0.55
        ax.fill_between(e[:-1], yb, yb + h, step="post", color=NULL, lw=0)
        ax.plot([o.loc[(k, md), "null_p95"]] * 2, [yb, yb + 0.6], color=GREY, lw=0.6)
        ax.plot(abs(o.loc[(k, md), "b"]), yb + 0.2, "D", color=c, ms=3.0, mec="white", mew=0.4, zorder=3)
        pv = o.loc[(k, md), "P_null"]
        ax.text(0.995 * xmax, yb + 0.72, f"{lab}   {fmt_p(pv)}", fontsize=5.0, va="center", ha="right", color=c)
    ax.set_ylim(-0.1, len(items) - 0.1); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(0, xmax); ax.set_xlabel("Absolute pooled OS log HR, group (i)")
    bs.title(ax, "Random-gene null, 1,000 sets")


def by_tumour(ax):
    mk = {"melanoma": "o", "kidney": "s", "lung": "^", "bladder": "D"}
    ax.axhline(0, color=GREY, lw=0.5, ls=(0, (2, 2))); ax.axvline(0, color=GREY, lw=0.5, ls=(0, (2, 2)))
    for k, c in [("CAF_ECM", STROMA), ("Immune", IMMUNE)]:
        t = tum[tum.score == k]
        for _, r in t.iterrows():
            ax.plot([r.prog_logHR - 1.96 * r.prog_se, r.prog_logHR + 1.96 * r.prog_se], [r.resp_logOR] * 2, color=c, lw=0.4, alpha=0.6)
            ax.plot([r.prog_logHR] * 2, [r.resp_logOR - 1.96 * r.resp_se, r.resp_logOR + 1.96 * r.resp_se], color=c, lw=0.4, alpha=0.6)
            ax.plot(r.prog_logHR, r.resp_logOR, mk[r.tumour], color=c, ms=3.2, mec="white", mew=0.4, zorder=3)
    for tname, m in mk.items():
        ax.plot([], [], m, color=GREY, ms=3, label=f"{tname} ({tum[tum.tumour == tname].n_icb.iloc[0]})")
    ax.plot([], [], "s", color=STROMA, ms=3.5, label="CAF_ECM"); ax.plot([], [], "s", color=IMMUNE, ms=3.5, label="Immune")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3, fontsize=4.9, columnspacing=0.6, handletextpad=0.2)
    ax.set_ylim(-1.1, 1.3); ax.set_xlim(-0.68, 0.34)
    ax.set_xlabel("TCGA OS log HR per SD (adjusted)"); ax.set_ylabel("Checkpoint log-OR per SD")
    bs.title(ax, "Within tumour type")


fig = plt.figure(figsize=(bs.W2, 160 * bs.MM))
gs = GridSpec(3, 2, figure=fig, width_ratios=[1.0, 1.05], height_ratios=[1.25, 0.6, 1.0], hspace=0.78, wspace=0.6)
ax_a = fig.add_subplot(gs[:, 0]); ax_b = fig.add_subplot(gs[0, 1]); ax_c = fig.add_subplot(gs[1, 1]); ax_d = fig.add_subplot(gs[2, 1])
forest(ax_a); scatter(ax_b); nulls(ax_c); by_tumour(ax_d)
for ax, L, dx, dy in [(ax_a, "a", -0.32, 1.02), (ax_b, "b", -0.2, 1.08), (ax_c, "c", -0.2, 1.16), (ax_d, "d", -0.2, 1.08)]:
    bs.panel(ax, L, dx=dx, dy=dy)
bs.save(fig, "Fig_E11_prognosis_vs_response")
print("saved Fig_E11_prognosis_vs_response")
