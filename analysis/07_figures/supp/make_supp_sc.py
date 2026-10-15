"""Supplementary Figs S1-S4 (single-cell): programme scores by cell type, composition by
clinical group, CAF vs suppressive-TAM fractions across patients, abundance-weighted attribution."""
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import kruskal, spearmanr
import bnstyle as bs, figstyle as fs
from bnstyle import INK, GREY, STROMA, MACRO, title, fmt_p

F = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures"; D = f"{F}/data/sc"; S = f"{F}/supp"
G = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
CT, ORDER = fs.CT, fs.CTORDER
CTN = {"Tcell": "T", "NK": "NK", "Bcell": "B", "Plasma": "Plasma", "Macro": "Macro", "DC": "DC", "Fibro": "Fibro",
       "Endo": "Endo", "Peri": "Peri", "Epi": "Epi"}
REPORT = {}


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(f"{S}/{stem}.{ext}", dpi=400 if ext == "png" else None)
    plt.close(fig)


def letter(ax, L, dx=-0.16, dy=1.12):
    bs.panel(ax, L, dx=dx, dy=dy)


# ---------------- S1 programme scores by cell type ----------------
obs = pd.read_parquet(f"{D}/cells_obs.parquet"); ct = obs.anno_level1.astype(str)
P = [("score_TLS", "TLS"), ("score_CYTO", "Cytotoxic"), ("score_EXHAUST", "Exhaustion"), ("score_IFN", "IFN"),
     ("score_CAF_ECM", "CAF/ECM"), ("score_TGFb", "TGF-β"), ("score_TAM_SUP", "Suppressive TAM")]
fig, axs = plt.subplots(2, 4, figsize=(bs.W2, 88 * bs.MM))
for k, (ax, (col, lab)) in enumerate(zip(axs.flat, P)):
    data = [obs.loc[ct == t, col].values for t in ORDER]
    bp = ax.boxplot(data, showfliers=False, widths=0.62, patch_artist=True, medianprops=dict(color=INK, lw=0.8),
                    whiskerprops=dict(lw=0.5, color=GREY), capprops=dict(lw=0.5, color=GREY), boxprops=dict(lw=0.4))
    for b, t in zip(bp["boxes"], ORDER):
        b.set_facecolor(CT[t]); b.set_edgecolor(CT[t])
    ax.set_xticks(range(1, 11)); ax.set_xticklabels([CTN[t] for t in ORDER], rotation=60, fontsize=5.0)
    ax.axhline(0, color=GREY, lw=0.4, ls=(0, (2, 2)))
    ax.set_ylabel("score" if k % 4 == 0 else ""); title(ax, lab)
    letter(ax, "abcdefg"[k], dx=-0.28)
axs.flat[-1].set_axis_off()
axs.flat[-1].text(0.02, 0.6, "Box, interquartile range;\nline, median; whiskers,\n1.5 × IQR; 126,559 cells.",
                  fontsize=5.4, color=GREY, va="top", transform=axs.flat[-1].transAxes)
fig.subplots_adjust(hspace=0.62, wspace=0.42)
save(fig, "FigS1")

# ---------------- S2 composition by clinical group ----------------
comp = pd.read_csv(f"{G}/out_scanpy_strict_GSE183904/composition_by_sample_anno_level1.csv").set_index("sample_id")
cl = pd.concat([pd.read_csv(f"{D}/sample_clinical_fibro.csv"), pd.read_csv(f"{D}/sample_clinical_macro.csv")]) \
       .drop_duplicates("sample_id").set_index("sample_id")
cl = cl[cl.sample_type == "Primary_Tumor"]
m = comp.join(cl, how="inner")
m["Lauren"] = m.Laurens.replace({"Mixed intestinal and diffuse": "Mixed", "N.A": np.nan})
m["TCGA_"] = m.TCGA.replace({"N.A": np.nan})
m["Stage_"] = m.Stage.replace({"DM": np.nan})
fig, axs = plt.subplots(1, 3, figsize=(bs.W2, 70 * bs.MM), gridspec_kw=dict(width_ratios=[4, 3, 3]))
minp = {}
for ax, (var, lab, order) in zip(axs, [("Stage_", "AJCC stage", ["I", "II", "III", "IV"]), ("Lauren", "Lauren type", ["Intestinal", "Diffuse", "Mixed"]),
                                         ("TCGA_", "TCGA subtype", ["CIN", "GS", "MSI"])]):
    order = [o for o in order if (m[var] == o).sum() > 0]
    mean = m.groupby(var)[ORDER].mean().loc[order]
    bottom = np.zeros(len(order))
    for t in ORDER:
        ax.bar(range(len(order)), mean[t], bottom=bottom, color=CT[t], width=0.72, edgecolor="white", lw=0.4)
        bottom += mean[t].values
    ns = m[var].value_counts()
    ax.set_xticks(range(len(order))); ax.set_xticklabels([f"{o}\n(n={ns.get(o, 0)})" for o in order], fontsize=5.4)
    ps = []
    for t in ORDER:
        groups = [m.loc[m[var] == o, t].values for o in order if (m[var] == o).sum() >= 2]
        if len(groups) >= 2: ps.append(kruskal(*groups).pvalue)
    minp[lab] = min(ps)
    ax.set_ylim(0, 1); ax.set_ylabel("mean fraction of cells" if var == "Stage_" else "")
    title(ax, f"{lab}\nKruskal–Wallis, smallest {fmt_p(min(ps))}")
    letter(ax, "abc"[list(axs).index(ax)], dx=-0.2)
from matplotlib.patches import Patch
axs[-1].legend(handles=[Patch(color=CT[t], label=CTN[t]) for t in ORDER[::-1]], loc="center left",
               bbox_to_anchor=(1.02, 0.5), fontsize=5.2)
fig.subplots_adjust(wspace=0.35)
save(fig, "FigS2"); REPORT["S2_min_KW_P"] = {k: round(v, 3) for k, v in minp.items()}; REPORT["S2_n_primary_samples"] = len(m)

# ---------------- S3 CAF vs suppressive-TAM fractions ----------------
co = pd.read_csv(f"{G}/results/Fig5_FibroMacro/FibroMacro_patient_coenrichment.csv")
ff = pd.read_csv(f"{G}/results/Fig5_Fibro/Fibro_patient_subtype_fractions.csv")
mf = pd.read_csv(f"{G}/results/Fig5_Macro/Macro_patient_state_fractions.csv")
j = ff.merge(mf, on="Patient_ID")
fig, axs = plt.subplots(1, 2, figsize=(bs.W2 * 0.72, 62 * bs.MM))
for ax, (x, y, xl, yl, L) in zip(axs, [(co.BarrierCAF_frac, co.SupTAM_frac, "barrier-type CAF fraction of fibroblasts",
                                        "suppressive TAM fraction of macrophages", "a"),
                                       (j.CAF_ECM, j.SPP1_TAM, "ECM-remodelling CAF fraction of fibroblasts", "SPP1⁺ TAM fraction of macrophages", "b")]):
    r, p = spearmanr(x, y)
    ax.scatter(x, y, s=12, color=INK, edgecolor="white", lw=0.4)
    ax.set_xlabel(xl); ax.set_ylabel(yl)
    title(ax, f"{len(x)} patients; Spearman ρ = {r:.2f}, {fmt_p(p)}".replace("-", "−"))
    letter(ax, L, dx=-0.22)
    REPORT[f"S3_{L}"] = (len(x), round(r, 3), round(p, 3))
fig.subplots_adjust(wspace=0.45)
save(fig, "FigS3")

# ---------------- S4 abundance-weighted attribution ----------------
c = pd.read_csv(f"{G}/validation/attribution/cell_attribution_contribution_pct.csv", index_col=0)
mods = list(c.columns)[::-1]
fig, ax = plt.subplots(figsize=(bs.W2 * 0.7, 58 * bs.MM))
left = np.zeros(len(mods))
for t in ORDER:
    v = c.loc[t, mods].values
    ax.barh(range(len(mods)), v, left=left, color=CT[t], height=0.7, edgecolor="white", lw=0.4, label=CTN[t])
    for i, (vv, ll) in enumerate(zip(v, left)):
        if vv >= 9: ax.text(ll + vv / 2, i, f"{vv:.0f}", ha="center", va="center", fontsize=5.0, color="white")
    left += v
ax.set_yticks(range(len(mods))); ax.set_yticklabels(mods)
for tl, mname in zip(ax.get_yticklabels(), mods):
    tl.set_color({"CAF_ECM": STROMA, "TAM_SUP": MACRO}.get(mname, bs.IMMUNE))
ax.set_xlim(0, 100); ax.set_xlabel("share of module signal (%), weighted by cell abundance")
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=5.2)
title(ax, "Which cells produce each module's bulk signal")
save(fig, "FigS4")
REPORT["S4_CAF_ECM_fibro_pct"] = round(c.loc["Fibro", "CAF_ECM"], 1)
REPORT["S4_T_CELL_lymphoid_pct"] = round(c.loc[["Tcell", "NK", "Bcell", "Plasma"], "T_CELL"].sum(), 1)
print(REPORT)
