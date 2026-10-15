"""Supplementary Fig. S11 (was main Figure 7): protein-level measurability in PDC000214."""
"""Figure 7 - protein-level measurability (PDC000214, iTRAQ tumour/adjacent log ratios).

Redrawn from the protein matrix with the barriernet 48-gene module definitions. The v13
panels compared Barrier-high with NetImmune-high tumours where both groups were defined
from the same protein scores (Cliff's delta about 1 by construction); those are replaced by
tests that the definitions do not fix: detection of each marker protein, within-module
co-variation against random protein sets, and correlation between different modules.
"""
import sys; sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import bnstyle as bs
from prot_stats import MOD, det, used, Z, score, coh, mm
import figstyle as fs

FAMCOL = {"CAF_ECM": bs.STROMA, "TAM_SUP": bs.MACRO, "APM": bs.IMMUNE, "IFN": bs.IMMUNE,
          "TLS_B": bs.IMMUNE, "T_CELL": bs.IMMUNE}
LAB = {"CAF_ECM": "CAF_ECM", "TAM_SUP": "TAM_SUP", "APM": "APM", "IFN": "IFN", "TLS_B": "TLS_B", "T_CELL": "T_CELL"}
N = len(Z)


def draw_detection(ax):
    x, xt, xl, bounds = 0, [], [], []
    for m, G in MOD.items():
        start = x
        for g in G:
            v = det[g]
            ax.bar(x, v, width=0.8, color=FAMCOL[m] if v >= 30 else bs.LIGHT, lw=0)
            if v == 0:
                ax.text(x, 3, "nd", ha="center", va="bottom", fontsize=5.0, color=bs.GREY, rotation=90)
            x += 1
        xt.append((start + x - 1) / 2); xl.append(LAB[m]); bounds.append(x - 0.5); x += 0.8
    genes = [g for G in MOD.values() for g in G]
    pos, x = [], 0.0
    for m, G in MOD.items():
        pos += [x + i for i in range(len(G))]; x += len(G) + 0.8
    ax.set_xticks(pos); ax.set_xticklabels(genes, rotation=90, fontsize=5.0)
    for p, m in zip(xt, xl):
        ax.text(p, N * 1.13, m, ha="center", va="bottom", fontsize=5.8, color=fs.INK)
    ax.axhline(30, color=bs.GREY, lw=0.5, ls=(0, (2, 2)))
    ax.text(pos[-1] + 0.8, 30, "≥30 samples:\nused in scores", fontsize=5.0, color=bs.GREY, va="center", ha="left")
    ax.set_ylim(0, N * 1.12); ax.set_xlim(-0.8, pos[-1] + 0.8)
    ax.set_ylabel(f"samples quantified (of {N})")
    ax.set_yticks([0, 30, 65, 130])
    bs.title(ax, "Detection of the 48 marker proteins", pad=12)


def draw_coherence(ax):
    y = np.arange(len(coh))[::-1]
    for yi, (_, r) in zip(y, coh.iterrows()):
        ax.plot([r.null_median, r.null_p95], [yi, yi], color=bs.LIGHT, lw=3, solid_capstyle="round")
        ax.plot(r.null_median, yi, "|", color=bs.GREY, ms=5, mew=0.8)
        ax.plot(r.rho_obs, yi, "o", color=FAMCOL[r.module], ms=4.2, mec="white", mew=0.5)
        ax.text(1.0, yi, f"{int(r.n_used)}/{int(r.n_defined)}", fontsize=5, va="center", ha="right", color=bs.GREY)
    ax.text(1.0, len(coh) - 0.35, "proteins\nused", fontsize=5.0, va="bottom", ha="right", color=bs.GREY)
    ax.set_yticks(y); ax.set_yticklabels([LAB[m] for m in coh.module])
    ax.set_xlim(-0.2, 1.02); ax.axvline(0, color=fs.INK, lw=0.4)
    ax.set_xlabel("median Spearman ρ between proteins of a module")
    ax.text(0.08, len(coh) - 0.35, "grey: random protein sets of equal size\n(median, range to 95th percentile)", fontsize=5.0, color=bs.GREY,
            va="bottom", ha="left")
    ax.set_ylim(-0.6, len(coh) + 0.2)
    bs.title(ax, "Within-module co-variation", pad=4)


def draw_modcorr(ax):
    mods = list(mm.columns)
    V = mm.values.copy(); np.fill_diagonal(V, np.nan)
    im = ax.imshow(V, cmap=fs.DIV, vmin=-1, vmax=1)
    for i in range(len(mods)):
        for j in range(len(mods)):
            if i != j:
                v = mm.values[i, j]
                ax.text(j, i, f"{v:.2f}".replace("-", "−"), ha="center", va="center", fontsize=5.0,
                        color="white" if abs(v) > 0.6 else fs.INK)
    ax.set_xticks(range(len(mods))); ax.set_xticklabels([LAB[m] for m in mods], rotation=45, ha="right")
    ax.set_yticks(range(len(mods))); ax.set_yticklabels([LAB[m] for m in mods])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03, ticks=[-1, 0, 1])
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5, length=1.5)
    cb.set_label("Spearman ρ", fontsize=5.4)
    bs.title(ax, "Correlation between module scores", pad=4)


def draw_heatmap(ax):
    genes, rowcol, seps = [], [], []
    for m, G in used.items():
        genes += G; rowcol += [FAMCOL[m]] * len(G); seps.append(len(genes))
    order = score[["CAF_ECM", "TAM_SUP"]].mean(1).sort_values(ascending=False).index
    M = Z.loc[order, genes].T.clip(-2.5, 2.5)
    im = ax.imshow(M.values, aspect="auto", cmap=fs.DIV, vmin=-2.5, vmax=2.5, interpolation="nearest")
    start = 0
    for (m, G), s_ in zip(used.items(), seps):
        if s_ < len(genes):
            ax.axhline(s_ - 0.5, color="white", lw=1.2)
        ax.plot([-3, -3], [start - 0.4, s_ - 0.6], color=FAMCOL[m], lw=2.2, clip_on=False, solid_capstyle="butt")
        ax.text(-5, (start + s_ - 1) / 2, f"{LAB[m]} ({len(G)})", ha="right", va="center", fontsize=5.4)
        start = s_
    ax.set_yticks([]); ax.set_xlim(-0.5, M.shape[1] - 0.5)
    ax.set_xticks([]); ax.set_xlabel(f"{N} tumour/adjacent pairs, ordered by protein Barrier score, high to low")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(axis="y", length=0)
    cb = plt.colorbar(im, ax=ax, fraction=0.015, pad=0.01, ticks=[-2, 0, 2])
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5, length=1.5)
    cb.set_label("protein z-score\n(missing = white)", fontsize=5.2)
    bs.title(ax, "Marker proteins across tumours", pad=4)


fig = plt.figure(figsize=(bs.W2, 185 * bs.MM))
gs = GridSpec(3, 2, figure=fig, height_ratios=[0.95, 1.0, 1.15], width_ratios=[1.15, 1], hspace=0.62, wspace=0.5)
sub = gs[2, :].subgridspec(1, 2, width_ratios=[0.1, 1], wspace=0.0)
axes = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1]), fig.add_subplot(sub[0, 1])]
DX = {"a": -0.075, "b": -0.2, "c": -0.42, "d": -0.14}
DY = {"a": 1.22, "b": 1.13, "c": 1.13, "d": 1.10}
for ax, fn, L in zip(axes, [draw_detection, draw_coherence, draw_modcorr, draw_heatmap], "abcd"):
    fn(ax); bs.panel(ax, L, dx=DX[L], dy=DY[L])
for ext in ("pdf", "png"):
    fig.savefig(f"{bs.OUT}/supp/FigS11.{ext}", dpi=400 if ext == "png" else None)
plt.close(fig)
print("saved Fig7")
