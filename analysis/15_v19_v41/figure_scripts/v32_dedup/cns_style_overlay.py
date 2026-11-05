"""House style for the v19 main figures (Nature/Cell/Science layout).

Physical size is final size: 180 mm double column (W2), 120 mm 1.5 column (W15), 88 mm single
column (W1); at most 210 mm tall. Type is Helvetica 5-7 pt, panel letters 8 pt bold lower case,
axes 0.5 pt, data lines 0.75-1 pt. Nothing is scaled after rendering, so what is drawn here is
what the journal prints.

Colour does one job per figure:
  identity  - module families: STROMA (CAF_ECM, Barrier), MACRO (TAM_SUP), IMMUNE (APM, IFN,
              TLS_B, T_CELL, Immune, NetImmune). Checked with ../check_palette.py (port of the
              dataviz validator): IMMUNE-STROMA and IMMUNE-MACRO pass; STROMA-MACRO sits at the
              CVD floor (7.1), so wherever both appear they also carry a direct label.
  polarity  - DIV: favorable (blue) <- neutral PAPER -> unfavorable (vermilion). "Favorable"
              means longer survival, response, or benefit from the checkpoint arm, so every
              heat map in the set reads the same way.
  magnitude - SEQ: one hue (blue), light -> dark.
Text is always INK or SLATE, never a series colour; colour sits on a mark beside the text.
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric")
import figstyle as fs
from figstyle import fmt_p  # noqa: F401  (Nature P-value format, shared with the rest of the paper)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

MM = 1 / 25.4
W2, W15, W1 = 180 * MM, 120 * MM, 88 * MM
HMAX = 210 * MM

ROOT = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2")
OUT = Path("/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-com-gastric-2/727be15e-8d67-407e-a66a-a54d26cc1ebb/scratchpad/v32fig/out")
RES = Path("/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-com-gastric-2/727be15e-8d67-407e-a66a-a54d26cc1ebb/scratchpad/v32fig/res")

INK, SLATE, GREY, LIGHT, PAPER = "#1A1F24", "#3F4C5A", "#8E9AA6", "#D8DEE4", "#F4F6F8"
IMMUNE, STROMA, MACRO = "#0F6FC5", "#0B6B47", "#C96F12"   # STROMA darkened: protan ΔE vs MACRO 7.1 -> 11.0
FITTED = "#3F4C5A"
FAV, UNFAV = "#0F6FC5", "#C1401A"
IMMUNE_LIGHT, STROMA_LIGHT, MACRO_LIGHT = "#E7F0FA", "#E6F3EE", "#FAF0E6"   # panel washes

MODCOL = {"CAF_ECM": STROMA, "Barrier": STROMA, "TAM_SUP": MACRO, "C1Q subset": MACRO,
          "SPP1 subset": MACRO, "CXCL9:SPP1": MACRO, "APM": IMMUNE, "IFN": IMMUNE, "TLS_B": IMMUNE,
          "T_CELL": IMMUNE, "Immune": IMMUNE, "NetImmune": IMMUNE}
MODULES = ["CAF_ECM", "TAM_SUP", "APM", "IFN", "TLS_B", "T_CELL"]

DIV = LinearSegmentedColormap.from_list(
    "fav_unfav", ["#0A4C8A", FAV, "#9CC3E6", PAPER, "#E9A890", UNFAV, "#8A2A10"])
SEQ = LinearSegmentedColormap.from_list("seq", [PAPER, "#9CC3E6", FAV, "#0A3D6E"])

# cell classes (shared by the single-cell and spatial figures)
CT = {**fs.CT, "Fibro": STROMA, "Endo": "#2F4858"}   # fibroblasts in the stroma green (as in Figs 1, 7); endothelium re-stepped to stay separable from it
CTORDER = fs.CTORDER

TUMCOL = {"melanoma": "#2E6FA8", "renal": "#C1401A", "kidney": "#C1401A", "NSCLC": "#1B7F79",
          "lung": "#1B7F79", "urothelial": "#4E6472", "bladder": "#4E6472", "gastric": "#E08A2E",
          "breast": "#B0457D", "HNSCC": "#7A8B99", "mixed": "#A9B6C1", "esophageal": "#C96F12"}


def apply():
    fs.apply()
    plt.rcParams.update({
        "font.size": 6.0, "axes.labelsize": 6.0, "axes.titlesize": 6.5,
        "xtick.labelsize": 5.6, "ytick.labelsize": 5.6, "legend.fontsize": 5.6,
        "axes.titlepad": 3.0, "axes.labelpad": 2.0,
        "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
        "lines.linewidth": 0.8, "lines.markersize": 3.0, "hatch.linewidth": 0.4,
        "savefig.bbox": "standard", "savefig.pad_inches": 0.0,
    })


apply()


def letter(fig, ax, s, dx_mm=6.0, dy_mm=1.5):
    """Panel letter a fixed physical distance left of / above the axes box, so letters in a
    row align whatever the width of each panel's tick labels."""
    bb = ax.get_position()
    W, H = fig.get_size_inches()
    fig.text(bb.x0 - dx_mm * MM / W, bb.y1 + dy_mm * MM / H, s, fontsize=8, fontweight="bold",
             va="bottom", ha="left", color=INK)


def letter_at(fig, x_mm, y_mm_from_top, s):
    """Panel letter at an absolute position (mm from the left edge and from the top)."""
    W, H = fig.get_size_inches()
    fig.text(x_mm * MM / W, 1 - y_mm_from_top * MM / H, s, fontsize=8, fontweight="bold",
             va="top", ha="left", color=INK)


def ax_mm(fig, x, y, w, h, **kw):
    """Axes placed in millimetres from the top-left corner of the figure (x, y, width, height)."""
    W, H = fig.get_size_inches()
    W, H = W / MM, H / MM
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H], **kw)


def title(ax, s, sub=None):
    """Short bold panel title, optional muted subtitle on the line below."""
    ax.set_title(s, loc="left", fontsize=6.5, fontweight="bold", color=INK,
                 pad=9 if sub else 3)
    if sub:
        ax.text(0, 1.0, sub, transform=ax.transAxes, fontsize=5.4, color=SLATE, va="bottom",
                ha="left")


def clean(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def rule(ax, x, **kw):
    """Reference line (null effect), dashed and recessive."""
    ax.axvline(x, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0, **kw)


def hrule(ax, y, **kw):
    ax.axhline(y, color=GREY, lw=0.5, ls=(0, (2, 2)), zorder=0, **kw)


def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.png", dpi=600)
    plt.close(fig)
    return OUT / f"{stem}.png"
