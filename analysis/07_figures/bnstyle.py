"""Shared style and assembly for the BarrierNet v14 figures.

Type, rules and sizes come from the companion manuscript's figstyle.py so the two papers
read alike. Category colours were re-checked with check_palette.py (a port of the dataviz
validator): immune/stroma/macrophage/fitted pass CVD and normal-vision separation; the
stroma–macrophage pair sits at the CVD floor (7.1), so wherever both appear they are also
separated by position and a direct label.
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric")
import figstyle as fs
from figstyle import MM, W1, W2, INK, LIGHT, GREY, PAPER, panel, title, fmt_p  # noqa: F401
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

IMMUNE = "#0F6FC5"   # NetImmune, immune modules
STROMA = "#138A5E"   # CAF_ECM, Barrier
MACRO = "#C96F12"    # TAM_SUP
FITTED = "#3F4C5A"   # learned / full-transcriptome models
NULL = LIGHT         # random-gene null distributions
OUT = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures"

fs.apply()


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/{stem}.{ext}", format=ext, dpi=400 if ext == "png" else None)
    plt.close(fig)


def load_image(path, dpi=400):
    """PNG/JPG as an array; PDF rasterised (first page) at `dpi`."""
    if path.lower().endswith(".pdf"):
        import fitz, numpy as np
        pg = fitz.open(path)[0]
        pix = pg.get_pixmap(dpi=dpi, alpha=False)
        return np.frombuffer(pix.samples, dtype="uint8").reshape(pix.height, pix.width, pix.n)
    return mpimg.imread(path)


def compose(stem, cells, grid, height_ratios=None, width_ratios=None, width=W2, height=None,
            hspace=0.12, wspace=0.08):
    """Assemble a multi-panel figure.

    cells: list of dicts with
      pos    - (row0, row1, col0, col1) slice in the grid (end-exclusive)
      letter - panel letter ('' for none)
      path   - image file to place (aspect preserved), or
      draw   - function(ax) that draws a data panel in the house style
    """
    nr, nc = grid
    if height is None:
        height = width * nr / nc * 0.8
    fig = plt.figure(figsize=(width, height))
    gs = GridSpec(nr, nc, figure=fig, height_ratios=height_ratios, width_ratios=width_ratios,
                  hspace=hspace, wspace=wspace)
    for c in cells:
        r0, r1, c0, c1 = c["pos"]
        ax = fig.add_subplot(gs[r0:r1, c0:c1])
        if "path" in c:
            ax.imshow(load_image(c["path"]), interpolation="lanczos")
            ax.set_axis_off()
            if c.get("letter"):
                ax.text(0.0, 1.0, c["letter"], transform=ax.transAxes, fontsize=8.2,
                        fontweight="bold", va="bottom", ha="left", color=INK)
        else:
            c["draw"](ax)
            if c.get("letter"):
                panel(ax, c["letter"])
    save(fig, stem)
    return f"{OUT}/{stem}.png"
