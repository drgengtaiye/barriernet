"""Row-packed figure layout for panels of mixed aspect ratio.

Each row is a list of items; every item has an `aspect` (width / height of its cell). A
row's height is chosen so that its cells fill the figure width, which keeps images at their
own aspect without the white margins a uniform grid would leave. An item is either an image
(`path`, optionally `strip=True` to remove an old baked-in panel letter) or a data panel
(`draw(ax)`, drawn inside the cell with `m` = (left, bottom, right, top) margins in mm for
tick labels). An image item may carry `crop` = (top, bottom) fractions to drop an old
panel title; its aspect is then taken from the cropped image unless given. Letters are
placed identically above every cell's top-left corner.
"""
import matplotlib.pyplot as plt
from bnstyle import MM, INK, load_image, save
from panelprep import strip_letter


def _image(it):
    img = strip_letter(it["path"]) if it.get("strip") else load_image(it["path"])
    t, b = it.get("crop", (0, 0))
    h = img.shape[0]
    return img[int(h * t): h - int(h * b)]


def pack(stem, rows, width_mm=180, hgap=3.0, vgap=2.5, letter_pad=3.2, row_scale=None):
    row_scale = row_scale or [1.0] * len(rows)
    for row in rows:
        for it in row:
            if "path" in it:
                it["_img"] = _image(it)
                it.setdefault("aspect", it["_img"].shape[1] / it["_img"].shape[0])
    geo, total = [], 0.0
    for row, sc in zip(rows, row_scale):
        avail = width_mm - hgap * (len(row) - 1)
        h = avail / sum(it["aspect"] for it in row) * sc
        geo.append(h)
        total += h + letter_pad
    total += vgap * (len(rows) - 1)
    fig = plt.figure(figsize=(width_mm * MM, total * MM))
    y = total
    for row, h in zip(rows, geo):
        y -= letter_pad
        x = 0.0
        for it in row:
            w = it["aspect"] * h
            if "path" in it:
                img = it["_img"]
                ax = fig.add_axes([x / width_mm, (y - h) / total, w / width_mm, h / total])
                ax.imshow(img, interpolation="lanczos", aspect="equal")
                ax.set_anchor("NW")
                ax.set_axis_off()
            else:
                l, b, r, t = it.get("m", (8, 7, 2, 3))
                ax = fig.add_axes([(x + l) / width_mm, (y - h + b) / total,
                                   (w - l - r) / width_mm, (h - b - t) / total])
                it["draw"](ax)
            if it.get("letter"):
                fig.text(x / width_mm, (y + 0.6) / total, it["letter"], fontsize=8.2,
                         fontweight="bold", ha="left", va="bottom", color=INK)
            x += w + hgap
        y -= h + vgap
    save(fig, stem)
    return total
