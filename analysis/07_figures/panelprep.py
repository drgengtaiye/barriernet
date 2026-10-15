"""Remove the bold panel letter baked into an old panel image (top-left corner).

The tallest near-black, near-neutral connected component in the top-left box (x < 20% of
width, y < 25% of height) that is letter-shaped (height >= 3% of the image, not a thin
rule) is taken as the old letter and painted white. Check the result visually. `strip_letter(path, show=True)` returns the box for visual checking.
"""
import numpy as np
from PIL import Image
from scipy import ndimage

def strip_letter(path, xfrac=0.2, yfrac=0.25, min_h=0.03, return_box=False):
    im = np.array(Image.open(path).convert("RGB")).copy()
    H, W, _ = im.shape
    reg = im[: int(H * yfrac), : int(W * xfrac)]
    reg16 = reg.astype(int)
    # near-neutral near-black only: letters are black, heatmap cells and dots are coloured
    dark = (reg16.mean(axis=2) < 70) & (reg16.max(axis=2) - reg16.min(axis=2) < 30)
    lab, n = ndimage.label(dark)
    best = None
    for sl in ndimage.find_objects(lab):
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        fill = (lab[sl] > 0).mean()
        if (h >= min_h * H and 0.005 * W <= w <= 0.06 * W and 0.7 <= h / w <= 4.5
                and fill > 0.15 and (best is None or h > best[0])):
            best = (h, sl)
    box = None
    if best:
        sl = best[1]
        pad = 4
        y0, y1 = max(0, sl[0].start - pad), sl[0].stop + pad
        x0, x1 = max(0, sl[1].start - pad), sl[1].stop + pad
        im[y0:y1, x0:x1] = 255
        box = (x0, y0, x1, y1)
    return (im, box) if return_box else im
