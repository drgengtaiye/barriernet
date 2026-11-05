"""Fig. 1b: "2,152 patients" -> "2,144 patients" (ten gastric cohorts after removal of eight shared patients).

The supplied BioRender image is edited at the pixel level, with no new font and no generated content: the digits "52"
are painted over with the panel's flat background (RGB 237, 244, 252) and the bold "4" of "1,496 patients" in the same
panel, at the same size, is copied twice into the freed digit slots (digit advance 36 px, baseline aligned). Only the
box x 1415-1487, y 1901-1951 changes. The authors should make the same change in the BioRender source.

Usage: python3 manuscript/patch_fig1_v33.py  (reads final_v32/figures/Figure1.png; writes figures/cns/v33_fig1/)
"""
import pathlib
import numpy as np
from PIL import Image

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "final_v32/figures/Figure1.png"
OUT = ROOT / "figures/cns/v33_fig1"
OUT.mkdir(parents=True, exist_ok=True)

im = Image.open(SRC)
dpi = im.info.get("dpi", (300, 300))
a = np.array(im.convert("RGB"))
bg = np.array([237, 244, 252], dtype=a.dtype)
assert (a[1880, 1420] == bg).all() and (a[1960, 1490] == bg).all()
a[1898:1953, 1415:1491] = bg                         # erase "52" of "2,152"
glyph = a[1898:1947, 2602:2636].copy()               # bold "4" of "1,496" (ink x 2604-2633, y 1902-1942) with margin
for ink_x in (1419, 1455):                           # ink start of a "4" in the two digit slots
    x0, y0 = ink_x - 2, 1898 + 3                     # +3: baseline of "2,152" is 3 px lower than that of "1,496"
    a[y0:y0 + glyph.shape[0], x0:x0 + glyph.shape[1]] = glyph
changed = (np.array(im.convert("RGB")) != a).any(2)
ys, xs = np.where(changed)
assert ys.min() >= 1898 and ys.max() <= 1952 and xs.min() >= 1415 and xs.max() <= 1490, (ys.min(), ys.max(), xs.min(), xs.max())
out = Image.fromarray(a)
out.save(OUT / "Figure1.png", dpi=dpi)
out.save(OUT / "Figure1.pdf", "PDF", resolution=float(dpi[0]))
print("changed pixels:", int(changed.sum()), "box", ys.min(), ys.max(), xs.min(), xs.max())
