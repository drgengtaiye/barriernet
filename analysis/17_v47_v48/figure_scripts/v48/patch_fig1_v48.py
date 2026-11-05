"""Fig. 1, v48: three panel titles rewritten so that panel c is not read as the three references of Fig. 2a and panels b and d no longer
describe the study as one fixed score. Only those title boxes change. The authors should make the same change in the BioRender source.
Reads figures/cns/v47/Figure1.png; writes figures/cns/v48/Figure1.png and .pdf."""
import pathlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont
H = pathlib.Path(__file__).resolve().parent; SRC = H.parent / "v47/Figure1.png"
im = Image.open(SRC).convert("RGB"); dpi = Image.open(SRC).info.get("dpi", (300, 300)); a = np.array(im); bg = np.array([237, 244, 252], dtype=a.dtype)
BOXES = [((60, 1362, 960, 1440), (70, 1370), "b. Datasets"), ((55, 2118, 600, 2190), (70, 2118), "c. Three kinds of evidence"), ((2135, 2120, 2820, 2192), (2147, 2120), "d. What was fixed, and when")]
for (x0, y0, x1, y1), _, _ in BOXES:
    assert (a[y0 - 3, x0] == bg).all() and (a[y1 + 3, x1] == bg).all(), (x0, y0, a[y0 - 3, x0], a[y1 + 3, x1]); a[y0:y1, x0:x1] = bg
im2 = Image.fromarray(a); d = ImageDraw.Draw(im2); f = ImageFont.truetype("/users/yc4229/.fonts/LiberationSans-Bold.ttf", 58)
for _, xy, t in BOXES:
    d.text(xy, t, font=f, fill=(0, 0, 0))
ch = (np.array(Image.open(SRC).convert("RGB")) != np.array(im2)).any(2)
mask = np.zeros(ch.shape, bool)
for (x0, y0, x1, y1), _, _ in BOXES:
    mask[y0:y1 + 10, x0:x1 + 400] = True
assert not (ch & ~mask).any(), "pixels changed outside the title boxes"
im2.save(H / "Figure1.png", dpi=dpi); im2.save(H / "Figure1.pdf", "PDF", resolution=float(dpi[0])); print("ok; changed pixels", int(ch.sum()))
