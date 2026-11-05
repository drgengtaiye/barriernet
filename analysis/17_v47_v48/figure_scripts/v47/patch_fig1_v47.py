"""Fig. 1d: the sentence below the arrows said the two extensions were "pre-registered", which they were not (dated, hash-logged
protocols, not registered publicly). The two text lines are painted over with the panel's flat background and rewritten in the
same size. Only that box changes. The authors should make the same change in the BioRender source.
Reads figures/cns/v33_fig1/Figure1.png; writes figures/cns/v47/Figure1.png and .pdf."""
import pathlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont
ROOT = pathlib.Path(__file__).resolve().parents[3]
SRC = ROOT / "figures/cns/v33_fig1/Figure1.png"; OUT = pathlib.Path(__file__).parent
im = Image.open(SRC).convert("RGB"); dpi = Image.open(SRC).info.get("dpi", (300, 300)); a = np.array(im)
Y0, Y1, X0, X1 = 2765, 2895, 2125, 3730
bg = np.array([237, 244, 252], dtype=a.dtype); assert (a[Y0 - 4, X0] == bg).all() and (a[Y1 + 4, X1] == bg).all()
a[Y0:Y1, X0:X1] = bg; im2 = Image.fromarray(a); d = ImageDraw.Draw(im2)
f = ImageFont.truetype("/users/yc4229/.fonts/LiberationSans-Regular.ttf", 50)
INK = (20, 20, 20)
d.text((2150, 2765), "Exceptions are stated in the text: two extensions followed dated", font=f, fill=INK)
d.text((2150, 2823), "protocols after P = 0.057; signature analyses are post hoc.", font=f, fill=INK)
changed = (np.array(Image.open(SRC).convert("RGB")) != np.array(im2)).any(2); ys, xs = np.where(changed)
assert ys.min() >= Y0 and ys.max() < Y1 and xs.min() >= X0 and xs.max() < X1, (ys.min(), ys.max(), xs.min(), xs.max())
im2.save(OUT / "Figure1.png", dpi=dpi); im2.save(OUT / "Figure1.pdf", "PDF", resolution=float(dpi[0]))
print("changed pixels", int(changed.sum()), "box", ys.min(), ys.max(), xs.min(), xs.max())
