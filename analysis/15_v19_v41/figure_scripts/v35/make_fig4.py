"""v35 Figure 4 - barriernet's immune and fibroblast components in the randomized comparisons: panels a-e of v34 Fig. 4
(= v33 Fig. 5; script figure_scripts/make_fig6.py, v31 numbering), cropped above panel f. Panels f (leave-one-out) and
g (random-gene reference) are in Supplementary Fig. S32, the full v34 figure. Writes only into figures/cns/v35/."""
from pathlib import Path
from PIL import Image
import matplotlib.pyplot as plt
ROOT = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2")
HERE = ROOT / "figures/cns/v35"
im = Image.open(ROOT / "final_v34/figures/Figure4.png").convert("RGB")
CUT = 3890                                         # inside the blank band 3860-3919 between panels c-e and f
crop = im.crop((0, 0, im.width, CUT)); crop.save(HERE / "Fig4.png", dpi=(600, 600))
w_mm = 180; h_mm = w_mm * crop.height / crop.width
fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4)); ax = fig.add_axes([0, 0, 1, 1]); ax.imshow(crop); ax.axis("off")
fig.savefig(HERE / "Fig4.pdf", dpi=600); print("Fig4", crop.size, f"{h_mm:.0f} mm tall")
