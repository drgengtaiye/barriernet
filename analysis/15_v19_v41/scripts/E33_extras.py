"""E33 add-on: breast-only pools of the barriernet and comparator scores, and per-comparison CD8A and coverage for the
main-text table of trials. Re-pools E32 per-comparison values with the manuscript rule (E28_base.pool_rule); no new fits."""
import re, numpy as np, pandas as pd
from scipy.stats import norm, t as tdist
from scipy.optimize import minimize_scalar
src = open("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/E28_base.py").read()
exec(src[src.index("def _reml_tau2"):src.index("# ------------------------------------------------------------------ data")])
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
PER = pd.read_csv(f"{O}/E32_per_comparison.csv"); PER = PER[PER.variant == "primary"]
BREAST = ["ISPY2_pembro", "ISPY2_durva", "NeoTRIP"]
rows = []
for s in ["NetImmune", "Immune", "CAF_ECM", "TAM_SUP", "Barrier", "ESTIMATE", "TMEscore_z", "CD8A"]:
    r = PER[(PER.signature == s) & PER.comparison.isin(BREAST)]
    q = pool_rule(r.y.values, r.se.values)
    rows.append(dict(set="breast three", name=s, est=q["est"], lo=q["lo"], hi=q["hi"], p=q["p"], k=len(r), k_pos=int((r.y > 0).sum())))
B = pd.DataFrame(rows); B.to_csv(f"{O}/E33_breast_extras.csv", index=False); print(B.round(3).to_string(index=False))
cov = pd.read_csv(f"{O}/E32_coverage.csv")
T = cov.groupby("comparison").scored.sum().rename("signatures_scored").to_frame()
c8 = PER[PER.signature == "CD8A"].set_index("comparison")[["y", "se", "p"]]
T = T.join(c8); T.to_csv(f"{O}/E33_trial_table.csv"); print(T.round(3))
