"""Check that the barriernet 0.3.0 package (barriernet.calibrate, bundled v0.3 weights) reproduces the analysis code:
barriernet 0.3 interactions per randomized comparison (E37_test_per_comparison.csv) and their pooled six-comparison
value (E37_test_pooled.csv). Setup as E37_test.py (E34 setup). Writes results/E37_package_check.csv."""
import sys
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E37p_step0_reproduction").replace("E34_coverage_check.csv", "E37p_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/repo_barriernet/barriernet")
import barriernet
from barriernet.calibrate import Trial, interaction, pool as bpool, weighted_score
from barriernet.genesets import load_weights
wt = load_weights()["weights"]; prev = pd.read_csv(f"{O}/E37_test_per_comparison.csv")
rows = []
for k in ALL6 + ["NeoPembrOV"]:
    d = TR[k]
    tr = Trial(k, d["Z"], d["arm"], time=d.get("t") if d["kind"] == "cox" else None, event=d.get("e") if d["kind"] == "cox" else None,
               response=d.get("y") if d["kind"] != "cox" else None, covariates=d.get("cov"), firth=d.get("firth", True))
    sc, g, w = weighted_score(d["Z"], wt)
    y, se = interaction(tr, sc)
    p = prev[(prev.comparison == k) & (prev.score == "BN03")].iloc[0]
    rows.append(dict(comparison=k, package=y, analysis=p.y, abs_diff=abs(y - p.y), se_package=se, se_analysis=p.se, genes=len(g)))
C = pd.DataFrame(rows)
q = bpool(C.package[:6], C.se_package[:6]); P = pd.read_csv(f"{O}/E37_test_pooled.csv")
ref = P[(P.set == "six (primary)") & (P.score == "barriernet 0.3")].iloc[0]
C.loc[len(C)] = dict(comparison="pooled six", package=q["est"], analysis=ref.est, abs_diff=abs(q["est"] - ref.est), se_package=np.nan, se_analysis=np.nan, genes=np.nan)
C.to_csv(f"{O}/E37_package_check.csv", index=False)
print("barriernet", barriernet.__version__); print(C.to_string(index=False))
print("max |difference|:", C.abs_diff.max()); assert C.abs_diff.max() < 1e-5
