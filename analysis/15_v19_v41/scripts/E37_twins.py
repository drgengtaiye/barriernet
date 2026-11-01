"""E37 step 2 - expression-twin check of the kidney and breast development cohorts against the randomized comparisons.

Protocol: analysis/E37_prereg.md. Rule as E12 (E12a_cohort_audit.R cross_cor): for each development sample, Spearman
correlation with every sample of a randomized comparison of the same cancer type over the shared genes (>= 200); the
sample is flagged when its best partner stands more than 0.03 clear of the second best. Reads no outcome.
Writes results/E37_twins.csv (one row per development cohort x comparison) and results/E37_twin_flags.csv.
"""
import hashlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h = hashlib.sha256(open(f"{A}/E37_prereg.md", "rb").read()).hexdigest()
assert _h == open(f"{A}/E37_prereg.sha256").read().split()[0], "E37 protocol hash mismatch"
_src = open(f"{A}/E28_base.py").read().split("# ------------------------------------------------------------------ data (copied from E5, E10, E21)")
exec(compile(_src[0], f"{A}/E28_base.py", "exec"))
_data = _src[1].split("def ok_mask(d, k):")[0]                         # data loading only (no outcome model is fitted)
exec(compile(_data, f"{A}/E28_base.py", "exec"))
from scipy.stats import rankdata
print("E37 protocol sha256 verified:", _h[:12])

PAIRS = {"ICB_Miao1": ["IMM150", "CM025", "JR101"], "ICB_Shiuan": ["IMM150", "CM025", "JR101"],
         "brca_iatlas_anders_2022": ["ISPY2_pembro", "ISPY2_durva", "NeoTRIP"]}


def ranks(M):
    return np.apply_along_axis(rankdata, 0, M)


rows, flags = [], []
for dev, cmps in PAIRS.items():
    Dm = pd.read_csv(f"{O}/E37_twin_{dev}.csv.gz", index_col=0)
    Dm = Dm[~Dm.index.duplicated()]
    for k in cmps:
        X = TR[k]["X"]; X = X[~X.index.duplicated()]
        g = Dm.index.intersection(X.index)
        Ra, Rb = ranks(Dm.loc[g].to_numpy(float)), ranks(X.loc[g].to_numpy(float))
        Ra = (Ra - Ra.mean(0)) / Ra.std(0); Rb = (Rb - Rb.mean(0)) / Rb.std(0)
        C = Ra.T @ Rb / len(g)                                        # Spearman = Pearson of ranks
        srt = -np.sort(-C, axis=1); best, gap = srt[:, 0], srt[:, 0] - srt[:, 1]
        rows.append(dict(development=dev, comparison=k, shared_genes=len(g), n_dev=C.shape[0], n_trial=C.shape[1],
                         median_best_r=np.median(best), max_best_r=best.max(), median_gap=np.median(gap),
                         max_gap=gap.max(), n_flagged=int((gap > 0.03).sum())))
        for i in np.where(gap > 0.03)[0]:
            flags.append(dict(development=dev, sample=Dm.columns[i], comparison=k, best_r=best[i], gap=gap[i],
                              partner=X.columns[int(np.argmax(C[i]))]))
R = pd.DataFrame(rows); R.to_csv(f"{O}/E37_twins.csv", index=False)
F = pd.DataFrame(flags, columns=["development", "sample", "comparison", "best_r", "gap", "partner"])
F.to_csv(f"{O}/E37_twin_flags.csv", index=False)
pd.set_option("display.width", 200)
print(R.round(3).to_string(index=False)); print(f"\nflagged samples: {len(F)}"); print(F.to_string(index=False))
