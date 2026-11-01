"""E30 Q1 step 0 (gate): each re-run with the filter switched off (results/E30_dedup_runs/off/) must reproduce the
published results files to |difference| < 1e-6 (numeric cells); text cells must be identical. Null draws (published as
.rds) are compared after export to csv by Rscript. Writes results/E30_dedup_gate.csv. Read-only on published files."""
import subprocess, sys, io
from pathlib import Path
import numpy as np, pandas as pd

A = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
OFF = A / "E30_dedup_runs" / "off"
TOL = 1e-6
FILES = {
    "E30_dedup_E15.R (E15_gastric_validation.R)": ["E15_per_cohort.csv", "E15_pooled.csv", "E15_null.csv", "E15_cindex.csv",
        "E15_deltaC.csv", "E15_cohort_table.csv", "E15_module_coverage.csv", "E15_overlap_check.csv", "E15_null_draws.rds"],
    "E30_dedup_E20.R (E20_extended_replication.R)": ["E20_per_cohort.csv", "E20_pooled.csv", "E20_estimates_used_for_pooling.csv",
        "E20_null.csv", "E20_cindex.csv", "E20_deltaC.csv", "E20_cohort_table.csv", "E20_duplicate_check.csv",
        "E20_duplicate_check_samples.csv", "E20_duplicate_flagged.csv", "E20_duplicate_exclusions.csv",
        "E20_module_coverage.csv", "E20_expectations.csv", "E20_null_draws.rds"],
    "E30_dedup_E29.R (E29_dissociation.R)": ["E29_per_cohort.csv", "E29_pooled.csv", "E29_step0_reproduction.csv", "E29_Q3_crossover.csv"],
    "E30_dedup_E29b.R (E29b_posthoc_adjusted.R)": ["E29b_posthoc_adjusted_per_cohort.csv", "E29b_posthoc_adjusted_pooled.csv"],
}


def read_pub(f):
    if f.endswith(".rds"):
        out = subprocess.run(["Rscript", "-e", f"write.csv(readRDS('{A / f}'), stdout(), row.names=FALSE)"],
                             capture_output=True, text=True, check=True).stdout
        return pd.read_csv(io.StringIO(out))
    return pd.read_csv(A / f)


def compare(pub, new):
    if pub.shape != new.shape or list(pub.columns) != list(new.columns):
        return dict(n_numeric=np.nan, max_abs_diff=np.nan, text_mismatch=np.nan,
                    note=f"shape/columns differ: published {pub.shape}, re-run {new.shape}")
    nmax, ntext, nnum = 0.0, 0, 0
    for c in pub.columns:
        a, b = pub[c], new[c]
        if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
            a, b = a.astype(float).values, b.astype(float).values
            both_na = np.isnan(a) & np.isnan(b)
            one_na = np.isnan(a) ^ np.isnan(b)
            d = np.abs(a - b)[~both_na & ~one_na]
            nnum += (~both_na).sum()
            if one_na.any():
                ntext += int(one_na.sum())
            if d.size:
                nmax = max(nmax, float(np.nanmax(d)))
        else:
            ntext += int((a.fillna("<NA>").astype(str) != b.fillna("<NA>").astype(str)).sum())
    return dict(n_numeric=int(nnum), max_abs_diff=nmax, text_mismatch=ntext, note="")


rows = []
for script, files in FILES.items():
    for f in files:
        newf = OFF / (f.replace(".rds", ".csv"))
        if not newf.exists():
            rows.append(dict(script=script, file=f, n_numeric=np.nan, max_abs_diff=np.nan, text_mismatch=np.nan,
                             note="re-run output missing", passed=False)); continue
        r = compare(read_pub(f), pd.read_csv(newf))
        r.update(script=script, file=f)
        r["passed"] = bool(r["note"] == "" and r["max_abs_diff"] < TOL and r["text_mismatch"] == 0)
        rows.append(r)
G = pd.DataFrame(rows)[["script", "file", "n_numeric", "max_abs_diff", "text_mismatch", "passed", "note"]]
G.to_csv(A / "E30_dedup_gate.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
print(G.to_string(index=False))
print("ALL PASSED" if G.passed.all() else "GATE FAILURES: " + ", ".join(G.file[~G.passed]))
