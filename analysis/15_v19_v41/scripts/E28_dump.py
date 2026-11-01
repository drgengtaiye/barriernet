"""E28 Q1b input: each comparison's log expression (gene symbol x sample) for the official ESTIMATE run in
E28_estimate.R. Missing values are kept; the R script handles them both ways."""
import sys
from E28_base import TR, ALL6
out = sys.argv[1]
for k in ALL6:
    X = TR[k]["X"]; X.index.name = "GeneSymbol"
    X.to_csv(f"{out}/{k}.tsv", sep="\t")
    print(k, X.shape, int(X.isna().sum().sum()), "missing")
