# E45 outcome: a reference gene frozen in single-arm data, a coexpression-matched random reference, the renal arms

Protocol E45_prereg.md (sha256 17135ddd…, logged 2026-11-05T15:39:08Z). Frozen genes results/E45_reference_genes.tsv (sha256 2ec89119…,
logged 2026-11-05T15:40:49Z). Eighteenth protocol written after the other results were known. Step 0 passed: 48 barriernet genes against
E37 (5.3e-15); CD8A per single-arm cohort against E40 (0); signatures minus CD8A against E32 (9.6e-17); CD8A per gastroesophageal
cohort against E44 (5.6e-17). 7 of 8 expectations held (A3 not counted: G2 is the same gene as G1).

## A. Reference gene

G1 = G2 = CXCL10: +0.470 (0.263 to 0.677) in the 20 single-arm cohorts; 1,996 of 2,082 candidates eligible. CD8A ranked 19th (+0.353).

| set | CXCL10 interaction | CD8A interaction | signatures above / below CXCL10 | CXCL10 minus CD8A |
|---|---|---|---|---|
| six | +0.171 (−0.044 to 0.386), P = 0.096; 95th percentile, rank 94 of 2,082 | +0.226 (0.050 to 0.403) | 8 / 2 | −0.134 (−0.253 to −0.015), P = 0.027 |
| renal three | +0.108 (−0.301 to 0.517) | +0.269 (−0.192 to 0.729) | 12 / 1 | −0.166 (−0.288 to −0.043), P = 0.008 |
| breast three | +0.451 (0.028 to 0.875), P = 0.043 | +0.060 (−0.331 to 0.451) | 0 / 9 | +0.368 (−0.123 to 0.858), P = 0.14 |

A2 failed: 8 of the 52 significantly exceeded CXCL10. "No signature significantly exceeded a single gene" holds for CD8A and most T-cell
genes, which were analyzed with knowledge of these data, and does not hold for the gene that single-arm data alone would have chosen.

Gastroesophageal cohorts (not used for the selection): CXCL10 +0.887 (0.387 to 1.386), P = 0.0005, fixed effect; main-text rule P = 0.12;
minus CD8A +0.331 (−0.280 to 0.941), P = 0.29. Two gastric cohorts: +0.877 unadjusted, +0.615 (−0.020 to 1.250), P = 0.058, adjusted for
MSI-H or EBV-positive status (CD8A +0.440, P = 0.19).

## B. Coexpression-matched random reference

All 287 signature × comparison tasks reached 1,000 sets at ±0.02; 234 of 287,000 fits failed. 16 of 52 exceeded the reference (10 of 20
inflamed, 6 of 24 other, 0 of 8 confounder), against 26 (15, 11, 0) with the size-matched reference. The median of the reference medians
moved from −0.059 to −0.017. No signature was gained. Ten were lost, among them STAT1_Care and Teff-IFNG_Fehrenbacher.

## C. Renal arms and vascular genes

CD8A, benefit-oriented slopes, checkpoint / control: IMmotion150 +0.289 / −0.238; CheckMate 025 +0.124 / −0.067; JAVELIN Renal 101 +0.157 /
−0.019; renal pool +0.185 (−0.017 to 0.388) / −0.065 (−0.292 to 0.162), P = 0.34. C1 held by direction only. In IMmotion150, the control-arm
slope accounts for a large part of the interaction; in the other two, most of it comes from the checkpoint arm.

Angio score, interaction: renal −0.267 (−0.560 to 0.026), P = 0.059; breast −0.012 (−0.750 to 0.725). Control-arm slopes +0.594 (IMmotion150),
+0.094 (CheckMate 025, everolimus), +0.320 (JAVELIN Renal 101): high angiogenesis goes with better outcome on sunitinib.

TUMOR-VASCULATURE-UP_Lu: 12 of its genes among the 2,082. Control-favoring 3 of 102 (COL5A3, PAMR1, PLXDC1), odds ratio 6.6, P = 0.018;
checkpoint-favoring 1 of 101 (ISG15), P = 0.45. The enrichment rests on three genes.

## Consequence for the manuscript

The headline changes from "no signature exceeded single T-cell genes" to: signatures did not exceed CD8A or most T-cell genes, but CD8A was
chosen with knowledge of the data, and 8 exceeded the gene chosen blind to the randomized data, in the renal comparisons and not the breast
comparisons. Which single gene matches the signatures cannot be read from single-arm data.

## Notes

E45_select.py was rerun once after the freeze log to correct a label in its log (it printed the hash variable of the E34 setup). The reference file and its hash did not change. E45_expectations.csv prints counts as 8.0 and 2.0 because the column held missing values; the workbook shows 8 and 2.
