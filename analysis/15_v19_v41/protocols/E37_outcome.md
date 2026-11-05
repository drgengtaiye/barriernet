# E37 outcome: barriernet 0.3

Protocol E37_prereg.md (sha256 89c5ce0f9c33…, logged 2026-10-31T02:59:12Z). Frozen weights
results/E37_frozen_weights.tsv (sha256 500264cc627b…, logged 2026-10-31T03:04:28Z, before E37_test.py was written).
Registration for controlled-access trials E37_registration.md (sha256 c24bb4141a73…, logged 2026-10-31T03:04:42Z).

## Development (single-arm cohorts only)

- E37_export.R: 20 cohorts (21 minus ICB_Braun), 1,198 patients, 369 responders; step 0 re-derived 78 published
  per-cohort log-odds (max difference 3.8e-15).
- E37_twins.py: 3 ICB_Miao1 samples were flagged against IMmotion150 (gap 0.033–0.040), all with the same partner
  (p77957) and with low correlations (0.66–0.74, against a median best of 0.86), so they are unlikely to be the same
  tumors; they were removed as the protocol requires. No breast or other kidney sample was flagged. Development set:
  1,195 patients.
- E37_develop.R, leave-one-cohort-out pooled log-odds of response per SD (20 held-out cohorts; NetImmune 19):
  NetImmune 0.2.1 +0.434; Immune axis +0.462; C2 module weights +0.488; C4 gene ridge +0.505; C3 gene elastic net
  +0.488. C4 was highest but within 0.02 of the simpler C2, so C2 was selected (protocol rule).
- Frozen barriernet 0.3 (C2 refitted on all 20; module totals): IFN +0.362, APM +0.225, TLS_B +0.126, T_CELL +0.104,
  CAF_ECM −0.160, TAM_SUP −0.023.
- D1 held (largest learned candidate − Immune axis = +0.043 ≤ 0.10); D2 held (CAF_ECM coefficient −0.106).

## Test in the randomized comparisons: 4 of 4 expectations held

| Q | Result | Held |
|---|---|---|
| Q1 | six-comparison interaction +0.252 (+0.013 to +0.490); P = 0.042 | Yes |
| Q2 | difference from CD8A −0.005 (−0.093 to +0.084); P = 0.92 | Yes |
| Q3 | difference from NetImmune 0.2.1 −0.016 (−0.104 to +0.071); P = 0.71 | Yes |
| Q4 | random reference median −0.045, 95th percentile +0.073; empirical P = 0.001 | Yes |

Q5 (descriptive): positive in all seven comparisons (+0.04 CheckMate 025 to +0.96 I-SPY2 pembrolizumab). Renal three
+0.19 (P = 0.18); breast three +0.57 (P = 0.069), above CD8A by +0.43 (+0.08 to +0.79; P = 0.017); seven comparisons
+0.24 (+0.05 to +0.43; P = 0.022); four chemotherapy-backbone comparisons +0.47 (P = 0.043), above CD8A by +0.36
(+0.06 to +0.67; P = 0.019). NeoPembrOv alone +0.24 (SE 0.34). Step 0: the NetImmune − CD8A paired difference
reproduced E32 (+0.025502).

## Addition beyond the protocol (written after the test results were seen)

E37_panel_extra.py: barriernet 0.3 against each single gene of the two panels by the paired bootstrap (the framework
compares every score with the whole panel; the protocol named only CD8A). Six comparisons: significantly above 3 of 27
genes (CD247, CXCL10, HAVCR2) and below none. Breast three: above 15 of 26 and below none; not significantly different
from the strongest breast genes (GZMB +0.06, P = 0.72; IDO1 +0.05, P = 0.80; CXCL10 +0.10, P = 0.59; STAT1 +0.15,
P = 0.28).

## Software

barriernet 0.3.0 (repo_barriernet): bundled weights (data/geneset_rna_v03.yaml; column barriernet03 in
`barriernet bulk`), module barriernet.calibrate (Efron Cox, logistic with Firth fallback, the conservative pooling
rule, size-matched random reference, within-arm paired bootstrap, detectable difference) and `barriernet calibrate`.
23 tests pass (Cox matches lifelines with ties to 1e-6). E37_package_check.py: the package reproduces the analysis
code's barriernet 0.3 interactions in all seven comparisons and the pooled value (largest difference 5.2e-6, from
lifelines' default precision in the analysis code).

## Reading

The frozen score, developed without randomized data, predicted benefit in them and was specific to its genes, but it
was no better than CD8A or than NetImmune 0.2.1 over the six comparisons; its gain was in the chemotherapy-backbone
comparisons, where CD8A is weak, and even there it did not exceed the best single genes. Its development gain in
single-arm cohorts (+0.03 to +0.05 log-odds over the fixed axes) did not carry into a measurable gain in the trials.

Note (2026-11-01): the step-0 table of E37_export.R (78 per-cohort log-odds) was saved as results/E37_step0_reproduction.csv and later overwritten by E37_test.py, which writes a file of the same name (now the trial step 0). The export step 0 was re-run unchanged into results/E37_export_step0_rerun.csv by E37_export_step0_rerun.R: 78 values, largest difference 3.8 × 10⁻¹⁵, as reported.
