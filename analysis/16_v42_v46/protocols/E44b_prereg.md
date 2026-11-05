# E44b addendum (post hoc): leave-one-out, power and heterogeneity for the gastric analyses

Written 2026-11-03, after the E44 results and after an independent fact-check of the v46 text, and before
E44b_checks.py is written or run. The hash of this file is logged in E44b_prereg.sha256 first. This addendum has no
expectations and is not counted in the tally of later protocols; the E44 results are reported as pre-specified.

## Why

The fact-check found three statements in the v46 draft that the E44 outputs do not support as written:

- "PRJEB25780 carries most of the weight": by inverse variance it has the smallest weight of the three cohorts; its
  influence comes from its large estimates.
- "Under the conservative rule nothing is significant": CD8A and barriernet 0.3 are not, but two signatures
  (CIN25_Carter, proliferation_Pabla) stay at P < 0.05 in the three-cohort set and CTLA4 in the two-gastric set.
- "The gastric trials could not have detected an interaction of CD8A's size": they had less than 80% power for it,
  which is not the same.

It also found that the file records of the project date the first PRJEB25780 analysis (2026-03-24) after barriernet
0.2.0 fixed the 48 genes (2026-02-11), against the E44 protocol's statement that PRJEB25780 was among the cohorts
explored when the genes were fixed. This addendum records the file dates; it does not settle the question, which the
author is asked to confirm.

## What is computed (descriptive, no expectations)

From results/E44_B_per_cohort.csv, E44_B_pooled.csv and E44_A_gastric_trials.csv, nothing refitted:

1. For CD8A, barriernet 0.3, the Immune axis and TMEscore: inverse-variance weight of each cohort in the three-cohort
   fixed-effect pool; its share of the pooled estimate (weight × estimate / sum); the fixed-effect pool leaving out each
   cohort in turn (estimate, 95% CI, P); Cochran's Q and I² over the three.
2. Signatures (of the 52) at P < 0.05 under the main-text rule (cons_p) in each set, and the composites and genes.
3. Power of each gastric trial, and of the two combined, for a per-SD interaction of +0.2264 (CD8A, six comparisons;
   E32_primary.csv), two-sided α = 0.05, with SE = detectable difference / 2.8016 as in E44.
4. File modification times (UTC) of barriernet_0.2.0.zip and of the earliest files named for PRJEB25780 or Kim2018
   under the gastric project directory.

Outputs: results/E44b_checks.csv, results/E44b.log.
