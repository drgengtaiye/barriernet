# E45 protocol: a reference gene frozen in single-arm data, a coexpression-matched random reference, and the renal arms

Written 2026-11-05, after all other results were known and before E45_export.R, E45_select.py or E45_test.py is written or
run. Eighteenth protocol written after the other results were known. All analyses are post hoc.

## Why

An external review of the manuscript raised three points that the existing analyses do not answer:

1. CD8A was the only single gene analyzed before the signature calibration, and the panel genes were compared in the same
   data. The conclusion "no signature significantly exceeded CD8A" therefore rests on a reference chosen with knowledge of
   the randomized data.
2. Size-matched random gene sets are not matched on coexpression, so any coherent immune module can exceed them.
3. In the renal trials the control arms contain anti-angiogenic drugs, and part of an interaction could come from the
   control arm.

## What is already known

Every result of E23–E44. While checking the review (2026-11-05) the author's assistant also looked at: the 50 genes with the
lowest and the 40 with the highest six-comparison interaction among those at P < 0.05 (results/E33_pooled.csv), which include
KDR, FLT1, ANGPT2, ESM1, EGFL7, EMCN and PDGFRB on the control side; the within-arm slopes of Tcell-inflamed_Ayers, CTL_Jiang
and Teff_McDermott per comparison (results/E40_randomized_per_comparison.csv); the per-comparison interaction of CD8A; and the
pools of results/E32_pooled_all_sets.csv, including CD8A adjusted for angiogenesis (+0.265, P = 0.029).

Not yet looked at: any single-arm result for genes outside the signature library; any coexpression-matched reference; the
within-arm slopes of CD8A; the enrichment test of part C.

## A. Reference gene frozen in single-arm data

**Selection set.** The 20 single-arm checkpoint cohorts used to develop barriernet 0.3 (E37, E40), with the three samples
flagged by the twin check removed (1,195 patients, 368 responders). No randomized comparison and no gastroesophageal cohort
is in this set.

**Candidates.** The 2,082 genes measured in all six randomized comparisons (results/E33_pooled.csv, `all_six`), measured in at
least 16 of the 20 cohorts. Export: E45_export.R, loading and z-scoring exactly as E40_export.R; step 0: the 48 barriernet
genes reproduce results/E37_dev_genes.csv (< 1e-10).

**Rule.** In each cohort, logistic regression of response on the standardized gene (Firth, as E40). Pooled by the main-text
rule (REML estimate; wider of the Hartung-Knapp and DerSimonian-Laird intervals with the larger P).

- G1 (primary): the candidate with the highest pooled log-odds per SD.
- G2 (secondary): the same among the 27 panel genes (13 T-cell, 14 interferon-γ/checkpoint).
- If G1 and G2 are the same gene, there is one reference.

E45_select.py writes results/E45_single_arm_gene_pooled.csv and results/E45_reference_genes.tsv. The hash of the second file
is logged in E45_frozen.sha256 before E45_test.py is written to read it and run. The rank of CD8A is reported.

**Test in the six comparisons (E45_test.py).** Setup as E42 (E34_trials.py up to its fitting stage). For G1, G2 and CD8A: the
interaction per SD by the main-text rule, in the six comparisons (primary), the three renal and the three breast comparisons.
For each of the 52 signatures, NetImmune and barriernet 0.3: score minus reference by the paired bootstrap (300 resamples;
the seeds of E32/E42; fixed-effect inverse variance of the bootstrap variances). Counted: significantly above, significantly
below (P < 0.05), interval within ±0.10 per SD, median detectable difference (2.8 × SE). Also the reference minus CD8A. Step 0:
the 52 signatures and NetImmune minus CD8A reproduce results/E32_primary.csv (< 1e-8).

**Gastroesophageal cohorts (independent of the selection).** PRJEB25780, PRJEB60680 and PERFECT from results/E44_Z_*.csv.gz.
Log-odds of response per SD for G1, G2 and CD8A, as E44: per cohort; fixed effect (primary, as E44) and the main-text rule over
the three; in the two gastric cohorts, also adjusted for MSI-H or EBV-positive status. Reference minus CD8A by the paired
bootstrap (1,000 resamples stratified by response; seeds [4045, cohort]; fixed effect). Step 0: CD8A per cohort reproduces
results/E44_B_per_cohort.csv (< 1e-6).

## B. Random reference matched on coexpression

For each signature and comparison with the signature scored in the size-matched reference: let n be the number of genes present
and c the mean pairwise Pearson correlation among them across the comparison's samples.

Draws come from the same universe as the size-matched reference (E32 `UNIV`). One attempt: a seed gene uniformly at random; a
neighborhood size Q uniformly at random from {n, 2n, 5n, 10n, 20n, 50n, 100n, all genes} (capped at the universe); n − 1 genes
uniformly from the Q genes most correlated with the seed. The attempt is accepted if its mean pairwise correlation is within
±0.02 of c. If 1,000 sets are not accepted within 500,000 attempts, the tolerance is widened to ±0.05 for the remaining draws
and this is recorded. The signature's weights are assigned in random order, and each set is scored and modeled exactly as the
signature (`wscore`, `m_int`). Seeds [4545, comparison, crc32 of the signature name]. B = 1,000; draws pooled across comparisons
draw by draw by the main-text rule.

Empirical P in the declared direction, by the formula of the main analysis. Counted: signatures below 0.05, overall and by
family, against the 26 of the size-matched reference. Reported per signature: c per comparison, acceptance, median and 95th
percentile of the reference. This is a sensitivity analysis; the size-matched reference stays as specified.

## C. Renal arms and vascular genes

1. Within-arm slopes (`m_within`, benefit-oriented) per comparison and pooled over the three renal comparisons, for CD8A, G1,
   G2, the Immune axis, the Angio score (McDermott 2018: VEGFA, KDR, ESM1, PECAM1, ANGPTL4, CD34; z-mean; at least three genes)
   and TUMOR-VASCULATURE-UP_Lu.
2. Interaction of the Angio score and of TUMOR-VASCULATURE-UP_Lu in the three renal and the three breast comparisons.
3. Fisher exact test: the genes of TUMOR-VASCULATURE-UP_Lu measured among the 2,082, against the 102 genes favoring the control
   arm and, separately, the 101 favoring the checkpoint arm (results/E33_pooled.csv, P < 0.05).

## Expectations (written before the run; reported whatever the result)

| id | expectation |
|---|---|
| A1 | G1 has a positive pooled interaction in the six comparisons |
| A2 | at most 2 of the 52 signatures are significantly above G1 in the six comparisons |
| A3 | at most 2 of the 52 signatures are significantly above G2 in the six comparisons |
| A4 | G1 has a positive pooled log-odds in the three gastroesophageal cohorts (fixed effect) |
| B1 | fewer than 26 of the 52 signatures exceed the coexpression-matched reference |
| B2 | at least 5 of the 20 inflamed signatures still exceed it |
| C1 | CD8A has a negative pooled control-arm slope over the three renal comparisons |
| C2 | the Angio score has a negative pooled interaction in the three renal comparisons |
| C3 | TUMOR-VASCULATURE-UP_Lu genes are enriched among the control-favoring genes (P < 0.05) and not among the checkpoint-favoring genes |

If G1 or G2 is exceeded by more signatures than expected, the manuscript's statement about single genes changes accordingly.

Outputs: results/E45_*.csv, results/E45_select.log, results/E45_test.log, E45_outcome.md.
