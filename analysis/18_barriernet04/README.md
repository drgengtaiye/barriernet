# barriernet 0.4.0 against the analysis code

`barriernet04_check.py` reruns the manuscript's analyses with the barriernet 0.4.0 package from the same inputs and
compares each result with the file the analysis code wrote; `barriernet04_check.csv` lists every comparison, its
tolerance and whether it held. No new analysis: nothing here relates new data to outcome.

- Benchmark: the six randomized comparisons as set up by `E32_trials.py` (up to its first stage), the 51 distinct
  library signatures, seeds 2032 (random reference) and 3032 (patient resamples, as E32 and E46). Compared with
  `E32_primary.csv` (pooled interactions, scoring, empirical P; Benjamini-Hochberg counts recomputed with scipy) and
  `E46_A_paired.csv` (paired differences from the two T-cell means).
- Single-arm cohorts: `E40_single_arm_genes.csv.gz` without the three twin-flagged samples, against
  `E40_single_arm_per_cohort.csv` and `E40_concordance.csv`.
- Published randomized trials: the table bundled with the package, against `E41_meta.csv`.
- Design: against `E39_required.csv`.

Result: 78 of 78 comparisons held. Pooled interactions agree within 1e-5 (the package's Cox solver against lifelines),
single-arm values within 1e-6 and published-trial pooling within 5e-5. Paired differences agree to 1e-4 for 50 of the
51 signatures; CRMA_Shukla differs because in JAVELIN Renal 101 its score takes 30 distinct values in 726 patients and a
few resamples give near-degenerate Cox fits; it is not significant either way and no count changes. Random-reference
draws are seeded differently from E32, so empirical P values differ by Monte Carlo error (at most 0.033) while every
classification at 0.05 is the same.

The script reads the lab's absolute paths, as the other analysis scripts here do (see `analysis/README_v48_map.md`).
