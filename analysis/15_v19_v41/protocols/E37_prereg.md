# E37 protocol: barriernet 0.3, developed without randomized data, frozen, then tested in the randomized comparisons

Written 2026-10-30, before any development or test analysis below. Post hoc: every randomized result of E25–E35 is known
to the author and the analyst (barriernet 0.2.1 NetImmune +0.32 per SD, P = 0.037, no score significantly above CD8A,
CAF_ECM pointing against benefit, best breast genes GZMB/IDO1/CXCL10). The hash of this file is logged in
E37_prereg.sha256 before E37_export.R runs. The frozen weights and their hash are logged in E37_frozen.sha256 before
E37_test.py runs; E37_test.py stops unless both hashes match.

## Why

The author asked for an improved algorithm. A version tuned on the six randomized comparisons could not be tested in
them. barriernet 0.3 is therefore developed only on single-arm checkpoint cohorts, frozen, and then put through the
calibration framework of the manuscript (treatment interaction, size-matched random genes, single genes) in the
randomized comparisons. The frozen version is registered for the controlled-access randomized trials before any of
them is requested (E37_registration.md).

## Development data (single-arm only)

The 21 single-arm checkpoint cohorts of the manuscript (11 from the companion compendium, 10 replication), response
R vs NR as recorded, minus ICB_Braun, which contains CheckMate 025 nivolumab-arm patients: 20 cohorts. Before any model
is fitted, development samples of kidney (ICB_Miao1, ICB_Shiuan) and breast (brca_iatlas_anders_2022) cohorts are
checked for expression twins among the samples of the randomized comparisons of the same cancer type (IMmotion150,
CheckMate 025, JAVELIN Renal 101; I-SPY2 ×2, NeoTRIP), by the E12 rule (Spearman correlation over shared genes; best
partner more than 0.03 clear of the second best); flagged samples are removed. This check reads no outcome.

Genes: the 48 barriernet 0.2.0 genes (fixed February 2026). Within each cohort, genes are z-scored across samples on the
log scale (as E5/E12); a gene not measured in a cohort is set to 0 there. Module scores as barriernet (mean z of the
genes present; 0 if fewer than three).

## Candidates

- C0: NetImmune 0.2.1 (fixed).
- C1: Immune axis 0.2.1 (fixed).
- C2: module weights — logistic regression of response on the six module scores with unpenalized cohort intercepts,
  ridge penalty (glmnet, alpha 0), lambda by leave-one-cohort-out deviance.
- C3: gene weights — the same on the 48 gene z-scores, elastic net (alpha 0.5).
- C4: gene weights — the same, ridge (alpha 0).

Each candidate is expressed as gene weights (C2: module coefficient divided by module size for each of its genes),
scaled to sum(|w|) = 1; zero weights dropped. Score = sum(w z) / sum(|w|) over the genes present (common engine), scored
when at least three genes and half of the non-zero-weight genes are present.

## Selection (development data only)

Leave-one-cohort-out: for each of the 20 cohorts, C2–C4 are fitted on the other 19 (lambda chosen by inner
leave-one-cohort-out over those 19, lambda.min) and scored in the held-out cohort; C0 and C1 are scored as they are.
In each held-out cohort the log-odds of response per SD of the score is estimated by logistic regression (as E12);
estimates are pooled by DerSimonian–Laird random effects (metafor). barriernet 0.3 is the candidate with the highest
pooled held-out log-odds; a candidate within 0.02 of a simpler one loses to it (order of simplicity C0, C1, C2, C4, C3).
The chosen candidate is refitted on all 20 cohorts (lambda by leave-one-cohort-out over the 20), and its weights are
frozen. If C0 is chosen, barriernet 0.3 equals NetImmune 0.2.1 and this is reported.

Development expectations: (D1) no learned candidate exceeds the Immune axis by more than 0.10 in pooled held-out
log-odds; (D2) the fitted CAF_ECM weight (C2) is at most zero.

## Test (randomized comparisons; opened only after the freeze)

Data and models as E33/E35 (E34 setup: six primary comparisons, I-SPY2 adjusted for hormone receptor and MammaPrint,
NeoTRIP for batch; NeoPembrOv as a sensitivity comparison; same within-arm bootstrap resamples: seeds [3032, comparison]
for the six, [3034, 7] for NeoPembrOv; 300 resamples).

| Q | Question | Expectation (written before the test) |
|---|---|---|
| Q1 (primary) | Pooled six-comparison interaction of barriernet 0.3 (manuscript pooling rule) | > 0 with P < 0.05 |
| Q2 | Paired difference from CD8A over the six | Not significant (P ≥ 0.05) |
| Q3 | Paired difference from NetImmune 0.2.1 over the six | Not significant (P ≥ 0.05) |
| Q4 | Against 1,000 size-matched random gene sets (its weights in random order; seed [2037, comparison, size]) | Above the 95th percentile |
| Q5 | Renal three, breast three, NeoPembrOv, seven comparisons | Descriptive |

## Registration

After the freeze, E37_registration.md states the frozen weights (with hash), and the primary analysis for each
controlled-access randomized dataset with expression from both arms (OAK and POPLAR, IMvigor010, IMmotion151,
CheckMate 649): score × arm interaction on the trial's primary survival endpoint, per SD, benefit-oriented, against CD8A
and NetImmune 0.2.1 by the paired bootstrap, and against size-matched random genes. Expectation for each: interaction
> 0; not significantly different from CD8A. Reported whatever the result.

## Reporting

Whatever D1–D2 and Q1–Q5 show is reported in the main text and in full in the Supplementary Information, with this
protocol, the development results of every candidate and the frozen weights.
