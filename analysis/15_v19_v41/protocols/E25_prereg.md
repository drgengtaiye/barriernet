# Pre-registration: E25, arm-specific and mutually adjusted interactions in the six randomized comparisons

Written 2026-10-18, before the analyses Q2–Q4 below were run. The scores, data, harmonization and pooling
rule are those of E5, E10, E21 and E21_meta.R, unchanged: *barriernet* 0.2.1 gene lists, row-z mean
scoring with at least 3 genes per module, scores standardized to SD 1 within each comparison, the six
primary comparisons (IMmotion150, I-SPY2 pembrolizumab, I-SPY2 durvalumab–olaparib, CheckMate 025,
JAVELIN Renal 101, NeoTRIP) with their primary endpoint and model (PFS Cox; pCR logistic, I-SPY2
adjusted for hormone-receptor status and NeoTRIP for batch). Benefit-oriented scale (−log HR, +log OR).
Pooling: REML estimate, the wider of the Hartung–Knapp and DerSimonian–Laird 95% intervals and the larger P.
All results are reported whatever their direction.

**What has already been seen.** The per-trial interactions and within-arm effects (E5, E10, E21) and the
cancer-specific pooled interactions in E21_meta.csv are published results. Q1 below was previewed on
2026-10-18 for five of the six comparisons with a DerSimonian–Laird approximation (CAF_ECM, checkpoint
arms −0.34, P = 0.006; control arms −0.07, P = 0.22); it is therefore confirmatory only in the sense of
adding IMmotion150 and applying the manuscript rule. Q2–Q4 have not been computed in any form.

## Q1. Is the fibroblast association confined to the checkpoint arm?

Within-arm effect per SD of CAF_ECM, Immune and NetImmune, pooled separately over the six checkpoint arms
and the six control arms.

Expected: CAF_ECM adverse in the checkpoint arms (pooled P < 0.05) and not in the control arms (P ≥ 0.05);
NetImmune favorable in the checkpoint arms (P < 0.05) and not in the control arms (P ≥ 0.05).

## Q2. Do the immune and fibroblast signals modify benefit independently, in opposite directions?

In each comparison, one model with both standardized scores, the arm and both score × arm products (plus
the comparison's covariates). The two mutually adjusted interactions are pooled over the six comparisons
(primary), and over the three breast and three renal comparisons (secondary). The Immune–CAF_ECM
correlation in each comparison is reported.

Expected: adjusted Immune × arm > 0 and adjusted CAF_ECM × arm < 0. The statement "immune and fibroblast
content modify checkpoint benefit in opposite directions" is supported if both pooled adjusted
interactions have P < 0.05 over the six comparisons; supported for breast cancer only if that holds in the
breast set only; not supported otherwise.

## Q3. Do the pooled interactions exceed random gene sets of the same structure?

2,000 random gene sets drawn once from the bulk gene pool measured in all six comparisons, the same genes
in every comparison, with the structure of NetImmune (each immune and barrier module at the size present
in that comparison), of Immune (the immune part alone) and of CAF_ECM (11 genes). For each draw the six
interactions, the six checkpoint-arm effects and the mutually adjusted interactions of Q2 are computed and
pooled with the same rule as the observed scores. Empirical P = (#|null| ≥ |observed| + 1) / 2,001.
The fraction of random gene sets reaching a nominal pooled P < 0.05 is reported for each set of comparisons.

Expected: NetImmune pooled interaction (six comparisons) empirical P < 0.05; CAF_ECM pooled interaction in
the breast comparisons empirical P < 0.05; in the renal comparisons, neither below 0.05.

## Q4. Does the nominal interaction test become anti-conservative as the trial grows?

JAVELIN Renal 101 (726 patients), subsampled without replacement within arm to 100, 200, 300, 400, 500, 600
and 726 patients (ten subsamples per size), with 500 random gene sets of the NetImmune and CAF_ECM
structures per subsample. Reported: the fraction of random sets reaching nominal interaction P < 0.05 and
the SD of the null interaction z.

Expected: the fraction rises with the number of patients, from about 5% at 100 patients.

## Decision rule for the manuscript

The reframing ("survival, single-arm response and randomized interaction see different parts of the tumor
microenvironment; immune and fibroblast content predict checkpoint benefit in opposite directions") is
adopted only if Q1 holds for CAF_ECM and either Q2 or the breast-set CAF_ECM part of Q3 holds. Otherwise the
manuscript keeps its present framing and E25 is reported as a sensitivity analysis.

## Deviations

(none yet)
