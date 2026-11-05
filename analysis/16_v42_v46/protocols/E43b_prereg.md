# E43b addendum (post hoc): bias of the cross-validated C-for-benefit, and a null-centered estimate

Written 2026-11-01 (local time), after the E43 results and after an independent audit of E43_benefit.py, and before
E43b_null.py is written or run. The hash of this file is logged in E43b_prereg.sha256 first. This addendum has no
expectations and is not counted in the tally of later protocols; the E43 results are reported as pre-specified.

## Why

The audit found that the pre-specified cross-validated predictions bias C-for-benefit and the benefit gradient
downward. Each patient's prediction comes from fits that leave that patient out, so its arm effect is shifted against
the patient's own outcome, and rank matching pairs exactly those patients. In the audit's null simulations (no effect
modification), the cross-validated C was 0.31 for a risk-only Cox model and 0.43 for a model with an interaction,
against about 0.50 in sample. The bias grows as the spread of predictions shrinks, so contrasts favor the model whose
predictions vary more. E43's M_RISK values in the survival comparisons (0.20 to 0.38) fit this pattern.

## Step 1: the bias in the real data

For each comparison, 20 joint permutations of the patients' four standardized scores (CD8A, Immune, CAF_ECM,
barriernet 0.3) within arm (seeds [4043, comparison index, permutation]); the arm, outcome and covariates stay in place.
Each permuted data set goes through the unchanged E43 pipeline (same folds and seeds, 20 repetitions), and C and the
gradient are computed on all patients. Reported: the mean null C and gradient per model and comparison.

## Step 2: null-centered apparent estimate

- Apparent predictions: each model fitted on all patients of a comparison, predicting all of them.
- Null: 200 joint within-arm permutations of the scores (seeds [5043, comparison index, permutation]), each refitted on
  all patients; C and the gradient on all patients.
- Null-centered C = apparent C − mean null C; same for the gradient. This removes the optimism that an in-sample fit
  has when there is no signal.
- Uncertainty: the E43 evaluation bootstrap (300 within-arm resamples, seeds [3043, comparison index]) on the
  apparent predictions, which are held fixed; the centering constants add no variance.
- Pooling: per model by the manuscript rule; differences between models by fixed-effect inverse variance; for the
  primary contrast M_2AX − M_CD8A, the difference between the survival three and the pCR three (Cochran Q, 1 df).
- Reading: in the survival comparisons, the pair outcome of C-for-benefit (who progressed first) is on the relative
  scale, so there the gradient (Kaplan–Meier difference at 12 months) is the absolute-scale measure.

Everything in this addendum is exploratory and post hoc.
