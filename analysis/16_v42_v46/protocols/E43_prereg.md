# E43 protocol: predicting absolute benefit with the immune and stroma axes (six randomized comparisons; exploratory)

Written 2026-11-01, after all other results were known and before E43_benefit.py is written or run. The hash of this
file is logged in E43_prereg.sha256 first. Sixteenth protocol written after the other results were known.

## Why

All earlier tests asked whether a score modifies the relative treatment effect (the interaction). A single gene can do
that as well as any signature. The framework has a second axis that a single gene does not: a stromal module that is
prognostic but not predictive (CAF_ECM; response −0.06 per SD; TCGA hazard ratio 1.13). On the absolute scale, a
patient's benefit depends on both the baseline risk and the relative effect, so a model with an immune axis for the
relative effect and a stromal axis for the baseline risk might rank patients by absolute benefit better than CD8A
alone. This protocol tests that with the patients of the six randomized comparisons.

## Status of the data

The outcomes of all six comparisons have been analyzed before (E21–E37), and the interactions of CD8A, the Immune axis,
CAF_ECM and barriernet 0.3 are known. Nothing below is fitted to choose a model: the model forms are fixed here. All
results are exploratory and can at most motivate a test in new trials.

## Data

E34 setup (E32 data, six primary comparisons), executed unchanged as in E37_test.py. Endpoints and covariates as in the
primary model of each comparison: progression-free survival (IMmotion150, CheckMate 025, JAVELIN Renal 101; times in
months); pathological complete response (I-SPY2 pembrolizumab and durvalumab–olaparib with hormone receptor and
MammaPrint class; NeoTRIP with batch). Scores: CD8A; the Immune axis and the CAF_ECM module (published barriernet
scores); barriernet 0.3 (frozen weights, hash checked). Each score is standardized within the comparison over all of its
patients (no outcome used).

## Models (fixed; all include the arm and the comparison's covariates)

| model | terms besides arm and covariates |
|---|---|
| M_CD8A (reference) | CD8A, arm × CD8A |
| **M_2AX (primary)** | Immune, CAF_ECM, arm × Immune |
| M_BN03 | barriernet 0.3, arm × barriernet 0.3 |
| M_CD8A_CAF | CD8A, CAF_ECM, arm × CD8A |
| M_2AXI | Immune, CAF_ECM, arm × Immune, arm × CAF_ECM |
| M_RISK | Immune, CAF_ECM (no interaction) |

Logistic models by maximum likelihood, with Firth's correction under the E28 separation rule. Cox models by lifelines
(Efron ties), Breslow baseline. Predicted absolute benefit of a patient: P(pCR | checkpoint) − P(pCR | control), or
S(12 months | checkpoint) − S(12 months | control).

## Cross-validation

Within each comparison, 5 folds stratified by arm and outcome (pCR; or progression event), 20 repetitions, seeds
[2043, comparison index, repetition]. A model is fitted on four folds and predicts the fifth. A patient's predicted
benefit is the mean of its 20 out-of-fold predictions (a fold whose fit fails is left out for that repetition).

## C-for-benefit (primary metric; van Klaveren et al., J Clin Epidemiol 2018)

- Matching: within a comparison, each arm is sorted by predicted benefit; the patient of rank j (0-based) in the smaller
  arm (n_s) is paired with the patient of rank floor((j + 0.5) n_l / n_s) in the larger arm (n_l); the other patients of
  the larger arm are left out. Pair prediction: the mean of the two predicted benefits.
- Observed pair benefit: pCR, y(checkpoint) − y(control). Survival: +1 if the checkpoint patient's time is longer and
  the control patient's time is an event; −1 if the control patient's time is longer and the checkpoint patient's time
  is an event; 0 if both times are equal events; otherwise undetermined and left out.
- C-for-benefit: over all pairs of pairs with different observed benefit, the proportion in which the pair with the
  larger observed benefit has the larger prediction (ties in prediction count one half).
- Uncertainty: 300 within-arm patient bootstrap resamples per comparison (seeds [3043, comparison index]) of the
  matching and C with the predictions fixed, common to all models; paired differences between models from the same
  resamples. The bootstrap does not repeat the cross-validation.
- Pooling over the six: C − 0.5 per model by the manuscript rule (wider of Hartung–Knapp and DerSimonian–Laird);
  differences between models by fixed-effect inverse variance. The minimum detectable difference (2.80 × SE) is
  reported for the primary contrast.

## Secondary metric: benefit gradient on the absolute scale

Patients of a comparison are split at the median of their predicted benefit. In each half, the observed absolute
benefit is the difference in pCR proportion (checkpoint − control) or in Kaplan–Meier progression-free survival at 12
months. Gradient = upper half − lower half. Same bootstrap and pooling as above (gradient per model by the manuscript
rule; differences between models by fixed-effect inverse variance).

## Expectations

(B1) M_CD8A's pooled C-for-benefit is above 0.5 with its 95% CI excluding 0.5.
(B2, primary) The difference M_2AX − M_CD8A in C-for-benefit has a 95% CI that includes 0.
(B3) The difference M_CD8A_CAF − M_CD8A is smaller than 0.02 in absolute value (point estimate).
(B4) M_RISK's pooled C-for-benefit has a 95% CI that includes 0.5.

## Reading rules

- If M_2AX exceeds M_CD8A with a 95% CI excluding 0: reported as an exploratory finding on previously analyzed data, to
  be tested in new trials; not a claim that the framework predicts benefit better than CD8A.
- Otherwise: the text says that ranking patients by absolute benefit did not separate the two-axis model from CD8A.
- Either way, the result is reported in the Supplementary Information, with the minimum detectable difference.
