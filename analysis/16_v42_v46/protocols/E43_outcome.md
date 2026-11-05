# E43 outcome: predicting absolute benefit with the immune and stroma axes (exploratory)

Protocol E43_prereg.md (sha256 e16d29823c6d…, logged 2026-11-02T01:24:04Z, before E43_benefit.py was written).
Sixteenth protocol written after the other results were known. Run: E43_benefit.py → results/E43.log,
E43_per_comparison.csv, E43_pooled.csv, E43_contrasts.csv, E43_coefficients.csv, E43_expectations.csv.
Post hoc addendum E43b_prereg.md (sha256 ec3fd91acb32…, logged 2026-11-02T01:43:49Z, after the E43 results and the
audit, before E43b_null.py was written): results/E43b.log, E43b_cv_null.csv, E43b_apparent.csv, E43b_pooled.csv,
E43b_contrasts.csv, E43b_heterogeneity.csv.

Step 0: the full-data arm × CD8A coefficient of M_CD8A equaled the manuscript's interaction in all six comparisons. No
fold fit failed. The protocol's background value for CAF_ECM in TCGA (hazard ratio 1.13) is the 24-type adjusted
overall-survival estimate of E11 (1.13, 1.06 to 1.21); the manuscript quotes the five checkpoint-type estimate (1.07).

## Pre-specified analysis (cross-validated): expectations 3 of 4 held, but the metric is biased

| | expectation | result | held |
|---|---|---|---|
| B1 | M_CD8A pooled C-for-benefit > 0.5, CI excluding 0.5 | 0.572 (0.476 to 0.668) | no |
| B2 (primary) | M_2AX − M_CD8A: CI includes 0 | −0.012 (−0.048 to +0.025); MDD 0.052 | yes |
| B3 | abs(M_CD8A_CAF − M_CD8A) < 0.02 | −0.011 | yes |
| B4 | M_RISK pooled C-for-benefit CI includes 0.5 | 0.441 (0.232 to 0.650) | yes |

An independent audit of the code found that the pre-specified cross-validated predictions bias C-for-benefit and the
gradient. Each patient's prediction comes from fits that leave the patient out, so the inherited arm effect is
shifted against the patient's own outcome, and rank matching pairs exactly those patients. The audit's null
simulations gave a cross-validated C of 0.31 for a risk-only Cox model (0.50 in sample). E43b step 1 confirmed this on
the real data: with the four scores permuted within arm, the cross-validated C averaged 0.26 to 0.39 for M_RISK in the
survival comparisons and 0.40 to 0.61 for the other models, depending on the comparison. The bias depends on the
spread of the predictions, so it also biases contrasts between models. The pooled survival-three contrasts of E43
(M_2AX − M_CD8A −0.089, M_BN03 − M_CD8A −0.108) are largely this artifact. The expectations are reported as computed;
their verdicts were the same in the null-centered analysis below.

The audit also noted that in the survival comparisons the pair outcome of C-for-benefit (who progressed first) is on
the relative scale, so it cannot reward absolute-benefit differences that come from baseline risk; there the
gradient (Kaplan–Meier difference at 12 months) is the absolute-scale measure.

## Post hoc null-centered analysis (E43b)

Apparent (full-data) predictions; each value minus the mean of 200 within-arm permutations of the scores, which keeps
the contribution of the arm and the covariates (in I-SPY2, hormone receptor and MammaPrint carry absolute-benefit
information of their own: null C 0.65 to 0.68). The centered value is what the scores add.

| six comparisons | null-centered C-for-benefit | null-centered gradient |
|---|---|---|
| M_CD8A | −0.008 (−0.050 to +0.034) | +0.052 (−0.099 to +0.202) |
| M_2AX | +0.036 (−0.036 to +0.108) | +0.041 (−0.085 to +0.167) |
| M_BN03 | +0.024 (−0.019 to +0.067) | +0.081 (−0.068 to +0.231) |
| M_CD8A_CAF | +0.025 (−0.053 to +0.103) | +0.059 (−0.098 to +0.216) |
| M_2AXI | +0.054 (+0.012 to +0.096) | +0.044 (−0.081 to +0.168) |
| M_RISK | −0.028 (−0.106 to +0.050) | −0.028 (−0.161 to +0.106) |

| contrast (six) | C-for-benefit | gradient |
|---|---|---|
| M_2AX − M_CD8A (primary) | +0.033 (−0.004 to +0.069), P = 0.078; MDD 0.052 | −0.023 (−0.115 to +0.068) |
| M_BN03 − M_CD8A | +0.033 (−0.007 to +0.073) | −0.004 (−0.098 to +0.090) |
| M_CD8A_CAF − M_CD8A | +0.015 (−0.004 to +0.035) | +0.001 (−0.010 to +0.013) |
| M_2AXI − M_CD8A | +0.059 (+0.014 to +0.105), P = 0.011 | −0.037 (−0.137 to +0.062) |
| M_RISK − M_CD8A | +0.015 (−0.037 to +0.066) | −0.038 (−0.158 to +0.083) |

- After centering, the survival three and pCR three no longer differ for the primary contrast (C-for-benefit
  Q = 0.20, P = 0.66; gradient Q = 2.69, P = 0.10).
- pCR three (post hoc subset): M_BN03 − M_CD8A +0.058 (+0.003 to +0.112) in C and +0.190 (+0.014 to +0.367) in the
  gradient, the breast pattern of E37 seen again in the same patients.

## Audit of E43b (independent agent, synthetic simulations only)

- Code faithful to the addendum; all quoted numbers reproduce from the CSVs (pools to 6 × 10⁻¹⁰).
- Null-centering is conservative for a single model: in-sample optimism shrinks as real signal grows, so a model with
  signal is centered about 0.02–0.03 too low. M_CD8A's −0.008 is therefore not evidence of no signal.
- For contrasts of two one-interaction models (M_2AX − M_CD8A) the centering shift is about 0.002, and the value can be
  biased toward the weaker model by up to about 0.02–0.03. For a model with one extra interaction (M_2AXI) the
  optimism gap is about +0.02 with or without signal, which the centering removes.
- The evaluation bootstrap overstated the spread across data sets (conservative standard errors).
- Centering the cross-validated C on its own permutation null is not valid (the artifact is not additive) and reverses
  the sign of the contrasts (M_2AX − M_CD8A −0.053, M_2AXI − M_CD8A −0.059, per-comparison means). M_2AXI's
  +0.059 rests on the choice of correction and is exploratory only.
- M_RISK's apparent C below its null in IMmotion150 (0.37) and JAVELIN Renal 101 (0.39) is unlikely to be an artifact:
  a risk-only ranking gives about 0.5 under proportional hazards, so it points to relative-scale heterogeneity running
  against baseline risk, or to non-proportional hazards.
- Within-arm permutation also breaks the correlation of the scores with the I-SPY2 covariates; a permutation within
  hormone receptor × MammaPrint strata was not run.

## Reading

Ranking patients by absolute benefit did not separate the two-axis model from CD8A. In the pre-specified analysis the
primary contrast was −0.012 (−0.048 to +0.025); that analysis is biased. In the null-centered analysis it was
+0.033 (−0.004 to +0.069) in C-for-benefit and −0.023 (−0.115 to +0.068) in the absolute-scale gradient, so the two
measures do not even agree in direction. Adding the stromal axis to CD8A changed nothing (+0.015; gradient +0.001).
The model with both interactions (M_2AXI) was the only one above chance and above CD8A in C-for-benefit (+0.059,
P = 0.011), but not in the gradient; it is one of about 30 post hoc contrasts, and its sign depends on the correction
(−0.030 in the pre-specified analysis). CD8A itself added nothing
measurable beyond the arm and the covariates (−0.008, a conservative value), although its interaction is +0.23 per SD: with six comparisons
of this size, C-for-benefit has little power. The premise that a prognostic stromal axis would improve absolute-benefit
prediction found no support here.
