# E29 pre-registration: direct test of Immune against CAF_ECM outside the randomized comparisons

Written 2026-10-21, before any E29 model is fitted, after E28 (E28_outcome.md) found no difference between
Immune and CAF_ECM in the six randomized comparisons (Δ = +0.07, −0.12 to 0.26, P = 0.47). Prompted by the same
external pre-submission review (point 3D: "one significant, the other not" does not show that two components
differ). The manuscript's dissociation claim ("fibroblast content tracks survival but not response; immune
programs track checkpoint response") now rests on the single-arm response cohorts and the gastric survival
cohorts, where it has so far been supported only by marginal estimates. E29 tests it directly in both.

Already known before this protocol (published marginal estimates, not re-tested): in the 21 single-arm
checkpoint cohorts Immune +0.43 and CAF_ECM about −0.06 log-OR of response per SD; in the ten gastric survival
cohorts Immune HR 0.84 (0.78–0.91) and CAF_ECM HR 1.30 (1.17–1.45) per SD (E20_pooled.csv, "all 10 gastric").
The correlation between the two scores within cohorts has not been examined and is not known.

## Data and scoring (no new data; scores exactly as published)

- **Response (single-arm):** the 21 cohorts of the published pooled set (E12 "all" = the 11 E5 cohorts loaded
  and scored as in E5_icb_split.R + the 10 "new" cohorts from E12_cohorts.rds scored as in
  E12_icb_expanded.R). Response = R/NR as published. Each score standardized to SD 1 within cohort.
- **Survival (gastric):** the ten cohorts of the published "all 10 gastric" set, unadjusted, published endpoint
  per cohort (OS; RFS only where it is the only outcome) and published duplicate exclusions: TCGA-STAD,
  GSE84437, GSE62254 (revision/results/B2_scores.rds), GSE15459, GSE34942, GSE26253 (E15_scores.rds), GSE13861,
  GSE26901, GSE26899, GSE29272 (E20_scores.rds). Cox model, score standardized to SD 1 within cohort.
- **Step 0 (gate).** Before anything new is computed, the script must reproduce the published marginal
  per-cohort estimates of Immune and CAF_ECM (E12_icb_per_cohort.csv; E20_estimates_used_for_pooling.csv,
  unadjusted, primary duplicate rule) to |difference| < 1e-6 in every cohort. If any cohort fails, E29 stops
  and the failure is reported; no E29 estimate is used until it is resolved.

## Estimand

Both components are oriented so that a positive value means "behaves as hypothesized": for response,
Immune +log OR and CAF_ECM −log OR; for survival, Immune −log HR and CAF_ECM +log HR. The dominance contrast is
D = strength(Immune) − strength(CAF_ECM):

- response: D_r = b(Immune) + b(CAF_ECM) (log OR per SD) — identical in form to E28's Δ;
- survival: D_s = −[β(Immune) + β(CAF_ECM)] (log HR per SD).

The dissociation hypothesis predicts D_r > 0 (immune dominates response) and D_s < 0 (fibroblast content
dominates survival).

- **(a) Mutually adjusted (primary):** one model per cohort with both scores; Var(D) from the model covariance.
- **(b) Marginal:** the two published single-score models; Cov from 2,000 bootstrap resamples within each
  cohort (stratified by response for (a)/(b) response; by event status for survival), seed 2029.

Pooling: the rule used for these two data sets in the manuscript — DerSimonian–Laird estimate; the wider of the
Hartung–Knapp and DerSimonian–Laird intervals; the larger P. Two-sided α = 0.05.

## Questions, expectations and decision rules

**Q1 (primary). Response, 21 cohorts, D_r (a).** Expected: D_r > 0 (from the marginal estimates, about +0.35).
*Decision rule.* If pooled D_r (a) > 0 with P < 0.05, the manuscript may state that in single-arm checkpoint
cohorts immune programs tracked response more strongly than fibroblast content, citing the direct test. If
P ≥ 0.05, it may report only the two marginal estimates and must not state or imply a difference between them
for response.

**Q2 (primary). Survival, ten gastric cohorts, D_s (a).** Expected: D_s < 0 (from the marginal estimates,
about −0.09; the direct test may be underpowered).
*Decision rule.* If pooled D_s (a) < 0 with P < 0.05, the manuscript may call fibroblast content the stronger
prognostic component. If P ≥ 0.05, "the strongest prognostic component" must be worded as the largest point
estimate, not as a demonstrated difference from the immune modules.

**Q3 (secondary, descriptive). Crossover.** D_r − D_s with variance Var(D_r) + Var(D_s) (independent data
sets). The two contrasts are on different scales (log OR and log HR per SD), so Q3 is reported as a direction
summary only and carries no decision. The word "dissociation" in the title and abstract is retained only if Q1
holds and D_s is negative in point estimate; otherwise the authors are told that the title over-states the
evidence (no automatic rewording).

**Q4 (descriptive). Collinearity.** Spearman correlation of Immune and CAF_ECM within each cohort (median and
range per data set). If the median |r| > 0.7 in a data set, its (a) estimates are flagged as unstable and (b)
is reported beside them.

## Sensitivity analyses (no decisions)

S1 marginal (b) for Q1 and Q2. S2 REML pooling. S3 response: E5 11 and new 10 separately; whole-transcriptome
cohorts only. S4 response: the two gastric checkpoint cohorts (PRJEB25780 as in E12; PRJEB60680 primary
batch-aware scores from E22_scores_PRJEB60680.csv), each reproduced against its published marginal estimates
first; reported beside the 21, not pooled with them. PERFECT is not included (patient-level scores were not
saved by E16a). S5 survival: benchmark 3 and locked 7 separately; adjusted for age and stage where recorded
(as in E20). S6 response and survival: each immune module (APM, IFN, TLS_B, T_CELL) in place of the Immune
axis, (a) only.

## Outputs

Scripts E29_dissociation.R; results E29_* in results/; outcome table E29_outcome.md written after the run
against this file (sha256 in E29_prereg.sha256, written before the run). The manuscript is not edited by E29.
