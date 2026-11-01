# E30 protocol: duplicate patients between GSE26253 and GSE62254, angiogenesis adjustment and endpoint strata in the randomized comparisons

Written 2026-10-21, after all other results were known and before any E30 model involving survival, response or
treatment arm is fitted. Prompted by a second external pre-submission review (points 1, 4 and 5) and by two
read-only checks of that review (session scratchpad `rev/agentA`, `rev/agentB`), which fitted no outcome model.
Like E28 and E29, E30 is *protocol-specified*: post hoc questions whose analyses are fixed here before they are run.

## Already known before this protocol (not re-tested)

- **Duplicates.** The 99.9th-percentile expression rule of E20 was never applied to the three locked cohorts
  against the benchmark cohorts. GSE26253 (Illumina WG-DASL, FFPE) and GSE62254 (Affymetrix U133 Plus 2, frozen)
  are both from Samsung Medical Center; cross-platform correlations top out at r = 0.36, below every
  within-cohort threshold (0.65–0.82), so that rule cannot detect shared patients between them. Among the 47
  mutual best cross-cohort expression matches between the two series (Pearson, all 11,520 shared genes,
  gene-centred within series, as in E20 `dup_pair()`), 11 agree on sex and on age within one year; under the same
  screen at most one such pair occurs between any other two of the ten gastric cohorts. Three of the 11 also agree
  to 0.1 month on recurrence/disease-free time and overall survival (GSM644726 = GSM1523990, GSM644801 =
  GSM1523747, GSM644802 = GSM1523975). The rule below was written after this candidate list was seen.
- **Pooling estimators.** Gastric pools use the DerSimonian–Laird (DL) τ² with Hartung–Knapp (HK) on that τ²;
  randomized pools use REML with REML-HK (E21/E25/E28); both report the wider of the HK and DL intervals and the
  larger P. For the first three locked gastric cohorts (NetImmune) P = 0.057 (DL-HK), 0.047 (REML-HK),
  3.3 × 10⁻⁵ (DL); for the first three randomized comparisons P = 0.057 (HK), 0.0020 (DL; τ² = 0).
- **CheckMate 025 with overall survival** instead of progression-free survival (E21 secondary endpoint): pooled
  six-comparison interactions NetImmune P = 0.059, ESTIMATE 0.107, TMEscore 0.095 / 0.129, CD8A 0.133, Immune
  0.25, IFN 0.18 (E21_meta.csv; workbook S16i).
- **I-SPY2.** The pembrolizumab (179 controls) and durvalumab–olaparib (29 controls) comparisons share no control
  patient (0 Research IDs in common; no expression match). Concurrency of controls cannot be established from the
  deposited data.
- **Angio score** (z-mean of VEGFA, KDR, ESM1, PECAM1, ANGPTL4, CD34; McDermott 2018) is scorable in all six
  comparisons (4/6 genes in IMmotion150, 5/6 in CheckMate 025, 6/6 elsewhere). Spearman ρ with NetImmune, Immune,
  CAF_ECM and Barrier is within ±0.14 in the three renal trials and 0.37–0.53 with CAF_ECM in the three breast
  comparisons. No outcome model with Angio has been fitted.

## Q1. Duplicates (primary for every estimate that uses GSE26253)

*Rule.* An array of GSE26253 is a duplicate of an array of GSE62254 if (i) the two are each other's best
cross-cohort expression match, (ii) sex agrees and age differs by at most one year, and (iii) their times do not
contradict each other: where GSE26253 records recurrence, GSE62254 disease-free time and overall survival are
each within 0.5 months of GSE26253 recurrence-free time and overall survival; where GSE26253 is censored,
GSE62254 disease-free time and overall survival are each at least the corresponding GSE26253 time minus 0.5 months.
Duplicates are removed from GSE26253 (the later-used cohort, as in E20); GSE62254 and every benchmark analysis
are unchanged. Expected: 8 arrays (the three above and GSM644593, GSM644733, GSM644810, GSM644871, GSM644940).
*Sensitivity.* Removing only the three arrays with matching event times.

*Re-run.* Every estimate in the manuscript, Supplementary Information or Supplementary Tables that uses GSE26253
is recomputed with the duplicates removed, by the published code path with only the sample filter changed:
per-cohort and pooled hazard ratios (locked three, locked seven, all ten; unadjusted and adjusted), the
structure-matched random-gene nulls that include GSE26253, the head-to-head ΔC, the E29 survival contrast and any
other pooled value that includes the cohort.
*Step 0 (gate).* With the filter switched off, each re-run reproduces its published values to |difference| < 1e-6.

*Expectation.* No reported pooled estimate changes by more than 0.01 on its reported scale; no conclusion changes.
*Decision rule.* The duplicate-free values become the primary values wherever they differ at reported precision
(text, tables; figures re-rendered when a plotted value moves by more than 0.01). The cohort counts become
GSE26253 424 patients and the corresponding totals. The previous values are kept in the Supplementary Tables as a
sensitivity analysis. The Methods state that the 99.9th-percentile rule cannot detect cross-platform duplicates
and describe this rule; the Supplementary Methods statement that expression matching found no duplicate between
GSE26253 and GSE62254 is withdrawn.

## Q2. Estimators (disclosure; no decision)

For every pooled estimate reported in the main text from k ≤ 7 cohorts or comparisons, a table gives the point
estimate and P under DL, REML, DL-HK and REML-HK, plus the manuscript rule. The Methods name the estimator used
for each data set and state that with three studies the HK interval uses t with two degrees of freedom
(critical value 4.30), so that the conservative rule can return P ≈ 0.06 for effects with DL P ≈ 10⁻³–10⁻⁵.

## Q3. Angiogenesis adjustment in the renal trials

*Model.* In IMmotion150, CheckMate 025 and JAVELIN Renal 101, the published primary model (same endpoint,
covariates and coding) with Angio and Angio × arm added: outcome ~ arm + score + Angio + score × arm + Angio × arm.
Angio is the z-mean of the available McDermott genes, standardized within trial. Scores: NetImmune (primary),
Immune, CAF_ECM, Barrier; ESTIMATE, TMEscore and CD8A (secondary). The breast comparisons are not adjusted: no arm
contains an anti-angiogenic drug and Angio tracks CAF_ECM there (ρ 0.37–0.53), so adjustment would remove part of
the signal being tested.
*Estimands.* (a) Adjusted score × arm interaction pooled over the three renal trials; (b) six-comparison pool
with the renal trials adjusted and the breast comparisons as published. Pooling as for the randomized comparisons
(REML; wider of REML-HK and DL intervals; larger P). Angio × arm itself is reported descriptively.
*Expectation.* Given |ρ| ≤ 0.14, each adjusted renal interaction is within 0.05 of the unadjusted one, and the
six-comparison NetImmune interaction keeps its sign with a 95% CI above zero.
*Decision rule.* If the six-comparison NetImmune interaction (b) keeps its sign with a CI excluding zero, the
manuscript states that the randomized result was not explained by angiogenesis. If the sign is kept but the CI
includes zero, it reports the attenuation and the main text describes the randomized evidence as not robust to
angiogenesis adjustment. If the sign reverses, the randomized result is reported as explained by angiogenesis.

## Q4. Endpoint strata

Each score's interaction pooled separately over the three hazard-ratio comparisons (renal) and the three
odds-ratio comparisons (breast), by the randomized pooling rule. *Decision rule (presentation).* Both strata are
reported side by side; the manuscript no longer compares their magnitudes, because cancer type and effect scale
are fully confounded. The six-comparison pool remains the pre-registered primary estimate.

## Q5. Weighting of the combined z

Stouffer's z weighted by √n (patients per comparison) is computed beside the unweighted version. *Presentation.*
The unweighted, weighted and calibrated (E28 parametric bootstrap) Stouffer P values move to the Supplementary
Tables; Table 4 keeps the random-effects estimates only.

## Q6. Disclosures without new computation

The full CheckMate 025 overall-survival sensitivity set (above) is reported in the Supplementary Tables and
summarized in the main text, including that every pooled interaction loses significance with it. The I-SPY2
statement (above) is added to the Methods.

## Outputs

Scripts `analysis/E30_dedup.R` (or `.py`) and `analysis/E30_trials.py`; results in `analysis/results/E30_*`;
outcome in `analysis/E30_outcome.md`, with each expectation marked held or not held.
