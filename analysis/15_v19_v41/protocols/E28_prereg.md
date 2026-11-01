# E28 pre-registration: comparator audit, interaction calibration and single-sample scoring

Written 2026-10-20, before any E28 model is fitted. Prompted by an external pre-submission review. All
analyses use the six primary randomized comparisons as loaded in E25_trials.py (step 0 must reproduce the
published interactions to 1e-6 before anything else runs) and the three gastric benchmark cohorts (TCGA-STAD,
GSE62254, GSE84437; revision/results/B1_*.rds). Pooling uses the manuscript rule (REML estimate; wider of the
Hartung–Knapp and DerSimonian–Laird interval; larger P) unless stated. Interactions are benefit-oriented
(−log HR for PFS, +log OR for pCR) per SD of the score.

Already known before this protocol (not re-tested): the official ssGSEA ESTIMATE immune score correlates
r = −0.50 with the mean z of PTPRC, CD2, CD3E and CD8A in GSE84437 (+0.88 in TCGA and GSE62254), whereas the
z-mean immune signature used in the manuscript correlates +0.87 (E9_estimate_check.R, rerun today). Renal-only,
breast-only and Stouffer pools are already in E21_meta.csv.

## Q1. ESTIMATE audit

1a. Diagnosis. In GSE84437, count missing values among the 141 immune and 141 stromal genes after
`filterCommonGenes`, rerun the official `estimateScore` after (i) removing genes with any missing value and
(ii) setting missing values to the gene's minimum, and report r with the four canonical immune genes.
Expected: the negative correlation is caused by missing values entering the within-sample ranking, and
disappears after (i) or (ii). Descriptive; no decision depends on it beyond the wording of the Methods.

1b. Official ESTIMATE in the six comparisons. Official contrast = z(ImmuneScore) − z(StromalScore) from
`estimate::estimateScore` (platform "affymetrix", as in B2b) on each comparison's log expression matrix.
Interaction per comparison and pooled.
Expected: correlation with the z-mean contrast ≥ 0.9 in every comparison, pooled interaction positive.
**Decision rule.** If the official pooled interaction is positive with P < 0.05, Table 4 keeps the z-mean
contrast and adds the official result in a footnote; otherwise both are shown and the text says that the
ESTIMATE result depends on its implementation.

## Q2. Direct test of Immune against CAF_ECM

CAF_ECM is reversed so that both scores are expected to favor the checkpoint arm (high Immune, low
CAF_ECM). Estimand Δ = θ(Immune) − θ(−CAF_ECM) = θ(Immune) + θ(CAF_ECM), benefit-oriented.
- Q2a, mutually adjusted: one model with both products; Var(Δ) from the model covariance.
- Q2b, marginal: the two published single-score models; Cov from 2,000 bootstrap resamples within each
  comparison (stratified by arm).
Pooled over six comparisons, and within renal and breast.
Expected: none (two-sided).
**Decision rule.** If the pooled Q2a Δ has P ≥ 0.05, the manuscript will not state or imply that Immune is a
stronger modifier of treatment benefit than CAF_ECM, and will report Δ. The asymmetry claimed for prognosis
and single-arm response is unaffected by this test.

## Q3. Calibration of the interaction test under no interaction

Parametric bootstrap from each comparison's fitted model **without** the score × arm term (score, arm and
covariates kept). Logistic: y* ~ Bernoulli(p̂). Cox: event times by inverting the Breslow cumulative baseline
hazard times exp(linear predictor); censoring times drawn from the reverse Kaplan–Meier estimate; follow-up
truncated at the largest observed time. B = 1,000 per comparison and score, for NetImmune, Immune, CAF_ECM,
ESTIMATE (z-mean) and CD8A. Firth is used when ML separates, as in the main analysis.
Reported: (i) per comparison, the share of null draws with nominal Wald P < 0.05; (ii) bootstrap-calibrated
P for the observed interaction; (iii) pooled: the six comparisons simulated jointly per draw and pooled with
the manuscript rule, giving the null share of pooled P < 0.05 and a calibrated pooled P.
Expected: the nominal test is close to calibrated (3–7% at α = 0.05) in every comparison, JAVELIN Renal 101
included.
**Decision rule.** If JAVELIN Renal 101's null rate is ≤ 7%, the 10–19% of random gene sets reaching P < 0.05
there (E25 Q4) is attributed to signal carried by random genes, not to an anti-conservative test: the
statement that nominal interaction tests were unreliable in the largest trial is withdrawn, and the
random-gene distribution is called a structure-matched reference distribution, not a null. If > 7%, the
statement stays and the bootstrap-calibrated P values become the primary P values in Table 4.

## Q4. Evidence tiers and multiplicity (no new models)

Tier 1, pre-registered replication: the three comparisons added in E21 (CheckMate 025, JAVELIN Renal 101,
NeoTRIP), pooled alone. Tier 2: all six (the first three motivated the extension). Tier 3: post hoc
(E23, E25, E28). For tier 2, Holm adjustment over the three axes (Barrier, Immune, NetImmune) and
Benjamini–Hochberg over the nine modules, axes and derived scores of Table 4.
**Decision rule.** The abstract and conclusions give the tier of every randomized-comparison P value they
cite; a claim resting only on tier 3 is labeled exploratory.

## Q5. Single-sample scoring

The manuscript z-scores each gene across the samples of a cohort. Three scorings that do not use the other
patients of the test cohort:
- (a) frozen reference: gene means and SDs from TCGA of the same tissue (KIRC for renal, BRCA for breast,
  Xena RSEM log2), never re-estimated in the trial;
- (b) reference-free: unit weights, the module mean of log2 expression;
- (c) within-sample rank: each gene's rank among all measured genes of the sample, divided by their number;
  the module score is the mean rank of its genes.
Axes are built as in the manuscript. Reported: Pearson r with the cohort-relative score; interaction per
comparison and pooled for NetImmune, Immune and CAF_ECM; in the gastric benchmark cohorts, HR per SD for
NetImmune and CAF_ECM with (b) and (c), and (a) with TCGA-STAD as the reference for the two array cohorts
(cross-platform, flagged).
Expected: r ≥ 0.9 and pooled NetImmune interaction within ±0.05 of +0.29 for (a) and (b); (c) may differ more.
**Decision rule.** If (a) and (b) meet the expectation, Methods state that the results do not depend on
cohort-relative standardization; otherwise the manuscript states that barriernet is a cohort-relative score
that needs a reference cohort for a single patient. The per-SD scaling of the reported effect is a reporting
unit in both cases.

## Q6. Increment over ESTIMATE in the six comparisons

(a) Paired difference θ(NetImmune) − θ(ESTIMATE z-mean) per comparison, SE from 2,000 arm-stratified
bootstrap resamples, pooled. (b) Both products in one model: adjusted interactions.
Expected: difference near zero. **Decision rule.** The manuscript makes no claim that barriernet predicts
treatment benefit better than ESTIMATE, and reports (a).

## Q7. Exploratory: PD-L1 in JAVELIN Renal 101

NetImmune × arm adjusted for PD-L1 status (PDL1FL) and PD-L1 × arm. Exploratory; no decision.
