# E44 protocol: what the calibration framework can say about gastric and gastroesophageal adenocarcinoma

Written 2026-11-02, after all other results were known and before E44_export.R or E44_gastric.py is written or run. The
hash of this file is logged in E44_prereg.sha256 first. Seventeenth protocol written after the other results were known.

## Why

The author asked for a gastric section without changing the main line of the manuscript. No randomized checkpoint trial
in gastric or gastroesophageal cancer has deposited expression from both arms, so the interaction test of the framework
cannot be run on gastric patient-level data. This protocol collects what the data at hand can say:

- A. the published signature-stratified results of randomized gastric trials;
- B. the published signatures, single genes and barriernet 0.3 in single-arm gastroesophageal checkpoint cohorts;
- C. prognosis without checkpoint inhibitors in gastric cancer;
- D. a score developed in gastric cancer (TMEscore) in the six randomized comparisons.

## What is known before writing

- A: the E41 values of the two gastric trials: CheckMate 649 (nivolumab + ipilimumab against chemotherapy; 10-gene
  inflammatory signature, hazard ratio 0.63 in the top tertile, 0.91 in the bottom) and KEYNOTE-061 (pembrolizumab
  against paclitaxel; T-cell-inflamed GEP, 0.88 non-low, 0.83 low). While planning this protocol the analyst combined
  the two by hand (about +0.10 per SD, P about 0.3). Part A is therefore descriptive, with no expectation.
- B: the barriernet 0.2.1 modules have been analyzed in all three cohorts. Immune axis log-odds of response per SD:
  PRJEB25780 +1.32 (SE 0.50; E12), PERFECT +0.40 (0.38; E16), PRJEB60680 +0.40 (0.34; E22). The analyst combined these
  by fixed effect by hand while planning (about +0.59, P about 0.01). NetImmune was null in PRJEB60680 (−0.06), and its
  three MSI-H or EBV-positive patients all responded. The companion manuscript analyzes the 52 signatures in PRJEB25780
  (its compendium's gastric cohort); the author has seen those results, the analyst does not recall them. The 52
  signatures, CD8A, the single-gene panel and barriernet 0.3 have not been analyzed in PRJEB60680 or PERFECT.
  None of the three cohorts was used to develop barriernet 0.3 (E37 used 20 cohorts, none gastric or oesophageal), but
  PRJEB25780 was among the cohorts explored when the 48 barriernet genes were fixed (barriernet 0.2.0), so barriernet 0.3
  is not independent of it at the level of its genes.
- C: the E20 pooled hazard ratios over the 10 gastric cohorts (primary duplicate rule of E30/E31): CD8A 0.90 per SD
  (0.84–0.96), Immune axis 0.84, ESTIMATE 0.73. Part C reports them; there is no new analysis.
- D: TMEscore's six-comparison interaction (E32: +0.31, 0.04–0.57), its random-gene percentile and its paired
  difference from CD8A (+0.06, P = 0.44) are known. Part D reports them.

## A. Published randomized gastric trials (descriptive)

Input: results/E41_inputs.csv, the rows with tumor_group gastric (CheckMate 649, KEYNOTE-061), with the values and
scales of E41 unchanged: y, se (log ratio of hazard ratios, low over high, positive = more benefit when the signature is
high) and y_sd, se_sd (per SD of the score).

- Per trial: ratio of hazard ratios and per-SD value with 95% CI; the per-SD interaction detectable with 80% power at
  two-sided 0.05, 2.80 × se_sd.
- The two trials combined by fixed effect (inverse variance) and by DerSimonian–Laird; Cochran's Q between them.
- The two against the other 16 E41 trials: each group combined by DerSimonian–Laird on the per-SD scale; difference with
  a z-test.
- The detectable interaction is set beside the six-comparison CD8A interaction (+0.23 per SD) and the 0.10-per-SD margin.

## B. Single-arm gastroesophageal checkpoint cohorts

Cohorts (pre-treatment samples; response labels and processing as used before):

1. PRJEB25780: metastatic gastric cancer, pembrolizumab (Kim 2018); 44 patients, 12 responders. Expression from
   results/E12_cohorts.rds (deposited log2 scale); genes z-scored over the cohort (E12).
2. PRJEB60680: gastric cancer, one cycle of 5-FU and platinum, then pembrolizumab added (NCT04249739); the 43 QC-pass
   patients, 29 with CR or PR (E22 primary). log2(TPM + 1) over nuclear protein-coding and IG/TR genes; genes z-scored
   within library batch, a gene constant within a batch set to 0 there (E22 primary).
3. PERFECT (GSE165252): oesophageal adenocarcinoma, neoadjuvant atezolizumab with chemoradiotherapy; 32 baseline
   biopsies, 12 responders (E16). log2 of the deposited vst matrix, Ensembl IDs to symbols as E16; genes z-scored over
   the cohort.

Primary set: all three (gastroesophageal adenocarcinoma, the population of CheckMate 649). Sensitivity: the two gastric
cohorts. GSE183924 (adjuvant durvalumab) is not used: it has relapse-free survival only and no comparator.

Export (E44_export.R): for each cohort, the z-scored matrix of all genes, the response, and MSI-H or EBV-positive status
where recorded (PRJEB25780: data_download/validation/kim2018_cps/kim2018_cps_clinical.csv, msi_h or ebv_pos; PRJEB60680:
Wright 2026 Table S1, MSI-H or EBV positive). Step 0: the Immune axis log-odds per SD computed from the exported
matrices (unpenalized logistic regression, as E12/E16/E22) reproduces the three values above within 1e-6.

Scores (E32 engine): the 52 library signatures (weighted mean z of the genes present, at least three genes and half of
the list), ESTIMATE, TMEscore, the barriernet 0.2.1 Immune axis and NetImmune, barriernet 0.3 (frozen weights, hash
checked), CD8A and the 27 genes of the T-cell and interferon/checkpoint panels. Raw-score orientation, as E40.

Per cohort: logistic regression of response on the standardized score (lfit of E28: maximum likelihood, Firth if the
data separate); log-odds per SD.

Pooling. Primary: fixed effect, inverse variance, over the three cohorts (three cohorts cannot estimate the
between-cohort variance). Sensitivity: the manuscript's conservative rule (pool_rule). A score enters a pooled set only
if it is scored in every cohort of that set.

Paired difference from CD8A. Within each cohort, 1,000 bootstrap resamples stratified by response (numpy seed [4044, i],
i = 0, 1, 2 for PRJEB25780, PRJEB60680, PERFECT). In each resample both models are refitted (lfit) and the difference
of log-odds is taken; the per-cohort SE is the SD of the differences over the resamples in which both fits succeeded
(failures are counted and reported). The observed differences with these SEs are pooled by fixed effect.

Adjusted sensitivity (the two gastric cohorts): response on the score and the MSI-H/EBV-positive indicator (lfit, Firth
if the data separate), pooled by fixed effect.

Agreement with the randomized interaction: Spearman correlation, over the 52 signatures pooled in the primary set,
between their pooled log-odds here and their six-comparison pooled interaction (results/E32_primary.csv, set "six
(primary)", raw orientation). The signatures are correlated, so its P is descriptive.

Expectations (written before any part of B is computed):

| Q | Expectation |
|---|---|
| G1 (primary) | CD8A pooled log-odds > 0 with P < 0.05 (fixed effect, three cohorts) |
| G2 | No library signature has a pooled paired difference from CD8A > 0 with P < 0.05 |
| G3 | Spearman correlation with the randomized interaction > 0 with P < 0.05 |
| G4 | barriernet 0.3 pooled log-odds > 0 with P < 0.05, and its paired difference from CD8A has P ≥ 0.05 |

Descriptive: the number of signatures with P < 0.05 here and, of these, the number with a same-direction randomized
interaction at P < 0.05 (as E40); the rank of CD8A among the 52 signatures and the 27 panel genes; the number of
signatures significantly below CD8A; the gastric-only set; the MSI/EBV-adjusted models; the conservative-rule pooling;
per-cohort values.

## C. Prognosis without checkpoint inhibitors (report only)

From results/E30_dedup_runs/primary/E20_pooled.csv, set "all 10 gastric", unadjusted, duplicate rule "pre-registered
exclusions (primary)": CD8A, the Immune axis, T_CELL, IFN, ESTIMATE.

## D. A score developed in gastric cancer in the randomized comparisons (report only)

From results/E32_primary.csv, set "six (primary)", TMEscore_z: pooled interaction, random-gene percentile (emp_p_declared)
and paired difference from CD8A.

## Reporting

Whatever A–D show is reported. If a gastric section is added to the manuscript, it reports the expectations above with
their results, and states that the gastric cohorts are single-arm and cannot separate a predictive from a prognostic
marker. With this protocol the later protocols number 17, with 67 expectations.

Outputs: results/E44_export.log, E44_A_gastric_trials.csv, E44_B_per_cohort.csv, E44_B_pooled.csv, E44_B_paired.csv,
E44_B_adjusted.csv, E44_C_prognosis.csv, E44_D_tmescore.csv, E44_expectations.csv, E44.log.
