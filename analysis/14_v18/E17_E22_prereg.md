# Pre-registration: validation analyses E17–E22

Written 2026-10-04, before any of the data below were downloaded or related to outcome. The scores are fixed:

- the *barriernet* 0.2.1 gene lists;
- row-z mean scoring with at least 3 genes per module;
- the axes Barrier, Immune and NetImmune;
- scores standardised to SD 1 within each cohort.

No gene, weight or threshold may be changed.

Rules carried over from the manuscript:

- **Pooling.** Random-effects meta-analysis, reporting the wider of the Hartung–Knapp and DerSimonian–Laird 95% intervals and the larger *P*.
- **Interactions.** Benefit-oriented: −log HR for time-to-event endpoints and +log OR for binary endpoints, so that a positive value means more benefit from the checkpoint-containing arm.
- **Expression scale.** Check the deposited scale of every matrix before scoring; apply log₂ once only if it is linear.
- **Reporting.** All results are reported whatever their direction. Any deviation from this plan is recorded at the end of this file with its date and reason.

## E17. Independent single-cell attribution

**Data.** Gastric cancer scRNA-seq atlases other than GSE183904, using tumour samples only. Candidates are GSE163558 and other atlases found by the data search.

**Method.** Cell types are annotated from canonical markers only, without the module genes; the annotation is fixed before module scores are computed.

**Expected:**

- CAF_ECM is highest in fibroblasts.
- TAM_SUP is highest in macrophages and monocytes.
- TLS_B is highest in B and plasma cells.
- T_CELL is highest in T and NK cells.
- Weighted by abundance:
  - fibroblasts contribute ≥ 50% of CAF_ECM signal;
  - T and NK cells contribute ≥ 50% of T_CELL signal;
  - macrophages contribute < 50% of TAM_SUP signal, replicating the 28% in the discovery atlas.

## E18. PD-L1 CPS in PRJEB25780

**Data.** Per-patient PD-L1 combined positive score, if it can be obtained from the original publication.

**Expected:**

- NetImmune correlates positively with CPS.
- In a Firth logistic model with MSI/EBV status and CPS, the NetImmune coefficient is positive. With 44 patients this is not expected to be significant.
- Reported: the AUC of CPS, the AUC of NetImmune, and the penalised likelihood-ratio *P* for NetImmune beyond MSI/EBV + CPS.

## E19. Robustness to gene loss and platform

**Data.** The three gastric benchmark cohorts (overall survival), PRJEB25780 (response) and the 21 checkpoint cohorts (response).

**Perturbations:**

- 1,000 random drops of 10%, 20% and 30% of the 48 genes, keeping at least 3 genes per module;
- leave-one-gene-out;
- the gene subset shared with the targeted-panel cohorts.

**Expected:**

- The median Spearman correlation with the full score is ≥ 0.9 at 30% dropout.
- The direction of the pooled NetImmune HR (< 1), of the pooled Immune log-OR (> 0) and of the pooled CAF_ECM log-OR (≤ 0.1) is unchanged in ≥ 95% of draws.

## E20. Extended locked gastric prognostic replication

**Data.** Gastric cohorts with recoverable survival that were not used so far. Candidates are GSE26901, GSE13861, GSE29272, GSE57303 and GSE51105, with outcomes taken from GEO or the original papers. They are de-duplicated against all cohorts already used.

**Status.** This extension is decided after the three-cohort locked replication gave *P* = 0.057, and is reported as such.

**Expected:**

- NetImmune has HR < 1 and CAF_ECM has HR > 1 in each cohort; overall survival is primary, and recurrence-free survival is used where it is the only outcome.
- The primary pooled estimate is NetImmune across all locked cohorts (the previous 3 plus the new ones): HR < 1.
- Also reported:
  - all gastric cohorts pooled;
  - the structure-matched random-gene null, applied as in E15.

## E21. Further randomised comparisons (treatment interaction)

**Data.** Any randomised trial of a checkpoint-containing arm against a non-checkpoint arm, with pre-treatment expression for both arms. The first candidate is CheckMate 025 (nivolumab against everolimus, renal cell carcinoma), if the everolimus arm is in the public Braun et al. 2020 data.

**Model.** score × arm interaction, Cox for time-to-event endpoints and logistic for response. The endpoint is progression-free survival where available, otherwise overall survival; the trial's own primary endpoint is noted.

**Expected interactions (benefit-oriented):**

- NetImmune, Immune and IFN: > 0.
- CAF_ECM: ≤ 0.

**Primary test.** The pooled benefit-oriented NetImmune interaction across all primary randomised comparisons: the 3 already reported plus the new ones. PRINCE stays a sensitivity analysis.

**Also reported.** CAF_ECM pooled, Stouffer *z*, and the structure-matched random-gene null for each new trial.

## E22. Gastric and gastro-oesophageal checkpoint cohorts

**Data.** Any gastric, gastro-oesophageal junction or oesophageal adenocarcinoma cohort treated with a PD-1/PD-L1 inhibitor, with or without chemotherapy, that has pre-treatment expression and response or survival.

**Expected:**

- Immune, IFN and NetImmune: log-OR of response > 0.
- CAF_ECM: ≤ 0.
- Pooled with PRJEB25780 (and with PERFECT as a sensitivity analysis).
- Where MSI/EBV status is available: in MSS/EBV-negative patients, the NetImmune log-OR is > 0.

## Deviations

Each entry is dated; all entries before "Analysis started" were written before any score was related to any of the new outcomes.

**2026-10-04: E19, written after E19 ran.** E19 results do not depend on new data. Two changes from the plan:
- **Panel subset.** In the 22 genes shared by all four panel cohorts, CAF_ECM keeps only one gene, so NetImmune cannot be scored. Each panel platform's own gene content was scored instead (post hoc).
- **Gene-retention rule.** The rule as implemented never dropped six genes (THY1, FN1, COL3A1, CD3D, CD8A, GZMB). The literal rule was added as a post-hoc check: at least 3 listed genes kept per module, and cohorts losing a module drop out.

Details are in results/E19_README.md.

**2026-10-04: E17 data.** GSE163558 has 3 primary tumours and no annotation, so it was replaced as planned by atlases found in the search:
- primary: NGDC OMIX001073 (Sun et al. 2022; 10 patients, 166,533 cells);
- secondary: GSE206785 (Kang et al. 2022; 24 tumours, EPCAM-depleted).

Cell types are taken from the authors' published annotation (major types). This annotation is fixed and independent of our scores, and replaces the planned marker-based re-annotation.
- **Sensitivity analysis.** The T-cell-enriched libraries are dropped. The deposited batch labels are used to identify them, and the discrepancy with the paper is recorded.
- The expectations are unchanged.

**2026-10-04: E18 data.** CPS was transcribed from Kim et al. 2018 (Table 2) and is available for 40 of the 44 patients.
- **Primary:** log(CPS + 1), as a continuous covariate.
- **Secondary:** CPS ≥ 1 as a binary covariate.
- The model is Firth logistic regression: response ~ MSI/EBV + CPS + NetImmune.
- **Primary test:** the penalised likelihood ratio for NetImmune.

**2026-10-04: E20 data.**
- **Primary (overall survival):**
  - GSE13861 (65 patients);
  - GSE26901 (109);
  - GSE29272 (126; GPL96 lacks C1QC, so TAM_SUP is scored on 10 genes);
  - GSE26899 (93), added here. It is the sister series of GSE26901, found during the data search and not listed in the original plan; its outcomes have not been related to any score.
- **Sensitivity only:** GSE57303 (68 patients, overall survival from a third-party re-curation).
- **Excluded:** GSE51105 (no per-sample outcome is public).
- **Duplicate check before pooling:** expression correlation of each new cohort against GSE84437, GSE62254, GSE26253, GSE15459 and GSE34942 where platforms allow. Samples are excluded where the maximum cross-cohort Pearson r exceeds the within-cohort 99.9th percentile.
- **Adjustment:** age and stage where recorded.
- **Primary pooled estimates:**
  - NetImmune across all locked cohorts, the previous 3 plus these 4;
  - all gastric cohorts together.

**2026-10-05: E21 data, written before any score was related to these outcomes.**

Four open randomised trials were found and downloaded (data_download/validation/; scouting report results/E21_E22_data_scout.md).

**Primary new comparisons.** Each is added to the 3 existing primary comparisons, for 6 in total.
- **CheckMate 025** (nivolumab against everolimus, renal cell carcinoma; 120 against 130 patients with RNA). Data are the Braun et al. 2020 supplementary matrix (log2) and clinical table.
  - Endpoint: progression-free survival, primary. Overall survival is secondary (it was the trial's primary endpoint). ORR is not used (5 everolimus responders).
  - Event coding: in this file `PFS_CNSR` and `OS_CNSR` = 1 means an event.
- **JAVELIN Renal 101** (avelumab + axitinib against sunitinib; 354 against 372). Endpoint: progression-free survival, the only one available.
- **NeoTRIP** (GSE319641; carboplatin + nab-paclitaxel ± atezolizumab, triple-negative breast cancer; 119 against 122). Endpoint: pCR.
  - Model: logistic, score × arm, with sequencing batch as a main-effect covariate (batches are uneven across arms). Without batch is a sensitivity analysis.
  - Expression is linear TPM; log2(TPM + 1) is applied once.

**Sensitivity only, like PRINCE: DUTRENEO** (GSE328930; durvalumab + tremelimumab against cisplatin-based chemotherapy, bladder; 28 against 19). It is a sensitivity analysis for three reasons: randomisation was only among immune-high tumours, the checkpoint group includes non-randomised run-in patients who cannot be identified, and many libraries are shallow. Raw counts are converted to log2(CPM + 1); endpoint: pathological response.

**Primary test (unchanged).** The pooled benefit-oriented NetImmune interaction across all primary comparisons, now 6, with the conservative rule. Also reported:
- CAF_ECM pooled, Stouffer z, and pooled estimates without each new trial;
- the structure-matched random-gene null for each new trial;
- whether JAVELIN Renal 101 is primary: the arms also differ in kinase inhibitor (axitinib against sunitinib), as IMmotion150's did (bevacizumab against sunitinib), so it is kept as primary.

**Check without outcomes.** For the CheckMate 025/010/009 nivolumab patients present in both sources, the score from the Braun supplementary matrix is correlated with the score from the compendium object used for the first 11 cohorts (ICB_Braun.rds). This checks the transform of the latter.

**2026-10-05: E22 data, written before any score was related to these outcomes.**
- **GSE183924** (adjuvant durvalumab after chemoradiation and surgery, oesophageal adenocarcinoma and GEJ; 37 patients; FPKM, log2 applied once). The endpoint is relapse-free survival, not response, so it is reported separately and not pooled with the response cohorts.
  - Expected: Immune, IFN and NetImmune HR < 1; CAF_ECM HR ≥ 1.
- **Awaiting the user's decision (FASTQ only):** PRJEB60680 (Samsung NCT04249739, 46 baseline gastric biopsies, about 150 GB), and the MSI-H-only gastric cohorts PRJEB40416 (Kwon) and DRA013565 (Chida).
  - If processed, they follow E22 as written.
  - The MSI-H-only cohorts are analysed within MSI-H patients, expecting NetImmune log-OR > 0.

**2026-10-05: controlled-access datasets, written before applying for access.** These datasets follow E21 and E22 as written. The endpoints are fixed here:
- **CheckMate 649** (nivolumab + chemotherapy against chemotherapy; gastric, GEJ and oesophageal adenocarcinoma; EGAD50000001105). It is analysed under both E21 and E22.
  - Primary: overall survival, the trial's primary endpoint. Secondary: progression-free survival and objective response.
  - Interaction models are stratified by the trial's randomisation factors where available. Analyses are repeated within PD-L1 CPS ≥ 5 and CPS < 5, and within MSS patients.
  - Expected: NetImmune, Immune and IFN interactions > 0 (benefit-oriented); CAF_ECM ≤ 0. Within the nivolumab arm, NetImmune is associated with longer overall survival.
- **IMvigor211** (atezolizumab against chemotherapy, urothelial; EGAS50000000497). Primary: overall survival, the trial's primary endpoint.
- **IMmotion151** (atezolizumab + bevacizumab against sunitinib, renal; EGAS00001004353). Primary: progression-free survival.
- **OAK and POPLAR** (atezolizumab against docetaxel, NSCLC; EGAS00001005013). Primary: overall survival, the trials' primary endpoint. The two trials are analysed separately.
- **IMvigor130:** only the atezolizumab-monotherapy and chemotherapy arms are compared, the checkpoint-without-chemotherapy contrast; overall survival.

Expectations for all of these are as in E21. Each enters the pooled primary test as a new primary comparison, and the pooled result is reported both with and without them.
