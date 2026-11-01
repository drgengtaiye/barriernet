# E28 outcome against the pre-registration (E28_prereg.md, sha256 unchanged)

Recorded 2026-10-21. Step 0 re-derived all 72 published interaction and within-arm values (largest
difference 2.2e-16); the gastric benchmark HRs per SD were reproduced (NetImmune 0.82 / 0.75 / 0.66, pooled
0.74, 0.58–0.95, P = 0.036; ESTIMATE z-mean 0.73, 0.57–0.92, P = 0.029; C as in Supplementary Table S3).
Scripts: E28_base.py, E28_dump.py, E28_estimate.R, E28_trials.py, E28_gastric.R, E28_tiers.py.

| Question | Pre-registered expectation | Result | Met? |
|---|---|---|---|
| Q1a cause of r = −0.50 | missing values in the ranking | The B2b input wrote missing values as "NA"; `estimateScore` reads with `na.strings = ""`, so every GSE84437 column (431/431; 35 genes with missing values, none in the ESTIMATE lists) was read as text and converted to factor codes. Removing those genes or setting them to the gene minimum: official immune score r = +0.89 with PTPRC/CD2/CD3E/CD8A, official contrast r = 0.98 with the z-mean contrast | Yes (mechanism more specific than expected) |
| Q1b official ESTIMATE, six comparisons | r ≥ 0.9 with z-mean; pooled positive | corrected official: r 0.93–0.96; pooled +0.33 (0.08 to 0.58), P = 0.020, 6/6 favorable (z-mean +0.30, P = 0.0089). As-is official, I-SPY2 broken by the same input problem (r 0.66, 0.78) | Yes → Table 4 keeps z-mean, footnote adds official |
| Q1 gastric prognosis (same script) | — | corrected official HR 0.73 (0.58–0.92), P = 0.029; C 0.565 / 0.594 / 0.610 (z-mean 0.569 / 0.586 / 0.611) | — |
| Q2a Immune vs CAF_ECM, adjusted | none | Δ = θ(Immune) + θ(CAF_ECM) = +0.07 (−0.12 to 0.26), P = 0.47; renal +0.11 (P = 0.28); breast −0.33 (P = 0.31) | Decision rule triggered: P ≥ 0.05 |
| Q2b marginal | none | Δ = +0.05 (−0.16 to 0.26), P = 0.64 | same |
| Q3 size of nominal test | 3–7% at α = 0.05 in every comparison | 2.4–6.8% over 30 score × comparison cells; JAVELIN Renal 101 4.5–5.5%; null SD of z 0.92–1.06 | Yes → "unreliable in the largest trial" withdrawn; random-gene distribution renamed a reference distribution |
| Q3 pooled rule | — | null rate of pooled P < 0.05: 1.8–3.3% (conservative). Calibrated pooled P: NetImmune 0.013 (nominal 0.034), Immune 0.009 (0.021), CAF_ECM 0.089 (0.14), ESTIMATE 0.006 (0.0089), CD8A 0.011 (0.035). Stouffer calibrated: 0.002, 0.003, 0.009, 0.001, 0.016 | nominal P stay primary (rule) |
| Q4 tiers | — | Tier 1 (three pre-registered comparisons alone): NetImmune +0.18 (−0.26 to 0.63), P = 0.22, 3/3 favorable, Stouffer P = 0.011; Immune +0.15, P = 0.094; CAF_ECM −0.09, P = 0.48. First three: NetImmune +0.47, P = 0.057. Tier 2 Holm over three axes: Immune 0.063, NetImmune 0.068 (nominal); 0.027, 0.027 (calibrated; Barrier nominal 0.47 used, it was not in Q3). BH over nine: smallest 0.10 (pooled P), 0.0012 for NetImmune (Stouffer) | — |
| Q5 single-sample, trials | r ≥ 0.9 and pooled NetImmune within ±0.05 of +0.29 for (a) and (b) | (a) frozen TCGA reference: r ≥ 0.98, +0.26 (0.02 to 0.50), P = 0.037; (b) unit weights: r ≥ 0.97, +0.28 (0.02 to 0.53), P = 0.041; (c) within-sample rank: r ≥ 0.86, +0.24 (0.07 to 0.41), P = 0.015. CAF_ECM −0.18 / −0.17 / −0.16, 0/6, 0/6, 1/6 favorable | Yes → Methods state independence from cohort-relative standardization |
| Q5 gastric | — | unit weights NetImmune HR 0.74 (0.58–0.95), P = 0.034; frozen TCGA-STAD reference in the two array cohorts 0.74 and 0.66 (cohort-relative 0.75, 0.66); within-sample rank 0.80 (0.55–1.15), P = 0.12, weakest in TCGA (0.92; r 0.87 with cohort-relative) | — |
| Q6 NetImmune vs ESTIMATE | difference near zero | paired difference −0.03 (−0.31 to 0.26), P = 0.81; joint model: NetImmune +0.13 (P = 0.57), ESTIMATE +0.22 (P = 0.30) adjusted for each other | Yes → no superiority claim |
| Q7 PD-L1, JAVELIN (exploratory) | — | NetImmune × arm +0.18 (P = 0.090) → +0.12 (P = 0.27) after PD-L1 and PD-L1 × arm (all 726 patients have PD-L1 status); ESTIMATE +0.32 (P = 0.0037) → +0.28 (P = 0.014) | — |

**Decision rules applied.** Q1b: official ESTIMATE as a Table 4 footnote; SI wording corrected (the negative
correlation belonged to our input to the official routine, not to the platform or to the z-mean contrast).
Q2: no statement that Immune modifies benefit more than CAF_ECM; Δ reported. Q3: the statement that nominal
interaction tests were unreliable in JAVELIN Renal 101 is withdrawn; random gene sets carry real
interaction signal there, and the random-gene distribution is called a structure-matched reference
distribution. Q4: abstract and conclusions give the tier of the randomized-comparison P values. Q5: Methods
state that results do not depend on cohort-relative standardization. Q6: no superiority claim over ESTIMATE.

Corrections (2026-10-21, found by an independent check while integrating E28 into the manuscript; no decision
changes): (1) Q5 unit-weights pooled NetImmune interaction is 0.27499 (results/E28_Q5_single_sample_pooled.csv), so
+0.27, not +0.28; (2) Q1b as-is official r in I-SPY2 is 0.664 and 0.775, so 0.66 and 0.77; (3) Q3 held in 27 of 30
score x comparison cells; three cells fell below 3% (I-SPY2 durvalumab-olaparib: CAF_ECM 2.6%, ESTIMATE 2.4%, CD8A
2.5%), i.e. conservative; JAVELIN Renal 101, the comparison the decision rule names, was 4.5-5.5% for every score.
