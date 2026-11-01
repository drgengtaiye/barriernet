# E30 outcome

Protocol: `E30_prereg.md` (sha256 in `E30_prereg.sha256`, verified before and after both runs). Scripts:
`E30_dedup*.R`, `E30_dedup_*.py` (Q1, Q2 gastric) and `E30_trials*.py|R` (Q2 randomized, Q3–Q6). Results in
`results/E30_*`. No published script or result file was changed.

## Gates

- Duplicate branch: with the filter off, all 29 published result files of E15, E20, E29 and E29b were reproduced
  exactly (largest difference 0, including null draws and bootstraps; `results/E30_dedup_gate.csv`).
- Randomized branch: 112 per-comparison interactions (largest difference 4.0 × 10⁻¹⁵), 80 pooled results (5.3 × 10⁻¹⁵)
  and all 13 rows of main-text Table 4 reproduced (`results/E30_trials_gate.csv`, `E30_trials_meta_gate.csv`).

## Q1. Duplicates — expectation NOT held

- The rule selected exactly the eight expected arrays of GSE26253 (GSM644593, GSM644726, GSM644733, GSM644801,
  GSM644802, GSM644810, GSM644871, GSM644940); the three-array sensitivity set is GSM644726, GSM644801, GSM644802.
  Of the 11 mutual best matches agreeing on sex and age, three failed only on time (GSM644782, GSM644832, GSM644878).
- Only E15, E20, E29 and E29b use GSE26253; each was re-run with a one-line sample filter.
- Expected: no reported pooled estimate changes by more than 0.01; no conclusion changes. **Not held:**
  - The first three locked cohorts, NetImmune: 0.73 (0.53–1.02; P = 0.057) → 0.72 (0.53–1.00; P = 0.049); the
    conservative interval now excludes 1 (three-array set: P = 0.048). The extension was decided on P = 0.057.
  - The first three locked cohorts, CAF_ECM: pooled HR 1.258 → 1.269 (plotted in Fig. 3d of v31).
  - Every point estimate quoted in the v31 main text moved by less than 0.01 (largest 0.0097).
- Main changed values: GSE26253 424 patients (174 recurrences), NetImmune HR 0.73 (0.63–0.86); seven locked cohorts
  0.70 (0.62–0.79; P = 2.9 × 10⁻⁴; I² = 1%); all ten 0.72 (0.67–0.78), CAF_ECM 1.31 (1.17–1.46); ΔC against CD8A
  +0.063 (0.034–0.093), against ESTIMATE +0.015 (−0.002 to 0.032); post hoc mutually adjusted Immune HR 0.82
  (0.75–0.90), CAF_ECM 1.33 (1.21–1.46). The E29 survival contrast (−0.08, P = 0.22) is unchanged.
- *Rule applied.* Duplicate-free values are primary in the main text and Table 1; Fig. 3 (v32 Fig. 2), Fig. 8
  (v32 Fig. 6) and Supplementary Fig. S10 were re-rendered from the duplicate-free results; the previous values are
  kept in the Supplementary Tables (S7, with the duplicate-free versions in S19b). The Methods describe the
  cross-platform screen; the Supplementary Methods statement that no duplicate was found between GSE26253 and
  GSE62254 is withdrawn. The text states that the extension was decided on the first result (P = 0.057) and that
  the duplicate-free value is P = 0.049.

## Q2. Estimators (disclosure)

| NetImmune pool | DL | REML | DL-HK | REML-HK | rule |
|---|---|---|---|---|---|
| Locked first three, with duplicates | 3.3 × 10⁻⁵ | — | 0.057 | 0.047 | 0.057 |
| Locked first three, without duplicates | 7.1 × 10⁻⁶ | — | 0.049 | 0.041 | 0.049 |
| Randomized, first three | 0.0020 | 0.0020 | 0.057 | 0.057 | 0.057 |
| Randomized, three added | 0.074 | 0.027 | 0.25 | 0.22 | 0.22 |
| Randomized, all six | 0.0033 | 0.0037 | 0.034 | 0.034 | 0.034 |

Full tables: `results/E30_estimators_gastric.csv`, `results/E30_trials_estimators.csv`.

## Q3. Angiogenesis adjustment in the renal trials

- Expectation 1 (each adjusted renal interaction within 0.05 of the unadjusted one): **not held** — 13 of 24
  score × trial pairs. NetImmune stayed within 0.05 in all three trials (+0.006, −0.013, −0.031). CAF_ECM moved
  from −0.19 to +0.01 (IMmotion150), −0.02 to +0.00 (CheckMate 025) and −0.06 to +0.05 (JAVELIN Renal 101);
  Immune from +0.39 to +0.55 in IMmotion150. The expectation rested on Spearman ρ ≤ 0.14 between Angio and the
  scores; within arms the Pearson correlation of Angio with CAF_ECM was 0.19 to 0.36 in IMmotion150 and JAVELIN
  Renal 101, and ESTIMATE and TMEscore correlated −0.20 to −0.30 with Angio.
- Expectation 2 (six-comparison NetImmune interaction keeps its sign with a CI above zero): **held** —
  +0.29 (0.01 to 0.56; P = 0.043), against +0.29 (0.03 to 0.55; P = 0.034) unadjusted.
- *Rule applied (branch 1):* the main text states that the randomized NetImmune result was not explained by
  angiogenesis. Not covered by the rule but reported: the fibroblast module's adverse renal interactions disappear
  after adjustment (renal pool −0.07 → +0.02), so its signal against checkpoint benefit rests on the breast trials.
- Other six-comparison pools with the renal trials adjusted: Immune +0.25 (P = 0.033), CAF_ECM −0.14 (P = 0.33),
  ESTIMATE +0.26 (P = 0.023), TMEscore +0.23 and +0.24 (P = 0.075 and 0.067), CD8A +0.25 (P = 0.050).
  Angio × arm pooled over the renal trials (Angio alone): −0.27 (−0.56 to 0.03; P = 0.059). IMmotion150
  atezolizumab monotherapy against sunitinib (descriptive): NetImmune +0.28 → +0.39.

## Q4. Endpoint strata (presentation)

NetImmune renal +0.20 (−0.17 to 0.57; P = 0.15), breast +0.63 (0.02 to 1.25; P = 0.048); CAF_ECM renal −0.07
(P = 0.33), breast −0.60 (−1.09 to −0.11; P = 0.034). Reported side by side; magnitudes not compared.

## Q5. Stouffer

NetImmune unweighted z 3.82 (P = 1.3 × 10⁻⁴), √n-weighted 3.86 (P = 1.1 × 10⁻⁴), calibrated P 0.002 for both.
Moved from Table 4 to the Supplementary Tables.

## Q6. CheckMate 025 with overall survival

No six-comparison pooled interaction remains significant (smallest P, NetImmune 0.059; Immune 0.25; IFN 0.18;
ESTIMATE 0.107; CD8A 0.133). Reported in the Supplementary Tables and in the main text.

## Count

Three expectations: one held (Q3, expectation 2), two not held (Q1; Q3, expectation 1). With E28 and E29, the
three later protocols carry ten expectations, of which six held fully.
