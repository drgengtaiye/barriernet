# E44 outcome: gastric and gastroesophageal applications of the calibration framework

Protocol E44_prereg.md (sha256 f6c0589648ed…, logged 2026-11-03T02:19:20Z, before E44_export.R and E44_gastric.py were
written). Seventeenth protocol written after the other results were known. Run: E44_export.R → results/E44_export.log,
E44_Z_*.csv.gz, E44_labels.csv; E44_gastric.py → results/E44.log, E44_A_gastric_trials.csv, E44_B_per_cohort.csv,
E44_B_pooled.csv, E44_B_paired.csv, E44_B_adjusted.csv, E44_C_prognosis.csv, E44_D_tmescore.csv, E44_expectations.csv.

Step 0: the Immune axis log-odds per SD from the exported matrices reproduced E12 (PRJEB25780), E22 (PRJEB60680) and
E16 (PERFECT) to 1.8 × 10⁻¹⁵. All 52 signatures, the 27 panel genes, ESTIMATE, TMEscore, the Immune axis, NetImmune and
barriernet 0.3 were scored in all three cohorts (barriernet 0.3: 48, 48 and 47 of its 48 genes). No bootstrap resample
failed (252 series × 1,000).

## Expectations: 4 of 4 held

| | expectation | result | held |
|---|---|---|---|
| G1 (primary) | CD8A pooled log-odds > 0, P < 0.05 (fixed effect, three cohorts) | +0.58 (0.12 to 1.04), P = 0.014 | yes |
| G2 | no signature significantly above CD8A | 0 above; 10 below | yes |
| G3 | Spearman with the randomized interaction > 0, P < 0.05 | 0.76 (P = 6 × 10⁻¹¹; 52 signatures) | yes |
| G4 | barriernet 0.3 > 0 with P < 0.05, not different from CD8A | +0.67 (0.21 to 1.14), P = 0.0045; minus CD8A −0.07 (−0.74 to 0.60), P = 0.84 | yes |

G1 and G4 were close to predictable: the Immune axis values of the three cohorts were known (combined by hand to about
+0.59 while planning), and CD8A and barriernet 0.3 track the Immune axis.

## A. Published randomized gastric trials (descriptive)

| | per SD (95% CI) | ratio of HRs, high vs low | detectable at 80% power |
|---|---|---|---|
| CheckMate 649 (nivolumab + ipilimumab vs chemotherapy; n = 402) | +0.17 (−0.07 to 0.40) | 0.69 (0.42–1.15) | 0.33 |
| KEYNOTE-061 (pembrolizumab vs paclitaxel; n = 274) | −0.04 (−0.37 to 0.30) | 1.06 (0.61–1.84) | 0.48 |
| both, fixed effect | +0.10 (−0.09 to 0.29), P = 0.30 | 0.84 (0.58–1.23) | 0.28 |
| other 16 E41 trials (DerSimonian–Laird) | +0.16 (0.06 to 0.26), P = 0.001 | | 0.14 |

- Between the two gastric trials: Q = 0.95, P = 0.33. Gastric minus the other 16: −0.06 (−0.27 to 0.16), P = 0.61.
- The detectable interaction of the gastric biomarker subsets (0.28–0.48 per SD) is larger than the six-comparison
  interaction of CD8A (+0.23), so these trials could not show even CD8A's interaction reliably, let alone a 0.10
  advantage over it.

## B. Single-arm gastroesophageal checkpoint cohorts

Per cohort (log-odds of response per SD):

| | PRJEB25780 (pembrolizumab; 44, 12 responders) | PRJEB60680 (chemotherapy then pembrolizumab; 43, 29) | PERFECT (atezolizumab + chemoradiotherapy; 32, 12) |
|---|---|---|---|
| CD8A | +1.56 (P = 0.005) | +0.49 (0.17) | +0.24 (0.53) |
| barriernet 0.3 | +1.66 (0.002) | +0.30 (0.37) | +0.65 (0.12) |
| Immune axis | +1.32 (0.009) | +0.40 (0.25) | +0.40 (0.30) |

Pooled, three cohorts (fixed effect):
- CD8A rank 20 of 79 (52 signatures and 27 panel genes). Highest: CXCL10 +0.89, CXCL13 +0.89, TMEscore +0.82,
  Teff_McDermott +0.75, IFN_Ayers +0.75.
- 21 of 52 signatures at P < 0.05 (all positive). Of these, 13 also had a same-direction randomized interaction at
  P < 0.05 (62%; E40 single-arm figure: half).
- Significantly below CD8A (paired): ADO_Sidders, Bcell_Budczies, Bcell_Helmink, CD8_Sade-Feldman, CYT_Davoli,
  EMT_Thompson, KDM5A_Wang, Myeloid-DC_Helmink, Response_Chen, TLS_Cabrita. Largest difference above CD8A:
  Chemokine_Messina +0.13 (−0.62 to 0.88).
- Under the manuscript's conservative rule (Hartung–Knapp with k = 3) CD8A (+0.65, P = 0.22) and barriernet 0.3
  (P = 0.18) are not significant; two signatures, CIN25_Carter (P = 0.028) and proliferation_Pabla (0.041), stay at
  P < 0.05 (their Q is small, so the rule falls back to DerSimonian–Laird), and CTLA4 in the two-gastric set (0.049).
  (Corrected after the v46 fact-check; the first version said nothing was significant.) The fixed-effect primary was chosen in the protocol because three cohorts cannot estimate
  between-cohort variance; the conservative values are reported beside it.
- Gastric only (two cohorts): CD8A +0.81 (0.21 to 1.40), P = 0.008, rank 9 of 79; barriernet 0.3 +0.68, P = 0.018,
  minus CD8A −0.12 (P = 0.74); 0 above, 11 below; 17 of 52 at P < 0.05; ρ = 0.79.
- MSI-H/EBV-positive status (two gastric cohorts; 11 of 44 and 3 of 43 patients) explains much of the signal: CD8A
  +0.81 → +0.44 (−0.22 to 1.10), P = 0.19; barriernet 0.3 +0.68 → +0.40 (P = 0.20); TMEscore +0.65 → +0.22
  (P = 0.55); signatures at P < 0.05: 17 → 1 of 52.
- PRJEB25780 has the smallest inverse-variance weight but the largest estimates, and supplies about half of each pooled
  estimate (see E44b below). The protocol states that it was among the cohorts explored when the 48 barriernet genes
  were fixed; the project's file records date its first analysis after those genes were fixed (see E44b).
  (Corrected after the v46 fact-check; the first version said it carried most of the weight.)

## C. Prognosis without checkpoint inhibitors (E20, 10 gastric cohorts, report only)

Hazard ratio per SD: CD8A 0.90 (0.84–0.96), P = 0.005; Immune axis 0.84 (0.78–0.91); T_CELL 0.86; IFN 0.86; ESTIMATE
0.73 (0.67–0.80). The immune scores that track response in single-arm gastric cohorts also track survival in patients
who never received a checkpoint inhibitor.

## D. TMEscore (developed in gastric cancer) in the six randomized comparisons (E32, report only)

Interaction +0.31 per SD (0.04 to 0.57), P = 0.030; above its size-matched random reference (empirical P = 0.001);
minus CD8A +0.06 (−0.09 to 0.20), P = 0.44.

## Reading

In gastric and gastroesophageal cancer the pattern of the randomized comparisons repeats: inflamed scores track
response, no published signature beats CD8A, and barriernet 0.3 equals CD8A. But none of this is randomized. The
gastric immune scores are prognostic without checkpoint inhibitors, much of their response signal is MSI-H/EBV status,
and the two randomized gastric trials that report signature subgroups disagree and were too small to detect an
interaction of CD8A's size. The gastric question is the one the framework says cannot yet be answered.

## E44b addendum (post hoc, after the independent fact-check of the v46 text)

Protocol E44b_prereg.md (sha256 090045882f3c…, logged 2026-11-03T15:51:32Z, before E44b_checks.py was written). No
expectations; not counted in the tally. Outputs results/E44b_checks.csv, E44b.log.

- Weights and leave-one-out (three cohorts, fixed effect). PRJEB25780: weight 18% (CD8A), 19% (barriernet 0.3), 20%
  (Immune axis), but 48%, 47% and 46% of the pooled estimates. Without it: CD8A +0.37 (−0.14 to 0.88), P = 0.16;
  barriernet 0.3 +0.44 (−0.07 to 0.96), P = 0.093; Immune axis +0.40, P = 0.12. Without either other cohort, both stay
  significant.
- Heterogeneity: I² 49% (CD8A, Q P = 0.14), 55% (barriernet 0.3, 0.11), 25% (Immune axis), 72% (TMEscore).
- Main-text rule: 2 of 52 significant in the three-cohort set (CIN25_Carter, proliferation_Pabla, both in the
  confounder family), none in the two-gastric set (CTLA4 among the genes).
- Power for CD8A's six-comparison interaction (+0.226 per SD): CheckMate 649 47%, KEYNOTE-061 26%, both 63%. "Could not
  have detected" overstated this; "less than 50% power each" is accurate.
- File dates (UTC): barriernet_0.2.0.zip 2026-02-12 01:43; earliest PRJEB25780/Kim2018 file 2026-03-25 02:23
  (GC_TME_Reprogrammer/validation). The protocol's statement that PRJEB25780 informed the gene choice was probably
  over-cautious; the author is asked to confirm.
- E44.log prints "MSI-H/EBV+ 0" for PERFECT because missing status was summed as 0; PERFECT has no recorded status.
