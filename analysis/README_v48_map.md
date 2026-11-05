# barriernet — map of the current manuscript to scripts and outputs

This file maps every table, figure panel and reported number of the current manuscript to the script that
produces it and to that script's output file. The manuscript text is `manuscript/v48_text.py`; the built package is
`final_v48/`. `README.md` in this folder maps the v14 paper; rows here point to it only where nothing has changed.
Rows for E38–E40 (Fig. 2g–i) are in section 8; the later protocols E41–E44 are in sections 4, 6 and 9, and E45–E46 in sections 4, 6 and 10.

**Counts.** The manuscript counts 51 distinct signatures, because the rows `TIS_Damotte` and `Tcell-inflamed_Ayers` hold the same 18 genes. The result files keep all 52 rows. `A/E46_E_distinct.csv` (from `E46_distinct.py`) gives every count in the text both ways. Rows in section 4 outside the abstract, the key messages, the gastric section and section 10 give the 52-row value of the file.

**Conventions**

- `A/` is the lab results folder `gastric/2/analysis/results/`. The scripts read and write absolute lab paths,
  so these must be changed for a local copy.
- Scripts are given by their repository path:
  - E23–E40 and `cmp_scores.py`: `analysis/15_v19_v41/scripts/`;
  - E41–E44: `analysis/16_v42_v46/scripts/`;
  - E45–E46: `analysis/17_v47_v48/scripts/`;
  - protocols and outcomes: `analysis/15_v19_v41/protocols/` (E41–E44: `analysis/16_v42_v46/protocols/`; E45–E46: `analysis/17_v47_v48/protocols/`);
  - figure scripts: `analysis/15_v19_v41/figure_scripts/`, with the same layout as `final_v41/figure_scripts/`; versions v42–v46 and the current `cns_style.py`: `analysis/16_v42_v46/figure_scripts/`; versions v47–v48: `analysis/17_v47_v48/figure_scripts/`;
  - earlier scripts: their stage folders (`analysis/05_v14/`, `09_v15/`, `10_v16/`, `11_v17/`, `14_v18/`, …).
- Build scripts are named by their lab path, `manuscript/<file>`. They contain the manuscript text, so they are provided
  with the submission (package folder `build_v48/`) and are not in this repository before publication.
- **How each mapping was checked.**
  - "Written by": the write call (`to_csv`, `write.csv`, `write.table`, `savez`) was found in that script.
  - "Input": the read call was found in the figure or build script.
  - Package figure files were matched to their lab renders by MD5.
- **"not established"**: no script on disk was found that writes the file or prints the number.

---

## 1. Run order of the randomized-trial analyses

The raw trial data are read by the loaders of E5, E10 and E21. Every later script either copies those loaders or
executes part of an earlier script with `exec(compile(...))`, after rewriting its step-0 and coverage file names.
Every executed copy re-runs the step-0 gate, which stops if any published interaction or within-arm value is not
reproduced. Running a script overwrites its `A/` outputs.

| Step | Script | Setup | Reads from `A/` | Writes to `A/` |
|---|---|---|---|---|
| 1 | `analysis/09_v15/E5_trials.py` | own loaders (IMmotion150, PRINCE); imports `analysis/15_v19_v41/scripts/cmp_scores.py` | `E23_comparator_genes.csv` (via `cmp_scores.py`) | `E5_trials.csv` |
| 2 | `analysis/10_v16/E10_ispy2_trials.py` | own loader (the two I-SPY2 comparisons); imports `cmp_scores.py` | as 1 | `E10_ispy2_interactions.csv`, `E10_ispy2_within_arm.csv`, `E10_ispy2_hr_strata.csv`, `E10_ispy2_null_draws.csv`, `E10_ispy2_null_summary.csv` |
| 3 | `analysis/14_v18/E21_trials.py` | own loaders (CheckMate 025, JAVELIN Renal 101, NeoTRIP, DUTRENEO); imports `cmp_scores.py` | as 1 | `E21_interactions.csv`, `E21_within_arm.csv`, `E21_scores.csv`, `E21_null_draws.csv`, `E21_null_summary.csv`, `E21_expectations_per_trial.csv` |
| 4 | `analysis/10_v16/E10_meta.R`; `analysis/14_v18/E21_meta.R` | pooling (metafor) | outputs of 1–3 | `E10_meta.csv`; `E21_meta_inputs.csv`, `E21_meta.csv`, `E21_expectations_pooled.csv` |
| 5 | `analysis/15_v19_v41/scripts/E25_trials.py` → `E25_meta.R` | copies the E5/E10/E21 loaders and models | `E21_meta.csv` (reproduction) | `E25_*.csv` |
| 6 | `analysis/15_v19_v41/scripts/E28_base.py` (never run alone) | copy of `E25_trials.py` lines 16–307 (data, scores, models, pooling rule `pool_rule`); step 0 | `E5_trials.csv`, `E10_ispy2_interactions.csv`, `E10_ispy2_within_arm.csv`, `E21_interactions.csv`, `E21_within_arm.csv` | step-0 table `E28_reproduction_check.csv` (each caller renames it) |
| 7 | `E28_trials.py` (`from E28_base import *`); `E28_dump.py` → `E28_estimate.R`; `E28_gastric.R`; `E28_tiers.py` | E28_base | `E21_meta.csv`, `E28_Q3_calibration_pooled.csv` (both read by `E28_tiers.py`) | `E28_Q*.csv`, `E28_bootstrap_draws.csv`, `E28_gastric_*.csv`, `E28_Q1_estimate_audit.csv`, `E28_Q4_tiers.csv` |
| 8 | `E30_trials.py` → `E30_trials_meta.R` → `E30_trials_report.py` | `exec` of the whole `E28_base.py` | — | `E30_trials_*.csv` |
| 9 | `E32_trials.py` | `exec` of the whole `E28_base.py` (step 0 → `E32_step0_reproduction.csv`); checks the protocol hash | `E23_comparator_genes.csv` | `E32_coverage.csv`, `E32_per_comparison.csv`, `E32_reference_draws.npz`, `E32_pooled_all_sets.csv`, `E32_primary.csv`, `E32_ispy2.csv`, `E32_expectations.csv`; stdout `E32.log` |
| 10 | `E32_native.py` (calls `Rscript E32_native_scores.R`) | `exec` of `E28_base.py` | — | `E32_native.csv`, `E32_native_per_comparison.csv` |
| 11 | `E32_locked_bionull.R` | standalone (gastric locked cohorts) | `E30_dedup_runs/primary/E20_estimates_used_for_pooling.csv`, `E8_biology_null_draws.rds`, E15/E20 score objects | `E32_locked_bionull.csv`, `E32_locked_bionull_draws.rds` |
| 12 | `E33_trials.py` | `exec` of `E32_trials.py` up to its `Q1:` section (step 0 and coverage → `E33_*`) | `E32_per_comparison.csv`, `E32_primary.csv`, `E32_reference_draws.npz` | `E33_per_comparison.csv`, `E33_pooled.csv`, `E33_paired.csv`, `E33_Q3_by_reference.csv`, `E33_pooled_reference_estimates.npz`, `E33_expectations.csv` |
| 13 | `E33_extras.py` | `exec` of the pooling functions of `E28_base.py` only; no new fits | `E32_per_comparison.csv`, `E32_coverage.csv` | `E33_breast_extras.csv`, `E33_trial_table.csv` |
| 14 | `E34_link.py` | standalone; uses expression values only | NeoPembrOv GEO counts and source data | `data_download/validation/NeoPembrOV_GSE227666/id_match_final.csv`, `id_match_checks.csv` (outside `A/`) |
| 15 | `E34_trials.py` | `exec` of `E32_trials.py` up to `Q1:`, plus NeoPembrOv | `id_match_final.csv`, `E32_per_comparison.csv`, `E33_per_comparison.csv`, `E33_paired.csv` | `E34_neopembrov_per_score.csv`, `E34_pooled.csv`, `E34_expectations.csv` |
| 16 | `E35_trials.py` | `exec` of `E34_trials.py` up to `def fit(task):` | `E33_paired.csv` (reproduction check) | `E35_paired.csv`, `E35_by_gene.csv`, `E35_expectations.csv` |
| 17 | `E36_sim.py` | `exec` of `E34_trials.py` up to `def fit(task):` | — | `E36_S0_S2.csv`, `E36_S3.csv`, `E36_S4.csv`, `E36_expectations.csv` |
| 18a | `E37_export.R` | own loaders (E5/E12 logic, without ICB_Braun) | single-arm compendium | `E37_dev_genes.csv`, `E37_twin_*.csv.gz`, `E37_step0_reproduction.csv` (later overwritten, see 18d) |
| 18b | `E37_twins.py` | `exec` of the data part of `E28_base.py` | `E37_twin_*.csv.gz` | `E37_twins.csv`, `E37_twin_flags.csv` |
| 18c | `E37_develop.R` | — | `E37_dev_genes.csv`, `E37_twin_flags.csv` | `E37_loco_per_cohort.csv`, `E37_loco_pooled.csv`, `E37_C2_module_coefficients.csv`, `E37_dev_expectations.csv`, `E37_frozen_weights.tsv` (hash in `analysis/15_v19_v41/protocols/E37_frozen.sha256`) |
| 18d | `E37_test.py`, then `E37_refdraws.py` | `exec` of `E34_trials.py` up to `def fit(task):` (step 0 and coverage → `E37_step0_reproduction.csv`, `E37_coverage_check.csv`); check the weights hash | `E37_frozen_weights.tsv`, `E32_primary.csv` | `E37_test_per_comparison.csv`, `E37_test_pooled.csv`, `E37_test_expectations.csv`; `E37_reference_draws.csv` |
| 18e | `E37_panel_extra.py` (added after the test); `E37_package_check.py` | `exec` of `E34_trials.py` (`E37x_*`, `E37p_*` step-0 names); the package check imports the repository's `barriernet` 0.3.0 | `E37_frozen_weights.tsv`, `E37_test_*.csv` | `E37_panel_extra.csv`; `E37_package_check.csv` |
| 19 | `E39_design.py`; `E40_export.R` → `E40_concordance.py` | `E39_design.py` and `E40_concordance.py` `exec` `E34_trials.py` up to `def fit(task):` | — | mapped separately (E38–E40) |
| 20 | `E45_export.R` → `E45_select.py` → `E45_test.py` | `E45_test.py` `exec`s `E34_trials.py` up to `def fit(task):`; checks the protocol, frozen-gene and frozen-weight hashes | `E45_reference_genes.tsv`, `E32_primary.csv`, `E33_pooled.csv`, `E37_frozen_weights.tsv`, `E44_Z_*`, `E44_labels.csv`, `E44_B_per_cohort.csv` | section 10 |
| 21 | `E46_test.py` → `E46_distinct.py`; `E46b_ispy2.py`; `E46c_pembro_pairs.py` | `exec` of `E34_trials.py` up to `def fit(task):`; each checks its protocol hash; `E46_distinct.py` fits nothing | `E32_primary.csv`, `E33_pooled.csv`, `E44_Z_*`, `E44_labels.csv`, `E44_B_per_cohort.csv`, `E46_D_ispy2_flags.csv` | section 10 |

**Gastric and single-arm results used in the current text.** The upstream steps are in `README.md`:

- `analysis/09_v15/E5_icb_split.R` and `analysis/10_v16/E12_icb_expanded.R`: the checkpoint cohorts;
- `analysis/10_v16/E11_pancancer_prognosis.R`: TCGA;
- `analysis/11_v17/E15_gastric_validation.R` → `analysis/14_v18/E20_extended_replication.R`: the locked gastric
  cohorts.

The duplicate-free values are re-run by the E30–E31 scripts below:

- `E30_dedup.R` (duplicate screen), then, for each `E30_MODE=off|primary|sens3`: `E30_dedup_E15.R` →
  `E30_dedup_E20.R` → `E30_dedup_E29.R` → `E30_dedup_E29b.R`. These write to `A/E30_dedup_runs/<mode>/` under the
  original file names.
- `E30_dedup_compare.py`, `E30_dedup_estimators.R`, `E30_dedup_gate.py` and `E30_dedup_summary.py` compare the runs.
- `E31_dedup_all_pairs.R` → `E31_rerun.R`.
- `E29_dissociation.R` → `E29b_posthoc_adjusted.R` hold the direct Immune-versus-CAF_ECM tests, as first run
  (`A/E29_*`).

---

## 2. Tables

Both tables are filled by `manuscript/apply_v48.py`. Its starting document is `final_v33/BarrierNet_Manuscript.docx`.

| Item | Columns | Script | Output read |
|---|---|---|---|
| Table 1, six primary rows, PRINCE and DUTRENEO | Trial, cancer, arms, endpoint, n | carried as text from v33 Table 3; built by `analysis/14_v18/make_v18_tables.py` (`TRIALS`, `trial_effect`) | n = `n_ici` / `n_ctrl` from `A/E5_trials.csv`, `A/E10_ispy2_interactions.csv`, `A/E21_interactions.csv` |
| Table 1, same rows | Signatures scored (the file's count minus 1: the same 18 genes count once); CD8A interaction, given as exp(∓y) | `analysis/15_v19_v41/scripts/E33_extras.py` | `A/E33_trial_table.csv` (from `E32_coverage.csv` and `E32_per_comparison.csv`) |
| Table 1, NeoPembrOv row | CD8A interaction | `analysis/15_v19_v41/scripts/E34_trials.py` | `A/E34_neopembrov_per_score.csv` (`kind == gene`, `CD8A`, `primary`) |
| Table 1, NeoPembrOv row | n "32 / 21"; "51" signatures | typed in `apply_v48.py` | consistent with `analysis/15_v19_v41/protocols/E34_outcome.md` |
| Table 2: ESTIMATE, TMEscore (z-mean), NetImmune, Immune axis, CAF_ECM | Six comparisons; favoring; reference P (`emp_p_declared`); minus CD8A | `analysis/15_v19_v41/scripts/E32_trials.py` | `A/E32_primary.csv` |
| Table 2, same rows | Breast comparisons | `E33_extras.py` | `A/E33_breast_extras.csv` |
| Table 2: barriernet 0.3 | Six and breast pools; reference P; minus CD8A | `analysis/15_v19_v41/scripts/E37_test.py` | `A/E37_test_pooled.csv` (`six (primary)`, `breast three`) |
| Table 2: barriernet 0.3 | Comparisons favoring | `E37_test.py` | `A/E37_test_per_comparison.csv` (`BN03`, NeoPembrOv excluded) |
| Table 2: CD8A, PRF1, NKG7, CD8B, GZMB, STAT1, CXCL13, CXCL9 | Six and breast pools; favoring | `analysis/15_v19_v41/scripts/E33_trials.py` | `A/E33_pooled.csv` (`six (primary)` and `breast three`, `kind == gene`) |
| Table 2: CD8A, PRF1, GZMA mean; 13 T-cell gene mean | Six and breast pools; favoring; percentile; minus CD8A | `analysis/17_v47_v48/scripts/E46_test.py` | `A/E46_A_reference_pooled.csv`, `A/E46_A_counts.csv` |
| Table 2: CXCL10 | Minus CD8A | `analysis/17_v47_v48/scripts/E45_test.py` | `A/E45_A_counts.csv` (`CXCL10 minus CD8A`); the other columns as the single genes |
| Table 2, single genes | Percentile among the 2,082 genes | computed in `apply_v48.py`: floor(100 × share of `all_six` genes with a lower estimate) | `A/E33_pooled.csv` |
| Supplementary Table S24 (SI document; cited in the Table 1 footnote) | NetImmune and CAF_ECM per comparison, I-SPY2 adjusted for hormone receptor only | `analysis/14_v18/make_v18_tables.py` `table3()`; carried as v33 Table 3 → v34 Table 3 → S24 | `A/E5_trials.csv`, `A/E10_ispy2_interactions.csv`, `A/E21_interactions.csv` |
| Supplementary Tables S22, S23 (SI document) | gastric survival; components | v33 Tables 1–2: `make_v18_tables.py` `table1_rows()` (E20 files), rows patched in `manuscript/apply_v32.py` from `manuscript/v32_fills.json` and in `apply_v33.py` | `v32_fills.json` holds typed E30 values; its writer, and the full source of S23, not established |

---

## 3. Figures

**Numbering.** Figs 6 and 7 of v46 are Supplementary Figs S37 and S38 here. Rows below and in later sections that cite panels of those figures read accordingly. All v48 figure scripts are written by `analysis/17_v47_v48/figure_scripts/v48/make_v48_scripts.py` from the earlier scripts named in each heading; it removes `TIS_Damotte` so that counts are of 51 signatures.

**Fig. 1** is a BioRender schematic. `analysis/15_v19_v41/figure_scripts/v33_fig1/patch_fig1_v33.py` changed
"2,152" to "2,144" on the supplied image, `analysis/17_v47_v48/figure_scripts/v47/patch_fig1_v47.py` reworded panel d and `analysis/17_v47_v48/figure_scripts/v48/patch_fig1_v48.py` reworded the titles of panels b, c and d. The patient counts in the legend are typed: 1,828 is in `A/E32.log`; the
source of the other legend counts is not established.

### Fig. 2a–f: `analysis/17_v47_v48/figure_scripts/v48/make_fig2.py` (from `v43/make_fig2.py`; panels b–f as in `v40/make_fig2.py`; g–i in section 8)

Lab render: `figures/cns/v48/Fig2.png`. Numbers file: `analysis/17_v47_v48/figure_scripts/v48/Fig2_numbers.md`. Inputs of panels b–f are
written by `analysis/15_v19_v41/scripts/E36_sim.py`; the protocol is `analysis/15_v19_v41/protocols/E36_prereg.md`.

| Panel | Content | Input |
|---|---|---|
| a | the three references applied to the T-cell-inflamed gene list | `A/E32_primary.csv`; `A/E40_concordance.csv`; `A/E33_pooled_reference_estimates.npz`; `A/E33_paired.csv`; `A/E46_A_paired.csv` (T-cell mean) |
| b | null and prognostic-only rejection | `A/E36_S0_S2.csv` (S0, S1) |
| c | power, CD8A truth | `A/E36_S0_S2.csv` (S2 β 0.1–0.3) |
| d | random-reference exceedance | `A/E36_S3.csv` |
| e | paired test, CD8A truth | `A/E36_S4.csv` (`S4 truth CD8A`) |
| f | paired-test power against normal theory | `A/E36_S4.csv`. The curves for SE 0.029, 0.062 and 0.076 are typed in the script: the first two are the `median_se` of T_CELL and APM in `E36_S4.csv`; 0.076 is `A/E33_expectations.csv` Q4. The predicted power, Φ(Δ/SE − 1.96), is computed in the script |

### Fig. 3a–f: `analysis/17_v47_v48/figure_scripts/v48/make_fig3.py` (from `v47/make_fig3.py`)

Lab render: `figures/cns/v48/Fig3.png`. Numbers file: `analysis/17_v47_v48/figure_scripts/v48/Fig3_numbers.md`. Panel b also reads `A/E45_B_coexpression.csv` and `A/E46_A_counts.csv`.

| Panel | Content | Input | Written by |
|---|---|---|---|
| a | six comparisons, patients per arm | typed in the script (`TRIALS`); the same counts are in `A/E32.log` | `E32_trials.py` (log) |
| b | counts over the 52 signatures | `A/E32_primary.csv` | `analysis/15_v19_v41/scripts/E32_trials.py` |
| c | pooled interaction by family; composites; CD8A band | `A/E32_primary.csv` | `E32_trials.py` |
| d | each signature against the 95th percentile of its random reference | `A/E32_primary.csv` (`ref_p95_declared`) | `E32_trials.py` |
| e | paired difference from CD8A; detectable-difference band | `A/E32_primary.csv`; `A/E33_expectations.csv` (Q4) | `E32_trials.py`; `analysis/15_v19_v41/scripts/E33_trials.py` |
| f | each comparison against random sets | `A/E32_per_comparison.csv`, `A/E32_reference_draws.npz` | `E32_trials.py` |

### Fig. 4a–f: `analysis/17_v47_v48/figure_scripts/v48/make_fig4.py` (from `v43/make_fig4.py`)

Lab render: `figures/cns/v48/Fig4.png`. Numbers file: `analysis/17_v47_v48/figure_scripts/v48/Fig4_numbers.md`. The counts of signatures above each gene are recomputed without `TIS_Damotte` from `A/E33_paired.csv` (panels b, d) and `A/E35_paired.csv` (panel f).

| Panel | Content | Input | Written by |
|---|---|---|---|
| a | 2,082 genes; panel genes and signatures marked | `A/E33_pooled.csv` (`six (primary)`, `all_six`); `A/E32_primary.csv` | `analysis/15_v19_v41/scripts/E33_trials.py`; `E32_trials.py` |
| b | 27 panel genes; signatures above each | `A/E33_pooled.csv`; `A/E33_Q3_by_reference.csv` | `E33_trials.py` |
| c | per comparison and breast pool | `A/E33_per_comparison.csv`; `A/E33_pooled.csv` (`breast three`) | `E33_trials.py` |
| d | signatures above each gene against the gene's interaction | `A/E33_Q3_by_reference.csv`; `A/E33_pooled.csv` | `E33_trials.py` |
| e | CD8A's rank among the 13 T-cell genes | `A/E33_per_comparison.csv`; `A/E33_pooled.csv` (six, five, breast) | `E33_trials.py` |
| f | breast: signatures above each gene; with NeoPembrOv | `A/E35_by_gene.csv` (`breast three`, `four chemotherapy-backbone`) | `analysis/15_v19_v41/scripts/E35_trials.py` |

### Fig. 5a–f: `analysis/17_v47_v48/figure_scripts/v48/make_fig5.py` (from `v44/make_fig5.py`)

Lab render: `figures/cns/v48/Fig5.png`. Numbers file: `analysis/17_v47_v48/figure_scripts/v48/Fig5_numbers.md`. All counts are recomputed without `TIS_Damotte` from the input files below.

| Panel | Content | Input | Written by |
|---|---|---|---|
| a | random sets reaching P < 0.05 against size; Spearman ρ computed in the script | `A/E32_primary.csv` (`ref_frac_nominal_p05`); `A/E32_coverage.csv` | `E32_trials.py` |
| b | volcano of the 2,082 genes | `A/E33_pooled.csv` | `E33_trials.py` |
| c | sensitivity counts; above CD8A | `A/E32_pooled_all_sets.csv`; `A/E32_native.csv`; `A/E33_pooled.csv`; `A/E33_paired.csv`; `A/E34_pooled.csv` | `E32_trials.py`; `analysis/15_v19_v41/scripts/E32_native.py`; `E33_trials.py`; `analysis/15_v19_v41/scripts/E34_trials.py` |
| d | breast pool against the six-comparison pool | `A/E33_pooled.csv`; `A/E33_paired.csv`; `A/E32_primary.csv`; `A/E33_breast_extras.csv` | `E33_trials.py`; `E32_trials.py`; `analysis/15_v19_v41/scripts/E33_extras.py` |
| e | before and after ESTIMATE-immune adjustment | `A/E32_primary.csv` (`est_estimate_adj`, `p_estimate_adj`) | `E32_trials.py` |
| f | inflamed signatures in 18 published randomized trials | `A/E41_inputs.csv`; `A/E41_meta.csv` | `analysis/16_v42_v46/scripts/E41_meta.R` |

### Fig. 6a–f: `analysis/17_v47_v48/figure_scripts/v48/make_fig6.py` (from `v47/make_fig6.py`)

Lab render: `figures/cns/v48/Fig6.png`. Numbers file: `analysis/17_v47_v48/figure_scripts/v48/Fig6_numbers.md`.

| Panel | Content | Input | Written by |
|---|---|---|---|
| a | prognosis, seven and ten gastric cohorts | `A/E30_dedup_runs/primary/E20_pooled.csv` | `E30_dedup_E20.R` |
| b | published gastric trials | `A/E44_A_gastric_trials.csv` | `E44_gastric.py` |
| c | single-arm cohorts, pooled, adjusted and MSS/EBV-negative | `A/E44_B_per_cohort.csv`, `A/E44_B_pooled.csv`, `A/E44_B_adjusted.csv`, `A/E46_C_pooled.csv` | `E44_gastric.py`, `E46_test.py` |
| d | responders by MSI-H or EBV-positive status | `A/E44_labels.csv` | `E44_export.R` |
| e | adjustment for MSI-H or EBV-positive status | `A/E44_B_adjusted.csv`, `A/E32_primary.csv` | `E44_gastric.py`, `E32_trials.py` |
| f | agreement with the randomized comparisons | `A/E44_B_pooled.csv`, `A/E32_primary.csv`, `A/E37_test_pooled.csv` | `E44_gastric.py`, `E32_trials.py`, `E37_test.py` |

### Supplementary Fig. S37a–f: `analysis/16_v42_v46/figure_scripts/v43/make_fig6.py` (Fig. 6 of v46; from `v40/make_fig6.py`, type one step larger)

Lab render: `figures/cns/v43/Fig6.png`. Numbers file: `analysis/16_v42_v46/figure_scripts/v43/Fig6_numbers.md`. The script imports
`barriernet.calibrate.pool` from the repository's package.

| Panel | Content | Input | Written by |
|---|---|---|---|
| a | NetImmune and CAF_ECM per comparison and pooled | per comparison: `A/E32_per_comparison.csv` (`primary`); six: `A/E32_primary.csv`; breast: `A/E33_breast_extras.csv`; renal three: pooled in the script from `E32_per_comparison.csv` | `E32_trials.py`; `E33_extras.py`; renal pool computed in the script |
| b | development: held-out log-odds of the five candidates | `A/E37_loco_pooled.csv` | `analysis/15_v19_v41/scripts/E37_develop.R` |
| c | frozen module weights | `A/E37_frozen_weights.tsv` (summed by module) | `E37_develop.R` |
| d | barriernet 0.3, NetImmune 0.2.1 and CD8A per comparison and pooled | `A/E37_test_per_comparison.csv`; `A/E37_test_pooled.csv` | `analysis/15_v19_v41/scripts/E37_test.py` |
| e | random reference (1,000 draws) | `A/E37_reference_draws.csv`; `A/E37_test_pooled.csv` | `analysis/15_v19_v41/scripts/E37_refdraws.py`; `E37_test.py` |
| f | paired difference from each panel gene (added after the test) | `A/E37_panel_extra.csv` | `analysis/15_v19_v41/scripts/E37_panel_extra.py` |

### Supplementary Fig. S38a–h: `analysis/16_v42_v46/figure_scripts/v43/make_fig7.py` (Fig. 7 of v46; from `v41/make_fig7.py`, type one step larger; panel c reads the duplicate-free pool)

Lab render: `figures/cns/v43/Fig7.png`. The script does not write a numbers file. Its seven data tables (`analysis/16_v42_v46/figure_scripts/v43/data/`)
are byte-identical to those of `v41/data/`.

| Panel | Content | Input | Written by |
|---|---|---|---|
| a | 21 checkpoint cohorts, response | `A/E12_icb_per_cohort.csv`; `A/E12_icb_pooled.csv` | `analysis/10_v16/E12_icb_expanded.R` |
| b | TCGA, 24 types, adjusted OS | `A/E11_per_type.csv`; `A/E11_pooled.csv` | `analysis/10_v16/E11_pancancer_prognosis.R` |
| c | prognosis against response | x: `A/E30_dedup_runs/primary/E20_pooled.csv` (`all 10 gastric`, unadjusted). y: `A/E12_icb_pooled.csv` (`all`) | `analysis/15_v19_v41/scripts/E30_dedup_E20.R`; `E12_icb_expanded.R` |
| d | the same within tumor type | `A/E11_dissociation_by_tumour.csv` | `E11_pancancer_prognosis.R` |
| e | leave-one-cohort-out range | `A/E12_icb_pooled.csv`; `A/E12_icb_loo.csv` | `E12_icb_expanded.R` |
| f | cohort subsets | `A/E12_icb_pooled.csv`; `A/E16_perfect_pooled.csv` | `E12_icb_expanded.R`; `analysis/11_v17/E16a_perfect.R` |
| g | size-matched null, first 11 cohorts | `A/E5_icb_null_summary.csv`; `A/E5_icb_null_draws.csv` | `analysis/09_v15/E5_icb_split.R` |
| h | random-gene null, TCGA survival | `A/E11_null.csv`; `A/E11_null_draws.csv` | `E11_pancancer_prognosis.R` |

---

## 4. Numbers in the text

Scripts are under `analysis/15_v19_v41/scripts/` unless another folder is given. Outcome files are
`analysis/15_v19_v41/protocols/E*_outcome.md`.

### Abstract

| Numbers | Source file | Script |
|---|---|---|
| six comparisons, 1,828 patients | `A/E32.log` | `E32_trials.py` |
| 51 distinct signatures; 16 favoring; 25 above size-matched random genes; none above CD8A | `A/E46_E_distinct.csv` (from `A/E32_primary.csv`) | `E46_distinct.py`; `E32_trials.py` |
| none favoring with CheckMate 025 overall survival | `A/E32_pooled_all_sets.csv` | `E32_trials.py` |
| 15 above coexpression-matched random genes | `A/E45_B_coexpression.csv`; `A/E46_E_distinct.csv` | `E45_test.py`; `E46_distinct.py` |
| T-cell mean +0.22 (0.07–0.37); detectable difference 0.20; none above it | `A/E46_A_reference_pooled.csv`; `A/E46_A_counts.csv` | `E46_test.py` |
| 7 above CXCL10 | `A/E45_A_paired.csv`; `A/E46_E_distinct.csv` | `E45_test.py`; `E46_distinct.py` |
| 18 other trials: ratio of hazard ratios 0.78 (0.67–0.92) | `A/E41_meta.csv` (`primary: survival, independent of the six`, as exp(−y)) | `E41_meta.R` |
| prognostic marker: interaction 1.1% | `A/E36_S0_S2.csv` (S1) | `E36_sim.py` |
| NetImmune HR 0.70 (0.62–0.79), seven cohorts; ESTIMATE | `A/E30_dedup_runs/primary/E20_pooled.csv` (`locked 7 (primary)`, unadjusted) | `E30_dedup_E20.R` |
| CD8A by the main rule, P = 0.22 | `A/E44_B_pooled.csv` (`cons_p`) | `E44_gastric.py` |
| 73 MSS, EBV-negative patients: no signature significant | `A/E46_C_pooled.csv`; `A/E46_C_counts.csv` | `E46_test.py` |
| each published trial under 50% power | `A/E44b_checks.csv` | `E44b_checks.py` |

### Key messages

| Numbers | Source file | Script |
|---|---|---|
| 7 of 51 above CXCL10; none above the T-cell mean or CD8A | as in the abstract | as in the abstract |
| 18 other randomized trials; about 1,200 to 2,800 events | `A/E41_meta.csv`; `A/E39_required.csv`, log `E39_design.log` | `E41_meta.R`; `E39_design.py` |

### Materials and Methods

| Numbers | Source file | Script |
|---|---|---|
| 52 signatures: 34 sensitive, 18 resistant; 20 inflamed, 8 confounder, 24 other | `A/E32_primary.csv` (`declared_dir`, `family`); `A/E32.log` | `E32_trials.py` |
| T-cell panel (13 genes) and IFN-γ/checkpoint panel (14 genes) | `T_PANEL`, `I_PANEL` in the script | `E33_trials.py` (repeated in `E34_trials.py`) |
| I-SPY2: all 69 pembrolizumab patients and 40 of 179 controls with the linear factor | `A/E32.log` (Q6) | `E32_trials.py` |
| first three comparisons, NetImmune P = 0.057 | `A/E10_meta.csv` (`primary: IMmotion150 + 2 I-SPY2`); `A/E21_meta.csv` (`reproduction: 3 published primary`) | `analysis/10_v16/E10_meta.R`; `analysis/14_v18/E21_meta.R` |
| first three locked gastric cohorts, P = 0.057 | `A/E15_pooled.csv` (`new cohorts`, NetImmune, unadjusted); the same value in `A/E30_dedup_runs/off/E20_pooled.csv` | `analysis/11_v17/E15_gastric_validation.R`; `E30_dedup_E20.R` |
| 1,000 random sets; 300 bootstrap resamples; smallest P 1/1,001 | `B_REF`, `B_BOOT` in the script | `E32_trials.py` |
| detectable difference = 2.8 × median paired SE | Q4 in the script | `E33_trials.py` |
| simulation design: 0.26 per SD; β 0.1, 0.2, 0.3; 1,000 and 200 replicates; 200 random sets | constants in the script | `E36_sim.py` |
| barriernet 0.3: 20 cohorts, 1,195 patients, 368 responders | `A/E37_dev_genes.csv` (1,198 / 369) minus the 3 samples in `A/E37_twin_flags.csv` | `E37_export.R`; `E37_twins.py`; `E37_develop.R` prints the number of samples removed |
| "within 0.02 of the best" selection rule | the selection loop in the script | `E37_develop.R` |
| 11 cohorts (796 patients, 259 responders); 10 cohorts (544, 149) | `A/E12_icb_pooled.csv` (`E5 11 (reference)`, `new only (replication)`) | `analysis/10_v16/E12_icb_expanded.R` |
| earlier unreported analysis: 107 lists, 40 of the 52, CD8A +0.22, composites +0.29 to +0.30 | declared in `analysis/15_v19_v41/protocols/E32_prereg.md` | the archived code and outputs of that analysis: `analysis/15_v19_v41/E24_archived/` (`e24_main.csv`, `e24_int_pooled.csv`) |
| seventeen further protocols | `analysis/15_v19_v41/protocols/E28_prereg.md` … `E40_prereg.md`; `analysis/16_v42_v46/protocols/E41_prereg.md` … `E44_prereg.md` (each with `.sha256`) | — |
| published trials: 629 records; 41 trials listed in advance; 1.60 for a median split | `analysis/16_v42_v46/protocols/E41_prereg.md`; `E41_outcome.md` | `E41_adjudicate.py`; `E41_meta.R` |
| barriernet 0.3.0 reproduced the analysis to within 5.2 × 10⁻⁶ | `A/E37_package_check.csv` (largest `abs_diff` 5.16 × 10⁻⁶) | `E37_package_check.py` |

### Results: Error rates and power in simulation

| Numbers | Source file | Script |
|---|---|---|
| null rejection 1.7%–2.5%; prognostic CAF_ECM 1.1% against 94.9% in the checkpoint arm alone | `A/E36_S0_S2.csv` (S0, S1) | `E36_sim.py` |
| power 17%, 59%, 86%; no bias (mean estimates 0.102, 0.200, 0.301) | `A/E36_S0_S2.csv` (S2) | `E36_sim.py` |
| random reference: 77%–99% (CD8A truth); 0.5%–5% (transcriptome-wide truth), with positive mean estimates for CD8A, T_CELL and IFN; 6%–10.5% (null); CAF_ECM 20.5% and 29.5% | `A/E36_S3.csv` (`exceed_rate`, `mean_est`) | `E36_sim.py` |
| paired test: a module above CD8A in at most 0.5%, below in 17%–48% | `A/E36_S4.csv` (`S4 truth CD8A`) | `E36_sim.py` |
| power matched normal theory; SE 0.029–0.062 | `A/E36_S4.csv` (`median_se`, `sig_above`); prediction computed in `figure_scripts/v40/make_fig2.py` (Fig2_numbers.md f) | `E36_sim.py` |
| three of six expectations held | `A/E36_expectations.csv`; `E36_outcome.md` | `E36_sim.py` |

### Results: Six randomized comparisons with expression from both arms

| Numbers | Source file | Script |
|---|---|---|
| 1,828 patients: 907 checkpoint, 921 control | `A/E32.log` (per-comparison counts); typed in Fig. 3a and Table 1 | `E32_trials.py` |
| all 52 scored in at least five comparisons; 25 below the coverage rule in IMmotion150 | `A/E32_primary.csv` (`k`); `A/E32_coverage.csv`; `A/E33_trial_table.csv` (IMM150: 27 scored) | `E32_trials.py`; `E33_extras.py` |
| IMmotion150 covers 2,475 genes | `A/E32.log` | `E32_trials.py` |
| the full T-cell panel and 2,082 genes measured in all six | `A/E33_pooled.csv` (`all_six`); `A/E33_trials.log` | `E33_trials.py` |

### Results: Published signatures in the randomized comparisons and against random genes

| Numbers | Source file | Script |
|---|---|---|
| 17 of 52 favoring at P < 0.05; IMPRES −0.23 (−0.41 to −0.06; P = 0.018); 19 of 20 inflamed toward benefit; significant estimates +0.15 to +0.24 | `A/E32_primary.csv`; `A/E32_expectations.csv` (Q1a, Q1b) | `E32_trials.py` |
| carried chiefly by IMmotion150 | `A/E32_per_comparison.csv` (Fig. 3f) | `E32_trials.py` |
| CYT and CD39–CD8 T cell labeled as resistance | `A/E32_primary.csv` (`declared_dir`) | `E32_trials.py` |
| audit: 3 contradicted, 12 unclear | `analysis/15_v19_v41/protocols/E33_direction_audit.csv` / `.md` | written by an AI assistant, not by a script |
| 26 above their reference: 15 of 20 inflamed, 11 of 24 other, 0 of 8 confounder | `A/E32_primary.csv` (`emp_p_declared`, `family`) | `E32_trials.py` |
| Fig. 2a sentence: T-cell-inflamed profile +0.20 (0.03 to 0.38; P = 0.029); empirical P = 0.002; +0.002 from CD8A (P = 0.96) | `A/E32_primary.csv`; `A/E33_pooled_reference_estimates.npz`; `A/E33_paired.csv` (analysis/16_v42_v46/figure_scripts/v43/Fig2_numbers.md) | `E32_trials.py`; `E33_trials.py` |
| 28 with the audited directions | `A/E33_expectations.csv` (Q7; recount from `A/E33_pooled_reference_estimates.npz`) | `E33_trials.py` |
| median of the reference medians −0.06; 11 of 26 not nominally significant; ImmuneCells negative | `A/E32_primary.csv` (`ref_median`, `p`, `est`) | `E32_trials.py` |

### Results: Comparison with single genes

| Numbers | Source file | Script |
|---|---|---|
| CD8A +0.23 (0.05 to 0.40; P = 0.022; 5 of 6) | `A/E32_primary.csv` | `E32_trials.py` |
| CD8A per trial: −0.17 (I-SPY2 durvalumab–olaparib) to +0.53 (IMmotion150) | `A/E33_trial_table.csv` (= `E32_per_comparison.csv`) | `E33_extras.py` |
| none above CD8A; 8 positive but not significant; 15 below, 7 of them resistance-labeled | `A/E32_primary.csv` (`diff_vs_CD8A`, `cd8a_p`) | `E32_trials.py` |
| none above CD8A over the five comparisons without IMmotion150 | `A/E33_paired.csv` (`five without IMmotion150`); `A/E33_expectations.csv` (Q5) | `E33_trials.py` |
| NetImmune +0.32, TMEscore +0.31, ESTIMATE +0.30; paired +0.03 (P = 0.69), +0.06 (0.44), +0.11 (0.098) | `A/E32_primary.csv` (Table 2) | `E32_trials.py` |
| CD8A at the 99th percentile, T-cell median at the 94th; CD8A second, behind PRF1 (+0.229 against +0.226); 7 of 13 at P < 0.05 | `A/E33_pooled.csv`; `A/E33_expectations.csv` (Q1, Q2) | `E33_trials.py` |
| CD8A's rank: first in CheckMate 025 to 12th in I-SPY2 pembrolizumab | `A/E33_per_comparison.csv` (ranked in `figure_scripts/v39/make_fig3.py`) | `E33_trials.py` |
| CXCL13 and STAT1 +0.28; CXCL9 +0.17 (−0.11 to 0.45), measured in five comparisons | `A/E33_pooled.csv` | `E33_trials.py` |
| for 11 of 13 T-cell genes at most two signatures above, and for 8 none; CD3G 3, CD247 8, CXCL10 8, HAVCR2 17 | `A/E33_Q3_by_reference.csv`; `A/E33_expectations.csv` (Q3) | `E33_trials.py` |
| 101 genes favor the checkpoint arm and 102 the control arm | `A/E33_pooled.csv` | `E33_trials.py` |
| median paired SE 0.076; detectable difference 0.21; upper limits up to +0.30 (median +0.05) | `A/E33_expectations.csv` (Q4) | `E33_trials.py` |

### Results: Sensitivity and subgroup analyses

| Numbers | Source file | Script |
|---|---|---|
| random sets at P < 0.05 in 0.4%–21% of draws (median 4.8%); ρ = 0.39 | `A/E32_primary.csv` (`ref_frac_nominal_p05`); ρ computed in `figure_scripts/v39/make_fig4.py` with `A/E32_coverage.csv` | `E32_trials.py` |
| per trial: 3 of 52 (NeoTRIP) to 28 of 52 (JAVELIN Renal 101) above random sets | computed in `figure_scripts/v34/make_fig3.py` from `A/E32_per_comparison.csv` and `A/E32_reference_draws.npz`; listed in `figure_scripts/v34/Fig3_numbers.md`. The current Fig. 3f script does not compute it | `E32_trials.py` |
| nominal test calibrated: 2.4%–6.8% (NetImmune, ESTIMATE, CD8A) | `A/E28_Q3_calibration_per_comparison.csv` | `E28_trials.py` |
| without IMmotion150: 9 signatures; CD8A +0.16 (0.007 to 0.32) | `A/E32_pooled_all_sets.csv` (`five without IMmotion150`) | `E32_trials.py` |
| CheckMate 025 overall survival: 1 against; CD8A +0.18 (−0.04 to 0.41; P = 0.088) | `A/E32_pooled_all_sets.csv` | `E32_trials.py` |
| + PRINCE and DUTRENEO 20 / 2; angiogenesis 18 / 1; hormone receptor only 14 / 1 | `A/E32_pooled_all_sets.csv` | `E32_trials.py` |
| published engines 23 / 1 | `A/E32_native.csv` | `E32_native.py` (+ `E32_native_scores.R`) |
| breast: IFN (Ayers) +0.47 (0.05 to 0.90); 2 against; 37 of 52 toward benefit | `A/E33_pooled.csv` (`breast three`) | `E33_trials.py` |
| breast: CD8A +0.06 (−0.33 to 0.45) | `A/E33_breast_extras.csv`; `A/E33_pooled.csv` | `E33_extras.py`; `E33_trials.py` |
| breast: six signatures above CD8A | `A/E33_paired.csv` (`breast three`); `A/E35_by_gene.csv` | `E33_trials.py`; `E35_trials.py` |
| "about one would be expected by chance" | 52 × 0.025 ≈ 1.3, stated in `E33_outcome.md` | not written by a script |
| breast: GZMB, IDO1 and CXCL10 at P < 0.05 | `A/E33_pooled.csv` (`breast three`) | `E33_trials.py` |
| none above GZMB, IDO1, CXCL10 or the 13 strongest breast genes; 3 to 15 signatures above each of nine T-cell genes | `A/E35_by_gene.csv`; `A/E35_expectations.csv` (Q1, Q2) | `E35_trials.py` |
| breast comparisons can detect only 0.6–0.7 per SD | 2.8 × the median SE by reference gene in `A/E35_paired.csv` (`breast three`; 0.54–0.76, GZMB 0.59, IDO1 0.73, CXCL10 0.71) | re-derived by `scripts/text_numbers_check.py` → `A/text_numbers_check.csv` |
| NeoPembrOv added: seven 18 / 1, none above CD8A | `A/E34_pooled.csv`; `A/E34_expectations.csv` | `E34_trials.py` |
| NeoPembrOv added: four-comparison pool CD8A +0.06 | `A/E34_pooled.csv` | `E34_trials.py` |
| NeoPembrOv added: none above GZMB, IDO1 or CXCL10 | `A/E35_expectations.csv` (Q4) | `E35_trials.py` |
| ESTIMATE-immune adjustment: the 18 grew; 2 kept P < 0.05 | `A/E32_primary.csv` (`est_estimate_adj`, `p_estimate_adj`); `A/E32_expectations.csv` (Q4) | `E32_trials.py` |
| MammaPrint: 15 → 18 significant; 40 controls enlarged the interactions | `A/E32_pooled_all_sets.csv` (`six, I-SPY2 hormone receptor only` against `six (primary)`); `A/E32_ispy2.csv` | `E32_trials.py` |

### Results: Immune and fibroblast components of barriernet

| Numbers | Source file | Script |
|---|---|---|
| NetImmune +0.32 (0.03 to 0.60; P = 0.037), 6 of 6; Immune +0.20 (0.05 to 0.35; P = 0.019); CAF_ECM −0.25 (−0.60 to 0.10; P = 0.12), 0 of 6 | `A/E32_primary.csv` | `E32_trials.py` |
| breast: NetImmune +0.71 (0.19 to 1.23; P = 0.027), CAF_ECM −0.71 (−1.25 to −0.18; P = 0.029) | `A/E33_breast_extras.csv` | `E33_extras.py` |
| CAF_ECM's renal interactions vanished after angiogenesis adjustment | `A/E30_trials_angio.csv` | `E30_trials.py` → `E30_trials_meta.R` |
| immune against reversed fibroblast in one model, P = 0.47 | `A/E28_Q2_Q6_pooled.csv` (Q2a, all six) | `E28_trials.py` |
| first three P = 0.057; three added P = 0.22 | `A/E21_meta.csv` (`reproduction: 3 published primary`, `3 new primary only`) | `analysis/14_v18/E21_meta.R` |

### Results: Development and test of barriernet 0.3

| Numbers | Source file | Script |
|---|---|---|
| held-out log-odds: learned candidates +0.49 to +0.50; NetImmune +0.43; Immune +0.46; module weights chosen | `A/E37_loco_pooled.csv`; `A/E37_dev_expectations.csv` | `E37_develop.R` |
| weights: IFN +0.36, APM +0.23, CAF_ECM −0.16 | `A/E37_frozen_weights.tsv` (summed by module) | `E37_develop.R` |
| +0.25 (0.01 to 0.49; P = 0.042), 6 of 6 and NeoPembrOv | `A/E37_test_pooled.csv`; `A/E37_test_per_comparison.csv` | `E37_test.py` |
| above random genes, empirical P = 0.001 | `A/E37_test_pooled.csv` (`emp_p`); draws in `A/E37_reference_draws.csv` | `E37_test.py`; `E37_refdraws.py` |
| against CD8A −0.005 (−0.09 to 0.08; P = 0.92); against NetImmune −0.02 (−0.10 to 0.07; P = 0.71) | `A/E37_test_pooled.csv` | `E37_test.py` |
| all four expectations held | `A/E37_test_expectations.csv`; `E37_outcome.md` | `E37_test.py` |
| breast: above CD8A by +0.43 (0.08 to 0.79; P = 0.017); four comparisons: by +0.36 (P = 0.019) | `A/E37_test_pooled.csv` | `E37_test.py` |
| whole panel: above 3 of 27 and below none (six); above 15 of 26 and below none (breast) | `A/E37_panel_extra.csv` | `E37_panel_extra.py` |
| GZMB, IDO1, CXCL10, STAT1: +0.05 to +0.15, P ≥ 0.28 | `A/E37_panel_extra.csv` | `E37_panel_extra.py` |
| against the 52: above 6, below none; CD8A above 4 | `A/E42_counts.csv` (`six (primary)`, `conservative`) | `E42_signatures.py` |
| interval within ±0.10 for CD8A and for 11 of the 15 signatures with the largest estimates | `A/E42_paired.csv`; `A/E42_counts.csv`; `E42_outcome.md` | `E42_signatures.py` |
| r = 0.81 with CD8A | `E42_outcome.md`; see section 7 | not written by a script |
| both expectations of this comparison held | `A/E42_expectations.csv` | `E42_signatures.py` |
| C-for-benefit, pre-specified: −0.012 (−0.048 to 0.025) | `A/E43_contrasts.csv` (`M_2AX - M_CD8A`) | `E43_benefit.py` |
| post hoc correction: +0.033 (−0.004 to 0.069) in C-for-benefit, −0.023 (−0.115 to 0.068) in gradient; fibroblast module added to CD8A | `A/E43b_contrasts.csv` (`M_2AX - M_CD8A`, `M_CD8A_CAF - M_CD8A`) | `E43b_null.py` |
| audit of the cross-validated estimate | `E43_outcome.md`; protocol `E43b_prereg.md` | `E43b_null.py` |
| three of the four expectations held | `A/E43_expectations.csv` | `E43_benefit.py` |

### Results: Prognostic and response associations of the barriernet components

| Numbers | Source file | Script |
|---|---|---|
| benchmark cohorts: CAF_ECM HR 1.29 (1.18–1.41; P = 0.0011), the largest of the six modules | `A/E20_pooled.csv` (`benchmark 3 (reference)`); identical in `A/E30_dedup_runs/primary/` | `analysis/14_v18/E20_extended_replication.R` |
| ten gastric cohorts: 1.31 (1.17–1.46) | `A/E30_dedup_runs/primary/E20_pooled.csv` (`all 10 gastric`, unadjusted) | `E30_dedup_E20.R` (`E30_MODE=primary`) |
| direct test against the Immune axis, ten cohorts: P = 0.22 | `A/E30_dedup_runs/primary/E29_pooled.csv` (`Q2 survival, 10 gastric`, Immune, method a); first run `A/E29_pooled.csv` (P = 0.2245) | `E30_dedup_E29.R`; `E29_dissociation.R` |
| 21 checkpoint cohorts: 1,340 patients, 408 responders; Immune +0.43 (0.24–0.61), positive in 19; CAF_ECM −0.06 (−0.21 to 0.09) | `A/E12_icb_pooled.csv` (`all`; wider of the DL and HK intervals) | `analysis/10_v16/E12_icb_expanded.R` |
| direct test, response: +0.26 (0.08 to 0.45; P = 0.0083) | `A/E29_pooled.csv` (`Q1 response, 21 cohorts`, Immune, method a) | `E29_dissociation.R` |
| TCGA, five checkpoint cancer types: CAF_ECM adjusted HR 1.07 (0.99–1.15) | `A/E11_pooled.csv` (`(i) SKCM KIRC LUAD LUSC BLCA`, OS, adjusted) | `analysis/10_v16/E11_pancancer_prognosis.R` |
| Immune axis not consistently prognostic; responders marked except in kidney cancer | `A/E11_dissociation_by_tumour.csv` (Fig. 7d) | `E11_pancancer_prognosis.R` |
| NetImmune in seven locked cohorts: HR 0.70 (0.62–0.79) | `A/E30_dedup_runs/primary/E20_pooled.csv` (`locked 7 (primary)`) | `E30_dedup_E20.R` |

### Results: Evidence in the source papers and in single-arm cohorts

See section 8 (E38 and E40).

### Results: Inflamed signatures in other randomized trials

| Numbers | Source file | Script |
|---|---|---|
| 18 comparisons in ten cancer types; overall survival 11, progression-free survival 7 | `A/E41_inputs.csv` | `E41_meta.R` |
| three more eligible trials not extractable; one formal interaction test (CheckMate 649) | `A/E41_adjudicated.csv`; `E41_outcome.md` | `E41_adjudicate.py` |
| ratio of hazard ratios 0.78 (0.67 to 0.92; P = 0.004; I² = 7%); about +0.15 per SD (0.06 to 0.24) | `A/E41_meta.csv` (`primary: survival, independent of the six`; `primary, per SD`) | `E41_meta.R` |
| leaving out one trial: 0.76 to 0.82 | `A/E41_meta.csv` (`leave one out`) | `E41_meta.R` |
| five trials with values in text or tables: 0.77 (0.51 to 1.17) | `A/E41_meta.csv` (`survival, excluding numbers read from figures`) | `E41_meta.R` |
| PD-L1 0.85 (0.65 to 1.10), signatures 0.79 (0.66 to 0.94), ten trials; tumor mutational burden 0.93 (0.64 to 1.36), seven | `A/E41_meta.csv` | `E41_meta.R` |
| all five expectations held | `A/E41_expectations.csv` | `E41_meta.R` |
| inflamed signatures in the six comparisons: median +0.19 | `A/E32_primary.csv` (`family`, median of `est` = 0.187); `E41_outcome.md`; see section 7 | not written by a script |
| POPLAR single genes: hazard ratios 0.39 to 0.46 and 0.97 to 1.13, against 0.43 and 1.10 | `E41_outcome.md`, read from the published figure; outside the protocol | none |

### Results: Gastric and gastroesophageal cancer

| Numbers | Source file | Script |
|---|---|---|
| CheckMate 649 +0.17 (−0.07 to 0.40); KEYNOTE-061 −0.04 (−0.37 to 0.30); detectable 0.33 and 0.48 | `A/E44_A_gastric_trials.csv` | `E44_gastric.py` |
| CD8A in the six comparisons +0.23 | `A/E32_primary.csv` | `E32_trials.py` |
| 119 patients in three cohorts | `A/E44.log` (per-cohort counts) | `E44_gastric.py` |
| CD8A +0.58 (0.12 to 1.04); barriernet 0.3 +0.67 (0.21 to 1.14), fixed effect | `A/E44_B_pooled.csv` | `E44_gastric.py` |
| barriernet 0.3 minus CD8A −0.07 (P = 0.84); none above CD8A, 10 below | `A/E44_B_paired.csv` | `E44_gastric.py` |
| ρ = 0.76; all four expectations held | `A/E44_expectations.csv`; `A/E44b_checks.csv` | `E44_gastric.py`; `E44b_checks.py` |
| main-text rule: CD8A +0.65 (−0.94 to 2.24; P = 0.22), barriernet 0.3 P = 0.18, CXCL10 P = 0.12; 2 of the 51 signatures, and 20 by fixed effect | `A/E44_B_pooled.csv` (`cons_est`, `cons_p`); `A/E45_A_gastric_pooled.csv`; `A/E46_E_distinct.csv` | `E44_gastric.py`; `E45_test.py`; `E46_distinct.py` |
| 13 of 14 MSI-H or EBV-positive patients responded, 28 of 73 others | `A/E44_labels.csv` (patient-level, not in the repository; recounted); see section 7 | `E44_export.R` |
| adjusted: CD8A +0.81 to +0.44 (P = 0.19); signatures at P < 0.05 from 16 to 1; 45 of 51 positive | `A/E44_B_adjusted.csv`; `A/E46_E_distinct.csv` | `E44_gastric.py`; `E46_distinct.py` |
| CXCL10 +0.89 (0.39 to 1.39); minus CD8A +0.33 (P = 0.29); adjusted +0.88 to +0.62 (P = 0.058) | `A/E45_A_gastric_pooled.csv` | `E45_test.py` |
| 73 other patients: CD8A +0.40 (P = 0.23), T-cell mean +0.42, barriernet 0.3 +0.39, Immune axis +0.40; 45 of 51 positive, none at P < 0.05; 28 responders, 26 in one cohort | `A/E46_C_pooled.csv`; `A/E46_C_counts.csv`; `A/E46_C_per_cohort.csv` | `E46_test.py` |
| CD8A prognostic in ten gastric survival cohorts: HR 0.90 (0.84–0.96) | `A/E44_C_prognosis.csv` | `E44_gastric.py` |
| seven cohorts: NetImmune 0.70 (0.62–0.79), fibroblast module 1.34 (1.09–1.65), CD8A 0.88 (0.80–0.97), ESTIMATE 0.73 (0.62–0.85); ten cohorts: 0.72, 1.31, 0.90 | `A/E30_dedup_runs/primary/E20_pooled.csv` | `E30_dedup_E20.R` |
| TMEscore minus CD8A +0.06 (P = 0.44) | `A/E44_D_tmescore.csv` | `E44_gastric.py` |

### Discussion

| Numbers | Source file | Script |
|---|---|---|
| random genes significant in up to a fifth of draws | as in "Sensitivity and subgroup analyses" (21%) | `E32_trials.py` |
| IMmotion150: 2,475 genes, 25 signatures unscored | `A/E32.log`; `A/E33_trial_table.csv` | `E32_trials.py`; `E33_extras.py` |
| 29 barriernet-validation expectations, 21 held (Supplementary Fig. S31b) | rows drawn by `figure_scripts/make_fig8.py` from `A/E17_expectations.csv`, `E18_cps.csv`, `E19_expectations.csv`, `E20_expectations.csv`, `E20_pooled.csv`, `E21_expectations_pooled.csv`, `E21_expectations_per_trial.csv`, `E22_expectations.csv` | the 21 / 29 tally is not written by any script: not established |
| First ten later protocols: 43 expectations, 33 held (later tally below) | E28–E31 6 of 11 and E32 6 of 8 (`E32_outcome.md`); E33 6 of 6, E34 3 of 3, E35 3 of 3, E36 3 of 6 (`A/E33…E36_expectations.csv`); E37 2 of 2 + 4 of 4 (`A/E37_dev_expectations.csv`, `A/E37_test_expectations.csv`) | the total is not written by any script |
| IMmotion151 823 tumors; OAK and POPLAR 891 | cited publications (refs 101, 102) | — |
| gastric trials: each below 50% power for an interaction of CD8A's size (47%, 26%) | `A/E44b_checks.csv` | `E44b_checks.py` |
| Limitations: 84 expectations in the nineteen later protocols, 68 held | `analysis/17_v47_v48/protocols/later_protocol_tally.md` | — |
| Limitations: without PRJEB25780 neither CD8A nor barriernet 0.3 was significant | `A/E44b_checks.csv` (`leave one out`) | `E44b_checks.py` |

---

## 5. Supplementary figures

Each package file `final_v48/figures/supplementary/FigureS<n>.png` was matched by MD5 to its lab render (S1–S33 are also byte-identical to `final_v41/`). Figure
scripts are under `analysis/15_v19_v41/figure_scripts/`. The supplementary scripts take a key (`python
make_supp_<x>.py <key>`), write `figures/cns/supp/FigS_<key>.png`, and write a numbers file
`FigS_<key>_numbers.md` beside it.

Short names:

- `R1` = `gastric/revision/results`;
- `R2` = `gastric/revision2/results`;
- `GC` = `gastric/GC_TME_Reprogrammer`;
- `FD` = `gastric/2/figures/data`.

| Fig. | Lab render | Script (key) | Main inputs → written by |
|---|---|---|---|
| S1 | `supp/FigS_sc_atlas` | `supp/make_supp_sc.py` (`sc_atlas`) | `FD/sc/cells_obs.parquet`, `cells_expr.parquet` → `analysis/07_figures/sc_extract.py` (the GSE183904 object it reads: not established; see `README.md`) |
| S2 | `supp/FigS_sc_programmes` | `supp/make_supp_sc.py` (`sc_programmes`) | `FD/sc/cells_obs.parquet` → `sc_extract.py` |
| S3 | `supp/FigS_sc_composition` | `supp/make_supp_sc.py` (`sc_composition`) | `GC/out_scanpy_strict_GSE183904/composition_by_sample_anno_level1.csv` → not established |
| S4 | `supp/FigS_sc_fractions` | `supp/make_supp_sc.py` (`sc_fractions`) | `GC/results/Fig5_Fibro/…`, `Fig5_Macro/…`, `Fig5_FibroMacro/…` → not established |
| S5 | `supp/FigS_sc_indep` | `supp/make_supp_sc.py` (`sc_indep`) | `A/E17_attribution_by_class.csv`, `E17_share_by_patient*.csv`, `E17_gene_expression_by_class.csv` → `analysis/14_v18/E17_singlecell_attribution.py` |
| S6 | `supp/FigS_meta_forest` | `supp/make_supp_gastric.py` (`meta_forest`) | `R1/B4_meta.csv`, `B4_univariate.csv` → `analysis/01_bulk_prognosis/B4_survival.R`; `A/E8_meta_DL_HK.csv` → `analysis/09_v15/E8_prognostic_hardening.R` |
| S7 | `supp/FigS_robust` | `supp/make_supp_gastric.py` (`robust`) | `A/E19_runs.csv`, `E19_leave_one_gene_out.csv`, `E19_panel_subset.csv` → `analysis/14_v18/E19_robustness.R` |
| S8 | `supp/FigS_fitted_models` | `supp/make_supp_gastric.py` (`fitted_models`) | `A/E6_model_classes.csv` → `analysis/09_v15/E6_models.R`; `E6_review_dC.csv`, `E6_within_cohort_cv.csv` → `analysis/09_v15/E6_review.R` |
| S9 | `supp/FigS_upgrades` | `supp/make_supp_gastric.py` (`upgrades`) | `R2/S3_state_loadings.csv`, `S5_loco_cindex.csv`, `S8_neigh_pred.csv`, `S10_prognosis.csv` → `analysis/03_ceiling_upgrades/S3_state_axes.py`, `S5_loco.R`, `S8_learn_arrangement.py`, `S10_decisive.R` |
| S10 | `v32_dedup/FigS_locked` | `supp/make_supp_gastric.py locked`, run with `v32_dedup/cns_style_overlay.py` (results overlay; `v32_dedup/README.md`) | the E15/E20 files of `A/E30_dedup_runs/primary/` → `E30_dedup_E15.R`, `E30_dedup_E20.R`; `A/E8_univariate.csv` → `E8_prognostic_hardening.R` |
| S11 | `supp/FigS_prjeb_detail` | `supp/make_supp_gastric.py` (`prjeb_detail`) | `A/E3_ici_ceiling.csv`, `E3_ici_matched_null.csv` → `analysis/05_v14/E3_ici_null_and_ceiling.py`; `E7_scores_log2.csv` → `analysis/09_v15/E7a_scores_log2.py`; `E7_auc_ci_log2.csv`, `E7_delong_log2.csv` → `E7b_roc_delong.R`; `E7_samepipeline_cv.csv` → `E7d_subgroup_samepipe.py` |
| S12 | `supp/FigS_gastric_icb` | `supp/make_supp_gastric.py` (`gastric_icb`) | `A/E22_PRJEB60680_response.csv`, `E22_pooled.csv` → `analysis/14_v18/E22_gastric_icb.R`; `E22_firth_cps.csv` → `analysis/14_v18/E22_firth_cps.py` |
| S13 | `supp/FigS_icb21_forest` | `supp/make_supp_pancancer.py` (`icb21_forest`) | `A/E12_icb_per_cohort.csv`, `E12_icb_pooled.csv` → `analysis/10_v16/E12_icb_expanded.R` |
| S14 | `supp/FigS_tcga` | `supp/make_supp_pancancer.py` (`tcga`) | `A/E11_per_type.csv`, `E11_pooled.csv`, `E11_dissociation_scores.csv` → `analysis/10_v16/E11_pancancer_prognosis.R` |
| S15 | `supp/FigS_perfect_metabric` | `supp/make_supp_pancancer.py` (`perfect_metabric`) | `A/E16_perfect_per_score.csv`, `E16_perfect_pooled.csv` → `analysis/11_v17/E16a_perfect.R`; `E16_metabric_cox.csv`, `E16_metabric_vs_tcga_brca.csv` → `analysis/11_v17/E16b_metabric.R` |
| S16 | `supp/FigS_icb11` | `supp/make_supp_pancancer.py` (`icb11`) | `A/E5_trials.csv` → `analysis/09_v15/E5_trials.py`; `E5_icb_pooled.csv` → `E5_icb_split.R`; `E5_gastric_arms_pooled.csv` → `E5_gastric_arms.R`; `E5_prjeb_arms.csv` → `E5_prjeb_arms.py` |
| S17 | `supp/FigS_trials_detail` | `supp/make_supp_pancancer.py` (`trials_detail`) | `A/E10_ispy2_interactions.csv`, `E10_ispy2_within_arm.csv`, `E10_ispy2_null_draws.csv`, `E10_ispy2_null_summary.csv` → `analysis/10_v16/E10_ispy2_trials.py` |
| S18 | `supp/FigS_trials_forest` | `supp/make_supp_pancancer.py` (`trials_forest`) | `A/E21_meta_inputs.csv`, `E21_meta.csv` → `analysis/14_v18/E21_meta.R` (inputs from `E21_trials.py`, `E5_trials.py`, `E10_ispy2_trials.py`) |
| S19 | `supp/FigS_visium` | `supp/make_supp_spatial.py` (`visium`) | `A/E27_rings_sections.csv`, `E27_spots_GSE251950.parquet` → `analysis/15_v19_v41/scripts/E27_visium_marker_tls.py` |
| S20 | `supp/FigS_visium_sens` | `supp/make_supp_spatial.py` (`visium_sens`) | `A/E27_rings_sections.csv`, `E27_rings_summary.csv` → `E27_visium_marker_tls.py` |
| S21 | `supp/FigS_barrierindex` | `supp/make_supp_spatial.py` (`barrierindex`) | `A/E27_rings_sections.csv` → `E27_visium_marker_tls.py`; the Spearman P is computed in the figure script |
| S22 | `supp/FigS_moran` | `supp/make_supp_spatial.py` (`moran`) | `GC/validation/bench/spatial_benchmark_results.csv` → `analysis/00_upstream/spatial_benchmark.py` |
| S23 | `supp/FigS_transport` | `supp/make_supp_spatial.py` (`transport`) | `R2/T4_sensitivity.csv`, `T6_marker_based.csv` → `analysis/04_transport/T4_sensitivity.py`, `T6_marker_based.py` |
| S24 | `supp/FigS_imc` | `supp/make_supp_spatial.py` (`imc`) | `A/E13_models.csv`, `E13_survival.csv` → `analysis/10_v16/E13_singlecell_barrier.py` |
| S25 | `supp/FigS_proteomics` | `supp/make_supp_pancancer.py` (`proteomics`) | `figures/Fig7_module_coherence.csv` → `analysis/07_figures/prot_stats.py`; `GC/expr.tsv` (upstream) |
| S26 | `supp/FigS_drugs` | `supp/make_supp_pancancer.py` (`drugs`) | `figures/Fig8_gdsc_points.csv` → `analysis/07_figures/gdsc_prep.py`. Not established: `GC/results/ai_recommendations_drug_top15_clean.csv`, `R1/C1_gdsc_fdr.csv` (as in `README.md`) and `gastric/validation/drug/drug_mechanism_table.csv` |
| S27 | `Fig2` (v31 Fig. 2) | `make_fig2.py` | `A/E17_attribution_by_class.csv`, `E17_gene_expression_by_class.csv` → `E17_singlecell_attribution.py`; `A/E17b_bcell_reassignment.csv` → `analysis/15_v19_v41/scripts/E17b_bcell_reassignment.py`; `FD/sc/*` → `sc_extract.py` |
| S28 | `Fig7` (v31 Fig. 7) | `make_fig7.py` | `A/E27_rings_sections.csv`, `E27_zone_means.csv`, `E27_rings_summary.csv` → `E27_visium_marker_tls.py`; `A/E13_*.csv` → `E13_singlecell_barrier.py`; `R2/T3_transport.csv`, `T6_marker_based.csv` → `analysis/04_transport/T3_transport.py`, `T6_marker_based.py`; `GC/results/ring_barrier_state_v2/BarrierIndex_by_slice_with_clinical.csv` → not established |
| S29 | `Fig4` (v31 Fig. 4) | `make_fig4.py` | `A/E7_scores_log2.csv`, `E7_auc_ci_log2.csv`, `E7_subgroup_auc.csv` → E7a/E7b/E7d; `A/E18_cps*.csv` → `analysis/14_v18/E18_cps.py`; `A/E22_*.csv` → `E22_gastric_icb.R`, `E22_firth_cps.py`; `A/E23_pooled_with_gastric.csv`, `E23_power_consistency.csv` → `analysis/15_v19_v41/scripts/E23_gastric_consistency.R`; `A/E3_ici_matched_null.csv`, `E3_ici_null_draws_log2.csv` → `E3_ici_null_and_ceiling.py` |
| S30 | `v34/Fig2` | `v34/make_fig2.py` | `A/E32_primary.csv`, `E32_per_comparison.csv` → `E32_trials.py` |
| S31 | `v34/Fig7` | `make_fig8.py`, run with `v32_dedup/cns_style_overlay.py` (the dedup overlay of `A/E30_dedup_runs/primary/`) and with the panel b title changed to "Pre-specified expectations" (`v34/README.md`). That edited copy is not saved, and the saved script still prints "Pre-registered expectations": the exact script is not established | (a) `A/E20_pooled.csv` (overlay), `E11_pooled.csv`, `E12_icb_pooled.csv`, `E22_meta_inputs.csv`, `E21_meta.csv`; (b) `A/E17_expectations.csv`, `E18_cps.csv`, `E19_expectations.csv`, `E20_expectations.csv`, `E21_expectations_*.csv`, `E22_expectations.csv` → their E11–E22 scripts |
| S32 | `Fig6` (v31 Fig. 6) | `make_fig6.py` | (a) n from `analysis/14_v18/make_v18_tables.py`; (b, f) `A/E21_meta_inputs.csv`, `E21_meta.csv` → `E21_meta.R`; (c) `A/E5_trials.csv`, `E10_ispy2_within_arm.csv`, `E21_within_arm.csv`; (d, e) `A/E21_scores.csv` → `E21_trials.py`, with I-SPY2 rescored in the script; (g) `A/E21_null_summary.csv` → `E21_trials.py` |
| S33 | `v34/Fig6` | `v33_fig2/make_fig3_v33.py` (its diff from `make_fig3.py` is in `v33_fig2/make_fig3_v33.diff`), run on the dedup overlay, with panel e extended by `A/E32_locked_bionull.csv` / `_draws.rds` and the label "pre-specified extension" (`v34/README.md`). That edited copy is not saved, and the saved script reads neither `E32_locked_bionull` file: the exact script is not established | sources listed in `v34/Fig6_numbers.md`: (a–c) `FD/km_*.csv` → `analysis/07_figures/export_fig4_data.R`, `A/E8_univariate.csv` → E8; (d, j) E15/E20 files of the overlay → `E30_dedup_E15.R`, `E30_dedup_E20.R`, and `A/E30_dedup_runs/off/E20_pooled.csv`; (e) `A/E2_null_joint.csv` → E2, `E8_biology_null.csv` → `E8b_biology_null.R`, overlay `E20_null*.csv`, `E32_locked_bionull.*` → `E32_locked_bionull.R`, `FD/null_draws.csv` → `export_fig4_data.R`, `FD/bio_null_draws.csv` → not established; (f) `A/E1_paired_dC.csv`, `E8_h2h_TMEscore.csv`, overlay `E15_deltaC.csv`, `E20_deltaC.csv`; (g) `A/E19_*summary.csv`; (h) `A/E6_model_classes.csv`, `E6_learning_*.csv`; (i) `R1/B4_cindex.csv`, `A/E8_lrt.csv` |
| S34 | `v44/FigS34` | `analysis/16_v42_v46/figure_scripts/v44/make_figS34.py` | `A/E32_ispy2.csv` → `E32_trials.py` (Fig. 5f until v43) |
| S35 | `v45/FigS35` | `analysis/16_v42_v46/figure_scripts/v45/make_figS35.py` | (a) `A/E42_paired.csv` → `E42_signatures.py`; (b, c) `A/E43b_contrasts.csv` → `E43b_null.py`, `A/E43_contrasts.csv` → `E43_benefit.py` |
| S36 | `v48/FigS36` | `analysis/17_v47_v48/figure_scripts/v48/make_figS36.py` | (a) `A/E45_single_arm_gene_pooled.csv`, `E33_pooled.csv`; (b) `A/E45_A_paired.csv`; (c) `A/E45_B_coexpression.csv`; (d) `A/E45_C_within_arm.csv`, `E45_C_pooled.csv` → `E45_select.py`, `E45_test.py`; (e) `A/E46_A_steadiness.csv`; (f) `A/E46_B_pooled.csv`, `E46_B_held_out.csv`, `E46_A_reference_pooled.csv`, `E32_primary.csv` → `E46_test.py` |
| S37 | `v43/Fig6` | `analysis/16_v42_v46/figure_scripts/v43/make_fig6.py` | Fig. 6 of v46; see section 3 |
| S38 | `v43/Fig7` | `analysis/16_v42_v46/figure_scripts/v43/make_fig7.py` | Fig. 7 of v46; see section 3 |

---

## 6. Workbook `BarrierNet_Supplementary_Tables.xlsx`

The workbook has 197 sheets. Its README sheet names a source script for each sheet; this table restates them
with repository paths.

- Sheets S1–S25 are copied from the v39 workbook.
- Sheets S27a–j (`manuscript/make_v41_workbook.py`), S28a–i (`make_v44_workbook.py`), S29a–c and S30a–i (`make_v45_workbook.py`) and S31a–i (`make_v46_workbook.py`) are copied from the files in sections 8 and 9, and S32a–l (`make_v47_workbook.py`) and S32m–w (`make_v48_workbook.py`) from the files in section 10.
- Sheets S26a–l are added by `manuscript/make_v40_workbook.py`, which copies them from `A/E36_*` and `A/E37_*`
  without recomputing anything.
- Tables S3, S5, S6, S10, S11 and S22–S24 are tables in the Supplementary Information document, not sheets
  (S24, S22 and S23: see section 2).

| Sheets | Content | Script (repository) | Main output in `A/` |
|---|---|---|---|
| S1a–S1g | single-cell attribution (S1g: B-cell reassignment) | `analysis/14_v18/E17_singlecell_attribution.py`; `analysis/15_v19_v41/scripts/E17b_bcell_reassignment.py` | `E17_*.csv`, `E17b_bcell_reassignment.csv` |
| S2a–S2b | paired ΔC; C by score | `analysis/05_v14/E1_paired_cindex.R` | `E1_*.csv` |
| S4a–S4e | gene loss | `analysis/14_v18/E19_robustness.R` | `E19_*.csv` |
| S7a–S7f | locked cohorts as published | `analysis/11_v17/E15_gastric_validation.R` | `E15_*.csv` |
| S7g–S7m | extension as published | `analysis/14_v18/E20_extended_replication.R` | `E20_*.csv` |
| S7n–S7o | direct Immune-versus-CAF_ECM test, gastric | `analysis/15_v19_v41/scripts/E29_dissociation.R`, `E29b_posthoc_adjusted.R` | `E29_*.csv`, `E29b_*.csv` |
| "Note on S7" | S7 holds values with the eight shared arrays; the duplicate-free values are in S19b | `analysis/15_v19_v41/scripts/E30_dedup.R` | — |
| S8 | gastric modules | `analysis/05_v14/E4_module_meta.R` | `E4_module_meta.csv` (README also names "figures/supp") |
| S9a–S9b | PRJEB25780 nulls and classifier | `analysis/05_v14/E3_ici_null_and_ceiling.py` | `E3_*.csv` |
| S12a–S12j | NCT04249739, GSE183924, post hoc pools | `analysis/14_v18/E22_gastric_icb.R`, `E22_firth_cps.py`; `analysis/15_v19_v41/scripts/E23_gastric_consistency.R` | `E22_*.csv`, `E23_*.csv` |
| S13a–S13b | checkpoint cohorts; duplicates | `analysis/10_v16/E12a_cohort_audit.R` | `E12_cohort_table.csv`, `E12_duplicate_check.csv` |
| S13c–S13e | per cohort, pooled, leave-one-out | `analysis/10_v16/E12_icb_expanded.R` | `E12_icb_*.csv` |
| S13f–S13g | direct test, response | `E29_dissociation.R`, `E29b_posthoc_adjusted.R` | `E29_*.csv` |
| S14a–S14b | PERFECT | `analysis/11_v17/E16a_perfect.R` | `E16_perfect_*.csv` |
| S14c | METABRIC | `analysis/11_v17/E16b_metabric.R` | `E16_metabric_*.csv` |
| S15a–S15d | TCGA | `analysis/10_v16/E11_pancancer_prognosis.R` | `E11_*.csv` |
| S16a–S16d | I-SPY2 | `analysis/10_v16/E10_ispy2_trials.py` | `E10_ispy2_*.csv` |
| S16e | trials pooled | `analysis/10_v16/E10_meta.R` | `E10_meta.csv` |
| S16f–S16k | added trials | `analysis/14_v18/E21_trials.py`, `E21_meta.R` | `E21_*.csv` |
| S16l–S16p | E25 | `analysis/15_v19_v41/scripts/E25_trials.py`, `E25_meta.R` | `E25_*.csv` |
| S16q–S16x | E28 | `analysis/15_v19_v41/scripts/E28_base.py`, `E28_trials.py`, `E28_estimate.R`, `E28_gastric.R`, `E28_tiers.py` | `E28_*.csv` |
| S17a–S17f | imaging mass cytometry | `analysis/10_v16/E13_singlecell_barrier.py` | `E13_*.csv` |
| S17g–S17j | Visium rings, TLS definition | `analysis/15_v19_v41/scripts/E27_visium_marker_tls.py` (S17j also `E26b_bridge_posthoc.py`) | `E27_*.csv`, `E26b_*.csv` |
| S17k–S17l | Visium sections and tests | `analysis/15_v19_v41/scripts/E26_gastric_visium.py` | `E26_*.csv` |
| S18 | GDSC2 | `analysis/00_upstream/gdsc_validation_fix.py` (r, P) | q values in `R1/C1_gdsc_fdr.csv`: writer not established |
| S19a | duplicates | `E30_dedup.R` | `E30_*` |
| S19b | without duplicates | `E30_dedup_compare.py` | `E30_*` |
| S19c | estimators | `E30_dedup_estimators.R`; `E30_trials_meta.R` | `E30_*` |
| S19d | angiogenesis | `E30_trials.py`; `E30_trials_meta.R` | `E30_trials_*.csv` |
| S19e | endpoint strata | `E30_trials_meta.R` | `E30_trials_*.csv` |
| S19f, S19g | Stouffer; CheckMate 025 overall survival | `E30_trials_report.py` | `E30_trials_*.csv` |
| S19h–S19j | all duplicate pairs; candidates; rerun | `analysis/15_v19_v41/scripts/E31_dedup_all_pairs.R`, `E31_rerun.R` | `E31_*` |
| S20a–S20d, S20f | library; per comparison; pooled; sensitivity sets; I-SPY2 checks | `analysis/15_v19_v41/scripts/E32_trials.py` (S20a also uses `icb/sig_methods_clean.csv`) | `E32_primary.csv`, `E32_per_comparison.csv`, `E32_pooled_all_sets.csv`, `E32_ispy2.csv` |
| S20e | published engines | `E32_native.py`, `E32_native_scores.R` | `E32_native*.csv` |
| S20g | locked immune/stromal null | `E32_locked_bionull.R` | `E32_locked_bionull.csv` |
| S21a | direction audit | `analysis/15_v19_v41/protocols/E33_direction_audit.md` / `.csv` (written by an AI assistant, not a script) | — |
| S21b–S21e | single genes | `E33_trials.py` (S21e also `E33_extras.py`) | `E33_pooled.csv`, `E33_paired.csv`, `E33_per_comparison.csv`, `E33_expectations.csv` |
| S25a | NeoPembrOv links | `E34_link.py` | `data_download/validation/NeoPembrOV_GSE227666/id_match_*.csv` |
| S25b–S25c | NeoPembrOv scores; seven and four pools | `E34_trials.py` | `E34_neopembrov_per_score.csv`, `E34_pooled.csv` |
| S25d–S25e | breast genes; breast paired | `E35_trials.py` | `E35_by_gene.csv`, `E35_paired.csv` |
| S25f | expectations | `E34_trials.py`; `E35_trials.py` | `E34_expectations.csv`, `E35_expectations.csv` |
| S26a, S26b, S26d | simulation: rejection; random reference; expectations | `E36_sim.py` | `E36_S0_S2.csv`, `E36_S3.csv`, `E36_expectations.csv` |
| S26c | simulation: paired test | `E36_sim.py` (README) | `E36_S4.csv`, merged with `E36_S4_power_vs_theory.csv` (computed after the E36 run; re-derived by `scripts/text_numbers_check.py`) |
| S26e, S26f, S26h | development per cohort; pooled; frozen weights | `E37_develop.R` | `E37_loco_per_cohort.csv`, `E37_loco_pooled.csv`, `E37_frozen_weights.tsv` |
| S26g | twin check | `E37_twins.py` | `E37_twins.csv` |
| S26i, S26j | test per comparison; pooled | `E37_test.py` | `E37_test_per_comparison.csv`, `E37_test_pooled.csv` |
| S26k | against the panel genes | `E37_panel_extra.py` | `E37_panel_extra.csv` |
| S26l | expectations | `E37_develop.R`; `E37_test.py` | `E37_dev_expectations.csv`, `E37_test_expectations.csv` |
| S27a–S27j | source papers; design analysis; single-arm evidence | see section 8 | see section 8 |
| S28a | published trials adjudicated, with citations | `E41_adjudicate.py`; `E41_citations.py` | `E41_adjudicated.csv`, `E41_citations.csv` |
| S28b, S28c | the two extraction passes | `E41_pass1/`, `E41_pass2/` (provided with the submission, not in this repository) | — |
| S28d, S28e | decisions; agreement | `E41_adjudicate.py` | `E41_decisions.csv`, `E41_agreement.csv` |
| S28f–S28i | inputs; pools; expectations; paired PD-L1 | `E41_meta.R` | `E41_inputs.csv`, `E41_meta.csv`, `E41_expectations.csv`, `E41_paired.csv` (post hoc pools: `E41_posthoc.csv`) |
| S29a–S29c | barriernet 0.3 against each signature | `E42_signatures.py` | `E42_paired.csv`, `E42_counts.csv`, `E42_expectations.csv` |
| S30a–S30e | absolute benefit, pre-specified | `E43_benefit.py` | `E43_per_comparison.csv`, `E43_pooled.csv`, `E43_contrasts.csv`, `E43_expectations.csv`, `E43_coefficients.csv` |
| S30f–S30i | absolute benefit, post hoc correction | `E43b_null.py` | `E43b_cv_null.csv`, `E43b_apparent.csv`, `E43b_pooled.csv`, `E43b_contrasts.csv` |
| S31a–S31h | gastric and gastroesophageal cancer | `E44_gastric.py` | `E44_A_gastric_trials.csv`, `E44_B_per_cohort.csv`, `E44_B_pooled.csv`, `E44_B_paired.csv`, `E44_B_adjusted.csv`, `E44_C_prognosis.csv`, `E44_D_tmescore.csv`, `E44_expectations.csv` |
| S31i | post hoc checks | `E44b_checks.py` | `E44b_checks.csv` |
| S32a–S32l | reference gene from single-arm data, coexpression-matched reference, renal arms | `E45_select.py`, `E45_test.py` | `E45_single_arm_gene_pooled.csv`, `E45_reference_genes.tsv`, `E45_A_reference_pooled.csv`, `E45_A_paired.csv`, `E45_A_counts.csv`, `E45_A_gastric_pooled.csv`, `E45_B_coexpression.csv`, `E45_B_counts.csv`, `E45_C_within_arm.csv`, `E45_C_pooled.csv`, `E45_C_enrichment.csv`, `E45_expectations.csv` |
| S32m–S32r | T-cell means, renal against breast pools, gene chosen in advance | `E46_test.py` | `E46_A_reference_pooled.csv`, `E46_A_paired.csv`, `E46_A_counts.csv`, `E46_A_steadiness.csv`, `E46_B_held_out.csv`, `E46_B_pooled.csv` |
| S32s–S32t | MSS, EBV-negative gastric patients | `E46_test.py` | `E46_C_per_cohort.csv`, `E46_C_pooled.csv` |
| S32u | I-SPY2 expression checks | `E46_test.py`, `E46b_ispy2.py`, `E46c_pembro_pairs.py` | `E46_D_ispy2_twins.csv`, `E46b_ispy2.csv`, `E46c_pembro_pairs.csv`, `E46c_sensitivity.csv` |
| S32v–S32w | counts as 52 rows and as 51 distinct signatures; expectations | `E46_distinct.py`, `E46_test.py` | `E46_E_distinct.csv`, `E46_expectations.csv` |

---

## 7. Not established

- `A/E23_comparator_genes.csv`: no script writes it. It holds the ESTIMATE and TMEscore gene lists that
  `cmp_scores.py` (and through it E5, E10, E21, E25, E28, E32 and later) and `E32_trials.py` read.
- `A/E36_S4_power_vs_theory.csv` (workbook S26c): not written by `E36_sim.py` or any other script. Fig. 2f
  recomputes the same prediction from `E36_S4.csv`.
- Supplementary Figs S31 and S33: the edited copies of `make_fig8.py` and `make_fig3_v33.py` that drew them are not
  saved.
- The source of the per-trial counts "3 of 52 to 28 of 52" is `figure_scripts/v34/make_fig3.py`, not a current
  figure script.
- The breast detectable difference "0.6 to 0.7 per SD", the "21 of 29" and "33 of 43" expectation tallies, and the
  E24 numbers (only in the package folder `analysis_E24_archived/`) are not printed or written by any script in the
  repository.
- `manuscript/v32_fills.json` (Supplementary Table S22 rows) was typed by hand; no script writes it.
- r = 0.81 between barriernet 0.3 and CD8A: stated in `E42_outcome.md`; no script in the repository prints it.
- The median +0.19 of the inflamed signatures and the counts 13 of 14 and 28 of 73: no script prints them. They were recounted from `A/E32_primary.csv` and `A/E44_labels.csv` when this map was updated, and they agree.
- Upstream objects listed in `README.md`: `FD/bio_null_draws.csv`, the GC single-cell composition and fraction tables
  (S3, S4), the BarrierIndex clinical table (S28), and the drug tables (S26).


## 8. Additions in the current version (source-paper review, design analysis, single-arm evidence) and closed gaps

| Item | Script | Output (`A/` = analysis results) |
|---|---|---|
| Protocols (eleventh to thirteenth) and hashes | — | `protocols/E38_prereg.md`, `E39_prereg.md`, `E40_prereg.md` (+ `.sha256`); outcomes `E38_outcome.md`, `E39_outcome.md`, `E40_outcome.md` |
| Fig. 2g, source papers: design of the evidence, single-gene and random-gene comparisons | two independent extraction passes (`results/E38_pass1.csv`, `E38_pass2.csv`; sources `protocols/E38_sources.csv`) → `scripts/E38_adjudicate.py` | `A/E38_adjudicated.csv`, `E38_agreement.csv`, `E38_crosstab.csv`, `E38_expectations.csv` |
| Fig. 2h, single-arm cohorts against the randomized interaction; checkpoint and control arms alone | `scripts/E40_export.R` → `scripts/E40_concordance.py` | `A/E40_single_arm_genes.csv.gz` (not in the repository: patient-level), `E40_single_arm_per_cohort.csv`, `E40_randomized_per_comparison.csv`, `E40_concordance.csv`, `E40_expectations.csv`, log `E40_concordance.log` |
| Fig. 2i, events needed to detect a difference from CD8A; ratio of bootstrap to normal-theory variance; detectable differences | `scripts/E39_design.py` | `A/E39_per_comparison.csv`, `E39_pooled.csv`, `E39_required.csv`, `E39_expectations.csv`, log `E39_design.log` |
| Analysis-plan addendum (about 1,190 / 1,400 events for 0.10 per SD) | from `A/E39_required.csv` and `E39_pooled.csv` | `protocols/E37_registration_addendum.md` (+ `.sha256`) |
| Text numbers: A 0, B 4, C 3, D 36, E 9; 18 / 2 / 17 papers; agreement 50/52, κ 0.92 | `scripts/E38_adjudicate.py` | `A/E38_adjudicated.csv`, `E38_agreement.csv` |
| Text numbers: ρ = 0.83; 17 of 18; 17 of 34; 0 of 52 in the control arms; 7 vs 18 | `scripts/E40_concordance.py` | `A/E40_concordance.csv`, log `E40_concordance.log` |
| Text numbers: ratios 1.45 / 1.17 / 1.76; Spearman 0.90–0.99; 2,825 and 1,186 events; 744 events; MDD 0.21, 0.13 | `scripts/E39_design.py` | `A/E39_*`, log `E39_design.log` |
| Tally at that stage: 13 later protocols, 39 of 52 expectations (current tally: section 9) | — | `analysis/15_v19_v41/protocols/later_protocol_tally.md` (v40 tally 33/43 + `A/E38_expectations.csv`, `E39_expectations.csv`, `E40_expectations.csv`) |
| Fig. 7c (prognosis against response) | `figure_scripts/v41/make_fig7.py` | x from `A/E30_dedup_runs/primary/E20_pooled.csv` (`all 10 gastric`, unadjusted) → `analysis/15_v19_v41/scripts/E30_dedup_E20.R`; y as before. Earlier renders read `A/E20_pooled.csv` as first written |
| Previously "not established": per-trial share above random sets (3/52 to 28/52), breast detectable difference 0.6–0.7 per SD, workbook S26c `predicted` column | `scripts/text_numbers_check.py` (re-derives them from saved results; nothing refitted) | `A/text_numbers_check.csv` |
| Previously overwritten: step 0 of the barriernet 0.3 export (78 per-cohort log-odds) | `scripts/E37_export_step0_rerun.R` (same code, output renamed) | `A/E37_export_step0_rerun.csv` (largest difference 3.8 × 10⁻¹⁵) |
| Workbook sheets S27a–S27j | `manuscript/make_v41_workbook.py` (copies the files above) | `BarrierNet_Supplementary_Tables.xlsx` |



## 9. Additions in v42–v46

v42 rewrote the prose and turned Fig. 2a into a worked example. v43 redrew Figs 2–7 with larger type and left the numbers
unchanged. The four protocols below were written after the tag v0.3.0 release. Their outputs in `A/` are not in the
repository, apart from the expectation tables (`analysis/16_v42_v46/results/`). The results are in the workbook, sheets
S28–S31. The extraction tables `E41_pass1/` and `E41_pass2/` and the patient-level `A/E44_labels.csv` and `A/E44_Z_*` are
not in the repository; the extraction tables are provided with the submission.

| Protocol | Files in `analysis/16_v42_v46/protocols/` | Scripts in `analysis/16_v42_v46/scripts/` | Output (`A/`) |
|---|---|---|---|
| fourteenth: published randomized trials | `E41_prereg.md`, `.sha256`, `E41_brief.md`, `E41_outcome.md` | `E41_adjudicate.py`, `E41_citations.py`, `E41_meta.R` | `E41_adjudicated.csv`, `E41_decisions.csv`, `E41_agreement.csv`, `E41_citations.csv`, `E41_inputs.csv`, `E41_meta.csv`, `E41_paired.csv`, `E41_posthoc.csv`, `E41_expectations.csv` |
| fifteenth: barriernet 0.3 against each signature | `E42_prereg.md`, `.sha256`, `E42_outcome.md` | `E42_signatures.py` | `E42_paired.csv`, `E42_counts.csv`, `E42_counts_by_family.csv`, `E42_expectations.csv`, `E42.log` |
| sixteenth: absolute benefit | `E43_prereg.md`, `.sha256`, `E43_outcome.md` | `E43_benefit.py` | `E43_per_comparison.csv`, `E43_pooled.csv`, `E43_contrasts.csv`, `E43_coefficients.csv`, `E43_expectations.csv`, `E43.log` |
| addendum to the sixteenth (post hoc) | `E43b_prereg.md`, `.sha256` | `E43b_null.py` | `E43b_cv_null.csv`, `E43b_apparent.csv`, `E43b_pooled.csv`, `E43b_contrasts.csv`, `E43b_heterogeneity.csv`, `E43b.log` |
| seventeenth: gastric and gastroesophageal cancer | `E44_prereg.md`, `.sha256`, `E44_outcome.md` | `E44_export.R`, `E44_gastric.py` | `E44_labels.csv`, `E44_A_gastric_trials.csv`, `E44_B_per_cohort.csv`, `E44_B_pooled.csv`, `E44_B_paired.csv`, `E44_B_adjusted.csv`, `E44_C_prognosis.csv`, `E44_D_tmescore.csv`, `E44_expectations.csv`, `E44.log` |
| addendum to the seventeenth (post hoc, descriptive) | `E44b_prereg.md`, `.sha256` | `E44b_checks.py` | `E44b_checks.csv`, `E44b.log` |
| Tally: seventeen later protocols, 67 expectations, 53 held | `later_protocol_tally.md` | — | — |
| Manuscript, SI, cover letter | — | `manuscript/make_v46_text.py` → `v46_text.py` → `apply_v46.py`, `apply_v46_si.py`, `make_v46_workbook.py`, `make_v46_package.py`, `scrub_v46_metadata.py` (provided with the submission) | `final_v46/` (not in the repository) |

---

## 10. Additions in v47–v48

v47 gave the paper two parts (framework and randomized comparisons, then gastric cancer), with a new Fig. 6. v48 counts 51 distinct signatures, adds the simple T-cell mean as a reference and corrects the empirical P equation to the one-sided form that `E32_trials.py` computes. The protocols below were written after the tag v0.3.0 release. Their outputs in `A/` are not in the repository, apart from the expectation tables and the frozen reference gene (`analysis/17_v47_v48/results/`). The results are in the workbook, sheet S32. `A/E45_single_arm_genes*` and the other patient-level exports are not in the repository.

| Protocol | Files in `analysis/17_v47_v48/protocols/` | Scripts in `analysis/17_v47_v48/scripts/` | Output (`A/`) |
|---|---|---|---|
| eighteenth: reference gene from single-arm data, coexpression-matched reference, renal arms | `E45_prereg.md`, `.sha256`, `E45_frozen.sha256`, `E45_outcome.md` | `E45_export.R`, `E45_select.py`, `E45_test.py` | `E45_single_arm_gene_pooled.csv`, `E45_reference_genes.tsv`, `E45_A_*.csv`, `E45_B_*.csv`, `E45_C_*.csv`, `E45_expectations.csv`, `E45_test.log` |
| nineteenth: T-cell means, gene chosen in advance, MSS/EBV-negative gastric patients, I-SPY2 controls, distinct counts | `E46_prereg.md`, `.sha256`, `E46_outcome.md` | `E46_test.py`, `E46_distinct.py` | `E46_A_*.csv`, `E46_B_*.csv`, `E46_C_*.csv`, `E46_D_*.csv`, `E46_E_distinct.csv`, `E46_expectations.csv`, `E46_test.log` |
| addendum: is the I-SPY2 flag a shared-patient signal (post hoc, no expectations) | `E46b_addendum.md`, `.sha256` | `E46b_ispy2.py` | `E46b_ispy2.csv`, `E46b_ispy2.log` |
| addendum: near-identical samples inside the I-SPY2 pembrolizumab comparison (post hoc, no expectations) | `E46c_addendum.md`, `.sha256` | `E46c_pembro_pairs.py` | `E46c_pembro_pairs.csv`, `E46c_sensitivity.csv` |
| Tally: nineteen later protocols, 84 expectations, 68 held | `later_protocol_tally.md` | — | — |
| Manuscript, SI, cover letter | — | `manuscript/make_v48_text.py` → `v48_text.py` → `apply_v48.py`, `apply_v48_si.py`, `make_v48_workbook.py`, `make_v48_package.py`, `scrub_v48_metadata.py` (provided with the submission) | `final_v48/` |

**Numbers new in the text of v47 and v48**

| Numbers | Source file | Script |
|---|---|---|
| 14 of 26 scored in all six; 2 of 25 not scored in IMmotion150; 18 of 19 inflamed positive | `A/E46_E_distinct.csv` | `E46_distinct.py` |
| 25 above size-matched (14 inflamed, 11 other); 27 with audited directions; 11 of 25 not nominal | `A/E46_E_distinct.csv`; `A/E33_Q7_audited_directions.csv` | `E46_distinct.py`; `E33_trials.py` |
| coexpression-matched: 15 (9 inflamed, 6 other), references centered at −0.02, 10 lost, none gained | `A/E45_B_coexpression.csv`; `A/E46_E_distinct.csv` | `E45_test.py`; `E46_distinct.py` |
| CXCL10: +0.47 (0.26 to 0.68) in 1,996 eligible genes; CD8A 19th, +0.35 | `A/E45_single_arm_gene_pooled.csv` | `E45_select.py` |
| CXCL10 in the randomized comparisons: +0.17; minus CD8A −0.13; renal +0.11, breast +0.45; above 7, 11, 0; below 2, 1, 9 | `A/E45_A_reference_pooled.csv`; `A/E45_A_counts.csv`; `A/E46_E_distinct.csv` | `E45_test.py`; `E46_distinct.py` |
| three-gene mean +0.22 (0.07 to 0.37), renal +0.24, breast +0.18; 0 above, 16 below; detectable 0.20 | `A/E46_A_reference_pooled.csv`; `A/E46_A_counts.csv` | `E46_test.py` |
| 13-gene mean +0.18 (0.04 to 0.32); one above (+0.06, P = 0.046), 11 below | `A/E46_A_reference_pooled.csv`; `A/E46_A_counts.csv`; `A/E46_A_paired.csv` | `E46_test.py` |
| breast comparisons: 1 above and 3 below the three-gene mean; 7 above the 13-gene mean; 6 above CD8A | `A/E46_A_counts.csv`; `A/E46_E_distinct.csv` | `E46_test.py`; `E46_distinct.py` |
| renal against breast: median 0.14 (panel genes), 0.09 (inflamed signatures), P = 0.27; 0.06 for each mean | `A/E46_A_steadiness.csv`; `A/E46_test.log` | `E46_test.py` |
| gene chosen in advance: +0.19 (0.04 to 0.34) against +0.28; chosen genes; −0.02 and −0.01 against the mean and the gene list | `A/E46_B_pooled.csv`; `A/E46_B_held_out.csv` | `E46_test.py` |
| I-SPY2: 15 of 29 flagged; 32 of 71 treated flagged; best correlation 0.60; duplicates at least 0.91; three pairs inside the pembrolizumab comparison | `A/E46_D_ispy2_twins.csv`; `A/E46b_ispy2.csv`; `A/E46c_pembro_pairs.csv`; `A/E46c_sensitivity.csv` | `E46_test.py`; `E46b_ispy2.py`; `E46c_pembro_pairs.py` |
| renal arms: CD8A slopes by arm; control-arm pool; angiogenesis score | `A/E45_C_within_arm.csv`; `A/E45_C_pooled.csv` | `E45_test.py` |
| Spearman ρ of random-gene significance with set size 0.41 on 51 signatures | recomputed from `A/E32_primary.csv` and `A/E32_coverage.csv` in `figure_scripts/v48/make_fig5.py` | `E32_trials.py` |
| sensitivity sets on 51: 19 and 2; 17 and 1; 13 and 1; 22 and 1 | `A/E46_E_distinct.csv` | `E46_distinct.py` |
| barriernet 0.3 minus the three-gene mean +0.006 (P = 0.89) | `A/E46_A_paired.csv` | `E46_test.py` |
