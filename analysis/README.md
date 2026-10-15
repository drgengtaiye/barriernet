# barriernet — analysis code map

Every script the manuscript depends on, grouped by stage, with the run order and the map from each
table, figure and reported number to the script that writes it. `python collect.py` re-copies the
scripts from their working locations (re-run it after any source change). The scripts keep the
absolute input/output paths they ran under on the lab filesystem
(`/groups/adv2105_gp/yichen/Yi/com/gastric/...`), so paths must be adjusted to a local copy of the
public datasets listed in the manuscript's data-availability statement. Two files exist only here:
`06_spatial/ring_peaks_sensitivity.py` (barriernet 0.2.1 rings under three TLS-peak definitions) and
`07_figures/export_fig4_data.R` (writes the tables `make_fig4.py` reads; output checked identical to
the files in use).

Short names used below: `R1` = `gastric/revision/results`, `R2` = `gastric/revision2/results`,
`A` = `gastric/2/analysis/results`, `F` = `gastric/2/figures`, `GC` = `gastric/GC_TME_Reprogrammer`.

## Run order

| Step | Scripts | Writes |
|---|---|---|
| 0 upstream | `00_upstream/benchmark_barriernet_official.py` (PRJEB25780 scores, GEP-18, IPS, TIDE, PD-L1; linear-TPM z-scores); `cell_attribution.py`; `spatial_benchmark.py` (Moran's I); `gdsc_validation_fix.py` (GDSC2 Spearman r, P) | `GC/validation/benchmark_scores.tsv`; `GC/validation/attribution/cell_attribution_*.csv`; `GC/validation/bench/spatial_benchmark_results.csv`; `gdsc_validation_results.csv` |
| 1 bulk | `B1_prep_cohorts.R` → `B2_scores.R` → `B3_baseline_sig.R` (signature-level ESTIMATE, `ESTsig_diff`) → `B2b_estimate_fix.R` (official ESTIMATE, **not used**: `EST_diff`) → `B4_survival.R` | `R1/B1_*.rds`, `R1/B2_scores.rds` + `B2_scores_*.csv`, `R1/B4_univariate.csv`, `B4_meta.csv`, `B4_multivariate.csv`, `B4_cindex.csv` |
| 2 checkpoint | `A2_roc_full.R`, `A3_cv_combo.py` (read `benchmark_scores.tsv`) | `R1/A2_auc_ci.csv`, `A2_delong.csv`, `A3_cv_combo.csv`, `A3_lrt.txt` |
| 3 ceiling/upgrades | `S1` → `S2` → `S3` → `S4` → `S5_loco.R`; `S6_ici_features.R`; `S7` → `S8` → `S9` → `S10_decisive.R`; `S11_ceiling.R`; `S13_fig.py` | `R2/S5_loco_cindex.csv`, `S6_ici_univariate.csv`, `S10_prognosis.csv`, `S10_ici_with_SAS.csv`, `S11_ceiling_prognosis.csv`, `revision2/figures/Fig_A1_ceiling.*`; S9 tuning in `revision2/logs_S9.txt` |
| 4 transport | `T3_transport.py`, `T4_sensitivity.py` (Tangram fractions), `T6_marker_based.py` (marker scores), `T5_fig.py` | `R2/T3_transport.csv`, `T4_sensitivity.csv`, `T6_marker_based.csv`, logs `revision2/logs_T4.txt`, `logs_T6.txt` |
| 5 v14 analyses | `E1_paired_cindex.R`, `E2_matched_null.R` (needs `A/bulk_gene_pool.txt`), `E3_ici_null_and_ceiling.py`, `E4_module_meta.R` | `A/E1_*.csv`, `A/E2_*`, `A/E3_*` + `E3.log`, `A/E4_module_meta.csv` |
| 6 spatial | `ring_recompute.py` (corrected lattice, stored TLS peaks; the numbers in the text); `ring_peaks_sensitivity.py` (0.2.1 package; stored vs kNN vs lattice peaks) | `F/ring_contrasts_recomputed*.csv`, `F/ring_hex_*.csv`; `06_spatial/ring_peaks_sensitivity*.csv` |
| 7 figures | `export_fig4_data.R`, `sc_extract.py`, `prot_stats.py`, `gdsc_prep.py`, `rank_baseline.py`, then `make_fig1.py` … `make_fig8.py` (all import `bnstyle.py` → `figstyle.py`); `supp/*.py` | `F/Fig1–8.{png,pdf}`, `F/supp/FigS*.{png,pdf}`, `F/Fig7_*.csv`, `F/Fig8_*.csv` |
| 8 manuscript | `08_manuscript/build.py` | `gastric/2/MANUSCRIPT_v14.{md,docx}` |

`99_superseded/` keeps scripts that produced earlier numbers but none in v14: `A1_ici_roc.R`
(read the stale PRJEB25780 score file), `v2_random_logrank.R` (drew a new gene set per cohort;
replaced by `E2`), `S12*` (replaced by `E3`), `T1`/`T2` (earlier transport versions),
`D1_figures.py`, `make_fig8_with_model.py` (ranking-model panels removed from Fig. 8).

## Manuscript item → script → output

| Item | Script(s) | Output |
|---|---|---|
| Cohort sizes (383/431/283; 1,097; 498 deaths) | `B1_prep_cohorts.R` | `R1/B1_*.rds` |
| Table 1 HRs, pooled HR, adjusted HRs | `B4_survival.R` | `B4_univariate.csv`, `B4_meta.csv`, `B4_multivariate.csv` |
| Table 1 C, null rank / empirical P | `E1_paired_cindex.R`, `E2_matched_null.R` | `E1_cindex_by_score.csv`, `E2_null_per_cohort.csv`, `E2_null_joint.csv` |
| Table 2 | `B4_survival.R` (pooled HRs), `E4_module_meta.R` (CAF_ECM HR), `E1_paired_cindex.R` (C, ΔC) | `B4_meta.csv`, `E4_module_meta.csv`, `E1_*.csv` |
| Table 3 (and ESTIMATE +0.009–0.022) | `B4_survival.R` | `B4_cindex.csv` |
| Table 4 | `A2_roc_full.R` | `A2_auc_ci.csv`, `A2_delong.csv` |
| Table 5, LRT χ² 3.47, P 0.062, OR 4.1 | `A3_cv_combo.py` | `A3_cv_combo.csv`, `A3_lrt.txt` |
| Table 6 prognosis | `S5_loco.R`, `S11_ceiling.R` | `S5_loco_cindex.csv`, `S11_ceiling_prognosis.csv` |
| Table 6 checkpoint | `E3_ici_null_and_ceiling.py` (log₂ AUCs, classifier 0.754), `S6_ici_features.R` (pericyte 0.786), `S10_decisive.R` (arrangement 0.503) | `E3_ici_matched_null.csv`, `E3_ici_ceiling.csv`, `S6_ici_univariate.csv`, `S10_ici_with_SAS.csv` |
| Fig. 1 | `make_fig1.py` | `F/Fig1.*` |
| Fig. 2 | `make_fig2.py` (+ `sc_extract.py`, `panelprep.py`); attribution from `cell_attribution.py` | `F/Fig2.*` |
| Fig. 3 | `make_fig3.py` (+ `ring_recompute.py` outputs) | `F/Fig3.*` |
| Fig. 4 | `export_fig4_data.R`, `make_fig4.py` | `F/Fig4.*` |
| Fig. 5 | `make_fig5.py` (reads `benchmark_scores.tsv`, `A2_*`, `E3_ici_null_draws_linear.csv`) | `F/Fig5.*` |
| Fig. 6 | `make_fig6.py` | `F/Fig6.*` |
| Fig. 7 and protein numbers (36/48 proteins; ρ 0.41–0.80; P 0.002–0.04; CAF–TAM ρ 0.77) | `prot_stats.py`, `make_fig7.py` | `F/Fig7_module_coherence.csv`, `F/Fig7_module_correlation.csv` |
| Fig. 8 (d,e: r −0.38/+0.39, P) | `gdsc_prep.py`, `make_fig8.py` | `F/Fig8_gdsc_points.csv` |
| Supplementary Figs S1–S7 | `07_figures/supp/make_supp_sc.py`, `make_supp_spatial.py`, `rings_export.py` (in progress by the figure work) | `F/supp/FigS*.{png,pdf}` |
| Supplementary Fig. S8 (transport) | `T5_fig.py` | `revision2/figures/Fig_A2_transport.*` |
| Supplementary Fig. S10 (upgrades) | `S13_fig.py` | `revision2/figures/Fig_A1_ceiling.*` |
| Supplementary Table S1 / S2 / S3 | `E1` / `E4` / `E3` | `E1_paired_dC.csv` / `E4_module_meta.csv` / `E3_ici_matched_null.csv` |
| Abstract: pooled HR 0.74; null P ≤ 0.001; ΔC +0.068, +0.046, +0.000 (−0.014 to +0.015); C 0.600 vs 0.589; upgrades < 0.01; AUC 0.797; null P 0.013; classifier 0.754; IFN 0.854; CAF_ECM HR 1.29; TAM_SUP not prognostic, higher in responders | `B4`, `E2`, `E1`, `S11`/`S5`, `S5`/`S10`, `A2`, `E3`, `E3`, `A2`, `E4`, `E4` + `A2` (TAM_SUP AUC 0.740) | as above |
| Mann–Whitney P 0.002; Cliff's δ +0.59 | `S6_ici_features.R` (P); `benchmark_barriernet_official.py` (δ, printed) | `S6_ici_univariate.csv`; console |
| CAF_ECM AUC 0.71 (low) / Barrier 0.52 | `make_fig5.py` (computed from `benchmark_scores.tsv`), `A2_roc_full.R` | panel 5d; `A2_auc_ci.csv` |
| Ring contrasts (P 0.014, 0.037, 0.32, 0.0098), BarrierIndex −3.1 to +3.2, +3.2 → −2.0 | `ring_recompute.py` | `F/ring_contrasts_recomputed*.csv` |
| Moran's I 0.45 / 0.41, 48/50 | `spatial_benchmark.py` | `GC/validation/bench/spatial_benchmark_results.csv` |
| Attribution z-scores; 80% / >70% | `cell_attribution.py` | `cell_attribution_zscore.csv`, `cell_attribution_contribution_pct.csv` |
| Upgrade A (289 profiles, 20 axes, r 0.98–1.00), C (16,091 neighbourhoods, r 0.22, shuffle 0.00, best of 9) | `S1`–`S3`, `S7`–`S9` | `R2/S1_*`, `S3_*`, `logs_S9.txt` |
| Transport (z −3.38, +0.74; 6/4 sections; P 0.19; 3.86 vs 3.46) | `T4_sensitivity.py`, `T6_marker_based.py` | `logs_T4.txt`, `logs_T6.txt`, `T6_marker_based.csv` |
| Drug scores (dasatinib 35.5, first in 10/10) and 13-compound table | `make_fig8.py` reads `GC/results/ai_recommendations_drug_top15_clean.csv` | see flags |

## Numbers without an identifiable producing script

- **Upstream objects**: the GSE183904 processing (`GC/out_scanpy_strict_GSE183904`), the
  Tangram mapping, and the spatial state scores and stored TLS peaks
  (`GC/results/ring_barrier_state_v2/*.with_state_and_ring.h5ad`). No script on disk writes them.
  The stored peaks were reverse-engineered: they are lattice maxima on the 0.2.0 broken graph.
- **Drug scores**: the pertpy per-section scores (`GC/results/ai_recommendations_drug_top15_clean.csv`).
- **GDSC q values**: the BH-corrected table `R1/C1_gdsc_fdr.csv` (q 0.24; Supplementary Table S4).
  `gdsc_validation_fix.py` produces r and P, but the BH step was not saved.
- **Single-cell statistics**:
  - Cohen's d ≈ 4.5 and 2.6 (Fig. 2e,g) are read off legacy panel images; no table exists.
  - The composition test (Kruskal–Wallis P > 0.2) has no script.
  - The CAF/TAM fraction correlation (r = 0.018, P = 0.93) has no script.
- **BarrierIndex vs age / pT / pN** (Spearman P 0.50, 0.44, 0.08): computed inline during
  drafting; no saved script. `supp/rings_export.py` may cover it once Supplementary Fig. S6
  is final.
- **Official ESTIMATE failure on GSE84437** (r = −0.40 with PTPRC, CD2, CD3E, CD8A):
  `B2b_estimate_fix.R` prints r with Immune and CD8A only; the four-gene figure is not reproduced.
- **PRJEB25780 log₂ values** (0.815, 0.831, 0.852) are from `E3`. The DeLong / combination
  models exist only on the linear scale.

## Software (conda env `r45`, used for the v14 analyses and figures)

Python 3.14.0: numpy 2.3.5, pandas 2.3.3, scipy 1.17.1, scikit-learn 1.9.0, statsmodels 0.14.6,
lifelines 0.30.3, matplotlib 3.11.0, anndata 0.12.17, scanpy 1.12.3, squidpy 1.8.3,
PyMuPDF 1.27.2.2, Pillow 12.2.0, PyYAML 6.0.3, pytest 9.1.0. R 4.5.2: survival 3.8.3,
metafor 5.0.1, glmnet 4.1.10, pROC 1.19.0.1, GSVA 2.4.4, estimate 1.0.13, matrixStats 1.5.0.
Tangram and pertpy are not installed in this environment. The upstream objects were
produced elsewhere, and the scripts from 3 September (revision/, revision2/) may have run
under other package versions.

---

## v15 additions (2026-09-23): `09_v15/`

These are snapshot copies. The originals live in `gastric/2/analysis/` and `gastric/2/figures/`.

| Script | Produces | Manuscript item |
|---|---|---|
| `E5_icb_split.R` | `results/E5_icb_*` | 11 external checkpoint cohorts: pooled log-OR, size-matched null, C1Q/SPP1 arms (Table 7; Fig 7a,b,d) |
| `E5_trials.py` | `results/E5_trials.csv` | IMmotion150 and PRINCE interactions (Table 7; Fig 7c) |
| `E5_gastric_arms.R`, `E5_prjeb_arms.py` | `results/E5_gastric_arms_*`, `E5_prjeb_arms.csv` | macrophage arms in the gastric cohorts and PRJEB25780 |
| `E6_common.R`, `E6_models.R`, `E6_learning.R`, `E6_review.R` | `results/E6_*` | eight model classes, ΔC, offset, within-cohort CV, inverse-normal, learning curves (Table 4; Fig 5) |
| `E7a`–`E7d` | `results/E7_*` (see `E7_SUMMARY.md`) | PRJEB25780 on log2: faithful IPS, DeLong, Firth, subgroups, same-pipeline CV (Tables 5–6; Fig 6) |
| `E8_prognostic_hardening.R` | `results/E8_*` | Hartung–Knapp, LR tests, TMEscore, subtype/Lauren (Tables 1–3; Supplementary Table S3) |
| `E8b_biology_null.R` | `results/E8_biology_null*` | immune/stromal-gene null (Table 1; Fig 4e) |
| `E9_estimate_check.R` | console | official ESTIMATE r = −0.50 on GSE84437 (Methods). This replaces the older "−0.40" entry above. |
| `make_v15_fig1.py`, `make_v15_fig4.py`–`make_v15_fig7.py` | `figures/v15_Fig*` | main Figures 1 and 4–7 (Figures 2–3 are `make_fig2.py`, `make_fig3.py`) |
| `build_v15.py` | `MANUSCRIPT_v15.{md,docx}` | reference numbering and export |

Primary spatial numbers come from `06_spatial/ring_peaks_sensitivity*.csv`, using kNN-9 peaks and Far ≥ 3 hops.

---

## v16 additions (2026-09-24): `10_v16/`

These are snapshot copies. The originals live in `gastric/2/analysis/`, `gastric/2/figures/` and `gastric/2/manuscript/`. Refresh them with `python collect.py`.

**PRJEB25780 scale correction.** The deposited matrix `mGC_Pembrolizumab_Kim2018.TPM` is already log2(TPM+1). v15 applied log2(x+1) to it a second time.

- `E3_ici_null_and_ceiling.py`, `E5_prjeb_arms.py`, `E7a_scores_log2.py` and `figures/supp/make_supp_tables.py` now use the matrix as provided. `E7a --twice` reproduces the v15 values.
- The v15 versions are saved as `*.v15bak` beside the originals, and the old outputs are in `analysis/results/v15bak/`.
- Every changed number is listed in `results/E14_prjeb_rescale_diff.csv`.
- This replaces the v14 note above that the DeLong and combination models "exist only on the linear scale".
- Legacy scripts `03_ceiling_upgrades/S6_ici_features.R` and `S10_decisive.R` still apply the log a second time. No v16 number comes from them.

| Script | Produces | Manuscript item |
|---|---|---|
| `E12a_cohort_audit.R` | `results/E12_cohort_table.csv`, `E12_duplicate_check.csv` | replication cohorts: inclusion, platform, module coverage, duplicate check (Methods; Supplementary Table S5a,b) |
| `E12_icb_expanded.R` | `results/E12_icb_per_cohort.csv`, `E12_icb_pooled.csv`, `E12_icb_loo.csv` | first 11 / replication 10 / all 21 pooled response effects and sensitivity sets (Table 7; Results; Supplementary Table S5c–e) |
| `E12_fig.py` | `figures/Fig_E12_icb_expanded.*` → `v16_Fig7.*` | Figure 7 |
| `E10_ispy2_trials.py` | `results/E10_ispy2_*` | I-SPY2 interactions, within-arm effects, hormone-receptor strata and random-gene null (Table 8; Supplementary Table S6a–d) |
| `E10_meta.R` | `results/E10_meta.csv` | benefit-oriented pooling across IMmotion150 and I-SPY2, ± PRINCE (Table 8 pooled column; Supplementary Table S6e) |
| `E10_fig_trials.py` | `figures/Fig_E10_trials.*` → `v16_Fig8.*` | Figure 8. The pooled panel uses the conservative interval rule. |
| `E11_pancancer_prognosis.R` | `results/E11_*` | TCGA 24-type prognosis, checkpoint-type pooling, random-gene null, prognosis–response table, TCGA-BRCA (Table 7; Supplementary Table S7) |
| `make_fig_E11.py` | `figures/Fig_E11_prognosis_vs_response.*` → `supp/FigS14.*` | Supplementary Figure S14 |
| `E13_singlecell_barrier.py` | `results/E13_*` | IMC continuity models, positive control, macrophage–aggregate contrasts, METABRIC survival (Results; Supplementary Table S8) |
| `E13_figure.py` | `figures/Fig_E13_singlecell_barrier.*` → `supp/FigS15.*` | Supplementary Figure S15 |
| `make_v16_fig1.py` | `figures/v16_Fig1.*` | Figure 1 |
| `figures/supp/make_supp_tables.py` | `figures/supp/Supplementary_Tables.xlsx` | Supplementary Tables S1–S8 |
| `build_v16.py` | `MANUSCRIPT_v16.{md,docx}` | reference numbering and export |

**Other figures.**

- `v16_Fig2`–`v16_Fig5` are copies of the v15 figures.
- `v16_Fig6` is `make_v15_fig6.py` rerun on the corrected scale.
- `supp/FigS13` is `v15_Fig7` (its panel d was rerun on the corrected scale).

**Conservative pooling rule.** Pooled estimates from the checkpoint cohorts, trials and TCGA use the wider of the Hartung–Knapp and DerSimonian–Laird intervals and the larger P value, as the gastric cohorts already did.

---

## v17 additions (2026-09-24): `11_v17/`

Data were downloaded on login1 with `data_download/get_gastric_geo.sh` into `data_download/geo/`: 11 gastric GEO series, GSE165252 and the METABRIC cBioPortal bundle. The GPL8432 probe table was saved to `data_download/geo/platforms/`.

| Script | Produces | Manuscript item |
|---|---|---|
| `E15_gastric_validation.R` | `results/E15_*` | locked replication in GSE15459, GSE34942 and GSE26253: per-cohort and pooled HRs, overlap check, random-gene null, ΔC (Results; Supplementary Table S9) |
| `E15_figure.py` | `figures/Fig_E15_gastric_validation.*` → `supp/FigS16.*` | Supplementary Figure S16 |
| `E16a_perfect.R` | `results/E16_perfect_*` | PERFECT (GSE165252) scoring and pooled checkpoint estimates with it added (Results; Supplementary Table S10a,b) |
| `E16b_metabric.R` | `results/E16_metabric_*` | METABRIC HRs by ER status, and a comparison with TCGA-BRCA (Results; Supplementary Table S10c) |
| `E16_figure.py` | `figures/Fig_E16_perfect_metabric.*` → `supp/FigS17.*` | Supplementary Figure S17 |
| `make_v17_fig1.py` | `figures/v17_Fig1.*` | Figure 1 |
| `build_v17.py` | `MANUSCRIPT_v17.{md,docx}` | reference numbering and export |

---

## Lean submission version (2026-09-28): `12_lean/`

- **`make_lean_supp.py`** builds `lean_supp_src.md` from `v17_src.md`: the Supplementary Methods, five former main tables and renumbered legends S1–S20.
- **`build_lean.py`** builds `MANUSCRIPT_lean` and `SUPPLEMENTARY_lean` (.md and .docx), numbering the references of each separately.
- **`make_lean_workbook.py`** builds `figures/lean/Supplementary_Tables_lean.xlsx` from the full workbook. It renames the sheets to the lean numbering, drops the doubly log-transformed rows, removes internal labels and adds the Hartung–Knapp columns to S6.
- **`make_lean_fig1.py`** redraws Figure 1 with the lean supplementary numbers.

Main-figure mapping to v17:

| Lean | v17 |
|---|---|
| Fig 2 | Fig 4 |
| Fig 3 | Fig 6 |
| Fig 4 | Fig 7 |
| Fig 5 | Fig 8 |

Supplementary-figure mapping, lean ← v17:

| Lean | v17 | Lean | v17 |
|---|---|---|---|
| S1 | Fig 2 | S11 | S14 |
| S2 | S1 | S12 | S17 |
| S3 | S2 | S13 | Fig 3 |
| S4 | S3 | S14 | S5 |
| S5 | S4 | S15 | S6 |
| S6 | S9 | S16 | S7 |
| S7 | Fig 5 | S17 | S8 |
| S8 | S10 | S18 | S15 |
| S9 | S16 | S19 | S11 |
| S10 | S13 | S20 | S12 |

## Final submission package (2026-10-03): `13_final/`

- **`apply_final_check.py`** applies the text fixes from `analysis/results/FINAL_CHECK.md` to `lean_src.md` and `lean_supp_src.md`. From here on `lean_supp_src.md` is edited directly; `make_lean_supp.py` was a one-time generator whose scratch input no longer exists.
- **`build_final.py`** (run after `build_lean.py`) writes `final/`: the manuscript and Supplementary Information as formatted .docx (US Letter, Times 12 pt, double spacing and line numbers in the main text, page numbers, References before Tables, table notes below each table), the workbook and the figures under submission names. The PDFs in `final/` were rendered from these .docx files with LibreOffice.
- Figure labels were corrected in `E10_fig_trials.py`, `E12_fig.py`, `E15_figure.py`, `E16_figure.py`, `E13_figure.py`, `make_lean_fig1.py`, `make_fig2.py`, `make_fig3.py` and `supp/make_supp_{spatial,models,sc}.py`, `supp/make_figS12.py`; each keeps its previous version as `<script>.prefinal`. Only labels and label positions changed, except Supplementary Fig. S14a, which now shows the nine sections other than the representative GC8 (it previously omitted GC2 and repeated GC8).

## v18: pre-registered validation (2026-10-05): `14_v18/`

- **Pre-registration:** `E17_E22_prereg.md`. The hash log `E17_E22_prereg.sha256` is kept next to it in `analysis/`.
- **Analyses:**
  - E17: independent single-cell attribution;
  - E18: PD-L1 CPS in PRJEB25780;
  - E19: robustness to gene loss;
  - E20: locked replication extended to seven cohorts;
  - E21: CheckMate 025, JAVELIN Renal 101 and NeoTRIP interactions, pooled over six comparisons;
  - E22: NCT04249739 (PRJEB60680, quantified from raw reads with fastp and salmon) and GSE183924.

  Each writes `analysis/results/E*_README.md`.
- **Manuscript build:** `refs.py`, then `apply_v18_edits.py`, `apply_v18_fixes.py` and `apply_v18_e22.py`, then `build_v18.py`. Supplementary numbering is generated by first citation into `figures/v18/numbering.json`. Tables 1, 3 and 4 come from `make_v18_tables.py`.
- **Figures and workbook:** `make_v18_fig1.py` and `make_v18_workbook.py`.
- **Final formatting:** `build_final.py v18` writes `final_v18/`.
