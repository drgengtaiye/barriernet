"""Copy (never move) every script the v14-v16 manuscripts depend on into this directory.

Re-run after any source script changes:  python collect.py
Sources stay where they are; the scripts keep their absolute input/output paths, so they
run from here exactly as from their original location. Missing sources are listed, not fatal.
"""
import glob, os, shutil

G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
GC = f"{G}/GC_TME_Reprogrammer"
HERE = os.path.dirname(os.path.abspath(__file__))

MANIFEST = {
    # upstream inputs produced before the revision (only those whose script exists)
    "00_upstream": [f"{GC}/validation/benchmark_barriernet_official.py",
                    f"{GC}/validation/attribution/cell_attribution.py",
                    f"{GC}/validation/bench/spatial_benchmark.py",
                    f"{GC}/validation/drug_vali/gdsc_validation_fix.py"],
    # bulk cohorts, scores, survival, C-index, matched null (per cohort)
    "01_bulk_prognosis": [f"{G}/revision/scripts/{s}" for s in
                          ["B1_prep_cohorts.R", "B2_scores.R", "B2b_estimate_fix.R",
                           "B3_baseline_sig.R", "B4_survival.R"]],
    # PRJEB25780 AUC, DeLong, combination models
    "02_checkpoint": [f"{G}/revision/scripts/{s}" for s in ["A2_roc_full.R", "A3_cv_combo.py"]],
    # ceiling and upgrade experiments A-C
    "03_ceiling_upgrades": [f"{G}/revision2/scripts/{s}" for s in
                            ["S1_pseudobulk.py", "S2_identity_state.py", "S3_state_axes.py",
                             "S4_bulk_features.R", "S5_loco.R", "S6_ici_features.R",  # S6, S10: PRJEB25780 logged twice; no v16 number uses them
                             "S7_neighborhoods.py", "S8_learn_arrangement.py",
                             "S9_tune_arrangement.py", "S10_decisive.R", "S11_ceiling.R",
                             "S13_fig.py"]],
    # least-cost transport
    "04_transport": [f"{G}/revision2/scripts/{s}" for s in
                     ["T3_transport.py", "T4_sensitivity.py", "T5_fig.py", "T6_marker_based.py"]],
    # analyses added in v14
    "05_v14": [f"{G}/2/analysis/{s}" for s in
               ["E1_paired_cindex.R", "E2_matched_null.R", "E3_ici_null_and_ceiling.py",
                "E4_module_meta.R"]],
    # corrected-lattice rings (stored peaks) - the numbers in the v14 text
    "06_spatial": [f"{G}/2/figures/ring_recompute.py"],
    # figures, style and figure-side statistics
    "07_figures": [f"{G}/figstyle.py"] + [f"{G}/2/figures/{s}" for s in
                  ["bnstyle.py", "check_palette.py", "panelprep.py", "rowpack.py", "sc_extract.py",
                   "prot_stats.py", "gdsc_prep.py", "rank_baseline.py"]]
                  + sorted(glob.glob(f"{G}/2/figures/make_fig[0-9].py")),
    "07_figures/supp": sorted(glob.glob(f"{G}/2/figures/supp/*.py")),
    "08_manuscript": [f"{G}/2/manuscript/build.py"],
    # analyses and figures added in v15 (E3/E5/E7 were corrected in v16 for the PRJEB25780 scale)
    "09_v15": [f"{G}/2/analysis/{s}" for s in
               ["E5_icb_split.R", "E5_trials.py", "E5_gastric_arms.R", "E5_prjeb_arms.py",
                "E6_common.R", "E6_models.R", "E6_learning.R", "E6_review.R",
                "E7a_scores_log2.py", "E7b_roc_delong.R", "E7c_firth_combo.py", "E7d_subgroup_samepipe.py",
                "E8_prognostic_hardening.R", "E8b_biology_null.R", "E9_estimate_check.R"]]
              + [f"{G}/2/figures/{s}" for s in
                 ["bnstyle.py", "make_fig_split_validation.py", "make_v15_fig1.py", "make_v15_fig4.py",
                  "make_v15_fig5.py", "make_v15_fig6.py", "make_v15_fig7.py"]]
              + [f"{G}/2/manuscript/build_v15.py"],
    # analyses and figures added in v16
    "10_v16": [f"{G}/2/analysis/{s}" for s in
               ["E10_ispy2_trials.py", "E10_meta.R", "E10_fig_trials.py",
                "E11_pancancer_prognosis.R", "E12a_cohort_audit.R", "E12_icb_expanded.R", "E12_fig.py",
                "E13_singlecell_barrier.py", "E13_figure.py"]]
              + [f"{G}/2/figures/{s}" for s in ["make_v16_fig1.py", "make_fig_E11.py"]]
              + [f"{G}/2/manuscript/build_v16.py"],
    # analyses and figures added in v17
    "11_v17": [f"{G}/2/analysis/{s}" for s in
               ["E15_gastric_validation.R", "E15_figure.py", "E16a_perfect.R", "E16b_metabric.R", "E16_figure.py"]]
              + [f"{G}/2/figures/make_v17_fig1.py", f"{G}/2/manuscript/build_v17.py",
                 f"{G}/2/data_download/get_gastric_geo.sh"],
    # lean (submission) version
    "12_lean": [f"{G}/2/manuscript/{s}" for s in ["build_lean.py", "make_lean_supp.py"]]
               + [f"{G}/2/figures/supp/make_lean_workbook.py", f"{G}/2/figures/make_lean_fig1.py",
                  f"{G}/2/analysis/E16_figure.py"],
    "13_final": [f"{G}/2/manuscript/{s}" for s in ["build_final.py", "apply_final_check.py"]],
    "14_v18": [f"{G}/2/manuscript/{s}" for s in ["refs.py", "make_v18_sources.py", "apply_v18_edits.py", "apply_v18_fixes.py",
                                                 "apply_v18_e22.py", "make_v18_tables.py", "build_v18.py"]]
              + [f"{G}/2/analysis/{s}" for s in ["E18_cps.py", "E17_singlecell_attribution.py", "E19_robustness.R", "E19_figure.py",
                                                "E20_extended_replication.R", "E20_figure.py", "E21_trials.py", "E21_meta.R",
                                                "E21_fig_trials.py", "E21_braun_check.R", "E21_braun_posthoc.R",
                                                "E22_gastric_icb.R", "E22_firth_cps.py", "E22_figure.py", "E17_E22_prereg.md"]]
              + [f"{G}/2/figures/make_v18_fig1.py", f"{G}/2/figures/supp/make_v18_workbook.py",
                 f"{G}/2/data_download/validation/bulk_PRJEB60680_Samsung_NCT04249739/quant/aggregate_qc.py",
                 f"{G}/2/data_download/validation/bulk_PRJEB60680_Samsung_NCT04249739/quant/quant_watch.py"],
    # kept for provenance, not used for any v14 number
    "99_superseded": [f"{G}/revision/scripts/A1_ici_roc.R",          # read the stale PRJEB25780 score file
                      f"{G}/revision/scripts/v2_random_logrank.R",   # new gene set per cohort; replaced by E2
                      f"{G}/revision/scripts/D1_figures.py",
                      f"{G}/revision2/scripts/S12_ici_ceiling.py", f"{G}/revision2/scripts/S12b.py",  # replaced by E3
                      f"{G}/revision2/scripts/T1_transport.py", f"{G}/revision2/scripts/T2_transport.py",
                      f"{G}/2/figures/make_fig8_with_model.py"],
}

missing = []
for sub, files in MANIFEST.items():
    d = os.path.join(HERE, sub)
    os.makedirs(d, exist_ok=True)
    for f in files:
        if os.path.exists(f):
            shutil.copy2(f, d)
        else:
            missing.append(f)
n = sum(len(v) for v in MANIFEST.values())
print(f"copied {n - len(missing)} of {n} scripts")
for m in missing:
    print("MISSING", m)
