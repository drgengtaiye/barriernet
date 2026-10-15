"""Supplementary Tables S1-S8 in one workbook (openpyxl), with a README sheet. S5-S8 (v16) are the
result tables of the replication checkpoint cohorts (E12), I-SPY2 (E10), TCGA pan-cancer prognosis (E11)
and single-cell imaging (E13), written as produced by those scripts (rounded to 4 d.p.)."""
import pandas as pd
from openpyxl.styles import Font, Alignment
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
R = "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
NM = {"NetImmune_ssGSEA": "NetImmune (ssGSEA scoring)", "Immune": "Immune (no Barrier term)", "Immune_minus_CAF": "Immune minus CAF_ECM (exploratory)",
      "ESTIMATE": "ESTIMATE immune-stromal (signature-level)", "CAF_ECM": "CAF_ECM module"}
COH = {"TCGA": "TCGA-STAD", "pooled": "Pooled (inverse-variance)"}
# S1
d = pd.read_csv(f"{A}/E1_paired_dC.csv")
s1 = pd.DataFrame({"Cohort": d.cohort.map(lambda c: COH.get(c, c)), "Score A": d.a.map(lambda x: NM.get(x, x)),
                   "Score B": d.b.map(lambda x: NM.get(x, x)), "Patients": d.n, "Delta C (A minus B)": d.dC.round(4),
                   "SE": d.se.round(4), "95% CI low": d.lo.round(4), "95% CI high": d.hi.round(4), "P": d.p.round(6)})
c = pd.read_csv(f"{A}/E1_cindex_by_score.csv")
s1c = pd.DataFrame({"Cohort": c.cohort.map(lambda x: COH.get(x, x)), "Score": c.score.map(lambda x: NM.get(x, x)),
                    "Harrell's C (direction fixed in advance)": c.C.round(4), "95% CI low": c.lo.round(4), "95% CI high": c.hi.round(4)})
# S2
m = pd.read_csv(f"{A}/E4_module_meta.csv"); pc = pd.read_csv("module_hr_by_cohort.csv")
s2 = pc.assign(cohort=pc.cohort.map(lambda x: COH.get(x, x))).rename(columns={"cohort": "Cohort", "module": "Module", "n": "Patients",
          "events": "Deaths", "HR": "HR per SD", "lo": "95% CI low", "hi": "95% CI high", "p": "P"})
s2m = m.rename(columns={"module": "Module", "HR": "Pooled HR per SD", "lo": "95% CI low", "hi": "95% CI high", "p": "P", "I2": "I2 (%)"})
s2 = pd.concat([s2, s2m.assign(Cohort="Pooled (DerSimonian-Laird)").rename(columns={"Pooled HR per SD": "HR per SD"})], ignore_index=True)
s2 = s2[["Module", "Cohort", "Patients", "Deaths", "HR per SD", "95% CI low", "95% CI high", "P", "I2 (%)"]].round(4)
# S3
n3 = pd.read_csv(f"{A}/E3_ici_matched_null.csv")
# the deposited matrix is already log2(TPM+1); "log2x2" re-applies log2(x+1) (the v15 error), kept for transparency
s3 = pd.DataFrame({"Expression scale": n3.scale.map({"log2": "log2(TPM+1) as deposited (primary)",
                                                    "log2x2": "log2 applied twice (v15 error; sensitivity)"}),
                   "Score": n3.score.map({"NetImmune": "NetImmune", "Immune": "Immune (no Barrier term)", "IFN": "IFN module"}),
                   "AUC observed": n3.AUC_obs.round(4), "Null median AUC (2,000 random axes, better direction)": n3.null_median.round(4),
                   "Null 95th percentile": n3.null_p95.round(4), "Empirical P": n3.p_emp.round(4)})
ce = pd.read_csv(f"{A}/E3_ici_ceiling.csv")
s3c = pd.DataFrame({"L2 penalty C": ce.C, "CV AUC (mean of 40 x 5-fold)": ce.cv_auc.round(4),
                    "2.5th percentile": ce.lo.round(4), "97.5th percentile": ce.hi.round(4)})
# S4
g = pd.read_csv(f"{R}/C1_gdsc_fdr.csv")
s4 = pd.DataFrame({"Compound": g.Drug, "GDSC2 name": g.GDSC_name, "Cell lines": g.n,
                   "Spearman r, Barrier vs ln(IC50)": g.Barrier_r.round(3), "P": g.Barrier_p.round(4), "BH q": g.Barrier_q.round(3),
                   "Spearman r, NetImmune vs ln(IC50)": g.NetImmune_r.round(3), "P (NetImmune)": g.NetImmune_p.round(4),
                   "BH q (NetImmune)": g.NetImmune_q.round(3)})
readme = pd.DataFrame({"Sheet": ["S1_head_to_head", "S1b_C_by_score", "S3_modules", "S2_PRJEB25780", "S2b_ICI_ceiling", "S4_GDSC2"],
    "Content": ["Paired differences in Harrell's C between scores, per cohort and pooled; 2,000 bootstrap resamples of patients.",
                "Harrell's C of each score per cohort, direction fixed in advance (higher immune / lower stromal = lower risk).",
                "Univariable Cox HR per SD for each of the six modules per cohort, and DerSimonian-Laird pooled estimates.",
                "PRJEB25780 (44 patients, 12 responders): observed AUC and structure-matched random-axis null on the deposited log2(TPM+1) scale, and with log2 applied twice (the v15 error).",
                "Full-transcriptome L2 logistic regression on the deposited log2(TPM+1), 40 x stratified 5-fold CV, by penalty.",
                "GDSC2, 28 gastric cancer cell lines: Spearman correlation of Barrier and NetImmune scores with ln(IC50), BH across 11 compounds."],
    "Source": ["analysis/E1_paired_cindex.R", "analysis/E1_paired_cindex.R", "analysis/E4_module_meta.R; figures/supp", "analysis/E3_ici_null_and_ceiling.py",
               "analysis/E3_ici_null_and_ceiling.py", "revision/scripts (C1_gdsc_fdr.csv)"]})

# S5-S8 (v16): result tables as written by E10-E13
def rd(f):
    d = pd.read_csv(f"{A}/{f}")
    return d.round({c: 4 for c in d.select_dtypes("float").columns})
def cons(d):
    """wider of the Hartung-Knapp and DerSimonian-Laird intervals, larger P (the manuscript's rule)"""
    hk = (d.hi_HK - d.lo_HK) >= (d.hi_DL - d.lo_DL)
    d = d.copy()
    d["lo_conservative"] = d.lo_HK.where(hk, d.lo_DL); d["hi_conservative"] = d.hi_HK.where(hk, d.hi_DL)
    d["p_conservative"] = d[["p_DL", "p_HK"]].max(axis=1)
    return d
v16 = [("S2c_PRJEB_v15_vs_v16", rd("E14_prjeb_rescale_diff.csv"), "Every PRJEB25780 number quoted in version 15 (matrix log-transformed twice) beside its corrected value on the deposited log2(TPM+1) scale, with its source file.", "analysis/E3, E5_prjeb_arms, E7a-d (see E14 note in code map)"),
       ("S5a_ICB_cohorts", rd("E12_cohort_table.csv"), "Checkpoint cohort audit: tumour, treatment, platform, n, responders, module gene coverage, inclusion status and reason.", "analysis/E12a_cohort_audit.R"),
       ("S5b_ICB_duplicates", rd("E12_duplicate_check.csv"), "Duplicate-patient checks between cohorts (identifier overlap and best-match expression correlation).", "analysis/E12a_cohort_audit.R"),
       ("S5c_ICB_per_cohort", rd("E12_icb_per_cohort.csv"), "Per-cohort log-odds of response per SD for every score.", "analysis/E12_icb_expanded.R"),
       ("S5d_ICB_pooled", cons(rd("E12_icb_pooled.csv")).round(4), "Pooled log-odds per SD for every set (first 11, replication 10, all 21 and sensitivity sets); DL, HK and the conservative interval.", "analysis/E12_icb_expanded.R"),
       ("S5e_ICB_leave_one_out", rd("E12_icb_loo.csv"), "Leave-one-cohort-out ranges of the pooled estimates.", "analysis/E12_icb_expanded.R"),
       ("S6a_ISPY2_interaction", rd("E10_ispy2_interactions.csv"), "I-SPY2 score x arm interaction ORs for pCR: overall, hormone-receptor adjusted (primary) and stratified; Firth estimates.", "analysis/E10_ispy2_trials.py"),
       ("S6b_ISPY2_within_arm", rd("E10_ispy2_within_arm.csv"), "Within-arm OR for pCR per SD, unadjusted and hormone-receptor adjusted.", "analysis/E10_ispy2_trials.py"),
       ("S6c_ISPY2_HR_strata", rd("E10_ispy2_hr_strata.csv"), "Interaction ORs within hormone-receptor strata.", "analysis/E10_ispy2_trials.py"),
       ("S6d_ISPY2_null", rd("E10_ispy2_null_summary.csv"), "Structure-matched random-gene null for the NetImmune and CAF_ECM interactions (2,000 sets).", "analysis/E10_ispy2_trials.py"),
       ("S6e_trials_pooled", rd("E10_meta.csv"), "Benefit-oriented pooling across randomised comparisons (REML + Hartung-Knapp, fixed effect, Stouffer), primary and sensitivity sets.", "analysis/E10_meta.R"),
       ("S7a_TCGA_per_type", rd("E11_per_type.csv"), "TCGA per cancer type: Cox HR per SD for every score, OS and PFI, unadjusted and adjusted.", "analysis/E11_pancancer_prognosis.R"),
       ("S7b_TCGA_pooled", rd("E11_pooled.csv"), "Pooled HRs for the pre-specified groups (DL, HK and conservative columns lo/hi/p).", "analysis/E11_pancancer_prognosis.R"),
       ("S7c_TCGA_null", rd("E11_null.csv"), "Random-gene null for pooled CAF_ECM and Immune effects in the five checkpoint cancer types.", "analysis/E11_pancancer_prognosis.R"),
       ("S7d_TCGA_BRCA", rd("E11_brca_substitute.csv"), "TCGA-BRCA (substitute for METABRIC bulk): HR per SD adjusted for age, stage and ER status.", "analysis/E11_pancancer_prognosis.R"),
       ("S8a_IMC_models", rd("E13_models.csv"), "Infiltration models: continuity, expected and excess continuity, positive-control fits.", "analysis/E13_singlecell_barrier.py"),
       ("S8b_IMC_patients", rd("E13_patient_metrics.csv"), "METABRIC per-patient metrics (continuity, densities, infiltration).", "analysis/E13_singlecell_barrier.py"),
       ("S8c_IMC_PDAC_regions", rd("E13_pdac_roi_model_data.csv"), "PDAC per-region metrics used in the model.", "analysis/E13_singlecell_barrier.py"),
       ("S8d_IMC_macrophages", rd("E13_macro_aggregates_patients.csv"), "Macrophage and other cell densities near versus far from lymphoid aggregates, per patient.", "analysis/E13_singlecell_barrier.py"),
       ("S8e_IMC_poscontrol", rd("E13_poscontrol.csv"), "Random-walker positive control by blocking probability and replicate.", "analysis/E13_singlecell_barrier.py"),
       ("S8f_IMC_survival", rd("E13_survival.csv"), "METABRIC disease-specific survival, HR per SD adjusted for age, NPI and ER.", "analysis/E13_singlecell_barrier.py"),
       ("S9a_locked_cohorts", rd("E15_cohort_table.csv"), "Locked gastric cohorts: platform, stored scale, n, events, endpoint, covariates, inclusion status.", "analysis/E15_gastric_validation.R"),
       ("S9b_locked_overlap", rd("E15_overlap_check.csv"), "Patient-overlap check between the locked cohorts and against the original three.", "analysis/E15_gastric_validation.R"),
       ("S9c_locked_coverage", rd("E15_module_coverage.csv"), "Module gene coverage per locked cohort.", "analysis/E15_gastric_validation.R"),
       ("S9d_locked_per_cohort", rd("E15_per_cohort.csv"), "HR per SD for every score per locked cohort, unadjusted and adjusted.", "analysis/E15_gastric_validation.R"),
       ("S9e_locked_pooled", rd("E15_pooled.csv"), "Pooled HRs (new cohorts, OS-only, new + original); DL, HK and the conservative choice.", "analysis/E15_gastric_validation.R"),
       ("S9f_locked_null_dC", pd.concat([rd("E15_null.csv").assign(table="random-gene null"), rd("E15_deltaC.csv").assign(table="delta C")], ignore_index=True), "Structure-matched random-gene null for NetImmune, and paired delta C against ESTIMATE and CD8A.", "analysis/E15_gastric_validation.R"),
       ("S10a_PERFECT", rd("E16_perfect_per_score.csv"), "PERFECT (GSE165252) per-score log-odds of response, primary and sensitivity scale.", "analysis/E16a_perfect.R"),
       ("S10b_ICB_with_PERFECT", rd("E16_perfect_pooled.csv"), "Pooled checkpoint estimates with PERFECT added; chemo-checkpoint and gastro-oesophageal sets.", "analysis/E16a_perfect.R"),
       ("S10c_METABRIC", rd("E16_metabric_cox.csv"), "METABRIC HR per SD, OS and disease-specific, unadjusted and adjusted, by ER stratum.", "analysis/E16b_metabric.R")]
readme = pd.concat([readme, pd.DataFrame({"Sheet": [v[0] for v in v16], "Content": [v[2] for v in v16], "Source": [v[3] for v in v16]})],
                   ignore_index=True)
import sys
out = sys.argv[1] if len(sys.argv) > 1 else "Supplementary_Tables.xlsx"
with pd.ExcelWriter(out, engine="openpyxl") as w:
    for name, df in [("README", readme), ("S1_head_to_head", s1), ("S1b_C_by_score", s1c), ("S2_PRJEB25780", s3),
                     ("S2b_ICI_ceiling", s3c), (v16[0][0], v16[0][1]), ("S3_modules", s2), ("S4_GDSC2", s4)] + [(v[0], v[1]) for v in v16[1:]]:
        df.to_excel(w, sheet_name=name, index=False)
        ws = w.sheets[name]
        for cell in ws[1]:
            cell.font = Font(bold=True); cell.alignment = Alignment(wrap_text=True, vertical="top")
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = min(60, max(12, max(len(str(c.value or "")) for c in col) * 0.9))
        ws.freeze_panes = "A2"
print("written", out)
