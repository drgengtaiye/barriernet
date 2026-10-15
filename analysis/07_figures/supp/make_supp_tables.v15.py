"""Supplementary Tables S1-S4 in one workbook (openpyxl), with a README sheet."""
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
import sys
out = sys.argv[1] if len(sys.argv) > 1 else "Supplementary_Tables.xlsx"
with pd.ExcelWriter(out, engine="openpyxl") as w:
    for name, df in [("README", readme), ("S1_head_to_head", s1), ("S1b_C_by_score", s1c), ("S2_PRJEB25780", s3),
                     ("S2b_ICI_ceiling", s3c), ("S3_modules", s2), ("S4_GDSC2", s4)]:
        df.to_excel(w, sheet_name=name, index=False)
        ws = w.sheets[name]
        for cell in ws[1]:
            cell.font = Font(bold=True); cell.alignment = Alignment(wrap_text=True, vertical="top")
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = min(60, max(12, max(len(str(c.value or "")) for c in col) * 0.9))
        ws.freeze_panes = "A2"
print("written", out)
