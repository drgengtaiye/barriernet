"""Supplementary_Tables_lean.xlsx: the full workbook's sheets renumbered for the lean manuscript
(S1, S5, S6, S7, S10-S15; S2-S4, S8, S9 are tables in SUPPLEMENTARY_lean). Rows of the PRJEB25780 null
table computed on the doubly log-transformed matrix are dropped."""
import pandas as pd
from openpyxl.styles import Font, Alignment
src = "Supplementary_Tables.xlsx"
x = pd.read_excel(src, sheet_name=None)
MAP = [("S1a_head_to_head", "S1_head_to_head"), ("S1b_C_by_score", "S1b_C_by_score"),
       ("S5a_locked_cohorts", "S9a_locked_cohorts"), ("S5b_locked_overlap", "S9b_locked_overlap"),
       ("S5c_locked_coverage", "S9c_locked_coverage"), ("S5d_locked_per_cohort", "S9d_locked_per_cohort"),
       ("S5e_locked_pooled", "S9e_locked_pooled"), ("S5f_locked_null_dC", "S9f_locked_null_dC"),
       ("S6_gastric_modules", "S3_modules"),
       ("S7a_PRJEB_nulls", "S2_PRJEB25780"), ("S7b_PRJEB_classifier", "S2b_ICI_ceiling"),
       ("S10a_ICB_cohorts", "S5a_ICB_cohorts"), ("S10b_ICB_duplicates", "S5b_ICB_duplicates"),
       ("S10c_ICB_per_cohort", "S5c_ICB_per_cohort"), ("S10d_ICB_pooled", "S5d_ICB_pooled"),
       ("S10e_ICB_leave_one_out", "S5e_ICB_leave_one_out"),
       ("S11a_PERFECT", "S10a_PERFECT"), ("S11b_ICB_with_PERFECT", "S10b_ICB_with_PERFECT"), ("S11c_METABRIC", "S10c_METABRIC"),
       ("S12a_TCGA_per_type", "S7a_TCGA_per_type"), ("S12b_TCGA_pooled", "S7b_TCGA_pooled"),
       ("S12c_TCGA_null", "S7c_TCGA_null"), ("S12d_TCGA_BRCA", "S7d_TCGA_BRCA"),
       ("S13a_ISPY2_interaction", "S6a_ISPY2_interaction"), ("S13b_ISPY2_within_arm", "S6b_ISPY2_within_arm"),
       ("S13c_ISPY2_receptor_strata", "S6c_ISPY2_HR_strata"), ("S13d_ISPY2_null", "S6d_ISPY2_null"),
       ("S13e_trials_pooled", "S6e_trials_pooled"),
       ("S14a_IMC_models", "S8a_IMC_models"), ("S14b_IMC_patients", "S8b_IMC_patients"),
       ("S14c_IMC_PDAC_regions", "S8c_IMC_PDAC_regions"), ("S14d_IMC_macrophages", "S8d_IMC_macrophages"),
       ("S14e_IMC_poscontrol", "S8e_IMC_poscontrol"), ("S14f_IMC_survival", "S8f_IMC_survival"),
       ("S15_GDSC2", "S4_GDSC2")]
old_readme = x["README"].set_index("Sheet")
rows, sheets = [], []
for new, old in MAP:
    df = x[old]
    if old == "S2_PRJEB25780":
        df = df[~df["Expression scale"].astype(str).str.contains("twice")]
    if old == "S5c_ICB_per_cohort":
        df = df[df["set"] != "gastric_sens"]                     # rows on the doubly log-transformed matrix
    if old == "S5a_ICB_cohorts":
        df = df.replace("analysed as a trial in v15", "analysed as a randomised comparison (sensitivity)", regex=True)
        df = df.replace("not rescuable offline", "not scorable from the compendium; rescored from GEO GSE165252 (Supplementary Table S11)", regex=True)
    if old in ("S5c_ICB_per_cohort", "S5d_ICB_pooled", "S5e_ICB_leave_one_out", "S10b_ICB_with_PERFECT"):
        df = df.replace({"E5 11 (reference)": "first 11", "all + chemo-ICB (E12)": "all + chemo-ICB", "new only, E5 entry rule (>= 25 labelled)": "new only, first-11 entry rule (>= 25 labelled)"})
        for c in df.columns:
            if df[c].dtype == object and c in ("set", "role"):
                df[c] = df[c].replace({"E5": "first 11"})
    if old == "S3_modules":                                         # add the Hartung-Knapp pooled columns
        hk = pd.read_csv("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results/E8_meta_DL_HK.csv").drop_duplicates("score").set_index("score")
        df = df.copy()
        pooled = df["Cohort"].astype(str).str.startswith("Pooled")
        for col, src in [("HK 95% CI low", "HK_lo"), ("HK 95% CI high", "HK_hi"), ("HK P", "HK_p")]:
            df[col] = None
            df.loc[pooled, col] = df.loc[pooled, "Module"].map(hk[src]).round(4)
    if old == "S9e_locked_pooled":                                  # "new/original" -> "locked/benchmark"
        df = df.replace({"original 3 (reference)": "benchmark 3 (reference)", "new + original 3": "locked + benchmark 3",
                         "new cohorts": "locked cohorts", "new cohorts, OS only": "locked cohorts, OS only"})
    if old in ("S5a_ICB_cohorts", "S5c_ICB_per_cohort", "S5d_ICB_pooled", "S5e_ICB_leave_one_out"):
        for c in ("role", "set"):
            if c in df.columns:
                df[c] = df[c].replace({"new": "replication", "new only (replication)": "replication only",
                                       "new, whole transcriptome": "replication, whole transcriptome",
                                       "new only, first-11 entry rule (>= 25 labelled)": "replication only, first-11 entry rule (>= 25 labelled)"})
    sheets.append((new, df))
    rows.append({"Sheet": new, "Content": old_readme.loc[old, "Content"] if old in old_readme.index else "",
                 "Source": old_readme.loc[old, "Source"] if old in old_readme.index else ""})
readme = pd.DataFrame(rows)
readme["Content"] = readme["Content"].str.replace(", and with log2 applied twice (the v15 error)", "", regex=False)
readme["Content"] = readme["Content"].str.replace("(substitute for METABRIC bulk)", "(for comparison with METABRIC)", regex=False)
readme["Content"] = readme["Content"].str.replace("against the original three", "against the three benchmark cohorts", regex=False)
readme["Content"] = readme["Content"].str.replace("(new cohorts, OS-only, new + original)", "(locked cohorts; overall-survival cohorts only; locked + benchmark cohorts)", regex=False)
readme.loc[readme.Sheet == "S15_GDSC2", "Source"] = "code/analysis/00_upstream/gdsc_validation_fix.py (Spearman r, P); Benjamini-Hochberg q values in C1_gdsc_fdr.csv"
readme.loc[readme.Sheet == "S6_gastric_modules", "Content"] = "Cox HR per SD for each module per gastric benchmark cohort, with DerSimonian-Laird pooled estimates and Hartung-Knapp pooled intervals (HK columns)."
out = "../lean/Supplementary_Tables_lean.xlsx"
with pd.ExcelWriter(out, engine="openpyxl") as w:
    for name, df in [("README", readme)] + sheets:
        df.to_excel(w, sheet_name=name, index=False)
        ws = w.sheets[name]
        for cell in ws[1]:
            cell.font = Font(bold=True); cell.alignment = Alignment(wrap_text=True, vertical="top")
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = min(60, max(12, max(len(str(c.value or "")) for c in col) * 0.9))
        ws.freeze_panes = "A2"
print("written", out, len(sheets), "sheets")
