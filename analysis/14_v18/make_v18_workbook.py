"""Supplementary_Tables_v18.xlsx: the lean workbook renumbered by figures/v18/numbering.json (written by
manuscript/build_v18.py), plus sheets for the pre-registered validation analyses (E17, E19, E20, E21; E22 when run).
Internal analysis codes are removed from cell values and README text."""
import json, re, pathlib
import pandas as pd
from openpyxl.styles import Font, Alignment

G = pathlib.Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2")
R = G / "analysis/results"
num = json.loads((G / "figures/v18/numbering.json").read_text())["ST"]
LEAN_KEY = {1: "headtohead_cohort", 5: "locked", 6: "modules", 7: "prjeb_nulls", 10: "icb", 11: "perfect_metabric",
            12: "tcga", 13: "trials", 14: "imc", 15: "gdsc"}
lean = pd.read_excel(G / "figures/lean/Supplementary_Tables_lean.xlsx", sheet_name=None)
lean_readme = lean.pop("README").set_index("Sheet")

sheets, readme = [], []
for name, df in lean.items():
    n = int(re.match(r"S(\d+)", name).group(1))
    key = LEAN_KEY[n]
    new = re.sub(r"^S\d+", f"S{num[key]}", name)
    sheets.append((num[key], new, df))
    readme.append((num[key], new, lean_readme.loc[name, "Content"], lean_readme.loc[name, "Source"]))

NEW = {
    "sc_indep": [
        ("a_celltype_mapping", "E17_celltype_mapping.csv", "Author cell annotations mapped to 11 classes (fixed before scoring)."),
        ("b_score_by_class", "E17_attribution_by_class.csv", "Module score and share of signal by cell class, discovery and two independent atlases."),
        ("c_top_class", "E17_top_class.csv", "Top class per module and its z, per atlas."),
        ("d_share_by_patient", "E17_share_by_patient_summary.csv", "Share of each module's signal in its intended class, per patient (median, IQR)."),
        ("e_expectations", "E17_expectations.csv", "Pre-specified expectations, held or not, per atlas."),
        ("f_sensitivity", "E17_sensitivity.csv", "Sensitivity analyses (enriched libraries dropped, scales, chemistries, reweighting)."),
    ],
    "robust": [
        ("a_dropout_summary", "E19_summary.csv", "Pooled estimates under random gene dropout (10%, 20%, 30%; 1,000 draws each)."),
        ("b_correlation", "E19_spearman_summary.csv", "Correlation of perturbed with full NetImmune, per cohort and overall."),
        ("c_leave_one_gene_out", "E19_leave_one_gene_out.csv", "Pooled estimates after removing each gene in turn."),
        ("d_panel_platforms", "E19_panel_subset.csv", "Scores restricted to each targeted-panel platform's genes."),
        ("e_expectations", "E19_expectations.csv", "Pre-specified expectations, held or not."),
    ],
    "locked": [
        ("g_ext_cohorts", "E20_cohort_table.csv", "Pre-registered extension: cohorts, outcome sources, patients and events."),
        ("h_ext_duplicates", "E20_duplicate_exclusions.csv", "Arrays excluded by the expression-correlation duplicate screen."),
        ("i_ext_per_cohort", "E20_per_cohort.csv", "Per-cohort HR per SD for every score, unadjusted and adjusted; OS and RFS."),
        ("j_ext_pooled", "E20_pooled.csv", "Pooled HRs: all locked, added, all ten gastric, and sensitivity sets; DL, HK and conservative."),
        ("k_ext_null", "E20_null.csv", "Structure-matched random-axis null per added cohort and pooled."),
        ("l_ext_deltaC", "E20_deltaC.csv", "ΔC of NetImmune against ESTIMATE and CD8A."),
        ("m_ext_expectations", "E20_expectations.csv", "Pre-specified expectations, held or not."),
    ],
    "trials": [
        ("f_added_interactions", "E21_interactions.csv", "CheckMate 025, JAVELIN Renal 101, NeoTRIP and DUTRENEO: score × arm interactions."),
        ("g_added_within_arm", "E21_within_arm.csv", "Within-arm effects of every score in the added trials."),
        ("h_added_null", "E21_null_summary.csv", "Structure-matched random-gene null for the added trials."),
        ("i_pooled_all_sets", "E21_meta.csv", "Pooled benefit-oriented interactions: six primary, leave-one-out, cancer type, sensitivity sets."),
        ("j_expectations", "E21_expectations_pooled.csv", "Pre-specified expectations for the pooled interactions."),
        ("k_checkmate_matrix", "E21_braun_consistency.csv", "CheckMate expression: compendium copy against the authors' matrix (scores only, no outcomes)."),
    ],
    "gastric_icb": [
        ("a_NCT04249739_qc", "../../data_download/validation/bulk_PRJEB60680_Samsung_NCT04249739/quant/PRJEB60680_baseline_qc.tsv", "NCT04249739 (PRJEB60680): per-sample quality control, batch, mapping and identity checks."),
        ("b_response", "E22_PRJEB60680_response.csv", "Response log-odds and AUC per score, all and MSS/EBV-negative, primary and pre-specified processing."),
        ("c_null", "E22_PRJEB60680_null.csv", "Structure-matched random-gene null for the AUCs."),
        ("d_survival", "E22_PRJEB60680_cox.csv", "Progression-free and overall survival HR per SD (secondary)."),
        ("e_firth_cps", "E22_firth_cps.csv", "Firth models with MSI/EBV status and CPS."),
        ("f_pooled", "E22_pooled.csv", "Pooled response estimates with PRJEB25780, and with PERFECT."),
        ("g_GSE183924", "E22_GSE183924_cox.csv", "Adjuvant durvalumab cohort: relapse-free survival HR per SD."),
        ("h_expectations", "E22_expectations.csv", "Pre-specified expectations, held or not."),
    ],
}
SRC = {"gastric_icb": "analysis/E22_gastric_icb.R, analysis/E22_firth_cps.py", "sc_indep": "analysis/E17_singlecell_attribution.py", "robust": "analysis/E19_robustness.R",
       "locked": "analysis/E20_extended_replication.R", "trials": "analysis/E21_trials.py, analysis/E21_meta.R"}
CLEAN = [(r"locked \(first 3\), E15 reference", "locked (first 3), reference"),
         (r"reproduction: 3 published primary", "first three primary comparisons"),
         (r"byproduct with the 7-cohort pool; E15_null\.csv is the reported value", "computed with the seven-cohort gene pool; the reported values for these cohorts are in sheet S7f_locked_null_dC"),
         (r"\bS4A\b", "authors' matrix"),
         (r"three new trials", "three added trials"), (r"3 new primary only", "3 added primary only"),
         (r"\bE1[0-9]\b", "")]


def clean(df):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == object:
            for a, b in CLEAN:
                df[c] = df[c].astype(str).str.replace(a, b, regex=True).where(df[c].notna(), None)
    return df


for key, items in NEW.items():
    n = num[key]
    for suffix, f, content in items:
        name = f"S{n}{suffix}"[:31]
        df = clean(pd.read_csv(R / f, sep="\t" if f.endswith(".tsv") else ","))
        sheets.append((n, name, df))
        readme.append((n, name, content, SRC[key]))

sheets.sort(key=lambda x: (x[0], x[1]))
readme.sort(key=lambda x: (x[0], x[1]))
rd = pd.DataFrame([r[1:] for r in readme], columns=["Sheet", "Content", "Source"])
rd["Content"] = rd["Content"].str.replace(r"\bE1[0-9]\b", "", regex=True)
out = G / "figures/v18/Supplementary_Tables_v18.xlsx"
with pd.ExcelWriter(out, engine="openpyxl") as w:
    for name, df in [("README", rd)] + [(s[1], s[2]) for s in sheets]:
        df.to_excel(w, sheet_name=name, index=False)
        ws = w.sheets[name]
        for cell in ws[1]:
            cell.font = Font(bold=True); cell.alignment = Alignment(wrap_text=True, vertical="top")
        for col in ws.columns:
            ws.column_dimensions[col[0].column_letter].width = min(60, max(12, max(len(str(c.value or "")) for c in col) * 0.9))
        ws.freeze_panes = "A2"
print("written", out, len(sheets), "sheets")
