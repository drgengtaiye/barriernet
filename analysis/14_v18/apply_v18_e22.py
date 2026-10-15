"""E22 (gastro-oesophageal checkpoint cohorts) into v18. Run after apply_v18_fixes.py:
apply_v18_edits.py -> apply_v18_fixes.py -> apply_v18_e22.py -> build_v18.py. Results: analysis/results/E22_README.md."""
import pathlib
from apply_v18_fixes import rep
HERE = pathlib.Path(__file__).resolve().parent


def main(s):
    # abstract (kept within 250 words)
    s = rep(s, "Checkpoint response was tested in a gastric pembrolizumab cohort and 21 non-gastric cohorts;",
            "Checkpoint response was tested in two gastric pembrolizumab cohorts and 21 non-gastric cohorts;")
    s = rep(s, "prognosis in 24 TCGA cancer types; and stromal geometry by imaging.", "prognosis in 24 TCGA cancer types; and stromal geometry.")
    s = rep(s, "In gastric cancer NetImmune predicted pembrolizumab response (AUC 0.80), not beyond microsatellite instability, Epstein–Barr virus or PD-L1 status.",
            "In gastric cancer NetImmune predicted pembrolizumab response in one cohort (AUC 0.80), not beyond microsatellite instability, Epstein–Barr virus or PD-L1 status, and not in a second (AUC 0.49).")
    s = rep(s, "CD8 exclusion tracked fibroblast number, not arrangement (underpowered test).", "CD8 exclusion tracked fibroblast number, not arrangement.")
    s = rep(s, "did not mark checkpoint responders (log-odds per SD −0.06, −0.21 to 0.09).", "did not mark checkpoint responders.")
    s = rep(s, "Gastric immune–stromal scores are usually validated by survival alone, which random gene sets also achieve.",
            "Gastric immune–stromal scores are validated by survival, which random gene sets also predict.")
    # methods
    s = rep(s, "- **Gastric pembrolizumab cohort:** PRJEB25780 (44 patients, 12 responders), with PD-L1 combined positive score (CPS) for 40[@kim2018].",
            "- **Gastric pembrolizumab cohorts:**\n"
            "  - PRJEB25780, pembrolizumab monotherapy (44 patients, 12 responders), with PD-L1 combined positive score (CPS) for 40[@kim2018];\n"
            "  - in the pre-registered extension, NCT04249739, first-line 5-fluorouracil and platinum followed by added pembrolizumab (43 patients with baseline RNA after quality control, 29 responders), processed from raw reads[@an2024; @wright2026].")
    s = rep(s, "  - Three chemotherapy–checkpoint cohorts, one of them in oesophageal adenocarcinoma (PERFECT), were used only as a sensitivity analysis[@hu2023; @thibaudin2023; @vandenende2021].",
            "  - Three chemotherapy–checkpoint cohorts, one of them in oesophageal adenocarcinoma (PERFECT), were used only as a sensitivity analysis[@hu2023; @thibaudin2023; @vandenende2021].\n"
            "  - An adjuvant durvalumab cohort in oesophageal and gastro-oesophageal junction adenocarcinoma, with relapse-free survival, was analysed separately[@mamdani2021].")
    s = rep(s, "independent single-cell atlases, PD-L1 CPS and robustness analyses were specified",
            "independent single-cell atlases, further gastro-oesophageal checkpoint cohorts, PD-L1 CPS and robustness analyses were specified")
    # results
    s = rep(s, "### In gastric cancer, fibroblast content tracks prognosis and immune programmes track response",
            "### In gastric cancer, fibroblast content tracks prognosis; the response association did not replicate")
    s = rep(s, "In the pembrolizumab cohort, NetImmune gave", "In the pembrolizumab monotherapy cohort, NetImmune gave")
    s = rep(s, "This cohort cannot separate the score from molecular subtype, so we turned to other cancers.",
            "This cohort cannot separate the score from molecular subtype.\n\n"
            "A second, pre-registered gastric cohort was treated first line with 5-fluorouracil and platinum followed by added pembrolizumab, and was mostly microsatellite-stable (43 patients, 29 responders; Supplementary Fig. {SF:gastric_icb}; Supplementary Table {ST:gastric_icb}).\n\n"
            "- **NetImmune** did not predict response (AUC 0.49; log-odds per SD −0.06, −0.70 to 0.59), in all patients or in the 40 MSS/EBV-negative patients (−0.01).\n"
            "- **The Immune axis and IFN module** pointed towards response, without reaching significance (+0.40 and +0.43; *P* = 0.25 and 0.22). In an analysis not pre-specified, the Immune axis was associated with longer progression-free survival (HR per SD 0.63; *P* = 0.009).\n"
            "- **CAF_ECM** also pointed towards response (+0.54; *P* = 0.14), opposite to expectation.\n"
            "- **PD-L1 CPS** was again the strongest single predictor (AUC 0.73).\n"
            "- **Processing.** The conclusions were unchanged with the pre-specified processing in place of the batch-aware processing adopted after quality control (NetImmune +0.02).\n"
            "- **Adjuvant durvalumab.** In the oesophageal adenocarcinoma cohort (37 patients, 20 relapses), no score was associated with relapse-free survival.\n\n"
            "The gastro-oesophageal checkpoint data therefore do not show that the score predicts response beyond MSI/EBV status and CPS, so we turned to other cancers.")
    # discussion and limitations
    s = rep(s, "In gastric cancer, PD-L1 CPS and MSI/EBV status had higher AUCs than NetImmune, which was not shown to add to them; the score's value lies less in selecting patients than in showing which component of the microenvironment carries which information.",
            "In gastric cancer the score did not hold up as a response marker. PD-L1 CPS and MSI/EBV status had higher AUCs, and NetImmune did not add to them. In a second, mostly microsatellite-stable cohort it did not predict response at all, although chemotherapy may account for much of the response in that cohort. The score's value lies less in selecting patients than in showing which component of the microenvironment carries which information.")
    i = s.index("- **Checkpoint cohorts.** Gastro-oesophageal"); j = s.index("\n", i)
    s = s[:i] + ("- **Checkpoint cohorts.** Gastro-oesophageal checkpoint data comprise two small, single-arm gastric cohorts whose results disagreed, and two small oesophageal cohorts. In the cohort treated with chemotherapy and pembrolizumab, response cannot be attributed to the checkpoint inhibitor.") + s[j:]
    # data availability
    s = rep(s, "Gastric pembrolizumab cohort: ENA PRJEB25780.",
            "Gastric pembrolizumab cohorts: ENA PRJEB25780 and PRJEB60680 (clinical data from the original publication). Adjuvant durvalumab cohort: GEO GSE183924.")
    return s


def supp(s):
    s = rep(s, "- Pembrolizumab cohort: PRJEB25780, a phase 2 trial of monotherapy in metastatic gastric cancer[@kim2018].",
            "- Pembrolizumab cohort: PRJEB25780, a phase 2 trial of monotherapy in metastatic gastric cancer[@kim2018]. A second cohort was added in the pre-registered extension (see Pre-registered validation analyses).")
    s = rep(s, "  Cell types were the authors' major annotations, mapped to 11 classes before scoring.",
            "  Cell types were the authors' major annotations, mapped to 11 classes before scoring.\n"
            "- **Gastro-oesophageal checkpoint cohorts.**\n"
            "  - **NCT04249739 (ENA PRJEB60680)**[@an2024; @wright2026]: first-line 5-fluorouracil and platinum for one cycle, then pembrolizumab added; 46 baseline biopsies taken before any treatment. Response was CR/PR against SD/PD; clinical data, MSI/EBV status and CPS are from the original publications.\n"
            "  - **Raw reads** (199 GB; md5-checked) were trimmed with fastp[@chen2018fastp] and quantified with salmon[@patro2017] against Ensembl 110 cDNA (k = 31; median mapping 76%).\n"
            "  - **Library chemistries.** Two were present (TruSeq-type, 18 samples; Nextera-type, 27, with a large mitochondrial fraction). Within-batch correlations were r = 0.81–0.87, against 0.52 between batches.\n"
            "  - **Processing for the primary analysis.** TPM were renormalised over nuclear protein-coding and immunoglobulin/T-cell-receptor genes, and gene z-scores were computed within batch.\n"
            "  - **Exclusions.** Three samples were excluded on identity or quality: mate files from different libraries, a sex mismatch, and a squamous-mucosa profile. 43 patients were analysed.\n"
            "  - **Timing of these decisions.** They were written down and hash-logged after technical quality control and before the full-cohort analysis. While the pipeline was built, test runs printed estimates for 7 and then 11 patients of one batch, and for the second cohort below; no decision used them.\n"
            "  - **Sensitivity analysis.** The pre-specified processing (all-transcript TPM, z-scores over the whole cohort) is reported alongside.\n"
            "  - **GSE183924**[@mamdani2021]: adjuvant durvalumab after chemoradiotherapy and surgery in oesophageal and gastro-oesophageal junction adenocarcinoma; 37 patients, 20 relapse events. FPKM were transformed to log₂ once. LTB is zero and TRAC and TRBC1 are absent in the deposited matrix, so TLS_B and T_CELL used five genes each. Cox models of relapse-free survival were reported separately and not pooled with the response cohorts.")
    s = rep(s, "  - a check of the CheckMate expression source.",
            "  - a check of the CheckMate expression source.\n"
            "- **Supplementary Table {ST:gastric_icb}.** Gastro-oesophageal checkpoint cohorts:\n"
            "  - NCT04249739 quality control;\n"
            "  - response log-odds, AUCs and random-gene null;\n"
            "  - progression-free and overall survival;\n"
            "  - Firth models with MSI/EBV status and CPS;\n"
            "  - pooled estimates with PRJEB25780 and PERFECT;\n"
            "  - GSE183924 relapse-free survival;\n"
            "  - pre-specified expectations.")
    k = s.index("## Supplementary references")
    k = s.rindex("\n---", 0, k)
    s = s[:k].rstrip() + "\n\n" + LEGEND + s[k:]
    return s


LEGEND = """![](FIG:gastric_icb)

**Supplementary Fig. {SF:gastric_icb}. Gastro-oesophageal checkpoint cohorts.**
(a) Log-odds of response per SD in two gastric cohorts, and pooled. The cohorts are NCT04249739 (PRJEB60680; 5-fluorouracil and platinum followed by added pembrolizumab; 43 patients, 29 CR/PR) and PRJEB25780 (pembrolizumab; 44 patients, 12 CR/PR). Pooled estimates use random effects with the conservative interval, without (filled) and with (open) PERFECT; ringed points, *P* < 0.05.
(b) ROC curves in NCT04249739 for NetImmune, IFN and CAF_ECM, with AUC (DeLong 95% CI) and the structure-matched null for NetImmune.
(c) Log-odds of response in the 40 MSS/EBV-negative patients (filled) and in all patients (open), with the penalised likelihood-ratio test for NetImmune beyond CPS.
(d) Relapse-free survival HR per SD in GSE183924 (adjuvant durvalumab, oesophageal and gastro-oesophageal junction adenocarcinoma; 37 patients, 20 events).
"""

if __name__ == "__main__":
    m = main((HERE / "v18_src.md").read_text())
    u = supp((HERE / "v18_supp_src.md").read_text())
    (HERE / "v18_src.md").write_text(m)
    (HERE / "v18_supp_src.md").write_text(u)
    print("E22 edits applied")
