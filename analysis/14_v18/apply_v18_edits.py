"""v18 text changes: the pre-registered validation analyses E17–E21 and E18 (analysis/E17_E22_prereg.md).
E22 (gastric checkpoint cohorts) is added separately once it has run. Each replacement must match once.
Backups: v18_src.pre_edits.md, v18_supp_src.pre_edits.md."""
import pathlib, shutil
HERE = pathlib.Path(__file__).resolve().parent

MAIN = [
    # ---- title, abstract, introduction ----
    ("21 checkpoint cohorts and three randomised comparisons\n", "21 checkpoint cohorts and six randomised comparisons\n"),
    ("**Methods.** The 48 genes were fixed before external validation. We benchmarked the axes in three gastric cohorts (1,097 patients) against random-gene nulls, established scores and fitted models, then tested them, locked, in three more (679 patients). Checkpoint response was tested in a gastric pembrolizumab cohort and 21 non-gastric cohorts; treatment interaction in three randomised comparisons; prognosis in 24 TCGA cancer types; and stromal geometry by imaging.",
     "**Methods.** The 48 genes were fixed before external validation. We benchmarked the axes in three gastric cohorts (1,097 patients) against random-gene nulls, established scores and fitted models, then tested them, locked, in seven more (1,055 patients). Checkpoint response was tested in a gastric pembrolizumab cohort and 21 non-gastric cohorts; treatment interaction in six randomised comparisons; prognosis in 24 TCGA cancer types; and stromal geometry by imaging. Later validations were pre-registered."),
    ("**Results.** NetImmune was prognostic in each benchmark cohort (pooled HR per SD 0.74, 95% CI 0.58–0.95), beyond stage and random axes, and pointed likewise in the locked cohorts (0.73, 0.53–1.02); it matched ESTIMATE and TMEscore, trailing the best fitted model by 0.024 in C-index. The fibroblast/matrix module, the strongest prognostic component (HR 1.29, 1.18–1.41), did not mark checkpoint responders (log-odds per SD −0.06, −0.21 to 0.09). Immune modules did, positively in all 10 replication cohorts (Immune +0.50, 0.28–0.71). NetImmune favoured the checkpoint arm in all three comparisons (pooled *P* = 0.057); the fibroblast module in none. In gastric cancer NetImmune predicted pembrolizumab response (AUC 0.80), not beyond microsatellite instability/Epstein–Barr virus status. CD8 exclusion tracked fibroblast number; arrangement had no detectable effect (underpowered test).",
     "**Results.** NetImmune was prognostic in the benchmark cohorts (pooled HR per SD 0.74, 95% CI 0.58–0.95) and the locked cohorts (0.70, 0.62–0.80), beyond stage and random axes; it matched ESTIMATE and TMEscore. The fibroblast/matrix module, the strongest prognostic component (HR 1.29, 1.18–1.41), did not mark checkpoint responders (log-odds per SD −0.06, −0.21 to 0.09). Immune modules did, positively in all 10 replication cohorts (Immune +0.50, 0.28–0.71). NetImmune favoured the checkpoint arm in all six randomised comparisons (pooled interaction +0.29, 0.03–0.55; *P* = 0.034); the fibroblast module in none. In gastric cancer NetImmune predicted pembrolizumab response (AUC 0.80), but not beyond microsatellite instability, Epstein–Barr virus or PD-L1 status. CD8 exclusion tracked fibroblast number; arrangement had no detectable effect (underpowered test)."),
    ("**Conclusions.** Fibroblast content tracked survival but not response; immune programmes tracked response and leaned towards checkpoint benefit. Predictive scores need validation against treatment interaction, not survival.",
     "**Conclusions.** Fibroblast content tracked survival but not response; immune programmes tracked response and, modestly, checkpoint benefit. Predictive scores need validation against treatment interaction, not survival."),
    ("judged in 21 checkpoint cohorts, three randomised comparisons and untreated tumours", "judged in 21 checkpoint cohorts, six randomised comparisons and untreated tumours"),
    # ---- data ----
    ("- **Gastric single-cell and spatial data:** GSE183904[@kumar2022] and 10 Visium sections from GSE251950[@lee2025].",
     "- **Gastric single-cell and spatial data:** GSE183904[@kumar2022] and 10 Visium sections from GSE251950[@lee2025]. Two further single-cell atlases, OMIX001073 and GSE206785[@sun2022; @kang2022], were used only to test cellular attribution."),
    ("- **Locked gastric cohorts,** analysed afterwards with every parameter fixed (679 patients): GSE15459 and GSE34942 with overall survival[@ooi2009; @lei2013], and GSE26253 with recurrence-free survival[@lee2014].",
     "- **Locked gastric cohorts,** analysed afterwards with every parameter fixed (1,055 patients):\n  - GSE15459 and GSE34942 with overall survival[@ooi2009; @lei2013], and GSE26253 with recurrence-free survival[@lee2014] (679 patients);\n  - in a pre-registered extension, GSE13861, GSE26899, GSE26901 and GSE29272 with overall survival[@cho2011; @oh2018; @wang2013] (376 patients after duplicate screening)."),
    ("- **Gastric pembrolizumab cohort:** PRJEB25780 (44 patients, 12 responders)[@kim2018].",
     "- **Gastric pembrolizumab cohort:** PRJEB25780 (44 patients, 12 responders), with PD-L1 combined positive score (CPS) for 40[@kim2018]."),
    ("  - Two I-SPY2 comparisons in HER2-negative breast cancer, with pathological complete response (pCR) as the endpoint[@wolf2022; @nanda2020; @pusztai2021]: paclitaxel with pembrolizumab (69 patients) against paclitaxel (179), and paclitaxel with durvalumab and olaparib (71) against paclitaxel (29).\n",
     "  - Two I-SPY2 comparisons in HER2-negative breast cancer, with pathological complete response (pCR) as the endpoint[@wolf2022; @nanda2020; @pusztai2021]: paclitaxel with pembrolizumab (69 patients) against paclitaxel (179), and paclitaxel with durvalumab and olaparib (71) against paclitaxel (29).\n"
     "  - In a pre-registered extension, every further randomised trial with open pre-treatment expression in both arms:\n"
     "    - CheckMate 025: nivolumab (120) against everolimus (130) in previously treated renal cell carcinoma, with progression-free survival[@motzer2015; @braun2020];\n"
     "    - JAVELIN Renal 101: avelumab plus axitinib (354) against sunitinib (372) in untreated renal cell carcinoma, with progression-free survival[@motzer2019; @motzer2020];\n"
     "    - NeoTRIP: carboplatin and nab-paclitaxel with (119) or without (122) atezolizumab in triple-negative breast cancer, with pCR[@gianni2022; @dugo2026].\n"),
    ("It was analysed in the same way as a sensitivity analysis[@padron2022].",
     "It was analysed in the same way as a sensitivity analysis[@padron2022], as was DUTRENEO (durvalumab and tremelimumab against cisplatin-based chemotherapy in bladder cancer, randomised only among immune-high tumours; 28 against 19)[@grande2026]."),
    ("The locked gastric cohorts, replication cohorts, PERFECT, I-SPY2, TCGA, METABRIC and IMC analyses were added later; for each, the expected direction of every primary score was written into the analysis code before it was run.",
     "The locked gastric cohorts, replication cohorts, PERFECT, I-SPY2, TCGA, METABRIC and IMC analyses were added later; for each, the expected direction of every primary score was written into the analysis code before it was run. The final validations (the locked-cohort extension, the further randomised trials, the independent single-cell atlases, PD-L1 CPS and robustness to gene loss) were specified in a dated, hash-logged protocol before any of their outcomes were examined; deviations are listed in the Supplementary Methods."),
    # ---- statistics ----
    ("  - Score × arm interaction per SD: Cox models for IMmotion150, and logistic models (odds ratios, ORs) for I-SPY2 adjusted for hormone-receptor status, which adaptive randomisation left unbalanced between arms.",
     "  - Score × arm interaction per SD: Cox models for IMmotion150, CheckMate 025 and JAVELIN Renal 101; logistic models (odds ratios, ORs) for I-SPY2, adjusted for hormone-receptor status, which adaptive randomisation left unbalanced between arms, and for NeoTRIP, adjusted for sequencing batch."),
    ("- **Fitted models.** Eight classes of full-transcriptome model, including penalised Cox regression[@zou2005; @simon2011], were trained leave-one-cohort-out.",
     "- **Fitted models.** Eight classes of full-transcriptome model, including penalised Cox regression[@zou2005; @simon2011], were trained leave-one-cohort-out.\n- **Robustness.** Scores were recomputed after removing 10%, 20% or 30% of the genes at random (1,000 draws each) and after removing each gene in turn."),
    ("For the analyses added later, the replication targets were stated in advance: the locked gastric HRs, the pooled response effects in the replication cohorts, the I-SPY2 interactions of NetImmune and CAF_ECM, the prognostic effect of CAF_ECM in the checkpoint cancer types and the IMC continuity coefficient.",
     "For the analyses added later, the replication targets were stated in advance: the locked gastric HRs (pooled across all seven locked cohorts in the extension), the pooled response effects in the replication cohorts, the I-SPY2 interactions of NetImmune and CAF_ECM, the pooled NetImmune interaction across all primary randomised comparisons, the prognostic effect of CAF_ECM in the checkpoint cancer types and the IMC continuity coefficient."),
    # ---- results: cell types ----
    ("but macrophages only 28% of TAM_SUP signal (Supplementary Fig. {SF:sc_share}).\n",
     "but macrophages only 28% of TAM_SUP signal (Supplementary Fig. {SF:sc_share}).\n\n"
     "In two independent gastric atlases (80,440 and 62,755 tumour cells), each module was again highest in its intended cell type (Supplementary Fig. {SF:sc_indep}; Supplementary Table {ST:sc_indep}). T and NK cells carried 97% and 93% of T_CELL signal and macrophages 15% and 39% of TAM_SUP signal. Fibroblasts carried 43% and 53% of CAF_ECM signal, below the pre-specified 50% in one atlas: fewer fibroblasts were captured than in the discovery atlas (1.8% and 4.9% of cells, against 10.9%), and pericytes and smooth-muscle cells carried about a quarter of the signal. CAF_ECM therefore marks perivascular as well as fibroblastic stroma.\n"),
    # ---- results: prognosis ----
    ("- **Upgrades to the score.** Separating cell abundance from cell state, learning the module weights, or adding a spatial arrangement score each changed C by less than 0.01 (Supplementary Fig. {SF:upgrades}).\n",
     "- **Upgrades to the score.** Separating cell abundance from cell state, learning the module weights, or adding a spatial arrangement score each changed C by less than 0.01 (Supplementary Fig. {SF:upgrades}).\n\n"
     "The score did not depend on particular genes. With 30% of the genes removed at random, NetImmune correlated with the full score at median ρ = 0.98, and every pooled estimate kept its direction in all 1,000 draws (Supplementary Fig. {SF:robust}; Supplementary Table {ST:robust}).\n"),
    ("In three further cohorts analysed with every parameter locked (Table 1; Supplementary Fig. {SF:locked}; Supplementary Table {ST:locked}):\n\n"
     "- **NetImmune.** The HR was below 1 in all three (0.65, 0.92 and 0.75). The pooled estimate was 0.73 (0.53–1.02; *P* = 0.057), not significant under the conservative rule, although it exceeded random axes (*P* ≤ 0.001).\n"
     "- **Adjustment.** After adjustment, NetImmune stayed protective in the two larger cohorts but not in GSE34942 (54 patients with stage; 1.38).\n"
     "- **CAF_ECM.** It was adverse in two of three (1.54 and 1.23; 0.95 in GSE34942; pooled 1.26, 0.74–2.14).\n"
     "- **Head to head.** In these cohorts NetImmune exceeded ESTIMATE in C (+0.036, 0.012–0.059).\n"
     "- **All six cohorts.** NetImmune gave 0.74 (0.67–0.82) and CAF_ECM 1.28 (1.17–1.41).\n",
     "In three further cohorts analysed with every parameter locked, NetImmune had an HR below 1 in each (0.65, 0.92 and 0.75), but the pooled estimate, 0.73 (0.53–1.02; *P* = 0.057), was not significant under the conservative rule. After this result we pre-registered an extension to four further cohorts with recoverable survival. Before any outcome model, 16 of their arrays were excluded by an expression-correlation duplicate screen; no excluded array matched a patient of another cohort on clinical data (Table 1; Supplementary Fig. {SF:locked}; Supplementary Table {ST:locked}).\n\n"
     "- **NetImmune.** The HR was below 1 in each added cohort (0.50–0.72; all *P* < 0.05). Across all seven locked cohorts it was 0.70 (0.62–0.80; *P* = 4.3 × 10⁻⁴; I² = 7%), beyond random axes (*P* ≤ 0.001). After adjustment for stage and age it was 0.75 (0.62–0.92).\n"
     "- **CAF_ECM.** The HR was above 1 in six of seven cohorts (pooled 1.34, 1.08–1.65).\n"
     "- **Head to head.** Across the seven cohorts NetImmune exceeded CD8A in C (+0.062, 0.032–0.091) but not ESTIMATE (+0.015, −0.002 to 0.031).\n"
     "- **All ten gastric cohorts.** NetImmune gave 0.72 (0.67–0.79) and CAF_ECM 1.30 (1.17–1.45). Estimates were unchanged without the duplicate exclusions, with overall survival only, or with a cohort whose survival data come from a third-party curation[@li2023].\n"),
    # ---- results: gastric pembrolizumab ----
    ("NetImmune added nothing significant to MSI/EBV status (Firth likelihood-ratio *P* = 0.097; Supplementary Table {ST:prjeb_combo}).",
     "NetImmune added nothing significant to MSI/EBV status (Firth likelihood-ratio *P* = 0.097; Supplementary Table {ST:prjeb_combo}). PD-L1 CPS, available for 40 patients, discriminated responders better than NetImmune (AUC 0.97 against 0.82 in the same patients); every responder with a CPS had CPS ≥ 1. NetImmune pointed the same way beyond MSI/EBV status and CPS together, but not significantly (OR per SD 4.6; *P* = 0.11)."),
    # ---- results: randomised comparisons ----
    ("### Randomised comparisons: immune scores lean towards the checkpoint arm, fibroblast content away from it\n",
     "### Randomised comparisons: immune scores favour the checkpoint arm, fibroblast content does not\n\nWe first analysed three randomised comparisons (Supplementary Fig. {SF:trials_detail}):\n\n"),
    ("(interaction HR per SD 0.68, 0.70 and 0.65; *P* = 0.035, 0.042 and 0.015). CAF_ECM did not (1.21, *P* = 0.25) (Fig. 5a; Table 3; Supplementary Fig. {SF:icb11}).",
     "(interaction HR per SD 0.68, 0.70 and 0.65; *P* = 0.035, 0.042 and 0.015). CAF_ECM did not (1.21, *P* = 0.25) (Table 3; Supplementary Fig. {SF:icb11})."),
    ("- **I-SPY2 pembrolizumab** (69 patients against 179 controls; Fig. 5b; Supplementary Table {ST:trials}):",
     "- **I-SPY2 pembrolizumab** (69 patients against 179 controls; Supplementary Table {ST:trials}):"),
    ("random-gene null *P* = 0.020; Fig. 5d).", "random-gene null *P* = 0.020)."),
    ("- **I-SPY2 durvalumab–olaparib** (71 against 29; Fig. 5c). No score modified the treatment effect.\n\n"
     "Across the three comparisons (Fig. 5a; Table 3):\n\n"
     "- **NetImmune** favoured the checkpoint arm in all three (pooled benefit-oriented log ratio +0.47, −0.04 to 0.97; *P* = 0.057, not significant; Stouffer *P* = 0.004).\n"
     "- **IFN** also favoured it in all three (*P* = 0.063).\n"
     "- **CAF_ECM** favoured it in none (−0.37, −1.16 to 0.42; *P* = 0.18; Stouffer *P* = 0.021, against benefit).\n\n"
     "With PRINCE added as a sensitivity analysis, NetImmune was +0.45 (0.16–0.73; *P* = 0.015). Among treated patients, fibroblast content did not distinguish responders; in the randomised comparisons it pointed against benefit each time, significantly in one.\n",
     "- **I-SPY2 durvalumab–olaparib** (71 against 29). No score modified the treatment effect.\n\n"
     "Pooled, NetImmune favoured the checkpoint arm in all three comparisons, but not significantly (benefit-oriented log ratio +0.47, −0.04 to 0.97; *P* = 0.057).\n\n"
     "In a pre-registered extension we added every further randomised trial with open expression data from both arms (Fig. 5; Tables 3 and 4):\n\n"
     "- **CheckMate 025** (nivolumab against everolimus; 120 against 130). Interactions were small: NetImmune HR 0.93 (0.70–1.23); CAF_ECM 1.02 (0.78–1.33).\n"
     "- **JAVELIN Renal 101** (avelumab plus axitinib against sunitinib; 354 against 372). NetImmune and the Immune axis leaned towards the checkpoint arm (HR 0.83, 0.67–1.03, *P* = 0.090; and 0.83, *P* = 0.081); CAF_ECM did not (1.06).\n"
     "- **NeoTRIP** (chemotherapy with or without atezolizumab; 119 against 122). NetImmune modified the benefit of atezolizumab (OR 1.90, 1.08–3.36; *P* = 0.027; random-gene null *P* = 0.019), and CAF_ECM did so in the opposite direction (0.58, 0.33–0.99; *P* = 0.048; null *P* = 0.001).\n\n"
     "Across the six primary comparisons (Table 4):\n\n"
     "- **NetImmune** favoured the checkpoint arm in all six (pooled +0.29, 0.03–0.55; *P* = 0.034; Stouffer *P* = 1.3 × 10⁻⁴). Leaving out any single trial, the estimate stayed positive (+0.23 to +0.37) but *P* exceeded 0.05 three times; with overall survival in place of progression-free survival in CheckMate 025 it was 0.059.\n"
     "- **The immune modules** also favoured it: Immune +0.20 (0.04–0.35; *P* = 0.021) and IFN +0.21 (0.04–0.39; *P* = 0.027).\n"
     "- **CAF_ECM** favoured it in none of the six (−0.17, −0.43 to 0.08; *P* = 0.14; Stouffer *P* = 0.0067, against benefit).\n"
     "- **By cancer type,** the effect was clearer in the three breast cancer comparisons (NetImmune +0.63, *P* = 0.048; CAF_ECM −0.60, *P* = 0.034) than in the three renal comparisons (+0.20, *P* = 0.15).\n\n"
     "With PRINCE and DUTRENEO added as sensitivity analyses, NetImmune was +0.27 (0.07–0.47; *P* = 0.014). Among treated patients, fibroblast content did not distinguish responders; across randomised comparisons it never favoured the checkpoint arm, and it pointed against benefit significantly in the two breast cancer trials.\n"),
    # ---- discussion ----
    ("In three locked cohorts it pointed the same way in each, although the pooled interval just included 1.",
     "In seven locked cohorts, four of them added in a pre-registered extension after the first three fell just short of significance, it was prognostic in each but one (pooled HR 0.70)."),
    ("Yet they marked responders in the first 11 cohorts and again in the 10 replication cohorts, and they leaned towards benefit in randomised comparisons.",
     "Yet they marked responders in the first 11 cohorts and again in the 10 replication cohorts, and they favoured the checkpoint arm across six randomised comparisons."),
    ("In the language of biomarker design[@ballman2015], immune engagement behaves as a predictive marker, although no pooled randomised interaction reached significance.",
     "In the language of biomarker design[@ballman2015], immune engagement behaves as a predictive marker. Its interaction with treatment was consistent in direction across six randomised comparisons and nominally significant when pooled, but it was small, clearest in breast cancer and sensitive to the removal of single trials."),
    ("In all three randomised comparisons fibroblast content pointed against benefit from the checkpoint arm, significantly in the I-SPY2 pembrolizumab comparison. That rests on one significant comparison, in breast cancer, with a pathological endpoint. We therefore read fibroblast content as prognostic, with a possible modest effect against checkpoint benefit that single-arm cohorts are not designed to detect.",
     "In none of six randomised comparisons did fibroblast content favour the checkpoint arm. It pointed against benefit significantly in the two breast cancer trials with pathological endpoints, I-SPY2 and NeoTRIP, and not in the renal trials. We therefore read fibroblast content as prognostic, with a possible effect against checkpoint benefit that single-arm cohorts are not designed to detect and that may depend on the setting."),
    ("A composite validated against survival can look strong while its predictive part is diluted.\n",
     "A composite validated against survival can look strong while its predictive part is diluted.\n\n"
     "In gastric cancer the score did not add to the clinical biomarkers. PD-L1 CPS and MSI/EBV status each discriminated pembrolizumab responders better than NetImmune. Its value lies less in selecting patients than in showing which component of the microenvironment carries which information.\n"),
    ("and the locked replication used only three cohorts, one small and one with a different endpoint, and was not significant under conservative pooling.",
     "The locked replication was extended from three to seven cohorts after the first three gave *P* = 0.057; the extension was pre-registered and every parameter stayed locked, but the decision to extend was taken in the light of that result."),
    ("- **Checkpoint cohorts.** Gastro-oesophageal checkpoint data are limited to one small gastric cohort dominated by MSI-H and EBV-positive responders, and one small oesophageal cohort treated with atezolizumab and chemoradiotherapy.",
     "- **Checkpoint cohorts.** Gastro-oesophageal checkpoint data are limited to one small gastric cohort dominated by MSI-H and EBV-positive responders, in which PD-L1 CPS outperformed the score, and one small oesophageal cohort treated with atezolizumab and chemoradiotherapy."),
    ("- **Randomised comparisons.** They differ in cancer type, comparator and endpoint. Their pooled interaction summarises direction, not effect size.",
     "- **Randomised comparisons.** They differ in cancer type, comparator and endpoint, and none is gastric. Their pooled interaction summarises direction, not effect size; it rests mainly on the breast cancer trials and loses significance when any of three trials is removed."),
    ("- **Cellular attribution.** It was assessed in the atlas used to refine the modules, so it confirms, but does not independently test, their intended cell sources.",
     "- **Cellular attribution.** Two independent atlases confirmed the intended cell types, but they captured few fibroblasts, and CAF_ECM is also expressed by perivascular cells."),
    ("A randomised gastric checkpoint cohort powered for the microsatellite-stable, EBV-negative majority would test the immune modules and CXCL9:SPP1 polarity where they would matter.",
     "A randomised gastric checkpoint trial with tumour expression, such as CheckMate 649, would test the immune modules and CXCL9:SPP1 polarity where they would matter, in the microsatellite-stable, EBV-negative majority."),
    ("whereas immune programmes track checkpoint response across cancers and lean towards benefit in randomised comparisons.",
     "whereas immune programmes track checkpoint response across cancers and modestly favour checkpoint benefit in randomised comparisons."),
    # ---- data availability ----
    ("locked cohorts GEO GSE15459, GSE34942 and GSE26253.",
     "locked cohorts GEO GSE15459, GSE34942, GSE26253, GSE13861, GSE26899, GSE26901 and GSE29272, with survival for the last four from the original publications. Further single-cell atlases: NGDC OMIX001073 and GEO GSE206785."),
    ("I-SPY2: GEO GSE194040 and GSE173839.",
     "I-SPY2: GEO GSE194040 and GSE173839. CheckMate 025 and JAVELIN Renal 101: supplementary data of the cited biomarker studies. NeoTRIP: GEO GSE319641 with the clinical annotation released by the authors. DUTRENEO: GEO GSE328930."),
]

SUPP = [
    ("Three further gastric cohorts with deposited outcomes were analysed afterwards with every parameter locked:",
     "Seven further gastric cohorts were analysed afterwards with every parameter locked. The first three had deposited outcomes:"),
    ("(Supplementary Table {ST:locked}).\n\nExternal checkpoint data",
     "(Supplementary Table {ST:locked}).\n\n"
     "The other four were added in a pre-registered extension (see Pre-registered validation analyses), with survival recovered from the original publications:\n\n"
     "- **GSE13861** (Yonsei) and **GSE26899** and **GSE26901** (Korea University Guro and Kosin), with overall survival from Oh et al.[@cho2011; @oh2018];\n"
     "- **GSE29272** (Shanxi, China), with overall survival from Wang et al.[@wang2013].\n\n"
     "GSE57303, whose survival data are available only from a third-party curation[@li2023], was a sensitivity analysis. GSE51105 has no public per-sample outcome.\n\n"
     "Before any outcome model, each sample of the added cohorts was correlated with every sample of the cohorts already used. Samples were excluded when their maximum cross-cohort correlation exceeded the cohort's within-cohort 99.9th percentile; this removed 16 arrays from the added cohorts and 4 from GSE57303. None of the flagged pairs matched on age, sex and survival, and results without the exclusions are also reported.\n\n"
     "External checkpoint data"),
    ("Prognosis outside gastric cancer was assessed",
     "Three further randomised trials with open pre-treatment expression in both arms were added in the pre-registered extension:\n\n"
     "- **CheckMate 025**[@motzer2015; @braun2020]: nivolumab (120 patients with RNA) against everolimus (130) in previously treated renal cell carcinoma; progression-free survival, with overall survival (the trial's primary endpoint) as a secondary endpoint. The authors' log₂ expression matrix was used.\n"
     "- **JAVELIN Renal 101**[@motzer2019; @motzer2020]: avelumab plus axitinib (354) against sunitinib (372) in untreated renal cell carcinoma; progression-free survival.\n"
     "- **NeoTRIP**[@gianni2022; @dugo2026]: carboplatin and nab-paclitaxel with (119) or without (122) atezolizumab in triple-negative breast cancer; pCR. Arm and pCR were taken from the clinical annotation released by the authors.\n\n"
     "DUTRENEO[@grande2026] (durvalumab and tremelimumab, 28, against cisplatin-based chemotherapy, 19, in bladder cancer) was a sensitivity analysis only. Randomisation was restricted to immune-high tumours, and the checkpoint group includes non-randomised run-in patients who cannot be identified in the data.\n\n"
     "Prognosis outside gastric cancer was assessed"),
    ("- **Pooling across trials.** Interaction effects were oriented so that positive values mean more benefit from the checkpoint-containing arm: −log interaction HR for IMmotion150 and +log interaction OR for I-SPY2.",
     "- **Added trials.** Cox models of progression-free survival for CheckMate 025 and JAVELIN Renal 101, and logistic models of pCR for NeoTRIP with sequencing batch as a covariate (batches were uneven across arms; results without it are also reported). Within-arm effects were reported for every trial.\n"
     "- **Pooling across trials.** Interaction effects were oriented so that positive values mean more benefit from the checkpoint-containing arm: −log interaction HR for time-to-event endpoints and +log interaction OR for pCR."),
    ("### Analysis hierarchy (as in the main text)",
     "### Pre-registered validation analyses\n\n"
     "The final validation analyses were specified in a protocol written before any of their data were related to outcome. It was hash-logged at creation and at every dated amendment, and is released with the code. The protocol covered:\n\n"
     "- cellular attribution in independent single-cell atlases;\n"
     "- PD-L1 CPS in PRJEB25780;\n"
     "- robustness to gene loss;\n"
     "- the extended locked gastric replication;\n"
     "- further randomised comparisons;\n"
     "- further gastro-oesophageal checkpoint cohorts.\n\n"
     "Amendments were made before the relevant outcomes were examined, as data sources were identified:\n\n"
     "- the choice of atlases and cohorts;\n"
     "- the CPS transformation;\n"
     "- the trial endpoints and the batch covariate for NeoTRIP.\n\n"
     "Deviations recorded after analysis:\n\n"
     "- **Robustness.** On the 22 genes shared by all four panel cohorts, CAF_ECM keeps only one gene, so each panel platform's own gene content was scored instead.\n"
     "- **Cellular attribution.** The authors' published cell-type annotations were used, in place of marker-based re-annotation.\n"
     "- **CheckMate 025 expression.** The compendium copy of the CheckMate expression data used for the first 11 cohorts agrees with the authors' matrix for NetImmune (ρ = 0.99) but less for single modules (ρ 0.85–0.95). Rescoring that cohort from the authors' matrix changed no pooled response estimate by more than 0.015.\n\n"
     "- **Cellular attribution.** Cells were scored and attributed as in the discovery atlas: the mean log-normalised expression of the module genes, class means standardised across classes, and the abundance-weighted share of signal. Two atlases were used:\n"
     "  - OMIX001073[@sun2022]: 10 patients, 80,440 tumour cells;\n"
     "  - GSE206785[@kang2022]: 24 patients, 62,755 tumour cells, depleted of epithelial cells.\n\n"
     "  Cell types were the authors' major annotations, mapped to 11 classes before scoring.\n"
     "- **PD-L1 CPS.** CPS (22C3) was transcribed from the original report[@kim2018] and matched to 40 of the 44 patients. It was entered as log(CPS + 1); a CPS ≥ 1 cut-off was a sensitivity analysis. NetImmune beyond MSI/EBV status and CPS was tested by penalised likelihood ratio in Firth logistic regression.\n"
     "- **Robustness.** Draws dropped 10%, 20% or 30% of the 48 genes, keeping at least three genes per scorable module. The same genes were dropped in every cohort of a draw.\n\n"
     "### Analysis hierarchy (as in the main text)"),
    ("- the pooled NetImmune and CAF_ECM hazard ratios in the three locked gastric cohorts;",
     "- the pooled NetImmune and CAF_ECM hazard ratios in the three locked gastric cohorts, and, in the pre-registered extension, across all seven;"),
    ("- the NetImmune and CAF_ECM interactions in the two I-SPY2 comparisons;",
     "- the NetImmune and CAF_ECM interactions in the two I-SPY2 comparisons, and the pooled NetImmune interaction across all primary randomised comparisons;"),
    # S9 (prjeb_combo) note
    ("Penalised likelihood-ratio test for NetImmune beyond MSI/EBV status: χ² = 2.75, 1 d.f., *P* = 0.097; OR 3.0 per SD (95% CI 0.76–11.7). LOOCV, leave-one-out cross-validation; CV, cross-validation.",
     "Penalised likelihood-ratio test for NetImmune beyond MSI/EBV status: χ² = 2.75, 1 d.f., *P* = 0.097; OR 3.0 per SD (95% CI 0.76–11.7).\n\n"
     "With PD-L1 CPS (40 patients with CPS; 11 responders):\n\n"
     "- **AUC:** log(CPS + 1) 0.970 (0.912–1.000), CPS ≥ 1 0.828, MSI/EBV status 0.937, NetImmune 0.825 (0.661–0.953).\n"
     "- **NetImmune beyond MSI/EBV status and log(CPS + 1):** OR 4.6 per SD, χ² = 2.58, *P* = 0.11.\n"
     "- **The same model with CPS ≥ 1:** OR 3.2, *P* = 0.18.\n\n"
     "LOOCV, leave-one-out cross-validation; CV, cross-validation."),
    # workbook guide
    ("- **Supplementary Table {ST:locked}.** Locked gastric cohorts: cohort table, overlap check, module coverage, per-cohort and pooled hazard ratios, random-gene null and ΔC.",
     "- **Supplementary Table {ST:locked}.** Locked gastric cohorts (first three and the pre-registered extension):\n"
     "  - cohort table;\n"
     "  - overlap and duplicate checks;\n"
     "  - module coverage;\n"
     "  - per-cohort and pooled hazard ratios, with sensitivity sets;\n"
     "  - random-gene null;\n"
     "  - ΔC.\n"
     "- **Supplementary Table {ST:sc_indep}.** Cellular attribution in two independent single-cell atlases:\n"
     "  - cell-class mapping;\n"
     "  - module scores by class;\n"
     "  - share of signal, pooled and per patient;\n"
     "  - pre-specified expectations;\n"
     "  - sensitivity analyses.\n"
     "- **Supplementary Table {ST:robust}.** Robustness to gene loss:\n"
     "  - random-dropout summaries;\n"
     "  - leave-one-gene-out estimates;\n"
     "  - panel-platform scores."),
    ("- **Supplementary Table {ST:trials}.** I-SPY2 interactions (unadjusted, adjusted, stratified), within-arm effects, random-gene null and pooled cross-trial estimates.",
     "- **Supplementary Table {ST:trials}.** Randomised comparisons:\n"
     "  - I-SPY2 interactions (unadjusted, adjusted, stratified);\n"
     "  - CheckMate 025, JAVELIN Renal 101, NeoTRIP and DUTRENEO interactions with within-arm effects;\n"
     "  - random-gene nulls;\n"
     "  - pooled estimates for every set, including leave-one-trial-out."),
]

LEGENDS = """
![](FIG:sc_indep)

**Supplementary Fig. {SF:sc_indep}. Cellular attribution in two independent gastric single-cell atlases.** (a) Module score by cell class (z across classes) in the discovery atlas (GSE183904) and in OMIX001073 and GSE206785; boxes mark each module's top class. (b) Share of each module's signal by cell class (cell fraction × class mean score), with per-patient medians (IQR); dashed line, 50%. (c) Expression of the module genes by class in OMIX001073; boxes mark the pre-specified intended class.

![](FIG:robust)

**Supplementary Fig. {SF:robust}. Robustness to gene loss.** (a) Spearman correlation of NetImmune with the full score as 10–30% of genes are removed at random (1,000 draws each; grey, per cohort; blue, median and 5th percentile over all cohorts). (b) Pooled NetImmune HR in the benchmark cohorts, NetImmune AUC in PRJEB25780, and pooled Immune and CAF_ECM log-odds of response across 21 checkpoint cohorts, with 2.5–97.5th percentiles over draws; filled, draws keeping at least three genes per module in every cohort; open, draws constrained only on the 48-gene list; squares, scores restricted to each panel platform's genes (n.e., not estimable); band, 95% CI of the full score. (c) Each pooled estimate after removing each gene in turn; open circles, removal leaves a module with fewer than three genes in at least one cohort.

![](FIG:trials_detail)

**Supplementary Fig. {SF:trials_detail}. The first three randomised comparisons.** (a) Interaction estimates per SD in IMmotion150 (progression-free survival), the I-SPY2 pembrolizumab and durvalumab–olaparib comparisons (pCR, adjusted for hormone-receptor status) and PRINCE (progression-free survival; sensitivity). Every axis is oriented so that more benefit from the checkpoint-containing arm lies to the right; ringed points, *P* < 0.05. The right panel pools the three comparisons (filled diamonds) and adds PRINCE (open diamonds). (b,c) Within-arm OR for pCR per SD in the two I-SPY2 comparisons, with the interaction OR and *P*. (d) Absolute interaction log-OR of 2,000 structure-matched random gene sets for NetImmune and CAF_ECM in I-SPY2, with the observed value and empirical *P*.
"""
TABLES34 = """**Table 3. Score × treatment interaction in each randomised comparison.**

{T3}

- **Interaction** per SD of the score: HR for progression-free survival (below 1 favours the checkpoint arm) or OR for pCR (above 1 favours the checkpoint arm).
- **Adjustment.** I-SPY2 is adjusted for hormone-receptor status and NeoTRIP for sequencing batch.
- **Sensitivity comparisons.** PRINCE and DUTRENEO are sensitivity analyses: in PRINCE the control arm also received an immunotherapy (sotigalimab), and DUTRENEO randomised only immune-high tumours.
- **Further results.** All scores and within-arm effects are in Supplementary Table {ST:trials}.

**Table 4. Pooled score × treatment interaction across randomised comparisons.**

{T4}

- **Pooled interaction:** benefit-oriented log ratio per SD (−log HR, +log OR), random effects, with the wider of the Hartung–Knapp and DerSimonian–Laird intervals and the larger *P*.
- **Favouring:** the number of primary comparisons in which the score favoured the checkpoint arm.
- **Stouffer *P*:** two-sided, from the combined *z*.
- **Robustness.** Leave-one-trial-out and cancer-type subsets are in Supplementary Table {ST:trials}.
"""
FIG5 = ("**Figure 5. Score × treatment interaction in six randomised comparisons.** (a) Benefit-oriented interaction per SD (−log HR for progression-free survival; log OR for pCR) in IMmotion150, CheckMate 025 and JAVELIN Renal 101 (renal cell carcinoma) and in the I-SPY2 pembrolizumab, I-SPY2 durvalumab–olaparib and NeoTRIP comparisons (breast cancer). Right of zero, more benefit from the checkpoint-containing arm; ringed points, *P* < 0.05. "
        "The pooled column shows random-effects estimates across the six primary comparisons (filled diamonds) and with PRINCE and DUTRENEO added (open diamonds), with the number favouring the checkpoint arm and the conservative *P*. (b) PRINCE and DUTRENEO (sensitivity). "
        "(c) Robustness of the pooled NetImmune and CAF_ECM interactions: all six, leaving out each trial, by cancer type, with the sensitivity trials, with overall survival for CheckMate 025 and without the NeoTRIP batch term. (d) Absolute interaction of 2,000 structure-matched random gene sets for NetImmune and CAF_ECM in the three added trials, with the observed value and empirical *P*; '(against)' marks an observed effect favouring the control arm.")
LOCKED_OLD_START = "**Supplementary Fig. {SF:locked}. "
LOCKED_NEW = ("**Supplementary Fig. {SF:locked}. Locked replication in seven gastric cohorts.** (a) HR per SD for NetImmune and CAF_ECM in the benchmark cohorts (circles), the first three locked cohorts and the four added in the pre-registered extension (squares; GSE26253 recurrence-free survival), a sensitivity cohort with third-party survival data (open square) and pooled estimates (diamonds; wider of the Hartung–Knapp and DerSimonian–Laird intervals). "
              "(b) Absolute log HR of 1,000 structure-matched random axes per added cohort and pooled, with NetImmune and empirical *P*. (c) Paired difference in Harrell's C, NetImmune minus ESTIMATE or CD8A, per added cohort and pooled.")


def apply(fn, reps):
    p = HERE / fn
    s = p.read_text()
    for a, b in reps:
        n = s.count(a)
        assert n == 1, (fn, n, a[:90])
        s = s.replace(a, b)
    return s


if __name__ == "__main__":
    for f in ("v18_src.md", "v18_supp_src.md"):
        bak = HERE / f.replace(".md", ".pre_edits.md")
        if not bak.exists(): shutil.copy(HERE / f, bak)
    m = apply("v18_src.pre_edits.md", MAIN)
    s = apply("v18_supp_src.pre_edits.md", SUPP)
    # replace the locked-cohort legend and append the new legends before the references
    i = s.index(LOCKED_OLD_START); j = s.index("\n", i)
    s = s[:i] + LOCKED_NEW + s[j:]
    k = s.index("## Supplementary references")
    k = s.rindex("\n---", 0, k)                        # keep the section rule after the last legend
    s = s[:k].rstrip() + "\n" + LEGENDS + s[k:]
    # tables generated from the analysis outputs
    import make_v18_tables as T
    old1 = ("| Pooled | 679 (299) | 0.73 (0.53–1.02) | 0.057 | 0.85 (0.38–1.87) | 0.46 | 1 (*P* ≤ 0.001) | — |\n"
            "| **All six cohorts, pooled** | 1,776 (797) | 0.74 (0.67–0.82) | 6.6 × 10⁻⁴ | — | — | — | — |")
    assert m.count(old1) == 1
    m = m.replace(old1, old1.split("\n")[0].replace("| Pooled |", "| Pooled, first three |") + "\n" + T.table1_rows())
    a = "| **Locked cohorts** | | | | | | | |"; assert m.count(a) == 1
    m = m.replace(a, "| **Locked cohorts, first three** | | | | | | | |")
    a = "- **Endpoint.** GSE26253 reports recurrence-free survival."; assert m.count(a) == 1
    m = m.replace(a, a + "\n- **Pre-registered extension.** The four added cohorts exclude 16 arrays flagged by the pre-specified duplicate screen. The extension was decided after the first three cohorts gave *P* = 0.057.")
    i = m.index("**Table 3. Score × treatment interaction in randomised comparisons.**")
    j = m.index("\n---", i)
    m = m[:i] + TABLES34.replace("{T3}", T.table3()).replace("{T4}", T.table4()) + m[j:]
    i = m.index("**Figure 5. Score × treatment interaction in randomised comparisons.**")
    j = m.index("\n", i)
    m = m[:i] + FIG5 + m[j:]
    (HERE / "v18_src.md").write_text(m)
    (HERE / "v18_supp_src.md").write_text(s)
    print("v18 sources written:", len(MAIN), "main and", len(SUPP), "supplement replacements")
