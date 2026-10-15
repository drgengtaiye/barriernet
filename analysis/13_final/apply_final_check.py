"""Apply the text fixes from analysis/results/FINAL_CHECK.md to lean_src.md and lean_supp_src.md
(2026-10-03). Backups: lean_src.pre_final.md, lean_supp_src.pre_final.md. Each replacement must match once."""
import pathlib
HERE = pathlib.Path(__file__).resolve().parent

MAIN = [
    # A7, B1: abstract
    ("*barriernet* scores six marker modules and two axes: **Barrier** (fibroblast/matrix and suppressive-macrophage modules) and **NetImmune** (four immune modules minus Barrier).",
     "*barriernet* scores six marker modules and three axes: **Barrier** (fibroblast/matrix and suppressive-macrophage modules), **Immune** (four immune modules) and **NetImmune** (Immune minus Barrier)."),
    ("**Methods.** The 48 genes were prespecified.", "**Methods.** The 48 genes were fixed before external validation."),
    ("a gastric pembrolizumab cohort and 21 cohorts from other cancers (10 new);", "a gastric pembrolizumab cohort and 21 cohorts from other cancers;"),
    ("(pooled HR per SD 0.74, 95% CI 0.58–0.95)", "(pooled hazard ratio [HR] per SD 0.74, 95% CI 0.58–0.95)"),
    ("trailed the best fitted model by 0.024 in C. The fibroblast/ECM module", "trailed the best fitted model by 0.024 in C-index. The fibroblast/matrix module"),
    ("Immune modules did, replicating in all 10 new cohorts (Immune +0.50, 0.28–0.71). NetImmune favoured the checkpoint arm in all three randomised comparisons (pooled *P* = 0.057)",
     "Immune modules did, positively in all 10 replication cohorts (Immune +0.50, 0.28–0.71). NetImmune favoured the checkpoint arm in all three comparisons (pooled *P* = 0.057)"),
    ("not beyond MSI/EBV. CD8 exclusion", "not beyond microsatellite instability/Epstein–Barr virus status. CD8 exclusion"),
    # B2, B3: introduction
    ("strongest predictors remain microsatellite instability and Epstein–Barr virus positivity[@kim2018]",
     "strongest predictors remain microsatellite instability (MSI) and Epstein–Barr virus (EBV) positivity[@kim2018]"),
    ("the ACRG molecular subtypes[@cristescu2015]", "the Asian Cancer Research Group (ACRG) molecular subtypes[@cristescu2015]"),
    ("six marker modules and two axes that can be scored unchanged", "six marker modules and three axes that can be scored unchanged"),
    # B3, B16: data and methods
    ("- **Gastric bulk cohorts with overall survival:** TCGA-STAD[@tcga2014]",
     "- **Gastric bulk cohorts with overall survival:** TCGA-STAD (The Cancer Genome Atlas stomach adenocarcinoma)[@tcga2014]"),
    ("  - The PRINCE trial, analysed in the same way as a sensitivity analysis[@padron2022].",
     "  - The PRINCE trial in metastatic pancreatic cancer: nivolumab-containing arms against sotigalimab without nivolumab, all with gemcitabine and nab-paclitaxel, with progression-free survival as the endpoint. It was analysed in the same way as a sensitivity analysis[@padron2022]."),
    ("- **CAF_ECM**: COL1A1", "- **CAF_ECM** (cancer-associated fibroblasts and extracellular matrix, ECM): COL1A1"),
    ("- **TAM_SUP**: SPP1", "- **TAM_SUP** (suppressive tumour-associated macrophages): SPP1"),
    ("- **APM**: B2M", "- **APM** (antigen-presentation machinery): B2M"),
    ("- **IFN**: STAT1", "- **IFN** (interferon-γ response): STAT1"),
    ("- **TLS_B**: CXCL13", "- **TLS_B** (TLS and B cells): CXCL13"),
    ("- **T_CELL**: CD3D", "- **T_CELL** (T cells and cytotoxicity): CD3D"),
    ("standardised to unit variance, so effects are per SD.", "standardised to unit variance, so effects are per standard deviation (SD)."),
    ("- **Survival.** Cox models, unadjusted and adjusted for age and stage.",
     "- **Survival.** Cox models giving hazard ratios (HRs), unadjusted and adjusted for age and stage."),
    ("In PRJEB25780 we used AUCs with DeLong tests[@delong1988] and Benjamini–Hochberg correction[@bh1995] against GEP-18[@ayers2017], PD-L1 and IPS[@charoentong2017].",
     "In PRJEB25780 we used areas under the receiver operating characteristic curve (AUCs) with DeLong tests[@delong1988] and Benjamini–Hochberg correction[@bh1995] against the 18-gene T-cell-inflamed expression profile (GEP-18)[@ayers2017], PD-L1 and the immunophenoscore (IPS)[@charoentong2017]."),
    ("and logistic models for I-SPY2 adjusted for hormone-receptor status", "and logistic models (odds ratios, ORs) for I-SPY2 adjusted for hormone-receptor status"),
    # B8
    ("- **NetImmune versus CAF_ECM.** NetImmune contained all of CAF_ECM's prognostic information: adding CAF_ECM never improved fit (*P* ≥ 0.076).",
     "- **NetImmune versus CAF_ECM.** Adding CAF_ECM to NetImmune did not significantly improve fit in any cohort (*P* ≥ 0.076)."),
    # B6
    ("- **CAF_ECM** was adverse in all five, but weakly: the adjusted HR was 1.07 (0.99–1.15), not beyond random gene sets.",
     "- **CAF_ECM** had an HR above 1 in all five, but weakly (adjusted HRs 1.00–1.18): the pooled adjusted HR was 1.07 (0.99–1.15), not beyond random gene sets."),
    # B9
    ("The immune half of the dissociation therefore holds within cancer types.", "The immune half of this prognostic–predictive dissociation therefore holds within cancer types."),
    # B21
    ("did not reveal an SPP1⁺ arm opposing a C1Q⁺ arm", "did not reveal an SPP1 arm opposing a C1Q arm"),
    # B10
    ("- SPP1⁺ macrophage states were lower next to TLS than in distant tissue (Ring1−Far lower in 9 of 10 sections; *P* = 0.0098).",
     "- SPP1⁺ macrophage states were lower next to TLS than in distant tissue (Ring1−Far lower in 9 of 10 sections; *P* = 0.0098), whereas an IL-10/STAT3 myeloid state was slightly higher (9 of 10; *P* = 0.0059)."),
    ("- A section-level barrier index varied widely", "- A section-level barrier index (BarrierIndex) varied widely"),
    # B11
    ("(breast −0.23 log₂ units per SD, −0.39 to −0.06).", "(breast −0.23 log₂ units per SD, −0.39 to −0.06; pancreas −0.34, −0.70 to 0.02)."),
    # B12
    ("We asked two things of an interpretable immune–stromal score for gastric cancer: how close it comes to what bulk transcriptomes can offer, and which of its parts carries which information.",
     "Two of our questions concern the score itself: how close it comes to what bulk transcriptomes can offer, and which of its parts carries which information."),
    # B13
    ("Among treated patients in 20 checkpoint cohorts it did not distinguish responders.",
     "Among treated patients in the 20 checkpoint cohorts in which it could be scored, it did not distinguish responders."),
    # B4
    ("and again in 10 new cohorts,", "and again in the 10 replication cohorts,"),
    # B7
    ("In lung and bladder cancer the immune modules were not prognostic in untreated tumours, yet they marked checkpoint responders.",
     "In lung and bladder cancer the Immune axis was not prognostic in untreated tumours, yet it marked checkpoint responders."),
    # B5
    ("with a pathological end point", "with a pathological endpoint"),
    ("one with a different end point,", "one with a different endpoint,"),
    ("comparator and end point.", "comparator and endpoint."),
    # B14
    ("Subtracting Barrier from the immune modules improved prognosis because it added the stromal signal. It did not help response in single-arm cohorts,",
     "Subtracting Barrier from the immune modules improved prognostic performance because it added the stromal signal. It did not strengthen the association with response in single-arm cohorts,"),
    # B15
    ("one small oesophageal cohort treated with chemoradiotherapy.", "one small oesophageal cohort treated with atezolizumab and chemoradiotherapy."),
    # B3
    ("powered for the MSS/EBV-negative majority", "powered for the microsatellite-stable, EBV-negative majority"),
    # B17
    ("External checkpoint cohorts and IMmotion150: as listed in the cited studies and the companion benchmark;",
     "External checkpoint cohorts: as listed in the cited studies and the companion benchmark; IMmotion150 and PRINCE: CRI iAtlas, through cBioPortal (rcc_iatlas_immotion150_2018 and paad_iatlas_prince_2022);"),
    # B18
    ("(d) HR per SD (95% CI), unadjusted (filled circles)", "(d) HR per SD (95% CI) for NetImmune (blue) and Barrier (green), unadjusted (filled circles)"),
    # B16, figure 5 legend
    ("and PRINCE (sensitivity). Every axis", "and PRINCE (progression-free survival; sensitivity). Every axis"),
]

SUPP = [
    # B16
    ("The PRINCE trial (nivolumab-containing versus sotigalimab-containing chemotherapy in pancreatic cancer)[@padron2022]",
     "The PRINCE trial in metastatic pancreatic cancer (nivolumab-containing arms versus sotigalimab without nivolumab, all with gemcitabine and nab-paclitaxel; progression-free survival)[@padron2022]"),
    # B23
    ("- **Programme:** a single-cell gene programme describing a cell state.",
     "- **Programme:** a gene programme describing a cell state; \"immune programmes\" in the title and main text refers to the immune modules collectively."),
    # B22
    ("(the published weights were not available offline)", "(an unweighted version of the published score)"),
    ("shown in place of TIDE, which could not be computed offline.", "shown in place of TIDE, which requires an external web service."),
    # B24
    ("Pancreatic regions were modelled with a patient random intercept and a patient-cluster bootstrap.",
     "Pancreatic regions were modelled with a patient random intercept and a patient-cluster bootstrap; 33 of the 34 regions entered the models, and 7 patients had lymphoid aggregates."),
    # B30
    ("**Supplementary Fig. S2. Programme scores by cell type** (box, interquartile range; 126,559 cells).",
     "**Supplementary Fig. S2. Programme scores by cell type** (box, interquartile range; 126,559 cells). (a–g) TLS, cytotoxic, exhaustion, IFN, CAF/ECM, TGF-β and suppressive-TAM programmes."),
    # A2
    ("(b) Elastic-net module weights.", "(b) Elastic-net coefficients on cell-abundance and cell-state features, by held-out cohort."),
    # B29
    ("(b) Feature comparison with Squidpy, BANKSY, SpatialLDA and CellCharter.",
     "(b) Feature comparison with Squidpy[@palla2022], BANKSY[@singhal2024], Spatial-LDA[@chen2020lda] and CellCharter[@varrone2024]."),
]


def apply(fn, reps):
    p = HERE / fn
    s = p.read_text()
    for a, b in reps:
        n = s.count(a)
        assert n == 1, (fn, n, a[:80])
        s = s.replace(a, b)
    p.write_text(s)
    print(fn, len(reps), "replacements")


if __name__ == "__main__":
    apply("lean_src.md", MAIN)
    apply("lean_supp_src.md", SUPP)
