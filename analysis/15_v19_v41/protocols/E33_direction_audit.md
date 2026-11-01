# E33: audit of the declared ICB direction for the 52 signatures

Companion to `E33_direction_audit.csv` (one row per signature, with a verbatim evidence quote and its source).

**Inputs.** The signature list and families come from `icb/v5_signature_clusters.csv`. The declared direction comes from `icb/sig_methods_clean.csv`, which is identical to `signature_information_io.csv` for all 52. All 52 are in the IO set, and the gene lists are in `SignatureSets-main/data-raw/IO/`.

**Method.**
- Each signature was checked against its original paper: first the PubMed abstract (via NCBI E-utilities), then the PMC full text when the abstract gave no direction.
- "audited" records the direction the paper gives for a HIGHER score: sensitivity, resistance, or not stated.
- "agrees" is marked `no` only when the paper clearly supports the opposite direction. It is marked `unclear` when the paper gives no ICB direction, when the direction is non-monotonic, or when the SignatureSets gene file or weights cannot reproduce the paper's score.

**Result: 37 yes, 3 no, 12 unclear.**

## Disagreements (3)

| signature | declared | audited | reasoning |
|---|---|---|---|
| **CD39-CD8Tcell_Chow** (Chow 2023, Immunity, PMID 36574773) | resistance | sensitivity | The abstract says: "Higher baseline frequency of CD39+ CD8+ T cells conferred improved clinical outcomes from ICB therapy" and "a gene signature of CD39+ CD8+ T cells predicted benefit from ICB, but not chemotherapy". The SignatureSets information text also calls it a response predictor, so the label contradicts the curator's own description. |
| **CYT_Davoli** (Davoli 2017, Science, PMID 28104840) | resistance | sensitivity (inferred) | The abstract says: "highly aneuploid tumors show reduced expression of markers of cytotoxic infiltrating immune cells, especially CD8+ T cells". In the title, aneuploidy correlates "with reduced response to immunotherapy". The cytotoxic score itself was not tested against ICB outcome (full text PMC5592794). However, the paper places low cytotoxic infiltrate in the tumours that respond poorly, and its 23-gene aneuploidy-up/immune-down signature was linked to worse survival on immunotherapy. The label appears to carry the aneuploidy direction over to the immune score, and no reading of the paper supports "high cytotoxic score = resistance". Under a strictly literal standard this would be `unclear`, but the inferred direction is unambiguous. |
| **IMS_Cui** (Cui 2021, npj Genom Med, PMID 33542239) | sensitive | resistance | IMS is the *immunosuppression* signature (FAP, PDGFRB, CD163, IL10, CCL2, AXL, VCAN ...). The full text (PMC7862369) says: "elevated expression levels of IMS genes were observed in patients from the nonresponder groups in the IFN-γ+ subgroup". The response predictor in this paper is the IFN-γ/IMS **ratio**, not IMS alone. The SignatureSets information text itself says the genes are "associated with negative response to ICB". This disagreement was not previously known. |

## Unclear (12)

### A. No ICB direction in the source paper (8)

| signature | declared | why unclear |
|---|---|---|
| CIN25_Carter, CIN70_Carter (Carter 2006) | resistance | These are pan-cancer prognostic chromosomal-instability signatures: "Net overexpression of this signature was predictive of poor clinical outcome". The paper has no immunotherapy data, so "resistance" is an extrapolation from later aneuploidy literature. |
| TUMOR-VASCULATURE-UP_Lu (Lu 2007) | resistance | The signature compares ovarian tumour endothelium with normal endothelium. The paper reports no immunotherapy or outcome association. |
| Chemokine_Messina (Messina 2012) | sensitive | A higher score marks lymph-node-like structures and "better overall survival". Immunotherapy is raised only as a question in the title. The label is consistent if anything, but it is not an ICB finding. |
| CTL_Jiang (Jiang 2018, TIDE) | sensitive | CTL level is the mean of CD8A, CD8B, GZMA, GZMB and PRF1, and is used inside TIDE. The paper does not test CTL alone against ICB response. Its survival benefit holds "only when TGFB1 has a low expression level" (PMC6487502). The label is plausible but not reported. |
| ImmunoScore_Hao (Hao 2018) | sensitive | The setting is ovarian cancer, using public non-ICB data. A higher score goes with better OS and chemosensitivity. The abstract says the score predicts "the response to immunotherapy", but no ICB cohort or direction was found. The label is consistent if anything. |
| Macrophage_Danaher (Danaher 2017) | sensitive | This is a descriptive marker paper. On-treatment, "most cell types" were higher in anti-PD1 responders, but macrophages are not singled out. In TCGA, macrophage enrichment "generally indicates poor prognosis". The label is weakly consistent. |
| Myeloid-DC_Helmink (Helmink 2020) | resistance | The gene list is MCP-counter's myeloid DC set. The full text (PMC8762581) names myeloid DCs only as an MCP-counter population and states no responder/non-responder direction. The SignatureSets claim ("more expressed in non-responder") may come from a figure panel, but it could not be verified. |

### B. Non-monotonic (1)

| signature | declared | why unclear |
|---|---|---|
| proliferation_Pabla (Pabla 2019) | resistance | Moderately proliferative tumours did best, and "a highly or poorly proliferative tumor microenvironment is associated with limited sensitivity to ICIs" (PMC6359802). No single direction applies to a higher score. "resistance" fits only the comparison between high and moderate proliferation. |

### C. The paper's direction is clear, but the SignatureSets gene file or weights do not reproduce the paper's score (3)

These are the most consequential for any analysis that scores the gene files directly.

| signature | declared | paper | problem |
|---|---|---|---|
| **PTEN-MITF_Cabrita** (Cabrita 2020, Cancers) | resistance | The MITFlow/PTENneg class has worse survival, including in ICB cohorts. | The paper uses a centroid classifier: "152 were downregulated and 127 were upregulated in MITFlow/PTENneg tumors". The SignatureSets file has 263 genes, all weighted +1, scored as a weighted mean. By gene identity, rows 1-150 (ITGB7..MLANA; e.g. CD8A, CXCL9, CXCL10, GZMB, HLA-B, STAT1, MITF, TYR, MLANA) are the immune/melanocytic block that is **down** in the poor-outcome class. Rows 151-263 (CADM1..SPRY1; adhesion/ECM/Wnt genes) are **up**. An unsigned mean dominated by the immune block may therefore track *benefit*, the opposite of the label. |
| **IRG_Yang** (Yang 2019, cervical cancer) | sensitive | A higher Cox risk score goes with poorer OS and a lower TCGA immunophenoscore, so the paper predicts a lower ICI response (in silico only, with no ICB cohort). | This is the opposite of "sensitive". However, the SignatureSets ±1 weights match the tumour-vs-normal "Regulation" column of Table 2, not the signs of the Cox coefficients: LEPR, PRL, NRP1 and ANGPTL5 have positive coefficients but weight −1, and IFI30 has a negative coefficient but weight +1. The implemented score is therefore not the paper's risk score. |
| **ImmuneCells_Xiong** (Xiong 2020) | sensitive | In the paper, the classifier z-score >0.19 predicts response. | The paper uses nearest-centroid classification, not a gene-set score. The 109 genes mix markers of TREM2hi macrophages and γδ T cells, which are enriched in non-responders (TREM2, SPP1, MMP9, CCL18), with markers of a B-cell subset enriched in responders. The paper does not determine the direction of an unsigned GSVA score. |

## Caveats on rows marked "yes"

- **C-ECM_Chakravarthy (scoring caveat).** The paper scores C-ECM-up (30 genes) and C-ECM-down (28 genes) separately, and the resistance direction belongs to C-ECM-up ("the C-ECM-up score was significantly higher in progressors"). The SignatureSets file merges both sets: rows 1-28 (TNXB..MFAP5) are the down set and rows 29-58 (MMP11..ADAM12) are the up set. All 58 genes are scored as one ssGSEA set. Both halves are stromal/ECM genes, so the label is kept as "yes", but this is not the published score.
- **IMPRES_Auslander.** IMPRES is implemented as a weighted mean of ±1 genes instead of the 15 pairwise comparisons. The weight directions are consistent with the paper. **Blood_Friedlander** is a whole-blood classifier ("9 predictors and six non-predictive enhancer variables"), so an unweighted score is an approximation.
- **Yes by indirect or weaker evidence.**
  - ImmuneScore_Roh: the score correlates with TCR clonality, a response marker; no responder vs non-responder test of the score itself.
  - TLS_Cabrita: direction comes from the title plus the abstract; the full text is paywalled.
  - VIGex_Hernando-Calvo: "Hot" is the highest-score group; the direction is implied.
  - GeparSixto_Karn: the endpoint is pCR on chemotherapy ± durvalumab, so the result is not ICB-specific.
  - MHC-I_Liu: non-significant trend only.
  - PDL1_Nishino: a review of PD-L1 IHC with "contradictory results"; no gene signature is defined.
  - STAT1_Care and KDM5A_Wang: the ICB evidence comes from mouse models.
  - Response_Chen and NonResponse_Chen: the direction holds by construction. However, "we did not observe significant associations of our derived gene signatures with response" in human ICB bulk cohorts, and the gene counts (443/258) differ from the 85/69 genes reported in the text.
- **Metadata errors that do not affect direction.**
  - Teff_McDermott: the regimen is atezolizumab + bevacizumab, not "+ sunitinib".
  - M1_Hwang: the cancer type is NSCLC, not melanoma.
  - IRG_Ayers: LAG3 is missing from the 18-gene expanded immune set.

## Limitations

- Quotes were extracted with a web-fetch tool, from PubMed abstracts (E-utilities) and PMC full text, in fragments of 125 characters or fewer. They were not checked by hand against the PDFs. Spot-check any quote before using it in the manuscript.
- No full text was accessible for Cabrita 2020 (Nature), Karn 2020 (Ann Oncol), Hao 2018 (Clin Cancer Res), Hernando-Calvo 2023 (Med) or Zemek 2019 (Sci Transl Med). These rows rest on the abstracts only.
- The gene-block assignments for PTEN-MITF and C-ECM are inferred from gene identity and block sizes. The supplementary tables were not accessed.
