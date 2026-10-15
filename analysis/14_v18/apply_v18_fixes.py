"""Fixes and cuts from analysis/results/VERIFY_v18.md, applied after apply_v18_edits.py
(order: apply_v18_edits.py -> apply_v18_fixes.py -> build_v18.py). Exact replacements and marker-delimited
block rewrites; every target must match once."""
import pathlib
HERE = pathlib.Path(__file__).resolve().parent


def rep(s, a, b):
    assert s.count(a) == 1, (s.count(a), a[:100])
    return s.replace(a, b)


def block(s, start, end, new):
    """Replace from the start marker to the end of the first end marker after it."""
    assert s.count(start) == 1, ("start", start[:80])
    i = s.index(start); j = s.index(end, i) + len(end)
    return s[:i] + new + s[j:]


def main(s):
    # A2 abstract
    s = rep(s, "then tested them, locked, in seven more (1,055 patients). Checkpoint",
            "then tested them, locked, in seven more (1,055 patients), four added by pre-registered extension after the first three gave *P* = 0.057. Checkpoint")
    s = rep(s, " Later validations were pre-registered.", "")
    s = rep(s, "NetImmune was prognostic in the benchmark cohorts (pooled HR per SD 0.74, 95% CI 0.58–0.95) and the locked cohorts (0.70, 0.62–0.80), beyond stage and random axes; it matched ESTIMATE and TMEscore.",
            "NetImmune was prognostic beyond stage and random axes in the benchmark (pooled HR per SD 0.74, 95% CI 0.58–0.95) and locked cohorts (0.70, 0.62–0.80), matching ESTIMATE and TMEscore.")
    s = rep(s, "CD8 exclusion tracked fibroblast number; arrangement had no detectable effect (underpowered test).",
            "CD8 exclusion tracked fibroblast number, not arrangement (underpowered test).")
    s = rep(s, "(AUC 0.80), but not beyond", "(AUC 0.80), not beyond")
    # O12 introduction
    s = rep(s, " Beyond these, the tumour microenvironment is the most studied source of predictive information.", "")
    # O14, C2, C1 methods data list
    s = rep(s, "[@tcga2014] (n = 383; 156 deaths), GSE84437[@yoon2020] (431; 207) and GSE62254 (ACRG)[@cristescu2015] (283; 135).",
            "[@tcga2014], GSE84437[@yoon2020] and GSE62254 (ACRG)[@cristescu2015] (1,097 patients).")
    s = block(s, "- **Locked gastric cohorts,**", "(376 patients after duplicate screening).",
              "- **Locked gastric cohorts,** analysed afterwards with every parameter fixed: GSE15459 and GSE34942 with overall survival[@ooi2009; @lei2013] and GSE26253 with recurrence-free survival[@lee2014] (679 patients), and, in a pre-registered extension, GSE13861, GSE26899, GSE26901 and GSE29272 with overall survival[@cho2011; @oh2018; @wang2013] (376 patients after duplicate screening).")
    s = block(s, "- **Randomised comparisons:**", "[@grande2026].",
              "- **Randomised comparisons,** with progression-free survival or pathological complete response (pCR) as the endpoint:\n"
              "  - IMmotion150 in renal cell carcinoma[@mcdermott2018], and two I-SPY2 comparisons in HER2-negative breast cancer[@wolf2022; @nanda2020; @pusztai2021];\n"
              "  - in a pre-registered extension, every further randomised trial with open pre-treatment expression in both arms: CheckMate 025[@motzer2015; @braun2020] and JAVELIN Renal 101[@motzer2019; @motzer2020] in renal cell carcinoma, and NeoTRIP in triple-negative breast cancer[@gianni2022; @dugo2026];\n"
              "  - as sensitivity analyses: PRINCE in pancreatic cancer, whose control arm also received an immunotherapy[@padron2022], and DUTRENEO in bladder cancer, randomised only among immune-high tumours[@grande2026].")
    # C10 / B5
    s = rep(s, "The final validations (the locked-cohort extension, the further randomised trials, the independent single-cell atlases, PD-L1 CPS and robustness to gene loss) were specified in a dated, hash-logged protocol before any of their outcomes were examined; deviations are listed in the Supplementary Methods.",
            "The locked-cohort extension, further randomised trials, independent single-cell atlases, PD-L1 CPS and robustness analyses were specified in a dated, hash-logged protocol before they were run; deviations are listed in the Supplementary Methods.")
    # C9, O9 statistics
    s = rep(s, "  - Score × arm interaction per SD: Cox models for IMmotion150, CheckMate 025 and JAVELIN Renal 101; logistic models (odds ratios, ORs) for I-SPY2, adjusted for hormone-receptor status, which adaptive randomisation left unbalanced between arms, and for NeoTRIP, adjusted for sequencing batch.",
            "  - Score × arm interaction per SD: Cox models for progression-free survival and logistic models (odds ratios, ORs) for pCR, adjusted for hormone-receptor status in I-SPY2 (left unbalanced by adaptive randomisation) and for sequencing batch in NeoTRIP.")
    s = block(s, "- **Spatial analysis.**", "with a simulated barrier as positive control.",
              "- **Spatial analysis.** On Visium, zones were defined around TLS peaks, with a least-cost transport analysis of stromal obstruction. In IMC, fibroblast-layer continuity at the tumour boundary, split into the part expected from fibroblast number and the excess over random placement, was related to CD8 infiltration, adjusting for fibroblast and CD8 density, with a simulated barrier as positive control.")
    # P7 software
    s = rep(s, " A code map links every table, figure and number to the script that produces it.", "")
    # B11 results, independent atlases
    s = rep(s, "T and NK cells carried 97% and 93% of T_CELL signal", "T and NK cells carried 97% and 93% of T_CELL signal (84% in the discovery atlas)")
    # O2, O3 prognosis
    s = rep(s, "(Fig. 2a–c; Table 1). The pooled HR per SD was 0.74 (95% CI 0.58–0.95; *P* = 0.036; Supplementary Fig. {SF:meta_forest}).",
            "(Fig. 2a–c; Table 1; pooled HR per SD 0.74, 95% CI 0.58–0.95; *P* = 0.036; Supplementary Fig. {SF:meta_forest}).")
    s = rep(s, "in each cohort (0.77, 0.80 and 0.69; Fig. 2d). It also remained prognostic after further adjustment",
            "in each cohort (Table 1; Fig. 2d) and after further adjustment")
    # O1 fitted models
    s = block(s, "Fitted models did little better:", "(Supplementary Fig. {SF:upgrades}).",
              "Fitted models did little better. The best, Cox regression on 20 principal components chosen after the fact from eight model classes, gained 0.024 in held-out C (0.000–0.048) and was still improving at the full training size (Supplementary Fig. {SF:fitted_models}; Supplementary Table {ST:fitted}). Separating cell abundance from cell state, learning the module weights or adding a spatial arrangement score each changed C by less than 0.01 (Supplementary Fig. {SF:upgrades}).")
    # C6 locked replication (with B1)
    s = block(s, "In three further cohorts analysed with every parameter locked", "third-party curation[@li2023].",
              "In three further cohorts analysed with every parameter locked, NetImmune had an HR below 1 in each (0.65, 0.92 and 0.75), but the pooled estimate, 0.73 (0.53–1.02; *P* = 0.057), was not significant under the conservative rule. After this result we pre-registered an extension to four further cohorts with recoverable survival; before any outcome model, an expression-correlation duplicate screen excluded 16 of their arrays, none of which matched a patient of another cohort on clinical data (Table 1; Supplementary Fig. {SF:locked}; Supplementary Table {ST:locked}). "
              "NetImmune's HR was below 1 in each added cohort (0.50–0.72; all *P* < 0.05) and was 0.70 across all seven locked cohorts (0.62–0.80; *P* = 4.3 × 10⁻⁴; I² = 7%), beyond random axes (*P* ≤ 0.001). CAF_ECM was adverse in six of the seven (pooled 1.34, 1.08–1.65). NetImmune exceeded CD8A in C (+0.062, 0.032–0.091) but not ESTIMATE (+0.015, −0.002 to 0.031). "
              "Across all ten gastric cohorts, NetImmune gave 0.72 (0.67–0.79) and CAF_ECM 1.30 (1.17–1.45); estimates were similar (0.68–0.72) without the duplicate exclusions, with overall survival only, or with a cohort whose survival data come from a third-party curation[@li2023].")
    # B12, P3 CPS
    s = rep(s, "discriminated responders better than NetImmune (AUC 0.97 against 0.82 in the same patients); every responder with a CPS had CPS ≥ 1.",
            "had a higher AUC than NetImmune (0.97 against 0.82 in the same patients).")
    # C15
    s = rep(s, "with chemotherapy–checkpoint cohorts added (24 cohorts: Immune +0.42, 0.25–0.59; CAF_ECM −0.07), and pointed",
            "with chemotherapy–checkpoint cohorts added, and pointed")
    # C3, C4, C5 randomised comparisons
    s = block(s, "We first analysed three randomised comparisons (Supplementary Fig. {SF:trials_detail}):", "benefit-oriented log ratio +0.47, −0.04 to 0.97; *P* = 0.057).",
              "We first analysed three randomised comparisons (Table 3; Supplementary Fig. {SF:trials_detail}). In IMmotion150, the Immune axis and NetImmune modified the benefit of atezolizumab with or without bevacizumab over sunitinib (interaction HR per SD 0.68 and 0.65; *P* = 0.035 and 0.015), but CAF_ECM did not (1.21; *P* = 0.25; Supplementary Fig. {SF:icb11}). "
              "In I-SPY2, NetImmune and CAF_ECM modified the effect of adding pembrolizumab to paclitaxel in opposite directions (*P* = 0.031 and 0.019; random-gene null *P* = 0.020 and 0.010), less clearly within hormone-receptor strata (Supplementary Table {ST:trials}); neither modified the effect of durvalumab–olaparib. "
              "Pooled over the three comparisons, NetImmune favoured the checkpoint arm in each, but not significantly (benefit-oriented log ratio +0.47, −0.04 to 0.97; *P* = 0.057).")
    s = block(s, "In a pre-registered extension we added every further randomised trial", "null *P* = 0.001).",
              "In a pre-registered extension, decided after this result, we added every further randomised trial with open expression data from both arms (Fig. 5; Tables 3 and 4). NetImmune leaned towards the checkpoint arm in CheckMate 025 and JAVELIN Renal 101 (*P* = 0.62 and 0.090). In NeoTRIP, NetImmune and CAF_ECM modified the benefit of atezolizumab in opposite directions (*P* = 0.027 and 0.048; random-gene null *P* = 0.019 and 0.001).")
    s = block(s, "Across the six primary comparisons (Table 4):", "significantly in the two breast cancer trials.",
              "Across the six primary comparisons (Table 4), NetImmune favoured the checkpoint arm in all six (pooled +0.29, 0.03–0.55; *P* = 0.034), and the pooled Immune and IFN interactions also favoured it (*P* = 0.021 and 0.027); CAF_ECM favoured it in none (*P* = 0.14). "
              "The NetImmune result was fragile: leaving out any single trial, it stayed positive (+0.23 to +0.37) but *P* exceeded 0.05 three times, and *P* was 0.059 with overall survival for CheckMate 025. "
              "In a subset analysis that was not pre-specified, effects were clearer in the three breast cancer comparisons (NetImmune +0.63, *P* = 0.048; CAF_ECM −0.60, *P* = 0.034) than in the three renal ones (NetImmune +0.20, *P* = 0.15).")
    # O5, O6 spatial results
    s = block(s, "On Visium, no stromal or myeloid barrier surrounded lymphoid aggregates", "(Supplementary Fig. {SF:transport}).",
              "On Visium, no stromal or myeloid barrier surrounded lymphoid aggregates (Supplementary Figs {SF:visium}, {SF:visium_sens}). SPP1⁺ macrophage states were lower next to TLS than in distant tissue (Ring1−Far lower in 9 of 10 sections; *P* = 0.0098), an IL-10/STAT3 myeloid state was slightly higher (9 of 10; *P* = 0.0059), and fibroblast states showed no consistent gradient. "
              "A section-level barrier index (BarrierIndex) varied widely and was not associated with age, pT or pN; stromal and myeloid states were spatially organised (Supplementary Figs {SF:barrierindex}, {SF:moran}); and a least-cost transport analysis gave no evidence in either direction once a compositional artefact of deconvolution had been removed (Supplementary Fig. {SF:transport}).")
    s = rep(s, "were denser, not sparser, than further away. The lower SPP1⁺ state next to TLS on Visium is therefore unlikely to reflect fewer macrophages; the imaging panels lack SPP1, so state itself was not measured.",
            "were denser, not sparser, so the lower SPP1⁺ state next to TLS on Visium is unlikely to reflect fewer macrophages (the imaging panels lack SPP1).")
    # C13, D2, C8, D1, C14, C7, O10, O7 discussion
    s = rep(s, "Two of our questions concern the score itself: how close it comes to what bulk transcriptomes can offer, and which of its parts carries which information.\n\nOn the first, NetImmune was prognostic",
            "NetImmune was prognostic")
    s = rep(s, "The second answer is the more useful one, and it concerns the components rather than the composite:",
            "The more useful findings concern the components rather than the composite:")
    s = rep(s, ", and they favoured the checkpoint arm across six randomised comparisons.", ".")
    s = rep(s, "Its interaction with treatment was consistent in direction across six randomised comparisons and nominally significant when pooled, but it was small, clearest in breast cancer and sensitive to the removal of single trials. This reading does not rest only on comparisons between cancer types. In lung and bladder cancer the Immune axis",
            "Its interaction with treatment was consistent in direction across six randomised comparisons and nominally significant when pooled, but small and clearest in breast cancer. Nor does this reading rest only on comparisons between cancer types: in lung and bladder cancer the Immune axis")
    s = rep(s, "In none of six randomised comparisons did fibroblast content favour the checkpoint arm. It pointed against benefit significantly in the two breast cancer trials with pathological endpoints, I-SPY2 and NeoTRIP, and not in the renal trials. We therefore read fibroblast content as prognostic, with a possible effect against checkpoint benefit that single-arm cohorts are not designed to detect and that may depend on the setting. Fibroblast-rich stroma is linked to checkpoint resistance in specific settings[@mariathasan2018; @kieffer2020].",
            "Fibroblast content never favoured the checkpoint arm, and it pointed significantly against benefit in the I-SPY2 pembrolizumab and NeoTRIP comparisons, both in breast cancer with pathological endpoints, but not in the renal trials. We therefore read it as prognostic, with a possible, setting-dependent effect against checkpoint benefit that single-arm cohorts are not designed to detect; fibroblast-rich stroma is linked to checkpoint resistance in specific settings[@mariathasan2018; @kieffer2020].")
    s = rep(s, " Across 20 cohorts here, however, an ECM-structural module did not identify non-responders.", "")
    s = rep(s, "In gastric cancer the score did not add to the clinical biomarkers. PD-L1 CPS and MSI/EBV status each discriminated pembrolizumab responders better than NetImmune. Its value lies less in selecting patients than in showing which component of the microenvironment carries which information.",
            "In gastric cancer, PD-L1 CPS and MSI/EBV status had higher AUCs than NetImmune, which was not shown to add to them; the score's value lies less in selecting patients than in showing which component of the microenvironment carries which information.")
    s = rep(s, "Two further points concern the macrophage module and the spatial data:\n\n", "")
    s = rep(s, "Neither the spatial nor the single-cell data supported it. CD8 exclusion followed how many fibroblasts lined the tumour, not whether they formed a continuous layer, although a modest physical effect cannot be excluded.",
            "Neither the spatial nor the single-cell data supported it: CD8 exclusion followed fibroblast number, not continuity, although a modest physical effect cannot be excluded.")
    # A1 / C11, B6, O11 limitations
    s = rep(s, "Exploratory analyses of the three benchmark cohorts preceded the fixing of the gene lists, The locked replication was extended from three to seven cohorts after the first three gave *P* = 0.057; the extension was pre-registered and every parameter stayed locked, but the decision to extend was taken in the light of that result.",
            "Exploratory analyses of the three benchmark cohorts preceded the fixing of the gene lists. The locked replication was extended from three to seven cohorts after the first three gave *P* = 0.057; the extension was pre-registered and fully locked, but the decision to extend followed that result.")
    s = rep(s, "it rests mainly on the breast cancer trials and loses significance when any of three trials is removed.",
            "it rests mainly on the breast cancer trials, loses significance when any of three trials is removed, and, like the locked replication, was extended after the first three comparisons gave *P* = 0.057.")
    s = rep(s, "- **Panel coverage.** Some replication cohorts used targeted panels on which not every module could be scored.\n", "")
    # B23
    s = rep(s, "such as CheckMate 649, would test", "such as CheckMate 649[@janjigian2021], would test")
    # B20 data availability
    s = rep(s, "with survival for the last four from the original publications.",
            "with survival for the last four from the original publications. Sensitivity cohort: GEO GSE57303, with survival from a third-party curation.")
    # B8 figure 1 legend, B9 table 3 note and figure 5 legend
    s = rep(s, "defined and characterised in a single-cell atlas and 10 Visium sections",
            "defined and characterised in a single-cell atlas (with two independent atlases for cell-type attribution) and 10 Visium sections")
    s = rep(s, "- **Further results.** All scores and within-arm effects are in Supplementary Table {ST:trials}.",
            "- **Further results.** All scores and within-arm effects are in Supplementary Table {ST:trials}.\n- **Abbreviations.** n, checkpoint / control arm; PFS, progression-free survival; TNBC, triple-negative breast cancer.")
    s = rep(s, "(−log HR for progression-free survival; log OR for pCR)", "(−log HR for progression-free survival, PFS; log OR for pCR)")
    s = rep(s, "in the three added trials, with the observed value and empirical *P*", "in each added trial and combined by inverse variance, with the observed value and empirical *P*")
    return s


def supp(s):
    s = rep(s, "21 checkpoint cohorts and three randomised comparisons**", "21 checkpoint cohorts and six randomised comparisons**")          # A4
    s = rep(s, "Supplementary Figures S1–S20;", "Supplementary Figures {SF:sc_atlas}–{SF:drugs};")                                          # A5
    s = rep(s, "Supplementary Tables {ST:locked}, {ST:modules}, {ST:prjeb_nulls} are sheets",
            "Supplementary Tables {ST:robust}, {ST:locked}, {ST:modules} and {ST:prjeb_nulls} are sheets")                                 # A6
    s = rep(s, "NetImmune 0.825 (0.661–0.953).", "NetImmune 0.824 (0.661–0.953).\n- **Correlation:** NetImmune with CPS, Spearman ρ = 0.42 (*P* = 0.006).\n- **Responders:** every responder with a CPS had CPS ≥ 1.")  # A8, A11, P3
    s = rep(s, "Eight other gastric series downloaded for this purpose had no deposited outcome and were not analysed.",
            "Eight other gastric series downloaded for this purpose had no deposited outcome and were not analysed at that stage; for four of them, outcomes were later recovered from publications (below).")  # A9
    s = rep(s, "The other four were added in a pre-registered extension (see Pre-registered validation analyses),",
            "The other four were added in a pre-registered extension, decided after the first three gave *P* = 0.057 (see Pre-registered validation analyses),")  # B7
    s = rep(s, "specified in a protocol written before any of their data were related to outcome.",
            "specified in a protocol written before any of these analyses was run and before any new data were related to outcome.")  # B5
    s = block(s, "Amendments were made before the relevant outcomes were examined, as data sources were identified:",
              "Rescoring that cohort from the authors' matrix changed no pooled response estimate by more than 0.015.\n",
              "Amendments were made before the relevant outcomes were examined, as data sources were identified:\n\n"
              "- the choice of atlases and cohorts, and the use of the authors' published cell-type annotations in place of marker-based re-annotation;\n"
              "- the CPS transformation;\n"
              "- the trial endpoints and the batch covariate for NeoTRIP.\n\n"
              "The locked-cohort extension and the trial extension were decided after the first three locked cohorts and the first three randomised comparisons each gave *P* = 0.057.\n\n"
              "Deviations recorded after analysis:\n\n"
              "- **Robustness, panel subset.** On the 22 genes shared by all four panel cohorts, CAF_ECM keeps only one gene, so each panel platform's own gene content was scored instead.\n"
              "- **Robustness, gene retention.** A draw was redrawn if any module scorable in a cohort would fall below three genes there, so six genes (THY1, FN1, COL3A1, CD3D, CD8A and GZMB) were never dropped. The literal rule (at least three listed genes per module; cohorts that lose a module leave the estimate) was run as a check and gave the same directions in every draw (Supplementary Fig. {SF:robust}, open symbols).\n"
              "- **CheckMate 025 expression.** The compendium copy of the CheckMate expression data used for the first 11 cohorts agrees with the authors' matrix for NetImmune (ρ = 0.99) but less for single modules (ρ 0.85–0.95). Rescoring that cohort from the authors' matrix changed no pooled estimate across the first 11 cohorts by more than 0.015; the 21-cohort estimates were not rescored.\n"
              "- **Analyses not pre-specified.** Pooled interactions by cancer type, and a DUTRENEO analysis restricted to libraries of at least one million counts.\n\n"
              "The analyses were run as follows:\n")  # A10
    s = rep(s, "per added cohort and pooled.", "per added cohort and pooled. The extension was decided after the first three locked cohorts gave *P* = 0.057 (row \"Pooled, locked (first 3)\").")  # A3
    s = rep(s, "boxes mark the pre-specified intended class.",
            "boxes mark the pre-specified intended class. Without LUM, DCN and POSTN (absent from the processed discovery atlas), CAF_ECM peaks in pericytes and smooth-muscle cells in both independent atlases.")  # B21
    s = rep(s, "  - random-dropout summaries;\n  - leave-one-gene-out estimates;\n  - panel-platform scores.",
            "  - random-dropout summaries;\n  - correlation with the full score;\n  - leave-one-gene-out estimates;\n  - panel-platform scores;\n  - pre-specified expectations.")  # B16
    s = rep(s, "  - random-gene null;\n  - ΔC.\n- **Supplementary Table {ST:sc_indep}.**", "  - random-gene null;\n  - ΔC;\n  - pre-specified expectations for the extension.\n- **Supplementary Table {ST:sc_indep}.**")
    s = rep(s, "  - pooled estimates for every set, including leave-one-trial-out.",
            "  - pooled estimates for every set, including leave-one-trial-out;\n  - pre-specified expectations;\n  - a check of the CheckMate expression source.")
    s = rep(s, "pooled by random effects (REML, Hartung–Knapp", "pooled by random effects (restricted maximum likelihood, REML, with Hartung–Knapp")  # B23
    return s


if __name__ == "__main__":
    m = main((HERE / "v18_src.md").read_text())
    u = supp((HERE / "v18_supp_src.md").read_text())
    (HERE / "v18_src.md").write_text(m)
    (HERE / "v18_supp_src.md").write_text(u)
    print("fixes applied")
