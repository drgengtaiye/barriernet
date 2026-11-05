# E41 protocol: inflamed expression signatures in published randomized checkpoint trials (aggregate data)

Written 2026-11-01, before any search or extraction is done. The hash of this file is logged in E41_prereg.sha256 first.

## Why

The manuscript's randomized evidence is six comparisons with expression from both arms (744 events in the three renal
trials). The main weakness a reader will raise is power. Many randomized trials of checkpoint inhibitors have published
exploratory biomarker analyses that give the treatment effect separately in patients with high and low values of a
T-cell-inflamed or interferon-γ expression signature. From these, the treatment × signature interaction can be estimated
as a ratio of hazard ratios without patient data. This analysis collects those results, pools them, and asks whether
they agree with the six-comparison estimate for the same class of signatures (Inflamed/IFN family, median +0.19 per SD;
T-cell-inflamed expression profile +0.20 per SD), and how the signatures compare with PD-L1 immunohistochemistry and
tumor mutational burden reported in the same trials. It uses published aggregate numbers only.

## Eligibility (per randomized comparison)

- Randomized (phase 2 or 3) comparison in which the arms differ by a PD-1, PD-L1 or CTLA-4 antibody: checkpoint
  regimen versus a control regimen without one (chemotherapy, targeted therapy, placebo or the same backbone alone).
  Comparisons in which both arms contain a checkpoint inhibitor are excluded.
- A tumor RNA expression signature of T-cell inflammation or interferon-γ signaling (for example T-effector,
  T-cell-inflamed gene expression profile, IFN-γ signature, CD8 T-cell or cytolytic signature, immune-infiltration
  score built from such genes). Single genes, PD-L1 mRNA alone, stromal, angiogenesis, myeloid, B-cell or proliferation
  signatures are not eligible as the primary signature.
- The treatment effect (checkpoint versus control) is reported with a 95% CI separately in signature-high and
  signature-low patients, or a treatment × signature interaction estimate with a 95% CI is reported.
- Peer-reviewed journal article (main text or supplement). Conference abstracts and press releases are excluded.
- The six comparisons analyzed with patient data in the manuscript (IMmotion150, CheckMate 025, JAVELIN Renal 101,
  I-SPY2 pembrolizumab, I-SPY2 durvalumab-olaparib, NeoTRIP) are recorded if found but excluded from the primary pool,
  so that the aggregate evidence is independent of the patient-level evidence.
- When several publications report the same comparison, the one with the most patients with signature data is used.

## Search

1. PubMed (NCBI E-utilities), all dates to 2026-10-31, query:
   (randomized OR randomised OR "phase 3" OR "phase III" OR "phase 2" OR "phase II") AND (atezolizumab OR avelumab
   OR durvalumab OR nivolumab OR pembrolizumab OR ipilimumab OR tremelimumab OR cemiplimab OR tislelizumab OR
   camrelizumab OR sintilimab OR toripalimab OR serplulimab OR sugemalimab OR penpulimab OR dostarlimab OR
   spartalizumab) AND ("gene expression" OR "gene signature" OR "expression signature" OR "RNA sequencing" OR
   "RNA-seq" OR "T-effector" OR "T-cell-inflamed" OR "gene expression profile" OR "interferon gamma" OR
   "biomarker analysis" OR "exploratory biomarker" OR "transcriptomic").
2. Seed list of trials known to have published exploratory biomarker work, checked individually whether or not the
   query finds them (being on the list does not make a trial eligible): POPLAR, OAK, IMpower150, IMpower131,
   IMpower133, IMmotion151, IMvigor211, IMvigor130, IMbrave150, IMpassion130, IMpassion031, IMspire150, JAVELIN
   Bladder 100, JAVELIN Lung 200, JAVELIN Gastric 100, JAVELIN Ovarian 100, KEYNOTE-010, KEYNOTE-045, KEYNOTE-048,
   KEYNOTE-061, KEYNOTE-062, KEYNOTE-181, KEYNOTE-426, CheckMate 026, CheckMate 141, CheckMate 214, CheckMate 227,
   CheckMate 649, CheckMate 9ER, DANUBE, MYSTIC, PACIFIC, CASPIAN, EMPOWER-Lung 1, JUPITER-02, JUPITER-06,
   ESCORT-1st, ORIENT-15, RATIONALE-306, GeparNuevo, ALICE.
3. Reference lists of included papers.

Screening and extraction may read the abstract, the PubMed Central full text and supplements where openly accessible,
and publisher full text where openly accessible; the source read is recorded.

## Items (per eligible comparison)

Trial, arms, tumor, line/setting, phase; endpoint used; signature name and number of genes; cutoff (median, tertile,
other) and the fraction classed high; checkpoint-versus-control HR (or OR) with 95% CI in the high and the low group,
or the interaction estimate with 95% CI; whether a formal interaction test is reported (and its P); the same items for
PD-L1 immunohistochemistry (the trial's own primary PD-L1 cutoff) and for tumor mutational burden when reported in the
same paper for the same endpoint; whether any number was read from a figure rather than a table or text; location
(table, figure, page) and a verbatim quote or the figure description supporting the numbers.

Rules: endpoint hierarchy OS, then PFS, then objective response or pathological complete response. If a paper reports
several eligible signatures, the primary is the one the paper names as primary; otherwise the T-effector or
T-cell-inflamed signature; otherwise the one with the most genes; others are recorded as secondary. For comparisons with
more than one checkpoint arm against one control, each checkpoint arm is a separate comparison only if reported
separately, and the pool then uses the arm with more patients.

## Procedure

Two independent search-and-extraction passes, each by a separate AI-assistant agent given only this protocol and the
instruction not to open the project's results or manuscripts; neither sees the other's output. Disagreements on
eligibility or numbers are adjudicated from the quoted evidence and, where needed, the paper; the adjudicated table is
primary. Agreement on eligibility is reported. Because the extraction is AI-assisted, the authors must check the
adjudicated table against the papers before submission, and the manuscript says so.

## Analyses

Effect measure: benefit-oriented log ratio of HRs, y = −[log HR(high) − log HR(low)], so that y > 0 means more benefit
from the checkpoint inhibitor in signature-high tumors; SE = sqrt(SE_high² + SE_low²), each SE from the reported 95% CI
on the log scale (the two subgroups are disjoint). A reported interaction HR is used directly. Response endpoints use
log odds ratios the same way.

1. Primary: random-effects meta-analysis (REML; Hartung-Knapp-Sidik-Jonkman 95% CI) of y over the eligible comparisons
   with a survival endpoint, excluding the six patient-level comparisons. Reported with τ², I² and the 95% prediction
   interval.
2. Per-SD scale: y divided by the expected difference in a standard-normal score between the high and low groups at
   the reported split (q = fraction low; φ(z_q)/(q(1−q)); 1.60 for a median split). This assumes a roughly normal score
   and a log-linear interaction and is stated as an approximation. Compared with the six-comparison patient-level
   estimates (Inflamed/IFN median +0.19; T-cell-inflamed expression profile +0.20 per SD).
3. Secondary: all endpoints (survival and response) together; PD-L1 immunohistochemistry and tumor mutational burden
   pooled the same way over the comparisons that report them; the paired, within-comparison difference between the
   signature's y and PD-L1's y, described by count (correlation between them is unknown, so no test); the number of
   comparisons reporting a formal interaction test.
4. Sensitivity: excluding numbers read from figures; leave-one-out; phase 3 only; by tumor type where at least three
   comparisons; adding the six patient-level comparisons' own published subgroup results if any were found.
   Small-study asymmetry (Egger) only if at least 10 comparisons.

## Expectations

- M1: at least 10 eligible comparisons with a survival endpoint.
- M2: the pooled y (primary) is positive with a 95% CI excluding 0.
- M3: the pooled per-SD estimate lies between +0.10 and +0.30.
- M4: fewer than half of the eligible comparisons report a formal interaction test.
- M5: over comparisons that report both, the pooled y of the signature and of PD-L1 immunohistochemistry differ by less
  than 0.20.

Reported whatever the result.
