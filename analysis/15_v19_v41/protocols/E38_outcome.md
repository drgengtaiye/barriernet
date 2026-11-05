# E38 outcome (protocol E38_prereg.md, sha256 f6e50c4b…, logged 2026-11-01T01:41:29Z)

## Runs
- **Extraction passes.** Two passes by separate AI-assistant agents (`E38_pass1.csv`, `E38_pass2.csv`), run after the hash. Each pass received only the signature names, the PMIDs (`E38_sources.csv`) and the rubric, and saw neither the other pass nor the results. Each checked by script that its quotes appear verbatim in the text it fetched.
- **Reconciliation.** `E38_adjudicate.py` resolves the disagreements in RESOLVE, with a reason for each, and writes:
  - `results/E38_adjudicated.csv`
  - `E38_agreement.csv`
  - `E38_crosstab.csv`
  - `E38_expectations.csv`

## Agreement

| item | agree | κ |
|---|---|---|
| D1 (design) | 50/52 | 0.92 |
| D2 (single gene) | 44/52 | 0.76 |
| D3 (random genes) | 50/52 | 0.90 |
| D4 (PD-L1 IHC or TMB) | 49/52 | 0.91 |

There were 15 disagreements across 8 signatures:
- **Pass 1 read the full text from the journal site, pass 2 only the abstract** (2): ADO_Sidders D2 and ImmunoScore_Hao D2 → yes.
- **Component-gene decompositions do not count as a single-gene comparison** (pass 1's rule; pass 2 had counted them as borderline): CRMA_Shukla, MHC-I_Liu and MHC-II_Liu D2 → no.
- **CRMA_Shukla D4 → no:** neoantigen load is neither TMB nor PD-L1 IHC.
- **APM_Thompson D2 → no:** the single-gene result is cited from an earlier study.
- **Response_Chen and NonResponse_Chen D1 → D:** a null test in ICB-treated bulk cohorts is reported, so D2–D4 → no.

## Expectations

| | expectation | result | held |
|---|---|---|---|
| L1 | at most 3 of 52 with D1 = A | 0 | yes |
| L2 | at least 38 of 52 with D1 = D or E | 45 | yes |
| L3 | at most 12 of 52 compared with a single gene | 18 | **no** |
| L4 | at most 3 of 52 compared with random gene sets | 2 | yes |

## Findings
- **Design of the evidence:** A 0, B 4, C 3, D 36, E 9.
  - The B papers are CD39-CD8 T cell (OAK), Teff McDermott (IMmotion150), Teff-IFNG Fehrenbacher (POPLAR) and GeparSixto Karn (GeparNuevo).
  - The last two were read as abstracts only, so an interaction test in them cannot be excluded. If both were A, L1 would still hold.
- **Comparators:**
  - Single gene: 18. More than expected; the count includes CD274 mRNA shown beside the signature.
  - Random gene sets: 2 (APM_Wang, IMPRES_Auslander).
  - PD-L1 IHC or TMB: 17.
- **Abstract-only papers:** six papers were read as abstracts only in both passes (GeparSixto_Karn, KDM5A_Wang, TLS_Cabrita, Teff-IFNG_Fehrenbacher, VIGex_Hernando-Calvo, WNT-beta-catenin_Takeuchi). Their D2–D4 = "no" means "not in the abstract".
- **Cross-tabulation (descriptive):**
  - All 4 signatures with randomized evidence favored the checkpoint arm at P < 0.05 in the six comparisons, against 13 of the other 48.
  - 3 of the 4 exceeded their random reference.
  - No signature exceeded CD8A, whatever its design.
- **Caveat:** the extraction was AI-assisted, so the authors must check the reconciled table against the papers before submission.

## Tally
Later-protocol tally: 13 protocols, 39 of 52 expectations held.

Correction and disclosure (2026-11-01, after the independent fact-check): the disagreements fall across 8 signatures, not 9. Pass 1 also read publisher or author-manuscript full texts where freely accessible (ADO_Sidders, ImmunoScore_Hao, STAT1_Care), which pass 2 had as abstracts only; the protocol named PubMed Central, so this is a deviation, and two reconciliations (ADO_Sidders D2, ImmunoScore_Hao D2) rest on it. The two abstract-only B papers (POPLAR, GeparNuevo) could be A; the text says 'none of the papers read in full'.
