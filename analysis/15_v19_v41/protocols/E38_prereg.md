# E38 protocol: what evidence the source papers give for the predictive claims of the 52 signatures

Written 2026-10-31, before any extraction is done. The hash of this file is logged in E38_prereg.sha256 first.

## Why

The manuscript argues that predictive claims for expression signatures need a treatment interaction in randomized
data, a size-matched random reference and a single-gene comparator. A reader will ask whether the published claims
were in fact made without these. This analysis records, for each of the 52 library signatures, the design of the
checkpoint-outcome evidence in its source paper and whether the paper compared the signature with a single gene or with
random gene sets. It is descriptive meta-research on the library already analyzed; it uses no patient data.

## Sources

One source paper per signature: the PMID in E33_direction_audit.csv (STAT1_Care: Zemek 2019, PMID 31316010, the paper
that tested it in checkpoint-treated patients; PDL1_Nishino: both PMIDs listed there are read). PubMed abstract via NCBI
E-utilities, and the PMC full text where the paper is in PMC; the source used is recorded. Supplementary files are not
read.

## Items (per signature; about the signature itself, not other scores in the paper)

- **D1, highest design of checkpoint-outcome evidence:**
  - A, randomized comparison of a checkpoint inhibitor with a non-checkpoint control, with a formal test of the
    signature × treatment interaction (or of a difference in the signature's effect between the arms);
  - B, randomized comparison with both arms analyzed, without a formal interaction test (within-arm associations or
    subgroup hazard ratios only);
  - C, checkpoint-treated patients, plus separately sourced, non-randomized patients without checkpoint treatment in
    whom the association is shown to be absent or weaker;
  - D, checkpoint-treated patients only (single-arm cohorts or the checkpoint arm of a trial);
  - E, no checkpoint-outcome data for the signature.
- **D2:** compared with a single-gene expression marker for checkpoint outcome (any gene, including CD274/PD-L1 mRNA,
  CD8A, CXCL9): yes / no / not applicable (D1 = E).
- **D3:** compared with random gene sets: yes / no / not applicable.
- **D4:** compared with PD-L1 immunohistochemistry or tumor mutational burden: yes / no / not applicable.
- A verbatim quote supporting D1, and the source read (abstract, full text, abstract only because no full text).

## Procedure

Two independent extraction passes, each by a separate AI-assistant agent given only the signature names, PMIDs, this
rubric and the instruction not to open the project's results or manuscripts; neither sees the other's output.
Disagreements are adjudicated from the quoted evidence and, where needed, the paper; the adjudicated table is primary.
Agreement is reported (percentage and Cohen's κ for D1). Because the extraction is AI-assisted, the authors must check
the adjudicated table before submission, and the manuscript says so.

## Analyses

1. Counts of D1 classes, D2, D3 and D4 (adjudicated).
2. Cross-tabulation of D1 with the framework's verdicts already reported (E32 primary: favored the checkpoint arm at
   P < 0.05; exceeded the size-matched random reference; E33: none differed significantly from CD8A in the checkpoint
   direction). Descriptive; no test is primary.

## Expectations

- L1: at most 3 of the 52 have D1 = A.
- L2: at least 38 of the 52 have D1 = D or E.
- L3: at most 12 of the 52 have D2 = yes.
- L4: at most 3 of the 52 have D3 = yes.

Reported whatever the result.
