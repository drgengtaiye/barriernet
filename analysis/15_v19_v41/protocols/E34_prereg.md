# E34 protocol: a seventh randomized comparison from public data (NeoPembrOV)

Written 2026-10-25, before any analysis involving outcomes in this dataset. Seventh protocol written after all other
results were known; post hoc. The hash of this file is logged in E34_prereg.sha256 before E34_trials.py is run.

## Why

The author asked to add randomized data that can be downloaded directly. A search of GEO, EGA, iAtlas, GSA and sponsor
portals (October 2026) found one new randomized checkpoint comparison with public pre-treatment expression from both
arms: NeoPembrOV (NCT03275506; GINECO randomized phase II; high-grade serous ovarian carcinoma; neoadjuvant
carboplatin–paclitaxel with or without pembrolizumab, randomized 2:1). It shares a chemotherapy backbone, like the
breast comparisons, and has no anti-angiogenic drug in either arm.

## Data (all public)

- Expression: GEO GSE227666, RNA-seq counts (GSE227666_neopembrov_counts.txt.gz), 109 libraries; pre-treatment 32
  NACT+P and 21 NACT (series-matrix titles).
- Outcome: Derbel et al., Nat Commun 2024 (PMC11252308), Source Data file (41467_2024_47000_MOESM4_ESM.xlsx), sheet
  "Figures 4 a,b,c,d": patient ID, arm, progression-free survival (del_pfs2, months; pfs2, event) and overall survival
  for 53 patients; "Table 1": progression status at 24 months (Pr/NPr) for 52.
- Linking: the Source Data use patient IDs, GEO uses library IDs. They were linked before this protocol using
  expression only (no outcome column read): NACT+P by 40 genes given per patient in the Source Data (all 32 assigned;
  correlation 0.79–0.99, second-best ≤ 0.59); NACT by 21 per-patient expression features (10 genes and 11 CD8B
  ratios; all 21 assigned; correlation 0.89–1.00), a method that reproduced all 32 NACT+P links. Links are in
  data_download/validation/NeoPembrOV_GSE227666/id_match_final.csv. Links with a margin over the second-best candidate
  below 0.10 are flagged and dropped in a sensitivity analysis.

## Already known (declared)

All E32/E33 results. From the source paper: pembrolizumab did not improve survival in this trial; within NACT+P, a
higher CD8B/FOXP3 ratio and lower endothelial (KDR/VEGFR2) signal went with response. No signature, single gene or
interaction in these data has been computed by us.

## Model

As E32/E33 (E32_trials.py setup; E33 engine): log2 CPM + 0.5 from the GEO counts (genes summed by symbol), genes
z-scored across the 53 pre-treatment samples of both arms, library scored by the common engine (≥ 3 genes and ≥ 50% of
the list), single genes as z-scores, barriernet scores by its published function. Primary endpoint: PFS, Cox model with
score, arm and score × arm (no covariates; Firth not needed for Cox), benefit-oriented interaction = −log HR per SD.
Sensitivity: progression status at 24 months (logistic, NPr = benefit, Firth if separation); dropping flagged links.
Paired difference from CD8A: 300 within-arm bootstrap resamples (seed [3034, 7]).

## Questions and expectations

| Q | Question | Expectation (written before running) |
|---|---|---|
| Q1 | Seven comparisons (six + NeoPembrOV): signatures at pooled P < 0.05 favoring / against the checkpoint arm | Within ±3 of the six-comparison counts (17 / 1) |
| Q2 | Seven comparisons: signatures significantly above CD8A (paired) | At most 2 of 52 |
| Q3 | Four chemotherapy-backbone comparisons (three breast + NeoPembrOV): signatures at P < 0.05 favoring the checkpoint arm; signatures above CD8A | Fewer than 10 favoring; reported descriptively for CD8A |
| Q4 | NeoPembrOV alone: interaction of CD8A, NetImmune, CAF_ECM and the T-cell panel | Descriptive (53 patients; no prior expectation of direction) |

## Reporting

NeoPembrOV is a sensitivity comparison, not part of the primary six-comparison pool. Whatever Q1–Q4 show is reported in
the main text (one paragraph and a panel or table row) and in full in the Supplementary Information.
