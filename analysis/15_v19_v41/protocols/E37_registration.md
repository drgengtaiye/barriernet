# Registration of barriernet 0.3 for controlled-access randomized trials

Written 2026-10-31, after barriernet 0.3 was frozen (results/E37_frozen_weights.tsv, sha256 500264cc627b3ef0240c0340296dc001bceb02748e3f3b0e83892b234b2f63bd) and before any of
the datasets below has been requested or obtained. Protocol: E37_prereg.md. The frozen version was developed only on
20 single-arm checkpoint cohorts (1,195 patients) and selected by leave-one-cohort-out validation.

## Frozen score

Genes z-scored across the samples of a dataset (log-scale expression); score = sum(w z) / sum(|w|) over the genes
present, scored when at least three genes and half of the 48 are present. Each gene of a module has the same weight:

| Module | Genes | Weight per gene | Module total |
|---|---|---|---|
| IFN | 5 | +0.07234 | +0.3617 |
| APM | 8 | +0.02814 | +0.2251 |
| TLS_B | 6 | +0.02102 | +0.1261 |
| T_CELL | 7 | +0.01488 | +0.1042 |
| CAF_ECM | 11 | -0.01452 | -0.1597 |
| TAM_SUP | 11 | -0.00211 | -0.0232 |

Gene lists are those of barriernet 0.2.0 (48 genes). barriernet 0.3 ships these weights (geneset_rna_v03.yaml).

## Datasets (controlled access; none obtained)

OAK and POPLAR (atezolizumab vs docetaxel, non-small-cell lung cancer; EGAS00001005013), IMvigor010 (adjuvant
atezolizumab vs observation, urothelial), IMmotion151 (atezolizumab–bevacizumab vs sunitinib, renal), CheckMate 649
(nivolumab–chemotherapy vs chemotherapy, gastric/esophageal; EGAS50000000747), and any other randomized checkpoint trial
with pre-treatment expression from both arms obtained later.

## Primary analysis (each dataset)

Cox model of the trial's primary survival endpoint (overall survival for OAK/POPLAR and CheckMate 649; disease-free
survival for IMvigor010; progression-free survival for IMmotion151) on score, arm and score × arm, score standardized,
benefit-oriented interaction = −log HR per SD, with the trial's stratification factors as covariates where available.
Comparators by the paired within-arm bootstrap (300 resamples): CD8A, NetImmune 0.2.1, and the published signatures of
the manuscript's library; against 1,000 size-matched random gene sets (weights in random order).

Expectations: interaction > 0; not significantly different from CD8A. Reported whatever the result, with this file.
