# E40 protocol: does the usual validation track the randomized interaction? (real-data counterpart of E36)

Written 2026-10-31, before E40_export.R and E40_concordance.py are written or run. The hash of this file is logged in
E40_prereg.sha256 first.

## Why

The simulation (E36) showed that a purely prognostic marker passed the interaction test in 1.1% of trials but appeared
predictive in 94.9% when the checkpoint arm was analyzed alone. This analysis asks the same question of the 52 published
signatures in real data: how often the checkpoint arm alone, and the usual validation in single-arm checkpoint cohorts,
point to a predictive signature, and how well they rank the signatures by their randomized interaction.

## Data

- Randomized: E34 setup (E32 data, six primary comparisons), unchanged. Per signature: the pooled interaction (E32
  primary, reproduced) and the pooled within-arm association in the checkpoint arm and in the control arm (E28 m_within;
  per SD, benefit-oriented; pooled over the comparisons by the manuscript rule, scored in ≥ 4 of the 6).
- Single-arm: the 20 single-arm checkpoint cohorts of the barriernet 0.3 development set (E37; ICB_Braun excluded; the
  three ICB_Miao1 samples removed by the twin check), loaded and z-scored exactly as E37_export.R, now for all genes of
  the 52 signatures. Each signature scored with the manuscript engine (sum(w z) / sum(|w|) over the genes present; ≥ 3
  genes and ≥ 50% of the list). Per cohort: logistic regression of response on the standardized score (Firth on
  separation); log-OR per SD pooled by DerSimonian–Laird with the manuscript rule (wider of Hartung–Knapp and DL),
  scored in ≥ 10 of the 20 cohorts.
- Step 0: the 48 barriernet gene z-scores of the new export must equal E37_dev_genes.csv (< 1e-10), and the E32 pooled
  interactions must be reproduced (< 1e-8).

All associations use the raw score (higher score), not the declared direction.

## Analyses

1. Counts of signatures at P < 0.05 (either direction): interaction; checkpoint arm alone; control arm alone;
   single-arm cohorts. Cross-tabulation of checkpoint-arm-alone against interaction and against the control arm.
2. Spearman correlation across the signatures between the interaction and (a) the checkpoint-arm-alone estimate,
   (b) the control-arm estimate, (c) the single-arm pooled log-OR. Bootstrap over signatures is not used (signatures are
   correlated); the correlations are descriptive.
3. Of the signatures significant in the single-arm cohorts, the number with a significant interaction in the same
   direction.

## Expectations

- C1: more signatures reach P < 0.05 in the checkpoint arm alone than in the interaction.
- C2: the Spearman correlation between the single-arm pooled log-OR and the interaction is positive and below 0.8.
- C3: fewer than 60% of the signatures significant in the single-arm cohorts have an interaction at P < 0.05 in the
  same direction.

Reported whatever the result.
