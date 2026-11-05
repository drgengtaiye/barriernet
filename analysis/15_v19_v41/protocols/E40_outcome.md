# E40 outcome (protocol E40_prereg.md, sha256 0a9f606e…, logged 2026-11-01T01:41:29Z)

## Runs (after the hash)
- `E40_export.R` (`results/E40_export.log`).
  - Step 0: the 48 barriernet gene z-scores equal E37_dev_genes.csv (max |diff| 5.3 × 10⁻¹⁵).
  - Deviation (implementation only): the first run stopped because the missing-value check used `identical()`, which also compared row names. It was replaced by an element-wise check, under which no element differed, and re-run.
- `E40_concordance.py` (`results/E40_concordance.log`).
  - Step 0: the pooled interactions reproduce E32_primary.csv (9.7 × 10⁻¹⁷).

## Data
- Single-arm cohorts: 20 cohorts, 1,195 patients, 368 responders (the three twin-flagged samples removed).
- All 52 signatures were pooled in ≥ 10 cohorts.
- Outputs: `results/E40_randomized_per_comparison.csv`, `E40_single_arm_per_cohort.csv`, `E40_concordance.csv`, `E40_expectations.csv`.

## Expectations

| | expectation | result | held |
|---|---|---|---|
| C1 | more signatures at P < 0.05 in the checkpoint arm alone than in the interaction | 7 vs 18 | **no** |
| C2 | Spearman(single-arm log-OR, interaction) > 0 and < 0.8 | 0.83 | **no** |
| C3 | < 60% of single-arm-significant signatures have a same-direction interaction at P < 0.05 | 17 of 34 (50%) | yes |

## Findings

### Checkpoint arm and control arm
- The checkpoint arm alone did not mislead for these signatures. No signature was associated with outcome in the control arms (0 of 52 at P < 0.05; Spearman of control-arm estimate with interaction 0.35).
- The failure mode shown by the simulation, a prognostic marker read as predictive, needs a marker that is prognostic in the control arm, and none of the 52 was.
- Pooled within-arm associations were less often significant than the interaction: 7 vs 18. Four of the 7 were not significant in the interaction, and none of those was significant in the control arm.
- Checkpoint-arm-alone estimates ranked the signatures as the interaction did (Spearman 0.81).

### Single-arm cohorts
- Single-arm cohorts ranked the signatures much as the randomized interaction did (Spearman 0.83). Seventeen of the 18 signatures with an interaction at P < 0.05 were also significant in single-arm cohorts.
- Single-arm significance overstated the evidence: half of the 34 signatures significant in single-arm cohorts (17) had no significant interaction. This partly reflects power: the single-arm set has 1,195 patients and the interaction is a difference.
- The one signature significant in the interaction but not in single-arm cohorts was IRG_Yang (+0.23, P = 0.022; single-arm +0.11, P = 0.21).
- IMPRES_Auslander was significant in both, in the same (adverse) direction.

### Reading
Single-arm validation ranks published immune signatures roughly as randomized interaction does, but it does not establish that a given signature is predictive. For these signatures the specific danger of prognostic confounding did not materialize in the control arms. The simulation's 94.9% applies to markers that are prognostic in the control arm (in the simulation, the fibroblast module made prognostic by construction), not to the immune signatures tested here. Both expectations that failed were more pessimistic about single-arm evidence than the data proved to be.

Later-protocol tally with E38 (3 of 4 held): 13 protocols, 39 of 52 expectations held.
