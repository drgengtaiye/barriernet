# E36 outcome: operating characteristics of the calibration framework

Protocol E36_prereg.md (sha256 4125147210c4…, logged 2026-10-31T02:59:12Z, before E36_sim.py was written).
Script E36_sim.py; results results/E36_S0_S2.csv, E36_S3.csv, E36_S4.csv, E36_expectations.csv, E36_S4_power_vs_theory.csv.

## Deviations from the protocol text (all declared)

1. Fast Cox check. Against lifelines at its default convergence precision the largest difference was 1.3e-6, above
   the protocol's 1e-6. The reference, not the threshold, was changed: with lifelines' precision set to 1e-12 the
   largest difference was 2.5e-8, and the check passed. Both values are printed by the script.
2. Missing fits. The first run stopped in S3 on a singular fit (random single genes as the reference of CD8A). As E32
   does for reference draws, a non-converging or singular fit is now recorded as missing and pooled over the remaining
   comparisons. S0–S2 were identical in both runs (no failed fit); in S3, 1,004, 958 and 903 of 240,000 random
   single-gene reference fits failed (0.4%), and none for the multi-gene markers.
3. The 200 random gene sets per comparison and size are drawn once and shared by all replicates of a scenario.

## Expectations: 3 of 6 held

| | Expectation | Result | Held |
|---|---|---|---|
| E1 | S0: interaction rejection 2.5%–7.5%, each marker | 1.7%–2.5% (CD8A 2.5%, NetImmune 1.7%) | No — the pooling rule is conservative |
| E2 | S1: CAF_ECM interaction ≤ 7.5%; checkpoint arm alone ≥ 50% | 1.1%; 94.9% | Yes |
| E3 | S2: CD8A power increasing, ≥ 60% at 0.3 | 16.7%, 58.9%, 85.9% at 0.1, 0.2, 0.3 | Yes |
| E4 | S3: null exceedance ≤ 10% each; CD8A truth: T_CELL ≥ 50% | null 6.0%–10.5% (T_CELL 10.5%); T_CELL 96.5% | No — T_CELL 10.5% |
| E5 | S4 truth CD8A: a module significantly above CD8A ≤ 5% | 0%–0.5% | Yes |
| E6 | S4: power ≥ 50% where oracle ≥ 0.20; ≤ 30% where ≤ 0.10 | no truth gave an oracle difference ≥ 0.20 (largest APM +0.112); TLS_B +0.095 had 36% | No |

## Other results (descriptive)

- S2: estimates were unbiased (mean +0.102, +0.200, +0.301 for true 0.1, 0.2, 0.3). Correlated markers recovered part
  of the CD8A truth (beta 0.3: T_CELL 79%, IFN 70%, NetImmune 57%, APM 47% power); CAF_ECM 2.7%.
- S3 (CD8A truth): every immune marker exceeded its random reference in 77%–98.5% of replicates; CAF_ECM, in its
  declared (adverse) direction, in 20.5%. Diffuse truth (first principal component, correlated with CD8A 0.0–0.5 by
  comparison): immune markers exceeded their reference in 0.5%–5% of replicates although their interactions were
  positive (+0.05 to +0.06 for CD8A, T_CELL, IFN) — the random reference declined to call a transcriptome-wide signal
  specific; CAF_ECM exceeded in 29.5%.
- S4 (truth CD8A): correlated modules were significantly below CD8A in 17%–48% of replicates and above in 0%–0.5%.
  With a module as the truth, the paired test's power matched the normal-theory prediction Φ(Δ/SE − 1.96) from the
  oracle difference and the median paired SE (post hoc comparison): T_CELL 10.0% vs 9.8%, IFN 25.5% vs 25.5%, APM
  44.5% vs 43.4%, TLS_B 36.0% vs 36.3%. The paired SE ranged from 0.029 (T_CELL, which contains CD8A) to 0.062 (APM),
  so the detectable difference depends on how closely a score tracks the reference gene: 0.08 to 0.17 per SD for these
  modules, against 0.21 for the median published signature.
