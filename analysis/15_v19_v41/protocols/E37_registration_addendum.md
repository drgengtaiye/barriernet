# Addendum to the controlled-access analysis plan for barriernet 0.3 (E37_registration.md)

Written 2026-10-31, after the design analysis (E39) and before any of the controlled-access datasets has been
requested or obtained. It adds a sample-size statement and two reporting rules to E37_registration.md (sha256
c24bb414…), whose primary analysis, comparators and expectations are unchanged.

## What the plan can detect

In the six public randomized comparisons barriernet 0.3 correlated with CD8A at r = 0.81 (within arm, mean over the
comparisons). By normal theory (E39), a 1 : 1 randomized survival comparison detects a difference of δ per SD between
the interactions of barriernet 0.3 and CD8A with 80% power at two-sided 0.05 when it has about
8 (1 − r) × 7.85 / δ² events: about 1,190 events for δ = 0.10, 530 for 0.15 and 300 for 0.20. In the public survival
comparisons the bootstrap variance of the paired difference was 1.17 times the normal-theory value (median), so the
plan assumes about 1,400, 620 and 350 events.

## Reporting rules (added)

1. The minimum detectable difference from CD8A (80% power, two-sided 0.05; 2.80 × the pooled standard error of the
   paired difference) is computed from the data actually obtained and reported with the result.
2. "No meaningful difference from CD8A" is stated only if the 95% confidence interval of the pooled paired difference
   lies within ±0.10 per SD; otherwise a non-significant difference is reported as inconclusive below the minimum
   detectable difference.
3. Pooling across the controlled-access datasets: interaction by the manuscript rule (wider of Hartung–Knapp and
   DerSimonian–Laird); paired differences by fixed-effect inverse variance of the bootstrap variances, as in the
   manuscript. The public six and the controlled-access datasets are reported separately and together.
