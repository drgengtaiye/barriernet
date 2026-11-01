# E31 outcome (protocol: E31_prereg.md, sha256 0e7c1008…, logged 2026-10-22T14:11:18Z)

Run 2026-10-22 (`E31_dedup_all_pairs.R`; a first run completed but its session was interrupted, and the rerun's log is
identical apart from the time stamp: `results/E31_interrupted_run1.log`, `results/E31.log`). Action re-run
`E31_rerun.R` (`results/E31_rerun.csv`, `E31_rerun.log`, `E31_removed_arrays.csv`).

## Expectation: NOT HELD

*No pair other than GSE26253–GSE62254 contains an array pair that meets rules i–iii.* 28 array pairs in 23 of the
other 44 cohort pairs met the rule (at most 2 per cohort pair).

## What the screen shows

| | GSE26253–GSE62254 | other 44 pairs |
|---|---|---|
| mutual best matches agreeing on sex and age (i–ii) | 11 | 64 (at most 7 per pair, GSE84437–GSE62254) |
| … of which the best match stands > 4 SD above its row | 9 | 2 |
| meeting rules i–iii | 9 (the eight E30 duplicates and GSM644832, which E30's stricter time rule excluded) | 28 |

- **GSE15459–GSE34942** (the reviewer's pair): 22 mutual best matches, none agreeing on sex and age; under the E20
  99.9th-percentile rule 8 GSE15459 arrays exceed GSE15459's threshold (0.574) against GSE34942 arrays, all matching
  a few hub arrays, and no GSE34942 array exceeds its own threshold (0.704). No duplicate.
- Post hoc reading (not part of the protocol): about one mutual best match per cohort pair agrees on sex and age, and
  rule iii (which lets any two censored records pass) keeps roughly half, so the rule finds chance matches when applied
  to unrelated cohorts. The GSE26253–GSE62254 matches differ in kind: 9 of 11 stand more than 4 SD above their row
  (2 of 64 elsewhere), and three have identical recurrence and survival times. This supports E30's removals and
  shows the rule is not specific outside a pair already known to share a hospital.
- Correction: the "already known" item of this protocol (taken from E30_prereg.md and repeated in the v32 SI) said
  that at most one mutual best match agreeing on sex and age occurs between any other two cohorts. That was wrong
  (up to 7); it holds only for matches standing more than 4 SD above their row. The SI text is corrected in v33.

## Action (as specified): remove and re-run, reported beside the current values

27 arrays removed (28 pairs; GSM662606 of GSE26901 was already excluded by the E20 rule): GSE13861 3, GSE15459 2,
GSE26253 2, GSE26899 6, GSE26901 4 (+1 already out), GSE29272 3, GSE34942 2, GSE62254 4, GSE84437 1; 2,144 → 2,117
patients. Within the benchmark group "later" = later in the protocol's cohort list (fixed before the re-run).

Implementation (`E31_rerun.R`): current patient-level scores reused, Cox refitted, manuscript pooling rule. **Deviation
from the protocol**: the protocol said "re-run with the E30 machinery", which recomputes the scores; here the scores
are not recomputed without the removed arrays. Step 0 reproduced 170 current per-cohort and 80 current pooled values
(largest difference 4.4 × 10⁻¹⁵). Null ranks and C not re-run.

Scope: the first version of the re-run covered 13 pooled estimates; the v33 fact-check found that the added-four pools,
the adjusted pools and the module pools of the main text were missing, and the re-run was extended to all 104
(benchmark 3, first three, added four, locked 7, all 10; NetImmune, CAF_ECM, TAM_SUP, Barrier, Immune, APM, IFN, TLS_B,
T_CELL, CD8A; unadjusted and, where every cohort has it, adjusted; Immune and CAF_ECM in one model; direct test).

- Largest change of a pooled HR: 0.036 (first three, TAM_SUP, not significant); 0.028 for NetImmune, CAF_ECM, Immune
  and Barrier; 0.017 for any NetImmune pool (added four: 0.67 → 0.65, P 0.012 → 0.020).
- Significance changed for 5 of 104: first three NetImmune (0.049 → 0.075), first three TLS_B (0.033 → 0.070), added
  four APM adjusted (0.058 → 0.037), added four TLS_B adjusted (0.064 → 0.041), locked 7 IFN adjusted (0.037 → 0.065).
  Only the first is reported in the main text.

| set | quantity | current | E31 removals |
|---|---|---|---|
| benchmark 3 | NetImmune | 0.74 (0.58–0.95); P = 0.036 | 0.74 (0.58–0.95); P = 0.036 |
| locked, first three | NetImmune | 0.72 (0.53–0.996); P = 0.049 | 0.73 (0.50–1.08); **P = 0.075** |
| locked, added four | NetImmune | 0.67 (0.53–0.85); P = 0.012 | 0.65 (0.49–0.88); P = 0.020 |
| locked 7 | NetImmune | 0.70 (0.62–0.79); P = 2.9 × 10⁻⁴ | 0.70 (0.60–0.81); P = 0.0011 |
| locked 7 | CAF_ECM | 1.34 (1.09–1.65); P = 0.014 | 1.35 (1.09–1.69); P = 0.015 |
| all 10 | NetImmune | 0.72 (0.67–0.78); P = 7.3 × 10⁻⁶ | 0.72 (0.66–0.79); P = 2.2 × 10⁻⁵ |
| all 10 | CAF_ECM | 1.31 (1.17–1.46) | 1.31 (1.17–1.47) |
| all 10 | Immune / CAF_ECM, one model (post hoc) | 0.82 / 1.33 | 0.82 / 1.33 |
| all 10 | direct test | −0.08 (−0.22 to 0.06); P = 0.22 | −0.08 (−0.22 to 0.05); P = 0.20 |

Full table: `results/E31_rerun.csv` (Supplementary Table S19j). The main text keeps the current values (the protocol
reports the removals beside them; most removed matches are chance, above) and states this in Methods and Table 1.

## Tally

Later protocols: 6 of 11 expectations held (E28–E30: 6 of 10; E31: 0 of 1).
