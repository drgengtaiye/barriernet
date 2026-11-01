# E39 outcome (protocol E39_prereg.md, sha256 1a8f51cf…, logged 2026-11-01T01:41:29Z)

Run: `E39_design.py` (`results/E39_design.log`), after the hash. Step 0: the 52 pooled paired differences from CD8A and
their standard errors reproduced E33_paired.csv (max |difference| 9.6 × 10⁻¹⁷ and 9.7 × 10⁻¹⁷).
Outputs: `results/E39_per_comparison.csv`, `E39_pooled.csv`, `E39_required.csv`, `E39_expectations.csv`.

| | expectation | result | held |
|---|---|---|---|
| X1 | median ratio of bootstrap to normal-theory variance in [0.67, 1.5] | 1.45 (IQR 1.16–1.85; 287 signature × comparison pairs) | yes |
| X2 | median over comparisons of Spearman(bootstrap, theory) ≥ 0.8 | 0.96 (range 0.90–0.99) | yes |

## Findings

- The formula ranks the uncertainty almost perfectly but underestimates it. The ratio of bootstrap to normal-theory
  variance was 1.17 in the survival comparisons and 1.76 in the pCR comparisons (median). In I-SPY2
  durvalumab–olaparib, with 29 controls, it was 9.9. X1 held only narrowly. Planning should multiply the theoretical
  events by about 1.2 for survival endpoints.
- Current minimum detectable difference from CD8A (80% power): median 0.21 per SD over the 52 signatures (10th–90th
  percentile 0.10–0.32). barriernet 0.3: 0.13. NetImmune 0.2.1: 0.18. The more a score correlates with CD8A, the
  smaller the difference that can be detected: barriernet 0.3 r = 0.81; NetImmune 0.68; median signature 0.55.
- Events needed in a 1 : 1 survival comparison for δ = 0.10 per SD (normal theory):
  - 2,825 at the median signature's correlation (0.55);
  - 1,186 at barriernet 0.3's correlation (0.81);
  - 4.5 times the information of the six current comparisons for the median signature, 1.6 times for barriernet 0.3.
  - For comparison, the three current survival comparisons have 744 events in all.
- Added to the analysis plan as `E37_registration_addendum.md` (sha256 c07e4e8c…, logged 2026-11-01T01:49:45Z, before
  any controlled-access dataset was requested). The addendum sets the target at about 1,400 events for δ = 0.10, and
  adds two rules: the minimum detectable difference is reported with the result, and "no meaningful difference" is
  claimed only if the 95% CI lies within ±0.10.

Later-protocol tally after E39: 11 protocols, 35 of 45 expectations held.
