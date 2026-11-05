# E42 protocol: barriernet 0.3 against each published signature (paired, six randomized comparisons)

Written 2026-11-01, after all other results were known and before E42_signatures.py is written or run. The hash of
this file is logged in E42_prereg.sha256 first. Fifteenth protocol written after the other results were known.

## Why

E37 compared the frozen barriernet 0.3 only with CD8A and NetImmune 0.2.1 (and, beyond its protocol, with the 27 panel
genes). The manuscript's question for every score is how it stands against the published alternatives. A reader needs
the paired comparison of barriernet 0.3 with each of the 52 library signatures, in the same six comparisons, with the
same bootstrap, so that "no better than CD8A" can be placed against how the published signatures fare.

## Already known (declared; not re-tested as discoveries)

- barriernet 0.3, six comparisons: +0.252 (+0.013 to +0.490); paired difference from CD8A −0.005 (−0.093 to +0.084),
  from NetImmune 0.2.1 −0.016 (E37). Its weights were frozen before any randomized data were used (E37).
- The 52 signatures' six-comparison pooled interactions (E32_primary.csv). The largest is Teff_McDermott, +0.242, so
  barriernet 0.3's point estimate is already known to be the largest; in E32, 15 of 52 were significantly below CD8A
  and none above, without orienting the signatures (E32 Q3 convention); 8 of 52 had a CI for the difference from CD8A
  within ±0.10.
- The ±0.10 rule for "no meaningful difference" (E37_registration_addendum.md) was written after the E37 test results
  were known, for the controlled-access data. Applying it to the public six is post hoc.
- Because barriernet 0.3 correlates with CD8A at r = 0.81 (E39), its paired differences from CD8A-like scores have
  smaller standard errors than those of a typical signature; a CI within ±0.10 is therefore easier for it to reach.

What is not known: the paired differences of barriernet 0.3 from each signature and their bootstrap intervals.

## Data and estimates

E34 setup (E32 data and library, six primary comparisons, NeoPembrOv as a seventh), executed unchanged as in
E37_test.py. barriernet 0.3 from the frozen weights (hash checked). Per signature and comparison, the interaction y
(benefit-oriented, per SD; primary model of each comparison) and 300 within-arm patient bootstrap resamples with the
E32 seeds [3032, comparison] (NeoPembrOv [3034, 7]), common to all scores. Paired difference per comparison:
y(barriernet 0.3) − y(signature); its variance from the bootstrap differences; pooled by fixed-effect inverse variance
over the comparisons where both are scored (at least 4 of the six), as in the manuscript.

Orientation of the signatures (the interaction is odd in the sign of the score, so y and its resamples flip exactly):

1. **Primary, conservative:** each signature in the direction that gave it the larger six-comparison pooled
   interaction in E32 (the sign of its E32 estimate). Every comparator thereby gets its better direction, chosen
   after the fact, which works against barriernet 0.3.
2. Declared direction (sensitive +1, resistance −1; icb/sig_methods_clean).
3. Unoriented (+1 for all; the E32 Q3 convention).

Step 0: the unoriented paired differences of NetImmune 0.2.1 and of all 52 signatures from CD8A over the six must
reproduce E32_primary.csv (diff_vs_CD8A) to < 1e-8, and barriernet 0.3 − CD8A must reproduce E37_test_pooled.csv.

## Reported

- Per signature: paired difference, 95% CI, P, number of comparisons, and the minimum detectable difference (2.80 × SE).
- Counts significantly above (P < 0.05, difference > 0) and below, per orientation; counts with the CI within ±0.10
  (no meaningful difference, the addendum rule, applied post hoc), and the same count for CD8A against each signature
  from E32 for reference.
- Seven comparisons (with NeoPembrOv), primary orientation, as a sensitivity analysis.
- The rank of barriernet 0.3's six-comparison pooled interaction among the 52 (already known; descriptive).

## Expectations

(A1, primary) In the conservative orientation, barriernet 0.3 is significantly above at least 5 of the 52 signatures
and significantly below none.
(A2) In the declared direction, it is significantly above at least 10 of the 52 and below none.

## Reading rules

- If A1 holds: the text may say that barriernet 0.3 was significantly better than N of the 52 published signatures,
  worse than none, and within ±0.10 of CD8A, stating that the comparators were oriented in their better direction, that
  the margin rule was applied post hoc, and that a score correlated with CD8A reaches a narrow interval more easily.
- If A1 fails: the count is reported as it is, and the text keeps "no better than CD8A".
- In no case does the text say that barriernet 0.3 is better than CD8A or than the best signatures.
