# E35 protocol: single genes as the reference in the breast comparisons

Written 2026-10-26, before any of the analyses below were run. Eighth protocol written after all other results were
known; post hoc. Hash logged in E35_prereg.sha256 before E35_trials.py is run.

## Why

A sixth pre-submission review: the clean-design half of the evidence (three breast comparisons, chemotherapy with or
without a checkpoint inhibitor) is where six signatures nominally exceeded CD8A, so the title may rest on the confounded
renal half. It asked, within the breast comparisons, how many signatures significantly exceed each single gene: if
signatures exceed only CD8A and not genes such as GZMB or CXCL10, the title holds there too; if they also exceed those,
the title must be limited to "across trials".

## Already known (declared)

All E32–E34 results, in particular for the three breast comparisons: pooled interaction of every panel gene (CD8A
+0.06; GZMB +0.49, IDO1 +0.47 and CXCL10 +0.45 at P < 0.05; STAT1 +0.45, PRF1 +0.35, NKG7 +0.30 not significant);
the six signatures significantly above CD8A (Bcell_Budczies, CTL_Jiang, IFN_Ayers, Teff-IFNG_Fehrenbacher,
Teff_McDermott, VIGex_Hernando-Calvo); and, with NeoPembrOV added (four chemotherapy-backbone comparisons), six above
CD8A. Paired differences between signatures and genes other than CD8A have not been computed in these comparisons.

## Method

As E33 Q3, restricted to the three breast comparisons (and, as a sensitivity, the four with NeoPembrOV): for each of
the 27 panel genes and each of the 52 signatures, the paired difference of interactions, its variance from 300
within-arm bootstrap resamples per comparison with the E32/E33 seeds (NeoPembrOV: the E34 seed), inverse-variance
pooled; "exceeds" = difference > 0 at P < 0.05. Genes must be measured in all comparisons of the set (CXCL9 is not
measured in I-SPY2 pembrolizumab and drops out of the breast set). Step 1 reproduces E33's breast differences from CD8A.

## Questions and expectations

| Q | Question | Expectation |
|---|---|---|
| Q1 (primary) | Breast: signatures significantly above GZMB, IDO1 and CXCL10 (the single genes at P < 0.05 there) | At most 2 of 52 for each of the three |
| Q2 | Breast: T-cell genes exceeded by more than 2 signatures | At least 5 of 13 (several T-cell genes are as weak as CD8A there) |
| Q3 | Breast: the gene exceeded by the fewest signatures and that number | Descriptive |
| Q4 | Four chemotherapy-backbone comparisons (with NeoPembrOV): Q1 repeated | At most 2 of 52 for each of GZMB, IDO1, CXCL10 |

## Rules for reporting

- If Q1 holds: the text says that in the breast comparisons signatures exceeded CD8A and other weak T-cell genes but
  not the strongest single genes there, so "no signature outperformed a single gene" holds in the clean designs when
  the gene is chosen from the panel, while the best gene differs between trial types.
- If Q1 fails for any of the three genes: the title and abstract are limited to "across trials" and the breast result is
  stated as signatures outperforming single genes there.
- Reported in the main text (Fig. 3f and one paragraph) whatever it shows.
