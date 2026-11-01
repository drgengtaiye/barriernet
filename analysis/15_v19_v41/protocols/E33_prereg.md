# E33 protocol: is CD8A a chosen gene, and where does the signal come from?

Written 2026-10-24, after a fifth pre-submission review, before any analysis below was run. Sixth protocol written after
all other results were known; its questions are post hoc and are reported as such. The hash of this file is logged in
E33_prereg.sha256 before E33_trials.py is run; the script stops if the file has changed.

## Prompted by the review

1. CD8A may look chosen: E32 was written after an earlier, unreported analysis (E24) in which CD8A (+0.22) was already
   close to the composite scores. Asked: a reference distribution of single genes (CD8B, CXCL9, GZMB, PRF1 and others)
   and the rank of CD8A in it.
2. "Did not significantly outperform" needs the upper bound of each signature's difference from CD8A and the size of
   difference the comparison can detect.
3. IMmotion150 is a 2,475-gene panel, so 25 signatures are not scored there and are compared with CD8A over different
   comparisons from the others.
4. The renal trials mix anti-angiogenic biology into the interaction; the breast trials (chemotherapy with or without a
   checkpoint inhibitor; I-SPY2 pembrolizumab, I-SPY2 durvalumab–olaparib, NeoTRIP) are the cleanest designs.
5. CYT_Davoli is labelled a resistance signature. Its genes (CD247, CD2, CD3E, GZMH, NKG7, PRF1, GZMK) are cytotoxic
   T-cell genes; CD39-CD8Tcell_Chow is also labelled resistance although its source (Chow 2023, Immunity) reports that
   it predicted benefit. Both favored the checkpoint arm in E32 (+0.17, +0.15).

## Already known before this protocol (declared)

- Every E32 result (analysis/E32_outcome.md), including all per-comparison interactions of the 52 signatures, ESTIMATE,
  TMEscore, CD8A and the barriernet scores (results/E32_per_comparison.csv). In particular CD8A per comparison:
  IMmotion150 +0.53, I-SPY2 pembrolizumab +0.01, I-SPY2 durvalumab–olaparib −0.17, CheckMate 025 +0.21, JAVELIN Renal
  101 +0.17, NeoTRIP +0.16 (so CD8A is near zero in the breast comparisons pooled); barriernet's NetImmune pooled over
  the three breast comparisons with hormone receptor only, +0.63 (P = 0.048).
- The barriernet T_CELL module (CD3D, CD3E, TRAC, TRBC1, CD8A, GZMB, NKG7) and IFN module (STAT1, IRF1, CXCL9, CXCL10,
  IFNG) and CXCL9−SPP1 results. No single gene other than CD8A has been analyzed in these comparisons in any earlier
  analysis (E21, E24, E25, E28, E32).
- The two direction labels above, from looking at E32 results; a full audit of the 52 declared directions against the
  source publications is run in parallel (E33_direction_audit.csv), by an agent that does not see any result.

## Data, models, pooling

As E32 (E32_trials.py, unchanged: common z-score engine, primary model per comparison, I-SPY2 adjusted for hormone
receptor and MammaPrint, pooling rule = REML with the wider of Hartung–Knapp and DerSimonian–Laird intervals; paired
difference from CD8A by 300 within-arm bootstrap resamples per comparison with the E32 seeds, inverse-variance pooled).
A single gene's score is its z-scored expression in the comparison. Genes are used when measured, finite and
non-constant in the comparison.

Panels, fixed here:
- T-cell panel (primary): CD8A, CD8B, CD3D, CD3E, CD3G, CD2, CD247, GZMA, GZMB, GZMH, GZMK, PRF1, NKG7.
- IFN-γ/checkpoint panel (secondary): IFNG, CXCL9, CXCL10, CXCL11, CXCL13, CCR5, STAT1, IDO1, CD274, PDCD1, LAG3,
  TIGIT, CTLA4, HAVCR2.
- All genes: every gene measured in all six comparisons (bounded by the IMmotion150 panel).

## Questions and expectations (written before running)

| Q | Question | Expectation |
|---|---|---|
| Q1 (primary) | Is CD8A an outlier among T-cell genes? Pooled six-comparison interaction of each T-cell panel gene. | CD8A is not ranked first among the 13, and at least 7 of the 13 reach pooled P < 0.05 favoring the checkpoint arm |
| Q2 | Where do T-cell genes sit among all genes? Percentile of CD8A and of the T-cell panel median in the distribution of pooled interactions of all genes measured in all six comparisons | CD8A and the T-cell panel median above the 90th percentile |
| Q3 | Is "no signature significantly exceeds it" specific to CD8A? For each T-cell panel gene, the number of the 52 signatures whose paired difference from that gene is > 0 at P < 0.05 | At most 2 of 52 for at least 7 of the 13 genes |
| Q4 | Upper bounds (descriptive) | Reported: the largest upper 95% limit of signature − CD8A among the 52, the median paired SE and the difference detectable with 80% power at α = 0.05 (2.8 × SE) |
| Q5 | Same comparisons for everyone: paired difference from CD8A over the five comparisons without IMmotion150 | At most 2 of 52 exceed CD8A at P < 0.05 |
| Q6 | Clean designs: the three breast comparisons pooled (chemotherapy with or without a checkpoint inhibitor) — pooled interaction of the 52, CD8A and the T-cell panel; paired difference from CD8A | Reported descriptively (CD8A is known to be near zero there); expectation: fewer than 10 of 52 signatures at P < 0.05 favoring the checkpoint arm |
| Q7 | Declared directions: Q2 of E32 (exceeding the size-matched random reference) recounted with the audited directions | The count changes by at most 3 |

Pooling minimum: a gene or signature is pooled when measured or scored in at least 4 of the 6 comparisons (Q1–Q5) or
in all 3 breast comparisons (Q6). Q5 and Q6 use the E32 random-reference P only where unchanged; they do not draw new
references.

## Rules for reporting

- Whatever Q1–Q3 show is reported. If CD8A ranks first among T-cell genes, the text says CD8A was the strongest T-cell
  gene here and that the comparison with signatures was against the best single gene; if not, the text may say that
  single T-cell genes in general were not outperformed, with the Q3 count.
- The breast-only result (Q6) is reported in the main text whatever it shows.
- Direction corrections (Q7) are a sensitivity analysis; the E32 primary result with SignatureSets labels is kept.
- Q4 values go in the abstract or title wording if the title says "did not significantly outperform".
