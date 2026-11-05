# E46 protocol: a simple T-cell reference score, choosing a gene in advance, the MSS/EBV-negative gastric patients, the I-SPY2 controls and distinct signature counts

Written 2026-11-05, after all other results were known and before E46_test.py or E46_distinct.py is written or run.
Nineteenth protocol written after the other results were known. All analyses are post hoc.

## Why

A second external review of the manuscript raised five points that the existing analyses do not answer:

1. The third reference: comparing a signature with every gene of a panel favors the panel, because some gene will be high by
   chance. The review proposes a simple equal-weight T-cell score as the reference (it names CD8A, PRF1 and GZMA as an example).
2. Whether a gene that matches the signatures can be identified in advance, and whether aggregated scores are steadier across
   tumor types than single genes.
3. Discussion of the MSS/EBV-negative gastric patients rests on 40 patients of one cohort and on barriernet scores only.
4. Whether the control arms of the two I-SPY2 comparisons overlap.
5. TIS and the T-cell-inflamed list are the same 18 genes and are counted twice.

## What is already known

Every result of E23-E45, including the interaction of every panel gene in the six comparisons (CD8A +0.23, PRF1, GZMA and the
others), so the genes of the reference score below were not chosen blind to the randomized data. CXCL13 and STAT1 had the
highest panel interactions (+0.28 each). TIS_Damotte and Tcell-inflamed_Ayers give identical paired differences in
E45_A_paired.csv. While checking the review (2026-11-05) the author's assistant also looked at: E20_pooled.csv (NetImmune in the
seven locked cohorts), E44_B_pooled.csv (CD8A by the two pooling rules), and the counts in E44_labels.csv (PRJEB25780: 33
MSS/EBV-negative patients, 2 responders; 11 MSI-H or EBV-positive, 10 responders. PRJEB60680: 40, 26 responders; 3, 3 responders).

Not yet looked at: any score built from several panel genes as a reference; any leave-one-comparison-out choice of gene; any score
within the MSS/EBV-negative patients other than the Immune axis and barriernet modules in PRJEB60680; any expression comparison of
the two I-SPY2 control arms.

## A. Simple T-cell reference scores

Setup as E45 (E34_trials.py up to its fitting stage). Two references, each the equal-weight mean z-score (`wscore`, weights 1) of
its genes within a comparison:

- R3 (primary): CD8A, PRF1, GZMA, as named by the review.
- R13 (secondary): the 13 genes of the T-cell panel, a list fixed before the panel analyses, with no gene selected by result.

For R3, R13, CD8A and CXCL10: interaction per SD by the main-text rule in the six comparisons (primary), the three renal and the
three breast comparisons, and the percentile of the six-comparison estimate among the 2,082 genes (results/E33_pooled.csv). For
each of the 52 signatures, NetImmune and barriernet 0.3: score minus reference by the paired bootstrap (300 resamples; the seeds
of E32/E42/E45; fixed-effect inverse variance of the bootstrap variances). Counted over the 51 distinct signatures (see E):
significantly above, significantly below (P < 0.05), interval within ±0.10 per SD, median detectable difference (2.8 × SE). Also
R3 minus CD8A, R3 minus CXCL10, R13 minus CD8A and R13 minus R3. Step 0: the 52 signatures and NetImmune minus CD8A reproduce
results/E32_primary.csv (< 1e-8).

Steadiness across tumor type (descriptive): D = |renal pool − breast pool| by the main-text rule, for the 25 panel genes measured
in all six comparisons, the 19 distinct inflamed signatures, R3 and R13. Median D of the panel genes against the inflamed
signatures, Mann-Whitney two-sided. The pools are on different scales (log HR and log OR), which applies to every score alike.

## B. Choosing a gene in advance (leave one comparison out)

For each comparison h in turn: among the 25 panel genes measured in all six, G*(h) is the gene with the highest pooled interaction
over the other five comparisons (main-text rule). Its interaction in h is recorded. The six held-out values are pooled by the
main-text rule and compared with the apparent best, the highest six-comparison pool among the same 25 genes. Paired: G*(h) minus
R3 and G*(h) minus the T-cell-inflamed list in comparison h, by the bootstrap of part A, fixed-effect inverse variance over h.
Secondary: the same with the choice restricted to the 13 T-cell genes.

## C. MSS/EBV-negative gastric patients

PRJEB25780 and PRJEB60680 from results/E44_Z_*.csv.gz and E44_labels.csv, restricted to msi_ebv = 0 (73 patients, 28 responders).
Scored as E44 (the whole-cohort z-scores of E44_Z), restandardized within the subgroup, Firth logistic regression per cohort as E44.
Scores: the 52 signatures, ESTIMATE, TMEscore_z, Immune, NetImmune, barriernet 0.3, the panel genes, R3 and R13. Pooled by fixed
effect (as E44) and by the main-text rule; PRJEB60680 alone is also reported, because PRJEB25780 contributes 2 responders. Counted
over the 51 distinct signatures: positive, P < 0.05 by each pooling. Step 0: CD8A in all patients of each cohort reproduces
results/E44_B_per_cohort.csv (< 1e-6).

## D. I-SPY2 control arms

Rule as E37 and E12: for each control of ISPY2_durva, Spearman correlation with every control of ISPY2_pembro over the shared genes
(at least 200); flagged when the best partner stands more than 0.03 clear of the second best. Also the same over all samples of the
two comparisons, and the overlap of sample identifiers. Reads no outcome. If any control is flagged, the flagged controls are
removed from ISPY2_durva and the six-comparison pools of CD8A, R3 and the T-cell-inflamed list are repeated as a sensitivity
analysis, and the statement in the Supplementary Methods that the comparisons shared no patients is corrected.

## E. Distinct signature counts

E46_distinct.py recomputes, from the existing results files and those of parts A-C, every count of signatures in the main text
with the row TIS_Damotte removed, so that the two rows with the same 18 genes count once (51 distinct signatures: 33 sensitivity,
18 resistance; 19 inflamed, 8 confounder, 24 other). The row kept is Tcell-inflamed_Ayers, whose name the main text already uses.
Checked first: the two rows have identical estimates in the six comparisons. No new model is fitted. The source-paper review
and the direction audit stay at 52 entries, because each entry has its own source paper and label.

## Expectations (written before the run; reported whatever the result)

| id | expectation |
|---|---|
| A1 | R3 favors the checkpoint arm in the six comparisons at P < 0.05 |
| A2 | at most 2 of the 51 distinct signatures are significantly above R3 in the six comparisons |
| A3 | R3 has a larger breast interaction than CD8A (+0.06) and a smaller renal interaction than CD8A |
| A4 | the median D of the inflamed signatures is below that of the panel genes |
| B1 | the pooled held-out interaction of the gene chosen in advance is below the apparent best of the 25 genes |
| B2 | the gene chosen in advance does not significantly exceed R3 (paired P >= 0.05) |
| C1 | CD8A has a positive fixed-effect log-odds in the 73 patients, with P >= 0.05 |
| C2 | fewer than 5 of the 51 distinct signatures reach P < 0.05 by fixed effect in the 73 patients |
| D1 | no I-SPY2 control is flagged as shared between the two comparisons |

If more than two signatures exceed R3, the statement that the signatures do not outperform simple T-cell scores changes accordingly.

Outputs: results/E46_*.csv, results/E46_test.log, E46_outcome.md.
