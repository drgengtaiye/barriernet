# E46b addendum to the E46 protocol: is the I-SPY2 flag a shared-patient signal?

Written 2026-11-05 after E46 part D flagged 15 of the 29 ISPY2_durva controls by the E37 rule, with 0 flagged in the reverse direction
and no shared identifiers. Post hoc; no expectations; reads no outcome.

The E37 rule compares whole-profile Spearman correlations, which gene means dominate; one sample with a high correlation to many others
can produce flags without any duplicate. Three checks:

1. Negative control. The durva-olaparib-treated samples of ISPY2_durva (GSE173839) cannot be patients of the pembrolizumab comparison.
   Flag rate of those samples against all samples of ISPY2_pembro, set beside the flag rate of the controls; and the number of distinct
   partners among flagged samples and of mutual best pairs.
2. Gene-centered profiles. Within each comparison each gene is centered on its mean over that comparison, and control profiles are
   correlated across the shared genes (Pearson). Independent tumors stay well below a duplicate, whose profiles should correlate above 0.8
   whatever the deposit's processing. A control is called shared if its best partner exceeds 0.8.
3. Positive control. The same statistic for each ISPY2_pembro sample against a copy of ISPY2_pembro processed differently: a monotone
   per-sample transform plus Gaussian noise with SD 0.3 per gene, to show that a duplicate would be detected.

Outputs: results/E46b_ispy2.csv, results/E46b_ispy2.log. If any control is called shared, the Supplementary Methods statement that the comparisons
shared no patients is corrected and the sensitivity analysis of E46 part D is the one reported.
