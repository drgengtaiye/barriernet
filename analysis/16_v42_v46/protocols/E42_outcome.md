# E42 outcome: barriernet 0.3 against each published signature

Protocol E42_prereg.md (sha256 1101490457b3…, logged 2026-11-02T01:24:04Z, before E42_signatures.py was written).
Fifteenth protocol written after the other results were known. Run: E42_signatures.py → results/E42.log,
E42_paired.csv, E42_counts.csv, E42_counts_by_family.csv, E42_expectations.csv.

Step 0: the unoriented paired differences from CD8A of the 52 signatures and NetImmune 0.2.1 reproduced E32 (largest
difference 9.6 × 10⁻¹⁷), and barriernet 0.3 − CD8A reproduced E37 (6.9 × 10⁻¹⁷). barriernet 0.3 was scored in all
seven comparisons; IMmotion150's 2,475-gene panel scored 27 of the 52 signatures.

## Expectations: 2 of 2 held

| | expectation | result | held |
|---|---|---|---|
| A1 (primary) | conservative orientation: significantly above ≥ 5 of the 52, below none | above 6, below 0 | yes |
| A2 | declared direction: significantly above ≥ 10, below none | above 16, below 0 | yes |

## Results (six comparisons)

| orientation | barriernet 0.3 above | below | CI within ±0.10 | CD8A above (same signatures) | CD8A within ±0.10 |
|---|---|---|---|---|---|
| conservative (primary) | 6 | 0 | 11 | 4 | 8 |
| declared | 16 | 0 | 11 | 14 | 7 |
| unoriented (E32 convention) | 16 | 0 | 11 | 15 | 8 |
| seven comparisons, conservative | 8 | 0 | 11 | 3 | 9 |

- Significantly above, conservative orientation: IFN_Ayers, APM_Wang, ImmuneScore_Roh, MHC-I_Liu, TLS_Cabrita, IMS_Cui.
  CD8A was above the last four of these.
- CI within ±0.10 (no meaningful difference by the addendum rule, applied post hoc): eleven signatures, all among the
  15 with the largest pooled estimates (two of them in the library's "Other" family; IMPRES, IRG_Yang, APM_Thompson
  and Bcell_Budczies of the top 15 were outside) (Teff_McDermott, CTL_Jiang, GeparSixto_Karn, TIS_Damotte, Tcell-inflamed_Ayers, VIGex,
  Inflammatory_Thompson, STAT1_Care, CCL5-CXCL9_Dangaj, Teff-IFNG_Fehrenbacher, IRG_Ayers), and CD8A
  (−0.005, −0.093 to +0.084).
- Median minimum detectable difference from a signature: 0.23 per SD (range 0.09–0.45); CD8A against the same
  signatures 0.22.
- barriernet 0.3's six-comparison pooled interaction (+0.252) ranks first of 53 with every signature in its better
  direction (known before the protocol); the nearest is Teff_McDermott, +0.242, paired difference −0.018 (−0.097 to
  +0.060).
- Margins: A1 held with one to spare, and one of the six, IMS_Cui, has an E32 estimate of +0.001, so its orientation is
  essentially arbitrary. The orientation by the sign of the E32 estimate is not always the one least favorable to
  barriernet 0.3 in the paired test (CRMA_Shukla, MHC-II_Liu); neither changes a count.
- The declared-direction counts include signatures whose library label is "resistance" although they are immune
  scores (for example CYT_Davoli and CD39-CD8Tcell_Chow); the conservative orientation avoids that artefact and is the
  one to quote.

## Reading

barriernet 0.3 was significantly better than 6 of the 52 published signatures (8 of 52 with NeoPembrOv), worse than
none, and within ±0.10 per SD of CD8A and of 11 of the 15 signatures with the largest pooled estimates, with every comparator in its
better direction. CD8A did almost as well against the same signatures (4 of 52 above). The paired comparison
therefore places barriernet 0.3 in the leading group of published scores, alongside CD8A; it does not separate it
from that group. Its narrower intervals partly reflect its correlation with CD8A-like scores (r = 0.81 with CD8A).

## Independent audit (2026-11-01)

The script was found faithful to the protocol; every number and both expectations were reproduced from the CSVs.
Additions beyond the protocol: counts by family, the median minimum detectable difference, barriernet 0.3 − NetImmune
0.2.1, and the CD8A reference recomputed in every orientation and in the seven-comparison set (E32 holds only the
unoriented six; those match). kmin = 4 was also used for the seven-comparison set; every signature there had k ≥ 6.
