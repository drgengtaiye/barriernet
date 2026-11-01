# E31 protocol: the E30 duplicate rule applied to every pair of the ten gastric survival cohorts

Written 2026-10-22, after all other results were known and before the rule below is run on any pair other than
GSE26253–GSE62254. Prompted by a third external pre-submission review (point 4b: GSE15459 and GSE34942 both come
from Singapore on the same platform, and the duplicate screens reported so far did not use one rule for every pair).
Like E28–E30, E31 is *protocol-specified*: a post hoc question whose analysis is fixed here before it is run.
No outcome model is fitted unless the rule finds a duplicate (Action, below).

## Already known before this protocol

- **GSE15459 vs GSE34942.** No shared patient identifier; the E15 best-match check (2,000 genes) gave a maximum
  cross-cohort r of 0.67, in the range of pairs from different countries (GSE15459 vs GSE62254, 0.66; Supplementary
  Table S7b). A read-only check during the second review (session scratchpad `rev/agentB`, no outcome model) applied
  the E20 99.9th-percentile rule to every pair of the ten cohorts: 8 GSE15459 arrays exceeded GSE15459's
  within-cohort threshold (0.574) against GSE34942 arrays, each best-matching one of a few GSE34942 arrays that are
  the best match of five to eight GSE15459 arrays; none agreed with its match on age, sex and survival time; no
  GSE34942 array exceeded its own threshold (0.704) against GSE15459. Matching on clinical data alone (age, sex,
  overall survival) gave 1 of 95 GSE15459 deaths with a match in GSE34942, against 3.3 expected by chance.
- **Mutual best matches with sex and age agreement** (rules i–ii below): 11 between GSE26253 and GSE62254 (eight
  of them are the E30 duplicates) and at most one between any other two of the ten cohorts. Whether those single
  pairs satisfy rule iii has not been examined.

## Rule (the E30 rule, written for any pair)

For cohorts A and B, an array of A and an array of B are taken to be the same patient if (i) they are each other's
best cross-cohort expression match (Pearson correlation over all shared genes with no missing value and non-zero
variance in both, genes z-scored within each cohort: `dup_pair()` of E20, verbatim), (ii) sex agrees and age
differs by at most one year, and (iii) their overall survival does not conflict: if both died, the two times differ by
at most 0.5 months; if one died at time t and the other is censored at time c, c ≤ t + 0.5 months; two censored
records never conflict. For GSE26253, age, sex, overall survival and death come from Oh et al. 2018 Supplementary
Data 2 (sheet SMC), linked as in E30; for GSE62254 from the same file (sheet ACRG); for the other cohorts from the
sources used in the manuscript. A pair with a missing value in (ii) or (iii) does not meet the rule. GSE26253–GSE62254
is re-run as a check and must return the eight E30 arrays (rule iii here uses overall survival only, so it can
only add pairs; any difference is reported, and E30 remains the rule for that pair).

All 45 unordered pairs of TCGA-STAD, GSE84437, GSE62254, GSE15459, GSE34942, GSE26253, GSE13861, GSE26901,
GSE26899 and GSE29272 are screened. The E20 99.9th-percentile rule is also reported for every pair (descriptive).

## Expectation

No pair other than GSE26253–GSE62254 contains an array pair that meets rules i–iii.

## Action

If any pair does, the array is removed from the cohort analyzed later (locked over benchmark; extension over the
first three locked cohorts), and every pooled estimate in the main text that uses that cohort is re-run with the E30
machinery (`E30_dedup_*` overlay) and reported beside the current value. If none does, nothing is re-run.

## Outputs

`analysis/E31_dedup_all_pairs.R` → `results/E31_pairs.csv` (one row per ordered pair: genes, thresholds, E20 flags,
mutual best matches, matches meeting ii, matches meeting i–iii), `results/E31_candidates.csv` (every mutual best match
meeting rule ii, with the rule iii columns), `results/E31.log`. Reported in Supplementary Table S19h and one sentence
of the Methods.
