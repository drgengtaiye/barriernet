# E29 outcome against the pre-registration (E29_prereg.md, sha256 unchanged)

Recorded 2026-10-21. Script E29_dissociation.R (log results/E29_log.txt). Step 0 reproduced every published
marginal estimate of Immune and CAF_ECM (21 + 2 response cohorts, 10 survival cohorts; largest difference
4.9e-15). GSE215868 has no CAF_ECM score (< 3 genes on its panel; as published, where the CAF_ECM pool has
k = 20), so the response contrasts use 20 of the 21 cohorts — a property of the data, not a choice.

| Question | Pre-registered expectation | Result | Met? / rule |
|---|---|---|---|
| Q1 response, D_r (a), primary | > 0 (about +0.35) | +0.26 (0.08 to 0.45; HK), P = 0.0083, 15/20 positive; REML +0.26, P = 0.0035 | Yes → the manuscript may state, citing the direct test, that in single-arm checkpoint cohorts immune programs tracked response more strongly than fibroblast content |
| Q2 survival, D_s (a), primary | < 0 (about −0.09; possibly underpowered) | −0.08 (−0.22 to 0.06; HK), P = 0.22, 3/10 positive; REML −0.08, P = 0.17 | Not significant → "the strongest prognostic component" must be worded as the largest point estimate, not as a demonstrated difference from the immune modules |
| Q3 crossover (descriptive) | — | D_r − D_s = +0.35 (0.11 to 0.58), P = 0.004 (different scales; no decision) | Title rule: Q1 holds and D_s < 0 in point estimate → "dissociation" may be retained |
| Q4 collinearity | — | Spearman ρ(Immune, CAF_ECM): response median 0.20 (−0.40 to 0.59); survival 0.01 (−0.10 to 0.24) | No data set flagged |

**Sensitivity (no decisions).**
- S1 marginal (b): response +0.26 (−0.05 to 0.56), P = 0.097, k = 19 (one cohort's bootstrap covariance was not
  estimable: separation in small cohorts gave unstable resampled logistic coefficients, so (b) is wide);
  survival −0.10 (−0.24 to 0.05), P = 0.18.
- S2 REML: as in the table.
- S3 response: E5 11 +0.31 (P = 0.066), new 10 (k = 9) +0.31 (P = 0.077), whole transcriptome +0.25 (P = 0.029).
  Same direction in each subset; neither half significant alone.
- S4 gastric checkpoint cohorts (not pooled): PRJEB25780 D_r +0.41 (SE 0.68); PRJEB60680 +0.71 (SE 0.45).
  In PRJEB60680 the adjusted CAF_ECM coefficient is positive (+0.45), so its positive D reflects CAF_ECM
  pointing the "wrong" way, not a large immune effect.
- S5 survival: benchmark 3 −0.07 (P = 0.41), locked 7 −0.11 (P = 0.31); **adjusted for age and stage +0.02
  (P = 0.66)** — after clinical adjustment fibroblast content no longer leads.
- S6 single immune modules against CAF_ECM, response: APM +0.21 (P = 0.024), IFN +0.31 (P = 0.0064), TLS_B +0.25
  (P = 0.053), T_CELL +0.23 (P = 0.013); survival: −0.10 to −0.12, only TLS_B significant (P = 0.022).

## Post hoc (E29b_posthoc_adjusted.R; not pre-registered, added after the run)

Pooled mutually adjusted coefficients (same pooling rule), to show what drives the contrasts:

| Data | Immune | CAF_ECM |
|---|---|---|
| Response (20 cohorts), log OR per SD | +0.46 (0.25 to 0.67), P = 0.0002; positive in 19/20 | **−0.17 (−0.32 to −0.01), P = 0.037**; negative in 14/20 |
| Survival (10 gastric), log HR per SD | −0.19 (−0.28 to −0.11), P = 0.0007 (HR 0.82); protective in 10/10 | +0.28 (0.18 to 0.38), P = 0.0001 (HR 1.32); adverse in 9/10 |

Two consequences for the manuscript's wording, both post hoc and to be labeled as such if used:
1. Immune programs are independently prognostic in every gastric cohort, of similar size to fibroblast
   content; the survival information is shared between the two components, not carried by fibroblasts alone.
2. "Fibroblast content does not track response" holds marginally (about −0.06) but not conditionally: at a
   given level of immune activity, more fibroblast content went with lower response odds (−0.17 per SD,
   P = 0.037; compare E28's randomized interaction for CAF_ECM, −0.17, −0.43 to 0.08). The accurate statement is
   that response information is concentrated in the immune programs, with at most a modest adverse contribution
   from fibroblast content.

## Summary for the manuscript (not edited by E29)

- Keep, with the direct test: immune programs track checkpoint response more strongly than fibroblast content
  (single-arm cohorts, D_r +0.26, P = 0.008). E28 found no such difference for randomized treatment interaction
  (P = 0.47); the two statements must not be merged.
- Change: "fibroblast content is the strongest prognostic component" → largest point estimate; the immune
  modules are also independently prognostic (post hoc) and the difference is not significant (P = 0.22; +0.02
  after age and stage).
- Change: "fibroblast content tracks survival but not response" → "...and contributes little to response"
  or equivalent; do not say it carries no response information.

## Run notes

Two coding errors were fixed between runs; neither changed an analysis choice: the first run stopped in step 0
on GSE215868's missing CAF_ECM score (now skipped and logged, as published), the second stopped after step 0
on a vector-naming error before any E29 estimate had been written. The third run produced all results above.

Correction (2026-10-21, after the outcome was first written): the Q2 P value was first transcribed as 0.23;
the stored value is 0.2245 (results/E29_pooled.csv), i.e. P = 0.22. No decision depends on the change.
Second correction (same day): S6 survival range is −0.10 to −0.12 (T_CELL −0.1017), first written as −0.11 to −0.12;
the Q1 P value is 0.0083 (0.008286), reported to two significant digits in the manuscript.
