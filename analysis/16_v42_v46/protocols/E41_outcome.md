# E41 outcome: inflamed expression signatures in published randomized checkpoint trials

Protocol E41_prereg.md (sha256 6227f3df…, logged 2026-11-01T12:46:04Z, before any search). Fourteenth protocol written
after the other results were known. Not yet public: the protocol was written after the tag v0.3.0 release.

## Search and extraction

Two blind AI-assistant passes (E41_pass1/, E41_pass2/; brief E41_brief.md), each running the pre-registered PubMed
query (629 records), the 41 seed trials and reference lists. Adjudication in E41_adjudicate.py (decisions in
results/E41_decisions.csv).

- Candidates rated eligible or unclear by either pass: 28. Eligibility verdicts agreed for 17. All 11 disagreements were
  comparisons whose full text pass 1 could not open (rated unclear); there was no direct yes/no conflict.
- Both passes extracted complete numbers for 11 comparisons; the numbers were identical for 10. IMvigor010 differed by
  model: pass 1 took covariate-adjusted HRs, pass 2 univariable ones. The univariable HRs were used.
- Seven comparisons rest on pass 2 alone (CheckMate 459, CLEAR, IMbassador250, IMbrave150, IMpower133, ORIENT-11,
  POPLAR). The coordinating session checked the POPLAR, CLEAR, CheckMate 459 and IMbrave150 values against the
  downloaded figures. IMpower133 was read from point and whisker positions (no printed values; about ±0.02).
- Eligible and estimable, independent of the six patient-level comparisons: 18, all with a survival endpoint (OS 11,
  PFS 7), in 10 cancer types. None had a response endpoint. A published subgroup result was found for one of the six
  (IMmotion150), and it was used only in the sensitivity analysis that adds it.
- Not eligible:
  - no 95% CI: AtezoTRIBE;
  - low-group HR not in the accessible text: IMmotion151, IMpower150;
  - within-arm associations or interaction P only: CheckMate 214, CASPIAN, JAVELIN Bladder 100, CheckMate 274,
    CheckMate 743, KEYNOTE-045;
  - no signature-stratified HRs found: RATIONALE-303, spartalizumab NPC phase 2, the pooled IMvigor urothelial
    analysis, GeparNuevo.

## Results (results/E41_meta.csv, E41_inputs.csv)

| analysis | k | estimate (95% CI) | P | I² |
|---|---|---|---|---|
| primary: −log(HR high / HR low), survival | 18 | +0.247 (0.088 to 0.406); ratio of HRs 0.78 (0.67 to 0.92) | 0.0044 | 7% |
| per SD (approximate) | 18 | +0.146 (0.057 to 0.235) | 0.0031 | 0% |
| 95% prediction interval, primary | | +0.016 to +0.477 | | |
| leave one out | 17 | range +0.202 to +0.275 | all ≤ 0.011 | |
| phase 3 only | 15 | +0.224 (0.060 to 0.388) | 0.011 | 0% |
| adding IMmotion150's published result | 19 | +0.273 (0.107 to 0.440) | 0.0029 | 17% |
| excluding all numbers read from figures | 5 | +0.260 (−0.155 to 0.676) | 0.16 | 44% |
| excluding the figure-position estimate only (post hoc) | 17 | +0.251 (0.079 to 0.424) | 0.0071 | 15% |
| breast | 3 | +0.583 (−0.372 to 1.539) | 0.12 | |
| non-small-cell lung | 3 | +0.596 (−0.546 to 1.738) | 0.15 | |
| signature, comparisons that also report PD-L1 | 10 | +0.238 (0.058 to 0.418) | 0.015 | 0% |
| PD-L1 immunohistochemistry, same 10 | 10 | +0.164 (−0.096 to 0.425) | 0.19 | 7% |
| tumor mutational burden | 7 | +0.071 (−0.304 to 0.447) | 0.66 | 0% |
| Egger test | 18 | z = 1.16 | 0.26 | |
| post hoc: without the three doubtful signature classes (OAK, ORIENT-11, RATIONALE-309) | 15 | +0.234 (0.044 to 0.424) | 0.019 | 21% |
| post hoc: without the covariate-adjusted IMbrave150 | 17 | +0.236 (0.078 to 0.394) | 0.0061 | 6% |

- Within the 10 comparisons that report both, the signature's value exceeded PD-L1's in 7 (no test; their correlation
  is unknown).
- A formal interaction test was reported in 1 of the 18 (CheckMate 649, likelihood-ratio P = 0.045, continuous score).
- The pooled per-SD estimate (+0.15) is close to the six-comparison patient-level estimates for the same family:
  Inflamed/IFN median +0.19, T-cell-inflamed expression profile +0.20, CD8A +0.23.

## Expectations (results/E41_expectations.csv): 5 of 5 held

| | expectation | observed | held |
|---|---|---|---|
| M1 | at least 10 eligible survival comparisons | 18 | yes |
| M2 | pooled y > 0, CI excluding 0 | +0.247 (0.088 to 0.406) | yes |
| M3 | per SD between +0.10 and +0.30 | +0.146 | yes |
| M4 | fewer than half report an interaction test | 1 of 18 | yes |
| M5 | signature and PD-L1 within 0.20 | +0.238 vs +0.164 | yes |

Tally after fourteen later protocols: 57 expectations, 44 held.

## Deviations and clarifications (decided before the analysis was run, after extraction)

1. **Per-SD conversion for top-versus-bottom splits.** The protocol's factor φ(z_q)/(q(1−q)) covers complementary
   splits. For top versus bottom tertiles, and ORIENT-11's high versus low clusters, the medium group was set aside and
   the factor generalized to φ(z_hi)/p_hi + φ(z_lo)/p_lo. This equals the protocol's factor for complementary splits.
2. **Adjusted versus unadjusted HRs.** Univariable or stratified HRs were preferred where both were given (IMvigor010).
   IMbrave150 published only covariate-adjusted HRs; a post hoc sensitivity analysis without it is reported.
3. **Signature class.** Three classes were judged eligible with doubt recorded: the OAK AP/T score (XCR1 and the
   IFN-γ-induced CXCL9, CXCL10, CXCL11) and the ORIENT-11 and RATIONALE-309 immune clusters. A post hoc sensitivity
   analysis without them is reported.
4. **Extra searches.** Pass 2 ran two extra PubMed searches with alternative wording, which are not in the protocol.
   They added one included comparison (ALICE), which pass 1 found through the protocol searches.
5. **Shared scratchpad.** Both passes wrote helper scripts to the shared session scratchpad, and pass 1 reported that
   its scripts were overwritten mid-run. It then moved its work to a private folder and states it did not open the
   other pass's files. Independence of the two passes therefore rests on their instructions and reports.
6. **Unverified extraction.** The extraction is AI-assisted. The authors must check results/E41_adjudicated.csv against
   the papers before submission.

## Outside the protocol (observation, not analyzed)

POPLAR's Fig. 4A also splits four single genes at the median. CD274, PDCD1LG2, PDCD1 and CD80 mRNA each separated
benefit as the Teff/IFN-γ signature did (HR 0.39–0.46 above the median, 0.97–1.13 below; signature 0.43 and 1.10).
Checked against the downloaded figure.

## Changes after the independent fact-check (v44)

- The cover-letter tally now reads fourteen protocols, 44 of 57.
- Code availability and the protocols paragraph now say that the fourteenth protocol is not part of the v0.3.0 release
  and is provided with the submission.
- The passes are described as blind to each other only, because the protocol they followed stated the patient-level
  estimates and the expectations.
- IMmotion151, IMpower150 and AtezoTRIBE are now described as eligible but not extractable from accessible material;
  the exclusion list reconciles to the 28 candidates.
- Covariate-adjusted HRs are named in the SI (IMbrave150, CLEAR, KEYNOTE-061, KEYNOTE-119) and mentioned in the
  limitations. Selective reporting is illustrated by CheckMate 649.
- "Identical for 10" now refers to the signature HRs only.
- The PD-L1, signature and TMB pools are given as ratios of HRs in the Results; the paired count (7 of 10) and the
  tumor-type pools are in the SI.
- IMpower133's position read is mentioned in the main text, and the 19 papers have full citations in S28a.
- The Results section moved after the sensitivity analyses so that figures are cited in order.
- The seed list is described as "listed in advance".
- Discussion and conclusion no longer imply specificity or a single-gene test from the aggregate data.
