# E33 outcome (protocol E33_prereg.md, sha256 611f13a8…, logged 2026-10-24T17:57:57Z)

Run: E33_trials.py (log results/E33_trials.log), E33_extras.py (breast pools of the barriernet and comparator scores,
per-comparison CD8A for the main-text table), Q7 recount from results/E33_pooled_reference_estimates.npz with the
direction audit (E33_direction_audit.csv / .md, run in parallel by an agent that saw no result). Step-1 reproduction:
library and CD8A per-comparison interactions equal E32 (max |diff| 2.2 × 10⁻¹⁶); paired differences from CD8A equal E32
(9.6 × 10⁻¹⁷); E32 random-reference P reproduced from the saved draws (1.1 × 10⁻¹⁶). 26 of 12,955 single-gene fits did
not converge (treated as missing).

| Q | Expectation | Result | Held |
|---|---|---|---|
| Q1 (primary) | CD8A not first among the 13 T-cell genes; ≥ 7 of 13 at P < 0.05 favoring the checkpoint arm | CD8A second (+0.226; PRF1 +0.229); 7 of 13 at P < 0.05 | yes |
| Q2 | CD8A and the T-cell median above the 90th percentile of the 2,082 genes measured in all six | CD8A 99.1%; T-cell median 94.2% | yes |
| Q3 | For ≥ 7 of 13 T-cell genes, at most 2 of 52 signatures exceed the gene (paired, P < 0.05) | 11 of 13 (none for 8 genes; CD3G 3, CD247 8) | yes |
| Q4 | descriptive | largest upper 95% limit of signature − CD8A +0.30 (IRG_Yang); median upper limit +0.05; median paired SE 0.076; difference detectable with 80% power 0.21 (CD8A interaction +0.23) | — |
| Q5 | ≤ 2 of 52 exceed CD8A over the five comparisons without IMmotion150 | 0 of 52 (11 below) | yes |
| Q6 | < 10 of 52 at P < 0.05 favoring the checkpoint arm in the three breast comparisons | 1 of 52 favoring (IFN_Ayers +0.47, P = 0.029), 2 against; 6 of 52 above CD8A there | yes |
| Q7 | Q2 count changes by ≤ 3 with audited directions | 26 → 28 (CD39-CD8Tcell_Chow and CYT_Davoli to sensitivity, IMS_Cui to resistance) | yes |

Six of six expectations held (Q4 descriptive). Across the six post hoc protocols, 18 of 25 expectations held.

Other results reported:
- IFN-γ/checkpoint genes: CXCL13 +0.28 (0.02 to 0.54), STAT1 +0.28 (0.14 to 0.43), LAG3 +0.26 highest of the 27; CXCL9
  +0.17 (−0.11 to 0.45; five comparisons, not measured in I-SPY2 pembrolizumab). Genes exceeded by > 2 signatures:
  CD3G 3, CD247 8, CXCL10 8, HAVCR2 17 (23 of 27 panel genes exceeded by at most 2).
- All 2,082 genes: 101 at P < 0.05 favoring the checkpoint arm and 102 against (9.8%); median
  −0.006; 95th percentile +0.165, 99th +0.224.
- CD8A in sensitivity sets: without IMmotion150 +0.16 (0.007 to 0.32; P = 0.040); CheckMate 025 on overall survival +0.19
  (−0.04 to 0.41; P = 0.088); renal adjusted for angiogenesis +0.27 (P = 0.029); I-SPY2 hormone receptor only +0.22
  (P = 0.035); breast only +0.06 (−0.33 to 0.45; P = 0.76).
- Breast only: signatures above CD8A at P < 0.05 — Bcell_Budczies, CTL_Jiang, IFN_Ayers, Teff-IFNG_Fehrenbacher,
  Teff_McDermott, VIGex_Hernando-Calvo (about 1.3 expected by chance in that direction); single genes at P < 0.05 there:
  GZMB +0.49, IDO1 +0.47, CXCL10 +0.45. barriernet: NetImmune +0.71 (0.19 to 1.23; P = 0.027), CAF_ECM −0.71 (−1.25 to
  −0.18; P = 0.029), Barrier −0.66 (P = 0.014), Immune +0.35 (P = 0.11); ESTIMATE +0.52 (P = 0.052).

Direction audit (52 source publications): 37 agree, 3 disagree (above), 12 unclear (8 report no checkpoint direction;
proliferation_Pabla non-monotonic; PTEN-MITF_Cabrita, IRG_Yang and ImmuneCells_Xiong cannot reproduce the published
score from the gene file). Quotes were taken through a web-fetch tool; the three disagreements were checked against
the abstracts by hand.

Rules applied: the main text says single T-cell genes in general were not outperformed, with the Q3 count, and that
CD8A ranked second of 13; reports the breast result in the main text; keeps the E32 SignatureSets labels as primary
with the audit as sensitivity; puts the detectable difference and the upper limits in the abstract.

Correction (2026-11-01): CD8A with CheckMate 025 on overall survival is +0.18 (file value 0.1845; the manuscript prints +0.18), not +0.19 as above.
