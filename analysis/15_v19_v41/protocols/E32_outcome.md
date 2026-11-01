# E32 outcome (protocol: E32_prereg.md, sha256 a33c1fd7…, logged 2026-10-22T18:01:38Z)

Runs (2026-10-22, all after the hash): `E32_trials.py` (`results/E32.log`), `E32_locked_bionull.R`
(`results/E32_locked_bionull.log`), `E32_native.py` + `E32_native_scores.R` (`results/E32_native.log`).
Step 0 re-derived all 72 published interaction / within-arm values (largest difference 2.2 × 10⁻¹⁶); the Q7 observed
per-cohort NetImmune log HRs equal the current values (4.4 × 10⁻¹⁶).

Implementation notes (none changes the protocol): the first `E32_trials.py` run stopped when one random-reference Cox fit
did not converge; failed reference fits are now recorded as missing (6 of 247,000) and that draw is pooled over its other
comparisons. The library loaded 52 signatures (Inflamed/IFN 20, Other 24, Confounder 8; 34 declared sensitive);
ImmunoScore_Hao had one missing weight, set to 1. Extra rows include the barriernet subset modules C1Q_arm and SPP1_arm
(present in the shared module list). Q6 is reported with z-scores recomputed within the subset (and with the published
z-scores; identical to the second decimal).

## Expectations: 6 of 8 held

| | Expectation | Result | Held |
|---|---|---|---|
| Q1a | ≥ 70% of Inflamed/IFN favor the checkpoint arm | 19/20 (95%) | yes |
| Q1b | 10–40% of the 52 reach P < 0.05, all favoring the checkpoint arm | 18/52 (35%); 17 favor, **IMPRES_Auslander against** (−0.23, P = 0.018) | **no** |
| **Q2 (primary)** | fewer than 20% of the 52 exceed their size-matched random reference | **26/52 (50%)**: Inflamed/IFN 15/20, Other 11/24, Confounder 0/8 | **no** |
| Q3 | no signature exceeds CD8A (at most 2 of 52) | 0/52 exceed; 15/52 significantly below CD8A | yes |
| Q4 | fewer than one-third of the Q1-significant keep P < 0.05 after ESTIMATE-immune adjustment | 2/18 (IRG_Yang, IMPRES) | yes, but see below |
| Q5 | NetImmune six-comparison pool keeps its sign with I-SPY2 adjusted for HR + MammaPrint | +0.32 (0.03 to 0.60; P = 0.037; 6/6) | yes |
| Q6 | NetImmune and Immune keep their sign in I-SPY2 pembrolizumab with Linear-adjusted controls only | NetImmune +1.47 (P = 0.015), Immune +1.18 (P = 0.046); both larger than with all controls | yes |
| Q7 | NetImmune locked-7 pooled ranks in the top 5% of immune/stromal random axes | rank 1 of 1,001 (P ≤ 0.001); per cohort P < 0.05 in 5 of 7 | yes |

Later protocols, cumulative: 12 of 19 expectations held (E28–E31 6 of 11; E32 6 of 8).

## What the results say

1. **Immune signatures carry real, shared predictive signal.** Of the 52 published signatures, 17 favored the checkpoint
   arm at P < 0.05 and 26 exceeded their size-matched random reference (15 of 20 Inflamed/IFN; no Confounder signature).
   The primary expectation (that few would exceed random genes) was wrong: against random genes the immune signatures
   are specific. Caveats: the random references were centered below zero (median of their medians −0.06), which makes
   "exceeding" easier; 11 of the 26 exceeded their reference without reaching nominal P < 0.05; and random gene sets
   themselves reached pooled nominal P < 0.05 in 0.4–21% of draws (median 4.8%; highest for large sets).
2. **No signature beats a single gene.** CD8A alone gave +0.23 (0.05 to 0.40; P = 0.022). No signature's interaction
   exceeded CD8A's (paired bootstrap); 15 were significantly smaller.
3. **The signal cannot be separated from general immune infiltration, but it is not explained away by it.** Adjusted for
   the ESTIMATE immune score and its interaction with arm, only 2 of the 18 kept P < 0.05, yet the adjusted point estimates
   were larger, not smaller (median ratio 1.75, range 0.93–2.53): the loss of significance comes from collinearity. Q4's
   wording ("explained by one immune axis") is therefore not supported even though its count criterion was met.
4. **The pooled evidence is fragile.** Without IMmotion150, 9 of 52 signatures reached P < 0.05 (all favoring) and
   NetImmune P = 0.095; with CheckMate 025 on overall survival, 1 of 52 (against the checkpoint arm) and NetImmune
   P = 0.058; adding PRINCE and DUTRENEO, 22 of 52 (20 favoring); I-SPY2 adjusted for hormone receptor only, 15 (14);
   renal trials adjusted for angiogenesis, 19 (18). Published engines (GSVA, ssGSEA, weighted mean): 24 of 52 at P < 0.05,
   IMPRES again against (`E32_native.csv`).
5. **barriernet.** NetImmune +0.32 (P = 0.037) with I-SPY2 adjusted for MammaPrint; the I-SPY2 pembrolizumab
   processing confound does not explain the breast interactions (they grow when only Linear-adjusted controls are used;
   six-comparison NetImmune with that subset +0.31, P = 0.057). In the seven locked gastric cohorts NetImmune ranked first
   of 1,001 immune/stromal random axes (review point 2), so its locked prognostic value is beyond immune/stromal random
   axes, not only beyond all-gene axes.

## Rules applied (as written in the protocol)

- Q2 headline: report the proportion (26/52, 50%) with its random-reference distribution; the "non-specific" wording
  planned in the protocol's reporting section does not apply, because Q2 did not come out as expected. The supported
  wording is: immune signatures beat random genes but not CD8A, and cannot be separated from general immune infiltration;
  the pooled evidence depends on IMmotion150 and on CheckMate 025's PFS endpoint.
- Q3: the text states that a single gene, CD8A, matched every published signature.
- Q4: the text states that signatures and general immune infiltration cannot be separated (collinearity), and does not
  say the signal is explained by it.
- Q5: the MammaPrint-adjusted NetImmune value (+0.32, 0.03 to 0.60; P = 0.037) replaces the hormone-receptor-only value.
- Q6: reported in the main text (the confound does not explain the breast interactions).
- Q7: the text says NetImmune's locked prognostic value exceeds immune/stromal random axes.
