# E32 protocol: calibrating signature × treatment interaction in randomized checkpoint trials

Written 2026-10-22 and confirmed by the authors; hash-logged before any analysis below that involves outcomes is run.
Prompted by a fourth external pre-submission
review, which judged the paper's real novelty to be the calibration of score × treatment interaction in randomized
trials (16% of structure-matched random gene sets reached P < 0.05 in JAVELIN Renal 101 while the nominal test was
calibrated; CD8A alone gives most of the composite's interaction), and by the authors' decision to rebuild the paper
around it. Like E28–E31, E32 is protocol-specified: post hoc questions whose analyses are fixed here before they are run.

## Already known before this protocol (declared; not re-tested as discoveries)

- **E24 (archived, never reported).** 107 SignatureSets lists (z-mean, ≥ 5 genes), score × arm in the same six
  comparisons, no random-gene reference and no protocol: 27 of 104 scorable lists had pooled P < 0.05, all favoring
  the checkpoint arm; 77 of 104 positive; of the 52 companion signatures, 40 scorable and 15 significant (Inflamed/IFN
  11/20, Other 4/15, Confounder 0/5). Signatures lay on one immune–stromal axis; TCGA prognostic and pooled predictive
  effects correlated (Spearman 0.63). The attempt to lead with "prognostic ≠ predictive" was abandoned because that
  claim did not hold. NetImmune +0.29 (P = 0.034), ESTIMATE +0.30 (0.009), TMEscore +0.30 (0.035), CD8A +0.22 (0.035).
- **Companion project, manuscript v5 (randomized analyses since removed from it).** 31-set trial library (mean-z);
  inflamed 12-gene score: three added trials one-sided P = 0.041, seven comparisons P = 0.0004; proliferation control
  favored the checkpoint arm in the added trials (P = 0.034; JAVELIN two-sided P = 0.015, the strongest interaction in
  that library). Against size-matched random sets: 8/17 immune signatures in IMmotion150, 5/17 in PRINCE, 0/17 in each
  of CheckMate 025, JAVELIN Renal 101 and NeoTRIP, 0/17 in both I-SPY2 comparisons; controls 0/3. In I-SPY2 with the
  52-signature library (logistic without hormone receptor): MHC-I_Liu P = 0.009, IMPRES P = 0.010 (opposite direction),
  NonResponse_Chen P = 0.040 (durvalumab).
- **barriernet (E10, E21, E25, E28, E30).** NetImmune six-comparison +0.29 (0.03–0.55; P = 0.034; 6/6), Immune and IFN
  significant, CAF_ECM 0/6; JAVELIN: 16% of 2,000 structure-matched random sets P < 0.05; the parametric bootstrap held
  nominal size (27/30 score × comparison pairs in 3–7%); common random reference: NetImmune empirical P = 0.006,
  CAF_ECM 0.031; Holm 0.068; CheckMate 025 overall survival: no pooled interaction significant.
- **Not yet known.** No one has scored the full companion library in all six comparisons against per-signature
  size-matched random gene sets; compared each signature with CD8A; adjusted signature interactions for general immune
  infiltration; adjusted I-SPY2 for MammaPrint; checked the I-SPY2 processing confound below; or ranked NetImmune in the
  locked gastric cohorts against immune/stromal random axes.
- **I-SPY2 pembrolizumab processing confound (found by a read-only inventory; counts only).** In GSE194040 the field
  "adjustment method" is aligned with arm: all 69 pembrolizumab patients are "Linear Adjustment Factor", controls 139
  "ComBat" and 40 "Linear". MammaPrint (0 = MP1, 1 = MP2): pembrolizumab 31/38, control 91/88; HR− MP1 has 2 patients
  in the pembrolizumab arm. Durvalumab–olaparib: MP1/MP2 24/47 and 9/20.

## Data, scores and model

- **Comparisons (primary):** IMmotion150 (2,475-gene panel), I-SPY2 pembrolizumab, I-SPY2 durvalumab–olaparib,
  CheckMate 025, JAVELIN Renal 101, NeoTRIP, loaded exactly as in E21/E25 (E28_base data code, imported only after
  hashing). Sensitivity: PRINCE, DUTRENEO.
- **Library:** the companion's 52 signatures (gene lists and weights from SignatureSets IO CSVs; family from
  v5_signature_clusters: Inflamed/IFN 20, Other 24, Confounder 8; declared direction: 34 sensitive, 18 resistance),
  plus ESTIMATE (immune − stromal), TMEscore (z-mean A − B), CD8A and the barriernet modules and axes (CAF_ECM,
  TAM_SUP, APM, IFN, TLS_B, T_CELL, Barrier, Immune, NetImmune). Gene aliases as in barriernet.
- **Primary scoring engine (common):** genes z-scored within each comparison across arms; score = weighted mean of the
  z-scores of the genes present (published weights and signs; weight 1 otherwise); scored in a comparison when ≥ 3
  genes and ≥ 50% of the list are present; a signature enters pooling when scored in ≥ 4 of the 6 comparisons, and in
  every pooled set it is reported with the comparisons it uses. Score standardized to SD 1 in each comparison.
  Sensitivity: the published engine (GSVA, ssGSEA or weighted mean) for the 52.
- **Model (as published):** Cox PFS ~ s + arm + s:arm (renal); logistic pCR ~ s + arm + s:arm + covariates (breast);
  Firth where maximum likelihood separates. **Changed for E32 (review point):** I-SPY2 covariates are hormone receptor
  and MammaPrint class (MP1/MP2) as main effects; NeoTRIP batch as published. Benefit-oriented interaction θ (−log HR,
  +log OR per SD), oriented additionally by declared direction (θ × −1 for resistance signatures) when counting
  "declared-direction" results.
- **Pooling:** the manuscript rule for randomized comparisons (REML estimate; wider of REML–Hartung–Knapp and
  DerSimonian–Laird intervals; larger P); also the two-sided sign test over comparisons.
- **Size-matched random reference (per signature):** for each signature and comparison, 1,000 random gene sets with the
  same number of genes present, drawn from that comparison's measured, non-constant genes (the IMmotion150 panel is
  immune-enriched, so its reference is harder), carrying the signature's weight vector in random order; scored and
  modeled exactly as the signature; draw b of all comparisons pooled together. Empirical P (declared direction) =
  (1 + #{θ_b ≥ θ_obs}) / 1,001. Seed 2032. Random genes carry real interaction signal (JAVELIN), so this is a reference
  distribution, not a null; its question is whether a signature's interaction is specific to its genes.

## Questions and expectations

| | Question | Expectation (written before running) | Rule applied |
|---|---|---|---|
| Q1 | How many signatures have a pooled interaction favoring the checkpoint arm, and how many reach nominal P < 0.05? | ≥ 70% of Inflamed/IFN signatures favor the checkpoint arm; 10–40% of all 52 reach P < 0.05, all of them in the benefit direction (E24 basis) | Reported as the extent of the shared signal |
| **Q2 (primary)** | How many signatures exceed the 95th percentile of their own size-matched random reference (empirical P < 0.05, declared direction)? | Fewer than 20% of the 52 | The headline: the proportion, whatever it is, with its random-reference distribution |
| Q3 | Does any signature's pooled interaction exceed CD8A alone? Paired difference per comparison by 300 patient-bootstrap resamples within arm (both models refitted per resample), pooled by inverse variance | No signature exceeds CD8A at P < 0.05 (at most 2 of 52 by chance) | The text says whether a single gene matches the published signatures |
| Q4 | Is the shared signal general immune infiltration? Each signature's interaction re-estimated with the ESTIMATE immune score and its interaction with arm in the same model | Fewer than one-third of the signatures significant in Q1 keep pooled P < 0.05 | The text states how much of the library's interaction is explained by one immune axis |
| Q5 | I-SPY2 with MammaPrint: does the barriernet six-comparison NetImmune interaction keep its sign with I-SPY2 adjusted for hormone receptor and MammaPrint? | Keeps its sign (CI excluding zero not required) | The adjusted value replaces the hormone-receptor-only value in the main text |
| Q6 | I-SPY2 processing confound: NetImmune, Immune, CAF_ECM, ESTIMATE and CD8A interactions in I-SPY2 pembrolizumab with only the 40 "Linear Adjustment Factor" controls | Same sign as with all controls for NetImmune and Immune | Reported in the main text whatever the result; if a sign flips, the breast evidence is described as confounded by processing |
| Q7 | Gastric locked cohorts (review point 2): NetImmune's pooled locked-7 absolute log HR ranked among 1,000 structure-matched axes drawn from the immune and stromal gene pool (as for the benchmark cohorts, restricted to genes measured in all seven) | Ranks in the top 5% (empirical P ≤ 0.05) | If not held, the text says its locked prognostic value is not beyond immune/stromal random axes |

Sensitivity (reported, no expectations): published engines; five comparisons without IMmotion150; + PRINCE and
DUTRENEO; I-SPY2 with hormone receptor only; CheckMate 025 overall survival; renal trials adjusted for the
angiogenesis score for Confounder/Other signatures (as E30).

## Reporting decided in advance (review point 3)

The six-comparison summaries are reported with the sign test as well as the random-effects pool, and described as
"direction consistent, significance unstable, not specific" if Q2–Q4 come out as expected.

## Outputs

`analysis/E32_*.py` / `.R`; `results/E32_*`: per-signature × comparison estimates, pooled table, random-reference
draws and empirical P, CD8A paired differences, ESTIMATE-adjusted interactions, I-SPY2 sensitivity, locked-cohort
immune/stromal null; `E32_outcome.md` with every expectation marked held / not held.
