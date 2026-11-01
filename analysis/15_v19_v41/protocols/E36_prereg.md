# E36 protocol: operating characteristics of the calibration framework (semi-synthetic simulation)

Written 2026-10-30, before E36_sim.py is written or run. The hash of this file is logged in E36_prereg.sha256 first.

## Why

The manuscript is reframed as a method: predictive claims for expression scores are tested against treatment
interaction, size-matched random genes and single genes. A method needs its operating characteristics. This simulation
keeps the real expression, arms and covariates of the six primary randomized comparisons and replaces the outcomes by
outcomes simulated under a known truth, so that error rates and power are measured at the trials' real sizes,
gene correlations and designs.

## Data and models

E34 setup (E32 data: six primary comparisons; I-SPY2 adjusted for hormone receptor and MammaPrint, NeoTRIP for batch).
Observed outcomes are used only for each comparison's event fraction (survival) or control-arm pCR rate.
Benefit-oriented linear predictor eta. Survival comparisons: event time exponential with rate exp(−eta); censoring
exponential with rate set so that the expected event fraction at eta = 0 equals the observed one. pCR comparisons:
logit P(pCR) = logit(observed control-arm rate) + eta. Analysis exactly as the manuscript: score standardized, model with
score, arm and their product (Cox; logistic with the comparison's covariates, Firth on separation), interaction pooled
over the six by the manuscript rule (wider of Hartung–Knapp and DerSimonian–Laird, larger P). For speed, the Cox model
is a Newton–Raphson partial-likelihood fit (no ties occur in simulated times), checked against lifelines to < 1e-6
before use.

Markers (scored as in the manuscript): CD8A; barriernet T_CELL, IFN, APM, TLS_B, CAF_ECM modules; NetImmune. P = the
CAF_ECM module (standardized); a = arm.

## Scenarios

- S0 (null): eta = 0. 1,000 replicates.
- S1 (prognostic only): eta = −0.26 P (hazard ratio 1.3 per SD of CAF_ECM, same in both arms). 1,000 replicates.
- S2 (predictive): eta = beta a T, T = CD8A, beta in {0.1, 0.2, 0.3}. 1,000 replicates each.
- S3 (random-gene reference): 200 replicates each of S0, S2 with beta = 0.3 (truth CD8A), and S3d (diffuse truth:
  T = first principal component of all complete genes of the comparison, sign set so that it correlates positively with
  CD8A, beta = 0.3). In each replicate every marker is compared with 200 size-matched random gene sets (pooled draw by
  draw), and "exceeds" means above the 95th percentile in the marker's declared direction (CAF_ECM −, others +).
- S4 (single-gene reference): 200 replicates each with truth T = CD8A, T_CELL, IFN, APM or TLS_B, beta = 0.3. In each
  replicate the paired difference of each module from CD8A (100 within-arm bootstrap resamples per comparison,
  inverse-variance pooled) is tested. The true difference of each truth is its mean estimated difference over the
  replicates (oracle).

Seeds: numpy default_rng([3036, scenario index, replicate]).

## Measures and expectations (written before running)

| | Measure | Expectation |
|---|---|---|
| E1 | S0: rejection rate of the pooled interaction test (P < 0.05), each marker | 2.5% to 7.5% |
| E2 | S1: CAF_ECM interaction rejection rate; CAF_ECM association in the checkpoint arms alone (single-arm analysis, pooled) | ≤ 7.5%; ≥ 50% |
| E3 | S2: power for CD8A at beta = 0.1, 0.2, 0.3 | increasing; ≥ 60% at 0.3 |
| E4 | S3: exceeding the random reference — S0, every marker; S2 (CD8A truth), T_CELL; S3d, every module | ≤ 10%; ≥ 50%; descriptive |
| E5 | S4, truth CD8A: a module significantly above CD8A (P < 0.05, positive difference), each module | ≤ 5% |
| E6 | S4, truth a module: power of the paired test against the oracle difference | ≥ 50% where the oracle difference ≥ 0.20; ≤ 30% where ≤ 0.10 |

## Reporting

All measures are reported in the main text (one figure) and in full in the Supplementary Information, whether or not
the expectations hold.
