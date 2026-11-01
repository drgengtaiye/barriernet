# E39 protocol: how many patients a signature-versus-gene comparison needs (design analysis)

Written 2026-10-31, before E39_design.py is written or run. The hash of this file is logged in E39_prereg.sha256 first.

## Why

The manuscript reports that differences of less than about 0.21 per SD between a signature's interaction and CD8A's
could not be reliably detected in the six randomized comparisons. A reader needs the constructive counterpart: how much
information, and how many events, a trial or a set of trials needs to detect a given advantage over a single gene, and
whether a simple normal-theory formula predicts the uncertainty actually observed, so that it can be used to plan the
controlled-access analysis (E37_registration.md) and future trials.

## Data and estimates

E34 setup (E32 data, six primary comparisons), executed unchanged as in E37_test.py. The 52 library signatures,
NetImmune 0.2.1 and barriernet 0.3 (frozen weights, hash checked) against CD8A. Per comparison i and score: the paired
difference of interactions and its within-arm bootstrap variance v_i (300 resamples, seeds [3032, comparison], as
E32/E33); the Pearson correlation ρ_i of the score with CD8A across the comparison's patients (within-arm, pooled);
events D_i and the checkpoint-arm fraction p_i (survival comparisons); arm sizes n_a and response rates π_a (pCR
comparisons).

Step 0: the pooled paired differences and standard errors of the 52 signatures over the six must reproduce
E33_paired.csv (six, primary) to < 1e-8.

## Normal theory

Survival (Cox, 1 : p allocation, standardized score): var(difference of interactions) ≈ 2(1 − ρ) / (D p (1 − p)).
pCR (logistic): ≈ 2(1 − ρ) × Σ_a 1 / (n_a π_a (1 − π_a)). Information I = Σ_i 1 / v_i; a difference δ is detected with
80% power at two-sided 0.05 when I ≥ (1.96 + 0.84)² / δ².

## Outputs

1. Ratio of bootstrap to theoretical variance for every score × comparison; median and interquartile range, by comparison.
2. Current information and minimum detectable difference per score (reproducing the median 0.21).
3. Events required in a 1 : 1 survival comparison for δ = 0.05, 0.10, 0.15 and 0.20 per SD, as a function of ρ, and at
   the median ρ of the 52 signatures and at barriernet 0.3's ρ; the same in multiples of the information in the six
   current comparisons.
4. Added to the controlled-access analysis plan as a dated, hash-logged addendum (E37_registration_addendum.md): the
   difference from CD8A detectable with the events of the datasets actually obtained, computed when they are obtained.

## Expectations

- X1: the median ratio of bootstrap to theoretical variance over all score × comparison pairs lies between 0.67 and 1.5.
- X2: within each comparison, the Spearman correlation between bootstrap and theoretical variance across scores is at
  least 0.8 (median over the six comparisons).

Outputs 2 and 3 are descriptive. Reported whatever the result.
