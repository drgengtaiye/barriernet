# E25 outcome against the pre-registration (E25_prereg.md, sha256 unchanged)

Recorded 2026-10-19. Step 0 re-derived all 72 published interaction and within-arm values (largest
difference 2.2e-16); the numpy pooling rule matched E21_meta.csv to 6.4e-5. Null gene pool: 1,899 genes
measured and non-constant in all six comparisons.

| Question | Pre-registered expectation | Result | Met? |
|---|---|---|---|
| Q1 CAF_ECM | checkpoint arms P < 0.05; control P ≥ 0.05 | checkpoint −0.29 (−0.59 to 0.01), P = 0.058, 0/6 favorable; control −0.04, P = 0.41 | **No** (checkpoint P = 0.058) |
| Q1 NetImmune | checkpoint P < 0.05; control P ≥ 0.05 | checkpoint +0.47 (0.08 to 0.87), P = 0.028; control +0.05, P = 0.35 | Yes |
| Q2 six comparisons | both adjusted interactions P < 0.05 | Immune +0.22 (0.07 to 0.37), P = 0.013; CAF_ECM −0.22 (−0.48 to 0.04), P = 0.085 | **No** |
| Q2 breast | both P < 0.05 | CAF_ECM −0.62 (−1.03 to −0.20), P = 0.0056; Immune +0.28, P = 0.20 | **No** |
| Q3 NetImmune, six | empirical P < 0.05 | +0.29, empirical P = 0.006 (nominal 0.034) | Yes |
| Q3 CAF_ECM, breast | empirical P < 0.05 | −0.60, empirical P = 0.0015 | Yes |
| Q3 renal | neither below 0.05 | NetImmune 0.067; CAF_ECM 0.39 | Yes |
| Q4 | nominal false-positive rate rises with n, from about 5% at 100 | NetImmune structure 10.2% (n 100) to 18.8% (726); CAF_ECM 9.9% to 19.2%, not monotonic | **No**: already about 10% at 100 patients, so the inflation is a property of JAVELIN Renal 101, not of trial size |

Not pre-specified but computed by the same run: CAF_ECM pooled interaction over six comparisons,
empirical P = 0.031 (nominal 0.14); within checkpoint arms, empirical P = 0.002.

**Decision rule.** The reframing required Q1 for CAF_ECM and either Q2 or the breast part of Q3. Q1 for
CAF_ECM was not met (P = 0.058), so the reframing is not adopted; the manuscript keeps its framing and E25 is
reported as a sensitivity analysis.
