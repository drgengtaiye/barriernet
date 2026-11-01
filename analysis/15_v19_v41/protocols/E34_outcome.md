# E34 outcome (protocol E34_prereg.md, sha256 327c8a1e…, logged 2026-10-25T16:36:59Z)

Data: NeoPembrOV (GSE227666 counts; outcomes from the Source Data of Derbel et al., Nat Commun 2024, PMC11252308),
53 pre-treatment samples (32 NACT+P, 21 NACT), PFS 45 events; 24-month status for 52. Patient-to-library links made
from expression only before the protocol (id_match_final.csv; 2 links flagged with margin < 0.10). All 52 signatures,
the composite scores and the 27 panel genes were scorable. Script E34_trials.py (log results/E34_trials.log); step 0
reproduced the published interactions (2.2 × 10⁻¹⁶).

| Q | Expectation | Result | Held |
|---|---|---|---|
| Q1 | seven-comparison counts within ±3 of 17 / 1 | 18 favoring / 1 against | yes |
| Q2 | ≤ 2 of 52 above CD8A over seven | 0 (14 below); largest upper limit +0.29 | yes |
| Q3 | < 10 favoring in the four chemotherapy-backbone comparisons | 1 favoring / 2 against; 6 above CD8A (CTL_Jiang, GeparSixto_Karn, IFN_Ayers, Teff-IFNG_Fehrenbacher, Teff_McDermott, VIGex) | yes |
| Q4 | descriptive | NeoPembrOV alone: no signature at P < 0.05; CD8A +0.06 (SE 0.30), NetImmune +0.14 (0.38), CAF_ECM −0.03 (0.31); T-cell genes +0.01 to +0.40, none significant | — |

Pools: seven comparisons — CD8A +0.22 (0.06 to 0.37; P = 0.013; 6/7), NetImmune +0.29 (P = 0.024), ESTIMATE +0.30
(P = 0.004), STAT1 +0.27 (P = 0.002). Four chemotherapy-backbone comparisons — CD8A +0.06 (−0.27 to 0.38; P = 0.72),
NetImmune +0.57 (P = 0.043), ESTIMATE +0.44 (P = 0.027), CXCL10 +0.41 (P = 0.021), CAF_ECM −0.52 (P = 0.084).
Sensitivity (CD8A in NeoPembrOV): 24-month status −0.04 (SE 0.61); dropping the 2 flagged links +0.25 (SE 0.31).
Three of three expectations held. Conclusion unchanged: adding the one directly downloadable randomized comparison (53
patients, about 3% more) moves no count by more than one; the breast pattern (CD8A near zero, six signatures above it)
persists in the chemotherapy-backbone pool.

Correction to the protocol text (written after logging, so the protocol file is unchanged): the source paper's first
author is Le Saux (Le Saux O, Ardin M, Berthet J, et al. Nat Commun 2024;15:5932; doi 10.1038/s41467-024-47000-5;
PMID 39013886), not "Derbel" as the protocol says.

Correction to the linking description (protocol text and the first v39 SI draft): the chemotherapy arm was linked by the
10 genes reported per patient for both arms, not by 21 features — the 11 CD8B ratios exist only for the pembrolizumab
arm. Pembrolizumab-arm second-best correlations reach 0.78 (smallest margin 0.17), not 0.59. analysis/E34_link.py now
reproduces both link files from the downloaded data (identical assignments) and shows that the 10 genes alone reproduce
all 32 pembrolizumab-arm links. The 2 flagged links are both in the chemotherapy arm.
Sensitivity detail: with 24-month progression one of 93 scores reached P < 0.05 (Tcell-exclusion_Arnon +1.77, P = 0.019);
without the two flagged links one did (TLS_B +0.74, P = 0.045); CD8A −0.04 (SE 0.61) and +0.25 (SE 0.31).
