# E26 outcome against the pre-registration (E26_prereg.md, sha256 unchanged)

Recorded 2026-10-19. GSE270678: 72 sections, 20 donors; 69 sections passed the pre-specified QC (three had
fewer than 300 spots with >= 500 UMIs), 104,716 spots.

| Test | Expectation | Result | Held? |
|---|---|---|---|
| P1 SPP1_arm Ring1 − Far | < 0 in more donors than not, one-sided P < 0.05 | 9 of 20 donors negative, P = 0.75 | **No** |
| P2 CAF_ECM Ring1 − Far | not consistently > 0 (one-sided P >= 0.05) | 14 of 20 positive, P = 0.058 | Yes (borderline) |
| P3 BarrierIndex vs stage | Spearman P >= 0.05 (sections) | rho −0.31, P = 0.010 | **No** |
| S1 fibroblast number → tumor T cells | coefficient < 0, P < 0.05 | +0.052 (0.013 to 0.091), P = 0.008 | **No** (opposite direction) |
| S2 continuity excess | not significantly negative | −0.016 (−0.046 to 0.014), P = 0.30 | Yes |

## Post hoc (E26b_bridge_posthoc.py; not pre-registered)

- The published ring result reproduces with the published inputs (Tangram states, peaks on TLS_niche_score):
  SPP1+ TAM state Ring1 < Far in 9 of 10 sections.
- With E26's inputs on the same 10 sections it weakens to 7 of 10 (P = 0.34). Swapping only the macrophage
  score (module SPP1_arm on the published peaks) keeps 9 of 10, so the difference lies in the TLS peaks.
- The published TLS_niche_score (Tangram B cell + plasma cell + DC fractions) is unrelated to lymphoid markers
  at spot level (median r: MS4A1 0.09, CD79A −0.02, CXCL13 0.01, CR2 0.08, CCL19 0.00, CCL21 0.01), and its
  plasma-cell component runs against MZB1 (median r −0.16). The marker TLS_B module tracks them (MS4A1 0.55,
  CD79A 0.57, CXCL13 0.50) including markers outside the module (CR2 0.27, CCL19 0.32, CCL21 0.25).
  The published TLS-centered rings were therefore not centered on lymphoid aggregates.
- BarrierIndex vs stage at donor level (stage is a donor property): rho −0.44, P = 0.050.

**Consequence.** Under the reporting rule set in advance, the published Visium reading must be qualified. The
published "SPP1+ macrophage states were lower next to TLS (9 of 10; P = 0.0098)" and the IL-10/STAT3 contrast
rest on peaks that do not mark lymphoid aggregates and should be recomputed with a marker-based TLS definition.
