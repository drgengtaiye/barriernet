# Pre-registration: E26, an independent gastric Visium atlas (GSE270678)

Written 2026-10-19 while the data were downloading, before any GSE270678 expression value, coordinate or
image was read. GSE270678 (Lee et al., spatially resolved gastric cancer atlas): 72 Visium sections from 20
donors, sampled at tumor core, edge and in between; stage and Lauren subtype per section.

Nothing is fitted to these data. The module genes are *barriernet* 0.2.1; rings, peaks and the BarrierIndex
are the package's `assign_rings` / `ring_summary` / `barrier_index` with their defaults (TLS peaks as local
maxima among the 9 nearest spots above the 95th percentile; Ring1 = 1 hop; Far >= 3 hops; Visium
doubled-column lattice).

## Processing (fixed)

- Spots: in-tissue spots of the filtered matrix with at least 500 UMIs; sections with fewer than 300 such
  spots are excluded.
- Expression: counts per spot normalized to 10,000 and log1p-transformed (natural log), per section.
- Spot module scores (within-section, for ring analyses): per-gene z-score across the section's spots, module
  = mean z of its genes present (at least 3), exactly as the bulk scoring. TLS score = TLS_B module.
- Absolute scores (for between-section analyses): module = mean log-normalized expression of its genes,
  without z-scoring. Epithelial score = mean of EPCAM, KRT8, KRT18, KRT19, CLDN4, CDH1 (not a barriernet module).
- Donor level: section values averaged within donor (20 donors) for sign tests.

## Primary: does the published Visium reading replicate in an independent gastric atlas?

The published analysis (10 sections) found SPP1+ macrophage states lower next to lymphoid aggregates than in
distant tissue (Ring1 − Far lower in 9 of 10 sections), no fibroblast or myeloid layer around them, and a
BarrierIndex unrelated to age, pT or pN.

- **P1.** SPP1_arm Ring1 − Far < 0 in more donors than not; one-sided sign test over donors, P < 0.05.
- **P2.** No fibroblast layer around lymphoid aggregates: CAF_ECM Ring1 − Far is not consistently positive
  (one-sided sign test for > 0 over donors, P >= 0.05).
- **P3.** BarrierIndex (package definition, across sections) not associated with stage (Spearman, sections,
  P >= 0.05).
- Reported without expectation: TAM_SUP, C1Q_arm, IFN and T_CELL Ring1 − Far; differences by sampling location
  and Lauren subtype.

## Secondary: fibroblast number against fibroblast arrangement in gastric tissue

- Tumor spots: epithelial score at or above the median over all spots of all sections. Fibroblast-high spots:
  absolute CAF_ECM at or above the 75th percentile over all spots of all sections.
- Per section: fibroblast number = fraction of in-tissue spots that are fibroblast-high; boundary = tumor spots
  with at least one non-tumor lattice neighbor; continuity = fraction of boundary spots with at least one
  fibroblast-high neighbor; expected continuity = mean over 200 random placements of the same number of
  fibroblast-high labels among the section's non-tumor spots (seed 2026); excess = observed − expected.
  T-cell infiltration = mean absolute T_CELL over tumor spots.
- Model: T-cell infiltration ~ z(fibroblast number) + z(continuity excess) + z(fraction tumor spots), donor
  random intercept (statsmodels MixedLM, REML), sections with at least 30 boundary spots.
- **S1.** Fibroblast number coefficient < 0 with P < 0.05.
- **S2.** Continuity excess coefficient not significantly negative (P >= 0.05 or > 0), i.e. arrangement adds
  nothing once number is known, as in the imaging mass cytometry data.

Limitations stated in advance: Visium spots (55 µm) hold several cells, so "arrangement" here is a coarse
proxy for the single-cell layer measured by imaging mass cytometry; 20 donors give limited power; replicate
sections of one donor are not independent (hence donor-level tests and random intercepts).

## Reporting

All results are reported whatever their direction. If P1–P3 hold, E26 is described as an independent
replication of the Visium findings; if not, the published Visium reading is qualified in the text.

## Deviations

(none yet)
