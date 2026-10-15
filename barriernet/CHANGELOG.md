# Changelog

## 0.2.1 (2026-09-23)

### Fixed
- **Visium neighbourhood.** `spatial.build_hex_graph` used axial hexagonal offsets
  (±1,0), (0,±1), (1,−1), (−1,1). Visium `array_row`/`array_col` are doubled-column
  coordinates (only row + col even is occupied; neighbours at (0,±2) and (±1,±1)), so on
  Visium data only the two diagonal steps ever matched: every spot had at most two
  neighbours (mean 1.95 on a real section, against 5.84 on the correct lattice) and
  TLS-centred rings grew along one diagonal. The Visium offsets are now the default;
  axial offsets remain available with `coords="axial"` for true axial coordinates. The
  graph builder warns when parity or mean degree show the wrong coordinate system.
- **Tissue edge.** Edge = fewer than six in-tissue neighbours. With the broken graph every
  spot was "edge", and a silent fallback (`far.sum() < 30`) then dropped edge exclusion
  from Far altogether. Edge exclusion is now explicit (`far_exclude_edge`, `edge_buffer`)
  and there is no silent fallback (a warning is raised when Far is small).
- **Far default.** `far_min_dist` default is now 3 hops (the definition used in the
  analyses); 0.2.0 defaulted to 4 in the CLI and API. Spots in tissue islands without a TLS
  peak (infinite distance) count as Far unless `far_include_unreachable=False`.
- **Packaging.** The built wheel omitted `data/*.yaml`, so an installed package could not
  load its default gene sets (only running from the source tree worked). Package data is
  now declared in `pyproject.toml`.
- **Orientation.** `load_expr_matrix(orientation="auto")` now decides samples×genes versus
  genes×samples by overlap with the marker genes; 0.2.0 used shape alone and misread, for
  example, 200 samples × 48 genes as genes × samples.

### Added
- `ring_summary` reports Ring1−Far as well as Ring2−Far contrasts; `barrier_index()`
  computes the section-level BarrierIndex = z(CAF_ECM Ring1−Far) + z(SPP1 TAM Ring1−Far)
  across sections (ddof 0); `barriernet spatial` writes `barriernet_sections.tsv` with it.
- `peak_method="lattice"` (maxima among the six lattice neighbours) beside the default
  `"knn"` (maxima among the `peak_k` nearest spots, the spot itself included); `peaks=` to
  pass a precomputed peak mask.
- `barriernet bulk --log2p1`, and a warning when input looks un-logged (max > 100):
  z-scoring linear TPM lets a few samples dominate each gene.
- Warning when an axis is built from fewer modules than defined (a module with < 3 genes).
- Tests (`tests/`): synthetic Visium lattice, real-section degree, exact regression against
  the corrected ring labels, gene-set and axis definitions, orientation.

### Effect on results produced with 0.2.0
The v13 manuscript's TLS-centred ring statistics and BarrierIndex were computed on the
broken graph, and so were the stored TLS peaks themselves: they are local maxima against at
most two diagonal neighbours (177 peaks in one section, against 142 with the documented
9-nearest-spot rule). Recomputed on the correct lattice, the direction of every ring result
is unchanged — SPP1⁺ TAM programmes lower next to TLS, no consistent CAF_ECM gradient,
IL-10/STAT3 higher — but the numbers change (for example SPP1_TAM Ring1−Far, stored peaks
P = 0.014, 8/10 sections; documented kNN peaks P = 0.0098, 9/10), BarrierIndex values change
(range −3.4 to +2.9 in v13; −2.6 to +3.5 with kNN peaks), and v13's edge-based robustness
checks were not valid. Bulk, proteomic and checkpoint scoring (`core.py`) is unaffected.

## 0.2.0
Initial public version.
