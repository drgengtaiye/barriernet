# Changelog

## 0.4.1 (2026-11-08)

### Added
- The web app also runs inside a web browser (Shinylive, Python compiled to WebAssembly) and is published on GitHub
  Pages at https://drgengtaiye.github.io/barriernet/; nothing a visitor uploads leaves their computer. A browser has no worker processes, so there the
  benchmark runs in the page on one core with lower limits (200 random sets, 100 resamples).
- A synthetic demonstration computed in advance (1,000 random sets, 300 resamples; `data/demo_benchmark/`,
  `demo.precomputed_benchmark`), shown at once in the app.
- `.github/workflows/pages.yml` builds the browser version with Shinylive and deploys it to GitHub Pages.

### Changed
- `multiprocessing` is optional; `shiny_app/make_bundle.sh` no longer needs rsync.

## 0.4.0 (2026-11-08)

### Added
- **Signature library** (`barriernet.library`, `data/signatures_io.tsv`): the 52 immuno-oncology entries of
  SignatureSets (BHKLAB, MIT License, `data/SignatureSets_LICENSE.txt`) with published weights (1 where none is given),
  declared direction and family (Inflamed/IFN 20, Other 24, Confounder 8); TIS_Damotte repeats the 18 genes of
  Tcell-inflamed_Ayers, so 51 are distinct. The 13-gene T-cell panel, the 14-gene interferon-γ and checkpoint panel,
  and two T-cell references (the 13-gene mean and the CD8A-PRF1-GZMA mean). `score_signature` skips missing values
  sample by sample.
- **`barriernet benchmark`** (`barriernet.benchmark.benchmark`): pooled score × arm interaction per signature, a
  minimum number of comparisons for pooling, size-matched random genes with the empirical P in the declared direction,
  paired differences from the T-cell references and single genes, Benjamini-Hochberg q values within each family of
  tests, and a summary of how many signatures pass each test. Multi-process bootstrap and random reference (`--jobs`).
- **`barriernet singlearm`** (`barriernet.singlearm`): log-odds of response per SD in single-arm cohorts, pooled
  by DerSimonian-Laird with the wider of the Hartung-Knapp and DerSimonian-Laird intervals, and the Spearman agreement
  with the randomized interaction.
- **`barriernet published`** (`barriernet.published`): pooling of published subgroup hazard ratios (signature high
  and low) by REML with the Knapp-Hartung adjustment, per-SD conversion from the reported split, leave-one-out, Egger
  test, PD-L1 and TMB over the same trials; `example_table()` holds the 19 comparisons extracted for the benchmark.
- **`barriernet design`** (`barriernet.design`): events needed, power and detectable difference for telling two
  correlated scores apart.
- **Report** (`barriernet.report`, `barriernet report`): figures and an HTML page and PDF; needs the `report` extra.
- **Web app** (`barriernet app`, `barriernet.app`, `shiny_app/`): needs the `app` extra (Shiny for Python). The benchmark
  runs in a separate process, so the page stays responsive and can be cancelled; one benchmark at a time, with limits on
  workers, draws and uploads for a shared server. `shiny_app/make_bundle.sh` makes a folder for Shiny Server.
- `barriernet demo` and `barriernet.demo`: synthetic randomized trials and single-arm cohorts with a known interaction.
- `barriernet library`: export of the bundled library.

### Changed
- `calibrate.bootstrap_indices` accepts a sequence of integers as the seed.
- `barriernet calibrate` and `benchmark` share the trial loader.

### Checked
Run on the six randomized comparisons of the benchmark, version 0.4.0 reproduces the analysis code: pooled interactions
within 1e-5 for all 51 signatures; the same counts of signatures favouring either arm, beyond random genes (25; 24 after
correction) and above or below the T-cell references and CD8A; paired differences equal to 1e-4 for 50 of 51 signatures
(CRMA_Shukla differs because its score takes 30 values in 726 patients of one trial and a few resamples give
near-degenerate Cox fits); single-arm log-odds and pooled values within 1e-6; published-trial pooling within 5e-5; and
the design formula.

## 0.3.0 (2026-10-31)

### Added
- **barriernet 0.3 weighted score** (`barriernet03`, `data/geneset_rna_v03.yaml`). The 48 genes of 0.2.0 with module
  weights (IFN +0.362, APM +0.225, TLS_B +0.126, T_CELL +0.104, CAF_ECM −0.160, TAM_SUP −0.023; every gene of a module
  weighted equally; sum of |w| = 1). Developed only on 20 single-arm checkpoint cohorts (ridge logistic regression with
  cohort intercepts, chosen by leave-one-cohort-out validation) and frozen before any randomized trial was analyzed.
  `barriernet bulk` now adds the `barriernet03` column (`--no-v03` to omit; `--weights` for another weight file);
  `core.score_weighted` and `genesets.load_weights` give the same in Python.
- **`barriernet.calibrate`** and **`barriernet calibrate --config trials.yaml`**: score × arm interaction in randomized
  comparisons (Cox with Efron ties, or logistic with a Firth fallback on separation), random-effects pooling by the
  conservative rule (wider of Hartung–Knapp and DerSimonian–Laird, larger P), size-matched random gene sets as the
  reference for specificity, a within-arm paired bootstrap against single genes, and the detectable difference.
- `calibrate.events_needed`: events a randomized survival comparison needs to detect a given difference between the
  treatment interactions of two correlated scores (normal theory, optional inflation).

### Changed
- `scipy` is now a dependency.

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
