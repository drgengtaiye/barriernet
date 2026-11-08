# barriernet v0.4.1

> 0.4.1 puts the web app in the browser at https://drgengtaiye.github.io/barriernet/ (CHANGELOG.md). 0.4.0 adds the analyses of a benchmark of published immune signatures in randomized checkpoint-inhibitor trials:
> `barriernet benchmark` (51 published signatures against simple T-cell scores, single genes and random genes),
> `singlearm`, `published` (pooling published subgroup hazard ratios), `design` (events needed), an HTML/PDF report and
> a web app (`barriernet app`); see sections 5 and 6 and CHANGELOG.md. 0.3.0 added the barriernet 0.3 weighted score and
> `barriernet calibrate`. 0.2.1 fixed the Visium neighbourhood used for TLS-centred rings; spatial results produced with
> 0.2.0 should be recomputed.

A portable, interpretable scoring toolkit that operationalizes **Barrier** and **NetImmune** axes across:

- **IHC/IF panels** (`barriernet ihc`)
- **Bulk RNA / bulk proteomics** (`barriernet bulk`)
- **Spatial platforms (Visium)** with ring/peak summaries when coordinates exist (`barriernet spatial`)

## Install (optional)
The report needs matplotlib (`pip install 'barriernet[report]'`) and the web app Shiny for Python
(`pip install 'barriernet[app]'`). You can run without installing:

```bash
python -m barriernet --help
```

Or install editable:

```bash
pip install -r requirements.txt
pip install -e .
```

## Quick start

### 1) Bulk RNA / proteomics (samples×genes or genes×samples)
```bash
barriernet bulk --expr expr.tsv --out out_bulk
```

Each module score is the mean of per-gene z-scores (across samples) over the module genes
present, if at least 3 are present; Barrier = mean(CAF_ECM, TAM_SUP), Immune = mean(APM, IFN,
TLS_B, T_CELL), NetImmune = Immune − Barrier. Supply log-scale expression (log2 TPM+1 or
normalised array values) or add `--log2p1`; a warning is printed if values look un-logged.

Outputs:
- `barriernet_bulk_scores.tsv` (modules + Barrier/Immune/NetImmune)
- `barriernet_bulk_used_genes.tsv`
- `barriernet_bulk_panel_matrix.tsv` (only genes actually used)

### 2) IHC / IF panel table (ROI/cell-level allowed; aggregated to sample)
```bash
barriernet ihc --table ihc_table.tsv --sample-col sample_id --out out_ihc
```

Outputs:
- `barriernet_ihc_scores.tsv`
- `barriernet_ihc_used_markers.tsv`

### 3) Spatial Visium (uses `obs['array_row']`, `obs['array_col']`, `obsm['spatial']`)
```bash
barriernet spatial --h5ad-glob "sections/*.h5ad" --out out_spatial            # primary: Far >= 3 hops
barriernet spatial --h5ad-glob "sections/*.h5ad" --out out_sens \
    --far-min-dist 4 --far-exclude-edge                                     # sensitivity
```

Neighbourhood: Visium array coordinates are doubled-column, so a spot's six neighbours
are at (0, ±2) and (±1, ±1) (`--coords visium`, default). Use `--coords axial` only for true
axial hexagonal coordinates. Tissue edge = fewer than six in-tissue neighbours.

Zones: Hotspot (TLS peak), Ring1 (1 hop), Ring2 (2 hops), Far (`--far-min-dist`, default 3).
TLS peaks: spots at or above the `--peak-q` quantile (0.95) that are local maxima among the
`--peak-k` nearest spots (9, the spot included; `--peak-method knn`) or among the lattice
neighbours (`--peak-method lattice`).

Outputs per input:
- `<stem>_spot_rings.tsv` (zone, hops to TLS peak and to tissue edge, lattice degree)
- `<stem>_ring_enrichment.tsv` (zone means; Ring1−Far and Ring2−Far contrasts)
- `<stem>.with_barriernet.h5ad` (ring metrics written to obs)

Across inputs (three or more): `barriernet_sections.tsv` with the contrasts and BarrierIndex
= z(CAF_ECM Ring1−Far) + z(SPP1 TAM Ring1−Far), standardised across the sections given.


### 4) Calibrating a predictive claim over randomized trials
```bash
barriernet calibrate --config trials.yaml --out out_cal
```
`trials.yaml` lists each randomized comparison:
```yaml
trials:
  - name: TRIAL1
    expr: trial1_expr.tsv          # log-scale expression, samples x genes or genes x samples
    clinical: trial1_clinical.tsv  # one row per sample
    sample_col: sample_id
    arm_col: arm                   # 1 = checkpoint-containing arm, 0 = control
    time_col: pfs_months           # survival endpoint ...
    event_col: pfs_event
  - name: TRIAL2
    expr: trial2_expr.tsv
    clinical: trial2_clinical.tsv
    arm_col: arm
    response_col: pcr              # ... or a binary benefit endpoint (1 = benefit)
    covariates: [hormone_receptor] # kept in the logistic model
n_random: 1000                     # size-matched random gene sets
n_boot: 300                        # within-arm bootstrap resamples
scores: {}                         # optional extra scores, name -> {gene: weight}
```
Genes are z-scored within each trial across both arms. Scored by default: `barriernet03`, `NetImmune` and the 27 single
genes of the T-cell and interferon-γ panels. Outputs: `calibrate_per_trial.tsv` (interaction per SD, benefit-oriented:
−log HR or log OR), `calibrate_pooled.tsv` (pooled interaction; random-gene reference median, 95th percentile,
empirical P) and `calibrate_paired.tsv` (difference of each score from each single gene, P and detectable difference).
The same is available in Python as `barriernet.calibrate.calibrate(trials, scores, genes)`.

Planning a trial: `barriernet.calibrate.events_needed(delta, rho)` gives the events a 1 : 1 randomized survival comparison
needs to detect a difference `delta` per SD between the interactions of two scores correlated at `rho` (80% power,
two-sided 0.05, normal theory; pass `inflation=1.17` for the excess variance observed in the manuscript's survival
comparisons). For example, `events_needed(0.10, 0.81)` ≈ 1,190 for barriernet 0.3 against CD8A.

### 5) Benchmarking published signatures in randomized trials (0.4)

`barriernet benchmark` runs the whole benchmark in one command: the 51 distinct published immuno-oncology signatures
bundled with the package (and any of your own) are scored in every randomized comparison, the score × arm interaction
is pooled over the comparisons, and each signature is set against the 13-gene T-cell mean, the CD8A-PRF1-GZMA mean,
the 27 genes of the T-cell and interferon-γ panels and size-matched random gene sets, with Benjamini-Hochberg q values
within each family of tests.

```bash
barriernet benchmark --config trials.yaml --out out_bm --jobs -1          # same YAML as calibrate
barriernet benchmark --config trials.yaml --signatures mine.tsv --out out_bm
```

The signature table is long: one row per gene with columns `signature`, `gene` and optionally `weight` (default 1),
`direction` (1 for a sensitivity signature, −1 for resistance; default 1) and `family`. A score is
sum(w z) / sum(|w|) over the listed genes present (missing values skipped sample by sample), computed when at least
three genes and half the list are present; a signature is pooled when scored in at least 60% of the comparisons
(`--min-trials`; four of six in the manuscript). The random reference draws, in each comparison, as many genes as the
signature uses there, with its weights in random order, and the empirical P is taken in the signature's declared
direction. Paired differences use within-arm patient resamples common to the two scores. Outputs:
`benchmark_pooled.tsv`, `benchmark_paired.tsv`, `benchmark_per_trial.tsv`, `benchmark_summary.tsv` (how many
signatures pass each test, nominally and after correction) and `report.html` / `report.pdf` with the figures.
`barriernet library --out signatures.tsv` writes the bundled library.

Single-arm checkpoint cohorts (response only, no control arm):

```bash
barriernet singlearm --config cohorts.yaml --randomized out_bm/benchmark_pooled.tsv --out out_sa
```

```yaml
cohorts:
  - name: COHORT1
    expr: cohort1_expr.tsv
    clinical: cohort1_clinical.tsv
    sample_col: sample_id
    response_col: responder        # 1 = responder, 0 = non-responder
```

Each score's log-odds of response per SD is estimated cohort by cohort (logistic regression, Firth's correction when the
fit separates) and pooled by DerSimonian-Laird with the wider of the Hartung-Knapp and DerSimonian-Laird intervals,
when scored in at least half the cohorts. With `--randomized`, the Spearman correlation between the single-arm and the
randomized rankings is reported. A single-arm association cannot separate a predictive score from a prognostic one.

Published subgroup results of randomized trials:

```bash
barriernet published --table my_trials.tsv --out out_pub
barriernet published --example --out out_pub      # the 19 comparisons extracted for the benchmark
```

One row per trial with `trial`, `endpoint` (OS or PFS are pooled) and the hazard ratio of the checkpoint arm against
control with its 95% CI in the signature-high (`hr_high`, `lo_high`, `hi_high`) and signature-low (`hr_low`, `lo_low`,
`hi_low`) groups; optional `frac_high` (share classed high, for the per-SD conversion), `from_figure`, `exclude`, and
PD-L1 or TMB subgroup columns. The log ratio of the two hazard ratios is pooled by REML with the Knapp-Hartung adjustment
(as R metafor `rma(method = "REML", test = "knha")`), with leave-one-out and Egger tests. The example values were
extracted from the trial reports; check them against the sources before relying on them.

Trial design:

```bash
barriernet design --delta 0.10 --rho 0.55 0.81 0.84 --events 800 --out out_design
```

gives the events a randomized survival comparison needs to tell two scores apart, 2(1 − r)(z₁₋α/₂ + z_power)² /
(δ² p(1 − p)), and the power and detectable difference at a given number of events.

`barriernet demo --out out_demo` runs every analysis on synthetic data and writes the report; the synthetic trials and
cohorts are simulated, not patient data. `barriernet report --dir out_dir` rebuilds the report from saved tables.

In Python:

```python
from barriernet.benchmark import benchmark
from barriernet.singlearm import Cohort, single_arm, concordance
from barriernet.published import pool_published, example_table
from barriernet.design import events_needed, power, detectable_difference
from barriernet.report import write_report

B = benchmark(trials, n_random=1000, n_boot=300, n_jobs=-1)      # trials: list of barriernet.calibrate.Trial
write_report("out", benchmark=B, published=pool_published(example_table()))
```

### 6) Web app

The app runs in the browser at https://drgengtaiye.github.io/barriernet/, with nothing to install; there the benchmark shows a demonstration
computed in advance and runs your own data in the page, more slowly than on a computer. To run it locally:

```bash
pip install 'barriernet[app]'
barriernet app                      # http://127.0.0.1:8000
```

The app has the benchmark (synthetic demonstration or a zip holding `trials.yaml` and the tables it names), single-arm
cohorts, an editable table of published trials whose pooled result updates as cells change, the trial-design
calculator and the signature library, with every result downloadable. Data stay on the computer running the app. To host it with Shiny Server (Open Source 1.5.21 or later),
`shiny_app/make_bundle.sh` makes a self-contained folder; see `shiny_app/README.md` for the server's Python and the
limits for a shared server.

## Default must-appear panels (YAML)
- RNA/spatial: `barriernet/data/geneset_rna.yaml` — the 48 marker genes of the six modules
  (plus an auxiliary STF3 module that enters no axis).
- IHC/protein: `barriernet/data/geneset_ihc.yaml` — a reduced protein panel (42 markers:
  without THY1, SOCS3, IFNG, LTB, TRBC1, NKG7); `barriernet ihc` defaults to `--min-genes 2`.

## Tests
```bash
python -m pytest -q tests
```

Override with your own:
```bash
barriernet bulk --expr expr.tsv --geneset my_geneset.yaml --out out_bulk
```
