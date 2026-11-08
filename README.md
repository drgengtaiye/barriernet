# barriernet

Marker-level Barrier / Immune / NetImmune scoring for spatial, single-cell, bulk and protein data, a framework that
calibrates predictive claims for expression scores in randomized trials (`barriernet calibrate`), the benchmark of
published immune signatures against simple T-cell scores with a web app (`barriernet benchmark`, `barriernet app`,
from version 0.4.0), and the analysis code behind the manuscript *Benchmarking published immune signatures against
simple T-cell scores in randomized checkpoint-inhibitor trials*.

## What is here

| Path | Contents |
|---|---|
| `barriernet/` | the Python package, version 0.4.0 (see `barriernet/CHANGELOG.md`): the library of 51 published signatures, the benchmark against T-cell references, single genes and random genes, single-arm cohorts, pooling of published subgroup results, trial design, an HTML/PDF report and the web app; the frozen barriernet 0.3 weights and the calibration framework of 0.3.0; 0.2.1 corrected the Visium neighborhood of 0.2.0, and every spatial result in the paper was computed on the corrected lattice |
| `barriernet/shiny_app/` | the web app (`app.py`, Shiny for Python), its requirements and `make_bundle.sh`, which makes a self-contained folder for Shiny Server |
| `barriernet-0.4.0-py3-none-any.whl`, `barriernet_0.4.0.zip` | built wheel and source archive of 0.4.0 (0.3.0 and 0.2.1 kept beside them) |
| `analysis/15_v19_v41/protocols/` | the later analysis protocols with their SHA-256 hashes and outcomes, the frozen barriernet 0.3 weights' hash, and the **analysis plan for the controlled-access randomized trials** (`E37_registration.md` and `E37_registration_addendum.md`) |
| `analysis/` | every analysis and figure script behind the manuscript, grouped by stage |
| `analysis/16_v42_v46/` | the protocols with their hashes, outcomes, scripts and figure scripts added after the v0.3.0 release (published randomized trials, barriernet 0.3 against each signature, absolute benefit, gastric and gastroesophageal cancer) |
| `analysis/17_v47_v48/` | the protocols with their hashes, outcomes, scripts and figure scripts of the current manuscript (a reference gene chosen in single-arm data, a coexpression-matched random reference, simple T-cell means, the MSS/EBV-negative gastric patients) |
| `analysis/18_barriernet04/` | the check that barriernet 0.4.0 reproduces the manuscript's analyses from the same inputs (78 of 78 comparisons with the analysis results) |
| `analysis/README_v48_map.md` | **the map from each table, figure panel and reported number of the current manuscript to the script that produces it** (`analysis/README.md`: the earlier version) |

## The score

Six marker modules of 48 genes in total:

| Module | Genes |
|---|---|
| `CAF_ECM` | COL1A1, COL1A2, LUM, DCN, POSTN, TAGLN, ACTA2, THY1, FAP, FN1, COL3A1 |
| `TAM_SUP` | SPP1, APOE, TREM2, C1QA, C1QB, C1QC, LGALS3, MARCO, IL10, SOCS3, STAT3 |
| `APM` | B2M, HLA-A, HLA-B, HLA-C, TAP1, TAP2, PSMB8, PSMB9 |
| `IFN` | STAT1, IRF1, CXCL9, CXCL10, IFNG |
| `TLS_B` | CXCL13, LTB, MS4A1, CD79A, CD74, MZB1 |
| `T_CELL` | CD3D, CD3E, TRAC, TRBC1, CD8A, GZMB, NKG7 |

A module's score is the mean of per-gene z-scores over the genes present in the dataset, provided at
least three are present. The axes are `Barrier = mean(CAF_ECM, TAM_SUP)`,
`Immune = mean(APM, IFN, TLS_B, T_CELL)` and `NetImmune = Immune − Barrier`. Genes are z-scored
within each dataset on a log2 scale and never across datasets. Nothing is fitted to outcome: the
gene lists, the weights within and across modules and the sign of the subtraction are all fixed.

## barriernet 0.3 and its analysis plan

barriernet 0.3 keeps the 48 genes and gives each module a weight learned only in 20 single-arm checkpoint cohorts
(IFN +0.362, APM +0.225, TLS_B +0.126, T_CELL +0.104, CAF_ECM −0.160, TAM_SUP −0.023; every gene of a module weighted
equally). The weights were frozen, and their SHA-256 hash logged, before any randomized trial was scored:

| File | SHA-256 |
|---|---|
| `analysis/15_v19_v41/results/E37_frozen_weights.tsv` (the same weights ship in `barriernet/barriernet/data/geneset_rna_v03.yaml`) | `500264cc627b3ef0240c0340296dc001bceb02748e3f3b0e83892b234b2f63bd` |
| `analysis/15_v19_v41/protocols/E37_registration.md` (analysis plan) | `c24bb4141a734c103f8e3f30681bfd39e2a191ac95bbe54713b028f9e962abc5` |
| `analysis/15_v19_v41/protocols/E37_registration_addendum.md` (events needed, reporting rules) | `c07e4e8cd9e663c249c43ed90b7a8a7de0afc1639198a26027aa248d3d7157f0` |

The analysis plan fixes the primary analysis of barriernet 0.3 in the larger randomized trials whose expression is held
under controlled access (OAK and POPLAR, IMvigor010, IMmotion151, CheckMate 649). It was made public here before any of
those datasets was requested or obtained; the result will be reported whatever it is.

## Install

```bash
pip install barriernet-0.3.0-py3-none-any.whl     # or, from the source tree:
cd barriernet && pip install -e .                 # add [spatial] for anndata and scanpy
```

It also runs without installing: `cd barriernet && python -m barriernet --help`.

## Use

The same three subcommands score bulk expression, IHC/IF panels and Visium sections:

```bash
barriernet bulk    --expr expr.tsv --out out_bulk       # samples x genes or genes x samples, log2 scale
barriernet ihc     --panel panel.tsv --out out_ihc
barriernet spatial --h5ad section.h5ad --out out_spatial   # adds TLS ring/peak summaries
```

Each writes the six module scores and the Barrier, Immune and NetImmune axes (and, for bulk, the
`barriernet03` score), together with the genes actually used. `barriernet calibrate --config trials.yaml` tests scores
against treatment interaction, size-matched random genes and single genes over randomized comparisons. `--log2p1` log-transforms linear input; a warning is printed when values look
un-logged. See `barriernet/README.md` for the full options and output files.

## Tests

```bash
cd barriernet && python -m pytest -q tests
```

## Reproducing the manuscript

`analysis/README.md` gives the run order and maps every item in the manuscript to its script. The
scripts keep the absolute input and output paths they ran under, so paths must be adjusted to a
local copy of the public datasets listed in the manuscript's data-availability statement.

## License

MIT; see `LICENSE`.
