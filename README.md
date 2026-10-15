# barriernet

Marker-level Barrier / Immune / NetImmune scoring for spatial, single-cell, bulk and protein
data, together with the analysis code behind the manuscript *Stroma tracks prognosis and immune
programs track checkpoint response: an interpretable tumor-microenvironment score tested in ten
gastric cohorts, 21 checkpoint-treated cohorts and six randomized comparisons*.

## What is here

| Path | Contents |
|---|---|
| `barriernet/` | the Python package, version 0.2.1 (see `barriernet/CHANGELOG.md`; 0.2.1 corrects the Visium neighborhood of 0.2.0, and every spatial result in the paper was computed on the corrected lattice) |
| `barriernet-0.2.1-py3-none-any.whl`, `barriernet_0.2.1.zip` | built wheel and source archive of that version |
| `analysis/` | every analysis and figure script behind the manuscript, grouped by stage |
| `analysis/README.md` | **the run order and the map from each table, figure and reported number to the script that produces it** |

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

## Install

```bash
pip install barriernet-0.2.1-py3-none-any.whl     # or, from the source tree:
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

Each writes the six module scores and the Barrier, Immune and NetImmune axes, together with the
genes actually used. `--log2p1` log-transforms linear input; a warning is printed when values look
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
