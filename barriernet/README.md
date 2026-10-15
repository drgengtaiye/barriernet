# barriernet v0.2.1

> 0.2.1 fixes the Visium neighbourhood used for TLS-centred rings (see CHANGELOG.md). Spatial results produced with 0.2.0 should be recomputed.

A portable, interpretable scoring toolkit that operationalizes **Barrier** and **NetImmune** axes across:

- **IHC/IF panels** (`barriernet ihc`)
- **Bulk RNA / bulk proteomics** (`barriernet bulk`)
- **Spatial platforms (Visium)** with ring/peak summaries when coordinates exist (`barriernet spatial`)

## Install (optional)
You can run without installing:

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
