"""Gene sets and scoring match the manuscript definitions."""
import numpy as np, pandas as pd, pytest
from barriernet.genesets import load_yaml, default_yaml
from barriernet.core import score_modules, score_axes

MANUSCRIPT = {
    "CAF_ECM": "COL1A1 COL1A2 LUM DCN POSTN TAGLN ACTA2 THY1 FAP FN1 COL3A1",
    "TAM_SUP": "SPP1 APOE TREM2 C1QA C1QB C1QC LGALS3 MARCO IL10 SOCS3 STAT3",
    "APM": "B2M HLA-A HLA-B HLA-C TAP1 TAP2 PSMB8 PSMB9",
    "IFN": "STAT1 IRF1 CXCL9 CXCL10 IFNG",
    "TLS_B": "CXCL13 LTB MS4A1 CD79A CD74 MZB1",
    "T_CELL": "CD3D CD3E TRAC TRBC1 CD8A GZMB NKG7",
}


def test_rna_geneset_is_the_48_manuscript_genes():
    gs = load_yaml(default_yaml("rna"))
    for m, g in MANUSCRIPT.items():
        assert gs.modules[m] == g.split(), m
    assert sum(len(g.split()) for g in MANUSCRIPT.values()) == 48


def test_axes_definitions():
    gs = load_yaml(default_yaml("rna"))
    rng = np.random.default_rng(0)
    genes = sorted({x for g in MANUSCRIPT.values() for x in g.split()})
    X = pd.DataFrame(rng.normal(size=(30, len(genes))), columns=genes)
    mods, used = score_modules(X, gs.modules)
    ax = score_axes(mods, gs.axes)
    Z = (X - X.mean()) / X.std(ddof=0)
    for m, g in MANUSCRIPT.items():
        assert np.allclose(mods[m], Z[g.split()].mean(axis=1), atol=1e-6)
    barrier = mods[["CAF_ECM", "TAM_SUP"]].mean(axis=1)
    immune = mods[["APM", "IFN", "TLS_B", "T_CELL"]].mean(axis=1)
    assert np.allclose(ax["Barrier"], barrier) and np.allclose(ax["Immune"], immune)
    assert np.allclose(ax["NetImmune"], immune - barrier)


def test_module_needs_three_genes():
    X = pd.DataFrame(np.random.default_rng(1).normal(size=(10, 2)), columns=["IFNG", "STAT1"])
    mods, _ = score_modules(X, {"IFN": MANUSCRIPT["IFN"].split()})
    assert mods["IFN"].isna().all()


def test_orientation_detected_by_gene_overlap(tmp_path):
    from barriernet.io import load_expr_matrix
    genes = MANUSCRIPT["APM"].split() + ["G%d" % i for i in range(40)]
    df = pd.DataFrame(np.random.default_rng(2).normal(size=(200, len(genes))), columns=genes)  # 200 samples x 48 genes
    df.index = ["S%d" % i for i in range(200)]
    p = tmp_path / "x.tsv"; df.to_csv(p, sep="\t")
    expr, orient = load_expr_matrix(str(p))
    assert orient == "samples_by_genes" and expr.shape == (200, len(genes))
