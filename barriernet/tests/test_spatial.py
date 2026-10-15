"""Lattice, ring and regression tests for barriernet.spatial (0.2.1)."""
import glob, os
import numpy as np, pandas as pd, pytest
from barriernet.spatial import (build_hex_graph, lattice_degree, assign_rings, AXIAL_STEPS,
                                VISIUM_STEPS, barrier_index)

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2"
FIG = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures"


def visium_lattice(nrow=10, ncol=20):
    rc = [(r, c) for r in range(nrow) for c in range(ncol) if (r + c) % 2 == 0]
    return pd.DataFrame({"array_row": [r for r, _ in rc], "array_col": [c for _, c in rc], "in_tissue": 1})


def test_synthetic_visium_interior_has_six_neighbours():
    obs = visium_lattice()
    nbrs, keep = build_hex_graph(obs)
    deg = lattice_degree(nbrs, keep)
    r, c = obs.array_row.values, obs.array_col.values
    interior = (r > 0) & (r < r.max()) & (c > 1) & (c < c.max() - 1)
    assert interior.sum() > 20
    assert (deg[interior] == 6).all()
    assert (deg[~interior] < 6).all()


def test_legacy_axial_offsets_break_on_visium():
    obs = visium_lattice()
    with pytest.warns(UserWarning, match="coordinate system is probably wrong"):
        nbrs, keep = build_hex_graph(obs, coords="axial")
    assert lattice_degree(nbrs, keep).max() <= 2


def test_axial_option_on_axial_coordinates():
    rc = [(q, r) for q in range(8) for r in range(8)]
    obs = pd.DataFrame({"array_row": [a for a, _ in rc], "array_col": [b for _, b in rc]})
    nbrs, keep = build_hex_graph(obs, coords="axial")
    assert lattice_degree(nbrs, keep).max() == 6


def test_rings_on_synthetic_lattice():
    obs = visium_lattice(12, 24)
    obs["TLS"] = 0.0
    centre = (obs.array_row == 6) & (obs.array_col == 12)
    rd = assign_rings(obs, np.c_[obs.array_row, obs.array_col].astype(float), "TLS",
                      peaks=centre.values, far_min_dist=3)
    assert (rd.ring == "Hotspot").sum() == 1
    assert (rd.ring == "Ring1").sum() == 6
    assert (rd.ring == "Ring2").sum() == 12


def test_barrier_index_is_sum_of_zscores():
    t = pd.DataFrame({"ECM_Ring1_minus_Far": [1.0, 2.0, 3.0], "TAM_Ring1_minus_Far": [0.0, 0.0, 3.0]})
    bi = barrier_index(t)
    z = lambda v: (v - v.mean()) / v.std(ddof=0)
    assert np.allclose(bi, z(t.iloc[:, 0]) + z(t.iloc[:, 1]))


real = sorted(glob.glob(f"{R}/*.with_state_and_ring.h5ad"))


@pytest.mark.skipif(not real, reason="section h5ad files not available")
def test_real_section_mean_degree():
    import anndata as ad
    o = ad.read_h5ad(real[1], backed="r").obs
    nbrs, keep = build_hex_graph(o)
    deg = lattice_degree(nbrs, keep)[keep]
    assert 5.6 < deg.mean() < 6.0


@pytest.mark.parametrize("sec", ["GSM7990474_21_00731_LI_SING", "GSM7990478_21_01251_LI_SING",
                                 "GSM7990479_21_01252_LI_SING"])
def test_regression_against_ring_recompute(sec):
    ref = f"{FIG}/ring_hex_{sec}.csv"
    h5 = f"{R}/{sec}.with_state_and_ring.h5ad"
    if not (os.path.exists(ref) and os.path.exists(h5)):
        pytest.skip("reference files not available")
    import anndata as ad
    a = ad.read_h5ad(h5, backed="r"); o = a.obs.copy()
    rd = assign_rings(o, np.asarray(a.obsm["spatial"], float), "TLS_niche_score",
                      peaks=o.TLS_peaks.astype(int).values.astype(bool), far_min_dist=3)
    exp = pd.read_csv(ref, index_col=0).ring_hex.reindex(o.index)
    assert (rd.ring.values == exp.values).all()


@pytest.mark.skipif(not real, reason="section h5ad files not available")
def test_stored_peaks_are_legacy_lattice_maxima():
    """Documents the second 0.2.0 artefact: stored TLS peaks = maxima against <= 2 diagonal neighbours."""
    import anndata as ad
    a = ad.read_h5ad(real[0], backed="r"); o = a.obs
    s = o.TLS_niche_score.values.astype(float); thr = np.quantile(s, 0.95)
    rc = list(zip(o.array_row.astype(int), o.array_col.astype(int))); pos = {p: i for i, p in enumerate(rc)}
    pk = np.array([s[i] >= thr and all(s[i] >= s[pos[(r + a_, c + b_)]] for a_, b_ in AXIAL_STEPS
                                        if (r + a_, c + b_) in pos) for i, (r, c) in enumerate(rc)])
    assert (pk == o.TLS_peaks.astype(int).values.astype(bool)).all()
