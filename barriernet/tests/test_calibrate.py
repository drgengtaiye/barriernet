"""barriernet 0.3 weights and the calibration module."""
import numpy as np, pandas as pd, pytest
from barriernet.genesets import load_yaml, default_yaml, load_weights
from barriernet.core import score_modules, score_weighted
from barriernet.calibrate import Trial, cox_ph, interaction, pool, mdd, events_needed, weighted_score, calibrate, logistic

TOTALS = {"IFN": 0.3617, "APM": 0.2251, "TLS_B": 0.1261, "T_CELL": 0.1042, "CAF_ECM": -0.1597, "TAM_SUP": -0.0232}


def test_v03_weights_are_the_frozen_module_weights():
    wt = load_weights(); gs = load_yaml(default_yaml("rna")); w = wt["weights"]
    assert len(w) == 48 and abs(sum(abs(v) for v in w.values()) - 1) < 1e-12
    for m, tot in TOTALS.items():
        vals = [w[g] for g in gs.modules[m]]
        assert np.allclose(vals, vals[0]) and abs(sum(vals) - tot) < 1e-4, m


def test_weighted_score_is_the_module_combination_when_all_genes_present():
    gs = load_yaml(default_yaml("rna")); wt = load_weights()
    genes = sorted(wt["weights"]); X = pd.DataFrame(np.random.default_rng(0).normal(size=(40, 48)), columns=genes)
    s, used = score_weighted(X, wt["weights"]); mods, _ = score_modules(X, gs.modules)
    comb = sum(sum(wt["weights"][g] for g in gs.modules[m]) * mods[m] for m in TOTALS)
    assert len(used) == 48 and np.allclose(s, comb)


def test_weighted_score_needs_half_the_genes():
    wt = load_weights(); genes = sorted(wt["weights"])[:20]
    X = pd.DataFrame(np.random.default_rng(1).normal(size=(30, 20)), columns=genes)
    s, _ = score_weighted(X, wt["weights"]); assert s.isna().all()


def test_cox_matches_lifelines_with_ties():
    lifelines = pytest.importorskip("lifelines")
    rng = np.random.default_rng(2); n = 300; X = rng.normal(size=(n, 3))
    t = np.round(rng.exponential(np.exp(-X @ [0.4, -0.2, 0.1])), 1) + 0.1; e = rng.random(n) < 0.7
    b, se = cox_ph(X, t, e)
    df = pd.DataFrame(X, columns=["a", "b", "c"]); df["T"] = t; df["E"] = e
    r = lifelines.CoxPHFitter().fit(df, "T", "E", fit_options={"precision": 1e-12}).summary
    assert np.allclose(b, r["coef"], atol=1e-6) and np.allclose(se, r["se(coef)"], atol=1e-6)


def _trial(n, beta, kind, seed):
    rng = np.random.default_rng(seed); arm = rng.integers(0, 2, n).astype(float)
    Z = pd.DataFrame(rng.normal(size=(60, n)), index=[f"G{i}" for i in range(60)])
    s = Z.loc["G0"].to_numpy(); eta = beta * arm * (s - s.mean()) / s.std(ddof=1)
    if kind == "cox":
        E = rng.exponential(np.exp(eta)); C = rng.exponential(2.0, n)      # benefit = lower hazard
        return Trial(f"T{seed}", Z, arm, time=np.minimum(E, C), event=(E <= C).astype(int))
    return Trial(f"T{seed}", Z, arm, response=(rng.random(n) < 1 / (1 + np.exp(-(eta - 0.5)))).astype(float))


@pytest.mark.parametrize("kind", ["cox", "logit"])
def test_interaction_recovers_the_simulated_effect(kind):
    tr = _trial(4000, 0.5, kind, 3); est, se = interaction(tr, tr.Z.loc["G0"].to_numpy())
    assert abs(est - 0.5) < 4 * se


def test_logistic_falls_back_to_firth_on_separation():
    X = np.column_stack([np.ones(20), np.r_[np.zeros(10), np.ones(10)]]); y = np.r_[np.zeros(10), np.ones(10)]
    b, se = logistic(X, y); assert np.isfinite(b).all() and np.isfinite(se).all()


def test_pool_and_mdd():
    q = pool([0.2, 0.2, 0.2], [0.1, 0.1, 0.1]); assert abs(q["est"] - 0.2) < 1e-9 and q["k"] == 3
    assert abs(mdd(0.076) - 0.213) < 0.001                                   # the manuscript's 0.21 per SD


def test_events_needed():
    assert abs(events_needed(0.10, 0.811) - 1186.8) < 1.0                    # barriernet 0.3 vs CD8A (E39)
    assert abs(events_needed(0.10, 0.811, inflation=1.17) - 1388.5) < 1.0
    assert events_needed(0.10, 0.5) > events_needed(0.10, 0.8) > events_needed(0.20, 0.8)


def test_calibrate_runs_end_to_end():
    trials = [_trial(300, 0.4, k, i) for i, k in enumerate(["cox", "logit", "cox"])]
    res = calibrate(trials, {"S": {"G0": 1.0, "G1": 0.5, "G2": 0.5}}, ["G0", "G5"], n_random=20, n_boot=20, seed=1)
    assert set(res) == {"per_trial", "pooled", "paired"}
    p = res["pooled"].set_index("score"); assert p.loc["S", "k"] == 3 and "emp_p" in p.columns
    assert set(res["paired"].gene) == {"G0", "G5"}
