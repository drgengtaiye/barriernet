"""Tests of the 0.4 modules: library, benchmark, single-arm, published trials, design, report, bundle and CLI."""
import io
import zipfile

import numpy as np
import pandas as pd
import pytest

from barriernet.benchmark import benchmark, bh
from barriernet.calibrate import events_needed
from barriernet.demo import demo_cohorts, demo_trials
from barriernet.design import curve, detectable_difference, power
from barriernet.library import (IFN_CHECKPOINT_PANEL, MAIN_REFERENCE, T_CELL_PANEL, load_library, read_signature_table,
                                reference_weights, score_signature, signature_meta, signature_weights)
from barriernet.published import example_table, per_sd_divisor, pool_published, rma
from barriernet.singlearm import concordance, pool_dl_hk, single_arm


def test_library_counts():
    L, A = load_library(), load_library(distinct=False)
    assert L.signature.nunique() == 51 and A.signature.nunique() == 52 and len(A) == 2112
    M = signature_meta(A)
    assert M.family.value_counts().to_dict() == {"Other": 24, "Inflamed/IFN": 20, "Confounder": 8}
    assert int((M.direction > 0).sum()) == 34
    W = signature_weights(A)
    assert set(W["TIS_Damotte"]) == set(W["Tcell-inflamed_Ayers"]) and len(W["TIS_Damotte"]) == 18
    assert len(T_CELL_PANEL) == 13 and len(IFN_CHECKPOINT_PANEL) == 14
    assert set(reference_weights()) == {MAIN_REFERENCE, "T-cell mean (CD8A, PRF1, GZMA)"}


def test_score_skips_missing_values_per_sample():
    Z = pd.DataFrame({"s1": [1.0, 2.0, np.nan, 0.0], "s2": [0.0, -1.0, 1.0, 0.0], "s3": [-1.0, -1.0, -1.0, 0.0]}, index=["A", "B", "C", "D"])
    sc, genes, w = score_signature(Z, {"A": 1.0, "B": -2.0, "C": 1.0, "D": 1.0})
    assert genes == ["A", "B", "C"]                                    # D is constant
    assert sc[0] == pytest.approx((1 - 4) / 3)                          # C missing in s1: divide by |1| + |-2|
    assert sc[1] == pytest.approx((0 + 2 + 1) / 4)
    assert score_signature(Z, {"A": 1, "X": 1, "Y": 1, "Z": 1})[0] is None   # fewer than three genes


def test_bh_matches_hand_values():
    q = bh([0.01, 0.04, 0.03, 0.20, np.nan])
    assert np.allclose(q[:4], [0.04, 0.16 / 3, 0.16 / 3, 0.20]) and np.isnan(q[4])


def test_read_signature_table(tmp_path):
    p = tmp_path / "mine.tsv"; p.write_text("signature\tgene\tweight\nMine\tcd8a\t1\nMine\tPRF1\t2\nMine\tGZMA\t\n")
    U = read_signature_table(p)
    assert list(U.gene) == ["CD8A", "PRF1", "GZMA"] and list(U.weight) == [1.0, 2.0, 1.0] and set(U.family) == {"User"}


@pytest.fixture(scope="module")
def small_benchmark():
    L = load_library(); keep = ["Tcell-inflamed_Ayers", "TGFB_Mariathasan", "CYT_Davoli", "IFN_Ayers"]
    return benchmark(demo_trials(), library=L[L.signature.isin(keep)], genes=("CD8A", "CXCL10"), n_random=19, n_boot=30, seed=3)


def test_benchmark_outputs(small_benchmark):
    B = small_benchmark
    assert set(B) == {"per_trial", "pooled", "paired", "summary", "meta", "settings"}
    P = B["pooled"].set_index("score")
    assert P.loc[MAIN_REFERENCE, "est"] > 0.1                             # the synthetic truth is +0.25 on the T-cell programme
    sig = P[P.kind == "signature"]
    assert sig.est.notna().sum() == 4 and sig.emp_p_declared.between(1 / 20, 1).all()
    assert B["settings"].min_trials.iloc[0] == 4 and B["settings"].boot_seed.iloc[0] == 1003
    PR = B["paired"]
    assert set(PR.reference) == {MAIN_REFERENCE, "T-cell mean (CD8A, PRF1, GZMA)", "CD8A", "CXCL10"}
    assert (PR.mdd > 0).all() and PR.q.notna().all()
    S = B["summary"].set_index("test")
    assert S.loc["favor the checkpoint arm", "of"] == 4


def test_benchmark_parallel_equals_serial():
    L = load_library(); L = L[L.signature.isin(["IFN_Ayers", "TGFB_Mariathasan"])]; T = demo_trials()[:3]
    a = benchmark(T, library=L, genes=("CD8A",), n_random=9, n_boot=12, min_trials=2, n_jobs=1)
    b = benchmark(T, library=L, genes=("CD8A",), n_random=9, n_boot=12, min_trials=2, n_jobs=2)
    pd.testing.assert_frame_equal(a["pooled"], b["pooled"]); pd.testing.assert_frame_equal(a["paired"], b["paired"])


def test_pool_dl_hk_hand_example():
    y, se = np.array([0.1, 0.3, 0.5]), np.array([0.1, 0.1, 0.1])
    r = pool_dl_hk(y, se)
    w0 = 1 / se ** 2; Q = np.sum(w0 * (y - 0.3) ** 2); t2 = (Q - 2) / (w0.sum() - np.sum(w0 ** 2) / w0.sum())
    assert r["est"] == pytest.approx(0.3) and r["tau2"] == pytest.approx(t2) and r["k"] == 3


def test_single_arm_and_concordance(small_benchmark):
    S = single_arm(demo_cohorts(), {**signature_weights(load_library()), **reference_weights()})
    assert S["per_cohort"].cohort.nunique() == 8 and len(S["pooled"]) > 40
    t = S["pooled"].set_index("score").loc[MAIN_REFERENCE]
    assert t.est > 0.2 and t.p < 0.05                                     # simulated log-odds +0.45 on the T-cell programme
    c = concordance(S["pooled"], small_benchmark["pooled"])
    assert c["n"] >= 4 and -1 <= c["rho"] <= 1


def test_published_example_reproduces_manuscript():
    T = example_table(); assert len(T) == 19 and (T.exclude == "yes").sum() == 1
    r = pool_published(T)
    assert r["primary"]["k"] == 18
    assert r["primary"]["est"] == pytest.approx(0.246745, abs=5e-5) and r["primary"]["p"] == pytest.approx(0.004436, abs=5e-5)
    assert r["per_sd"]["est"] == pytest.approx(0.145889, abs=5e-5)
    assert r["egger"]["p"] == pytest.approx(0.2623, abs=5e-4)


def test_rma_without_heterogeneity_is_fixed_effect():
    y, se = np.array([0.2, 0.2, 0.2]), np.array([0.1, 0.2, 0.3]); m = rma(y, se)
    assert m["tau2"] == pytest.approx(0, abs=1e-8) and m["est"] == pytest.approx(0.2)


def test_per_sd_divisor_median_split():
    assert float(per_sd_divisor(np.array([0.5]), np.array([np.nan]))[0]) == pytest.approx(1.5958, abs=1e-4)


def test_design_numbers():
    assert events_needed(0.10, 0.55) == pytest.approx(2825.6, abs=0.1)
    assert events_needed(0.10, 0.84) == pytest.approx(1004.7, abs=0.1)
    n = events_needed(0.15, 0.7)
    assert power(n, 0.15, 0.7) == pytest.approx(0.8, abs=1e-3)
    assert detectable_difference(n, 0.7) == pytest.approx(0.15)
    C = curve(deltas=(0.1, 0.2)); assert set(C.delta) == {0.1, 0.2} and (C.groupby("delta").events.apply(lambda s: s.is_monotonic_decreasing)).all()


def test_report_writes_html_and_pdf(tmp_path, small_benchmark):
    pytest.importorskip("matplotlib")
    from barriernet.report import write_report
    r = write_report(tmp_path, benchmark=small_benchmark, published=pool_published(example_table()), design=curve(), title="test")
    html = r["html"].read_text()
    assert "data:image/png;base64" in html and "Published subgroup results" in html and r["pdf"].stat().st_size > 10000
    assert (tmp_path / "figures" / "benchmark.png").exists() and (tmp_path / "figures" / "design.pdf").exists()


def test_bundle_rejects_paths_outside(tmp_path):
    from barriernet.bundle import extract
    z = tmp_path / "bad.zip"
    with zipfile.ZipFile(z, "w") as f:
        f.writestr("../evil.txt", "x")
    with pytest.raises(ValueError):
        extract(z, tmp_path / "out")


def test_bundle_loads_trials(tmp_path):
    from barriernet.bundle import load_trials_zip
    tr = demo_trials()[0]; Z = tr.Z.iloc[:60]
    buf = io.StringIO(); Z.T.to_csv(buf, sep="\t")
    clin = pd.DataFrame({"sample_id": Z.columns, "arm": tr.arm, "time": tr.time, "event": tr.event})
    yml = "trials:\n  - name: one\n    expr: data/expr.tsv\n    clinical: data/clin.tsv\n    arm_col: arm\n    time_col: time\n    event_col: event\n"
    z = tmp_path / "t.zip"
    with zipfile.ZipFile(z, "w") as f:
        f.writestr("trials.yaml", yml); f.writestr("data/expr.tsv", buf.getvalue()); f.writestr("data/clin.tsv", clin.to_csv(sep="\t", index=False))
    trials, cfg = load_trials_zip(z)
    assert len(trials) == 1 and trials[0].Z.shape == (60, len(tr.arm)) and trials[0].kind == "cox"


def test_cli_design_and_library(tmp_path, capsys):
    from barriernet.cli import main
    assert main(["design", "--delta", "0.1", "--rho", "0.55", "--no-report", "--out", str(tmp_path / "d")]) == 0
    assert "2825" in capsys.readouterr().out
    assert main(["library", "--out", str(tmp_path / "lib.tsv")]) == 0
    assert pd.read_csv(tmp_path / "lib.tsv", sep="\t").signature.nunique() == 51


def test_user_signature_name_clash():
    from barriernet.library import combine
    U = pd.DataFrame({"signature": ["IFN_Ayers"], "gene": ["CD8A"], "weight": [1.0]})
    with pytest.raises(ValueError):
        combine(load_library(), U)


def test_benchmark_with_nothing_poolable():
    L = load_library(); L = L[L.signature == "IFN_Ayers"]
    B = benchmark(demo_trials()[:2], library=L, genes=("CD8A",), n_random=0, n_boot=5, min_trials=3)
    assert B["pooled"].set_index("score").loc["IFN_Ayers", "est"] != B["pooled"].set_index("score").loc["IFN_Ayers", "est"]   # NaN


def test_job_runs_in_its_own_process_and_cancels():
    import time
    from barriernet.jobs import Job
    spec = dict(source="demo", keep=["IFN_Ayers", "TGFB_Mariathasan", "CYT_Davoli", "Tcell-inflamed_Ayers"], genes=["CD8A"],
                n_random=5, n_boot=20, seed=1, n_jobs=1)
    j = Job(spec)
    t0 = time.time()
    while not j.done() and time.time() - t0 < 600:
        time.sleep(0.5)
    B = j.result(); j.cleanup()
    assert B["pooled"].set_index("score").loc["IFN_Ayers", "est"] == B["pooled"].set_index("score").loc["IFN_Ayers", "est"]
    k = Job(dict(spec, n_random=500, n_boot=300)); time.sleep(3); k.cancel()
    assert k.done(); k.cleanup()


def test_precomputed_demonstration():
    from barriernet.demo import DEMO_SETTINGS, precomputed_benchmark
    B = precomputed_benchmark()
    st = B["settings"].iloc[0]
    assert int(st.n_random) == DEMO_SETTINGS["n_random"] and int(st.n_boot) == DEMO_SETTINGS["n_boot"]
    assert B["pooled"].query("kind == 'signature'").est.notna().sum() == 51 and len(B["summary"]) == 11
