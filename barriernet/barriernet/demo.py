"""Synthetic data for trying the package. Nothing here is patient data: expression is simulated from a few latent
cell programmes, and outcomes are simulated with a known score x arm interaction on the T-cell programme, so the
benchmark has a true answer (signatures rich in T-cell genes carry the interaction; random genes do not)."""
from __future__ import annotations
from typing import List
import numpy as np
import pandas as pd
from .calibrate import Trial
from .library import IFN_CHECKPOINT_PANEL, T_CELL_PANEL, load_library
from .singlearm import Cohort

TRUE_INTERACTION = 0.25          # benefit-oriented, per SD of the T-cell programme


def _genes_and_loadings(seed: int = 0, n_filler: int = 600):
    rng = np.random.default_rng(seed); L = load_library(distinct=False)
    fam = L.groupby("gene")["family"].agg(lambda s: "Inflamed/IFN" if (s == "Inflamed/IFN").any() else s.iloc[0])
    genes = sorted(set(L.gene) | set(T_CELL_PANEL) | set(IFN_CHECKPOINT_PANEL)) + [f"FILLER{i:04d}" for i in range(n_filler)]
    lo = {}
    for g in genes:
        if g in T_CELL_PANEL:
            t = rng.uniform(0.65, 0.85)
        elif g in IFN_CHECKPOINT_PANEL:
            t = rng.uniform(0.45, 0.70)
        elif fam.get(g) == "Inflamed/IFN":
            t = rng.uniform(0.15, 0.55)
        else:
            t = rng.uniform(0.0, 0.25)
        s = rng.uniform(0.0, 0.6) if fam.get(g) in ("Other", "Confounder") else rng.uniform(0.0, 0.2)
        lo[g] = (t, min(s, np.sqrt(max(0.0, 0.95 - t * t))))
    return genes, lo


def _expression(n, genes, lo, rng):
    T = rng.standard_normal(n); S = 0.3 * T + np.sqrt(1 - 0.09) * rng.standard_normal(n)
    X = np.empty((len(genes), n))
    for i, g in enumerate(genes):
        t, s = lo[g]; X[i] = t * T + s * S + np.sqrt(max(1e-6, 1 - t * t - s * s)) * rng.standard_normal(n)
    Z = pd.DataFrame(X, index=genes, columns=[f"S{j:04d}" for j in range(n)])
    Z = Z.sub(Z.mean(1), axis=0).div(Z.std(1), axis=0)
    return Z, T


def demo_trials(seed: int = 1, sizes=(260, 240, 300, 220, 360, 200)) -> List[Trial]:
    """Six synthetic randomized comparisons: four with survival, two with response."""
    genes, lo = _genes_and_loadings(seed); rng = np.random.default_rng(seed + 100); out = []
    for i, n in enumerate(sizes):
        Z, T = _expression(n, genes, lo, rng); arm = rng.permutation(np.r_[np.ones(n // 2), np.zeros(n - n // 2)])
        Z.columns = [f"T{i + 1}_{c}" for c in Z.columns]
        if i < 4:
            lp = -0.25 * arm - 0.05 * T - TRUE_INTERACTION * T * arm
            t = rng.exponential(1 / (0.08 * np.exp(lp))); c = rng.uniform(6, 36, n)
            out.append(Trial(f"Synthetic trial {i + 1}", Z, arm, time=np.minimum(t, c), event=(t <= c).astype(int)))
        else:
            lp = -0.8 + 0.5 * arm + 0.05 * T + TRUE_INTERACTION * T * arm
            out.append(Trial(f"Synthetic trial {i + 1}", Z, arm, response=(rng.random(n) < 1 / (1 + np.exp(-lp))).astype(float)))
    return out


def demo_cohorts(seed: int = 2, sizes=(48, 60, 72, 40, 110, 56, 90, 64)) -> List[Cohort]:
    """Eight synthetic single-arm checkpoint cohorts with response."""
    genes, lo = _genes_and_loadings(seed - 1); rng = np.random.default_rng(seed + 200); out = []
    for i, n in enumerate(sizes):
        Z, T = _expression(n, genes, lo, rng); Z.columns = [f"C{i + 1}_{c}" for c in Z.columns]
        lp = -0.9 + 0.45 * T
        out.append(Cohort(f"Synthetic cohort {i + 1}", Z, (rng.random(n) < 1 / (1 + np.exp(-lp))).astype(float)))
    return out


DEMO_SETTINGS = dict(n_random=1000, n_boot=300, seed=0)


def precomputed_benchmark() -> dict:
    """The benchmark of the synthetic trials computed in advance with DEMO_SETTINGS (data/demo_benchmark/), so that a
    page can show it at once; ``write_demo_benchmark`` recomputes it."""
    from importlib import resources
    out = {}
    base = resources.files("barriernet").joinpath("data/demo_benchmark")
    for k in ("per_trial", "pooled", "paired", "summary", "meta", "settings"):
        with base.joinpath(f"benchmark_{k}.tsv").open("r", encoding="utf-8") as fh:
            out[k] = pd.read_csv(fh, sep="\t", dtype={"pmid": str} if k == "meta" else None, keep_default_na=True)
    out["summary"]["fdr"] = out["summary"]["fdr"].astype(float)
    return out


def write_demo_benchmark(folder, n_jobs: int = 1) -> None:
    """Recompute the precomputed demonstration (used when the package is built)."""
    import pathlib
    from .benchmark import benchmark
    folder = pathlib.Path(folder); folder.mkdir(parents=True, exist_ok=True)
    B = benchmark(demo_trials(), n_jobs=n_jobs, **DEMO_SETTINGS)
    for k, df in B.items():
        df.to_csv(folder / f"benchmark_{k}.tsv", sep="\t", index=False)
