"""Benchmark published signatures in randomized comparisons, as in the manuscript.

For each signature: the pooled score x arm interaction (conservative random-effects rule); a size-matched random-gene
reference in the signature's declared direction; paired differences from the T-cell references (the 13-gene T-cell
mean first, then the CD8A-PRF1-GZMA mean) and from single genes, with the difference detectable at 80% power; and
Benjamini-Hochberg correction within each family of tests over the signatures. A signature is pooled when it can be
scored in at least ``min_trials`` comparisons (default: 60% of them, at least two; four of six in the manuscript)."""
from __future__ import annotations
import math
try:
    import multiprocessing as mp
except ImportError:                  # e.g. Python in a web browser (Pyodide)
    mp = None
import os
from typing import Callable, Dict, Optional, Sequence
import numpy as np
import pandas as pd
from .calibrate import Trial, bootstrap_indices, bootstrap_interactions, interaction, paired_difference, pool, random_reference
from .library import MAIN_REFERENCE, PANEL_GENES, load_library, reference_weights, score_signature, signature_meta, signature_weights


_STATE: dict = {}


def _boot_task(key):
    S, TR, idx = _STATE["S"], _STATE["TR"], _STATE["idx"]
    return bootstrap_interactions(TR[key[1]], S[key], idx[key[1]])


def _ref_task(args):
    ks, n_random, seed = args[1], _STATE["n_random"], _STATE["seed"]
    TR, USED = _STATE["TR"], _STATE["USED"]
    return random_reference([TR[t] for t in ks], {t: USED[(args[0], t)] for t in ks}, n_random, seed)


def _map(fn, items, n_jobs, say, label):
    items = list(items); out = []
    if n_jobs > 1 and len(items) > 1 and mp is not None and "fork" in mp.get_all_start_methods():
        with mp.get_context("fork").Pool(n_jobs) as P:
            for i, r in enumerate(P.imap(fn, items, chunksize=1)):
                out.append(r); say(label, i + 1, len(items))
    else:
        for i, it in enumerate(items):
            out.append(fn(it)); say(label, i + 1, len(items))
    return out


def bh(p) -> np.ndarray:
    """Benjamini-Hochberg q values (NaN kept)."""
    p = np.asarray(p, float); q = np.full(p.shape, np.nan); ok = np.isfinite(p); m = int(ok.sum())
    if m:
        o = np.argsort(p[ok]); r = p[ok][o] * m / np.arange(1, m + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]; out = np.empty(m); out[o] = np.minimum(r, 1.0); q[ok] = out
    return q


def benchmark(trials: Sequence[Trial], library: Optional[pd.DataFrame] = None, references: Optional[Dict[str, list]] = None,
              genes: Sequence[str] = PANEL_GENES, extra_scores: Optional[Dict[str, dict]] = None, n_random: int = 1000,
              n_boot: int = 300, seed: int = 0, min_trials: Optional[int] = None, n_jobs: int = 1, boot_seed: Optional[int] = None, fdr_genes: Sequence[str] = ("CD8A", "CXCL10"),
              progress: Optional[Callable[[str, int, int], None]] = None) -> Dict[str, pd.DataFrame]:
    """``extra_scores``: name -> dict(weights={gene: w}, direction=+1/-1, family=str) benchmarked with the library.
    ``boot_seed``: seed of the patient resamples, drawn as default_rng([boot_seed, trial index]) (default seed + 1000;
    the manuscript used 3032). ``n_jobs``: worker processes for the bootstrap and the random reference (fork start method; -1 = all cores).
    Returns per_trial, pooled, paired, summary, meta and settings."""
    say = progress or (lambda *a: None)
    L = load_library() if library is None else library
    meta = signature_meta(L); W = signature_weights(L)
    for nm, v in (extra_scores or {}).items():
        W[nm] = {str(g).upper(): float(w) for g, w in v["weights"].items()}
        meta.loc[len(meta)] = dict(signature=nm, family=v.get("family", "User"), direction=int(v.get("direction", 1)),
                                   published_engine="", pmid="", n_genes=len(v["weights"]))
    REF = reference_weights(references)
    k = len(trials); min_trials = min_trials or max(2, math.ceil(0.6 * k))
    TR = {t.name: t for t in trials}
    S, USED = {}, {}
    for nm, w in list(W.items()) + list(REF.items()):
        for t in trials:
            sc, _, wu = score_signature(t.Z, w)
            if sc is not None:
                S[(nm, t.name)], USED[(nm, t.name)] = sc, wu
    for g in genes:
        for t in trials:
            sc, _, wu = score_signature(t.Z, {g: 1.0}, min_genes=1)
            if sc is not None:
                S[(g, t.name)] = sc
    R = {}
    for i, (key, sc) in enumerate(S.items()):
        try:
            R[key] = interaction(TR[key[1]], sc)
        except Exception:
            R[key] = (np.nan, np.nan)
        say("interactions", i + 1, len(S))
    n_jobs = (os.cpu_count() or 1) if n_jobs == -1 else max(1, int(n_jobs))
    bs = seed + 1000 if boot_seed is None else boot_seed
    idx = {t.name: bootstrap_indices(t, n_boot, [bs, i]) for i, t in enumerate(trials)}
    _STATE.update(S=S, TR=TR, idx=idx, USED=USED, n_random=n_random, seed=seed)
    keys = [key for key in S if np.isfinite(R[key][0])]
    BT = dict(zip(keys, _map(_boot_task, keys, n_jobs, say, "bootstrap"))) if n_boot > 1 else {}
    per = pd.DataFrame([dict(score=a, trial=b, est=v[0], se=v[1]) for (a, b), v in R.items()])

    kinds = {**{s: "signature" for s in W}, **{r: "T-cell reference" for r in REF}, **{g: "single gene" for g in genes}}
    M = meta.set_index("signature")
    sigs = list(W)
    need_ref = []
    for nm in sigs:
        ks = [t.name for t in trials if (nm, t.name) in R and np.isfinite(R[(nm, t.name)][0])]
        if len(ks) >= min_trials and n_random > 0:
            need_ref.append((nm, ks))
    REFD = dict(zip([a for a, _ in need_ref], _map(_ref_task, need_ref, n_jobs, say, "random reference")))
    rows = []
    for nm in sigs + list(REF) + list(genes):
        ks = [t.name for t in trials if (nm, t.name) in R and np.isfinite(R[(nm, t.name)][0])]
        need = min_trials if kinds[nm] == "signature" else 2
        row = dict(score=nm, kind=kinds[nm], family=M.loc[nm, "family"] if nm in M.index else kinds[nm],
                   direction=int(M.loc[nm, "direction"]) if nm in M.index else 1, trials_scored=len(ks))
        if len(ks) >= need:
            row.update(pool([R[(nm, t)][0] for t in ks], [R[(nm, t)][1] for t in ks]))
            if nm in REFD:
                ref = REFD[nm]; d = row["direction"]; ok = np.isfinite(ref)
                row.update(ref_median=float(np.nanmedian(ref)), ref_p95_declared=float(np.nanpercentile(d * ref, 95) * d),
                           emp_p_declared=float((1 + np.sum(d * ref[ok] >= d * row["est"])) / (ok.sum() + 1)))
        rows.append(row)
    P = pd.DataFrame(rows)
    for c in ("est", "lo", "hi", "p", "k"):
        if c not in P:
            P[c] = np.nan
    sel = (P["kind"] == "signature") & P["est"].notna()
    P["q"] = np.nan; P.loc[sel, "q"] = bh(P.loc[sel, "p"])
    if "emp_p_declared" in P:
        P["q_reference"] = np.nan; P.loc[sel, "q_reference"] = bh(P.loc[sel, "emp_p_declared"])

    refs_for_pairs = list(REF) + list(genes)
    prs = []
    for nm in [s for s in sigs if s in set(P.loc[sel, "score"])]:
        for rf in refs_for_pairs:
            ks = [t.name for t in trials if (nm, t.name) in BT and (rf, t.name) in BT]
            if len(ks) < 2:
                continue
            q = paired_difference([R[(nm, t)][0] - R[(rf, t)][0] for t in ks], [BT[(nm, t)] for t in ks], [BT[(rf, t)] for t in ks])
            prs.append(dict(score=nm, reference=rf, reference_kind=kinds[rf], **q))
    PR = pd.DataFrame(prs)
    if len(PR):
        PR["q"] = np.nan
        for rf in list(REF) + [g for g in fdr_genes if g in genes]:
            m = PR["reference"] == rf; PR.loc[m, "q"] = bh(PR.loc[m, "p"])

    def count(df, cond, qcol="q"):
        return int(cond.sum()), (int((cond & (df[qcol] < 0.05)).sum()) if qcol in df else None)
    S_ = []
    Ps = P[sel]
    for lab, cond in (("favor the checkpoint arm", (Ps["est"] > 0) & (Ps["p"] < 0.05)), ("favor the control arm", (Ps["est"] < 0) & (Ps["p"] < 0.05))):
        n, nq = count(Ps, cond); S_.append(dict(test=lab, nominal=n, fdr=nq, of=len(Ps)))
    if "emp_p_declared" in Ps:
        n, nq = count(Ps, Ps["emp_p_declared"] < 0.05, "q_reference"); S_.append(dict(test="beyond size-matched random genes", nominal=n, fdr=nq, of=len(Ps)))
    for rf in list(REF) + [g for g in fdr_genes if g in genes]:
        D = PR[PR["reference"] == rf] if len(PR) else PR
        if not len(D):
            continue
        for lab, cond in ((f"above {rf}", (D["diff"] > 0) & (D["p"] < 0.05)), (f"below {rf}", (D["diff"] < 0) & (D["p"] < 0.05))):
            n, nq = count(D, cond); S_.append(dict(test=lab, nominal=n, fdr=nq, of=len(D)))
    return dict(per_trial=per, pooled=P, paired=PR, summary=pd.DataFrame(S_), meta=meta,
                settings=pd.DataFrame([dict(trials=k, min_trials=min_trials, n_random=n_random, n_boot=n_boot, seed=seed, boot_seed=bs, main_reference=MAIN_REFERENCE)]))
