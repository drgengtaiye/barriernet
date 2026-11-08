from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

from .genesets import load_yaml, default_yaml, load_weights
from .io import load_expr_matrix, read_table, looks_linear
from .core import score_modules, score_axes, summarise_used_genes, score_weighted
from .spatial import infer_spatial_keys, assign_rings, ring_summary, barrier_index

def _ensure_out(out: str) -> Path:
    p = Path(out)
    p.mkdir(parents=True, exist_ok=True)
    return p

def cmd_bulk(args) -> int:
    out = _ensure_out(args.out)
    gs = load_yaml(args.geneset or default_yaml("rna"))

    expr, _ = load_expr_matrix(args.expr, gene_col=args.gene_col, orientation=args.orientation)
    if args.log2p1:
        expr = np.log2(expr.clip(lower=0) + 1)
    elif looks_linear(expr):
        print("[barriernet bulk] WARNING: values look un-logged (max > 100). Z-scoring linear TPM/counts "
              "lets a few highly expressed samples dominate each gene; consider --log2p1.")
    scores, used = score_modules(expr, gs.modules, min_genes=args.min_genes, zscore=not args.no_zscore)
    axes = score_axes(scores, gs.axes)

    out_df = pd.concat([scores, axes], axis=1)
    if not args.no_v03:
        wt = load_weights(args.weights)
        out_df["barriernet03"], _ = score_weighted(expr, wt["weights"], wt["min_genes"], wt["min_frac"], zscore=not args.no_zscore)
    out_df.index.name = "sample_id"
    out_df.to_csv(out/"barriernet_bulk_scores.tsv", sep="\t")

    summarise_used_genes(used).to_csv(out/"barriernet_bulk_used_genes.tsv", sep="\t", index=False)
    panel_genes = sorted({g for gl in used.values() for g in gl})
    if panel_genes:
        expr.loc[:, [g for g in panel_genes if g in expr.columns]].to_csv(out/"barriernet_bulk_panel_matrix.tsv", sep="\t")

    print(f"[barriernet bulk] wrote outputs to: {out}")
    return 0

def _load_trials(cfg):
    """Trials from the YAML layout of 'calibrate' (see README)."""
    from .calibrate import Trial
    trials = []
    for t in cfg["trials"]:
        expr, _ = load_expr_matrix(t["expr"], gene_col=t.get("gene_col"), orientation=t.get("orientation", "auto"))
        if t.get("log2p1"):
            expr = np.log2(expr.clip(lower=0) + 1)
        clin = read_table(t["clinical"]).set_index(t.get("sample_col", "sample_id"))
        ids = [i for i in expr.index if i in clin.index]
        if len(ids) < 10:
            raise SystemExit(f"{t['name']}: fewer than 10 samples shared by expression and clinical tables")
        X = expr.loc[ids].T; X = X[X.std(axis=1) > 0]
        Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1), axis=0)
        c = clin.loc[ids]; cov = t.get("covariates") or []
        trials.append(Trial(t["name"], Z, c[t["arm_col"]].to_numpy(float),
                            time=c[t["time_col"]].to_numpy(float) if t.get("time_col") else None,
                            event=c[t["event_col"]].to_numpy(int) if t.get("event_col") else None,
                            response=c[t["response_col"]].to_numpy(float) if t.get("response_col") else None,
                            covariates=c[cov].to_numpy(float) if cov else None))
    return trials


def _load_cohorts(cfg):
    """Single-arm cohorts: name, expr, clinical, sample_col, response_col (1 = responder), optional covariates."""
    from .singlearm import Cohort
    out = []
    for t in cfg["cohorts"]:
        expr, _ = load_expr_matrix(t["expr"], gene_col=t.get("gene_col"), orientation=t.get("orientation", "auto"))
        if t.get("log2p1"):
            expr = np.log2(expr.clip(lower=0) + 1)
        clin = read_table(t["clinical"]).set_index(t.get("sample_col", "sample_id"))
        ids = [i for i in expr.index if i in clin.index and pd.notna(clin.loc[i, t["response_col"]])]
        X = expr.loc[ids].T; X = X[X.std(axis=1) > 0]; Z = X.sub(X.mean(axis=1), axis=0).div(X.std(axis=1), axis=0)
        cov = t.get("covariates") or []
        out.append(Cohort(t["name"], Z, clin.loc[ids, t["response_col"]].to_numpy(float), clin.loc[ids, cov].to_numpy(float) if cov else None))
    return out


def _save(res: dict, out: Path, prefix: str) -> None:
    for k, df in res.items():
        if isinstance(df, pd.DataFrame):
            df.to_csv(out / f"{prefix}_{k}.tsv", sep="\t", index=False)


def cmd_calibrate(args) -> int:
    """Calibrate scores over randomized trials described in a YAML file (see README, 'calibrate')."""
    import yaml
    from .calibrate import calibrate
    out = _ensure_out(args.out)
    cfg = yaml.safe_load(open(args.config))
    trials = _load_trials(cfg)
    gs = load_yaml(default_yaml("rna")); wt = load_weights(cfg.get("weights"))
    scores = {"barriernet03": wt["weights"]}
    net = {g: 0.25 / len(gs.modules[m]) for m in ("APM", "IFN", "TLS_B", "T_CELL") for g in gs.modules[m]}
    net.update({g: -0.5 / len(gs.modules[m]) for m in ("CAF_ECM", "TAM_SUP") for g in gs.modules[m]})
    scores["NetImmune"] = net
    for nm, w in (cfg.get("scores") or {}).items():
        scores[nm] = {str(k).upper(): float(v) for k, v in w.items()}
    from .library import PANEL_GENES
    genes = cfg.get("genes") or list(PANEL_GENES)
    res = calibrate(trials, scores, genes, n_random=int(cfg.get("n_random", 1000)), n_boot=int(cfg.get("n_boot", 300)),
                    seed=int(cfg.get("seed", 0)))
    for k, df in res.items():
        df.to_csv(out / f"calibrate_{k}.tsv", sep="\t", index=False)
    print(res["pooled"].round(3).to_string(index=False))
    print(f"[barriernet calibrate] wrote outputs to: {out}")
    return 0


# ------------------------------------------------------------------ 0.4: the manuscript's analyses
def _progress(stage, i, n):
    if n and (i == n or i % max(1, n // 10) == 0):
        print(f"  {stage}: {i}/{n}", flush=True)


def cmd_library(args) -> int:
    from .library import load_library
    L = load_library(distinct=not args.all); L.to_csv(args.out, sep="\t", index=False)
    print(f"[barriernet library] {L.signature.nunique()} signatures, {len(L)} gene rows -> {args.out}"); return 0


def cmd_benchmark(args) -> int:
    import yaml
    from .benchmark import benchmark
    from .library import combine, load_library, read_signature_table
    out = _ensure_out(args.out); cfg = yaml.safe_load(open(args.config)); trials = _load_trials(cfg)
    lib = load_library()
    if args.signatures:
        U = read_signature_table(args.signatures); lib = U if args.only_user else combine(lib, U)
    B = benchmark(trials, library=lib, n_random=args.n_random, n_boot=args.n_boot, seed=args.seed, boot_seed=args.boot_seed,
                  min_trials=args.min_trials, n_jobs=args.jobs, progress=_progress)
    _save(B, out, "benchmark"); print(B["summary"].to_string(index=False))
    if not args.no_report:
        from .report import write_report
        r = write_report(out, benchmark=B, title="Benchmark of published signatures"); print(f"[barriernet benchmark] report: {r['html']}")
    print(f"[barriernet benchmark] wrote outputs to: {out}"); return 0


def cmd_singlearm(args) -> int:
    import yaml
    from .library import load_library, read_signature_table, reference_weights, signature_weights
    from .singlearm import concordance, single_arm
    out = _ensure_out(args.out); cfg = yaml.safe_load(open(args.config)); cohorts = _load_cohorts(cfg)
    lib = read_signature_table(args.signatures) if args.signatures else load_library()
    scores = {**signature_weights(lib), **reference_weights()}
    S = single_arm(cohorts, scores, min_cohorts=args.min_cohorts); _save(S, out, "singlearm")
    conc = None
    if args.randomized:
        conc = concordance(S["pooled"], pd.read_csv(args.randomized, sep=None, engine="python"))
        conc["table"].to_csv(out / "singlearm_concordance.tsv", sep="\t", index=False)
        print(f"Spearman with the randomized interaction: {conc['rho']:.3f} over {conc['n']} scores")
    print(S["pooled"].sort_values("est", ascending=False).head(20).round(3).to_string(index=False))
    if not args.no_report:
        from .report import write_report
        write_report(out, single_arm=S, concordance=conc, title="Single-arm cohorts")
    print(f"[barriernet singlearm] wrote outputs to: {out}"); return 0


def _published_rows(res):
    return pd.DataFrame([dict(analysis=k, **{c: v.get(c) for c in ("k", "est", "lo", "hi", "p", "tau2", "I2", "pi_lo", "pi_hi", "t", "trials")})
                         for k, v in res.items() if isinstance(v, dict)])


def cmd_published(args) -> int:
    from .published import example_table, pool_published
    out = _ensure_out(args.out)
    T = example_table() if args.example else pd.read_csv(args.table, sep=None, engine="python")
    res = pool_published(T); A = _published_rows(res); A.to_csv(out / "published_analyses.tsv", sep="\t", index=False)
    res["per_trial"].to_csv(out / "published_per_trial.tsv", sep="\t", index=False)
    print(A.drop(columns=["trials"]).round(4).to_string(index=False))
    if not args.no_report:
        from .report import write_report
        write_report(out, published=res, title="Published randomized trials")
    print(f"[barriernet published] wrote outputs to: {out}"); return 0


def cmd_design(args) -> int:
    from .design import curve, detectable_difference, events_needed, power
    rows = []
    for r in args.rho:
        for d in args.delta:
            row = dict(rho=r, delta=d, events_needed=events_needed(d, r, args.allocation, args.power, args.alpha, args.inflation))
            if args.events:
                row.update(events=args.events, power_at_events=power(args.events, d, r, args.allocation, args.alpha, args.inflation),
                           detectable_at_events=detectable_difference(args.events, r, args.allocation, args.power, args.alpha, args.inflation))
            rows.append(row)
    T = pd.DataFrame(rows); print(T.round(3).to_string(index=False))
    if args.out:
        out = _ensure_out(args.out); T.to_csv(out / "design_table.tsv", sep="\t", index=False)
        C = curve(deltas=tuple(args.delta), allocation=args.allocation, power_=args.power, alpha=args.alpha, inflation=args.inflation)
        C.to_csv(out / "design_curve.tsv", sep="\t", index=False)
        if not args.no_report:
            from .report import write_report
            write_report(out, design=C, design_marks={f"r = {r:g}": r for r in args.rho}, title="Trial design")
        print(f"[barriernet design] wrote outputs to: {out}")
    return 0


def cmd_report(args) -> int:
    """Rebuild report.html and report.pdf from the tables written by the other commands in one folder."""
    from .report import write_report
    d = Path(args.dir)

    def rd(n):
        return pd.read_csv(d / n, sep="\t") if (d / n).exists() else None
    B = None
    if (d / "benchmark_pooled.tsv").exists():
        B = {k: rd(f"benchmark_{k}.tsv") for k in ("per_trial", "pooled", "paired", "summary", "meta", "settings")}
        B["paired"] = B["paired"] if B["paired"] is not None else pd.DataFrame()
    S = {k: rd(f"singlearm_{k}.tsv") for k in ("per_cohort", "pooled")} if (d / "singlearm_pooled.tsv").exists() else None
    conc = None
    if (d / "singlearm_concordance.tsv").exists():
        from scipy.stats import spearmanr
        T = rd("singlearm_concordance.tsv"); r = spearmanr(T.interaction, T.single_arm)
        conc = dict(rho=float(r[0]), p=float(r[1]), n=len(T), table=T)
    pub = None
    if (d / "published_per_trial.tsv").exists():
        I = rd("published_per_trial.tsv"); A = rd("published_analyses.tsv")
        pub = {r.analysis: r.dropna().to_dict() for _, r in A.iterrows()}; pub["per_trial"] = I
    C = rd("design_curve.tsv")
    r = write_report(d, benchmark=B, single_arm=S, concordance=conc, published=pub, design=C, title=args.title)
    print(f"[barriernet report] {r['html']} and {r['pdf']}"); return 0


def cmd_demo(args) -> int:
    """Every analysis on synthetic data (not patient data), with the report."""
    from .benchmark import benchmark
    from .demo import demo_cohorts, demo_trials
    from .design import curve
    from .library import load_library, reference_weights, signature_weights
    from .published import example_table, pool_published
    from .report import write_report
    from .singlearm import concordance, single_arm
    out = _ensure_out(args.out); n_r, n_b = (100, 100) if args.quick else (1000, 300)
    print("synthetic randomized trials ...")
    B = benchmark(demo_trials(), n_random=n_r, n_boot=n_b, n_jobs=args.jobs, progress=_progress); _save(B, out, "benchmark")
    print("synthetic single-arm cohorts ...")
    S = single_arm(demo_cohorts(), {**signature_weights(load_library()), **reference_weights()}); _save(S, out, "singlearm")
    conc = concordance(S["pooled"], B["pooled"]); conc["table"].to_csv(out / "singlearm_concordance.tsv", sep="\t", index=False)
    pub = pool_published(example_table()); _published_rows(pub).to_csv(out / "published_analyses.tsv", sep="\t", index=False)
    pub["per_trial"].to_csv(out / "published_per_trial.tsv", sep="\t", index=False)
    C = curve(); C.to_csv(out / "design_curve.tsv", sep="\t", index=False)
    r = write_report(out, benchmark=B, single_arm=S, concordance=conc, published=pub, design=C,
                     design_marks={"r = 0.55": 0.55, "r = 0.81": 0.81, "r = 0.84": 0.84}, title="barriernet demonstration",
                     note="The randomized trials and single-arm cohorts here are synthetic, simulated with a known interaction on a T-cell "
                          "programme; they are not patient data. The published-trials section uses subgroup results extracted from "
                          "published reports; check them against the sources before relying on them.")
    print(B["summary"].to_string(index=False)); print(f"[barriernet demo] report: {r['html']}"); return 0


def cmd_app(args) -> int:
    try:
        import shiny  # noqa: F401
    except ImportError:
        raise SystemExit("the web app needs Shiny for Python: pip install 'barriernet[app]'")
    import uvicorn
    from .app import app
    print(f"barriernet app on http://{args.host}:{args.port}")
    uvicorn.run(app, host=args.host, port=args.port, log_level="warning"); return 0


def cmd_ihc(args) -> int:
    out = _ensure_out(args.out)
    gs = load_yaml(args.geneset or default_yaml("ihc"))

    df = read_table(args.table)
    if args.sample_col not in df.columns:
        raise SystemExit(f"sample_col '{args.sample_col}' not found in columns")
    sample_col = args.sample_col

    drop_cols = set([sample_col] + (args.drop_cols or []))
    marker_cols = [c for c in df.columns if c not in drop_cols]

    X = df[marker_cols].copy()
    X.columns = [str(c).upper() for c in X.columns]
    for c in X.columns:
        X[c] = pd.to_numeric(X[c], errors="coerce")

    if args.agg == "median":
        Xg = pd.concat([df[[sample_col]], X], axis=1).groupby(sample_col).median(numeric_only=True)
    else:
        Xg = pd.concat([df[[sample_col]], X], axis=1).groupby(sample_col).mean(numeric_only=True)

    Xg.index = Xg.index.astype(str)
    scores, used = score_modules(Xg, gs.modules, min_genes=args.min_genes, zscore=not args.no_zscore)
    axes = score_axes(scores, gs.axes)

    out_df = pd.concat([scores, axes], axis=1)
    out_df.index.name = "sample_id"
    out_df.to_csv(out/"barriernet_ihc_scores.tsv", sep="\t")
    summarise_used_genes(used).to_csv(out/"barriernet_ihc_used_markers.tsv", sep="\t", index=False)

    print(f"[barriernet ihc] wrote outputs to: {out}")
    return 0

def _read_h5ad_any(path: str):
    try:
        import anndata as ad
        return ad.read_h5ad(path)
    except Exception:
        import scanpy as sc
        return sc.read_h5ad(path)

def cmd_spatial(args) -> int:
    out = _ensure_out(args.out)
    paths = []
    if args.h5ad:
        paths = [args.h5ad]
    else:
        import glob
        paths = sorted(glob.glob(args.h5ad_glob))

    if not paths:
        raise SystemExit("No h5ad files found.")

    sections = []
    for p in paths:
        adata = _read_h5ad_any(p)
        if "spatial" not in adata.obsm:
            raise SystemExit(f"{p}: missing obsm['spatial'].")

        xy = np.asarray(adata.obsm["spatial"]).astype(float)
        keys = infer_spatial_keys(adata.obs)

        ring_df = assign_rings(adata.obs, xy, keys.tls, peak_q=args.peak_q, peak_k=args.peak_k,
                               far_min_dist=args.far_min_dist, far_exclude_edge=args.far_exclude_edge,
                               edge_buffer=args.edge_buffer, peak_method=args.peak_method, coords=args.coords)

        stem = Path(p).name.replace(".h5ad","")
        spot = pd.concat([ring_df], axis=1)
        # Keep a minimal per-spot export that matches your Figure3 expectations
        spot.to_csv(out/f"{stem}_spot_rings.tsv", sep="\t")

        rs = ring_summary(adata.obs, ring_df, keys)
        rs.to_csv(out/f"{stem}_ring_enrichment.tsv", sep="\t", index=False)
        sections.append(dict(section=stem, **{c: rs[c].iloc[0] for c in rs.columns if "_minus_Far" in c}))

        # Write back to h5ad
        for c in ring_df.columns:
            adata.obs[c] = ring_df[c].values
        out_h5 = out/f"{stem}.with_barriernet.h5ad"
        try:
            adata.write(out_h5)
        except Exception:
            import scanpy as sc
            sc.write(out_h5, adata)

        print(f"[barriernet spatial] processed: {stem}")

    if len(sections) >= 3:
        tab = pd.DataFrame(sections).set_index("section")
        tab["BarrierIndex"] = barrier_index(tab)
        tab.to_csv(out/"barriernet_sections.tsv", sep="\t")
    print(f"[barriernet spatial] wrote outputs to: {out}")
    return 0

def _spatial_args(p):
    p.add_argument("--peak-q", type=float, default=0.95, help="TLS peak quantile (default 0.95).")
    p.add_argument("--peak-k", type=int, default=9, help="Local-maximum window: k nearest spots incl. itself (default 9).")
    p.add_argument("--far-min-dist", type=int, default=3, help="Far = at least this many hops from a TLS peak (default 3; 0.2.0 used 4).")
    p.add_argument("--far-exclude-edge", action="store_true", help="Exclude tissue-edge spots (<6 in-tissue neighbours) from Far.")
    p.add_argument("--edge-buffer", type=int, default=0, help="With --far-exclude-edge, also exclude spots within this many hops of the edge.")
    p.add_argument("--peak-method", default="knn", choices=["knn", "lattice"], help="Local maxima among k nearest spots (default) or among lattice neighbours.")
    p.add_argument("--coords", default="visium", choices=["visium", "axial"], help="Coordinate system of array_row/array_col.")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="barriernet",
        description="Benchmarking immune signatures in randomized checkpoint trials; Barrier/NetImmune scoring (bulk, IHC/IF) and Visium ring/peak summaries."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    p_bulk = sub.add_parser("bulk", help="Compute Barrier/NetImmune axes from bulk RNA/proteomics expression matrix.")
    p_bulk.add_argument("--expr", required=True, help="Expression matrix TSV/CSV (samples×genes or genes×samples).")
    p_bulk.add_argument("--gene-col", default=None, help="Optional column name containing gene symbols.")
    p_bulk.add_argument("--orientation", default="auto", choices=["auto","samples_by_genes","genes_by_samples"])
    p_bulk.add_argument("--geneset", default=None, help="Path to geneset YAML (default: bundled RNA geneset).")
    p_bulk.add_argument("--min-genes", type=int, default=3)
    p_bulk.add_argument("--no-zscore", action="store_true")
    p_bulk.add_argument("--log2p1", action="store_true", help="Apply log2(x + 1) before z-scoring (use for linear TPM/counts).")
    p_bulk.add_argument("--weights", default=None, help="Weighted-score YAML (default: bundled barriernet 0.3 weights).")
    p_bulk.add_argument("--no-v03", action="store_true", help="Do not add the barriernet03 weighted score column.")
    p_bulk.add_argument("--out", required=True)
    p_bulk.set_defaults(func=cmd_bulk)

    p_cal = sub.add_parser("calibrate", help="Calibrate scores against treatment interaction, random genes and single genes "
                                             "over randomized trials (YAML config).")
    p_cal.add_argument("--config", required=True, help="YAML listing the trials (expression, clinical table, columns).")
    p_cal.add_argument("--out", required=True)
    p_cal.set_defaults(func=cmd_calibrate)

    p_lib = sub.add_parser("library", help="Export the bundled library of published immune signatures (genes, weights, family, direction).")
    p_lib.add_argument("--out", required=True)
    p_lib.add_argument("--all", action="store_true", help="Keep TIS_Damotte, a duplicate of Tcell-inflamed_Ayers.")
    p_lib.set_defaults(func=cmd_library)

    p_bm = sub.add_parser("benchmark", help="Benchmark the published signatures (and your own) in randomized trials against "
                                            "T-cell references, random genes and single genes.")
    p_bm.add_argument("--config", required=True, help="YAML listing the trials (same layout as 'calibrate').")
    p_bm.add_argument("--signatures", help="Your signatures as a long table: signature, gene, optional weight, direction, family.")
    p_bm.add_argument("--only-user", action="store_true", help="Benchmark only the signatures in --signatures.")
    p_bm.add_argument("--n-random", type=int, default=1000)
    p_bm.add_argument("--n-boot", type=int, default=300)
    p_bm.add_argument("--seed", type=int, default=0)
    p_bm.add_argument("--boot-seed", type=int, default=None)
    p_bm.add_argument("--min-trials", type=int, default=None, help="Pool a signature scored in at least this many trials (default 60%%, at least 2).")
    p_bm.add_argument("--jobs", type=int, default=1, help="Worker processes (-1: all cores).")
    p_bm.add_argument("--no-report", action="store_true")
    p_bm.add_argument("--out", required=True)
    p_bm.set_defaults(func=cmd_benchmark)

    p_sa = sub.add_parser("singlearm", help="Single-arm checkpoint cohorts: pooled log-odds of response per SD and agreement with "
                                            "the randomized interaction.")
    p_sa.add_argument("--config", required=True, help="YAML listing the cohorts (name, expr, clinical, sample_col, response_col).")
    p_sa.add_argument("--signatures")
    p_sa.add_argument("--min-cohorts", type=int, default=None)
    p_sa.add_argument("--randomized", help="benchmark_pooled.tsv from 'benchmark', for the concordance.")
    p_sa.add_argument("--no-report", action="store_true")
    p_sa.add_argument("--out", required=True)
    p_sa.set_defaults(func=cmd_singlearm)

    p_pub = sub.add_parser("published", help="Pool published subgroup hazard ratios (signature high and low) of randomized trials.")
    g3 = p_pub.add_mutually_exclusive_group(required=True)
    g3.add_argument("--table", help="One row per trial: trial, endpoint, hr_high, lo_high, hi_high, hr_low, lo_low, hi_low, ...")
    g3.add_argument("--example", action="store_true", help="The 19 comparisons extracted from published reports for the manuscript.")
    p_pub.add_argument("--no-report", action="store_true")
    p_pub.add_argument("--out", required=True)
    p_pub.set_defaults(func=cmd_published)

    p_des = sub.add_parser("design", help="Events needed to tell two correlated scores apart in a randomized trial.")
    p_des.add_argument("--delta", type=float, nargs="+", default=[0.10], help="Difference in interaction per SD.")
    p_des.add_argument("--rho", type=float, nargs="+", default=[0.55, 0.81, 0.84], help="Correlation of the two scores.")
    p_des.add_argument("--allocation", type=float, default=0.5)
    p_des.add_argument("--power", type=float, default=0.8)
    p_des.add_argument("--alpha", type=float, default=0.05)
    p_des.add_argument("--inflation", type=float, default=1.0)
    p_des.add_argument("--events", type=float, default=None, help="Also give power and detectable difference at this many events.")
    p_des.add_argument("--no-report", action="store_true")
    p_des.add_argument("--out", default=None)
    p_des.set_defaults(func=cmd_design)

    p_rep = sub.add_parser("report", help="Rebuild report.html and report.pdf from the tables in a results folder.")
    p_rep.add_argument("--dir", required=True)
    p_rep.add_argument("--title", default="barriernet report")
    p_rep.set_defaults(func=cmd_report)

    p_demo = sub.add_parser("demo", help="Run every analysis on synthetic data and write the report.")
    p_demo.add_argument("--out", required=True)
    p_demo.add_argument("--quick", action="store_true", help="100 random draws and 100 resamples.")
    p_demo.add_argument("--jobs", type=int, default=1)
    p_demo.set_defaults(func=cmd_demo)

    p_app = sub.add_parser("app", help="Open the web app (Shiny for Python).")
    p_app.add_argument("--host", default="127.0.0.1")
    p_app.add_argument("--port", type=int, default=8000)
    p_app.set_defaults(func=cmd_app)

    p_ihc = sub.add_parser("ihc", help="Compute Barrier/NetImmune axes from IHC/IF/protein panel table.")
    p_ihc.add_argument("--table", required=True, help="Long table with marker columns + sample_id column (ROI/cell rows allowed).")
    p_ihc.add_argument("--sample-col", required=True, help="Column name for sample/slide/patient ID.")
    p_ihc.add_argument("--drop-cols", nargs="*", default=None, help="Optional metadata columns to drop (besides sample_col).")
    p_ihc.add_argument("--agg", default="mean", choices=["mean","median"], help="Aggregate ROI/cell values to sample.")
    p_ihc.add_argument("--geneset", default=None, help="Path to geneset YAML (default: bundled IHC geneset).")
    p_ihc.add_argument("--min-genes", type=int, default=2)
    p_ihc.add_argument("--no-zscore", action="store_true")
    p_ihc.add_argument("--out", required=True)
    p_ihc.set_defaults(func=cmd_ihc)

    p_sp = sub.add_parser("spatial", help="Visium spatial: assign TLS-centered rings using hex-grid adjacency.")
    g = p_sp.add_mutually_exclusive_group(required=True)
    g.add_argument("--h5ad", help="Single h5ad path.")
    g.add_argument("--h5ad-glob", help="Glob for multiple h5ad files.")
    _spatial_args(p_sp)
    p_sp.add_argument("--out", required=True)
    p_sp.set_defaults(func=cmd_spatial)

    # Backward compatible alias of spatial
    p_run = sub.add_parser("run", help="Alias of 'spatial' for backward compatibility.")
    g2 = p_run.add_mutually_exclusive_group(required=True)
    g2.add_argument("--h5ad", help="Single h5ad path.")
    g2.add_argument("--h5ad-glob", help="Glob for multiple h5ad files.")
    _spatial_args(p_run)
    p_run.add_argument("--out", required=True)
    p_run.set_defaults(func=cmd_spatial)

    return p

def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
