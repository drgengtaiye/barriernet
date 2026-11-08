"""barriernet web app (Shiny for Python): benchmark published signatures in randomized trials, single-arm cohorts,
published subgroup results and trial design. Every statistic calls the package functions; data stay on the computer
running the app. Start with ``barriernet app`` (needs ``pip install 'barriernet[app]'``)."""
from __future__ import annotations
import os
import pathlib
import shutil
import tempfile
import zipfile
import numpy as np
import pandas as pd
from shiny import App, reactive, render, ui
from shiny.types import SafeException
from . import __version__
from .benchmark import benchmark
from .design import curve, detectable_difference, events_needed, power
from .library import MAIN_REFERENCE, REFERENCES, load_library, read_signature_table, reference_weights, signature_meta, signature_weights
from .published import example_table, pool_published
from .report import fig_benchmark, fig_concordance, fig_counts, fig_design, fig_paired, fig_published, fig_reference, benchmark_text, write_report
from .singlearm import concordance, single_arm

def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, default))
    except ValueError:
        return default


# Limits for a shared server; each can be changed with an environment variable.
CPU = os.cpu_count() or 2
MAX_WORKERS = max(1, _env_int("BARRIERNET_MAX_WORKERS", max(1, min(4, CPU - 1))))   # worker processes per benchmark
MAX_JOBS = max(1, _env_int("BARRIERNET_MAX_JOBS", 1))                              # benchmarks running at once
MAX_RANDOM = max(0, _env_int("BARRIERNET_MAX_RANDOM", 1000))
MAX_BOOT = max(20, _env_int("BARRIERNET_MAX_BOOT", 300))
MAX_UPLOAD_MB = max(1, _env_int("BARRIERNET_MAX_UPLOAD_MB", 200))
JOBS = MAX_WORKERS
_ACTIVE: set = set()                                                                 # benchmarks running in this process
PUB_COLS = ["trial", "endpoint", "hr_high", "lo_high", "hi_high", "hr_low", "lo_low", "hi_low", "frac_high", "from_figure", "exclude"]
PUB_NUM = {"hr_high", "lo_high", "hi_high", "hr_low", "lo_low", "hi_low", "frac_high"}
PAPER_RHO = {"r = 0.55": 0.55, "r = 0.81": 0.81, "r = 0.84": 0.84}

CSS = """
body{font-family:-apple-system,Helvetica,Arial,sans-serif}
.bn-lead{font-size:1rem;margin:.25rem 0 1rem;max-width:62rem}
.bn-note{color:#5d6577;font-size:.88rem;max-width:62rem}
.bn-num{font-size:1.9rem;font-weight:600;color:#3C5488;line-height:1.1}
.bn-numlab{color:#5d6577;font-size:.85rem;margin-bottom:.8rem}
.bn-synthetic{border-left:3px solid #e7213a;padding-left:.6rem}
"""


def _sidebar_bm():
    return ui.sidebar(
        ui.input_radio_buttons("bm_src", "Randomized trials", {"demo": "Synthetic demonstration (not patient data)", "zip": "Upload a zip"}),
        ui.panel_conditional("input.bm_src === 'zip'",
                             ui.input_file("bm_zip", "Zip with trials.yaml and the tables it names", accept=[".zip"]),
                             ui.p("Same YAML layout as the command line (README, 'benchmark').", class_="bn-note")),
        ui.input_file("bm_sigs", "Your signatures (optional)", accept=[".tsv", ".csv", ".txt"]),
        ui.p("Long table: signature, gene, and optionally weight, direction (1 or -1), family.", class_="bn-note"),
        ui.input_numeric("bm_nrand", "Random gene sets per signature", min(200, MAX_RANDOM), min=0, max=MAX_RANDOM, step=100),
        ui.input_numeric("bm_nboot", "Patient resamples per trial", min(200, MAX_BOOT), min=20, max=MAX_BOOT, step=50),
        ui.input_numeric("bm_seed", "Seed", 0),
        ui.input_numeric("bm_jobs", "Worker processes", JOBS, min=1, max=MAX_WORKERS),
        ui.input_action_button("bm_run", "Run the benchmark", class_="btn-primary"),
        ui.p("The manuscript used 1,000 random sets and 300 resamples; fewer run faster with more Monte Carlo error.", class_="bn-note mt-2"),
        width=330)


app_ui = ui.page_navbar(
    ui.nav_panel(
        "Benchmark",
        ui.layout_sidebar(
            _sidebar_bm(),
            ui.p("Each signature is scored in every randomized comparison, its score × arm interaction is pooled, and the result is set "
                 "against the 13-gene T-cell mean, the CD8A-PRF1-GZMA mean, single genes and size-matched random genes.", class_="bn-lead"),
            ui.output_ui("bm_progress"),
            ui.output_ui("bm_text"),
            ui.navset_underline(
                ui.nav_panel("Tests passed", ui.output_plot("bm_counts", height="420px")),
                ui.nav_panel("Signatures", ui.output_plot("bm_forest", height="980px")),
                ui.nav_panel("Against a reference", ui.input_select("bm_ref", None, choices=list(REFERENCES) + ["CD8A"], selected=MAIN_REFERENCE),
                             ui.output_plot("bm_paired", height="980px")),
                ui.nav_panel("Against random genes", ui.output_plot("bm_reference", height="480px")),
                ui.nav_panel("Table", ui.output_data_frame("bm_table")),
                id="bm_tabs"),
            ui.output_ui("bm_dl_ui"))),
    ui.nav_panel(
        "Single-arm cohorts",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_radio_buttons("sa_src", "Cohorts", {"demo": "Synthetic demonstration (not patient data)", "zip": "Upload a zip"}),
                ui.panel_conditional("input.sa_src === 'zip'", ui.input_file("sa_zip", "Zip with cohorts.yaml and the tables it names", accept=[".zip"])),
                ui.input_numeric("sa_min", "Pool when scored in at least (cohorts; blank: half)", None, min=2),
                ui.input_action_button("sa_run", "Pool the cohorts", class_="btn-primary"), width=330),
            ui.p("Log-odds of response per SD of each score, cohort by cohort (logistic regression with Firth's correction when the fit "
                 "separates), pooled by random effects. A single-arm association cannot tell a predictive score from a prognostic one; "
                 "the comparison with the randomized interaction shows how far the two rankings agree.", class_="bn-lead"),
            ui.output_ui("sa_text"),
            ui.output_plot("sa_conc", height="440px"),
            ui.output_data_frame("sa_table"))),
    ui.nav_panel(
        "Published trials",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_action_button("pub_example", "Load the 19 published comparisons"),
                ui.input_file("pub_file", "Or upload a table", accept=[".tsv", ".csv", ".txt"]),
                ui.input_action_button("pub_add", "Add an empty row"),
                ui.download_button("pub_dl_table", "Download the table"),
                ui.download_button("pub_dl_report", "Download the report (zip)"),
                ui.p("Hazard ratio, checkpoint arm against control, with its 95% CI, separately for patients with a high and a low signature. "
                     "frac_high: share classed high (blank: half). from_figure: yes if read off a figure. exclude: yes to leave a row out.",
                     class_="bn-note mt-2"), width=330),
            ui.p("The benefit-oriented log ratio of the two hazard ratios is pooled by REML with the Knapp-Hartung adjustment. "
                 "Edit any cell; the pooled result updates at once.", class_="bn-lead"),
            ui.output_ui("pub_text"),
            ui.output_data_frame("pub_grid"),
            ui.output_plot("pub_forest", height="620px"),
            ui.output_data_frame("pub_analyses"))),
    ui.nav_panel(
        "Trial design",
        ui.layout_sidebar(
            ui.sidebar(
                ui.input_slider("d_delta", "Difference between the two interactions (per SD)", 0.02, 0.40, 0.10, step=0.01),
                ui.input_slider("d_rho", "Correlation between the two scores", 0.0, 0.99, 0.81, step=0.01),
                ui.input_slider("d_alloc", "Share randomized to the checkpoint arm", 0.2, 0.8, 0.5, step=0.05),
                ui.input_slider("d_power", "Power", 0.5, 0.95, 0.8, step=0.05),
                ui.input_numeric("d_alpha", "Two-sided alpha", 0.05, min=0.001, max=0.2, step=0.005),
                ui.input_numeric("d_events", "Events available (optional)", None, min=10),
                width=330),
            ui.p("How many events a randomized survival comparison needs before it can tell two correlated scores apart: "
                 "2(1 − r)(z₁₋α/₂ + z_power)² / (δ² p(1 − p)), by normal theory.", class_="bn-lead"),
            ui.output_ui("d_numbers"),
            ui.output_plot("d_curve", height="440px"))),
    ui.nav_panel(
        "Signature library",
        ui.p("The 51 distinct published immuno-oncology signatures benchmarked (SignatureSets, BHKLAB, MIT License), with their published "
             "weights, declared direction and family. Select a row to see its genes.", class_="bn-lead"),
        ui.layout_columns(ui.output_data_frame("lib_table"), ui.output_data_frame("lib_genes"), col_widths=(8, 4)),
        ui.download_button("lib_dl", "Download the library")),
    ui.nav_panel(
        "About",
        ui.div(
            ui.h4(f"barriernet {__version__}"),
            ui.p("These are the analyses of a benchmark of published immune signatures in randomized checkpoint-inhibitor trials. "
                 "A score predicts benefit from a checkpoint inhibitor only if its association with outcome differs between the arms of a "
                 "randomized comparison, so the benchmark tests the score × arm interaction and then asks whether a signature does better than "
                 "simple T-cell scores, than single genes and than random genes of the same size."),
            ui.p("Scores are sum(w z) / sum(|w|) over the listed genes present (genes z-scored within each comparison across arms; missing "
                 "values skipped sample by sample), scored when at least three genes and half the list are present. Interactions are "
                 "benefit-oriented: −log HR per SD for survival, log OR per SD for response. Pooling uses the wider of the Hartung-Knapp and "
                 "DerSimonian-Laird intervals and the larger P. Paired differences use within-arm patient resamples common to both scores. "
                 "q values are Benjamini-Hochberg within each family of tests."),
            ui.p("The synthetic data are simulated from a few latent cell programmes with a known interaction on the T-cell programme; they are "
                 "not patient data. The 19 published comparisons were extracted from the trial reports; check them against the sources before "
                 "relying on them.", class_="bn-note"),
            style="max-width:60rem")),
    title=f"barriernet {__version__}", id="tabs", header=ui.tags.style(CSS), window_title="barriernet")


def _read_any(path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python")


def _too_big(files) -> bool:
    if files and any(int(f.get("size") or 0) > MAX_UPLOAD_MB * 1024 * 1024 for f in files):
        ui.notification_show(f"Files above {MAX_UPLOAD_MB} MB are not accepted here; run the command line instead.", type="error")
        return True
    return False


def _clamp(v, lo, hi, default):
    try:
        return int(min(hi, max(lo, int(v))))
    except (TypeError, ValueError):
        return default


def _opt(input, name):
    """An optional input, or None while the browser has not sent it (reading it would silently cancel the caller)."""
    v = getattr(input, name)
    return v() if v.is_set() else None


def server(input, output, session):
    bm = reactive.value(None)          # dict(B=..., synthetic=bool)
    sa = reactive.value(None)
    pub = reactive.value(example_table())

    # ---------------------------------------------------------------- benchmark
    job = reactive.value(None)
    prog = reactive.value(None)
    _mine: list = []                                                                  # this session's job, for cleanup
    _tmp: list = []                                                                   # this session's download folders

    def _tmpdir(prefix: str) -> pathlib.Path:
        d = pathlib.Path(tempfile.mkdtemp(prefix=prefix)); _tmp.append(d); return d

    def _end():
        for j in _mine:
            if not j.done():
                j.cancel()
            _ACTIVE.discard(j); j.cleanup()
        for d in _tmp:
            shutil.rmtree(d, ignore_errors=True)
    session.on_ended(_end)

    @reactive.effect
    @reactive.event(input.bm_run)
    def _bm():
        if job() is not None:
            ui.notification_show("A benchmark is already running.", type="warning"); return
        for j in [j for j in _ACTIVE if j.done()]:
            _ACTIVE.discard(j)
        if len(_ACTIVE) >= MAX_JOBS:
            ui.notification_show("Another benchmark is running on this server. Please try again in a few minutes.", type="warning"); return
        synthetic = input.bm_src() == "demo"; z = _opt(input, "bm_zip"); sg = _opt(input, "bm_sigs")
        if not synthetic and not z:
            ui.notification_show("Upload a zip first.", type="warning"); return
        if _too_big((z if not synthetic else None) or []) or _too_big(sg or []):
            return
        spec = dict(source="demo" if synthetic else "zip", zip=z[0]["datapath"] if z and not synthetic else None,
                    signatures=sg[0]["datapath"] if sg else None, n_random=_clamp(input.bm_nrand(), 0, MAX_RANDOM, 200),
                    n_boot=_clamp(input.bm_nboot(), 20, MAX_BOOT, 200), seed=_clamp(input.bm_seed(), 0, 2 ** 31 - 1, 0),
                    n_jobs=_clamp(input.bm_jobs(), 1, MAX_WORKERS, 1),
                    synthetic=synthetic)
        from .jobs import Job
        j = Job(spec); _ACTIVE.add(j); _mine.clear(); _mine.append(j)
        bm.set(None); prog.set(("starting", 0, 1)); job.set(j)

    @reactive.effect
    def _poll():
        j = job()
        if j is None:
            return
        reactive.invalidate_later(1)
        if not j.done():
            prog.set(j.progress()); return
        try:
            bm.set(dict(B=j.result(), synthetic=j.spec["synthetic"]))
        except Exception as e:
            ui.notification_show(f"The benchmark failed: {e}", type="error", duration=None)
        _ACTIVE.discard(j); j.cleanup(); job.set(None); prog.set(None)

    @reactive.effect
    @reactive.event(input.bm_cancel)
    def _cancel():
        j = job()
        if j is not None:
            j.cancel(); _ACTIVE.discard(j); j.cleanup(); job.set(None); prog.set(None); ui.notification_show("Benchmark cancelled.")

    @render.ui
    def bm_progress():
        p = prog()
        if job() is None or p is None:
            return None
        stage, i, n = p
        return ui.div(ui.tags.progress(value=str(i), max=str(max(n, 1)), style="width:100%;height:10px"),
                      ui.p(f"{stage.capitalize()}: {i} of {n}" if n > 1 else stage.capitalize(), class_="bn-note"),
                      ui.input_action_button("bm_cancel", "Cancel", class_="btn-sm btn-outline-secondary"), class_="my-2")

    @render.ui
    def bm_text():
        v = bm()
        if v is None:
            return None if job() is not None else ui.p("Choose the data and press Run. With the settings shown, the synthetic "
                                                       "demonstration takes a few minutes.", class_="bn-note")
        out = [ui.p(benchmark_text(v["B"]))]
        if v["synthetic"]:
            out.insert(0, ui.p("Synthetic demonstration: simulated trials, not patient data.", class_="bn-note bn-synthetic"))
        return ui.div(*out)

    def _need():
        v = bm()
        if v is None:
            raise SafeException("Run the benchmark first.")
        return v["B"]

    @render.plot
    def bm_counts():
        return fig_counts(_need())

    @render.plot
    def bm_forest():
        return fig_benchmark(_need())

    @render.plot
    def bm_paired():
        return fig_paired(_need(), input.bm_ref())

    @render.plot
    def bm_reference():
        return fig_reference(_need())

    @render.data_frame
    def bm_table():
        P = _need()["pooled"]; cols = [c for c in ("score", "kind", "family", "direction", "trials_scored", "est", "lo", "hi", "p", "q", "emp_p_declared", "q_reference") if c in P]
        return render.DataGrid(P[cols].sort_values("est", ascending=False).round(4), width="100%", filters=True, summary=False)

    @render.ui
    def bm_dl_ui():
        return ui.download_button("bm_dl", "Download the results and report (zip)", class_="mt-3") if bm() is not None else None

    @render.download_button(filename="barriernet_benchmark.zip")
    def bm_dl():
        v = bm()
        if v is None:
            raise SafeException("Run the benchmark first.")
        d = _tmpdir("barriernet_bm_")
        for k, df in v["B"].items():
            df.to_csv(d / f"benchmark_{k}.tsv", sep="\t", index=False)
        write_report(d, benchmark=v["B"], title="Benchmark of published signatures",
                     note="Synthetic demonstration: simulated trials, not patient data." if v["synthetic"] else "")
        z = _tmpdir("barriernet_zip_") / (d.name + ".zip")
        with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as f:
            for p in d.rglob("*"):
                if p.is_file():
                    f.write(p, p.relative_to(d))
        return str(z)

    # ---------------------------------------------------------------- single-arm
    @reactive.effect
    @reactive.event(input.sa_run)
    def _sa():
        synthetic = input.sa_src() == "demo"
        with ui.Progress(min=0, max=1) as p:
            p.set(0.2, message="Scoring the cohorts")
            try:
                if synthetic:
                    from .demo import demo_cohorts
                    cohorts = demo_cohorts()
                else:
                    if not _opt(input, "sa_zip"):
                        ui.notification_show("Upload a zip first.", type="warning"); return
                    if _too_big(_opt(input, "sa_zip")):
                        return
                    from .bundle import load_cohorts_zip
                    cohorts, _ = load_cohorts_zip(_opt(input, "sa_zip")[0]["datapath"])
            except Exception as e:
                ui.notification_show(f"Could not read the cohorts: {e}", type="error"); return
            S = single_arm(cohorts, {**signature_weights(load_library()), **reference_weights()},
                           min_cohorts=int(_opt(input, "sa_min")) if _opt(input, "sa_min") else None)
        sa.set(dict(S=S, synthetic=synthetic, n=len(cohorts)))

    @render.ui
    def sa_text():
        v = sa()
        if v is None:
            return ui.p("Choose the cohorts and press Pool.", class_="bn-note")
        P = v["S"]["pooled"]; out = []
        if v["synthetic"]:
            out.append(ui.p("Synthetic demonstration: simulated cohorts, not patient data.", class_="bn-note bn-synthetic"))
        out.append(ui.p(f"{len(P)} scores pooled over {v['n']} cohorts; {int((P.p < 0.05).sum())} at P < 0.05."))
        if bm() is None:
            out.append(ui.p("Run the benchmark as well to compare the single-arm ranking with the randomized interaction.", class_="bn-note"))
        return ui.div(*out)

    @render.plot
    def sa_conc():
        v, b = sa(), bm()
        if v is None or b is None:
            raise SafeException("Needs both the single-arm pooling and the benchmark.")
        c = concordance(v["S"]["pooled"], b["B"]["pooled"])
        if c["n"] < 3:
            raise SafeException("Fewer than three scores pooled in both.")
        return fig_concordance(c)

    @render.data_frame
    def sa_table():
        v = sa()
        if v is None:
            return None
        return render.DataGrid(v["S"]["pooled"][["score", "k", "est", "lo", "hi", "p", "tau2"]].sort_values("est", ascending=False).round(4),
                               width="100%", filters=True, summary=False)

    # ---------------------------------------------------------------- published trials
    @reactive.effect
    @reactive.event(input.pub_example)
    def _pub_ex():
        pub.set(example_table())

    @reactive.effect
    @reactive.event(input.pub_file)
    def _pub_file():
        try:
            f = _opt(input, "pub_file")
            if not f or _too_big(f):
                return
            T = _read_any(f[0]["datapath"]); pool_published(T); pub.set(T)
        except Exception as e:
            ui.notification_show(f"Could not use the table: {e}", type="error")

    @reactive.effect
    @reactive.event(input.pub_add)
    def _pub_add():
        T = pub().copy(); T.loc[len(T)] = {c: ("OS" if c == "endpoint" else ("no" if c in ("from_figure", "exclude") else (f"Trial {len(T) + 1}" if c == "trial" else np.nan))) for c in T.columns}
        pub.set(T)

    def _pub_res():
        T = pub()
        try:
            return pool_published(T)
        except Exception:
            return None

    @render.data_frame
    def pub_grid():
        T = pub()
        for c in PUB_COLS:
            if c not in T:
                T[c] = np.nan
        return render.DataGrid(T[PUB_COLS], width="100%", editable=True, summary=False, height="360px")

    @pub_grid.set_patch_fn
    def _pub_patch(*, patch):
        col = PUB_COLS[patch["column_index"]]; i = patch["row_index"]; v = patch["value"]
        if col in PUB_NUM:
            if v in ("", None):
                v = np.nan
            else:
                try:
                    v = float(str(v).strip())
                except ValueError:
                    raise SafeException("Enter a number, for example 0.78.")
                if v <= 0 or (col == "frac_high" and v >= 1):
                    raise SafeException("Hazard ratios and limits must be above 0; frac_high between 0 and 1.")
        elif col == "endpoint" and str(v).strip() not in ("OS", "PFS", "EFS", "DFS", "pCR", "ORR"):
            raise SafeException("Endpoint: OS or PFS are pooled; other endpoints are kept but not pooled.")
        T = pub().copy(); T.loc[T.index[i], col] = v; pub.set(T)
        return "" if (isinstance(v, float) and np.isnan(v)) else v

    @render.ui
    def pub_text():
        r = _pub_res()
        if r is None or "primary" not in r:
            return ui.p("Fewer than two complete survival rows: nothing to pool yet.", class_="bn-note")
        m = r["primary"]; e, lo, hi = np.exp(-m["est"]), np.exp(-m["hi"]), np.exp(-m["lo"])
        s = (f"Over {m['k']} trials the ratio of hazard ratios, signature high over low, was {e:.2f} ({lo:.2f} to {hi:.2f}), P = {m['p']:.4f}; "
             f"I² = {m['I2']:.0f}%.")
        if "per_sd" in r:
            s += f" Per SD of the signature: {np.exp(-r['per_sd']['est']):.2f}."
        if "leave_one_out" in r:
            s += f" Leaving out one trial at a time gave log ratios from {r['leave_one_out']['lo']:.3f} to {r['leave_one_out']['hi']:.3f}."
        return ui.p(s)

    @render.plot
    def pub_forest():
        r = _pub_res()
        if r is None or "primary" not in r:
            raise SafeException("Nothing to plot yet.")
        return fig_published(r)

    @render.data_frame
    def pub_analyses():
        r = _pub_res()
        if r is None:
            return None
        rows = [dict(analysis=k, k=v.get("k"), est=v.get("est"), lo=v.get("lo"), hi=v.get("hi"), p=v.get("p"), I2=v.get("I2"))
                for k, v in r.items() if isinstance(v, dict)]
        return render.DataGrid(pd.DataFrame(rows).round(4), width="100%", summary=False)

    @render.download_button(filename="published_trials.tsv")
    def pub_dl_table():
        d = _tmpdir("barriernet_pub_") / "published_trials.tsv"; pub().to_csv(d, sep="\t", index=False); return str(d)

    @render.download_button(filename="barriernet_published.zip")
    def pub_dl_report():
        r = _pub_res()
        if r is None or "primary" not in r:
            raise SafeException("Nothing to pool yet.")
        d = _tmpdir("barriernet_pub_"); pub().to_csv(d / "published_trials.tsv", sep="\t", index=False)
        r["per_trial"].to_csv(d / "published_per_trial.tsv", sep="\t", index=False); write_report(d, published=r, title="Published randomized trials")
        z = _tmpdir("barriernet_zip_") / (d.name + ".zip")
        with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as f:
            for p in d.rglob("*"):
                if p.is_file():
                    f.write(p, p.relative_to(d))
        return str(z)

    # ---------------------------------------------------------------- design
    @render.ui
    def d_numbers():
        d, r, a, pw, al = input.d_delta(), input.d_rho(), input.d_alloc(), input.d_power(), float(input.d_alpha() or 0.05)
        n = events_needed(d, r, a, pw, al)
        items = [ui.div(ui.div(f"{n:,.0f}", class_="bn-num"), ui.div(f"events to detect {d:.2f} per SD at r = {r:.2f} with {pw:.0%} power", class_="bn-numlab"))]
        ev = _opt(input, "d_events")
        if ev:
            items.append(ui.div(ui.div(f"{power(ev, d, r, a, al):.0%}", class_="bn-num"), ui.div(f"power with {ev:,.0f} events", class_="bn-numlab")))
            items.append(ui.div(ui.div(f"{detectable_difference(ev, r, a, pw, al):.2f}", class_="bn-num"),
                                ui.div(f"difference per SD detectable with {ev:,.0f} events", class_="bn-numlab")))
        return ui.layout_columns(*items)

    @render.plot
    def d_curve():
        d = input.d_delta(); deltas = sorted({round(d, 2), 0.10, 0.15, 0.20})
        C = curve(deltas=tuple(deltas), allocation=input.d_alloc(), power_=input.d_power(), alpha=float(input.d_alpha() or 0.05))
        return fig_design(C, {**PAPER_RHO, "chosen": input.d_rho()}, letter="")

    # ---------------------------------------------------------------- library
    LIB = load_library(); META = signature_meta(LIB)

    @render.data_frame
    def lib_table():
        M = META.assign(direction=META.direction.map({1: "sensitivity", -1: "resistance"}))
        return render.DataGrid(M[["signature", "family", "direction", "n_genes", "pmid"]], width="100%", filters=True, selection_mode="row", summary=False, height="520px")

    @render.data_frame
    def lib_genes():
        sel = lib_table.cell_selection()["rows"] if lib_table.cell_selection() else ()
        if not sel:
            return None
        s = lib_table.data_view(selected=True)["signature"].iloc[0]
        return render.DataGrid(LIB[LIB.signature == s][["gene", "weight"]], width="100%", summary=False, height="520px")

    @render.download_button(filename="barriernet_signatures.tsv")
    def lib_dl():
        d = _tmpdir("barriernet_lib_") / "signatures.tsv"; LIB.to_csv(d, sep="\t", index=False); return str(d)


app = App(app_ui, server)
