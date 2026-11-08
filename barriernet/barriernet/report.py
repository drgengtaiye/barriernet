"""Figures in the manuscript's style and a report (HTML with the figures embedded, and a PDF of the same pages).

Needs matplotlib (``pip install barriernet[report]``). Each ``fig_*`` function takes the output of one analysis and
returns a matplotlib Figure; ``write_report`` puts whichever analyses it is given into report.html and report.pdf."""
from __future__ import annotations
import base64
import datetime as _dt
import html
import io
import logging
import pathlib
import warnings
from typing import Dict, Optional
import numpy as np
import pandas as pd

NAVY, RED, TEAL, ORANGE = "#3C5488", "#e7213a", "#008573", "#E18727"
INK, MUTED, GREY, LIGHT = "#1a1a1a", "#595959", "#a6a6a6", "#e0e0e0"
FAMILY_COLOURS = {"Inflamed/IFN": NAVY, "Other": GREY, "Confounder": ORANGE, "User": RED}
MM = 1 / 25.4


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR); logging.getLogger("fontTools").setLevel(logging.ERROR)
    warnings.filterwarnings("ignore", message="Glyph .* missing")
    plt.rcParams.update({"font.family": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7, "axes.titlesize": 7.5,
                         "axes.labelsize": 7, "xtick.labelsize": 6.6, "ytick.labelsize": 6.6, "legend.fontsize": 6.6,
                         "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2,
                         "ytick.major.size": 2, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                         "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK})
    return plt


def _logticks(ax, axis, ticks):
    from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
    a = ax.xaxis if axis == "x" else ax.yaxis
    a.set_major_locator(FixedLocator(ticks)); a.set_minor_locator(NullLocator()); a.set_minor_formatter(NullFormatter())
    a.set_ticklabels([f"{t:,.0f}" if t >= 10 else f"{t:g}" for t in ticks])


def _title(ax, letter, text):
    ax.set_title(f"{letter}  {text}" if letter else text, loc="left", fontweight="bold", pad=6)


# ------------------------------------------------------------------ benchmark
def fig_benchmark(B: Dict[str, pd.DataFrame], letter: str = "a"):
    """Pooled score x arm interaction per signature, with the T-cell references marked."""
    plt = _plt(); P = B["pooled"]; S = P[(P.kind == "signature") & P.est.notna()].sort_values("est")
    refs = P[(P.kind == "T-cell reference") & P.est.notna()]
    fig, ax = plt.subplots(figsize=(88 * MM, max(60, 2.6 * len(S) + 22) * MM))
    y = np.arange(len(S))
    for i, (_, r) in enumerate(S.iterrows()):
        c = FAMILY_COLOURS.get(r.family, GREY)
        ax.plot([r.lo, r.hi], [i, i], color=c, lw=0.6)
        ax.plot(r.est, i, "o", ms=2.6, color=c, mfc=c if r.p < 0.05 else "white", mew=0.6)
    for (_, r), ls in zip(refs.iterrows(), ("-", "--", ":")):
        ax.axvline(r.est, color=RED, lw=0.6, ls=ls, label=f"{r.score}: {r.est:+.2f}")
    ax.axvline(0, color=INK, lw=0.5)
    ax.set_yticks(y); ax.set_yticklabels(S.score, fontsize=6.6); ax.set_ylim(-1, len(S))
    ax.set_xlabel("Pooled interaction per SD (positive: more checkpoint benefit)")
    for fam, c in FAMILY_COLOURS.items():
        if fam in set(S.family):
            ax.plot([], [], "o", ms=3, color=c, label=fam)
    ax.plot([], [], "o", ms=3, color=GREY, mfc="white", label="open: P ≥ 0.05")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.4, -0.065 * 60 / max(60, len(S) + 10)), ncol=2, handlelength=1.6)
    _title(ax, letter, "Signatures in the randomized comparisons"); fig.tight_layout(); return fig


def fig_paired(B: Dict[str, pd.DataFrame], reference: Optional[str] = None, letter: str = "b"):
    """Paired difference of each signature from a reference, with its 95% CI and the difference detectable at 80% power."""
    from .library import MAIN_REFERENCE
    plt = _plt(); reference = reference or MAIN_REFERENCE; D = B["paired"]
    D = D[D.reference == reference].sort_values("diff") if len(D) else D
    fig, ax = plt.subplots(figsize=(88 * MM, max(60, 2.6 * len(D) + 22) * MM))
    if not len(D):
        ax.text(0.5, 0.5, "no paired comparisons", ha="center", transform=ax.transAxes); return fig
    md = float(np.nanmedian(D["mdd"])); ax.axvspan(-md, md, color=LIGHT, lw=0)
    for i, (_, r) in enumerate(D.iterrows()):
        c = NAVY if r["diff"] > 0 else RED
        ax.plot([r.lo, r.hi], [i, i], color=c if r.p < 0.05 else GREY, lw=0.6)
        ax.plot(r["diff"], i, "o", ms=2.6, color=c if r.p < 0.05 else GREY, mfc="white" if r.p >= 0.05 else None, mew=0.6)
    ax.axvline(0, color=INK, lw=0.5); ax.set_yticks(range(len(D))); ax.set_yticklabels(D.score, fontsize=6.6); ax.set_ylim(-1, len(D))
    ax.set_xlabel(f"Difference in interaction per SD\ngrey: detectable at 80% power (median ±{md:.2f})", loc="right")
    _title(ax, letter, f"Against {reference.replace('T-cell mean', 'the T-cell mean')}"); fig.tight_layout(); return fig


def fig_reference(B: Dict[str, pd.DataFrame], letter: str = "c"):
    """Observed pooled interaction against the 95th percentile of size-matched random gene sets (declared direction)."""
    plt = _plt(); P = B["pooled"]
    S = P[(P.kind == "signature") & P.get("ref_p95_declared", pd.Series(np.nan, index=P.index)).notna()] if "ref_p95_declared" in P else P.iloc[:0]
    fig, ax = plt.subplots(figsize=(72 * MM, 72 * MM))
    if not len(S):
        ax.text(0.5, 0.5, "random reference not run", ha="center", transform=ax.transAxes); return fig
    obs, ref = S.est * S.direction, S.ref_p95_declared * S.direction
    hit = S.emp_p_declared < 0.05
    ax.scatter(ref[~hit], obs[~hit], s=9, facecolors="white", edgecolors=GREY, lw=0.6, label="within random sets")
    ax.scatter(ref[hit], obs[hit], s=9, color=NAVY, lw=0, label="beyond, empirical P < 0.05")
    lim = [min(ref.min(), obs.min()) - 0.03, max(ref.max(), obs.max()) + 0.03]
    ax.plot(lim, lim, color=INK, lw=0.5); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("95th percentile, size-matched random genes"); ax.set_ylabel("Observed, declared direction")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, handletextpad=0.2, columnspacing=1.0)
    _title(ax, letter, "Against random genes"); fig.tight_layout(); return fig


def fig_counts(B: Dict[str, pd.DataFrame], letter: str = "d"):
    """How many signatures pass each test, nominally and after Benjamini-Hochberg correction."""
    plt = _plt(); S = B["summary"].iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(88 * MM, (4.2 * len(S) + 18) * MM))
    for i, r in S.iterrows():
        ax.plot([0, r.nominal], [i, i], color=LIGHT, lw=2.2, solid_capstyle="butt")
        ax.plot(r.nominal, i, "o", ms=3.6, color=NAVY, mfc="white", mew=0.7)
        if r.fdr is not None and np.isfinite(r.fdr):
            ax.plot(r.fdr, i, "o", ms=3.6, color=NAVY)
        ax.text(r.nominal + 0.022 * (int(S["of"].max()) + 2) + 0.4, i, f"{int(r.nominal)}", va="center", fontsize=6.6, color=MUTED)
    ax.set_yticks(range(len(S))); ax.set_yticklabels(S.test); ax.set_xlim(0, max(1, int(S["of"].max())) + 2)
    ax.set_xlabel(f"Signatures (of {int(S['of'].max())})")
    ax.plot([], [], "o", color=NAVY, mfc="white", ms=3.6, label="P < 0.05"); ax.plot([], [], "o", color=NAVY, ms=3.6, label="q < 0.05")
    ax.legend(frameon=False, loc="lower right"); _title(ax, letter, "Tests passed"); fig.tight_layout(); return fig


# ------------------------------------------------------------------ single-arm, published trials, design
def fig_concordance(conc: dict, letter: str = "e"):
    """Pooled single-arm log-odds against the pooled randomized interaction (one dot per score)."""
    plt = _plt(); T = conc["table"]; fig, ax = plt.subplots(figsize=(72 * MM, 72 * MM))
    both = (T.interaction_p < 0.05) & (T.single_arm_p < 0.05)
    ax.scatter(T.interaction[~both], T.single_arm[~both], s=9, facecolors="white", edgecolors=GREY, lw=0.6, label="not both P < 0.05")
    ax.scatter(T.interaction[both], T.single_arm[both], s=9, color=NAVY, lw=0, label="both P < 0.05")
    ax.axhline(0, color=INK, lw=0.4); ax.axvline(0, color=INK, lw=0.4)
    ax.set_xlabel("Randomized interaction per SD"); ax.set_ylabel("Single-arm log-odds of response per SD")
    ax.text(0.98, 0.04, f"Spearman {conc['rho']:.2f} ({conc['n']} scores)", transform=ax.transAxes, ha="right", fontsize=6.6, color=MUTED)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, handletextpad=0.2, columnspacing=1.0)
    _title(ax, letter, "Single-arm against randomized"); fig.tight_layout(); return fig


def fig_published(res: dict, letter: str = "f"):
    """Forest plot of the published subgroup results as ratios of hazard ratios, with the pooled estimate."""
    plt = _plt(); I = res["per_trial"]; P = I[I["survival"] & ~I["exclude"] & np.isfinite(I["y"])].sort_values("y")
    z = 1.959963984540054; fig, ax = plt.subplots(figsize=(88 * MM, (3.4 * len(P) + 26) * MM))
    for i, (_, r) in enumerate(P.iterrows()):
        ax.plot(np.exp(-np.array([r.y + z * r.se, r.y - z * r.se])), [i, i], color=NAVY, lw=0.6)
        ax.plot(np.exp(-r.y), i, "s", ms=2.8, color=NAVY, mfc="white" if r.figure else NAVY, mew=0.6)
    k = len(P); lab = "pooled"
    if "primary" in res:
        m = res["primary"]; e, l, h = np.exp(-m["est"]), np.exp(-m["hi"]), np.exp(-m["lo"]); yy = k + 0.6
        ax.fill([l, e, h, e], [yy, yy + 0.35, yy, yy - 0.35], color=RED)
        lab = f"pooled, {int(m['k'])} trials"
        ax.set_xlabel(f"Ratio of hazard ratios, signature high over low\npooled {e:.2f} ({l:.2f}–{h:.2f}), P = {m['p']:.3f}")
    else:
        ax.set_xlabel("Ratio of hazard ratios, signature high over low")
    ax.axvline(1, color=INK, lw=0.5); ax.set_xscale("log")
    xl, xh = float(np.exp(-(P.y + z * P.se)).min()) * 0.85, float(np.exp(-(P.y - z * P.se)).max()) * 1.15
    ax.set_xlim(min(xl, 0.5), max(xh, 2)); _logticks(ax, "x", [t for t in (0.1, 0.25, 0.5, 1, 2, 4, 8) if min(xl, 0.5) <= t <= max(xh, 2)])
    ax.set_yticks(list(range(k)) + [k + 0.6]); ax.set_yticklabels(list(P.trial) + [lab]); ax.set_ylim(k + 1.4, -0.8)
    ax.plot([], [], "s", color=NAVY, mfc="white", ms=2.8, label="read from a figure"); ax.plot([], [], "s", color=NAVY, ms=2.8, label="text or table")
    ax.legend(frameon=False, loc="upper left", borderaxespad=0.1, handletextpad=0.3)
    _title(ax, letter, "Published randomized trials"); fig.tight_layout(); return fig


def fig_design(curve: pd.DataFrame, marks: Optional[Dict[str, float]] = None, letter: str = "g"):
    """Events needed against the correlation of the two scores, for each difference."""
    plt = _plt(); fig, ax = plt.subplots(figsize=(80 * MM, 64 * MM)); cols = [NAVY, RED, TEAL, ORANGE, GREY]
    for c, (d, g) in zip(cols, curve.groupby("delta")):
        ax.plot(g.rho, g.events, color=c, lw=0.9, label=f"difference {d:.2f} per SD")
    ax.set_yscale("log"); lo_, hi_ = curve.events.min(), curve.events.max()
    ticks = [t for t in (50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000) if lo_ * 0.8 <= t <= hi_ * 1.25]
    _logticks(ax, "y", ticks or [round(lo_), round(hi_)])
    for lab, r in (marks or {}).items():
        ax.axvline(r, color=GREY, lw=0.5, ls="--"); ax.text(r, 1.0, lab, transform=ax.get_xaxis_transform(), rotation=90, ha="right", va="top", fontsize=6.6, color=MUTED)
    ax.set_xlabel("Correlation between the two scores"); ax.set_ylabel("Events needed (80% power, two-sided 0.05)")
    ax.legend(frameon=False, loc="lower left"); _title(ax, letter, "Trial size to separate two scores"); fig.tight_layout(); return fig


# ------------------------------------------------------------------ report
CSS = """body{font:14px/1.55 -apple-system,Helvetica,Arial,sans-serif;color:#1a1a1a;background:#fff;margin:0}
main{max-width:860px;margin:0 auto;padding:24px 16px 64px}h1{font-size:22px;margin:0 0 4px}h2{font-size:17px;margin:32px 0 6px;border-bottom:1px solid #e0e0e0;padding-bottom:4px}
p.meta{color:#595959;margin:0 0 16px}img{max-width:100%;height:auto;display:block;margin:10px 0}
table{border-collapse:collapse;font-size:12px;margin:8px 0;width:100%}th,td{padding:3px 6px;border-bottom:1px solid #eee;text-align:right}th:first-child,td:first-child{text-align:left}
.note{color:#595959;font-size:13px}@media (prefers-color-scheme:dark){body{background:#fff}}"""


def _png(fig) -> str:
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=200, bbox_inches="tight"); return base64.b64encode(buf.getvalue()).decode()


def _table(df: pd.DataFrame, n: int = 60) -> str:
    return df.head(n).to_html(index=False, float_format=lambda v: f"{v:.3g}", border=0, na_rep="")


def benchmark_text(B) -> str:
    S = B["summary"].set_index("test"); st = B["settings"].iloc[0]; P = B["pooled"]; n = int((P.kind == "signature").sum())
    k = int(P.loc[P.kind == "signature", "est"].notna().sum())
    fav = S.loc["favor the checkpoint arm"] if "favor the checkpoint arm" in S.index else None
    ref = P[P.score == st.main_reference]
    s = (f"{k} of {n} signatures could be scored in at least {int(st.min_trials)} of the {int(st.trials)} comparisons and were pooled. ")
    if fav is not None:
        s += f"{int(fav.nominal)} favoured the checkpoint arm at P < 0.05 ({int(fav.fdr)} after Benjamini-Hochberg correction). "
    if len(ref) and np.isfinite(ref.est.iloc[0]):
        r = ref.iloc[0]; s += f"The {st.main_reference} gave {r.est:+.2f} per SD ({r.lo:+.2f} to {r.hi:+.2f}). "
    a, b = f"above {st.main_reference}", f"below {st.main_reference}"
    if a in S.index:
        na, nb = int(S.loc[a, "nominal"]), int(S.loc[b, "nominal"])
        s += f"{na} signature{'s' if na != 1 else ''} lay above it and {nb} below it at P < 0.05."
    return s


def write_report(out_dir, benchmark: Optional[dict] = None, single_arm: Optional[dict] = None, concordance: Optional[dict] = None,
                 published: Optional[dict] = None, design: Optional[pd.DataFrame] = None, design_marks: Optional[Dict[str, float]] = None,
                 title: str = "barriernet report", note: str = "") -> Dict[str, pathlib.Path]:
    """Write report.html, report.pdf and each figure as PNG and PDF under ``out_dir``."""
    from matplotlib.backends.backend_pdf import PdfPages
    from . import __version__
    out = pathlib.Path(out_dir); (out / "figures").mkdir(parents=True, exist_ok=True)
    parts = [f"<h1>{html.escape(title)}</h1><p class='meta'>barriernet {__version__} · {_dt.date.today().isoformat()}</p>"]
    if note:
        parts.append(f"<p class='note'>{html.escape(note)}</p>")
    figs = []
    if benchmark is not None:
        parts.append("<h2>Published signatures in randomized comparisons</h2><p>" + html.escape(benchmark_text(benchmark)) + "</p>")
        for nm, f in (("benchmark", fig_benchmark(benchmark)), ("paired", fig_paired(benchmark)), ("reference", fig_reference(benchmark)), ("counts", fig_counts(benchmark))):
            figs.append((nm, f)); parts.append(f"<img alt='{nm}' src='data:image/png;base64,{_png(f)}'>")
        P = benchmark["pooled"]; cols = [c for c in ("score", "kind", "family", "trials_scored", "est", "lo", "hi", "p", "q", "emp_p_declared", "q_reference") if c in P]
        parts.append("<p class='note'>Pooled interactions (positive: more checkpoint benefit per SD; q: Benjamini-Hochberg over the signatures).</p>" + _table(P[cols].sort_values("est", ascending=False), 200))
    if single_arm is not None:
        PS = single_arm["pooled"]
        parts.append(f"<h2>Single-arm cohorts</h2><p>{len(PS)} scores were pooled over the single-arm cohorts (DerSimonian-Laird, with the wider of "
                     "the Hartung-Knapp and DerSimonian-Laird intervals). A single-arm association cannot separate a predictive score from a prognostic one.</p>")
        if concordance is not None and concordance["n"] > 2:
            f = fig_concordance(concordance); figs.append(("concordance", f)); parts.append(f"<img alt='concordance' src='data:image/png;base64,{_png(f)}'>")
        parts.append(_table(PS.sort_values("est", ascending=False)[["score", "k", "est", "lo", "hi", "p", "tau2"]], 200))
    if published is not None and "primary" in published:
        m = published["primary"]
        parts.append(f"<h2>Published subgroup results</h2><p>Over {m['k']} trials the benefit-oriented log ratio of hazard ratios was {m['est']:.3f} "
                     f"({m['lo']:.3f} to {m['hi']:.3f}), P = {m['p']:.4f}; τ² = {m['tau2']:.4f}, I² = {m['I2']:.0f}% (REML with the Knapp-Hartung adjustment).</p>")
        f = fig_published(published); figs.append(("published", f)); parts.append(f"<img alt='published' src='data:image/png;base64,{_png(f)}'>")
        rows = [dict(analysis=k, k=v.get("k"), est=v.get("est"), lo=v.get("lo"), hi=v.get("hi"), p=v.get("p")) for k, v in published.items() if isinstance(v, dict)]
        parts.append(_table(pd.DataFrame(rows)))
    if design is not None:
        parts.append("<h2>Trial design</h2><p>Events needed to detect a given difference between two scores' interactions, by their correlation: "
                     "2(1 − r)(z<sub>α/2</sub> + z<sub>β</sub>)² / (δ² p(1 − p)).</p>")
        f = fig_design(design, design_marks); figs.append(("design", f)); parts.append(f"<img alt='design' src='data:image/png;base64,{_png(f)}'>")
    page = f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>{html.escape(title)}</title><style>{CSS}</style></head><body><main>{''.join(parts)}</main></body></html>"
    (out / "report.html").write_text(page, encoding="utf-8")
    with PdfPages(out / "report.pdf") as pdf:
        for nm, f in figs:
            pdf.savefig(f, bbox_inches="tight"); f.savefig(out / "figures" / f"{nm}.pdf", bbox_inches="tight"); f.savefig(out / "figures" / f"{nm}.png", dpi=300, bbox_inches="tight")
    import matplotlib.pyplot as plt
    for _, f in figs:
        plt.close(f)
    return dict(html=out / "report.html", pdf=out / "report.pdf", figures=out / "figures")
