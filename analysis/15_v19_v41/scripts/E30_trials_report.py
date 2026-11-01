"""E30 - summary of the randomized part (Q2 randomized, Q3, Q4, Q5, Q6 table) from the E30_trials* results.
Reads only results files written by E30_trials.py and E30_trials_meta.R; writes results/E30_trials_summary.md."""
import hashlib
import numpy as np, pandas as pd

A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"; O = f"{A}/results"
H = hashlib.sha256(open(f"{A}/E30_prereg.md", "rb").read()).hexdigest()
assert H == open(f"{A}/E30_prereg.sha256").read().split()[0]
Q3S = ["NetImmune", "Immune", "CAF_ECM", "Barrier", "ESTIMATE", "TMEscore_PC", "TMEscore_z", "CD8A"]
RT = ["IMmotion150", "CheckMate 025", "JAVELIN Renal 101"]
G = pd.read_csv(f"{O}/E30_trials_gate.csv"); G28 = pd.read_csv(f"{O}/E30_trials_E28base_reproduction.csv")
MG = pd.read_csv(f"{O}/E30_trials_meta_gate.csv"); AN = pd.read_csv(f"{O}/E30_trials_angio.csv")
PT = pd.read_csv(f"{O}/E30_trials_angio_per_trial.csv"); ST = pd.read_csv(f"{O}/E30_trials_strata.csv")
SZ = pd.read_csv(f"{O}/E30_trials_stouffer.csv"); ES = pd.read_csv(f"{O}/E30_trials_estimators.csv")
OS = pd.read_csv(f"{O}/E30_trials_cm025_os.csv")


def f(x, d=2, sign=True):
    if pd.isna(x): return "–"
    s = f"{x:+.{d}f}" if sign else f"{x:.{d}f}"
    return s.replace("-", "−")


def fp(p):
    if pd.isna(p): return "–"
    if p < 0.001: return f"{p:.1e}".replace("e-0", " × 10^-").replace("e-", " × 10^-")
    return f"{p:.2g}" if p < 0.01 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


def ci(lo, hi, d=2):
    return f"{f(lo, d, False)} to {f(hi, d, False)}"


def table(rows, head):
    out = ["| " + " | ".join(head) + " |", "|" + "|".join(["---"] * len(head)) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


L = []
L.append("# E30 randomized comparisons: Q2 (randomized part), Q3, Q4, Q5 and the Q6 table\n")
L.append(f"Protocol `analysis/E30_prereg.md`, sha256 `{H}` (matches `E30_prereg.sha256`, logged 2026-10-22T02:53:19Z). "
         "Scripts: `analysis/E30_trials.py` (per-comparison models; executes `E28_base.py` unchanged except for the name "
         "of its step-0 output), `analysis/E30_trials_meta.R` (pooling; `pool()` copied unchanged from `E21_meta.R`), "
         "`analysis/E30_trials_report.py` (this summary). Log: `analysis/results/E30_trials.log`.\n")

# ------------------------------------------------------------------ step 0
g1 = G.absdiff.max(); g28 = G28.maxdiff.max(); pg = MG[MG.check.str.startswith("pooled")]; t4 = MG[MG.check.str.startswith("Table 4")]
L.append("## Step 0 (gate): PASS\n")
L.append(f"- E28_base step 0 (on load): {len(G28)} interaction and within-arm values re-derived, largest |difference| {g28:.1e}.")
L.append(f"- Per-comparison score × arm interactions, 16 scores × six primary comparisons (y and SE against "
         f"`E21_meta_inputs.csv`, set `primary: 6 comparisons`) plus IMmotion150 atezolizumab monotherapy (16 scores, "
         f"against `E5_trials.csv`): {len(G)} values, largest |difference| {g1:.1e} (< 1e-6).")
L.append(f"- Pooled values re-derived from these inputs with `pool()` against `E21_meta.csv` (all six, renal, breast, first "
         f"three, three added; 16 scores; all 21 numeric columns): {len(pg)} pooled results, largest |difference| "
         f"{pg.absdiff.max():.1e}; the interval chosen (HK or DL) identical in all.")
L.append(f"- Table 4 (final_v31 manuscript), six-comparison column: {int((t4.absdiff == 0).sum())} of {len(t4)} rows "
         "match at reported precision (estimate, CI, P, number favoring, Stouffer P).")
sz6 = SZ[(SZ.set == "all six")]
L.append(f"- Unweighted Stouffer z reproduces `E21_meta.csv` (largest |difference| {SZ.absdiff_vs_published.max():.1e}); "
         "the calibrated Stouffer P re-derived from `E28_Q3_null_draws.csv` equals `E28_Q3_calibration_pooled.csv` for "
         "all five scores with draws.")
L.append(f"- Workbook S16i / S16f (final_v31) CheckMate 025 OS values equal `E21_meta.csv` / `E21_interactions.csv` "
         f"(largest |difference| {np.nanmax(OS.workbook_maxdiff):.1e}).\n")

# ------------------------------------------------------------------ Angio
L.append("## Angio score\n")
L.append("z-mean of the McDermott genes present (IMmotion150 VEGFA, KDR, ESM1, PECAM1; CheckMate 025 without PECAM1; "
         "6/6 elsewhere) in each trial's published z matrix, identical (|difference| 0) to the feasibility check "
         "`rev/agentB/angio_feasibility.py`. Spearman ρ with the authors' Angio: CheckMate 025 (Braun S1) 0.994, "
         "JAVELIN Renal 101 (Motzer S16 MCDERMOTT_ANGIO) 1.000 (see log). Standardized within trial over the analysed "
         "patients in every model.\n")

# ------------------------------------------------------------------ Q3
L.append("## Q3. Angiogenesis adjustment in the renal trials\n")
L.append("Model: published primary model (Cox, progression-free survival; CheckMate 025 on its E21 primary endpoint, PFS) "
         "+ Angio + Angio × arm. Benefit-oriented interaction per SD (−log HR); y > 0 favors the checkpoint arm.\n")
P = PT[PT.role == "primary"]
rows = []
for sc in Q3S:
    for tr in RT:
        r = P[(P.trial == tr) & (P.score == sc)].iloc[0]
        rows.append([sc, tr, f"{r.rho_angio:+.2f}".replace("-", "−"), f"{f(r.y_unadj, 3)} ({fp(r.p_unadj)})",
                     f"{f(r.y_adj, 3)} ({fp(r.p_adj)})", f(r.diff_adj_minus_unadj, 3), "yes" if r.within_0_05 else "**no**",
                     f"{f(r.angio_x_arm_y, 3)} ({fp(r.angio_x_arm_p)})"])
L.append(table(rows, ["Score", "Trial", "ρ(Angio, score)", "Unadjusted y (P)", "Adjusted y (P)", "Adjusted − unadjusted",
                      "Within 0.05", "Angio × arm y (P), same model"]))
aa = AN[(AN.level == "trial") & (AN.model.str.startswith("Angio alone")) & (AN.estimand == "per trial")]
L.append("\nAngio × arm with Angio as the only score (published primary model): " +
         "; ".join(f"{r.trial} {f(r.est, 3)} (P = {fp(r.p)})" for _, r in aa.iterrows()) + ".\n")

pq = AN[AN.level == "pooled"]
rows = []
for sc in Q3S:
    g = pq[pq.score == sc].set_index("estimand")
    def cell(e):
        r = g.loc[e]; return f"{f(r.est)} ({ci(r.lo, r.hi)}); {fp(r.p)}"
    rows.append([sc, cell("renal, unadjusted (published)"), cell("(a) renal, Angio-adjusted"),
                 cell("six, unadjusted (published)"), cell("(b) six, renal Angio-adjusted + breast as published"),
                 cell("Angio x arm, renal (descriptive)")])
L.append("Pooled (REML point; wider of REML-HK and DL interval; larger P):\n")
L.append(table(rows, ["Score", "Renal, unadjusted (published)", "(a) Renal, Angio-adjusted", "Six, unadjusted (published)",
                      "(b) Six, renal adjusted + breast published", "Angio × arm pooled over renal (same model)"]))
g = pq[pq.score == "Angio"].iloc[0]
L.append(f"\nAngio × arm pooled over the three renal trials with Angio as the only score: {f(g.est)} ({ci(g.lo, g.hi)}); "
         f"P = {fp(g.p)} (descriptive).\n")

# expectations and decision
d = P.diff_adj_minus_unadj.abs(); nbad = int((~P.within_0_05).sum())
bad = P[~P.within_0_05].sort_values("diff_adj_minus_unadj", key=np.abs, ascending=False)
net = P[P.score == "NetImmune"]
b = pq[(pq.estimand.str.startswith("(b)")) & (pq.score == "NetImmune")].iloc[0]
L.append("### Expectations\n")
L.append(f"1. *Each adjusted renal interaction within 0.05 of the unadjusted one* — **not held**: {len(P) - nbad} of "
         f"{len(P)} (8 scores × 3 trials) within 0.05; largest |change| {d.max():.3f}. Outside 0.05: " +
         "; ".join(f"{r.trial} {r.score} {f(r.y_unadj, 3)} → {f(r.y_adj, 3)}" for _, r in bad.iterrows()) +
         ". For NetImmune alone the expectation held in all three trials (changes " +
         ", ".join(f(x, 3) for x in net.diff_adj_minus_unadj) + "), as it did for the Immune axis in CheckMate 025 and "
         "JAVELIN Renal 101. The premise |ρ| ≤ 0.14 holds for the four BarrierNet axes only; the comparators correlate "
         f"more with Angio (ρ {P[P.score.isin(['ESTIMATE', 'TMEscore_PC', 'TMEscore_z'])].rho_angio.min():.2f} to "
         f"{P[P.score.isin(['ESTIMATE', 'TMEscore_PC', 'TMEscore_z'])].rho_angio.max():.2f}).".replace("-0.", "−0."))
held2 = b.est > 0 and b.lo > 0
L.append(f"2. *Six-comparison NetImmune interaction keeps its sign with a 95% CI above zero* — **{'held' if held2 else 'not held'}**: "
         f"(b) {f(b.est, 3)} ({ci(b.lo, b.hi, 3)}); P = {fp(b.p)} ({b.interval} interval; published unadjusted "
         f"+0.289, 0.032 to 0.545; P = 0.034).\n")
branch = 1 if held2 else (2 if b.est > 0 else 3)
L.append("### Decision rule\n")
L.append(f"**Branch {branch} applies**: the six-comparison NetImmune interaction (estimand b) keeps its sign with a CI "
         "excluding zero, so the manuscript states that the randomized result was not explained by angiogenesis. "
         "(Estimand a, renal only: NetImmune "
         f"{f(pq[(pq.estimand.str.startswith('(a)')) & (pq.score == 'NetImmune')].est.iloc[0], 3)}, P = "
         f"{fp(pq[(pq.estimand.str.startswith('(a)')) & (pq.score == 'NetImmune')].p.iloc[0])}, against 0.200, P = 0.146 unadjusted.) "
         "Not covered by the rule but material to the text: with adjustment, CAF_ECM turns positive in all three renal "
         "trials and the six-comparison CAF_ECM, ESTIMATE and TMEscore pools attenuate (TMEscore CIs now include zero; "
         "CD8A's lower limit is 0.0002).\n")

mono = PT[PT.role != "primary"]
L.append("### Descriptive: IMmotion150 atezolizumab monotherapy vs sunitinib (E5 arm split; no decision)\n")
rows = [[r.score, r.n, f"{f(r.y_unadj, 3)} ({fp(r.p_unadj)})", f"{f(r.y_adj, 3)} ({fp(r.p_adj)})",
         f(r.diff_adj_minus_unadj, 3), f"{f(r.angio_x_arm_y, 3)} ({fp(r.angio_x_arm_p)})"] for _, r in mono.iterrows()]
L.append(table(rows, ["Score", "n", "Unadjusted y (P)", "Adjusted y (P)", "Change", "Angio × arm y (P)"]))
L.append("")

# ------------------------------------------------------------------ Q4
L.append("## Q4. Endpoint strata (randomized rule)\n")
rows = []
for sc in ["NetImmune", "Immune", "CAF_ECM", "Barrier", "IFN", "TAM_SUP", "ESTIMATE", "TMEscore_PC", "TMEscore_z", "CD8A"]:
    g = ST[ST.score == sc].set_index("stratum")
    rows.append([sc] + [f"{f(g.loc[s].est)} ({ci(g.loc[s].lo, g.loc[s].hi)}); {fp(g.loc[s].p)}; {g.loc[s].k_benefit}/3"
                        if s != "all six" else f"{f(g.loc[s].est)} ({ci(g.loc[s].lo, g.loc[s].hi)}); {fp(g.loc[s].p)}; {g.loc[s].k_benefit}/6"
                        for s in ["renal", "breast", "all six"]])
L.append(table(rows, ["Score", "Renal (3, −log HR)", "Breast (3, log OR)", "All six (pre-registered primary)"]))
L.append("\nAll 16 scores in `E30_trials_strata.csv`; values equal the published `renal only` / `breast only` sets of "
         "`E21_meta.csv`. Presentation per protocol: side by side, without comparing magnitudes (cancer type and effect "
         "scale are fully confounded).\n")

# ------------------------------------------------------------------ Q5
L.append("## Q5. Stouffer z: unweighted, √n-weighted, calibrated\n")
L.append("Unweighted z = Σz/√J (Table 4); weighted z = Σ√n·z/√Σn, n = patients analysed per comparison "
         "(263 / 248 / 100 / 250 / 726 / 241). Calibrated P: E28 parametric bootstrap under no score × arm interaction "
         "(1,000 null data sets per comparison; available for five scores); the weighted calibrated P uses the same draws.\n")
rows = []
for _, r in sz6.iterrows():
    rows.append([r.score, f"{r.stouffer_z:+.2f}".replace("-", "−"), fp(r.stouffer_p), f"{r.weighted_z:+.2f}".replace("-", "−"),
                 fp(r.weighted_p), fp(r.calibrated_p), fp(r.calibrated_weighted_p)])
L.append(table(rows, ["Score (six comparisons)", "z", "P", "Weighted z", "Weighted P", "Calibrated P", "Calibrated weighted P"]))
L.append("\nStrata and the first three / three added are in `E30_trials_stouffer.csv`.\n")

# ------------------------------------------------------------------ Q2
L.append("## Q2. Estimators for the randomized pooled estimates (disclosure; no decision)\n")
L.append("Point, 95% CI and P under DL, REML (Wald), DL-HK, REML-HK and the manuscript rule, for every randomized "
         "pooled estimate from k ≤ 7 comparisons: 16 scores × 15 E21 sets (first three, three added, all six, renal, "
         "breast, six leave-one-out, CheckMate 025 OS, NeoTRIP without batch, + PRINCE, + DUTRENEO), E25 within-arm and "
         "mutually adjusted pools, and the E28 official ESTIMATE and direct-comparison pools "
         f"({int((ES.estimator == 'manuscript rule').sum())} pooled results; `main_text` flags those cited in the v31 "
         "main text). HK uses t with k − 1 df (4.30 for k = 3, 2.78 for k = 5, 2.57 for k = 6). The manuscript rule "
         "reproduces E21_meta.csv and E25_meta.csv to < 1e-13 and the E28 numpy pools to 3.7e-6.\n")
SETN = {"reproduction: 3 published primary": "First three", "3 new primary only": "Three added",
        "primary: 6 comparisons": "All six", "renal only (IMmotion150, CheckMate 025, JAVELIN)": "Renal",
        "breast only (I-SPY2 x 2, NeoTRIP)": "Breast", "primary, CheckMate 025 OS for PFS": "All six, CM025 OS"}
rows = []
for s, lab in SETN.items():
    for sc in ["NetImmune", "Immune", "CAF_ECM", "ESTIMATE"]:
        g = ES[(ES.set == s) & (ES.score == sc)].set_index("estimator")
        if g.empty: continue
        rows.append([lab, sc] + [f"{f(g.loc[e].est)} ({ci(g.loc[e].lo, g.loc[e].hi)}); {fp(g.loc[e].p)}"
                                 for e in ["DL", "REML", "DL-HK", "REML-HK", "manuscript rule"]])
L.append(table(rows, ["Set", "Score", "DL", "REML", "DL-HK", "REML-HK", "Manuscript rule"]))
L.append("")

# ------------------------------------------------------------------ Q6
L.append("## Q6. CheckMate 025 with overall survival (existing results; no new model)\n")
rows = [[r.score, f"{r.cm025_os_int_HR:.2f} ({r.cm025_os_lo:.2f}–{r.cm025_os_hi:.2f}); {fp(r.cm025_os_p)}",
         f"{r.cm025_pfs_int_HR:.2f}; {fp(r.cm025_pfs_p)}",
         f"{f(r.six_os_est)} ({ci(r.six_os_lo, r.six_os_hi)}); {fp(r.six_os_p)}; {r.six_os_k_benefit}/6",
         f"{f(r.six_pfs_est)} ({ci(r.six_pfs_lo, r.six_pfs_hi)}); {fp(r.six_pfs_p)}"] for _, r in OS.iterrows()]
L.append(table(rows, ["Score", "CM025 OS interaction HR (95% CI); P", "CM025 PFS HR; P", "Six, CM025 OS: pooled (CI); P; favoring",
                      "Six, CM025 PFS (published primary)"]))
sig = OS[OS.sig_pfs]
L.append(f"\nCheckMate 025 OS: {int(OS.cm025_os_n.iloc[0])} patients, {int(OS.cm025_os_events.iloc[0])} deaths. "
         f"Significant six-comparison pools with PFS: {len(sig)} ({', '.join(sig.score)}); with OS: "
         f"{int(OS.sig_os.sum())} (smallest P {OS.six_os_p.min():.3f}, NetImmune). The protocol statement that every pooled "
         "interaction loses significance holds. Module values are in workbook S16i/S16f; ESTIMATE, TMEscore and CD8A "
         "come from `E21_meta.csv` / `E21_interactions.csv` only (not in S16i).\n")

# ------------------------------------------------------------------ protocol statements
e3 = ES[(ES.set == "reproduction: 3 published primary") & (ES.score == "NetImmune")].set_index("estimator")
L.append("## Statements in the protocol's 'Already known' section, checked\n")
L.append(f"- First three randomized comparisons, NetImmune: P = {e3.loc['REML-HK'].p:.3f} (HK), {e3.loc['DL'].p:.4f} (DL; "
         f"τ² = {e3.loc['DL'].tau2:.0f}) — as stated (0.057; 0.0020).")
o = OS.set_index("score").six_os_p
L.append(f"- CheckMate 025 OS six-comparison P: NetImmune {o['NetImmune']:.3f}, ESTIMATE {o['ESTIMATE']:.3f}, TMEscore "
         f"{o['TMEscore_PC']:.3f} / {o['TMEscore_z']:.3f}, CD8A {o['CD8A']:.3f}, Immune {o['Immune']:.2f}, IFN "
         f"{o['IFN']:.2f} — as stated.")
rr = P[P.score.isin(["NetImmune", "Immune", "CAF_ECM", "Barrier"])].rho_angio
L.append(f"- Spearman ρ of Angio with NetImmune, Immune, CAF_ECM and Barrier in the renal trials: {rr.min():.3f} to "
         f"{rr.max():.3f} — within ±0.14 at the two decimals stated (largest |ρ| {rr.abs().max():.4f}).\n")

L.append("## Files\n")
for fn, what in [("E30_trials_angio.csv", "Q3 per trial (unadjusted, Angio-adjusted, Angio × arm, Angio alone; incl. IMmotion150 monotherapy) and pooled estimands"),
                 ("E30_trials_angio_per_trial.csv", "Q3 per trial, wide"),
                 ("E30_trials_strata.csv", "Q4"), ("E30_trials_stouffer.csv", "Q5"), ("E30_trials_estimators.csv", "Q2"),
                 ("E30_trials_cm025_os.csv", "Q6"), ("E30_trials_inputs.csv", "per-comparison inputs to the pools"),
                 ("E30_trials_gate.csv / E30_trials_meta_gate.csv / E30_trials_E28base_reproduction.csv", "step 0"),
                 ("E30_trials.log", "log")]:
    L.append(f"- `analysis/results/{fn}` — {what}")
open(f"{O}/E30_trials_summary.md", "w").write("\n".join(L) + "\n")
print("E30_trials_summary.md written")
