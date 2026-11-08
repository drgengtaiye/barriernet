"""Check that barriernet 0.4.0 reproduces the manuscript's analyses from the same inputs. No new analysis: every value
is compared with a result file already written by the analysis code.

- benchmark: E32 setup (E32_trials.py up to Q1, as in E34/E37/E40), the six comparisons, the 51 distinct library
  signatures, seeds 2032 (random reference) and 3032 (patient resamples, as E32 and E46). Against E32_primary.csv
  (with Benjamini-Hochberg counts recomputed by scipy) and E46_A_paired.csv.
- single-arm: results/E40_single_arm_genes.csv.gz minus the three twin-flagged samples, as in E40_concordance.py.
  Against E40_single_arm_per_cohort.csv and E40_concordance.csv.
- published trials: the bundled example table (E41_adjudicated.csv, eligible rows). Against E41_meta.csv.
- design: against E39_required.csv.
Writes results/barriernet04_check.csv."""
import sys, time
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_s32 = open(f"{A}/E32_trials.py").read().split("# ------------------------------------------------------------------ Q1:")[0]
_s32 = _s32.replace("E32_step0_reproduction", "bn04_step0_reproduction").replace("E32_coverage.csv", "bn04_coverage_check.csv")
exec(compile(_s32, f"{A}/E32_trials.py", "exec"))
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/repo_barriernet/barriernet")
import barriernet
from barriernet.calibrate import Trial
from barriernet.benchmark import benchmark
from barriernet.singlearm import Cohort, single_arm, concordance
from barriernet.published import example_table, pool_published
from barriernet.design import events_needed
from barriernet.library import load_library, signature_weights
CHK = []


def chk(item, pkg, art, tol, note=""):
    d = abs(float(pkg) - float(art)) if np.isfinite(pkg) and np.isfinite(art) else np.nan
    CHK.append(dict(item=item, package=float(pkg), article=float(art), abs_diff=d, tolerance=tol, ok=bool(d <= tol) if np.isfinite(d) else False, note=note))


# ------------------------------------------------------------------ benchmark
trials = [Trial(k, TR[k]["Z"], TR[k]["arm"], time=TR[k].get("t") if TR[k]["kind"] == "cox" else None,
                event=TR[k].get("e") if TR[k]["kind"] == "cox" else None, response=TR[k].get("y") if TR[k]["kind"] != "cox" else None,
                covariates=TR[k].get("cov"), firth=TR[k].get("firth", True)) for k in ALL6]
t0 = time.time()
B = benchmark(trials, n_random=1000, n_boot=300, seed=2032, boot_seed=3032, n_jobs=NPROC)
print(f"benchmark {time.time() - t0:.0f} s")
O_ = f"{A}/results"
E = pd.read_csv(f"{O_}/E32_primary.csv"); E = E[(E.set == "six (primary)") & (E.kind == "library") & (E.signature != "TIS_Damotte")].set_index("signature")
P = B["pooled"].set_index("score"); S = P[P.kind == "signature"]
j = S.join(E[["est", "p", "k", "emp_p_declared", "ref_median"]], rsuffix="_E32", how="inner")
chk("signatures pooled (of 51)", S.est.notna().sum(), E.est.notna().sum(), 0)
chk("pooled interaction, max |difference| over 51", np.abs(j.est - j.est_E32).max(), 0, 1e-5, "package Cox solver against lifelines")
chk("pooled interaction P, max |difference|", np.abs(j.p - j.p_E32).max(), 0, 1e-4)
chk("comparisons scored, mismatches", int((j.trials_scored != j.k).sum()), 0, 0)
chk("empirical P < 0.05 classification, disagreements", int(((j.emp_p_declared < 0.05) != (j.emp_p_declared_E32 < 0.05)).sum()), 0, 0,
    "random draws differ from E32 (seeding by trial subset)")
chk("empirical P, max |difference|", np.abs(j.emp_p_declared - j.emp_p_declared_E32).max(), 0, 0.05, "Monte Carlo error with 1,000 draws")
from scipy.stats import false_discovery_control as _bh
SM = B["summary"].set_index("test"); qa, qr = _bh(E.p.values), _bh(E.emp_p_declared.values)
chk("favor checkpoint arm, P < 0.05", SM.loc["favor the checkpoint arm", "nominal"], int(((E.est > 0) & (E.p < 0.05)).sum()), 0)
chk("favor control arm, P < 0.05", SM.loc["favor the control arm", "nominal"], int(((E.est < 0) & (E.p < 0.05)).sum()), 0)
chk("favor checkpoint arm, q < 0.05", SM.loc["favor the checkpoint arm", "fdr"], int(((E.est > 0) & (qa < 0.05)).sum()), 0)
chk("beyond size-matched random genes, P < 0.05", SM.loc["beyond size-matched random genes", "nominal"], int((E.emp_p_declared < 0.05).sum()), 0)
chk("beyond size-matched random genes, q < 0.05", SM.loc["beyond size-matched random genes", "fdr"], int((qr < 0.05).sum()), 0)
AP = pd.read_csv(f"{O_}/E46_A_paired.csv"); PR = B["paired"]
for r, nm in {"R3": "T-cell mean (CD8A, PRF1, GZMA)", "R13": "T-cell mean (13 genes)"}.items():
    a = AP[(AP.set == "six (primary)") & (AP.reference == r) & (AP.score != "TIS_Damotte")].set_index("score"); b = PR[PR.reference == nm].set_index("score")
    jj = b.join(a[["diff", "p", "above", "below"]], rsuffix="_E46", how="inner")
    chk(f"{r}: signatures above, P < 0.05", int(((jj["diff"] > 0) & (jj.p < 0.05)).sum()), int(jj.above.sum()), 0)
    chk(f"{r}: signatures below, P < 0.05", int(((jj["diff"] < 0) & (jj.p < 0.05)).sum()), int(jj.below.sum()), 0)
    ok = (jj["diff"] - jj.diff_E46).abs() < 1e-4
    chk(f"{r}: paired differences equal to 1e-4 (of {len(jj)})", int(ok.sum()), len(jj) - 1, 0,
        "CRMA_Shukla differs: in JAVELIN Renal 101 its score has 30 distinct values in 726 patients and a few resamples give near-degenerate Cox fits")
E2 = pd.read_csv(f"{O_}/E32_primary.csv"); E2 = E2[(E2.set == "six (primary)") & (E2.kind == "library") & (E2.signature != "TIS_Damotte")].set_index("signature")
b = PR[PR.reference == "CD8A"].set_index("score").join(E2[["diff_vs_CD8A", "cd8a_p"]], how="inner")
chk("CD8A: signatures below, P < 0.05", int(((b["diff"] < 0) & (b.p < 0.05)).sum()), int(((b.diff_vs_CD8A < 0) & (b.cd8a_p < 0.05)).sum()), 0)
chk("CD8A: signatures above, P < 0.05", int(((b["diff"] > 0) & (b.p < 0.05)).sum()), int(((b.diff_vs_CD8A > 0) & (b.cd8a_p < 0.05)).sum()), 0)
for nm, art in (("T-cell mean (13 genes)", 0.18), ("T-cell mean (CD8A, PRF1, GZMA)", 0.22)):
    chk(f"{nm}: pooled interaction (manuscript, 2 decimals)", round(P.loc[nm, "est"], 2), art, 0)

# ------------------------------------------------------------------ single-arm
SA = pd.read_csv(f"{O_}/E40_single_arm_genes.csv.gz"); fl = pd.read_csv(f"{O_}/E37_twin_flags.csv")
SA = SA[~SA.set_index(["cohort", "sample"]).index.isin(list(zip(fl.development, fl["sample"])))]
GC = [c for c in SA.columns if c not in ("cohort", "set", "sample", "response")]
cohorts = []
for c, g in SA.groupby("cohort", sort=False):
    Z = g[GC].T; Z.columns = g["sample"].values; Z = Z.dropna(how="all")
    cohorts.append(Cohort(c, Z, g.response.to_numpy(float)))
W = signature_weights(load_library(distinct=False))
SAR = single_arm(cohorts, W, min_cohorts=10)
PC = pd.read_csv(f"{O_}/E40_single_arm_per_cohort.csv")
m = SAR["per_cohort"].merge(PC, left_on=["score", "cohort"], right_on=["signature", "cohort"], suffixes=("", "_E40"))
chk("single-arm per-cohort rows matched", len(m), len(PC), 0)
chk("single-arm per-cohort log-odds, max |difference|", np.abs(m.logOR - m.logOR_E40).max(), 0, 1e-6)
C40 = pd.read_csv(f"{O_}/E40_concordance.csv").set_index("signature"); sp = SAR["pooled"].set_index("score").join(C40[["sa_est", "sa_p"]], how="inner")
chk("single-arm pooled, max |difference|", np.abs(sp.est - sp.sa_est).max(), 0, 1e-6)
chk("single-arm pooled P, max |difference|", np.abs(sp.p - sp.sa_p).max(), 0, 1e-6)
RP = C40[["int_est", "int_p"]].rename(columns={"int_est": "est", "int_p": "p"}).reset_index().rename(columns={"signature": "score"})
cc = concordance(SAR["pooled"], RP)
chk("Spearman, single-arm against randomized (52 rows)", round(cc["rho"], 2), 0.83, 0, "E40_outcome.md")

# ------------------------------------------------------------------ published trials
M = pd.read_csv(f"{O_}/E41_meta.csv").set_index("analysis"); pub = pool_published(example_table())
for key, lab in (("primary", "primary: survival, independent of the six"), ("per_sd", "primary, per SD"),
                 ("pdl1", "PD-L1 IHC, survival, same comparisons"), ("tmb", "TMB, survival, same comparisons"),
                 ("signature_where_pdl1", "signature, survival, comparisons that also report PD-L1")):
    for c in ("est", "lo", "hi", "p"):
        chk(f"published {key}: {c}", pub[key][c], M.loc[lab, c], 5e-5)
    chk(f"published {key}: k", pub[key]["k"], M.loc[lab, "k"], 0)
lo = M.loc["leave one out, survival (lo/hi = range of estimates)"]
chk("published leave one out: lowest", pub["leave_one_out"]["lo"], lo["lo"], 5e-5); chk("published leave one out: highest", pub["leave_one_out"]["hi"], lo["hi"], 5e-5)
chk("published Egger P", pub["egger"]["p"], M.loc["Egger test, survival", "p"], 5e-4)

# ------------------------------------------------------------------ design
D39 = pd.read_csv(f"{O_}/E39_required.csv")
for r in (0.5, 0.6, 0.7, 0.8, 0.9, 0.95):
    for _, row in D39.iterrows():
        chk(f"events, delta {row.delta}, r {r} (ratio to E39)", events_needed(row.delta, r) / row[f"events_rho{r}"], 1, 1e-6,
            "E39 used rounded normal quantiles")

C = pd.DataFrame(CHK); C.to_csv(f"{O_}/barriernet04_check.csv", index=False)
print("barriernet", barriernet.__version__); print(C[~C.item.str.startswith("events,")].to_string(index=False))
print(f"design rows: {int(C.item.str.startswith('events,').sum())}, all within tolerance: {bool(C[C.item.str.startswith('events,')].ok.all())}")
print(f"checks passed: {int(C.ok.sum())} of {len(C)}")
