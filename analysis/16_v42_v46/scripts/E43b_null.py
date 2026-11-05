"""E43b (post hoc addendum) - bias of the cross-validated C-for-benefit, and a null-centered apparent estimate.

Protocol: analysis/E43b_prereg.md (hash checked). Code: the setup, models, cross-validation and evaluation of
E43_benefit.py, executed unchanged (step-0 and coverage tables written under E43b names). Within-arm joint permutations
of the four standardized scores: step 1, 20 per comparison through the E43 cross-validation (seeds [4043, k, j]);
step 2, 200 per comparison with full-data fits (seeds [5043, k, j]).
"""
import hashlib, sys
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_hb = hashlib.sha256(open(f"{A}/E43b_prereg.md", "rb").read()).hexdigest()
assert _hb == open(f"{A}/E43b_prereg.sha256").read().split()[0], "E43b protocol hash mismatch"
_src = open(f"{A}/E43_benefit.py").read()
_s0, _r = _src.split("# ------------------------------------------------------------------ cross-validated predictions")
_s1, _r = _r.split("# ------------------------------------------------------------------ evaluation")
_s1 = _s1.split("NP = int(sys.argv")[0]; _s2 = _r.split("etasks =")[0]
_s0 = _s0.replace("E43_step0_reproduction", "E43b_step0_reproduction").replace("E43_coverage_check.csv", "E43b_coverage_check.csv")
for _s in (_s0, _s1, _s2):
    exec(compile(_s, f"{A}/E43_benefit.py", "exec"))
print("E43b protocol sha256 verified:", _hb[:12])
from scipy.stats import chi2

NP = int(sys.argv[1]) if len(sys.argv) > 1 else NPROC
SS0 = {k: dict(v) for k, v in SS.items()}


def perm_idx(k, seed):
    a = TR[k]["arm"]; rng = np.random.default_rng(seed); p = np.arange(len(a))
    for v in (0, 1):
        ix = np.where(a == v)[0]; p[ix] = ix[rng.permutation(len(ix))]
    return p


def with_perm(k, p, f):
    SS[k] = {s: SS0[k][s][p] for s in SS0[k]}
    try:
        return f()
    finally:
        SS[k] = SS0[k]

# ------------------------------------------------------------------ step 1: cross-validated C under the permutation null


def s1_task(args):
    k, m, j = args; n = len(TR[k]["arm"])
    def go():
        prs = [cv_task((k, m, r))[0] for r in range(NREP)]
        return evaluate(k, np.nanmean(np.vstack(prs), 0), np.arange(n))[:2]
    return with_perm(k, perm_idx(k, [4043, KX[k], j]), go)

t1 = [(k, m, j) for k in ALL6 for m in MODELS for j in range(20)]
with mp.get_context("fork").Pool(NP) as pool:
    r1 = pool.map(s1_task, t1, chunksize=1)
N1 = pd.DataFrame([dict(comparison=k, model=m, perm=j, C=c, gradient=g) for (k, m, j), (c, g) in zip(t1, r1)])
E43 = pd.read_csv(f"{O}/E43_per_comparison.csv").set_index(["comparison", "model"])
CVN = N1.groupby(["comparison", "model"], sort=False).agg(cv_null_C=("C", "mean"), cv_null_C_sd=("C", "std"),
                                                          cv_null_g=("gradient", "mean")).reset_index()
CVN["cv_C_E43"] = [E43.loc[(k, m), "C"] for k, m in zip(CVN.comparison, CVN.model)]
CVN["cv_g_E43"] = [E43.loc[(k, m), "gradient"] for k, m in zip(CVN.comparison, CVN.model)]
CVN.to_csv(f"{O}/E43b_cv_null.csv", index=False)

# ------------------------------------------------------------------ step 2: apparent estimate, null-centered
APP = {}
for k in ALL6:
    n = len(TR[k]["arm"]); allp = np.arange(n)
    for m in MODELS:
        pred = fit_predict(k, m, allp, allp)[0]
        APP[(k, m)] = (evaluate(k, pred, allp), np.array([evaluate(k, pred, ix)[:2] for ix in BOOTI[k]]), np.std(pred))


def s2_task(args):
    k, m, j = args; n = len(TR[k]["arm"]); allp = np.arange(n)
    def go():
        try:
            return evaluate(k, fit_predict(k, m, allp, allp)[0], allp)[:2]
        except Exception:
            return (np.nan, np.nan)
    return with_perm(k, perm_idx(k, [5043, KX[k], j]), go)

t2 = [(k, m, j) for k in ALL6 for m in MODELS for j in range(200)]
with mp.get_context("fork").Pool(NP) as pool:
    r2 = pool.map(s2_task, t2, chunksize=8)
N2 = pd.DataFrame([dict(comparison=k, model=m, perm=j, C=c, gradient=g) for (k, m, j), (c, g) in zip(t2, r2)])
rows = []
for (k, m), (o, b, sd) in APP.items():
    nn = N2[(N2.comparison == k) & (N2.model == m)]
    rows.append(dict(comparison=k, model=m, kind=TR[k]["kind"], C_app=o[0], null_C=nn.C.mean(), null_C_p95=nn.C.quantile(0.95),
                     C_ctr=o[0] - nn.C.mean(), C_se=np.nanstd(b[:, 0], ddof=1), emp_p=(1 + np.sum(nn.C >= o[0])) / (1 + nn.C.notna().sum()),
                     g_app=o[1], null_g=nn.gradient.mean(), g_ctr=o[1] - nn.gradient.mean(), g_se=np.nanstd(b[:, 1], ddof=1),
                     pred_sd=sd, null_failed=int(nn.C.isna().sum())))
AP = pd.DataFrame(rows); AP.to_csv(f"{O}/E43b_apparent.csv", index=False)
CTR = {(r.comparison, r.model): (r.null_C, r.null_g) for r in AP.itertuples()}


def contrast(m1, m0, j, keys):
    dd, vv = [], []
    for k in keys:
        (o1, b1, _), (o0, b0, _) = APP[(k, m1)], APP[(k, m0)]; db = b1[:, j] - b0[:, j]; ok = np.isfinite(db)
        dd.append((o1[j] - CTR[(k, m1)][j]) - (o0[j] - CTR[(k, m0)][j])); vv.append(np.var(db[ok], ddof=1))
    dd, vv = np.array(dd), np.array(vv); w = 1 / vv; D_ = np.sum(w * dd) / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(diff=D_, se=se, lo=D_ - 1.96 * se, hi=D_ + 1.96 * se, p=2 * norm.sf(abs(D_ / se)), mdd=2.80 * se, k=len(keys))

SETS = [("six", ALL6), ("survival three", RENAL), ("pCR three", BREAST)]
POOL, CON = [], []
for lab, keys in SETS:
    for m in MODELS:
        s = AP[(AP.model == m) & AP.comparison.isin(keys)]
        c = pool_rule(s.C_ctr.values, s.C_se.values); g = pool_rule(s.g_ctr.values, s.g_se.values)
        POOL.append(dict(set=lab, model=m, C_ctr=c["est"], C_lo=c["lo"], C_hi=c["hi"], C_p=c["p"], g_ctr=g["est"], g_lo=g["lo"],
                         g_hi=g["hi"], g_p=g["p"], k=c["k"]))
        if m != "M_CD8A":
            for j, metric in ((0, "C-for-benefit"), (1, "gradient")):
                CON.append(dict(set=lab, contrast=f"{m} - M_CD8A", metric=metric, **contrast(m, "M_CD8A", j, keys)))
POOL, CON = pd.DataFrame(POOL), pd.DataFrame(CON)
HET = []
for metric in ("C-for-benefit", "gradient"):
    a = CON[(CON.set == "survival three") & (CON.contrast == "M_2AX - M_CD8A") & (CON.metric == metric)].iloc[0]
    b = CON[(CON.set == "pCR three") & (CON.contrast == "M_2AX - M_CD8A") & (CON.metric == metric)].iloc[0]
    q = (a["diff"] - b["diff"]) ** 2 / (a.se ** 2 + b.se ** 2); HET.append(dict(contrast="M_2AX - M_CD8A", metric=metric, Q=q, p=chi2.sf(q, 1)))
HET = pd.DataFrame(HET)
POOL.to_csv(f"{O}/E43b_pooled.csv", index=False); CON.to_csv(f"{O}/E43b_contrasts.csv", index=False)
HET.to_csv(f"{O}/E43b_heterogeneity.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print("\n== step 1: cross-validated C and gradient under the permutation null (20 permutations)\n", CVN.round(3).to_string(index=False))
print("\n== step 2: apparent, null-centered (200 permutations)\n", AP.round(3).to_string(index=False))
print("\n== pooled (null-centered)\n", POOL.round(3).to_string(index=False))
print("\n== contrasts (null-centered)\n", CON.round(3).to_string(index=False))
print("\n== survival three vs pCR three, primary contrast\n", HET.round(4).to_string(index=False))
