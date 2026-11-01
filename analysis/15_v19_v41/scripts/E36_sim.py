"""E36 - operating characteristics of the calibration framework (semi-synthetic simulation).

Protocol: analysis/E36_prereg.md (sha256 in E36_prereg.sha256, logged 2026-10-31T02:59:12Z, before this script was
written). Setup: E34_trials.py up to its fitting stage (E32 data of the six primary comparisons), step-0 and coverage
tables written under E36 names. Real expression, arms and covariates; outcomes simulated under a known truth.
Survival: event time ~ Exp(rate exp(-eta)), censoring ~ Exp(rate (1 - f) / f), f = observed event fraction.
pCR: logit P = logit(observed control-arm rate) + eta. Analysis models as the manuscript; Cox by Newton-Raphson
(checked against lifelines on simulated data before use). The random reference uses 200 size-matched random gene sets
per comparison and size, drawn once (seeded) and shared by all replicates of a scenario.
"""
import hashlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h36 = hashlib.sha256(open(f"{A}/E36_prereg.md", "rb").read()).hexdigest()
assert _h36 == open(f"{A}/E36_prereg.sha256").read().split()[0], "E36 protocol hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E36_step0_reproduction").replace("E34_coverage_check.csv", "E36_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
print("E36 protocol sha256 verified:", _h36[:12])
import sys

NP = int(sys.argv[1]) if len(sys.argv) > 1 else 96
MODS = ["T_CELL", "IFN", "APM", "TLS_B", "CAF_ECM"]
MARK = ["CD8A"] + MODS + ["NetImmune"]
DIRN = {m: (-1 if m == "CAF_ECM" else 1) for m in MARK}


# ------------------------------------------------------------------ fast Cox (no ties) and the analysis models
def cox_fast(X, t, e, it=60, tol=1e-10):
    o = np.argsort(-t, kind="stable"); X = X[o]; e = e[o].astype(bool); p = X.shape[1]; b = np.zeros(p)

    def parts(b):
        eta = X @ b; c = eta.max(); w = np.exp(eta - c)
        S0 = np.cumsum(w); S1 = np.cumsum(w[:, None] * X, 0); S2 = np.cumsum(w[:, None, None] * (X[:, :, None] * X[:, None, :]), 0)
        ll = np.sum((eta - c - np.log(S0))[e])
        m1 = S1[e] / S0[e, None]
        g = np.sum(X[e] - m1, 0)
        H = np.sum(S2[e] / S0[e, None, None] - m1[:, :, None] * m1[:, None, :], 0)
        return ll, g, H
    ll, g, H = parts(b)
    for _ in range(it):
        step = np.linalg.solve(H, g); s_ = 1.0
        while True:
            ll2, g2, H2 = parts(b + s_ * step)
            if ll2 >= ll - 1e-12 or s_ < 1e-6:
                break
            s_ /= 2
        b = b + s_ * step; ll, g, H = ll2, g2, H2
        if np.max(np.abs(s_ * step)) < tol:
            break
    return b, np.sqrt(np.diag(np.linalg.inv(H)))


def fit_int(d, s, out):
    s = std(s); a = d["arm"]
    if d["kind"] == "cox":
        b, se = cox_fast(np.column_stack([s, a, s * a]), out[0], out[1]); return -b[2], se[2]
    X = np.column_stack([np.ones_like(s), s, a, s * a] + ([d["cov"]] if d["cov"] is not None else []))
    return lfit(X, out, [3], d["firth"])[0]


def fit_within(d, s, out):
    s = std(s); m = d["arm"] == 1
    if d["kind"] == "cox":
        b, se = cox_fast(s[m][:, None], out[0][m], out[1][m]); return -b[0], se[0]
    cov = d["cov"][m] if d["cov"] is not None else None
    X = np.column_stack([np.ones(m.sum()), s[m]] + ([cov] if cov is not None else []))
    return lfit(X, out[m], [1], d["firth"])[0]


def sfit(f, *args):
    """A non-converging or singular fit is recorded as missing (as E32 does for reference draws)."""
    try:
        return f(*args)
    except Exception:
        return (np.nan, np.nan)


def sub_out(d, out, idx):
    return (out[0][idx], out[1][idx]) if d["kind"] == "cox" else out[idx]


# ------------------------------------------------------------------ markers, truths, baseline rates
NC, PRESW = {}, {}
for k in ALL6:
    Z = TR[k]["Z"]; nc = set(nonconst(Z)); NC[k] = nc
    PRESW[(k, "CD8A")] = OBS[("CD8A", k)][1:]; PRESW[(k, "NetImmune")] = OBS[("NetImmune", k)][1:]
    for m in MODS:
        g = [x for x in MOD[m] if x in nc]; PRESW[(k, m)] = (g, np.full(len(g), 1.0))
SCO = {}
for k in ALL6:
    for m in MARK:
        g, w = PRESW[(k, m)]; SCO[(k, m)] = wscore(TR[k]["Z"], list(g), np.asarray(w, float))
        assert np.isfinite(SCO[(k, m)]).all(), (k, m)
UNIV = {}
for k in ALL6:
    Z = TR[k]["Z"]; g = [x for x in nonconst(Z) if np.isfinite(Z.loc[x].to_numpy(float)).all()]
    UNIV[k] = Z.loc[g].to_numpy(float)
PC1 = {}
for k in ALL6:
    U = UNIV[k] - UNIV[k].mean(1, keepdims=True); u, sv, vt = np.linalg.svd(U, full_matrices=False)
    pc = vt[0] * np.sign(np.corrcoef(vt[0], SCO[(k, "CD8A")])[0, 1]); PC1[k] = std(pc)
BASE = {}
for k in ALL6:
    d = TR[k]
    BASE[k] = d["e"].mean() if d["kind"] == "cox" else d["y"][d["arm"] == 0].mean()
print("baseline (event fraction or control pCR):", {k: round(v, 3) for k, v in BASE.items()})
print("corr(PC1, CD8A):", {k: round(np.corrcoef(PC1[k], SCO[(k, 'CD8A')])[0, 1], 2) for k in ALL6})


def simulate(k, eta, rng):
    d = TR[k]
    if d["kind"] == "cox":
        f = BASE[k]; E = rng.exponential(1 / np.exp(-eta)); C = rng.exponential(1 / ((1 - f) / f), len(eta))
        return (np.minimum(E, C), (E <= C).astype(int))
    p = 1 / (1 + np.exp(-(np.log(BASE[k] / (1 - BASE[k])) + eta)))
    return rng.binomial(1, p).astype(float)


def truth(k, sc):
    """eta for comparison k under scenario sc = (kind, param)."""
    a = TR[k]["arm"].astype(float); kind, par = sc
    if kind == "null":
        return np.zeros(len(a))
    if kind == "prog":
        return -0.26 * std(SCO[(k, "CAF_ECM")])
    T = PC1[k] if par[0] == "PC1" else std(SCO[(k, par[0])])
    return par[1] * a * T


# ------------------------------------------------------------------ check the fast Cox against lifelines
rng0 = np.random.default_rng([3036, 99, 0]); worst = worst0 = 0.0
for k in ("IMM150", "CM025", "JR101"):
    d = TR[k]; out = simulate(k, truth(k, ("pred", ("CD8A", 0.3))), rng0); s = std(SCO[(k, "CD8A")]); a = d["arm"]
    b1, se1 = cox_fast(np.column_stack([s, a, s * a]), out[0], out[1])
    df = pd.DataFrame({"s": s, "arm": a, "sx": s * a, "T": out[0], "E": out[1]})
    r0 = CoxPHFitter().fit(df, "T", "E").summary                                       # lifelines default precision
    r1 = CoxPHFitter().fit(df, "T", "E", fit_options={"precision": 1e-12}).summary      # tightened convergence
    worst0 = max(worst0, np.max(np.abs(b1 - r0["coef"].values)), np.max(np.abs(se1 - r0["se(coef)"].values)))
    worst = max(worst, np.max(np.abs(b1 - r1["coef"].values)), np.max(np.abs(se1 - r1["se(coef)"].values)))
print(f"fast Cox vs lifelines: max |difference| {worst0:.2e} at lifelines' default precision, {worst:.2e} at precision 1e-12")
assert worst < 1e-6


# ------------------------------------------------------------------ S0-S2: interaction and single-arm tests
def rep_basic(args):
    si, sc, r = args
    rng = np.random.default_rng([3036, si, r]); per = {}
    for k in ALL6:
        out = simulate(k, truth(k, sc), rng); d = TR[k]
        for m in MARK:
            per[(k, m, "int")] = sfit(fit_int, d, SCO[(k, m)], out)
            if sc[0] == "prog" and m == "CAF_ECM":
                per[(k, m, "within")] = sfit(fit_within, d, SCO[(k, m)], out)
    res = {}
    for m in MARK:
        q = pool_rule([per[(k, m, "int")][0] for k in ALL6], [per[(k, m, "int")][1] for k in ALL6]); res[(m, "int")] = (q["est"], q["p"])
    if sc[0] == "prog":
        q = pool_rule([per[(k, "CAF_ECM", "within")][0] for k in ALL6], [per[(k, "CAF_ECM", "within")][1] for k in ALL6])
        res[("CAF_ECM", "within")] = (q["est"], q["p"])
    return res

BASIC = [("S0", 0, ("null", None)), ("S1", 1, ("prog", None)), ("S2 beta 0.1", 2, ("pred", ("CD8A", 0.1))),
         ("S2 beta 0.2", 3, ("pred", ("CD8A", 0.2))), ("S2 beta 0.3", 4, ("pred", ("CD8A", 0.3)))]
rows = []
with mp.get_context("fork").Pool(NP) as pool:
    for lab, si, sc in BASIC:
        rs = pool.map(rep_basic, [(si, sc, r) for r in range(1000)], chunksize=4)
        for key in rs[0]:
            est = np.array([x[key][0] for x in rs]); p = np.array([x[key][1] for x in rs])
            rows.append(dict(scenario=lab, marker=key[0], test="interaction" if key[1] == "int" else "checkpoint arm alone",
                             replicates=len(rs), mean_est=np.nanmean(est), reject_p05=np.nanmean(p < 0.05),
                             reject_p05_benefit=np.nanmean((p < 0.05) & (est > 0)), failed=int(np.isnan(p).sum())))
        print(lab, "done", flush=True)
B = pd.DataFrame(rows); B.to_csv(f"{O}/E36_S0_S2.csv", index=False)
print(B[B.marker.isin(["CD8A", "CAF_ECM", "NetImmune"])].round(3).to_string(index=False))

# ------------------------------------------------------------------ S3: random-gene reference
NREF = 200
REFS = {}
for k in ALL6:
    for m in MARK:
        g, w = PRESW[(k, m)]; n = len(g); w = np.asarray(w, float)
        rng = np.random.default_rng([3036, 50, ALL6.index(k), n]); U = UNIV[k]; L = []
        for b in range(NREF):
            idx = rng.choice(U.shape[0], n, replace=False); ww = w[rng.permutation(n)]
            L.append(U[idx].T @ ww / np.abs(ww).sum())
        REFS[(k, m)] = np.array(L)


def rep_ref(args):
    si, sc, r = args
    rng = np.random.default_rng([3036, si, r]); obs, ref = {}, {}
    for k in ALL6:
        out = simulate(k, truth(k, sc), rng); d = TR[k]
        for m in MARK:
            obs[(k, m)] = sfit(fit_int, d, SCO[(k, m)], out)
            ref[(k, m)] = [sfit(fit_int, d, REFS[(k, m)][b], out) for b in range(NREF)]
    res = {}
    for m in MARK:
        o = pool_rule([obs[(k, m)][0] for k in ALL6], [obs[(k, m)][1] for k in ALL6])["est"]
        eb = np.array([pool_rule([ref[(k, m)][b][0] for k in ALL6], [ref[(k, m)][b][1] for k in ALL6])["est"] for b in range(NREF)])
        dr = DIRN[m]; nf = sum(int(not np.isfinite(ref[(k, m)][b][0])) for k in ALL6 for b in range(NREF))
        res[m] = (dr * o > np.nanpercentile(dr * eb, 95), o, np.nanmedian(eb), nf)
    return res

SREF = [("S3 null", 10, ("null", None)), ("S3 CD8A truth", 11, ("pred", ("CD8A", 0.3))), ("S3d diffuse (PC1) truth", 12, ("pred", ("PC1", 0.3)))]
rows = []
with mp.get_context("fork").Pool(NP) as pool:
    for lab, si, sc in SREF:
        rs = pool.map(rep_ref, [(si, sc, r) for r in range(200)], chunksize=1)
        for m in MARK:
            rows.append(dict(scenario=lab, marker=m, declared=DIRN[m], replicates=len(rs), exceed_rate=np.mean([x[m][0] for x in rs]),
                             mean_est=np.nanmean([x[m][1] for x in rs]), mean_ref_median=np.nanmean([x[m][2] for x in rs]),
                             failed_reference_fits=int(sum(x[m][3] for x in rs))))
        print(lab, "done", flush=True)
S3 = pd.DataFrame(rows); S3.to_csv(f"{O}/E36_S3.csv", index=False); print(S3.round(3).to_string(index=False))

# ------------------------------------------------------------------ S4: single-gene reference (paired bootstrap)
NB = 100
PM = ["T_CELL", "IFN", "APM", "TLS_B"]


def rep_pair(args):
    si, sc, r = args
    rng = np.random.default_rng([3036, si, r]); dd, vv = {m: [] for m in PM}, {m: [] for m in PM}
    for k in ALL6:
        out = simulate(k, truth(k, sc), rng); d = TR[k]; a = d["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
        idxs = [np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]) for _ in range(NB)]
        y0 = sfit(fit_int, d, SCO[(k, "CD8A")], out)[0]
        b0 = np.array([sfit(fit_int, take(d, ix), SCO[(k, "CD8A")][ix], sub_out(d, out, ix))[0] for ix in idxs])
        for m in PM:
            y1 = sfit(fit_int, d, SCO[(k, m)], out)[0]
            b1 = np.array([sfit(fit_int, take(d, ix), SCO[(k, m)][ix], sub_out(d, out, ix))[0] for ix in idxs])
            db = b1 - b0; ok = np.isfinite(db); dd[m].append(y1 - y0); vv[m].append(np.var(db[ok], ddof=1))
    res = {}
    for m in PM:
        dv, v = np.array(dd[m]), np.array(vv[m]); w = 1 / v; D_ = np.sum(w * dv) / w.sum(); se = np.sqrt(1 / w.sum())
        res[m] = (D_, se, 2 * norm.sf(abs(D_ / se)))
    return res

SPAIR = [("S4 truth CD8A", 20, ("pred", ("CD8A", 0.3)))] + [(f"S4 truth {m}", 21 + i, ("pred", (m, 0.3))) for i, m in enumerate(PM)]
rows = []
with mp.get_context("fork").Pool(NP) as pool:
    for lab, si, sc in SPAIR:
        rs = pool.map(rep_pair, [(si, sc, r) for r in range(200)], chunksize=1)
        for m in PM:
            D_ = np.array([x[m][0] for x in rs]); se = np.array([x[m][1] for x in rs]); p = np.array([x[m][2] for x in rs])
            rows.append(dict(scenario=lab, truth=sc[1][0], module=m, replicates=len(rs), oracle_diff=np.nanmean(D_),
                             median_se=np.nanmedian(se), sig_above=np.nanmean((p < 0.05) & (D_ > 0)), sig_below=np.nanmean((p < 0.05) & (D_ < 0))))
        print(lab, "done", flush=True)
S4 = pd.DataFrame(rows); S4.to_csv(f"{O}/E36_S4.csv", index=False); print(S4.round(3).to_string(index=False))

# ------------------------------------------------------------------ expectations
g = lambda sc, m, t="interaction": B[(B.scenario == sc) & (B.marker == m) & (B.test == t)].iloc[0]
e1 = {m: g("S0", m).reject_p05 for m in MARK}
e2a, e2b = g("S1", "CAF_ECM").reject_p05, g("S1", "CAF_ECM", "checkpoint arm alone").reject_p05
e3 = [g(f"S2 beta {b}", "CD8A").reject_p05 for b in (0.1, 0.2, 0.3)]
e4a = {m: S3[(S3.scenario == "S3 null") & (S3.marker == m)].exceed_rate.iloc[0] for m in MARK}
e4b = S3[(S3.scenario == "S3 CD8A truth") & (S3.marker == "T_CELL")].exceed_rate.iloc[0]
e5 = {m: S4[(S4.truth == "CD8A") & (S4.module == m)].sig_above.iloc[0] for m in PM}
own = S4[S4.truth == S4.module]
hi_, lo_ = own[own.oracle_diff >= 0.20], own[own.oracle_diff <= 0.10]
EXP = pd.DataFrame([
    dict(e="E1", expectation="S0 interaction rejection 2.5%-7.5%, each marker", result=str({m: round(v, 3) for m, v in e1.items()}),
         held=all(0.025 <= v <= 0.075 for v in e1.values())),
    dict(e="E2", expectation="S1: CAF_ECM interaction <= 7.5%; checkpoint arm alone >= 50%", result=f"{e2a:.3f}; {e2b:.3f}",
         held=bool(e2a <= 0.075 and e2b >= 0.5)),
    dict(e="E3", expectation="S2 CD8A power increasing; >= 60% at 0.3", result=str([round(x, 3) for x in e3]),
         held=bool(e3[0] < e3[1] < e3[2] and e3[2] >= 0.6)),
    dict(e="E4", expectation="S3 null: every marker exceeds <= 10%; CD8A truth: T_CELL >= 50%; diffuse: descriptive",
         result=f"null {({m: round(v, 3) for m, v in e4a.items()})}; T_CELL {e4b:.3f}", held=bool(max(e4a.values()) <= 0.10 and e4b >= 0.5)),
    dict(e="E5", expectation="S4 truth CD8A: module significantly above CD8A <= 5%, each", result=str({m: round(v, 3) for m, v in e5.items()}),
         held=max(e5.values()) <= 0.05),
    dict(e="E6", expectation="S4 truth module: power >= 50% where oracle >= 0.20; <= 30% where <= 0.10",
         result="; ".join(f"{r.module}: oracle {r.oracle_diff:+.3f}, power {r.sig_above:.3f}" for r in own.itertuples()),
         held=bool((hi_.sig_above >= 0.5).all() and (lo_.sig_above <= 0.3).all()) if len(hi_) + len(lo_) else None)])
EXP.to_csv(f"{O}/E36_expectations.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 200)
print("\n== expectations\n", EXP.to_string(index=False))
