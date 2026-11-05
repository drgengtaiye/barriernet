"""E43 - predicting absolute benefit with the immune and stroma axes in the six randomized comparisons (exploratory).

Protocol: analysis/E43_prereg.md (hash checked); frozen weights results/E37_frozen_weights.tsv (hash checked). Setup:
E34_trials.py up to its fitting stage, executed unchanged as in E37_test.py except that its step-0 and coverage tables
are written under E43 names. Six fixed models; 5-fold cross-validation stratified by arm and outcome, 20 repetitions
(seeds [2043, comparison, repetition]); C-for-benefit by rank matching and a benefit gradient at the median, with 300
within-arm bootstrap resamples of the evaluation (seeds [3043, comparison]).
Step 0: the full-data arm x CD8A coefficient of M_CD8A equals m_int (the manuscript's interaction) in every comparison.
"""
import hashlib, sys, warnings
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h43 = hashlib.sha256(open(f"{A}/E43_prereg.md", "rb").read()).hexdigest()
assert _h43 == open(f"{A}/E43_prereg.sha256").read().split()[0], "E43 protocol hash mismatch"
_hw = hashlib.sha256(open(f"{A}/results/E37_frozen_weights.tsv", "rb").read()).hexdigest()
assert _hw == open(f"{A}/E37_frozen.sha256").read().split()[0], "frozen weights hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E43_step0_reproduction").replace("E34_coverage_check.csv", "E43_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
print("E43 protocol sha256 verified:", _h43[:12], "| frozen weights sha256 verified:", _hw[:12])
warnings.filterwarnings("ignore")

FW = pd.read_csv(f"{O}/E37_frozen_weights.tsv", sep="\t"); WG = dict(zip(FW.gene, FW.weight))
T_STAR, NFOLD, NREP, NB = 12.0, 5, 20, 300
MODELS = {"M_CD8A": (["CD8A"], ["CD8A"]), "M_2AX": (["Immune", "CAF_ECM"], ["Immune"]), "M_BN03": (["BN03"], ["BN03"]),
          "M_CD8A_CAF": (["CD8A", "CAF_ECM"], ["CD8A"]), "M_2AXI": (["Immune", "CAF_ECM"], ["Immune", "CAF_ECM"]),
          "M_RISK": (["Immune", "CAF_ECM"], [])}
SS = {}
for k in ALL6:
    Z = TR[k]["Z"]; nc = set(nonconst(Z))
    pres = [g for g in WG if g in nc and np.isfinite(Z.loc[g].to_numpy(float)).all()]
    SS[k] = {"BN03": std(wscore(Z, pres, np.array([WG[g] for g in pres]))),
             **{s: std(OBS[(s, k)][0]) for s in ("CD8A", "Immune", "CAF_ECM")}}
    d = TR[k]; assert all(np.isfinite(v).all() for v in SS[k].values())
    if d["kind"] == "cox":
        assert np.nanmax(d["t"]) < 200, "survival times are expected in months"
    print(f"{k:13s} n {len(d['arm'])}  checkpoint {int(d['arm'].sum())}  "
          + (f"events {int(d['e'].sum())}" if d["kind"] == "cox" else f"pCR {int(d['y'].sum())}")
          + f"  barriernet 0.3 genes {len(pres)}")


def design(k, m, idx, arm):
    """Columns (name, vector) for patients idx with the arm set to `arm` (array or scalar)."""
    mains, inter = MODELS[m]; a = np.broadcast_to(np.asarray(arm, float), (len(idx),)).astype(float)
    cols = [(s, SS[k][s][idx]) for s in mains] + [("arm", a)] + [(f"x_{s}", SS[k][s][idx] * a) for s in inter]
    d = TR[k]
    if d["kind"] == "logit" and d["cov"] is not None:
        cv = np.asarray(d["cov"], float); cv = cv[:, None] if cv.ndim == 1 else cv
        cols += [(f"c{j}", cv[idx, j]) for j in range(cv.shape[1])]
    return cols


def fit_logit(X, y):
    b, V_, ok = ml(X, y)
    sep = (not ok) or np.any(np.abs(b) > 10) or np.any(~np.isfinite(np.diag(V_))) or np.any(np.sqrt(np.abs(np.diag(V_))) > 10)
    if sep:
        b, V_ = firth(X, y)
    if not np.all(np.isfinite(b)):
        raise RuntimeError("logistic fit failed")
    return b, V_


def fit_predict(k, m, tr, te):
    """Fit model m on patients tr; predicted absolute benefit for patients te; also the fitted coefficients."""
    d = TR[k]
    if d["kind"] == "logit":
        cols = design(k, m, tr, d["arm"][tr]); X = np.column_stack([np.ones(len(tr))] + [v for _, v in cols])
        b, V_ = fit_logit(X, d["y"][tr])
        X1 = np.column_stack([np.ones(len(te))] + [v for _, v in design(k, m, te, 1.0)])
        X0 = np.column_stack([np.ones(len(te))] + [v for _, v in design(k, m, te, 0.0)])
        ex = lambda Xm: 1 / (1 + np.exp(-np.clip(Xm @ b, -30, 30)))
        return ex(X1) - ex(X0), dict(zip(["const"] + [n for n, _ in cols], b)), V_
    cols = design(k, m, tr, d["arm"][tr]); df = pd.DataFrame(dict(cols)); df["T"] = d["t"][tr]; df["E"] = d["e"][tr]
    cph = CoxPHFitter().fit(df, "T", "E")
    s1 = cph.predict_survival_function(pd.DataFrame(dict(design(k, m, te, 1.0))), times=[T_STAR]).to_numpy()[0]
    s0 = cph.predict_survival_function(pd.DataFrame(dict(design(k, m, te, 0.0))), times=[T_STAR]).to_numpy()[0]
    return s1 - s0, cph.params_.to_dict(), cph.summary

# ------------------------------------------------------------------ step 0
for k in ALL6:
    n = len(TR[k]["arm"]); allp = np.arange(n); _, par, _ = fit_predict(k, "M_CD8A", allp, allp)
    yb = par["x_CD8A"] * (1 if TR[k]["kind"] == "logit" else -1); y0 = m_int(TR[k], OBS[("CD8A", k)][0])[0]
    print(f"step 0 {k}: M_CD8A interaction {yb:+.6f}  m_int {y0:+.6f}"); assert abs(yb - y0) < 1e-6

# ------------------------------------------------------------------ cross-validated predictions
KX = {k: KIDX[k] for k in ALL6}


def folds(k, r):
    d = TR[k]; out = d["y"] if d["kind"] == "logit" else d["e"]
    strata = d["arm"].astype(int) * 2 + np.asarray(out).astype(int)
    rng = np.random.default_rng([2043, KX[k], r]); f = np.empty(len(strata), int)
    for s in np.unique(strata):
        ix = np.where(strata == s)[0]; ix = ix[rng.permutation(len(ix))]; f[ix] = np.arange(len(ix)) % NFOLD
    return f


def cv_task(args):
    k, m, r = args; f = folds(k, r); pr = np.full(len(f), np.nan); nfail = 0
    for j in range(NFOLD):
        tr, te = np.where(f != j)[0], np.where(f == j)[0]
        try:
            pr[te] = fit_predict(k, m, tr, te)[0]
        except Exception:
            nfail += 1
    return pr, nfail

NP = int(sys.argv[1]) if len(sys.argv) > 1 else NPROC
ctasks = [(k, m, r) for k in ALL6 for m in MODELS for r in range(NREP)]
with mp.get_context("fork").Pool(NP) as pool:
    cres = pool.map(cv_task, ctasks, chunksize=1)
PRED, FAIL = {}, {}
for (k, m, r), (pr, nf) in zip(ctasks, cres):
    PRED.setdefault((k, m), []).append(pr); FAIL[(k, m)] = FAIL.get((k, m), 0) + nf
PRED = {key: np.nanmean(np.vstack(v), 0) for key, v in PRED.items()}
print("failed fold fits:", {f"{k}|{m}": v for (k, m), v in FAIL.items() if v} or "none")
assert all(np.isfinite(v).all() for v in PRED.values())

# ------------------------------------------------------------------ evaluation


def km_at(t, e, t0):
    if len(t) == 0:
        return np.nan
    s = 1.0
    for u in np.unique(t[(e == 1) & (t <= t0)]):
        at = np.sum(t >= u); dd = np.sum((t == u) & (e == 1)); s *= 1 - dd / at
    return s


def evaluate(k, pred, idx):
    """C-for-benefit and benefit gradient on the patients idx (possibly with repeats)."""
    d = TR[k]; a = d["arm"][idx]; p = pred[idx]
    i1, i0 = idx[a == 1], idx[a == 0]
    i1 = i1[np.argsort(pred[i1], kind="stable")]; i0 = i0[np.argsort(pred[i0], kind="stable")]
    small, large = (i1, i0) if len(i1) <= len(i0) else (i0, i1)
    ns, nl = len(small), len(large); partner = large[np.floor((np.arange(ns) + 0.5) * nl / ns).astype(int)]
    t1, t0 = (small, partner) if small is i1 else (partner, small)
    pp = (pred[t1] + pred[t0]) / 2
    if d["kind"] == "logit":
        ob = d["y"][t1] - d["y"][t0]
    else:
        T1, T0, E1, E0 = d["t"][t1], d["t"][t0], d["e"][t1], d["e"][t0]
        ob = np.where((T1 > T0) & (E0 == 1), 1.0, np.where((T0 > T1) & (E1 == 1), -1.0,
                      np.where((T1 == T0) & (E1 == 1) & (E0 == 1), 0.0, np.nan)))
    ok = np.isfinite(ob); pp, ob = pp[ok], ob[ok]
    gt = ob[:, None] > ob[None, :]
    num = np.sum(gt & (pp[:, None] > pp[None, :])) + 0.5 * np.sum(gt & (pp[:, None] == pp[None, :]))
    C = num / gt.sum() if gt.sum() else np.nan
    up = p > np.median(p); g = []
    for h in (up, ~up):
        ih = idx[h]; ah = d["arm"][ih]
        if d["kind"] == "logit":
            g.append(d["y"][ih][ah == 1].mean() - d["y"][ih][ah == 0].mean())
        else:
            g.append(km_at(d["t"][ih][ah == 1], d["e"][ih][ah == 1], T_STAR) - km_at(d["t"][ih][ah == 0], d["e"][ih][ah == 0], T_STAR))
    return C, g[0] - g[1], int(ok.sum())

BOOTI = {}
for k in ALL6:
    rng = np.random.default_rng([3043, KX[k]]); a = TR[k]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
    BOOTI[k] = [np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]) for _ in range(NB)]


def ev_task(args):
    k, m = args; pred = PRED[(k, m)]; n = len(pred)
    obs = evaluate(k, pred, np.arange(n)); bs = np.array([evaluate(k, pred, ix)[:2] for ix in BOOTI[k]])
    return obs, bs

etasks = [(k, m) for k in ALL6 for m in MODELS]
with mp.get_context("fork").Pool(min(NP, len(etasks))) as pool:
    eres = pool.map(ev_task, etasks, chunksize=1)
EV = dict(zip(etasks, eres))
PER = pd.DataFrame([dict(comparison=k, model=m, kind=TR[k]["kind"], C=o[0], C_se=np.nanstd(b[:, 0], ddof=1), gradient=o[1],
                         g_se=np.nanstd(b[:, 1], ddof=1), pairs_used=o[2], pred_sd=np.std(PRED[(k, m)]),
                         pred_mean=np.mean(PRED[(k, m)]))
                    for (k, m), (o, b) in EV.items()])
PER.to_csv(f"{O}/E43_per_comparison.csv", index=False)


def contrast(m1, m0, j, keys):
    dd, vv = [], []
    for k in keys:
        (o1, b1), (o0, b0) = EV[(k, m1)], EV[(k, m0)]; db = b1[:, j] - b0[:, j]; ok = np.isfinite(db)
        dd.append(o1[j] - o0[j]); vv.append(np.var(db[ok], ddof=1))
    dd, vv = np.array(dd), np.array(vv); w = 1 / vv; D_ = np.sum(w * dd) / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(diff=D_, se=se, lo=D_ - 1.96 * se, hi=D_ + 1.96 * se, p=2 * norm.sf(abs(D_ / se)), mdd=2.80 * se, k=len(keys))

SETS = [("six (primary)", ALL6), ("survival three", RENAL), ("pCR three", BREAST)]
POOL, CON = [], []
for lab, keys in SETS:
    for m in MODELS:
        s = PER[(PER.model == m) & PER.comparison.isin(keys)]
        c = pool_rule(s.C.values - 0.5, s.C_se.values); g = pool_rule(s.gradient.values, s.g_se.values)
        POOL.append(dict(set=lab, model=m, C=c["est"] + 0.5, C_lo=c["lo"] + 0.5, C_hi=c["hi"] + 0.5, C_p=c["p"],
                         gradient=g["est"], g_lo=g["lo"], g_hi=g["hi"], g_p=g["p"], k=c["k"]))
        if m != "M_CD8A":
            for j, metric in ((0, "C-for-benefit"), (1, "gradient")):
                CON.append(dict(set=lab, contrast=f"{m} - M_CD8A", metric=metric, **contrast(m, "M_CD8A", j, keys)))
    for j, metric in ((0, "C-for-benefit"), (1, "gradient")):
        CON.append(dict(set=lab, contrast="M_2AX - M_BN03", metric=metric, **contrast("M_2AX", "M_BN03", j, keys)))
POOL, CON = pd.DataFrame(POOL), pd.DataFrame(CON)
POOL.to_csv(f"{O}/E43_pooled.csv", index=False); CON.to_csv(f"{O}/E43_contrasts.csv", index=False)

# full-data coefficients of M_2AX and M_CD8A_CAF (descriptive: is CAF_ECM prognostic in these trials?)
COEF = []
for k in ALL6:
    n = len(TR[k]["arm"]); allp = np.arange(n)
    for m in ("M_2AX", "M_CD8A_CAF"):
        _, par, _ = fit_predict(k, m, allp, allp)
        COEF.append(dict(comparison=k, model=m, **{x: v for x, v in par.items() if x in ("CAF_ECM", "Immune", "CD8A", "x_Immune", "x_CD8A")}))
pd.DataFrame(COEF).to_csv(f"{O}/E43_coefficients.csv", index=False)

six = POOL[POOL.set == "six (primary)"].set_index("model"); c6 = CON[(CON.set == "six (primary)") & (CON.metric == "C-for-benefit")].set_index("contrast")
EXP = pd.DataFrame([
    dict(q="B1", expectation="M_CD8A pooled C-for-benefit > 0.5, CI excluding 0.5",
         result=f"{six.loc['M_CD8A', 'C']:.3f} ({six.loc['M_CD8A', 'C_lo']:.3f} to {six.loc['M_CD8A', 'C_hi']:.3f})",
         held=bool(six.loc["M_CD8A", "C_lo"] > 0.5)),
    dict(q="B2 (primary)", expectation="M_2AX - M_CD8A: CI includes 0",
         result=f"{c6.loc['M_2AX - M_CD8A', 'diff']:+.3f} ({c6.loc['M_2AX - M_CD8A', 'lo']:+.3f} to "
                f"{c6.loc['M_2AX - M_CD8A', 'hi']:+.3f}); MDD {c6.loc['M_2AX - M_CD8A', 'mdd']:.3f}",
         held=bool(c6.loc["M_2AX - M_CD8A", "lo"] <= 0 <= c6.loc["M_2AX - M_CD8A", "hi"])),
    dict(q="B3", expectation="|M_CD8A_CAF - M_CD8A| < 0.02",
         result=f"{c6.loc['M_CD8A_CAF - M_CD8A', 'diff']:+.3f}", held=bool(abs(c6.loc["M_CD8A_CAF - M_CD8A", "diff"]) < 0.02)),
    dict(q="B4", expectation="M_RISK pooled C-for-benefit CI includes 0.5",
         result=f"{six.loc['M_RISK', 'C']:.3f} ({six.loc['M_RISK', 'C_lo']:.3f} to {six.loc['M_RISK', 'C_hi']:.3f})",
         held=bool(six.loc["M_RISK", "C_lo"] <= 0.5 <= six.loc["M_RISK", "C_hi"]))])
EXP.to_csv(f"{O}/E43_expectations.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print("\n== expectations\n", EXP.to_string(index=False))
print("\n== pooled\n", POOL.round(3).to_string(index=False))
print("\n== contrasts\n", CON.round(3).to_string(index=False))
print("\n== per comparison\n", PER.round(3).to_string(index=False))
print("\n== full-data coefficients\n", pd.DataFrame(COEF).round(3).to_string(index=False))
