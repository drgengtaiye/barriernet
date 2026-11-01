"""E28 - comparator audit, interaction calibration and single-sample scoring in the six primary randomized
comparisons. Pre-registered in E28_prereg.md (2026-10-20; sha256 in E28_prereg.sha256) before any model below
was fitted. Data, scores, models and the pooling rule come from E28_base.py (E25 code; step 0 passes on import).

  Q1b  official ESTIMATE (E28_estimate.R output) against the z-mean contrast
  Q2   Immune against CAF_ECM: difference of benefit-oriented interactions, adjusted and marginal
  Q3   parametric bootstrap under no score x arm interaction: size of the nominal test, calibrated P
  Q5   single-sample scorings: frozen TCGA reference, unit weights, within-sample rank
  Q6   NetImmune against ESTIMATE: paired difference and joint model
  Q7   PD-L1 status in JAVELIN Renal 101 (exploratory)
Q4 needs no model and is computed in E28_tiers.py.  Usage: python3 E28_trials.py <E28_estimate.R workdir>
"""
import sys
from E28_base import *
from E28_base import _cmp

W = sys.argv[1]
NBOOT, NNULL, NPROC = 2000, 1000, 64
SETS = {"all six": ALL6, "renal": RENAL, "breast": BREAST}
TISSUE = {"IMM150": "KIRC", "CM025": "KIRC", "JR101": "KIRC", "ISPY2_pembro": "BRCA", "ISPY2_durva": "BRCA",
          "NeoTRIP": "BRCA"}
P5E = "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan/expr"


def pool_table(rows, keys, label):
    """rows: dict key -> (y, se). Manuscript rule plus Stouffer, for each set."""
    out = []
    for sn, ks in SETS.items():
        y = [rows[k][0] for k in ks]; se = [rows[k][1] for k in ks]
        r = pool_rule(y, se); z = np.nansum(np.array(y) / np.array(se)) / np.sqrt(len(ks))
        out.append(dict(label, set=sn, k=r["k"], est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"],
                        n_positive=int(np.sum(np.array(y) > 0)), stouffer_z=z, stouffer_p=2 * norm.sf(abs(z))))
    return out


def lfitV(X, y, firth_ok=True):
    """Full coefficient vector and covariance, with lfit's separation rule."""
    b, V_, ok = ml(X, y)
    sep = (not ok) or np.any(np.abs(b) > 10) or np.any(~np.isfinite(np.diag(V_))) or np.any(np.sqrt(np.diag(V_)) > 10)
    if sep and firth_ok:
        b, V_ = firth(X, y)
    elif sep:
        return None, None
    return b, V_


def joint_cov(d, s1, s2):
    """Both products in one model: benefit-oriented (y1, y2), variances and covariance."""
    s1, s2 = std(s1), std(s2); a = d["arm"]
    if d["kind"] == "cox":
        df = pd.DataFrame({"s1": s1, "s2": s2, "arm": a, "x1": s1 * a, "x2": s2 * a, "T": d["t"], "E": d["e"]})
        f = CoxPHFitter().fit(df, "T", "E"); V_ = f.variance_matrix_
        return -f.params_["x1"], -f.params_["x2"], V_.loc["x1", "x1"], V_.loc["x2", "x2"], V_.loc["x1", "x2"]
    X = np.column_stack([np.ones_like(s1), s1, s2, a, s1 * a, s2 * a] + ([d["cov"]] if d["cov"] is not None else []))
    b, V_ = lfitV(X, d["y"], d["firth"])
    return b[4], b[5], V_[4, 4], V_[5, 5], V_[4, 5]


# ------------------------------------------------------------------ Q1b: official ESTIMATE
print("\n" + "=" * 90 + "\nQ1b official ESTIMATE in the six comparisons")
q1, q1rows = [], {}
for k in ALL6:
    d = TR[k]; o = pd.read_csv(f"{W}/{k}_official_scores.csv"); o["sample"] = o["sample"].astype(str)
    o = o.set_index("sample").loc[d["S"].index.astype(str)]
    for v in ["as_is", "na_min"]:
        s = std(o[f"imm_{v}"]) - std(o[f"str_{v}"])
        y, se = m_int(d, s); q1rows[(v, k)] = (y, se)
        q1.append(dict(trial=LABEL[k], variant=f"official, {v}", r_with_zmean=np.corrcoef(s, d["S"]["ESTIMATE"])[0, 1],
                       y=y, se=se, p=2 * norm.sf(abs(y / se))))
    y, se = m_int(d, d["S"]["ESTIMATE"].values); q1rows[("zmean", k)] = (y, se)
    q1.append(dict(trial=LABEL[k], variant="z-mean (manuscript)", r_with_zmean=1.0, y=y, se=se, p=2 * norm.sf(abs(y / se))))
Q1 = pd.DataFrame(q1); Q1.to_csv(f"{O}/E28_Q1_trials_per_comparison.csv", index=False)
Q1P = pd.DataFrame(sum([pool_table({k: q1rows[(v, k)] for k in ALL6}, ALL6, {"variant": v})
                        for v in ["zmean", "as_is", "na_min"]], []))
Q1P.to_csv(f"{O}/E28_Q1_trials_pooled.csv", index=False)
print(Q1.round(3).to_string(index=False)); print(Q1P.round(4).to_string(index=False))


# ------------------------------------------------------------------ Q2a, Q6b: joint models with covariance
print("\n" + "=" * 90 + "\nQ2a / Q6b joint models")
JROWS, jrec = {}, []
for k in ALL6:
    d = TR[k]
    for tag, a_, b_ in [("Q2a", "Immune", "CAF_ECM"), ("Q6b", "NetImmune", "ESTIMATE")]:
        y1, y2, v1, v2, c12 = joint_cov(d, d["S"][a_].values, d["S"][b_].values)
        if tag == "Q2a":   # same estimates as E25 Q2
            e = m_joint(d, d["S"]["Immune"].values, d["S"]["CAF_ECM"].values)
            assert abs(e[0] - y1) < 1e-8 and abs(e[2] - y2) < 1e-8 and abs(e[1] - np.sqrt(v1)) < 1e-8
            dlt, vd = y1 + y2, v1 + v2 + 2 * c12          # theta(Immune) - theta(-CAF_ECM)
        else:
            dlt, vd = y1 - y2, v1 + v2 - 2 * c12
        JROWS[(tag, k)] = (dlt, np.sqrt(vd)); JROWS[(tag + "_1", k)] = (y1, np.sqrt(v1)); JROWS[(tag + "_2", k)] = (y2, np.sqrt(v2))
        jrec.append(dict(test=tag, trial=LABEL[k], score1=a_, score2=b_, y1=y1, se1=np.sqrt(v1), y2=y2, se2=np.sqrt(v2),
                         corr_estimates=c12 / np.sqrt(v1 * v2), delta=dlt, delta_se=np.sqrt(vd),
                         delta_p=2 * norm.sf(abs(dlt / np.sqrt(vd)))))
JR = pd.DataFrame(jrec); JR.to_csv(f"{O}/E28_Q2a_Q6b_joint.csv", index=False); print(JR.round(3).to_string(index=False))


# ------------------------------------------------------------------ Q2b, Q6a: bootstrap covariance of marginal interactions
BS = ["Immune", "CAF_ECM", "NetImmune", "ESTIMATE"]


def one_boot(job):
    k, b = job; d = TR[k]; r = np.random.default_rng([2028, ALL6.index(k), b])
    ix = np.concatenate([r.choice(np.where(d["arm"] == a)[0], int((d["arm"] == a).sum()), replace=True) for a in (0, 1)])
    dd = sub(d, ix); out = dict(key=k, b=b)
    for sc in BS:
        try:
            out[sc] = m_int(dd, d["S"][sc].values[ix])[0]
        except Exception:
            out[sc] = np.nan
    return out


print(f"\nQ2b / Q6a: {NBOOT} arm-stratified bootstrap resamples x 6 comparisons ...")
with mp.get_context("fork").Pool(NPROC) as pp:
    BT = pd.DataFrame(pp.map(one_boot, [(k, b) for k in ALL6 for b in range(NBOOT)], chunksize=20))
BT.to_csv(f"{O}/E28_bootstrap_draws.csv", index=False)
mrec = []
for k in ALL6:
    d = TR[k]; bt = BT[BT.key == k]
    th = {sc: m_int(d, d["S"][sc].values) for sc in BS}
    for tag, a_, b_, sg in [("Q2b", "Immune", "CAF_ECM", +1), ("Q6a", "NetImmune", "ESTIMATE", -1)]:
        ok = bt[[a_, b_]].dropna(); rho = np.corrcoef(ok[a_], ok[b_])[0, 1]
        (y1, s1), (y2, s2) = th[a_], th[b_]
        dlt = y1 + sg * y2; se = np.sqrt(s1 ** 2 + s2 ** 2 + 2 * sg * rho * s1 * s2)
        se_boot = np.std(ok[a_] + sg * ok[b_], ddof=1)
        JROWS[(tag, k)] = (dlt, se)
        mrec.append(dict(test=tag, trial=LABEL[k], score1=a_, score2=b_, y1=y1, se1=s1, y2=y2, se2=s2, boot_corr=rho,
                         n_boot=len(ok), delta=dlt, delta_se=se, delta_se_pure_bootstrap=se_boot,
                         delta_p=2 * norm.sf(abs(dlt / se))))
MR = pd.DataFrame(mrec); MR.to_csv(f"{O}/E28_Q2b_Q6a_marginal.csv", index=False); print(MR.round(3).to_string(index=False))
DP = pd.DataFrame(sum([pool_table({k: JROWS[(t, k)] for k in ALL6}, ALL6, {"test": t})
                       for t in ["Q2a", "Q2b", "Q6a", "Q6b", "Q6b_1", "Q6b_2"]], []))
DP.to_csv(f"{O}/E28_Q2_Q6_pooled.csv", index=False); print(DP.round(4).to_string(index=False))


# ------------------------------------------------------------------ Q3: parametric bootstrap under no interaction
NS3 = ["NetImmune", "Immune", "CAF_ECM", "ESTIMATE", "CD8A"]


def breslow(t, e, lp):
    o = np.argsort(t); t, e, w = t[o], e[o], np.exp(lp[o])
    risk = np.cumsum(w[::-1])[::-1]
    ut = np.unique(t[e == 1]); H = np.empty(len(ut)); acc = 0.0
    for j, u in enumerate(ut):
        i0 = np.searchsorted(t, u, "left"); acc += np.sum((t == u) & (e == 1)) / risk[i0]; H[j] = acc
    return np.r_[0.0, ut], np.r_[0.0, H]


def reverse_km(t, e):
    """Censoring survival G at its jump times (censorings as events; ties: censoring after event)."""
    c = 1 - e; ut = np.unique(t[c == 1]); G, g = [], 1.0
    for u in ut:
        n_at = np.sum(t >= u) - np.sum((t == u) & (e == 1)); g *= 1 - np.sum((t == u) & (c == 1)) / n_at; G.append(g)
    return ut, np.array(G)


def null_model(d, s):
    """Fitted no-interaction model, ready for simulation."""
    s = std(s); a = d["arm"]
    if d["kind"] == "cox":
        f = CoxPHFitter().fit(pd.DataFrame({"s": s, "arm": a, "T": d["t"], "E": d["e"]}), "T", "E")
        lp = f.params_["s"] * s + f.params_["arm"] * a
        tg, H = breslow(d["t"], d["e"], lp); ct, G = reverse_km(d["t"], d["e"])
        return dict(lp=lp, tg=tg, H=H, ct=ct, G=G, tmax=d["t"].max())
    X0 = np.column_stack([np.ones_like(s), s, a] + ([d["cov"]] if d["cov"] is not None else []))
    b, _ = lfitV(X0, d["y"], True)
    return dict(p=1 / (1 + np.exp(-X0 @ b)))


def simulate(d, nm, r):
    dd = dict(d)
    if d["kind"] == "cox":
        n = len(nm["lp"]); target = -np.log(r.random(n)) / np.exp(nm["lp"])
        T = np.where(target <= nm["H"][-1], np.interp(target, nm["H"], nm["tg"]), np.inf)
        u = r.random(n); j = np.searchsorted(-nm["G"], -u, "left")          # first jump with G <= u
        C = np.where(j < len(nm["ct"]), nm["ct"][np.minimum(j, len(nm["ct"]) - 1)], np.inf)
        obs = np.minimum(np.minimum(T, C), nm["tmax"]); ev = (T <= np.minimum(C, nm["tmax"])).astype(int)
        dd["t"], dd["e"] = np.maximum(obs, 1e-6), ev
    else:
        dd["y"] = (r.random(len(nm["p"])) < nm["p"]).astype(float)
    return dd


def one_null(job):
    k, sc, b0, nb = job; d = TR[k]; s = d["S"][sc].values; m = np.isfinite(s); d = sub(d, m) if not m.all() else d
    s = s[m]; nm = null_model(d, s); out = []
    for b in range(b0, b0 + nb):
        r = np.random.default_rng([2029, ALL6.index(k), NS3.index(sc), b])
        try:
            y, se = m_int(simulate(d, nm, r), s)
        except Exception:
            y, se = np.nan, np.nan
        out.append(dict(key=k, score=sc, b=b, y=y, se=se))
    return out


print(f"\nQ3: parametric bootstrap, {NNULL} null data sets x {len(NS3)} scores x 6 comparisons ...")
jobs = [(k, sc, b0, 25) for k in ALL6 for sc in NS3 for b0 in range(0, NNULL, 25)]
with mp.get_context("fork").Pool(NPROC) as pp:
    NB = pd.DataFrame([r for blk in pp.map(one_null, jobs, chunksize=1) for r in blk])
NB.to_csv(f"{O}/E28_Q3_null_draws.csv", index=False)
q3, q3p = [], []
for sc in NS3:
    obs = {}
    for k in ALL6:
        d = TR[k]; s = d["S"][sc].values; m = np.isfinite(s)
        y, se = m_int(sub(d, m) if not m.all() else d, s[m]); obs[k] = (y, se)
        nb = NB[(NB.key == k) & (NB.score == sc)].dropna(); z = nb.y / nb.se
        q3.append(dict(score=sc, trial=LABEL[k], obs_y=y, obs_se=se, obs_p_nominal=2 * norm.sf(abs(y / se)), n_null=len(nb),
                       null_mean_y=nb.y.mean(), null_sd_z=z.std(), pct_null_nominal_p05=100 * np.mean(np.abs(z) > norm.ppf(0.975)),
                       pct_null_nominal_p01=100 * np.mean(np.abs(z) > norm.ppf(0.995)),
                       p_calibrated=(1 + np.sum(np.abs(z) >= abs(y / se))) / (1 + len(nb))))
    po = pool_rule([obs[k][0] for k in ALL6], [obs[k][1] for k in ALL6])
    zo = sum(obs[k][0] / obs[k][1] for k in ALL6) / np.sqrt(6)
    W6 = NB[NB.score == sc].pivot(index="b", columns="key", values=["y", "se"]).dropna()
    pb = np.array([pool_rule(W6.loc[b, "y"][ALL6].values, W6.loc[b, "se"][ALL6].values)["p"] for b in W6.index])
    zb = (W6["y"][ALL6].values / W6["se"][ALL6].values).sum(1) / np.sqrt(6)
    q3p.append(dict(score=sc, obs_est=po["est"], obs_p_nominal=po["p"], n_null=len(pb),
                    pct_null_pooled_p05=100 * np.mean(pb < 0.05), p_calibrated=(1 + np.sum(pb <= po["p"])) / (1 + len(pb)),
                    obs_stouffer_p=2 * norm.sf(abs(zo)), pct_null_stouffer_p05=100 * np.mean(np.abs(zb) > norm.ppf(0.975)),
                    stouffer_p_calibrated=(1 + np.sum(np.abs(zb) >= abs(zo))) / (1 + len(zb))))
Q3 = pd.DataFrame(q3); Q3.to_csv(f"{O}/E28_Q3_calibration_per_comparison.csv", index=False)
Q3P = pd.DataFrame(q3p); Q3P.to_csv(f"{O}/E28_Q3_calibration_pooled.csv", index=False)
print(Q3.round(4).to_string(index=False)); print(Q3P.round(4).to_string(index=False))


# ------------------------------------------------------------------ Q5: single-sample scorings
print("\n" + "=" * 90 + "\nQ5 single-sample scorings")
GENES48 = sorted({g for m in ["CAF_ECM", "TAM_SUP", "APM", "IFN", "TLS_B", "T_CELL"] for g in MOD[m]})
REF = {}
for tis in ["KIRC", "BRCA"]:
    T_ = pd.read_csv(f"{P5E}/{tis}.gz", sep="\t", index_col=0)
    T_ = T_.loc[:, [c for c in T_.columns if int(c[13:15]) < 10]]                       # primary tumors
    T_ = T_.loc[T_.index.intersection(GENES48)]
    REF[tis] = (T_.mean(1), T_.std(1), T_.shape[1])
    print(f"  TCGA-{tis}: {T_.shape[1]} tumors, {T_.shape[0]} of 48 genes")


def axes(mod):
    S_ = pd.DataFrame(mod); S_["Barrier"] = S_[["CAF_ECM", "TAM_SUP"]].mean(1)
    S_["Immune"] = S_[["APM", "IFN", "TLS_B", "T_CELL"]].mean(1); S_["NetImmune"] = S_.Immune - S_.Barrier
    return S_


def modules(V, genes_ok):
    out = {}
    for m in ["CAF_ECM", "TAM_SUP", "APM", "IFN", "TLS_B", "T_CELL"]:
        g = [x for x in MOD[m] if x in genes_ok]
        out[m] = V.loc[g].mean(0) if len(g) >= 3 else pd.Series(np.nan, index=V.columns)
    return out


def single_sample(k):
    X = TR[k]["X"].astype(float); X = X[~X.index.duplicated()]
    mu, sd, _ = REF[TISSUE[k]]; ga = [g for g in X.index if g in mu.index and sd[g] > 0]
    Za = ((X.loc[ga].T - mu[ga]) / sd[ga]).T
    R_ = X.rank(axis=0, na_option="keep") / X.notna().sum(0)
    return {"frozen TCGA reference": axes(modules(Za, set(ga))), "unit weights": axes(modules(X, set(X.index))),
            "within-sample rank": axes(modules(R_, set(X.index)))}


q5, q5rows = [], {}
for k in ALL6:
    d = TR[k]; ss = single_sample(k)
    for v, S_ in ss.items():
        for sc in ["NetImmune", "Immune", "CAF_ECM"]:
            s = S_[sc].values; ref = d["S"][sc].values; m = np.isfinite(s) & np.isfinite(ref)
            y, se = m_int(sub(d, m) if not m.all() else d, s[m]); q5rows[(v, sc, k)] = (y, se)
            q5.append(dict(scoring=v, score=sc, trial=LABEL[k], r_with_cohort_relative=np.corrcoef(s[m], ref[m])[0, 1],
                           y=y, se=se, y_cohort_relative=m_int(d, ref)[0]))
Q5 = pd.DataFrame(q5); Q5.to_csv(f"{O}/E28_Q5_single_sample_per_comparison.csv", index=False)
Q5P = pd.DataFrame(sum([pool_table({k: q5rows[(v, sc, k)] for k in ALL6}, ALL6, {"scoring": v, "score": sc})
                        for v in ["frozen TCGA reference", "unit weights", "within-sample rank"]
                        for sc in ["NetImmune", "Immune", "CAF_ECM"]], []))
Q5P.to_csv(f"{O}/E28_Q5_single_sample_pooled.csv", index=False)
print(Q5.groupby(["scoring", "score"]).r_with_cohort_relative.agg(["min", "median"]).round(3).to_string())
print(Q5P[Q5P.set == "all six"].round(4).to_string(index=False))


# ------------------------------------------------------------------ Q7: PD-L1 in JAVELIN Renal 101 (exploratory)
print("\n" + "=" * 90 + "\nQ7 JAVELIN Renal 101, PD-L1")
d = TR["JR101"]; pd1 = d["C"]["PDL1FL"].map({"Y": 1.0, "N": 0.0}).values; m = np.isfinite(pd1)
q7 = []
for sc in ["NetImmune", "Immune", "CAF_ECM", "ESTIMATE"]:
    s = std(d["S"][sc].values[m]); a = d["arm"][m]; p_ = pd1[m]
    for adj in [False, True]:
        cols = {"s": s, "arm": a, "sx": s * a}
        if adj:
            cols.update({"pdl1": p_, "px": p_ * a})
        df = pd.DataFrame(cols); df["T"] = d["t"][m]; df["E"] = d["e"][m]
        r = CoxPHFitter().fit(df, "T", "E").summary
        q7.append(dict(score=sc, model="+ PD-L1 and PD-L1 x arm" if adj else "same patients, unadjusted", n=int(m.sum()),
                       y=-r.loc["sx", "coef"], se=r.loc["sx", "se(coef)"], p=r.loc["sx", "p"],
                       pdl1_int_y=-r.loc["px", "coef"] if adj else np.nan, pdl1_int_p=r.loc["px", "p"] if adj else np.nan))
Q7 = pd.DataFrame(q7); Q7.to_csv(f"{O}/E28_Q7_jr101_pdl1.csv", index=False); print(Q7.round(4).to_string(index=False))
print("\nE28 trials done")
