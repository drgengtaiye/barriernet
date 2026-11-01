"""E40 step 2 - checkpoint arm alone, control arm and single-arm cohorts against the randomized interaction (52 signatures).

Protocol: analysis/E40_prereg.md (hash checked). Randomized: E34_trials.py up to its fitting stage (E32 setup), executed
unchanged except that its step-0 and coverage tables are written under E40 names, as in E37_test.py. Single-arm:
results/E40_single_arm_genes.csv.gz (E40_export.R; 48-gene step 0 against E37 passed there), minus the three samples
removed by the E37 twin check. Raw-score orientation throughout (benefit-oriented: higher score, more benefit).
Step 0: pooled interactions reproduce E32_primary.csv (< 1e-8).
"""
import hashlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h40 = hashlib.sha256(open(f"{A}/E40_prereg.md", "rb").read()).hexdigest()
assert _h40 == open(f"{A}/E40_prereg.sha256").read().split()[0], "E40 protocol hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E40_step0_reproduction").replace("E34_coverage_check.csv", "E40_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
from scipy.stats import spearmanr, t as tdist_
print("E40 protocol sha256 verified:", _h40[:12])

# ------------------------------------------------------------------ randomized: interaction and within-arm, per comparison
tasks = [(s, k) for s in LIB for k in ALL6 if (s, k) in OBS]


def fit3(task):
    s, k = task; sc = OBS[(s, k)][0]; d = TR[k]; out = []
    for f in (lambda: m_int(d, sc), lambda: m_within(d, sc, 1), lambda: m_within(d, sc, 0)):
        try:
            out.append(f())
        except Exception:
            out.append((np.nan, np.nan))
    return out

with mp.get_context("fork").Pool(NPROC) as pool:
    res = pool.map(fit3, tasks, chunksize=1)
PER = pd.DataFrame([dict(signature=s, comparison=k, int_y=r[0][0], int_se=r[0][1], icb_y=r[1][0], icb_se=r[1][1],
                         ctl_y=r[2][0], ctl_se=r[2][1]) for (s, k), r in zip(tasks, res)])
PER.to_csv(f"{O}/E40_randomized_per_comparison.csv", index=False)
RP = []
for s, g in PER.groupby("signature", sort=False):
    row = dict(signature=s, k=len(g))
    for c in ("int", "icb", "ctl"):
        p_ = pool_rule(g[f"{c}_y"], g[f"{c}_se"]) if len(g) >= 4 else dict(est=np.nan, lo=np.nan, hi=np.nan, p=np.nan, k=len(g))
        row.update({f"{c}_est": p_["est"], f"{c}_lo": p_["lo"], f"{c}_hi": p_["hi"], f"{c}_p": p_["p"]})
    RP.append(row)
RP = pd.DataFrame(RP).set_index("signature")
E32P = pd.read_csv(f"{O}/E32_primary.csv").set_index("signature")
d0 = float(np.nanmax(np.abs(RP["int_est"] - E32P.loc[RP.index, "est"]))); p0 = float(np.nanmax(np.abs(RP["int_p"] - E32P.loc[RP.index, "p"])))
print(f"step 0: pooled interactions vs E32_primary, max |est diff| {d0:.2e}, max |p diff| {p0:.2e}"); assert d0 < 1e-8 and p0 < 1e-8

# ------------------------------------------------------------------ single-arm cohorts
SA = pd.read_csv(f"{O}/E40_single_arm_genes.csv.gz")
fl = pd.read_csv(f"{O}/E37_twin_flags.csv")
drop = SA.set_index(["cohort", "sample"]).index.isin(list(zip(fl.development, fl["sample"])))
print("removed as possible twins:", int(drop.sum())); assert drop.sum() == 3
SA = SA[~drop]; GC = [c for c in SA.columns if c not in ("cohort", "set", "sample", "response")]
print(f"single-arm: {SA.cohort.nunique()} cohorts, {len(SA)} patients, {int(SA.response.sum())} responders")


def dl_hk(y, se):
    """DerSimonian-Laird; CI and P from the wider of Hartung-Knapp (DL tau2) and DL, as E12 (metafor DL / knha)."""
    y, se = np.asarray(y, float), np.asarray(se, float); ok = np.isfinite(y) & np.isfinite(se) & (se > 0); y, v = y[ok], se[ok] ** 2
    k = len(y); w0 = 1 / v; m0 = np.sum(w0 * y) / w0.sum(); Q = np.sum(w0 * (y - m0) ** 2)
    t2 = max(0.0, (Q - (k - 1)) / (w0.sum() - np.sum(w0 ** 2) / w0.sum())); w = 1 / (v + t2); mu = np.sum(w * y) / w.sum()
    sd = np.sqrt(1 / w.sum()); dl = (mu - 1.959964 * sd, mu + 1.959964 * sd, 2 * norm.sf(abs(mu / sd)))
    sh = np.sqrt(np.sum(w * (y - mu) ** 2) / ((k - 1) * w.sum())); q = tdist_.ppf(0.975, k - 1)
    hk = (mu - q * sh, mu + q * sh, 2 * tdist_.sf(abs(mu / sh), k - 1)) if sh > 0 else dl
    lo, hi = (hk[0], hk[1]) if (hk[1] - hk[0]) >= (dl[1] - dl[0]) else (dl[0], dl[1])
    return dict(est=mu, lo=lo, hi=hi, p=max(dl[2], hk[2]), k=k, tau2=t2)


SROWS = []
for c, g in SA.groupby("cohort", sort=False):
    Z = g[GC].T; Z.columns = g["sample"].values; nc = set(nonconst(Z)); y = g.response.to_numpy(float)
    for s, v in LIB.items():
        pres = [x for x in v["genes"] if x in nc]
        if len(pres) < 3 or len(pres) < 0.5 * len(v["genes"]):
            continue
        sc = wscore(Z, pres, v["w"][[v["genes"].index(x) for x in pres]])
        if not np.isfinite(sc).all():
            continue
        X = np.column_stack([np.ones(len(y)), std(sc)]); b, se = lfit(X, y, [1], True)[0]
        SROWS.append(dict(signature=s, cohort=c, n=len(y), responders=int(y.sum()), genes=len(pres), logOR=b, se=se))
SPC = pd.DataFrame(SROWS); SPC.to_csv(f"{O}/E40_single_arm_per_cohort.csv", index=False)
SP = {}
for s, g in SPC.groupby("signature", sort=False):
    if len(g) >= 10:
        SP[s] = dl_hk(g.logOR, g.se)
SP = pd.DataFrame(SP).T.add_prefix("sa_")
T = RP.join(SP, how="left"); T["family"] = [LIB[s]["family"] for s in T.index]; T["declared"] = [LIB[s]["dir"] for s in T.index]
T.to_csv(f"{O}/E40_concordance.csv")
print(f"signatures pooled in single-arm cohorts (>= 10 of 20): {SP.shape[0]} of {len(LIB)}; not pooled:",
      sorted(set(LIB) - set(SP.index)))

# ------------------------------------------------------------------ analyses
sig = lambda col: T[f"{col}_p"] < 0.05
n_int, n_icb, n_ctl = int(sig("int").sum()), int(sig("icb").sum()), int(sig("ctl").sum())
n_sa = int(sig("sa").sum())
print(f"\n1. P < 0.05 (either direction): interaction {n_int}; checkpoint arm alone {n_icb}; control arm alone {n_ctl}; "
      f"single-arm cohorts {n_sa} (of {SP.shape[0]})")
print("checkpoint arm alone x interaction:"); print(pd.crosstab(sig("icb").rename("checkpoint arm P<0.05"), sig("int").rename("interaction P<0.05")))
print("checkpoint arm alone x control arm:"); print(pd.crosstab(sig("icb").rename("checkpoint arm P<0.05"), sig("ctl").rename("control arm P<0.05")))
print("single-arm x interaction:"); print(pd.crosstab(sig("sa").rename("single-arm P<0.05"), sig("int").rename("interaction P<0.05")))
icb_only = T[sig("icb") & ~sig("int")]
print(f"significant in the checkpoint arm alone but not in the interaction: {len(icb_only)}; of these, control arm P < 0.05 "
      f"in the same direction: {int(((icb_only.ctl_p < 0.05) & (np.sign(icb_only.ctl_est) == np.sign(icb_only.icb_est))).sum())}")
TT = T.dropna(subset=["sa_est"])
rho = dict(icb=spearmanr(T.int_est, T.icb_est)[0], ctl=spearmanr(T.int_est, T.ctl_est)[0], sa=spearmanr(TT.int_est, TT.sa_est)[0],
           sa_icb=spearmanr(TT.icb_est, TT.sa_est)[0])
print("\n2. Spearman with the interaction: checkpoint arm alone {icb:.3f}; control arm {ctl:.3f}; single-arm log-OR {sa:.3f} "
      "(single-arm vs checkpoint arm alone {sa_icb:.3f})".format(**rho))
sa_sig = TT[TT.sa_p < 0.05]; same = sa_sig[(sa_sig.int_p < 0.05) & (np.sign(sa_sig.int_est) == np.sign(sa_sig.sa_est))]
frac = len(same) / len(sa_sig) if len(sa_sig) else np.nan
print(f"\n3. significant in single-arm cohorts: {len(sa_sig)}; with an interaction at P < 0.05 in the same direction: {len(same)} ({frac:.0%})")
print(T[["family", "int_est", "int_p", "icb_est", "icb_p", "ctl_est", "ctl_p", "sa_est", "sa_p", "sa_k"]].sort_values("int_est", ascending=False).round(3).to_string())

EXP = pd.DataFrame([
    dict(q="C1", expectation="more signatures at P < 0.05 in the checkpoint arm alone than in the interaction",
         result=f"{n_icb} vs {n_int}", held="yes" if n_icb > n_int else "no"),
    dict(q="C2", expectation="Spearman(single-arm log-OR, interaction) > 0 and < 0.8", result=f"{rho['sa']:.3f}",
         held="yes" if 0 < rho["sa"] < 0.8 else "no"),
    dict(q="C3", expectation="< 60% of single-arm-significant signatures have a same-direction interaction at P < 0.05",
         result=f"{len(same)}/{len(sa_sig)} ({frac:.0%})", held="yes" if frac < 0.6 else "no")])
EXP.to_csv(f"{O}/E40_expectations.csv", index=False); print(EXP.to_string(index=False))
