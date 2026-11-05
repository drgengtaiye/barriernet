"""E45 step 3 - the frozen reference gene in the randomized comparisons and the gastroesophageal cohorts; the coexpression-matched
random reference; the renal arms and vascular genes. Protocol: analysis/E45_prereg.md; frozen genes: results/E45_reference_genes.tsv
(both hashes checked). Setup: E34_trials.py up to its fitting stage, as E42. Step 0: the 52 signatures and NetImmune minus CD8A
reproduce E32_primary.csv; CD8A per gastroesophageal cohort reproduces E44_B_per_cohort.csv.
Writes results/E45_A_*.csv, E45_B_*.csv, E45_C_*.csv, E45_expectations.csv and results/E45_test.log."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import hashlib, sys, zlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h45 = hashlib.sha256(open(f"{A}/E45_prereg.md", "rb").read()).hexdigest()
assert _h45 == open(f"{A}/E45_prereg.sha256").read().split()[0], "E45 protocol hash mismatch"
_hf = hashlib.sha256(open(f"{A}/results/E45_reference_genes.tsv", "rb").read()).hexdigest()
assert _hf == open(f"{A}/E45_frozen.sha256").read().split()[0], "E45 frozen genes hash mismatch"
_hw = hashlib.sha256(open(f"{A}/results/E37_frozen_weights.tsv", "rb").read()).hexdigest()
assert _hw == open(f"{A}/E37_frozen.sha256").read().split()[0], "frozen weights hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E45_step0_reproduction").replace("E34_coverage_check.csv", "E45_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
from scipy.stats import fisher_exact


class Tee:
    def __init__(self, *f): self.f = f
    def write(self, s): [x.write(s) for x in self.f]
    def flush(self): [x.flush() for x in self.f]


sys.stdout = Tee(sys.__stdout__, open(f"{O}/E45_test.log", "w"))
print("E45 protocol sha256 verified:", _h45[:12], "| frozen genes:", _hf[:12], "| frozen weights:", _hw[:12])
NP = int(sys.argv[1]) if len(sys.argv) > 1 else 90
REFG = pd.read_csv(f"{O}/E45_reference_genes.tsv", sep="\t").set_index("reference").gene
G1, G2 = REFG["G1"], REFG["G2"]; REFS = list(dict.fromkeys([G1, G2])); print("frozen references:", REFS)
RENAL3 = ["IMM150", "CM025", "JR101"]; BREAST3 = ["ISPY2_pembro", "ISPY2_durva", "NeoTRIP"]
SETS = {"six (primary)": (ALL6, 4), "renal three": (RENAL3, 2), "breast three": (BREAST3, 2)}
E32P = pd.read_csv(f"{O}/E32_primary.csv"); E32P = E32P[E32P.set == "six (primary)"].set_index("signature")
ANGIO = ["VEGFA", "KDR", "ESM1", "PECAM1", "ANGPTL4", "CD34"]; VASC = "TUMOR-VASCULATURE-UP_Lu"
ZQ = norm.ppf(0.975)


def zm(Z, genes):
    g = [x for x in genes if x in Z.index]
    return Z.loc[g].mean(0).to_numpy(float) if len(g) >= 3 else None


# ================================================================== A. reference gene in the six comparisons
FW = pd.read_csv(f"{O}/E37_frozen_weights.tsv", sep="\t"); WG = dict(zip(FW.gene, FW.weight))
SC = {}
for k in ALL6:
    Z = TR[k]["Z"]; nc = set(nonconst(Z))
    for s in list(LIB) + ["CD8A", "NetImmune"]:
        if (s, k) in OBS:
            SC[(s, k)] = OBS[(s, k)][0]
    for g in REFS:
        assert g in nc, (g, k); SC[(g, k)] = Z.loc[g].to_numpy(float)
    pres = [g for g in WG if g in nc and np.isfinite(Z.loc[g].to_numpy(float)).all()]
    if len(pres) >= 3 and len(pres) >= 0.5 * len(WG):
        SC[("BN03", k)] = wscore(Z, pres, np.array([WG[g] for g in pres]))
BIDX = {}
for k in ALL6:
    rng_ = np.random.default_rng([3032, KIDX[k]]); a = TR[k]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
    BIDX[k] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]


def run(task):
    name, k = task; sc = SC[(name, k)]; d = TR[k]; y, se = m_int(d, sc); bs = []
    for idx in BIDX[k]:
        try:
            bs.append(m_int(take(d, idx), sc[idx])[0])
        except Exception:
            bs.append(np.nan)
    return y, se, np.array(bs)


tasks = list(SC)
with mp.get_context("fork").Pool(NP) as pool:
    res = pool.map(run, tasks, chunksize=1)
R = dict(zip(tasks, res))


def paired(a, b, keys, kmin=4):
    ks = [k for k in keys if (a, k) in R and (b, k) in R]
    if len(ks) < kmin:
        return None
    dd, vv = [], []
    for k in ks:
        ya, _, ba = R[(a, k)]; yb, _, bb = R[(b, k)]; db = ba - bb; ok = np.isfinite(db)
        dd.append(ya - yb); vv.append(np.var(db[ok], ddof=1))
    dd, vv = np.array(dd), np.array(vv); w = 1 / vv; D_ = np.sum(w * dd) / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(k=len(ks), diff=D_, se=se, lo=D_ - 1.96 * se, hi=D_ + 1.96 * se, p=2 * norm.sf(abs(D_ / se)), mdd=2.80 * se)


gap = max(abs(paired(s, "CD8A", ALL6)["diff"] - E32P.loc[s, "diff_vs_CD8A"]) for s in list(LIB) + ["NetImmune"])
print(f"step 0: 52 signatures + NetImmune minus CD8A reproduce E32 (max |diff| {gap:.1e})"); assert gap < 1e-8
P33 = pd.read_csv(f"{O}/E33_pooled.csv"); G6 = P33[(P33.set == "six (primary)") & (P33.kind == "gene") & (P33.all_six.astype(str) == "True")].copy()
RP = []
for st, (keys, kmin) in SETS.items():
    for g in ["CD8A"] + REFS:
        ks = [k for k in keys if (g, k) in R]; r = pool_rule([R[(g, k)][0] for k in ks], [R[(g, k)][1] for k in ks])
        row = dict(set=st, gene=g, k=len(ks), est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"], k_pos=sum(R[(g, k)][0] > 0 for k in ks))
        if st.startswith("six"):
            row["percentile_of_2082"] = 100 * (G6.est < r["est"]).mean(); row["rank_of_2082"] = int((G6.est > r["est"]).sum() + 1)
        RP.append(row)
RP = pd.DataFrame(RP); RP.to_csv(f"{O}/E45_A_reference_pooled.csv", index=False); print(RP.round(3).to_string(index=False))
perc = []
for g in ["CD8A"] + REFS:
    for k in ALL6:
        perc.append(dict(gene=g, comparison=k, y=R[(g, k)][0], se=R[(g, k)][1]))
pd.DataFrame(perc).to_csv(f"{O}/E45_A_reference_per_comparison.csv", index=False)
PA, CT = [], []
for st, (keys, kmin) in SETS.items():
    for ref in ["CD8A"] + REFS:
        rows = []
        for s in list(LIB) + ["NetImmune", "BN03"]:
            q = paired(s, ref, keys, kmin)
            if q is None:
                continue
            rows.append(dict(set=st, reference=ref, score="barriernet 0.3" if s == "BN03" else s, family=LIB[s]["family"] if s in LIB else "composite", **q))
        T = pd.DataFrame(rows); T["above"] = (T["diff"] > 0) & (T.p < 0.05); T["below"] = (T["diff"] < 0) & (T.p < 0.05); T["within010"] = (T.lo > -0.10) & (T.hi < 0.10); PA.append(T)
        L_ = T[T.score.isin(LIB)]
        CT.append(dict(set=st, reference=ref, n=len(L_), above=int(L_.above.sum()), below=int(L_.below.sum()), within010=int(L_.within010.sum()), median_mdd=L_.mdd.median(),
                       median_hi=L_.hi.median(), max_hi=L_.hi.max(), above_names="; ".join(L_[L_.above].sort_values("diff", ascending=False).score)))
    for ref in REFS:
        q = paired(ref, "CD8A", keys, kmin); CT.append(dict(set=st, reference=f"{ref} minus CD8A", **{x: q[x] for x in ("k", "diff", "lo", "hi", "p", "mdd")}))
PA = pd.concat(PA); PA.to_csv(f"{O}/E45_A_paired.csv", index=False); CT = pd.DataFrame(CT); CT.to_csv(f"{O}/E45_A_counts.csv", index=False)
print(CT.drop(columns=["above_names"]).round(3).to_string(index=False)); print(CT[CT.above_names.notna() & (CT.above > 0)][["set", "reference", "above_names"]].to_string(index=False))

# ================================================================== A. gastroesophageal cohorts
COH = ["PRJEB25780", "PRJEB60680", "PERFECT"]; GASTRIC = COH[:2]
LB = pd.read_csv(f"{O}/E44_labels.csv"); PC44 = pd.read_csv(f"{O}/E44_B_per_cohort.csv")


def logor(y, sc, x=None):
    cols = [np.ones(len(y))] + ([x] if x is not None else []) + [std(sc)]
    try:
        return lfit(np.column_stack(cols), y, [len(cols) - 1], True)[0]
    except Exception:
        return (np.nan, np.nan)


def fe(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float); w = 1 / se ** 2; m = np.sum(w * y) / w.sum(); s = np.sqrt(1 / w.sum())
    return dict(est=m, lo=m - ZQ * s, hi=m + ZQ * s, p=2 * norm.sf(abs(m / s)), se=s)


GRows, GB = [], {}
for i, c in enumerate(COH):
    Z = pd.read_csv(f"{O}/E44_Z_{c}.csv.gz", index_col=0); lb = LB[LB.cohort == c]; assert list(lb["sample"]) == list(Z.columns)
    y = lb.response.to_numpy(float); me = lb.msi_ebv.to_numpy(float)
    rng_ = np.random.default_rng([4045, i]); i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
    BI = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(1000)]
    for g in ["CD8A"] + REFS:
        sc = Z.loc[g].to_numpy(float); b, se = logor(y, sc); row = dict(cohort=c, gene=g, n=len(y), responders=int(y.sum()), logOR=b, se=se, p=2 * norm.sf(abs(b / se)))
        if c in GASTRIC:
            ba, sa = logor(y, sc, me); row.update(logOR_msiebv=ba, se_msiebv=sa)
        GRows.append(row); GB[(c, g)] = np.array([logor(y[idx], sc[idx])[0] for idx in BI])
GP_ = pd.DataFrame(GRows); GP_.to_csv(f"{O}/E45_A_gastric_per_cohort.csv", index=False); print(GP_.round(3).to_string(index=False))
ck = PC44[(PC44.name == "CD8A") & (PC44.kind == "gene")].set_index("cohort").logOR
d0 = max(abs(GP_[(GP_.cohort == c) & (GP_.gene == "CD8A")].logOR.iloc[0] - ck[c]) for c in COH); print(f"step 0: CD8A per cohort reproduces E44 (max |diff| {d0:.1e})"); assert d0 < 1e-6
GPOOL = []
for lab, cs in [("three gastroesophageal (primary)", COH), ("two gastric", GASTRIC)]:
    for g in ["CD8A"] + REFS:
        q = GP_[(GP_.gene == g) & GP_.cohort.isin(cs)]; f = fe(q.logOR, q.se); r = pool_rule(q.logOR.values, q.se.values)
        row = dict(set=lab, gene=g, est=f["est"], lo=f["lo"], hi=f["hi"], p=f["p"], k_pos=int((q.logOR > 0).sum()), cons_est=r["est"], cons_lo=r["lo"], cons_hi=r["hi"], cons_p=r["p"])
        if cs == GASTRIC:
            fa = fe(q.logOR_msiebv, q.se_msiebv); row.update(est_msiebv=fa["est"], lo_msiebv=fa["lo"], hi_msiebv=fa["hi"], p_msiebv=fa["p"])
        if g != "CD8A":
            dd, ss = [], []
            for c in cs:
                db = GB[(c, g)] - GB[(c, "CD8A")]; ok = np.isfinite(db)
                dd.append(q[q.cohort == c].logOR.iloc[0] - GP_[(GP_.cohort == c) & (GP_.gene == "CD8A")].logOR.iloc[0]); ss.append(np.std(db[ok], ddof=1))
            fd = fe(dd, ss); row.update(diff_vs_CD8A=fd["est"], diff_lo=fd["lo"], diff_hi=fd["hi"], diff_p=fd["p"])
        GPOOL.append(row)
GPOOL = pd.DataFrame(GPOOL); GPOOL.to_csv(f"{O}/E45_A_gastric_pooled.csv", index=False); print(GPOOL.round(3).to_string(index=False))

# ================================================================== B. coexpression-matched reference
print("\nB. coexpression-matched random reference")
UNIV, CORR = {}, {}
for k in ALL6:
    Z = TR[k]["Z"]; g = [x for x in nonconst(Z) if np.isfinite(Z.loc[x].to_numpy(float)).all()]; U = Z.loc[g].to_numpy(float)
    UNIV[k] = U; X = (U - U.mean(1, keepdims=True)); X = X / np.sqrt((X ** 2).sum(1, keepdims=True))
    CORR[k] = (X @ X.T).astype(np.float32); print(f"  {k}: universe {U.shape[0]} genes, {U.shape[1]} samples", flush=True)
B_X = 1000; MAXATT = 500000


def coex_task(task):
    s, k = task; _, pres, w = OBS[(s, k)]; n = len(pres); Zk = TR[k]["Z"].loc[pres].T.corr().to_numpy(float)
    c0 = (np.nansum(Zk) - np.sum(np.isfinite(np.diag(Zk)))) / (np.sum(np.isfinite(Zk)) - n)
    U, C, d = UNIV[k], CORR[k], TR[k]; G = U.shape[0]
    rng = np.random.default_rng([4545, KIDX[k], zlib.crc32(s.encode())]); w = np.array(w, float)
    grid = sorted(set(min(G, q) for q in (n, 2 * n, 5 * n, 10 * n, 20 * n, 50 * n, 100 * n, G)))
    tol, att, acc, att2, out, cs = 0.02, 0, 0, 0, [], []
    while acc < B_X and att < 2 * MAXATT:
        if tol == 0.02 and att >= MAXATT:
            tol = 0.05; att2 = acc
        att += 1; seed = rng.integers(G); Q = grid[rng.integers(len(grid))]
        if Q < G:
            nb = np.argpartition(-C[seed], Q - 1)[:Q]; nb = nb[nb != seed]
        else:
            nb = np.delete(np.arange(G), seed)
        if len(nb) < n - 1:
            continue
        idx = np.concatenate([[seed], rng.choice(nb, n - 1, replace=False)]); cc = (C[np.ix_(idx, idx)].sum() - n) / (n * (n - 1))
        if abs(cc - c0) > tol:
            continue
        acc += 1; cs.append(cc); ww = w[rng.permutation(n)]
        try:
            out.append(m_int(d, U[idx].T @ ww))
        except Exception:
            out.append((np.nan, np.nan))
    while len(out) < B_X:
        out.append((np.nan, np.nan))
    return dict(signature=s, comparison=k, n=n, c_obs=c0, attempts=att, accepted=acc, accepted_at_002=(acc if att2 == 0 and tol == 0.02 else att2), tol_final=tol, c_ref_mean=float(np.mean(cs)) if cs else np.nan), np.array(out[:B_X])


xt = [(s, k) for s in LIB for k in ALL6 if (s, k) in OBS]
with mp.get_context("fork").Pool(NP) as pool:
    xres = pool.map(coex_task, xt, chunksize=1)
XPC = pd.DataFrame([r[0] for r in xres]); XPC.to_csv(f"{O}/E45_B_per_comparison.csv", index=False); XD = {t: r[1] for t, r in zip(xt, xres)}
print(f"  draws: {len(xt)} signature x comparison; at +/-0.02 throughout {(XPC.tol_final == 0.02).sum()}; widened {(XPC.tol_final > 0.02).sum()}; below 1,000 accepted {(XPC.accepted < B_X).sum()}; "
      f"failed fits {sum(int((~np.isfinite(v[:, 0])).sum()) for v in XD.values())}")


def coex_pool(s):
    ks = [k for k in ALL6 if (s, k) in XD]; D = np.stack([XD[(s, k)] for k in ks], 1); dr = float(LIB[s]["dir"])
    est = np.array([pool_rule(D[b, :, 0], D[b, :, 1])["est"] for b in range(B_X)]); obs = E32P.loc[s, "est"]; ok = np.isfinite(est)
    q = XPC[XPC.signature == s]
    return dict(signature=s, family=LIB[s]["family"], declared_dir=dr, k=len(ks), est=obs, emp_p_size=E32P.loc[s, "emp_p_declared"], ref_median_size=E32P.loc[s, "ref_median"],
                emp_p_coex=(1 + np.sum(dr * est[ok] >= dr * obs)) / (ok.sum() + 1), ref_median_coex=np.median(est[ok]), ref_p95_coex=np.percentile(dr * est[ok], 95) * dr, c_obs_mean=q.c_obs.mean(), draws=int(ok.sum()), widened=int((q.tol_final > 0.02).sum()))


with mp.get_context("fork").Pool(NP) as pool:
    XR = pd.DataFrame(pool.map(coex_pool, list(LIB), chunksize=1))
XR.to_csv(f"{O}/E45_B_coexpression.csv", index=False)
assert int((XR.emp_p_size < 0.05).sum()) == 26
XC = []
for fam, g in list(XR.groupby("family")) + [("all", XR)]:
    XC.append(dict(family=fam, n=len(g), exceed_size=int((g.emp_p_size < 0.05).sum()), exceed_coex=int((g.emp_p_coex < 0.05).sum()), median_ref_size=g.ref_median_size.median(), median_ref_coex=g.ref_median_coex.median()))
XC = pd.DataFrame(XC); XC.to_csv(f"{O}/E45_B_counts.csv", index=False); print(XC.round(3).to_string(index=False))
print(XR[XR.emp_p_coex < 0.05].sort_values("emp_p_coex")[["signature", "family", "est", "emp_p_size", "emp_p_coex", "ref_median_coex", "c_obs_mean", "widened"]].round(3).to_string(index=False))

# ================================================================== C. renal arms and vascular genes
print("\nC. renal arms and vascular genes")
def csc(name, k):
    Z = TR[k]["Z"]
    if name == "Angio":
        return zm(Z, ANGIO)
    if name == "Immune":
        return TR[k]["S"]["Immune"].to_numpy(float)
    if (name, k) in OBS:
        return OBS[(name, k)][0]
    return Z.loc[name].to_numpy(float) if name in Z.index else None


CW = []
names = list(dict.fromkeys(["CD8A"] + REFS + ["Immune", "Angio", VASC]))
for nm in names:
    for k in ALL6:
        sc = csc(nm, k)
        if sc is None or not np.isfinite(sc).all():
            continue
        d = TR[k]; iy, ise = m_int(d, sc); ay, ase = m_within(d, sc, 1); cy, cse = m_within(d, sc, 0)
        CW.append(dict(score=nm, comparison=k, int_y=iy, int_se=ise, icb_y=ay, icb_se=ase, ctl_y=cy, ctl_se=cse))
CW = pd.DataFrame(CW); CW.to_csv(f"{O}/E45_C_within_arm.csv", index=False); print(CW[CW.comparison.isin(RENAL3)].round(3).to_string(index=False))
CP = []
for st, keys in [("renal three", RENAL3), ("breast three", BREAST3), ("six (primary)", ALL6)]:
    for nm in names:
        q = CW[(CW.score == nm) & CW.comparison.isin(keys)]
        row = dict(set=st, score=nm, k=len(q))
        if len(q) >= 2:
            for c in ("int", "icb", "ctl"):
                r = pool_rule(q[f"{c}_y"].values, q[f"{c}_se"].values); row.update({f"{c}_est": r["est"], f"{c}_lo": r["lo"], f"{c}_hi": r["hi"], f"{c}_p": r["p"]})
        CP.append(row)
CP = pd.DataFrame(CP); CP.to_csv(f"{O}/E45_C_pooled.csv", index=False); print(CP.round(3).to_string(index=False))
VG = set(LIB[VASC]["genes"]); inlist = G6.name.isin(VG); neg = (G6.p < 0.05) & (G6.est < 0); pos = (G6.p < 0.05) & (G6.est > 0); assert neg.sum() == 102 and pos.sum() == 101
EN = []
for lab, m in [("control-favoring", neg), ("checkpoint-favoring", pos)]:
    a, b, c_, d_ = int((m & inlist).sum()), int((m & ~inlist).sum()), int((~m & inlist).sum()), int((~m & ~inlist).sum()); o_, p_ = fisher_exact([[a, b], [c_, d_]], alternative="two-sided")
    EN.append(dict(side=lab, vascular_in=a, other_in=b, vascular_out=c_, other_out=d_, odds_ratio=o_, p=p_, genes=" ".join(sorted(G6[m & inlist].name))))
EN = pd.DataFrame(EN); EN.to_csv(f"{O}/E45_C_enrichment.csv", index=False); print(f"vascular genes measured among 2,082: {int(inlist.sum())}"); print(EN.round(4).to_string(index=False))

# ================================================================== expectations
six = RP[RP.set.str.startswith("six")].set_index("gene"); ctG1 = CT[(CT.set.str.startswith("six")) & (CT.reference == G1)].iloc[0]; ctG2 = CT[(CT.set.str.startswith("six")) & (CT.reference == G2)].iloc[0]
gG1 = GPOOL[(GPOOL.set.str.startswith("three")) & (GPOOL.gene == G1)].iloc[0]; xa = XC[XC.family == "all"].iloc[0]; xi = XC[XC.family.str.startswith("Inflamed")].iloc[0]
cdr = CP[(CP.set == "renal three") & (CP.score == "CD8A")].iloc[0]; agr = CP[(CP.set == "renal three") & (CP.score == "Angio")].iloc[0]; en = EN.set_index("side")
X = [("A1", "G1 has a positive pooled interaction in the six comparisons", f"{G1} {six.loc[G1, 'est']:+.3f} ({six.loc[G1, 'lo']:+.3f} to {six.loc[G1, 'hi']:+.3f}), P = {six.loc[G1, 'p']:.3g}", six.loc[G1, "est"] > 0),
     ("A2", "at most 2 of the 52 signatures are significantly above G1", f"{ctG1.above} above, {ctG1.below} below", ctG1.above <= 2)]
if G2 != G1:
    X.append(("A3", "at most 2 of the 52 signatures are significantly above G2", f"{G2}: {ctG2.above} above, {ctG2.below} below", ctG2.above <= 2))
X += [("A4", "G1 has a positive pooled log-odds in the three gastroesophageal cohorts (fixed effect)", f"{gG1.est:+.3f} ({gG1.lo:+.3f} to {gG1.hi:+.3f}), P = {gG1.p:.3g}", gG1.est > 0),
      ("B1", "fewer than 26 of the 52 exceed the coexpression-matched reference", f"{xa.exceed_coex} of 52", xa.exceed_coex < 26),
      ("B2", "at least 5 of the 20 inflamed signatures still exceed it", f"{xi.exceed_coex} of 20", xi.exceed_coex >= 5),
      ("C1", "CD8A has a negative pooled control-arm slope over the three renal comparisons", f"{cdr.ctl_est:+.3f} ({cdr.ctl_lo:+.3f} to {cdr.ctl_hi:+.3f}), P = {cdr.ctl_p:.3g}", cdr.ctl_est < 0),
      ("C2", "the Angio score has a negative pooled interaction in the three renal comparisons", f"{agr.int_est:+.3f} ({agr.int_lo:+.3f} to {agr.int_hi:+.3f}), P = {agr.int_p:.3g}", agr.int_est < 0),
      ("C3", "TUMOR-VASCULATURE-UP_Lu genes enriched among control-favoring (P < 0.05), not among checkpoint-favoring", f"control OR {en.loc['control-favoring', 'odds_ratio']:.2f}, P = {en.loc['control-favoring', 'p']:.3g}; checkpoint OR {en.loc['checkpoint-favoring', 'odds_ratio']:.2f}, P = {en.loc['checkpoint-favoring', 'p']:.3g}",
       (en.loc['control-favoring', 'p'] < 0.05) and (en.loc['control-favoring', 'odds_ratio'] > 1) and not ((en.loc['checkpoint-favoring', 'p'] < 0.05) and (en.loc['checkpoint-favoring', 'odds_ratio'] > 1)))]
XP = pd.DataFrame([dict(q=q, expectation=e, result=r, held="yes" if h else "no") for q, e, r, h in X]); XP.to_csv(f"{O}/E45_expectations.csv", index=False); print(XP.to_string(index=False)); print("expectations held:", (XP.held == "yes").sum(), "of", len(XP))
