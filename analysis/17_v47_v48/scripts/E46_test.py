"""E46 - simple T-cell reference scores; choosing a gene in advance (leave one comparison out); the MSS/EBV-negative gastric
patients; the I-SPY2 control arms. Protocol: analysis/E46_prereg.md (hash checked). Setup: E34_trials.py up to its fitting stage,
as E42 and E45. Step 0: the 52 signatures and NetImmune minus CD8A reproduce E32_primary.csv; CD8A per gastric cohort reproduces
E44_B_per_cohort.csv. Writes results/E46_A_*.csv, E46_B_*.csv, E46_C_*.csv, E46_D_*.csv, E46_expectations.csv, results/E46_test.log."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import hashlib, sys
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h46 = hashlib.sha256(open(f"{A}/E46_prereg.md", "rb").read()).hexdigest()
assert _h46 == open(f"{A}/E46_prereg.sha256").read().split()[0], "E46 protocol hash mismatch"
_hw = hashlib.sha256(open(f"{A}/results/E37_frozen_weights.tsv", "rb").read()).hexdigest()
assert _hw == open(f"{A}/E37_frozen.sha256").read().split()[0], "frozen weights hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E46_step0_reproduction").replace("E34_coverage_check.csv", "E46_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
from scipy.stats import mannwhitneyu, rankdata


class Tee:
    def __init__(self, *f): self.f = f
    def write(self, s): [x.write(s) for x in self.f]
    def flush(self): [x.flush() for x in self.f]


sys.stdout = Tee(sys.__stdout__, open(f"{O}/E46_test.log", "w"))
print("E46 protocol sha256 verified:", _h46[:12], "| frozen weights:", _hw[:12])
NP = int(sys.argv[1]) if len(sys.argv) > 1 else 90
RENAL3 = ["IMM150", "CM025", "JR101"]; BREAST3 = ["ISPY2_pembro", "ISPY2_durva", "NeoTRIP"]
SETS = {"six (primary)": (ALL6, 4), "renal three": (RENAL3, 2), "breast three": (BREAST3, 2)}
E32P = pd.read_csv(f"{O}/E32_primary.csv"); E32P = E32P[E32P.set == "six (primary)"].set_index("signature")
R3G = ["CD8A", "PRF1", "GZMA"]; R13G = list(T_PANEL); DUP, AYERS = "TIS_Damotte", "Tcell-inflamed_Ayers"
ZQ = norm.ppf(0.975); PANEL = list(T_PANEL) + list(I_PANEL)

# ================================================================== A. reference scores
FW = pd.read_csv(f"{O}/E37_frozen_weights.tsv", sep="\t"); WG = dict(zip(FW.gene, FW.weight))
SC = {}
for k in ALL6:
    Z = TR[k]["Z"]; nc = set(nonconst(Z))
    for s in list(LIB) + ["CD8A", "NetImmune"]:
        if (s, k) in OBS:
            SC[(s, k)] = OBS[(s, k)][0]
    for g in PANEL:
        if g in nc and np.isfinite(Z.loc[g].to_numpy(float)).all() and g != "CD8A":
            SC[(g, k)] = Z.loc[g].to_numpy(float)
    for nm, gs in [("R3", R3G), ("R13", R13G)]:
        assert all(g in nc for g in gs), (nm, k)
        SC[(nm, k)] = wscore(Z, gs, np.ones(len(gs)))
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
print(len(tasks), "score x comparison tasks", flush=True)
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
dup = max(abs(R[(DUP, k)][0] - R[(AYERS, k)][0]) for k in ALL6); print(f"TIS_Damotte and Tcell-inflamed_Ayers: max per-comparison |diff| {dup:.1e}"); assert dup < 1e-10


def pooled(name, keys):
    ks = [k for k in keys if (name, k) in R]; r = pool_rule([R[(name, k)][0] for k in ks], [R[(name, k)][1] for k in ks])
    return dict(k=len(ks), est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"], k_pos=sum(R[(name, k)][0] > 0 for k in ks))


P33 = pd.read_csv(f"{O}/E33_pooled.csv"); G6 = P33[(P33.set == "six (primary)") & (P33.kind == "gene") & (P33.all_six.astype(str) == "True")].copy()
RP = []
for st, (keys, kmin) in SETS.items():
    for nm in ["R3", "R13", "CD8A", "CXCL10"]:
        row = dict(set=st, score=nm, **pooled(nm, keys))
        if st.startswith("six"):
            row["percentile_of_2082"] = 100 * (G6.est < row["est"]).mean(); row["rank_of_2082"] = int((G6.est > row["est"]).sum() + 1)
        RP.append(row)
RP = pd.DataFrame(RP); RP.to_csv(f"{O}/E46_A_reference_pooled.csv", index=False); print(RP.round(3).to_string(index=False))
pd.DataFrame([dict(score=nm, comparison=k, y=R[(nm, k)][0], se=R[(nm, k)][1]) for nm in ["R3", "R13", "CD8A", "CXCL10"] for k in ALL6]).to_csv(f"{O}/E46_A_reference_per_comparison.csv", index=False)
PA, CT = [], []
for st, (keys, kmin) in SETS.items():
    for ref in ["R3", "R13"]:
        rows = []
        for s in list(LIB) + ["NetImmune", "BN03"]:
            q = paired(s, ref, keys, kmin)
            if q is None:
                continue
            rows.append(dict(set=st, reference=ref, score="barriernet 0.3" if s == "BN03" else s, family=LIB[s]["family"] if s in LIB else "composite", **q))
        T = pd.DataFrame(rows); T["above"] = (T["diff"] > 0) & (T.p < 0.05); T["below"] = (T["diff"] < 0) & (T.p < 0.05); T["within010"] = (T.lo > -0.10) & (T.hi < 0.10); PA.append(T)
        L_ = T[T.score.isin(LIB) & (T.score != DUP)]
        CT.append(dict(set=st, reference=ref, n=len(L_), above=int(L_.above.sum()), below=int(L_.below.sum()), within010=int(L_.within010.sum()), median_mdd=L_.mdd.median(),
                       median_hi=L_.hi.median(), max_hi=L_.hi.max(), above_names="; ".join(L_[L_.above].sort_values("diff", ascending=False).score), below_names="; ".join(L_[L_.below].sort_values("diff").score)))
    for a_, b_ in [("R3", "CD8A"), ("R3", "CXCL10"), ("R13", "CD8A"), ("R13", "R3")]:
        q = paired(a_, b_, keys, kmin); CT.append(dict(set=st, reference=f"{a_} minus {b_}", **{x: q[x] for x in ("k", "diff", "lo", "hi", "p", "mdd")}))
PA = pd.concat(PA); PA.to_csv(f"{O}/E46_A_paired.csv", index=False); CT = pd.DataFrame(CT); CT.to_csv(f"{O}/E46_A_counts.csv", index=False)
print(CT.drop(columns=["above_names", "below_names"]).round(3).to_string(index=False)); print(CT[CT.above > 0][["set", "reference", "above_names"]].to_string(index=False)); print(CT[CT.below > 0][["set", "reference", "below_names"]].to_string(index=False))

# steadiness across tumor type
P25 = [g for g in PANEL if all((g, k) in R for k in ALL6)]; assert len(P25) == 25, len(P25)
INF = [s for s, v in LIB.items() if v["family"] == E32P.loc[AYERS, "family"] and s != DUP]; print("distinct inflamed signatures:", len(INF))
ST = []
for lab, names in [("panel gene", P25), ("inflamed signature", INF), ("reference", ["R3", "R13"])]:
    for nm in names:
        a_, b_ = pooled(nm, RENAL3), pooled(nm, BREAST3); ST.append(dict(group=lab, score=nm, renal=a_["est"], renal_k=a_["k"], breast=b_["est"], breast_k=b_["k"], D=abs(a_["est"] - b_["est"])))
ST = pd.DataFrame(ST); ST.to_csv(f"{O}/E46_A_steadiness.csv", index=False)
mg, ms = ST[ST.group == "panel gene"].D, ST[ST.group == "inflamed signature"].D; mw = mannwhitneyu(mg, ms, alternative="two-sided")
print(f"steadiness: median D panel genes {mg.median():.3f} (n {len(mg)}), inflamed signatures {ms.median():.3f} (n {len(ms)}), Mann-Whitney P {mw.pvalue:.3f}"); print(ST[ST.group == "reference"].round(3).to_string(index=False)); print(ST[ST.score.isin(["CD8A", "CXCL10", AYERS])].round(3).to_string(index=False))

# ================================================================== B. choosing a gene in advance
BR = []; BH = []
APP = {g: pooled(g, ALL6)["est"] for g in P25}; TOP = max(APP, key=APP.get); print(f"apparent best of the 25: {TOP} {APP[TOP]:+.3f}; of the 13 T-cell genes: {max(R13G, key=lambda g: APP[g])} {max(APP[g] for g in R13G):+.3f}")
for lab, cand in [("25 panel genes (primary)", P25), ("13 T-cell genes", R13G)]:
    ch = {}
    for h in ALL6:
        oth = [k for k in ALL6 if k != h]; est = {g: pooled(g, oth)["est"] for g in cand}; g_ = max(est, key=est.get); ch[h] = g_
        BH.append(dict(candidates=lab, held_out=h, gene=g_, est_other_five=est[g_], y_held_out=R[(g_, h)][0], se_held_out=R[(g_, h)][1], y_R3=R[("R3", h)][0], y_inflamed=R[(AYERS, h)][0]))
    r = pool_rule([R[(ch[h], h)][0] for h in ALL6], [R[(ch[h], h)][1] for h in ALL6]); ap = max(APP[g] for g in cand)
    row = dict(candidates=lab, est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"], apparent_best=ap, shrinkage=ap - r["est"])
    for nm, other in [("R3", "R3"), ("inflamed", AYERS)]:
        dd, vv = [], []
        for h in ALL6:
            db = R[(ch[h], h)][2] - R[(other, h)][2]; ok = np.isfinite(db); dd.append(R[(ch[h], h)][0] - R[(other, h)][0]); vv.append(np.var(db[ok], ddof=1))
        dd, vv = np.array(dd), np.array(vv); w = 1 / vv; D_ = np.sum(w * dd) / w.sum(); s_ = np.sqrt(1 / w.sum())
        row.update({f"minus_{nm}": D_, f"minus_{nm}_lo": D_ - 1.96 * s_, f"minus_{nm}_hi": D_ + 1.96 * s_, f"minus_{nm}_p": 2 * norm.sf(abs(D_ / s_))})
    BR.append(row)
BH = pd.DataFrame(BH); BR = pd.DataFrame(BR); BH.to_csv(f"{O}/E46_B_held_out.csv", index=False); BR.to_csv(f"{O}/E46_B_pooled.csv", index=False)
print(BH.round(3).to_string(index=False)); print(BR.round(3).T.to_string())

# ================================================================== C. MSS/EBV-negative gastric patients
GASTRIC = ["PRJEB25780", "PRJEB60680"]; LB = pd.read_csv(f"{O}/E44_labels.csv"); PC44 = pd.read_csv(f"{O}/E44_B_per_cohort.csv"); RAL = {"WARS1": "WARS"}


def logor(y, sc):
    try:
        return lfit(np.column_stack([np.ones(len(y)), std(sc)]), y, [1], True)[0]
    except Exception:
        return (np.nan, np.nan)


def fe(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float); w = 1 / se ** 2; m = np.sum(w * y) / w.sum(); s = np.sqrt(1 / w.sum())
    return dict(est=m, lo=m - ZQ * s, hi=m + ZQ * s, p=2 * norm.sf(abs(m / s)))


CR = []; d0 = 0
for c in GASTRIC:
    Z = pd.read_csv(f"{O}/E44_Z_{c}.csv.gz", index_col=0); Z = Z.rename(index={a: b for a, b in RAL.items() if a in Z.index and b not in Z.index})
    lb = LB[LB.cohort == c]; assert list(lb["sample"]) == list(Z.columns); y = lb.response.to_numpy(float); me = lb.msi_ebv.to_numpy(float); sub = me == 0
    nc = set(nonconst(Z)); S = scores(Z); OB = {}
    for s, v in LIB.items():
        pres = [g for g in v["genes"] if g in nc]
        if len(pres) >= 3 and len(pres) >= 0.5 * len(v["genes"]):
            sc = wscore(Z, pres, v["w"][[v["genes"].index(g) for g in pres]])
            if np.isfinite(sc).all():
                OB[("sig", s)] = sc
    for s in ("ESTIMATE", "TMEscore_z", "Immune", "NetImmune"):
        if s in S and np.isfinite(S[s].to_numpy(float)).all():
            OB[("sig", s)] = S[s].to_numpy(float)
    pres = [g for g in WG if g in nc and np.isfinite(Z.loc[g].to_numpy(float)).all()]
    OB[("sig", "barriernet 0.3")] = wscore(Z, pres, np.array([WG[g] for g in pres]))
    for g in PANEL:
        if g in nc:
            OB[("gene", g)] = Z.loc[g].to_numpy(float)
    for nm, gs in [("R3", R3G), ("R13", R13G)]:
        g_ = [g for g in gs if g in nc]; assert len(g_) == len(gs), (c, nm); OB[("ref", nm)] = wscore(Z, g_, np.ones(len(g_)))
    d0 = max(d0, abs(logor(y, OB[("gene", "CD8A")])[0] - PC44[(PC44.cohort == c) & (PC44.name == "CD8A") & (PC44.kind == "gene")].logOR.iloc[0]))
    print(f"{c}: MSS/EBV-negative n {int(sub.sum())}, responders {int(y[sub].sum())}; scored {sum(1 for x in OB if x[0] == 'sig' and x[1] in LIB)}/52")
    for kk, sc in OB.items():
        b, se = logor(y[sub], sc[sub]); CR.append(dict(cohort=c, kind=kk[0], name=kk[1], n=int(sub.sum()), responders=int(y[sub].sum()), logOR=b, se=se, p=2 * norm.sf(abs(b / se))))
print(f"step 0: CD8A in all patients reproduces E44 (max |diff| {d0:.1e})"); assert d0 < 1e-6
CR = pd.DataFrame(CR); CR.to_csv(f"{O}/E46_C_per_cohort.csv", index=False)
CP = []
for (kd, nm), q in CR.groupby(["kind", "name"], sort=False):
    q = q[np.isfinite(q.logOR) & np.isfinite(q.se)]
    if len(q) < 2:
        continue
    f = fe(q.logOR, q.se); r = pool_rule(q.logOR.values, q.se.values); one = q[q.cohort == "PRJEB60680"].iloc[0]
    CP.append(dict(kind=kd, name=nm, k=len(q), est=f["est"], lo=f["lo"], hi=f["hi"], p=f["p"], rule_est=r["est"], rule_lo=r["lo"], rule_hi=r["hi"], rule_p=r["p"], est_PRJEB60680=one.logOR, se_PRJEB60680=one.se, p_PRJEB60680=one.p))
CP = pd.DataFrame(CP); CP.to_csv(f"{O}/E46_C_pooled.csv", index=False)
SG = CP[(CP.kind == "sig") & CP.name.isin(LIB) & (CP.name != DUP)]
CC = dict(n_distinct=len(SG), positive=int((SG.est > 0).sum()), p05_fixed=int((SG.p < 0.05).sum()), p05_fixed_positive=int(((SG.p < 0.05) & (SG.est > 0)).sum()), p05_rule=int((SG.rule_p < 0.05).sum()),
          p05_PRJEB60680=int((SG.p_PRJEB60680 < 0.05).sum()), positive_PRJEB60680=int((SG.est_PRJEB60680 > 0).sum()), names_fixed="; ".join(SG[SG.p < 0.05].sort_values("est", ascending=False).name))
pd.DataFrame([CC]).to_csv(f"{O}/E46_C_counts.csv", index=False); print(CC)
print(CP[CP.name.isin(["CD8A", "CXCL10", "CXCL13", "barriernet 0.3", "Immune", "NetImmune", "ESTIMATE", "TMEscore_z", "R3", "R13", AYERS])].round(3).to_string(index=False))

# ================================================================== D. I-SPY2 control arms (reads no outcome)
assert "X" in TR["ISPY2_durva"] and "X" in TR["ISPY2_pembro"]


def ranks(M):
    return np.apply_along_axis(rankdata, 0, M)


Xd, Xp = TR["ISPY2_durva"]["X"], TR["ISPY2_pembro"]["X"]; Xd = Xd[~Xd.index.duplicated()]; Xp = Xp[~Xp.index.duplicated()]; g_ = Xd.index.intersection(Xp.index); assert len(g_) >= 200
ad, ap = TR["ISPY2_durva"]["arm"], TR["ISPY2_pembro"]["arm"]
idd, idp = Xd.columns.astype(str), Xp.columns.astype(str); DR, FL = [], []
Ra = ranks(np.nan_to_num(Xd.loc[g_].to_numpy(float))); Rb = ranks(np.nan_to_num(Xp.loc[g_].to_numpy(float))); Ra = (Ra - Ra.mean(0)) / Ra.std(0); Rb = (Rb - Rb.mean(0)) / Rb.std(0); C = Ra.T @ Rb / len(g_)
for lab, ia, ib in [("controls of durva vs controls of pembro", ad == 0, ap == 0), ("all of durva vs all of pembro", np.ones(len(ad), bool), np.ones(len(ap), bool))]:
    for direction, CC_, ida, idb in [("durva to pembro", C[np.ix_(ia, ib)], idd[ia], idp[ib]), ("pembro to durva", C[np.ix_(ia, ib)].T, idp[ib], idd[ia])]:
        srt = -np.sort(-CC_, axis=1); best, gp = srt[:, 0], srt[:, 0] - srt[:, 1]
        DR.append(dict(set=lab, direction=direction, shared_genes=len(g_), n_query=CC_.shape[0], n_target=CC_.shape[1], id_overlap=len(set(ida) & set(idb)), median_best_r=np.median(best), max_best_r=best.max(), median_gap=np.median(gp), max_gap=gp.max(), n_flagged=int((gp > 0.03).sum())))
        for i in np.where(gp > 0.03)[0]:
            FL.append(dict(set=lab, direction=direction, sample=ida[i], best_r=best[i], gap=gp[i], partner=idb[int(np.argmax(CC_[i]))]))
DR = pd.DataFrame(DR); DR.to_csv(f"{O}/E46_D_ispy2_twins.csv", index=False); pd.DataFrame(FL, columns=["set", "direction", "sample", "best_r", "gap", "partner"]).to_csv(f"{O}/E46_D_ispy2_flags.csv", index=False); print(DR.round(3).to_string(index=False)); print("flagged:", len(FL))
nctrl = int(DR[DR.set.str.startswith("controls")].n_flagged.sum())
if nctrl:
    fl = {r["sample"] for r in FL if r["set"].startswith("controls") and r["direction"] == "durva to pembro"}; keep = np.array([c not in fl for c in idd]); dq = take(TR["ISPY2_durva"], np.where(keep)[0]); SENS = []
    for nm in ["CD8A", "R3", AYERS]:
        ys = [R[(nm, k)][0] if k != "ISPY2_durva" else m_int(dq, SC[(nm, k)][keep])[0] for k in ALL6]; ss = [R[(nm, k)][1] if k != "ISPY2_durva" else m_int(dq, SC[(nm, k)][keep])[1] for k in ALL6]; r = pool_rule(ys, ss); SENS.append(dict(score=nm, removed=int((~keep).sum()), est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"]))
    SENS = pd.DataFrame(SENS); SENS.to_csv(f"{O}/E46_D_sensitivity.csv", index=False); print(SENS.round(3).to_string(index=False))

# ================================================================== expectations
r3 = RP[(RP.score == "R3")].set_index("set"); c8 = RP[RP.score == "CD8A"].set_index("set"); a2 = CT[(CT.set == "six (primary)") & (CT.reference == "R3")].iloc[0]; b0 = BR.iloc[0]; c8g = CP[CP.name == "CD8A"].iloc[0]
EX = [("A1", "R3 favors the checkpoint arm in the six comparisons at P < 0.05", r3.loc["six (primary)", "est"] > 0 and r3.loc["six (primary)", "p"] < 0.05, f"{r3.loc['six (primary)', 'est']:+.3f}, P {r3.loc['six (primary)', 'p']:.3f}"),
      ("A2", "at most 2 of the 51 distinct signatures significantly above R3", a2.above <= 2, f"{int(a2.above)} above, {int(a2.below)} below"),
      ("A3", "R3 larger breast and smaller renal interaction than CD8A", r3.loc["breast three", "est"] > c8.loc["breast three", "est"] and r3.loc["renal three", "est"] < c8.loc["renal three", "est"], f"breast {r3.loc['breast three', 'est']:+.3f} vs {c8.loc['breast three', 'est']:+.3f}; renal {r3.loc['renal three', 'est']:+.3f} vs {c8.loc['renal three', 'est']:+.3f}"),
      ("A4", "median D of inflamed signatures below that of panel genes", ms.median() < mg.median(), f"{ms.median():.3f} vs {mg.median():.3f}; Mann-Whitney P {mw.pvalue:.3f}"),
      ("B1", "held-out interaction of the gene chosen in advance below the apparent best", b0.est < b0.apparent_best, f"{b0.est:+.3f} vs {b0.apparent_best:+.3f}"),
      ("B2", "gene chosen in advance does not significantly exceed R3", not (b0.minus_R3 > 0 and b0.minus_R3_p < 0.05), f"{b0.minus_R3:+.3f}, P {b0.minus_R3_p:.3f}"),
      ("C1", "CD8A positive fixed-effect log-odds in the 73 patients, P >= 0.05", c8g.est > 0 and c8g.p >= 0.05, f"{c8g.est:+.3f}, P {c8g.p:.3f}"),
      ("C2", "fewer than 5 of 51 distinct signatures at P < 0.05 by fixed effect in the 73", CC["p05_fixed"] < 5, f"{CC['p05_fixed']}"),
      ("D1", "no I-SPY2 control flagged as shared", nctrl == 0, f"{nctrl} flagged")]
EX = pd.DataFrame([dict(id=i, expectation=e, held=bool(h), result=v) for i, e, h, v in EX]); EX.to_csv(f"{O}/E46_expectations.csv", index=False); print(EX.to_string(index=False)); print(f"{int(EX.held.sum())} of {len(EX)} expectations held")
