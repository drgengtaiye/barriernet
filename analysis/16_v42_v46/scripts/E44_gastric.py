"""E44 step 2 - gastric and gastroesophageal applications of the calibration framework.

Protocol: analysis/E44_prereg.md (hash checked); frozen weights results/E37_frozen_weights.tsv (hash checked). Setup:
E34_trials.py up to its fitting stage (E32 library, engine and helpers), executed unchanged as in E37_test.py except that
its step-0 and coverage tables are written under E44 names. Cohort matrices from E44_export.R (step 0 passed there).
A: published gastric randomized trials (E41 inputs). B: the three single-arm cohorts (expectations G1-G4).
C, D: earlier results, reported only. Bootstrap seeds [4044, i], i = 0, 1, 2 (PRJEB25780, PRJEB60680, PERFECT).
"""
import hashlib, sys
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h44 = hashlib.sha256(open(f"{A}/E44_prereg.md", "rb").read()).hexdigest()
assert _h44 == open(f"{A}/E44_prereg.sha256").read().split()[0], "E44 protocol hash mismatch"
_hw = hashlib.sha256(open(f"{A}/results/E37_frozen_weights.tsv", "rb").read()).hexdigest()
assert _hw == open(f"{A}/E37_frozen.sha256").read().split()[0], "frozen weights hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E44_step0_reproduction").replace("E34_coverage_check.csv", "E44_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
from scipy.stats import spearmanr, chi2


class Tee:
    def __init__(self, *f): self.f = f
    def write(self, s): [x.write(s) for x in self.f]
    def flush(self): [x.flush() for x in self.f]


sys.stdout = Tee(sys.__stdout__, open(f"{O}/E44.log", "w"))
print("E44 protocol sha256 verified:", _h44[:12], "| frozen weights sha256 verified:", _hw[:12])
ZQ = norm.ppf(0.975); MDD_K = norm.ppf(0.975) + norm.ppf(0.80)


def fe(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float); w = 1 / se ** 2
    m = np.sum(w * y) / w.sum(); s = np.sqrt(1 / w.sum())
    return dict(est=m, se=s, lo=m - ZQ * s, hi=m + ZQ * s, p=2 * norm.sf(abs(m / s)), k=len(y),
                Q=float(np.sum(w * (y - m) ** 2)))


def dl(y, se):
    y, v = np.asarray(y, float), np.asarray(se, float) ** 2; w0 = 1 / v; k = len(y)
    m0 = np.sum(w0 * y) / w0.sum(); Q = np.sum(w0 * (y - m0) ** 2)
    t2 = max(0.0, (Q - (k - 1)) / (w0.sum() - np.sum(w0 ** 2) / w0.sum())); w = 1 / (v + t2)
    m = np.sum(w * y) / w.sum(); s = np.sqrt(1 / w.sum())
    return dict(est=m, se=s, lo=m - ZQ * s, hi=m + ZQ * s, p=2 * norm.sf(abs(m / s)), k=k, tau2=t2, Q=Q)


# ================================================================== A. published randomized gastric trials
print("\n" + "=" * 100 + "\nA. Published randomized gastric trials (E41 inputs, descriptive)")
E41 = pd.read_csv(f"{O}/E41_inputs.csv")
E41 = E41[(E41.eligible == "yes") & (E41.in_ipd_six == "no")]; assert len(E41) == 18
GA = E41[E41.tumor_group == "gastric"]; OT = E41[E41.tumor_group != "gastric"]
assert sorted(GA.trial) == ["CheckMate 649", "KEYNOTE-061"] and len(OT) == 16
RA = []
for _, r in GA.iterrows():
    RA.append(dict(row=r.trial, k=1, rhr=np.exp(-r.y), rhr_lo=np.exp(-(r.y + ZQ * r.se)), rhr_hi=np.exp(-(r.y - ZQ * r.se)),
                   per_sd=r.y_sd, per_sd_lo=r.y_sd - ZQ * r.se_sd, per_sd_hi=r.y_sd + ZQ * r.se_sd,
                   p=2 * norm.sf(abs(r.y_sd / r.se_sd)), mdd80=MDD_K * r.se_sd, n_biomarker=r.n_biomarker))
for lab, f in [("gastric two, fixed effect", fe), ("gastric two, DerSimonian-Laird", dl)]:
    q = f(GA.y_sd, GA.se_sd); qr = f(GA.y, GA.se)
    RA.append(dict(row=lab, k=2, rhr=np.exp(-qr["est"]), rhr_lo=np.exp(-qr["hi"]), rhr_hi=np.exp(-qr["lo"]), per_sd=q["est"],
                   per_sd_lo=q["lo"], per_sd_hi=q["hi"], p=q["p"], mdd80=MDD_K * q["se"], Q=q["Q"],
                   Q_p=float(chi2.sf(q["Q"], 1))))
qo = dl(OT.y_sd, OT.se_sd); qg = dl(GA.y_sd, GA.se_sd)
RA.append(dict(row="other 16 trials, DerSimonian-Laird", k=16, per_sd=qo["est"], per_sd_lo=qo["lo"], per_sd_hi=qo["hi"],
               p=qo["p"], mdd80=MDD_K * qo["se"]))
dz = (qg["est"] - qo["est"]) / np.sqrt(qg["se"] ** 2 + qo["se"] ** 2)
RA.append(dict(row="gastric minus other 16 (DerSimonian-Laird each)", k=18, per_sd=qg["est"] - qo["est"],
               per_sd_lo=qg["est"] - qo["est"] - ZQ * np.sqrt(qg["se"] ** 2 + qo["se"] ** 2),
               per_sd_hi=qg["est"] - qo["est"] + ZQ * np.sqrt(qg["se"] ** 2 + qo["se"] ** 2), p=2 * norm.sf(abs(dz))))
RA = pd.DataFrame(RA); RA.to_csv(f"{O}/E44_A_gastric_trials.csv", index=False)
print(RA.round(4).to_string(index=False))
print(f"reference: six-comparison CD8A interaction +0.226 per SD; margin 0.10 per SD")

# ================================================================== B. single-arm gastroesophageal checkpoint cohorts
print("\n" + "=" * 100 + "\nB. Single-arm gastroesophageal checkpoint cohorts")
COH = ["PRJEB25780", "PRJEB60680", "PERFECT"]; GASTRIC = ["PRJEB25780", "PRJEB60680"]
LB = pd.read_csv(f"{O}/E44_labels.csv")
FW = pd.read_csv(f"{O}/E37_frozen_weights.tsv", sep="\t"); WG = dict(zip(FW.gene, FW.weight))
RAL = {"WARS1": "WARS"}                                          # library alias (E32 ALIAS_); matrices keep symbols otherwise
OBG, YG, MEG, NG = {}, {}, {}, {}
for c in COH:
    Z = pd.read_csv(f"{O}/E44_Z_{c}.csv.gz", index_col=0)
    Z = Z.rename(index={a: b for a, b in RAL.items() if a in Z.index and b not in Z.index})
    lb = LB[LB.cohort == c]; assert list(lb["sample"]) == list(Z.columns)
    YG[c] = lb.response.to_numpy(float); MEG[c] = lb.msi_ebv.to_numpy(float)
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
    if len(pres) >= 3 and len(pres) >= 0.5 * len(WG):
        OB[("sig", "barriernet 0.3")] = wscore(Z, pres, np.array([WG[g] for g in pres]))
    for g in T_PANEL + I_PANEL:
        if g in nc:
            OB[("gene", g)] = Z.loc[g].to_numpy(float)
    OBG[c] = OB; NG[c] = dict(bn03_genes=len(pres))
    print(f"{c}: n {len(YG[c])}, responders {int(YG[c].sum())}, MSI-H/EBV+ {int(np.nansum(MEG[c]))}; genes {Z.shape[0]}; "
          f"scored {sum(1 for kk in OB if kk[0] == 'sig' and kk[1] in LIB)}/52 signatures, extras "
          f"{[kk[1] for kk in OB if kk[0] == 'sig' and kk[1] not in LIB]}, panel genes {sum(1 for kk in OB if kk[0] == 'gene')}/27, "
          f"barriernet 0.3 genes {len(pres)}/{len(WG)}")
assert ("gene", "CD8A") in OBG["PRJEB25780"] and all(("gene", "CD8A") in OBG[c] for c in COH)


def logor(y, sc, x=None):
    cols = [np.ones(len(y))] + ([x] if x is not None else []) + [std(sc)]
    try:
        return lfit(np.column_stack(cols), y, [len(cols) - 1], True)[0]
    except Exception:
        return (np.nan, np.nan)


PC = []
for c in COH:
    for kk, sc in OBG[c].items():
        b, se = logor(YG[c], sc)
        row = dict(cohort=c, kind=kk[0], name=kk[1], n=len(YG[c]), responders=int(YG[c].sum()), logOR=b, se=se)
        if c in GASTRIC:
            ba, sa = logor(YG[c], sc, MEG[c]); row.update(logOR_msiebv=ba, se_msiebv=sa)
        PC.append(row)
PC = pd.DataFrame(PC); PC["p"] = 2 * norm.sf(np.abs(PC.logOR / PC.se)); PC.to_csv(f"{O}/E44_B_per_cohort.csv", index=False)

# ------------------------------------------------------------------ paired bootstrap against CD8A
BI = {}
for i, c in enumerate(COH):
    rng_ = np.random.default_rng([4044, i]); y = YG[c]; i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
    BI[c] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(1000)]


def boot(task):
    c, kk = task; sc = OBG[c][kk]; y = YG[c]
    return np.array([logor(y[idx], sc[idx])[0] for idx in BI[c]])


btasks = [(c, kk) for c in COH for kk in OBG[c]]
with mp.get_context("fork").Pool(NPROC) as pool:
    bres = pool.map(boot, btasks, chunksize=2)
BT = dict(zip(btasks, bres))
PD_ = []
for c in COH:
    b0 = BT[(c, ("gene", "CD8A"))]; y0 = PC[(PC.cohort == c) & (PC.name == "CD8A") & (PC.kind == "gene")].logOR.iloc[0]
    for kk in OBG[c]:
        if kk == ("gene", "CD8A"):
            continue
        db = BT[(c, kk)] - b0; ok = np.isfinite(db)
        y1 = PC[(PC.cohort == c) & (PC.name == kk[1]) & (PC.kind == kk[0])].logOR.iloc[0]
        PD_.append(dict(cohort=c, kind=kk[0], name=kk[1], diff=y1 - y0, se=np.std(db[ok], ddof=1), boot_ok=int(ok.sum()),
                        boot_failed=int((~ok).sum())))
PD_ = pd.DataFrame(PD_)
print(f"\nbootstrap: {len(btasks)} score x cohort series x 1,000; failed resamples (either fit): total "
      f"{int(PD_.boot_failed.sum())}, max per series {int(PD_.boot_failed.max())}")

# ------------------------------------------------------------------ pooling
POOL, PAIR, ADJ = [], [], []
for setname, cs in [("three gastroesophageal (primary)", COH), ("two gastric", GASTRIC)]:
    for (kind, name), g in PC.groupby(["kind", "name"], sort=False):
        g = g[g.cohort.isin(cs)]
        if len(g) < len(cs) or not np.isfinite(g.logOR).all() or not np.isfinite(g.se).all():
            continue
        q = fe(g.logOR, g.se); r = pool_rule(g.logOR.values, g.se.values)
        POOL.append(dict(set=setname, kind=kind, name=name, k=len(g), est=q["est"], lo=q["lo"], hi=q["hi"], p=q["p"], Q=q["Q"],
                         k_pos=int((g.logOR > 0).sum()), cons_est=r["est"], cons_lo=r["lo"], cons_hi=r["hi"], cons_p=r["p"]))
        if cs == GASTRIC:
            qa = fe(g.logOR_msiebv, g.se_msiebv)
            ADJ.append(dict(set=setname, kind=kind, name=name, est=q["est"], p=q["p"], est_msiebv=qa["est"], lo_msiebv=qa["lo"],
                            hi_msiebv=qa["hi"], p_msiebv=qa["p"]))
        if (kind, name) != ("gene", "CD8A"):
            d = PD_[(PD_.kind == kind) & (PD_.name == name) & PD_.cohort.isin(cs)]
            if len(d) == len(cs):
                qd = fe(d["diff"], d.se)
                PAIR.append(dict(set=setname, kind=kind, name=name, diff=qd["est"], lo=qd["lo"], hi=qd["hi"], p=qd["p"]))
POOL, PAIR, ADJ = pd.DataFrame(POOL), pd.DataFrame(PAIR), pd.DataFrame(ADJ)
POOL.to_csv(f"{O}/E44_B_pooled.csv", index=False); PAIR.to_csv(f"{O}/E44_B_paired.csv", index=False)
ADJ.to_csv(f"{O}/E44_B_adjusted.csv", index=False)

E32P = pd.read_csv(f"{O}/E32_primary.csv"); E32P = E32P[E32P.set == "six (primary)"].set_index("signature")
SUM = {}
for setname in ["three gastroesophageal (primary)", "two gastric"]:
    P = POOL[POOL.set == setname]; PR_ = PAIR[PAIR.set == setname]
    L = P[(P.kind == "sig") & P.name.isin(list(LIB))].set_index("name")
    pr = PR_[(PR_.kind == "sig") & PR_.name.isin(list(LIB))].set_index("name")
    cd = P[(P.kind == "gene") & (P.name == "CD8A")].iloc[0]
    bn = P[(P.kind == "sig") & (P.name == "barriernet 0.3")]
    bnd = PR_[(PR_.kind == "sig") & (PR_.name == "barriernet 0.3")]
    rho, rp = spearmanr(L.est, E32P.loc[L.index, "est"])
    sig = L[L.p < 0.05]
    same = sig[(E32P.loc[sig.index, "p"] < 0.05) & (np.sign(E32P.loc[sig.index, "est"]) == np.sign(sig.est))]
    allr = pd.concat([L.est, P[(P.kind == "gene")].set_index("name").est])
    SUM[setname] = dict(n_lib=len(L), cd8a=cd, bn=bn, bnd=bnd, rho=rho, rho_p=rp, above=int(((pr["diff"] > 0) & (pr.p < 0.05)).sum()),
                        below=int(((pr["diff"] < 0) & (pr.p < 0.05)).sum()), n_sig=len(sig), n_same=len(same),
                        rank=int((allr > cd.est).sum()) + 1, n_rank=len(allr), above_names=list(pr[(pr["diff"] > 0) & (pr.p < 0.05)].index))
    s_ = SUM[setname]
    print(f"\n--- {setname}: {s_['n_lib']} signatures pooled")
    print(f"CD8A log-OR per SD {cd.est:+.3f} ({cd.lo:+.3f} to {cd.hi:+.3f}), P = {cd.p:.2g}; conservative rule {cd.cons_est:+.3f} "
          f"({cd.cons_lo:+.3f} to {cd.cons_hi:+.3f}), P = {cd.cons_p:.2g}; rank {s_['rank']} of {s_['n_rank']} (52 signatures + panel genes)")
    if len(bn):
        b_ = bn.iloc[0]; d_ = bnd.iloc[0]
        print(f"barriernet 0.3 {b_.est:+.3f} ({b_.lo:+.3f} to {b_.hi:+.3f}), P = {b_.p:.2g}; minus CD8A {d_['diff']:+.3f} "
              f"({d_.lo:+.3f} to {d_.hi:+.3f}), P = {d_.p:.2g}")
    print(f"signatures significantly above CD8A: {s_['above']} {s_['above_names']}; below: {s_['below']}")
    print(f"signatures with P < 0.05: {s_['n_sig']} of {s_['n_lib']}; of these with a same-direction randomized interaction at "
          f"P < 0.05: {s_['n_same']}")
    print(f"Spearman with the six-comparison interaction: {rho:.3f} (P = {rp:.2g})")
    show = P[P.name.isin(["CD8A", "barriernet 0.3", "Immune", "NetImmune", "ESTIMATE", "TMEscore_z", "CXCL9", "PRF1", "GZMB",
                          "IFNG", "CD274"])]
    print(show[["kind", "name", "est", "lo", "hi", "p", "k_pos", "cons_p"]].round(3).to_string(index=False))
print("\nper cohort, CD8A / barriernet 0.3 / Immune:")
print(PC[PC.name.isin(["CD8A", "barriernet 0.3", "Immune"])][["cohort", "name", "n", "responders", "logOR", "se", "p"]].round(3)
      .to_string(index=False))
print("\nMSI-H/EBV-adjusted (two gastric):")
print(ADJ[ADJ.name.isin(["CD8A", "barriernet 0.3", "Immune", "TMEscore_z"])].round(3).to_string(index=False))
al = ADJ[ADJ.kind.eq("sig") & ADJ.name.isin(list(LIB))]
print(f"library signatures with P < 0.05 in the two gastric cohorts: unadjusted {int((al.p < 0.05).sum())}, "
      f"MSI-H/EBV-adjusted {int((al.p_msiebv < 0.05).sum())} (of {len(al)})")

# ================================================================== C, D: earlier results, reported only
print("\n" + "=" * 100 + "\nC. Prognosis without checkpoint inhibitors, 10 gastric cohorts (E20, primary duplicate rule)")
E20 = pd.read_csv(f"{O}/E30_dedup_runs/primary/E20_pooled.csv")
C_ = E20[(E20.set == "all 10 gastric") & (E20.model == "unadjusted") & (E20.duplicate_rule == "pre-registered exclusions (primary)")
         & E20.score.isin(["CD8A", "Immune", "T_CELL", "IFN", "ESTIMATE"])]
assert len(C_) == 5
C_[["set", "score", "k", "HR", "lo", "hi", "p", "I2"]].to_csv(f"{O}/E44_C_prognosis.csv", index=False)
print(C_[["score", "k", "HR", "lo", "hi", "p", "I2"]].round(3).to_string(index=False))
print("\nD. TMEscore (developed in gastric cancer) in the six randomized comparisons (E32)")
D_ = E32P.loc[["TMEscore_z"], ["est", "lo", "hi", "p", "k", "emp_p_declared", "diff_vs_CD8A", "cd8a_lo", "cd8a_hi", "cd8a_p"]]
D_.to_csv(f"{O}/E44_D_tmescore.csv"); print(D_.round(4).to_string())

# ================================================================== expectations
s_ = SUM["three gastroesophageal (primary)"]; cd = s_["cd8a"]; b_ = s_["bn"].iloc[0]; d_ = s_["bnd"].iloc[0]
EXP = pd.DataFrame([
    dict(q="G1 (primary)", expectation="CD8A pooled log-odds > 0 with P < 0.05 (fixed effect, three cohorts)",
         result=f"{cd.est:+.3f} ({cd.lo:+.3f} to {cd.hi:+.3f}), P = {cd.p:.2g}", held="yes" if cd.est > 0 and cd.p < 0.05 else "no"),
    dict(q="G2", expectation="no library signature significantly above CD8A (pooled paired difference > 0, P < 0.05)",
         result=f"{s_['above']} above ({', '.join(s_['above_names']) or 'none'}); {s_['below']} below",
         held="yes" if s_["above"] == 0 else "no"),
    dict(q="G3", expectation="Spearman correlation with the six-comparison randomized interaction > 0 with P < 0.05",
         result=f"{s_['rho']:.3f} (P = {s_['rho_p']:.2g}; {s_['n_lib']} signatures)", held="yes" if s_["rho"] > 0 and s_["rho_p"] < 0.05 else "no"),
    dict(q="G4", expectation="barriernet 0.3 pooled log-odds > 0 with P < 0.05, and its paired difference from CD8A P >= 0.05",
         result=f"{b_.est:+.3f} ({b_.lo:+.3f} to {b_.hi:+.3f}), P = {b_.p:.2g}; minus CD8A {d_['diff']:+.3f} ({d_.lo:+.3f} to {d_.hi:+.3f}), "
                f"P = {d_.p:.2g}", held="yes" if (b_.est > 0 and b_.p < 0.05 and d_.p >= 0.05) else "no")])
EXP.to_csv(f"{O}/E44_expectations.csv", index=False)
print("\n" + EXP.to_string(index=False))
