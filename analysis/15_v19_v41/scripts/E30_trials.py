"""E30 - randomized comparisons: Q3 (Angio adjustment in the renal trials), inputs for Q2 / Q4 (pooled in
E30_trials_meta.R), Q5 (unweighted, sqrt(n)-weighted and calibrated Stouffer z) and the Q6 table (CheckMate 025
overall survival, from existing results). Protocol: E30_prereg.md (2026-10-21; sha256 in E30_prereg.sha256),
written before any E30 model involving survival, response or arm was fitted.

Data, scoring, models: E28_base.py (the E25 code; E21 / E10 / E5 data loading) is executed unchanged except that
its step-0 table is written to results/E30_trials_E28base_reproduction.csv instead of E28_reproduction_check.csv,
so no published file is touched. Its own step 0 (interaction and within-arm effects < 1e-6) runs on load.

Step 0 (gate, extended here): every per-comparison interaction (y, se) of the 16 published scores in the six
primary comparisons against E21_meta_inputs.csv, and IMmotion150 atezolizumab monotherapy against E5_trials.csv,
|difference| < 1e-6; the pooled values are checked in E30_trials_meta.R. The script stops if any check fails.

Angio (McDermott 2018): z-mean of VEGFA, KDR, ESM1, PECAM1, ANGPTL4, CD34 present in the trial's z matrix
(>= 3 genes, the module rule), as in the read-only feasibility check (rev/agentB/angio_feasibility.py);
standardized within trial over the analysed patients, like every score.
Q3 model (renal trials, published primary model = Cox PFS, Efron ties):
    PFS ~ arm + s + Angio + s:arm + Angio:arm         benefit-oriented y = -coef(s:arm)
Usage: python3 E30_trials.py
"""
import hashlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h = hashlib.sha256(open(f"{A}/E30_prereg.md", "rb").read()).hexdigest()
assert _h == open(f"{A}/E30_prereg.sha256").read().split()[0], "E30 protocol hash mismatch"
print(f"E30_prereg.md sha256 {_h} (matches E30_prereg.sha256)")
_src = open(f"{A}/E28_base.py").read()
print(f"E28_base.py sha256 {hashlib.sha256(_src.encode()).hexdigest()}")
_OLD = 'REP.to_csv(f"{O}/E28_reproduction_check.csv", index=False)'
assert _src.count(_OLD) == 1 and _src.count(".to_csv(") == 1
_src = _src.replace(_OLD, 'REP.to_csv(f"{O}/E30_trials_E28base_reproduction.csv", index=False)')
exec(compile(_src, f"{A}/E28_base.py", "exec"))          # defines TR, scores, m_int, cox, std, sub, ... ; step 0

from scipy.stats import spearmanr

ORDER = ["CAF_ECM", "TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1", "IFN", "APM", "TLS_B", "T_CELL",
         "Barrier", "Immune", "NetImmune", "ESTIMATE", "TMEscore_PC", "TMEscore_z", "CD8A"]
Q3S = ["NetImmune", "Immune", "CAF_ECM", "Barrier", "ESTIMATE", "TMEscore_PC", "TMEscore_z", "CD8A"]
ANGIO = ["VEGFA", "KDR", "ESM1", "PECAM1", "ANGPTL4", "CD34"]
FIRST3 = ["IMM150", "ISPY2_pembro", "ISPY2_durva"]; ADDED3 = ["CM025", "JR101", "NeoTRIP"]
LABEL["IMM150_mono"] = "IMmotion150, atezolizumab monotherapy vs sunitinib"
V21 = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/validation"
FEAS = ("/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-com-gastric-2/"
        "727be15e-8d67-407e-a66a-a54d26cc1ebb/scratchpad/rev/agentB/angio_scores_six_trials.csv.gz")
WB = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/final_v31/BarrierNet_Supplementary_Tables.xlsx"
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
cancer = lambda k: "breast" if k in BREAST else "renal"

# ------------------------------------------------------------------ IMmotion150 arm split (E5_trials.py code path)
print("=" * 90, "\nIMmotion150 arm split (E5_trials.py: atezolizumab mono vs sunitinib)")
C = pd.read_csv(f"{P5}/v4_immotion150_clinical_api.csv", index_col=0)
M = pd.read_csv(f"{P5}/ici/rcc_iatlas_immotion150_2018.csv", index_col=0)
keep = [s for s in M.columns if s in C.index]; C, M = C.loc[keep], M[keep]
Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
S = scores(Z)
t = pd.to_numeric(C["PFS_MONTHS"], errors="coerce")
e = pd.to_numeric(C["PFS_STATUS"].astype(str).str.extract(r"^\s*([01])\s*:")[0], errors="coerce").fillna(0).astype(int)
atezo = C["ICI_RX"].astype(str) == "Atezolizumab"
suni = C["NON_ICI_RX"].astype(str).str.contains("Sunitinib", na=False)
bev = C["NON_ICI_RX"].astype(str).str.contains("Bevacizumab", na=False)
for key, trt in [("IMM150_check", atezo), ("IMM150_mono", atezo & ~bev)]:
    sel = ((trt | suni) & t.notna() & (t > 0)).values
    TR[key] = dict(kind="cox", Z=Z.loc[:, sel], S=S.loc[sel], arm=trt.values[sel].astype(int),
                   t=t.values[sel].astype(float), e=e.values[sel].astype(int), firth=True)
a_, b_ = TR.pop("IMM150_check"), TR["IMM150"]
assert (a_["arm"] == b_["arm"]).all() and np.allclose(a_["t"], b_["t"]) and (a_["e"] == b_["e"]).all()
assert np.nanmax(np.abs(a_["S"][ORDER].values - b_["S"][ORDER].values)) == 0
d = TR["IMM150_mono"]
print(f"  atezolizumab monotherapy {int(d['arm'].sum())}  sunitinib {int((1 - d['arm']).sum())}  events {int(d['e'].sum())}")


# ------------------------------------------------------------------ step 0 (gate), per comparison
def fit_int(d, sc):
    m = np.isfinite(d["S"][sc].values); dd = sub(d, m)
    y, se = m_int(dd, d["S"][sc].values[m])
    return y, se, int(m.sum()), dd


MI = pd.read_csv(f"{O}/E21_meta_inputs.csv"); MI6 = MI[MI.set == "primary: 6 comparisons"]
E5 = pd.read_csv(f"{O}/E5_trials.csv")
gate, inputs = [], []
for k in ALL6 + ["IMM150_mono"]:
    d = TR[k]
    for sc in ORDER:
        y, se, n, dd = fit_int(d, sc)
        if k == "IMM150_mono":
            r = E5[(E5.trial == "IMmotion150: atezolizumab mono vs sunitinib") & (E5.score == sc)].iloc[0]
            yp, sp, src = -np.log(r.int_HR), (np.log(r.int_hi) - np.log(r.int_lo)) / (2 * norm.ppf(0.975)), "E5_trials.csv (mono)"
        else:
            r = MI6[(MI6.trial == LABEL[k]) & (MI6.score == sc)].iloc[0]
            yp, sp, src = r.y, r.se, "E21_meta_inputs.csv (primary: 6 comparisons)"
        gate.append(dict(check="per-comparison interaction (y, se)", source=src, trial=LABEL[k], score=sc, y_e30=y,
                         y_pub=yp, se_e30=se, se_pub=sp, absdiff=max(abs(y - yp), abs(se - sp))))
        ev = int(dd["e"].sum()) if d["kind"] == "cox" else int(dd["y"].sum())
        inputs.append(dict(analysis="unadjusted", key=k, trial=LABEL[k], cancer=cancer(k), score=sc, n=n, events=ev,
                           y=y, se=se, p=2 * norm.sf(abs(y / se))))
G = pd.DataFrame(gate)
print(f"\nstep 0 (E30): {len(G)} per-comparison interactions (y and se) of 16 scores in the six primary comparisons "
      f"and IMmotion150 monotherapy; largest |difference| from the published values {G.absdiff.max():.2e}")
G.to_csv(f"{O}/E30_trials_gate.csv", index=False)
assert G.absdiff.max() < 1e-6, G.sort_values("absdiff").tail()
print("step 0 per-comparison gate: PASS")

# ------------------------------------------------------------------ Angio score and its checks
print("=" * 90, "\nAngio (VEGFA, KDR, ESM1, PECAM1, ANGPTL4, CD34)")
feas = pd.read_csv(FEAS, index_col=0); feas.index = feas.index.astype(str)
ang_rows = []
for k in ALL6 + ["IMM150_mono"]:
    d = TR[k]; g = [x for x in ANGIO if x in d["Z"].index]
    assert len(g) >= 3
    d["Angio"] = d["Z"].loc[g].mean(0).values; d["angio_genes"] = " ".join(g)
    row = dict(trial=LABEL[k], genes=" ".join(g), n_genes=len(g))
    if k != "IMM150_mono":
        f = feas[feas.trial == k]; f = f.loc[d["S"].index.astype(str)]
        row["maxdiff_vs_feasibility"] = float(np.max(np.abs(f.Angio.values - d["Angio"])))
    for sc in ["NetImmune", "Immune", "CAF_ECM", "Barrier"]:
        row[f"rho_{sc}"] = spearmanr(d["Angio"], d["S"][sc].values)[0]
    ang_rows.append(row)
c = pd.read_csv(f"{V21}/Braun2020_CheckMate/S1_clinical.tsv", sep="\t"); c = c[c.Cohort == "CM-025"].set_index("RNA_ID")
ref = c.loc[TR["CM025"]["S"].index, "Angio"].values; ok = np.isfinite(ref)
rho_braun = spearmanr(TR["CM025"]["Angio"][ok], ref[ok])[0]
j = pd.read_excel(f"{V21}/JAVELIN_Renal101/41591_2020_1044_MOESM3_ESM.xlsx", sheet_name="S16_Pathway_scores", header=1)
j = j.set_index(j.row_id.astype(str))
ids = [s for s in TR["JR101"]["S"].index.astype(str)]
jm = j.reindex(ids)["MCDERMOTT_ANGIO"].values; ok2 = np.isfinite(jm)
rho_jav = spearmanr(TR["JR101"]["Angio"][ok2], jm[ok2])[0]
AN = pd.DataFrame(ang_rows)
print(AN.round(3).to_string(index=False))
print(f"Spearman rho with the published Angio: CheckMate 025 (Braun S1 'Angio', n = {ok.sum()}) {rho_braun:.3f}; "
      f"JAVELIN Renal 101 (Motzer S16 MCDERMOTT_ANGIO, n = {ok2.sum()}) {rho_jav:.3f}")
assert np.nanmax(AN.maxdiff_vs_feasibility.values) < 1e-12 and ok2.sum() == len(ids)


# ------------------------------------------------------------------ Q3: Angio-adjusted interactions (renal)
def m_adj(d, s, g):
    """Published primary model + Angio + Angio x arm: (score x arm y, se), (Angio x arm y, se), benefit-oriented."""
    s, g = std(s), std(g); a = d["arm"]
    (bs, ss), (bg, sg) = cox({"s": s, "arm": a, "g": g, "sx": s * a, "gx": g * a}, d["t"], d["e"], ["sx", "gx"])
    return -bs, ss, -bg, sg


print("=" * 90, "\nQ3 per renal trial: published primary model + Angio + Angio x arm")
q3 = []
for k in RENAL + ["IMM150_mono"]:
    d = TR[k]; g = d["Angio"]
    ya, sa, n, _ = fit_int(dict(d, S=pd.DataFrame({"Angio": g}, index=d["S"].index)), "Angio")
    inputs.append(dict(analysis="Angio alone", key=k, trial=LABEL[k], cancer="renal", score="Angio", n=n,
                       events=int(d["e"].sum()), y=ya, se=sa, p=2 * norm.sf(abs(ya / sa))))
    for sc in Q3S:
        s = d["S"][sc].values; assert np.isfinite(s).all() and np.isfinite(g).all()
        yu, su = m_int(d, s)
        y1, s1, y2, s2 = m_adj(d, s, g)
        for an, yy, ss in [("Angio-adjusted", y1, s1), ("Angio x arm (in score model)", y2, s2)]:
            inputs.append(dict(analysis=an, key=k, trial=LABEL[k], cancer="renal", score=sc, n=len(s),
                               events=int(d["e"].sum()), y=yy, se=ss, p=2 * norm.sf(abs(yy / ss))))
        q3.append(dict(trial=LABEL[k], key=k, role="descriptive (no decision)" if k == "IMM150_mono" else "primary",
                       score=sc, n=len(s), events=int(d["e"].sum()), angio_genes=d["angio_genes"],
                       rho_angio=spearmanr(g, s)[0], y_unadj=yu, se_unadj=su, p_unadj=2 * norm.sf(abs(yu / su)),
                       y_adj=y1, se_adj=s1, p_adj=2 * norm.sf(abs(y1 / s1)), diff_adj_minus_unadj=y1 - yu,
                       within_0_05=bool(abs(y1 - yu) <= 0.05), angio_x_arm_y=y2, angio_x_arm_se=s2,
                       angio_x_arm_p=2 * norm.sf(abs(y2 / s2)), angio_alone_y=ya, angio_alone_se=sa,
                       angio_alone_p=2 * norm.sf(abs(ya / sa))))
Q3 = pd.DataFrame(q3)
print(Q3[["trial", "score", "n", "rho_angio", "y_unadj", "y_adj", "diff_adj_minus_unadj", "p_unadj", "p_adj",
          "angio_x_arm_y", "angio_x_arm_p"]].round(4).to_string(index=False))
Q3.to_csv(f"{O}/E30_trials_angio_per_trial.csv", index=False)
IN = pd.DataFrame(inputs); IN.to_csv(f"{O}/E30_trials_inputs.csv", index=False)

# ------------------------------------------------------------------ Q5: Stouffer, unweighted / sqrt(n)-weighted / calibrated
print("=" * 90, "\nQ5 Stouffer z: unweighted (Table 4), sqrt(n)-weighted, calibrated on the E28 parametric bootstrap")
SETS5 = {"all six": ALL6, "renal": RENAL, "breast": BREAST, "first three": FIRST3, "three added": ADDED3}
E21SET = {"all six": "primary: 6 comparisons", "renal": "renal only (IMmotion150, CheckMate 025, JAVELIN)",
          "breast": "breast only (I-SPY2 x 2, NeoTRIP)", "first three": "reproduction: 3 published primary",
          "three added": "3 new primary only"}
M21 = pd.read_csv(f"{O}/E21_meta.csv"); NB = pd.read_csv(f"{O}/E28_Q3_null_draws.csv")
CAL = pd.read_csv(f"{O}/E28_Q3_calibration_pooled.csv").set_index("score")
U = IN[IN.analysis == "unadjusted"].set_index(["key", "score"])
q5 = []
for sc in ORDER:
    W6 = None
    if sc in NB.score.unique():
        W6 = NB[NB.score == sc].pivot(index="b", columns="key", values=["y", "se"]).dropna()      # as E28_trials.py
    for sn, ks in SETS5.items():
        y = np.array([U.loc[(k, sc), "y"] for k in ks]); se = np.array([U.loc[(k, sc), "se"] for k in ks])
        n = np.array([U.loc[(k, sc), "n"] for k in ks], float); z = y / se; J = len(ks)
        zu = z.sum() / np.sqrt(J); zw = np.sum(np.sqrt(n) * z) / np.sqrt(n.sum())
        pub = M21[(M21.set == E21SET[sn]) & (M21.score == sc)].iloc[0]
        row = dict(score=sc, set=sn, k=J, n_total=int(n.sum()), n_per_comparison=" / ".join(str(int(v)) for v in n),
                   stouffer_z=zu, stouffer_p=2 * norm.sf(abs(zu)), weighted_z=zw, weighted_p=2 * norm.sf(abs(zw)),
                   published_stouffer_z=pub.stouffer_z, published_stouffer_p=pub.stouffer_p,
                   absdiff_vs_published=abs(zu - pub.stouffer_z))
        if W6 is not None:
            zb = (W6["y"][ks].values / W6["se"][ks].values)
            zbu = zb.sum(1) / np.sqrt(J); zbw = (np.sqrt(n) * zb).sum(1) / np.sqrt(n.sum())
            row.update(n_null=len(zb), calibrated_p=(1 + np.sum(np.abs(zbu) >= abs(zu))) / (1 + len(zb)),
                       calibrated_weighted_p=(1 + np.sum(np.abs(zbw) >= abs(zw))) / (1 + len(zb)),
                       pct_null_p05=100 * np.mean(np.abs(zbu) > norm.ppf(0.975)),
                       pct_null_weighted_p05=100 * np.mean(np.abs(zbw) > norm.ppf(0.975)))
            if sn == "all six":
                row["published_calibrated_p_E28"] = CAL.loc[sc, "stouffer_p_calibrated"]
        q5.append(row)
Q5 = pd.DataFrame(q5)
assert Q5.absdiff_vs_published.max() < 1e-6, Q5.sort_values("absdiff_vs_published").tail()
chk = Q5.dropna(subset=["published_calibrated_p_E28"])
assert np.allclose(chk.calibrated_p, chk.published_calibrated_p_E28, atol=1e-12, rtol=0), chk
print(f"unweighted Stouffer z reproduces E21_meta.csv (largest |difference| {Q5.absdiff_vs_published.max():.2e}); "
      f"calibrated P for the six comparisons reproduces E28_Q3_calibration_pooled.csv exactly ({len(chk)} scores)")
Q5.to_csv(f"{O}/E30_trials_stouffer.csv", index=False)
print(Q5[Q5.set == "all six"][["score", "k", "n_total", "stouffer_z", "stouffer_p", "weighted_z", "weighted_p",
                               "calibrated_p", "calibrated_weighted_p"]].to_string(index=False))

# ------------------------------------------------------------------ Q6: CheckMate 025 overall survival (existing results)
print("=" * 90, "\nQ6 CheckMate 025 with overall survival (E21_interactions.csv, E21_meta.csv; workbook S16i / S16f)")
I21 = pd.read_csv(f"{O}/E21_interactions.csv")
os_ = I21[(I21.key == "CM025") & (I21.endpoint == "OS") & (I21.model == "Cox")].set_index("score")
pfs = I21[(I21.key == "CM025") & (I21.endpoint == "PFS") & (I21.model == "Cox")].set_index("score")
mo = M21[M21.set == "primary, CheckMate 025 OS for PFS"].set_index("score")
mp6 = M21[M21.set == "primary: 6 comparisons"].set_index("score")
s16i = pd.read_excel(WB, sheet_name="S16i_pooled_all_sets"); s16i = s16i[s16i.set == "primary, CheckMate 025 OS for PFS"].set_index("score")
s16f = pd.read_excel(WB, sheet_name="S16f_added_interactions")
s16f = s16f[(s16f.key == "CM025") & (s16f.endpoint == "OS")].set_index("score")
q6 = []
for sc in ORDER:
    o, p_, m, m6 = os_.loc[sc], pfs.loc[sc], mo.loc[sc], mp6.loc[sc]
    row = dict(score=sc, cm025_os_n=int(o.n_ici + o.n_ctrl), cm025_os_events=int(o.ev_ici + o.ev_ctrl),
               cm025_os_int_HR=o.int_est, cm025_os_lo=o.int_lo, cm025_os_hi=o.int_hi, cm025_os_p=o.int_p,
               cm025_os_y=o.y, cm025_os_se=o.int_se, cm025_pfs_int_HR=p_.int_est, cm025_pfs_y=p_.y, cm025_pfs_p=p_.int_p,
               six_os_est=m.est, six_os_lo=m.lo, six_os_hi=m.hi, six_os_p=m.p, six_os_interval=m.interval,
               six_os_k_benefit=int(m.k_benefit), six_os_stouffer_p=m.stouffer_p,
               six_pfs_est=m6.est, six_pfs_lo=m6.lo, six_pfs_hi=m6.hi, six_pfs_p=m6.p,
               sig_pfs=bool(m6.p < 0.05), sig_os=bool(m.p < 0.05),
               source="E21_interactions.csv; E21_meta.csv (set 'primary, CheckMate 025 OS for PFS')")
    if sc in s16i.index:
        row["workbook_maxdiff"] = max(np.max(np.abs(s16i.loc[sc, ["est", "lo", "hi", "p"]].values.astype(float)
                                                    - np.array([m.est, m.lo, m.hi, m.p]))),
                                      np.max(np.abs(s16f.loc[sc, ["y", "int_se", "int_p"]].values.astype(float)
                                                    - np.array([o.y, o.int_se, o.int_p]))))
        row["source"] += "; workbook S16i / S16f"
    else:
        row["source"] += "; not in workbook S16i (comparator)"
    q6.append(row)
Q6 = pd.DataFrame(q6); Q6.to_csv(f"{O}/E30_trials_cm025_os.csv", index=False)
print(Q6[["score", "cm025_os_int_HR", "cm025_os_p", "six_os_est", "six_os_lo", "six_os_hi", "six_os_p", "six_pfs_p",
          "workbook_maxdiff"]].round(4).to_string(index=False))
print(f"workbook S16i / S16f against E21 results files: largest |difference| {np.nanmax(Q6.workbook_maxdiff):.2e}")
print("\nE30_trials.py done")
