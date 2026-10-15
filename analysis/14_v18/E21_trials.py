"""E21 - score x arm interaction in further randomised checkpoint trials (pre-registered in
E17_E22_prereg.md, section E21 and the Deviations entry of 2026-10-05). Steps 2 and 3 of E21:
per-trial interactions with within-arm effects, and the structure-matched random-gene null.
Step 1 (Braun consistency check, no outcomes) is E21_braun_check.R; step 4 (pooling) is
E21_meta.R; the figure is E21_fig_trials.py.

Trials (data_download/validation/<folder>; see each README.txt):
  Primary, added to the 3 existing primary comparisons (IMmotion150, I-SPY2 x 2):
    CheckMate 025      nivolumab (1) vs everolimus (0), renal; Braun 2020 S1 + S4A (log2 already:
                       UQ-TPM, ComBat, +21; not re-logged). PFS primary; OS secondary (the trial's
                       own primary endpoint). PFS_CNSR / OS_CNSR = 1 means EVENT in this file.
    JAVELIN Renal 101  avelumab + axitinib (1) vs sunitinib (0), renal; S11 + S13 (log2 TPM floored
                       at 0.01; not re-logged). PFS only. PFS_P_CNSR follows ADaM: 0 = EVENT
                       (verified below from arm-level medians against Supplementary Table 2).
    NeoTRIP            carboplatin + nab-paclitaxel + atezolizumab (1) vs the same chemotherapy (0),
                       TNBC; GSE319641 baseline TPM (linear, ComBat; log2(max(x,0) + 1) once) with arm
                       and pCR from the authors' GitHub table. pCR; sequencing batch as a main-effect
                       covariate (primary), without batch (sensitivity).
  Sensitivity only (like PRINCE):
    DUTRENEO           durvalumab + tremelimumab (1) vs cisplatin-based chemotherapy (0), TIS-high
                       bladder; GSE328930 raw counts -> log2(CPM + 1). pCR (no surgery = no pCR).
                       Post-hoc further sensitivity: libraries >= 1M counts only.

Scoring exactly as in E5_trials.py / E10_ispy2_trials.py: duplicate symbols collapsed by mean;
genes constant across the trial dropped (no z; only DUTRENEO TAP2, all zero in 3' QuantSeq);
genes row-z within trial (across both arms; CheckMate 025 = its 250 randomised patients with RNA,
not CM-009/010); module = mean z of its genes present (>= 3); Barrier = mean(CAF_ECM, TAM_SUP);
Immune = mean(APM, IFN, TLS_B, T_CELL); NetImmune = Immune - Barrier; C1Q arm, SPP1 arm,
CXCL9:SPP1 = z(CXCL9) - z(SPP1). Each score standardised to SD 1 over the analysed patients.

Models (as E5 / E10):
  Cox    PFS ~ s + arm + s:arm (lifelines, Efron ties); interaction HR = exp(coef of s:arm).
  Logit  pCR ~ s + arm + s:arm [+ batch]; interaction OR; Firth (Jeffreys) where ML separates
         (non-convergence, |coef| > 10 or SE > 10), Firth estimates also written as a check.
  Within-arm HR / OR per SD from separate fits in each arm (NeoTRIP: batch-adjusted and not).
  Benefit-oriented log ratio y = -log(interaction HR) or +log(interaction OR): y > 0 means the
  score predicts more benefit from the checkpoint-containing arm.

Structure-matched random-gene null (as E10): 2,000 random gene sets from the bulk gene pool
(results/bulk_gene_pool.txt) measured, non-constant, in all three new primary trials; the same
genes in each trial. NetImmune: four "immune" modules minus two "barrier" modules, all distinct,
sized as the real modules in that trial (CheckMate 025 8/5/6/6 - 11/11, JAVELIN 8/5/6/5 - 11/11,
NeoTRIP 8/5/6/7 - 11/11; slots filled in the same order so a shorter module uses the first genes
of the slot); CAF_ECM: 11 genes. Statistic: each trial's primary-model y, and their inverse-
variance combination; random sets credited with their better direction,
P = (#|null| >= |obs| + 1) / (N + 1).

Pre-registered expectations (benefit-oriented): NetImmune, Immune, IFN > 0; CAF_ECM <= 0.
"Held" = point estimate on the expected side; "supported" = also P < 0.05.
"""
import warnings, functools; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy.stats import norm
from lifelines import CoxPHFitter, KaplanMeierFitter

V = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/validation"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
NRAND = 2000
rng = np.random.default_rng(2026)
MOD = dict(
    CAF_ECM=["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
    TAM_SUP=["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
    APM=["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    IFN=["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
    TLS_B=["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
    T_CELL=["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"],
    C1Q_arm=["C1QA", "C1QB", "C1QC"],
    SPP1_arm=["SPP1", "TREM2", "MARCO", "APOE", "LGALS3"])
ORDER = ["CAF_ECM", "TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1", "IFN", "APM", "TLS_B", "T_CELL",
         "Barrier", "Immune", "NetImmune"]
EXPECT = {"NetImmune": ">0", "Immune": ">0", "IFN": ">0", "CAF_ECM": "<=0"}


# ------------------------------------------------------------------ scoring (E5 / E10, unchanged)
def scores(Z):
    S = {}
    for k, G in MOD.items():
        g = [x for x in G if x in Z.index]
        S[k] = Z.loc[g].mean(0) if len(g) >= 3 else pd.Series(np.nan, index=Z.columns)
    S = pd.DataFrame(S)
    S["Barrier"] = S[["CAF_ECM", "TAM_SUP"]].mean(1)
    S["Immune"] = S[["APM", "IFN", "TLS_B", "T_CELL"]].mean(1)
    S["NetImmune"] = S.Immune - S.Barrier
    S["CXCL9_SPP1"] = Z.loc["CXCL9"] - Z.loc["SPP1"]
    return S


def rowz(X):
    """Row-z as E5/E10; genes constant across the trial carry no z and are dropped first (as in
    E5_icb_split.R), so they do not count towards a module's >= 3 genes (only DUTRENEO TAP2)."""
    X = X[X.std(1) > 0]
    return ((X.T - X.mean(1)) / X.std(1)).T


def std(s):
    s = np.asarray(s, float)
    return (s - s.mean()) / s.std(ddof=1)


# ------------------------------------------------------------------ logistic ML / Firth (E10, unchanged)
def ml(X, y, it=100, tol=1e-10):
    b = np.zeros(X.shape[1])
    for _ in range(it):
        eta = np.clip(X @ b, -30, 30); p = 1 / (1 + np.exp(-eta)); W = p * (1 - p)
        I = X.T @ (X * W[:, None])
        try:
            step = np.linalg.solve(I, X.T @ (y - p))
        except np.linalg.LinAlgError:
            return b, np.full((len(b), len(b)), np.nan), False
        b = b + step
        if np.max(np.abs(step)) < tol:
            p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
            return b, np.linalg.inv(X.T @ (X * W[:, None])), True
    return b, np.full((len(b), len(b)), np.nan), False


def firth(X, y, tol=1e-9, it=200):
    k = X.shape[1]; b = np.zeros(k)
    def pll(b):
        p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
        return np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)) + 0.5 * np.linalg.slogdet(X.T @ (X * W[:, None]))[1]
    for _ in range(it):
        p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
        I = X.T @ (X * W[:, None]); Ii = np.linalg.inv(I)
        h = np.einsum("ij,jk,ik->i", X * np.sqrt(W)[:, None], Ii, X * np.sqrt(W)[:, None])
        step = np.linalg.solve(I, X.T @ (y - p + h * (0.5 - p)))
        f0, t = pll(b), 1.0
        while t > 1e-6 and pll(b + t * step) < f0 - 1e-12:
            t /= 2
        b = b + t * step
        if np.max(np.abs(t * step)) < tol:
            break
    p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
    return b, np.linalg.inv(X.T @ (X * W[:, None]))


def lfit(X, y, j):
    """Coefficient j: (estimate, SE, method, Firth estimate, Firth SE); Firth if ML separates."""
    b, V_, ok = ml(X, y)
    sep = (not ok) or np.any(np.abs(b) > 10) or np.any(~np.isfinite(np.diag(V_))) or np.any(np.sqrt(np.diag(V_)) > 10)
    bf, Vf = firth(X, y)
    if sep:
        return bf[j], np.sqrt(Vf[j, j]), "Firth", bf[j], np.sqrt(Vf[j, j])
    return b[j], np.sqrt(V_[j, j]), "ML", bf[j], np.sqrt(Vf[j, j])


def ci(b, se):
    return dict(log=b, se=se, est=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se),
                p=2 * norm.sf(abs(b / se)))


# ------------------------------------------------------------------ Cox (E5 fit(), plus SE)
def cox_int(s, arm, t, e):
    df = pd.DataFrame({"T": t, "E": e, "s": s, "arm": arm}); df["sx"] = df.s * df.arm
    r = CoxPHFitter().fit(df, "T", "E").summary.loc["sx"]
    return r["coef"], r["se(coef)"]


def cox_within(s, t, e):
    r = CoxPHFitter().fit(pd.DataFrame({"T": t, "E": e, "s": s}), "T", "E").summary.loc["s"]
    return r["coef"], r["se(coef)"]


def design(s, arm, cov, interaction=True):
    one = np.ones_like(s)
    cols = [one, s, arm, s * arm] if interaction else [one, s]
    if cov is not None and cov.shape[1]:
        cols += [cov[:, j] for j in range(cov.shape[1])]
    return np.column_stack(cols)


# ------------------------------------------------------------------ data loading + scale checks
scale_rows, cover_rows = [], []


def scale_check(name, X, deposited, transform):
    v = X.values.ravel()
    scale_rows.append(dict(trial=name, deposited_scale=deposited, transform_applied=transform, genes=X.shape[0],
                           samples=X.shape[1], min=v.min(), median=np.median(v), max=v.max(), frac_negative=np.mean(v < 0)))


def coverage(name, X):
    for k, G in MOD.items():
        g = [x for x in G if x in X.index]
        zero = [x for x in g if X.loc[x].std() == 0]
        cover_rows.append(dict(trial=name, module=k, genes_listed=len(G), genes_present=len(g), genes_scored=len(g) - len(zero),
                               missing=";".join(x for x in G if x not in X.index), zero_variance=";".join(zero)))


def collapse(X):
    nd = int(X.index.duplicated().sum())
    return (X.groupby(level=0).mean() if nd else X), nd


def km_median(t, e):
    return KaplanMeierFitter().fit(t, e).median_survival_time_


TR = {}   # key -> dict(label, Z, S, arm, outcome..., kind)

print("=" * 100, "\nCheckMate 025 (Braun 2020)")
C = pd.read_csv(f"{V}/Braun2020_CheckMate/S1_clinical.tsv", sep="\t")
X = pd.read_csv(f"{V}/Braun2020_CheckMate/S4A_RNA_expression.tsv.gz", sep="\t", index_col=0)
X, nd = collapse(X)
c25 = C[C.Cohort == "CM-025"]
print("arm-level check of event coding (all CM-025 patients in S1, _CNSR = 1 as event; NEJM 2015: PFS 4.6 vs 4.4 mo, OS 25.0 vs 19.6 mo):")
for a in ["NIVOLUMAB", "EVEROLIMUS"]:
    d = c25[c25.Arm == a]
    print(f"   {a:10s} n={len(d)}  median PFS {km_median(d.PFS, d.PFS_CNSR):.1f}  median OS {km_median(d.OS, d.OS_CNSR):.1f}")
c = c25[c25.RNA_ID.isin(X.columns)].set_index("RNA_ID")
Xc = X[c.index]
scale_check("CheckMate 025", Xc, "log2 (UQ-TPM, ComBat, +21)", "none (already log2)")
print(f"duplicate symbols collapsed: {nd}; RNA patients: {c.Arm.value_counts().to_dict()}; scale min/median/max "
      f"{Xc.values.min():.2f}/{np.median(Xc.values):.2f}/{Xc.values.max():.2f}")
coverage("CheckMate 025", Xc)
Z = rowz(Xc); S = scores(Z)
TR["CM025"] = dict(label="CheckMate 025", Z=Z, S=S, clin=c, arm=(c.Arm == "NIVOLUMAB").astype(int).values)

print("=" * 100, "\nJAVELIN Renal 101")
C = pd.read_csv(f"{V}/JAVELIN_Renal101/JR101_S11_clinical.tsv", sep="\t")
X = pd.read_csv(f"{V}/JAVELIN_Renal101/JR101_S13_expression_log2TPM.tsv.gz", sep="\t", index_col=0)
X, nd = collapse(X)
Ca = C[C.TRT01P.notna()]
print("arm-level check of PFS_P_CNSR coding (PD-L1+ patients; Supplementary Table 2: 13.34 vs 8.21 mo):")
for code in [0, 1]:
    out = []
    for a in ["Avelumab+Axitinib", "Sunitinib"]:
        d = Ca[(Ca.TRT01P == a) & (Ca.PDL1FL == "Y")]
        out.append(f"{a} n={len(d)} median {km_median(d.PFS_P, (d.PFS_P_CNSR == code).astype(int)):.2f}")
    print(f"   event = (PFS_P_CNSR == {code}): " + " | ".join(out))
c = C.set_index("ID").loc[X.columns]
scale_check("JAVELIN Renal 101", X, "log2 TPM, values < 0.01 set to 0.01", "none (already log2)")
print(f"duplicate symbols collapsed: {nd}; RNA patients: {c.TRT01P.value_counts().to_dict()}")
coverage("JAVELIN Renal 101", X)
Z = rowz(X); S = scores(Z)
TR["JR101"] = dict(label="JAVELIN Renal 101", Z=Z, S=S, clin=c, arm=(c.TRT01P == "Avelumab+Axitinib").astype(int).values)

print("=" * 100, "\nNeoTRIP (GSE319641)")
P = pd.read_csv(f"{V}/NeoTRIP_GSE319641/github_Early_GEP_NeoTRIP/data/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt", sep="\t")
X = pd.read_csv(f"{V}/NeoTRIP_GSE319641/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz", sep="\t")
X = X.drop(columns=["HGNC_id", "entrez_gene_id"]).set_index("HGNC_approved_symbol")
c = P[P.Timepoint == "Baseline"].set_index("Sample_ID")
Xb = X[c.index]
v = Xb.values.ravel()
print(f"deposited baseline TPM: min {v.min():.3f}, median {np.median(v):.2f}, max {v.max():.3g}, frac < 0 {np.mean(v < 0):.3f} -> linear")
Xb = np.log2(Xb.clip(lower=0) + 1)
Xb, nd = collapse(Xb)
scale_check("NeoTRIP", Xb, "linear TPM after ComBat (7% slightly negative)", "log2(max(x, 0) + 1), once")
print(f"duplicate symbols collapsed: {nd}; baseline patients: {c.Arm.value_counts().to_dict()}")
print("arm x batch:\n", pd.crosstab(c.Arm, c.Batch_correlation).to_string())
coverage("NeoTRIP", Xb)
Z = rowz(Xb); S = scores(Z)
bat = pd.get_dummies(c.Batch_correlation)[["Batch 1", "Batch 2"]].astype(float).values   # reference: Batch 3
TR["NeoTRIP"] = dict(label="NeoTRIP", Z=Z, S=S, clin=c, arm=(c.Arm == "CT/A").astype(int).values,
                     y=c["pCR.num"].astype(float).values, batch=bat)

print("=" * 100, "\nDUTRENEO (GSE328930) - sensitivity only")
M = pd.read_csv(f"{V}/DUTRENEO_GSE328930/GSE328930_metadata.tsv.gz", sep="\t", index_col=0)
X = pd.read_csv(f"{V}/DUTRENEO_GSE328930/GSE328930_counts.tsv.gz", sep="\t", index_col=0)
c = M[M.Group.isin(["HOT-ICI", "HOT-NAC"])]
lib = X.sum(0)
print("response categories:", c.groupby(["Group", "Pathological_response"]).size().to_dict())
for key, lab, keep in [("DUTRENEO", "DUTRENEO", c.index),
                       ("DUTRENEO_1M", "DUTRENEO, libraries >= 1M", c.index[lib[c.index].values >= 1e6])]:
    cc = c.loc[keep]
    cpm = X[keep] / lib[keep] * 1e6
    Xd = np.log2(cpm + 1)
    Xd, nd = collapse(Xd)
    scale_check(lab, Xd, "raw UMI counts", "log2(CPM + 1), CPM on all genes' library size")
    print(f"{lab}: n={len(cc)} {cc.Group.value_counts().to_dict()}; library size min {lib[keep].min():.0f}, median {lib[keep].median():.0f}; dup collapsed {nd}")
    coverage(lab, Xd)
    Z = rowz(Xd); S = scores(Z)
    TR[key] = dict(label=lab, Z=Z, S=S, clin=cc, arm=(cc.Group == "HOT-ICI").astype(int).values,
                   y=(cc.Pathological_response == "Pathological complete response").astype(float).values, batch=None)

SCL = pd.DataFrame(scale_rows); COV = pd.DataFrame(cover_rows)
SCL.to_csv(f"{O}/E21_scale_checks.csv", index=False); COV.to_csv(f"{O}/E21_gene_coverage.csv", index=False)
print("\nmodule genes scored (present and non-constant):\n", COV.pivot(index="trial", columns="module", values="genes_scored").to_string())
print("missing / zero-variance:", COV[(COV.missing != "") | (COV.zero_variance != "")][["trial", "module", "missing", "zero_variance"]].to_string(index=False))

# per-patient scores (standardised within trial over all scored patients)
sc_rows = []
for key, d in TR.items():
    S = d["S"].copy(); S.insert(0, "arm", d["arm"]); S.insert(0, "trial", d["label"]); S.index.name = "sample"
    sc_rows.append(S.reset_index())
pd.concat(sc_rows).to_csv(f"{O}/E21_scores.csv", index=False)

# ------------------------------------------------------------------ step 2: interactions + within-arm
rows, arms = [], []


def survival(key, tcol, ecol, endpoint, role, event_is):
    d = TR[key]; c = d["clin"]
    t = pd.to_numeric(c[tcol], errors="coerce").values
    e = (pd.to_numeric(c[ecol], errors="coerce").values == event_is).astype(int)
    sel = np.isfinite(t) & (t > 0)
    arm = d["arm"][sel]; t, e = t[sel], e[sel]
    print(f"\n{d['label']} {endpoint}: checkpoint n={arm.sum()} (events {e[arm == 1].sum()}), control n={(1 - arm).sum()} "
          f"(events {e[arm == 0].sum()}); excluded (missing or non-positive time): {(~sel).sum()}")
    for k in ORDER:
        s = std(d["S"][k].values[sel])
        b, se = cox_int(s, arm, t, e)
        r = ci(b, se)
        rows.append(dict(trial=d["label"], key=key, endpoint=endpoint, model="Cox", role=role, score=k,
                         n_ici=int(arm.sum()), n_ctrl=int((1 - arm).sum()), ev_ici=int(e[arm == 1].sum()),
                         ev_ctrl=int(e[arm == 0].sum()), int_est=r["est"], int_lo=r["lo"], int_hi=r["hi"], int_p=r["p"],
                         int_log=b, int_se=se, y=-b, method="Cox"))
        for a, lab in [(1, "checkpoint"), (0, "control")]:
            m = arm == a
            bw, sw = cox_within(s[m], t[m], e[m]); w = ci(bw, sw)
            arms.append(dict(trial=d["label"], key=key, endpoint=endpoint, model="unadjusted", score=k, arm=lab,
                             n=int(m.sum()), events=int(e[m].sum()), est=w["est"], lo=w["lo"], hi=w["hi"], p=w["p"],
                             log=bw, se=sw, method="Cox"))


def binary(key, role, adjust, model_tag):
    d = TR[key]; y, arm = d["y"], d["arm"].astype(float)
    cov = d["batch"] if adjust else None
    print(f"\n{d['label']} pCR [{model_tag}]: checkpoint n={int(arm.sum())} (pCR {int(y[arm == 1].sum())}), "
          f"control n={int((1 - arm).sum())} (pCR {int(y[arm == 0].sum())})")
    for k in ORDER:
        s = std(d["S"][k].values)
        b, se, meth, bf, sf = lfit(design(s, arm, cov), y, 3)
        r = ci(b, se)
        rows.append(dict(trial=d["label"], key=key, endpoint="pCR", model=model_tag, role=role, score=k,
                         n_ici=int(arm.sum()), n_ctrl=int((1 - arm).sum()), ev_ici=int(y[arm == 1].sum()),
                         ev_ctrl=int(y[arm == 0].sum()), int_est=r["est"], int_lo=r["lo"], int_hi=r["hi"], int_p=r["p"],
                         int_log=b, int_se=se, y=b, method=meth, firth_est=np.exp(bf), firth_p=2 * norm.sf(abs(bf / sf))))
        for a, lab in [(1, "checkpoint"), (0, "control")]:
            m = arm == a
            cv = cov[m] if cov is not None else None
            bw, sw, mw, _, _ = lfit(design(s[m], arm[m], cv, interaction=False), y[m], 1)
            w = ci(bw, sw)
            arms.append(dict(trial=d["label"], key=key, endpoint="pCR", model=model_tag, score=k, arm=lab,
                             n=int(m.sum()), events=int(y[m].sum()), est=w["est"], lo=w["lo"], hi=w["hi"], p=w["p"],
                             log=bw, se=sw, method=mw))


survival("CM025", "PFS", "PFS_CNSR", "PFS", "primary", 1)
survival("CM025", "OS", "OS_CNSR", "OS", "secondary", 1)
survival("JR101", "PFS_P", "PFS_P_CNSR", "PFS", "primary", 0)
binary("NeoTRIP", "primary", True, "logistic + batch")
binary("NeoTRIP", "sensitivity", False, "logistic")
binary("DUTRENEO", "sensitivity", False, "logistic")
binary("DUTRENEO_1M", "sensitivity (post hoc)", False, "logistic")

R = pd.DataFrame(rows); A = pd.DataFrame(arms)
R.to_csv(f"{O}/E21_interactions.csv", index=False); A.to_csv(f"{O}/E21_within_arm.csv", index=False)
pd.set_option("display.width", 250)
for (tr, ep, mo), g in R.groupby(["trial", "endpoint", "model"], sort=False):
    w = A[(A.trial == tr) & (A.endpoint == ep) & (A.model == ("unadjusted" if mo == "Cox" else mo))]
    wi = w.pivot(index="score", columns="arm", values=["est", "p"]).loc[ORDER]
    g = g.set_index("score").loc[ORDER]
    out = pd.DataFrame({"int": g.int_est, "lo": g.int_lo, "hi": g.int_hi, "P": g.int_p, "y_benefit": g.y, "method": g.method,
                        "ckpt": wi[("est", "checkpoint")], "ckpt_P": wi[("p", "checkpoint")],
                        "ctrl": wi[("est", "control")], "ctrl_P": wi[("p", "control")]})
    print(f"\n=== {tr} | {ep} | {mo} | interaction {'HR' if mo == 'Cox' else 'OR'} per SD; within-arm {'HR' if mo == 'Cox' else 'OR'} ===")
    print(out.round(3).to_string())
print("\nmethods:", R.method.value_counts().to_dict(), "| within-arm:", A.method.value_counts().to_dict())

# expectations, per new trial (primary model; also the sensitivity models)
ex_rows = []
for _, r in R[R.score.isin(EXPECT)].iterrows():
    want = EXPECT[r.score]
    held = r.y > 0 if want == ">0" else r.y <= 0
    ex_rows.append(dict(trial=r.trial, endpoint=r.endpoint, model=r.model, role=r.role, score=r.score, expected=f"y {want}",
                        y=r.y, p=r.int_p, held=bool(held), supported=bool(held and r.int_p < 0.05)))
EX = pd.DataFrame(ex_rows); EX.to_csv(f"{O}/E21_expectations_per_trial.csv", index=False)
print("\n=== expectations per trial (y = benefit-oriented log ratio) ===")
print(EX.round(3).to_string(index=False))

# ------------------------------------------------------------------ step 3: structure-matched null
NEW = ["CM025", "JR101", "NeoTRIP"]
pool_all = open(f"{O}/bulk_gene_pool.txt").read().split()
pool = [g for g in pool_all if all(g in TR[k]["Z"].index for k in NEW)]
pool = [g for g in pool if all(TR[k]["Z"].loc[g].notna().all() and TR[k]["Z"].loc[g].std() > 0 for k in NEW)]
print(f"\nnull gene pool (bulk pool, measured and non-constant in all three new primary trials): {len(pool)} of {len(pool_all)}")
SLOT = {"APM": 8, "IFN": 5, "TLS_B": 6, "T_CELL": 7, "CAF_ECM": 11, "TAM_SUP": 11}
STRUCT = {k: {m: sum(g in TR[k]["Z"].index for g in MOD[m]) for m in SLOT} for k in NEW}
for k in NEW:
    print(f"  {TR[k]['label']}: NetImmune null structure immune",
          [STRUCT[k][m] for m in ["APM", "IFN", "TLS_B", "T_CELL"]], "minus barrier", [STRUCT[k][m] for m in ["CAF_ECM", "TAM_SUP"]],
          "| CAF_ECM null size", STRUCT[k]["CAF_ECM"])

# outcome arrays for the primary models
NUL = {}
for k in NEW:
    d = TR[k]; c = d["clin"]
    if k == "CM025":
        t = c.PFS.values.astype(float); e = (c.PFS_CNSR.values == 1).astype(int)
    elif k == "JR101":
        t = c.PFS_P.values.astype(float); e = (c.PFS_P_CNSR.values == 0).astype(int)
    if k in ("CM025", "JR101"):
        sel = np.isfinite(t) & (t > 0)
        NUL[k] = dict(kind="cox", sel=sel, t=t[sel], e=e[sel], arm=d["arm"][sel], ZP=d["Z"].loc[pool].values[:, sel])
    else:
        NUL[k] = dict(kind="logit", y=d["y"], arm=d["arm"].astype(float), cov=d["batch"], ZP=d["Z"].loc[pool].values)


def null_stat(k, s):
    n = NUL[k]; s = std(s)
    if n["kind"] == "cox":
        b, se = cox_int(s, n["arm"], n["t"], n["e"]); return -b, se
    b, V_, ok = ml(design(s, n["arm"], n["cov"]), n["y"])
    return b[3], np.sqrt(V_[3, 3])


def modscore(Zp, idx, sizes):
    return np.mean([Zp[ix].mean(0) for ix in (idx[m][:sizes[m]] for m in sizes)], 0)


DRAWS = []
for i in range(NRAND):
    g = rng.choice(len(pool), sum(SLOT.values()), replace=False)
    g11 = rng.choice(len(pool), 11, replace=False)
    idx, j = {}, 0
    for m, n in SLOT.items():
        idx[m] = g[j:j + n]; j += n
    DRAWS.append((idx, g11))


def one_draw(i):
    idx, g11 = DRAWS[i]; row = {}
    for k in NEW:
        Zp = NUL[k]["ZP"]; st = STRUCT[k]
        imm = modscore(Zp, idx, {m: st[m] for m in ["APM", "IFN", "TLS_B", "T_CELL"]})
        bar = modscore(Zp, idx, {m: st[m] for m in ["CAF_ECM", "TAM_SUP"]})
        row[f"NetImmune_{k}"], row[f"NetImmune_{k}_se"] = null_stat(k, imm - bar)
        row[f"CAF_ECM_{k}"], row[f"CAF_ECM_{k}_se"] = null_stat(k, Zp[g11[:st["CAF_ECM"]]].mean(0))
    return row


def ivw(bs, ses):
    w = 1 / np.asarray(ses) ** 2
    return np.sum(w * np.asarray(bs)) / w.sum(), np.sqrt(1 / w.sum())


with mp.get_context("fork").Pool(24) as pool_:
    draws = pool_.map(one_draw, range(NRAND), chunksize=20)
ND = pd.DataFrame(draws)
for k in ["NetImmune", "CAF_ECM"]:
    cb = [ivw([r[f"{k}_{t}"] for t in NEW], [r[f"{k}_{t}_se"] for t in NEW]) for _, r in ND.iterrows()]
    ND[f"{k}_combined"] = [x[0] for x in cb]; ND[f"{k}_combined_se"] = [x[1] for x in cb]
ND.to_csv(f"{O}/E21_null_draws.csv", index=False)
prim = {"CM025": ("PFS", "Cox"), "JR101": ("PFS", "Cox"), "NeoTRIP": ("pCR", "logistic + batch")}
nrows = []
for k in ["NetImmune", "CAF_ECM"]:
    vals = {}
    for t in NEW:
        o = R[(R.key == t) & (R.endpoint == prim[t][0]) & (R.model == prim[t][1]) & (R.score == k)].iloc[0]
        vals[t] = (o.y, o.int_se)
    vals["combined"] = ivw([v[0] for v in vals.values()], [v[1] for v in vals.values()])
    for t, (b, se) in vals.items():
        nd = ND[f"{k}_{t}"].values
        nrows.append(dict(score=k, contrast=t, label=TR[t]["label"] if t in TR else "three new trials, inverse-variance",
                          obs_y=b, obs_se=se, null_median_abs=np.median(np.abs(nd)), null_p95_abs=np.quantile(np.abs(nd), 0.95),
                          p_emp=((np.abs(nd) >= abs(b)).sum() + 1) / (NRAND + 1), frac_null_ge_obs_signed=(nd >= b).mean()))
NS = pd.DataFrame(nrows); NS.to_csv(f"{O}/E21_null_summary.csv", index=False)
print("\n=== structure-matched null (primary-model benefit-oriented interaction y; direction credited) ===")
print(NS.round(4).to_string(index=False))
