"""E10 - score x arm interaction for pathological complete response in two randomised I-SPY2
comparisons (HER2-negative breast cancer, neoadjuvant), to add to IMmotion150 (E5_trials).

Data (built by submission/validation/ispy2_prep.py; no new downloads):
  ISPY2_pembro  GSE194040 [@wolf2022; @nanda2020]: paclitaxel + pembrolizumab (n 69) vs the
                HER2-negative paclitaxel controls (n 179); fresh-frozen Agilent gene-level log2.
  ISPY2_durva   GSE173839 [@pusztai2021]: durvalumab + olaparib + paclitaxel (n 71) vs its
                paclitaxel controls (n 29), evaluable pCR only; FFPE Agilent gene-level log2.
  No research ID is shared between the two files (the 29 durvalumab-comparison controls are not
  among the 988 GSE194040 patients), so the two comparisons have no patients in common.

Scoring exactly as in the manuscript: genes row-z within cohort (each contrast file), module =
mean z of its genes present (>= 3), Barrier = mean(CAF_ECM, TAM_SUP), Immune = mean(APM, IFN,
TLS_B, T_CELL), NetImmune = Immune - Barrier, C1Q / SPP1 macrophage arms, CXCL9:SPP1 =
z(CXCL9) - z(SPP1). TRAC and TRBC1 are not on the Agilent array, so T_CELL has 5 genes here.
Each score is standardised to SD 1 over the patients of the contrast. CXCL9 or SPP1 is missing
in 3 patients of the pembrolizumab contrast, who are left out of the CXCL9:SPP1 fits only.

Model: logistic pCR ~ s + arm + s:arm (arm 1 = checkpoint-containing). Interaction OR =
exp(coef of s:arm); > 1 means the score predicts pCR more strongly with the checkpoint
inhibitor (benefit direction). Fitted three ways:
  overall        - as above;
  HR-adjusted    - + HR status (primary: I-SPY2 randomisation is adaptive by subtype, and HR+
                   patients are over-represented in both experimental arms);
  HR-stratified  - the model fitted within HR- and HR+ separately, the two interaction
                   log-ORs combined by inverse variance (sensitivity).
Within-arm OR per SD from separate fits in each arm (unadjusted and HR-adjusted). Firth
penalisation (Jeffreys prior, as in E7c) replaces ML wherever ML shows separation
(non-convergence, |coef| > 10 or SE > 10); Firth estimates are also written for every fit as
a check.

Structure-matched random-gene null (NetImmune and CAF_ECM interactions): 2,000 random gene
sets drawn from the bulk gene pool measured (no missing values) in both contrasts, the same
set applied to both. NetImmune: four "immune" modules of 8/5/6/5 genes minus two "barrier"
modules of 11, all distinct; CAF_ECM: 11 genes. Statistic: the HR-adjusted interaction log-OR
in each contrast and their inverse-variance combination; random sets credited with their
better direction, P = (#|null| >= |obs| + 1) / (N + 1).

Pre-stated expectation (written before running), from the gastric, 11-cohort and IMmotion150
results. "Held" = point estimate on the expected side; "supported" = also P < 0.05.
  IFN, Immune, NetImmune - benefit-direction interaction with the checkpoint arm (OR > 1)
  CAF_ECM                - no interaction (OR near 1; 95% CI including 1)
  TAM_SUP                - leans towards the immune side (OR > 1, weaker than the immune modules)
  CXCL9_SPP1             - behaves as the immune modules do (OR > 1)
  (C1Q arm, SPP1 arm, APM, TLS_B, T_CELL, Barrier are reported without a stated expectation
  beyond the E5 ones: C1Q arm OR > 1; SPP1 arm no benefit-direction interaction.)
"""
import warnings, functools; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from scipy.stats import norm

V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation/cohorts"
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
TRIALS = {"ISPY2_pembro": "I-SPY2: pembrolizumab + paclitaxel vs paclitaxel",
          "ISPY2_durva": "I-SPY2: durvalumab + olaparib + paclitaxel vs paclitaxel"}


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


def ml(X, y, it=100, tol=1e-10):
    """Newton-Raphson logistic ML; returns b, cov, converged."""
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
    """Firth (Jeffreys-penalised) logistic regression, as in E7c_firth_combo.py."""
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


def fit(X, y, j):
    """Coefficient j: (estimate, SE, method), Firth if ML separates."""
    b, V, ok = ml(X, y)
    sep = (not ok) or np.any(np.abs(b) > 10) or np.any(~np.isfinite(np.diag(V))) or np.any(np.sqrt(np.diag(V)) > 10)
    bf, Vf = firth(X, y)
    if sep:
        return bf[j], np.sqrt(Vf[j, j]), "Firth", bf[j], np.sqrt(Vf[j, j])
    return b[j], np.sqrt(V[j, j]), "ML", bf[j], np.sqrt(Vf[j, j])


def est(b, se):
    return dict(logOR=b, se=se, OR=np.exp(b), lo=np.exp(b - 1.96 * se), hi=np.exp(b + 1.96 * se),
                p=2 * norm.sf(abs(b / se)))


def interaction(s, arm, y, hr, adjust):
    one = np.ones_like(s)
    X = np.c_[one, s, arm, s * arm] if not adjust else np.c_[one, s, arm, s * arm, hr]
    return fit(X, y, 3)


def within(s, y, hr, adjust):
    X = np.c_[np.ones_like(s), s] if not adjust else np.c_[np.ones_like(s), s, hr]
    return fit(X, y, 1)


def ivw(bs, ses):
    w = 1 / np.asarray(ses) ** 2
    return np.sum(w * bs) / w.sum(), np.sqrt(1 / w.sum())


D, rows, arms_rows, strat_rows = {}, [], [], []
for key, label in TRIALS.items():
    X = pd.read_csv(f"{V}/{key}.expr.tsv.gz", sep="\t", index_col=0); X.index.name = None
    X.columns = X.columns.astype(str)
    C = pd.read_csv(f"{V}/{key}.clin.tsv", sep="\t"); C["sample"] = C["sample"].astype(str)
    C = C.set_index("sample").loc[X.columns]
    Z = ((X.T - X.mean(1)) / X.std(1).replace(0, 1)).T
    S = scores(Z)
    print(f"\n=== {label} ===")
    print("module genes present:", {k: sum(g in Z.index for g in G) for k, G in MOD.items()})
    print("n / pCR by arm:", {a: (len(g), int(g.y.sum())) for a, g in C.groupby("arm")},
          "| HR+ by arm:", {a: int(g.hr.sum()) for a, g in C.groupby("arm")})
    print("pCR by arm x HR:", C.groupby(["arm", "hr"]).y.agg(["size", "sum"]).to_dict("index"))
    y0, arm0, hr0 = C.y.values.astype(float), C.arm.values.astype(float), C.hr.values.astype(float)
    D[key] = dict(Z=Z, y=y0, arm=arm0, hr=hr0)
    for k in ORDER:
        ok = np.isfinite(S[k].values)       # CXCL9 or SPP1 missing in 3 pembrolizumab-contrast patients
        y, arm, hr = y0[ok], arm0[ok], hr0[ok]
        s = S[k].values[ok]; s = (s - s.mean()) / s.std(ddof=1)
        r = dict(trial=label, cohort=key, score=k, n_ici=int(arm.sum()), n_ctrl=int((1 - arm).sum()),
                 pcr_ici=int(y[arm == 1].sum()), pcr_ctrl=int(y[arm == 0].sum()))
        for adj, tag in [(False, "overall"), (True, "HRadj")]:
            b, se, m, bf, sef = interaction(s, arm, y, hr, adj)
            r.update({f"{tag}_{a}": v for a, v in est(b, se).items()}); r[f"{tag}_method"] = m
            r[f"{tag}_firth_OR"] = np.exp(bf); r[f"{tag}_firth_p"] = 2 * norm.sf(abs(bf / sef))
        # HR-stratified: interaction within HR- and HR+, combined by inverse variance
        sb, ss = [], []
        for h in (0, 1):
            m_ = hr == h
            b, se, m, bf, sef = interaction(s[m_], arm[m_], y[m_], hr[m_], False)
            sb.append(b); ss.append(se)
            strat_rows.append(dict(trial=label, cohort=key, score=k, stratum=f"HR{'+' if h else '-'}",
                                   n=int(m_.sum()), pcr=int(y[m_].sum()), method=m, **est(b, se)))
        b, se = ivw(sb, ss)
        r.update({f"strat_{a}": v for a, v in est(b, se).items()})
        rows.append(r)
        for a, lab in [(1, "checkpoint"), (0, "control")]:
            m_ = arm == a
            for adj, tag in [(False, "unadjusted"), (True, "HRadj")]:
                b, se, m, bf, sef = within(s[m_], y[m_], hr[m_], adj)
                arms_rows.append(dict(trial=label, cohort=key, score=k, arm=lab, model=tag, n=int(m_.sum()),
                                      pcr=int(y[m_].sum()), method=m, **est(b, se)))

R = pd.DataFrame(rows); A = pd.DataFrame(arms_rows); ST = pd.DataFrame(strat_rows)
R.to_csv(f"{O}/E10_ispy2_interactions.csv", index=False)
A.to_csv(f"{O}/E10_ispy2_within_arm.csv", index=False)
ST.to_csv(f"{O}/E10_ispy2_hr_strata.csv", index=False)
pd.set_option("display.width", 250)
for key, g in R.groupby("cohort", sort=False):
    print(f"\n--- {key}: interaction OR per SD (overall | HR-adjusted | HR-stratified) ---")
    print(g[["score", "overall_OR", "overall_p", "HRadj_OR", "HRadj_lo", "HRadj_hi", "HRadj_p", "HRadj_method",
             "HRadj_firth_OR", "strat_OR", "strat_p"]].round(3).to_string(index=False))
    a = A[(A.cohort == key) & (A.model == "HRadj")].pivot(index="score", columns="arm", values=["OR", "p"]).loc[ORDER]
    print("within-arm OR per SD (HR-adjusted):\n", a.round(3).to_string())
print("\nmethods used:", R.overall_method.value_counts().to_dict(), R.HRadj_method.value_counts().to_dict(),
      ST.method.value_counts().to_dict(), A.method.value_counts().to_dict())

# ------------------------------------------------------------------ structure-matched random-gene null
pool_all = open(f"{O}/bulk_gene_pool.txt").read().split()
pool = [g for g in pool_all if all(g in d["Z"].index and d["Z"].loc[g].notna().all() for d in D.values())]
print(f"\nnull gene pool (bulk pool, measured without gaps in both contrasts): {len(pool)}")
ZP = {k: d["Z"].loc[pool].values for k, d in D.items()}
SI = [sum(g in D["ISPY2_pembro"]["Z"].index for g in MOD[m]) for m in ["APM", "IFN", "TLS_B", "T_CELL"]]
SB = [11, 11]
print("NetImmune null structure: immune", SI, "minus barrier", SB)


def modscore(Zp, idx, sizes):
    out, k = [], 0
    for s_ in sizes:
        out.append(Zp[idx[k:k + s_]].mean(0)); k += s_
    return np.mean(out, 0)


def int_hradj(s, d):
    s = (s - s.mean()) / s.std(ddof=1)
    X = np.c_[np.ones_like(s), s, d["arm"], s * d["arm"], d["hr"]]
    b, V, ok = ml(X, d["y"])
    return b[3], np.sqrt(V[3, 3])


draws = []
for i in range(NRAND):
    g = rng.choice(len(pool), sum(SI) + sum(SB), replace=False)
    g11 = rng.choice(len(pool), 11, replace=False)
    row = {}
    for key, d in D.items():
        net = modscore(ZP[key], g[:sum(SI)], SI) - modscore(ZP[key], g[sum(SI):], SB)
        caf = ZP[key][g11].mean(0)
        row[f"NetImmune_{key}"], row[f"NetImmune_{key}_se"] = int_hradj(net, d)
        row[f"CAF_ECM_{key}"], row[f"CAF_ECM_{key}_se"] = int_hradj(caf, d)
    for k in ["NetImmune", "CAF_ECM"]:
        row[f"{k}_combined"] = ivw([row[f"{k}_ISPY2_pembro"], row[f"{k}_ISPY2_durva"]],
                                   [row[f"{k}_ISPY2_pembro_se"], row[f"{k}_ISPY2_durva_se"]])[0]
    draws.append(row)
ND = pd.DataFrame(draws)
ND.to_csv(f"{O}/E10_ispy2_null_draws.csv", index=False)
nrows = []
for k in ["NetImmune", "CAF_ECM"]:
    obs = {key: R[(R.cohort == key) & (R.score == k)].iloc[0] for key in TRIALS}
    vals = {key: (o.HRadj_logOR, o.HRadj_se) for key, o in obs.items()}
    vals["combined"] = ivw([v[0] for v in vals.values()], [v[1] for v in vals.values()])
    for key, (b, se) in vals.items():
        nd = ND[f"{k}_{key}"].values
        nrows.append(dict(score=k, contrast=key, obs_logOR=b, obs_se=se, null_median_abs=np.median(np.abs(nd)),
                          null_p95_abs=np.quantile(np.abs(nd), 0.95),
                          p_emp=((np.abs(nd) >= abs(b)).sum() + 1) / (NRAND + 1),
                          frac_null_ge_obs_signed=(nd >= b).mean()))
NS = pd.DataFrame(nrows)
NS.to_csv(f"{O}/E10_ispy2_null_summary.csv", index=False)
print("\n=== structure-matched null (HR-adjusted interaction log-OR; direction credited) ===")
print(NS.round(4).to_string(index=False))
