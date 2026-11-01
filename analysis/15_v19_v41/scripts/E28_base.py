"""Shared base for E28 (E28_prereg.md, 2026-10-20): imports, modules, scoring, models, pooling rule and the
six primary randomized comparisons, copied from E25_trials.py lines 16-307. Two changes: each comparison also
keeps its log expression matrix (X, same samples as Z), and the step-0 reproduction table is written as
E28_reproduction_check.csv. Step 0 still stops if any published interaction or within-arm value differs.
"""
import warnings, functools; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy.stats import norm, t as tdist
from scipy.optimize import minimize_scalar
from lifelines import CoxPHFitter
import sys as _sys; _sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis")
from cmp_scores import comparators as _cmp

O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
P5 = "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"                              # E5
V10 = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation/cohorts"                 # E10
V21 = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/validation"         # E21
rng = np.random.default_rng(2025)
MOD = dict(
    CAF_ECM=["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
    TAM_SUP=["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
    APM=["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    IFN=["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
    TLS_B=["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
    T_CELL=["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"],
    C1Q_arm=["C1QA", "C1QB", "C1QC"],
    SPP1_arm=["SPP1", "TREM2", "MARCO", "APOE", "LGALS3"])
CHECK = ["CAF_ECM", "TAM_SUP", "IFN", "APM", "TLS_B", "T_CELL", "Barrier", "Immune", "NetImmune", "ESTIMATE", "CD8A"]
Q1 = ["CAF_ECM", "TAM_SUP", "Barrier", "IFN", "APM", "TLS_B", "T_CELL", "Immune", "NetImmune", "ESTIMATE",
      "TMEscore_PC", "CD8A"]
BREAST = ["ISPY2_pembro", "ISPY2_durva", "NeoTRIP"]; RENAL = ["IMM150", "CM025", "JR101"]
ALL6 = ["IMM150", "ISPY2_pembro", "ISPY2_durva", "CM025", "JR101", "NeoTRIP"]
LABEL = {"IMM150": "IMmotion150", "ISPY2_pembro": "I-SPY2 pembrolizumab", "ISPY2_durva": "I-SPY2 durvalumab-olaparib",
         "CM025": "CheckMate 025", "JR101": "JAVELIN Renal 101", "NeoTRIP": "NeoTRIP"}


# ------------------------------------------------------------------ scoring and models (E21, unchanged)
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
    return pd.concat([S, _cmp(Z)], axis=1)


def rowz(X):
    X = X[X.std(1) > 0]
    return ((X.T - X.mean(1)) / X.std(1)).T


def std(s):
    s = np.asarray(s, float)
    return (s - s.mean()) / s.std(ddof=1)


def collapse(X):
    nd = int(X.index.duplicated().sum())
    return (X.groupby(level=0).mean() if nd else X), nd


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


def lfit(X, y, js, firth_ok=True):
    """Coefficients js: list of (estimate, SE); Firth if ML separates (observed data only)."""
    b, V_, ok = ml(X, y)
    sep = (not ok) or np.any(np.abs(b) > 10) or np.any(~np.isfinite(np.diag(V_))) or np.any(np.sqrt(np.diag(V_)) > 10)
    if sep and firth_ok:
        b, V_ = firth(X, y)
    elif sep:
        return [(np.nan, np.nan) for _ in js]
    return [(b[j], np.sqrt(V_[j, j])) for j in js]


def cox(cols, t, e, names):
    df = pd.DataFrame(cols); df["T"] = t; df["E"] = e
    r = CoxPHFitter().fit(df, "T", "E").summary
    return [(r.loc[n, "coef"], r.loc[n, "se(coef)"]) for n in names]


# ------------------------------------------------------------------ trial models (benefit-oriented y)
def m_int(d, s):
    """Interaction y, se for standardized s with the comparison's primary model."""
    s = std(s); a = d["arm"]
    if d["kind"] == "cox":
        b, se = cox({"s": s, "arm": a, "sx": s * a}, d["t"], d["e"], ["sx"])[0]; return -b, se
    X = np.column_stack([np.ones_like(s), s, a, s * a] + ([d["cov"]] if d["cov"] is not None else []))
    return lfit(X, d["y"], [3], d["firth"])[0]


def m_within(d, s, arm):
    """Within-arm y, se per SD (s standardized over the whole comparison, as in E5/E10/E21)."""
    s = std(s); m = d["arm"] == arm
    if d["kind"] == "cox":
        b, se = cox({"s": s[m]}, d["t"][m], d["e"][m], ["s"])[0]; return -b, se
    cov = d["cov"][m] if d["cov"] is not None else None
    X = np.column_stack([np.ones(m.sum()), s[m]] + ([cov] if cov is not None else []))
    return lfit(X, d["y"][m], [1], d["firth"])[0]


def m_joint(d, si, sc):
    """Immune and CAF_ECM in one model with both products: adjusted interaction y, se for each."""
    si, sc = std(si), std(sc); a = d["arm"]
    if d["kind"] == "cox":
        (bi, sei), (bc, sec) = cox({"si": si, "sc": sc, "arm": a, "xi": si * a, "xc": sc * a}, d["t"], d["e"], ["xi", "xc"])
        return -bi, sei, -bc, sec
    X = np.column_stack([np.ones_like(si), si, sc, a, si * a, sc * a] + ([d["cov"]] if d["cov"] is not None else []))
    (bi, sei), (bc, sec) = lfit(X, d["y"], [4, 5], d["firth"])
    return bi, sei, bc, sec


# ------------------------------------------------------------------ pooling: manuscript rule in numpy
def _reml_tau2(y, v):
    def nll(t2):
        w = 1 / (v + t2); mu = np.sum(w * y) / w.sum()
        return 0.5 * (np.sum(np.log(v + t2)) + np.log(w.sum()) + np.sum(w * (y - mu) ** 2))
    hi = max(10 * np.var(y, ddof=1), 10 * v.max(), 1e-4)
    r = minimize_scalar(nll, bounds=(0, hi), method="bounded", options=dict(xatol=1e-12))
    return r.x if nll(r.x) < nll(0.0) else 0.0


def pool_rule(y, se):
    y, se = np.asarray(y, float), np.asarray(se, float)
    ok = np.isfinite(y) & np.isfinite(se) & (se > 0); y, se = y[ok], se[ok]; k = len(y)
    if k < 2:
        return dict(est=np.nan, lo=np.nan, hi=np.nan, p=np.nan, k=k)
    v = se ** 2
    t2 = _reml_tau2(y, v); w = 1 / (v + t2); mu = np.sum(w * y) / w.sum()
    se_hk = np.sqrt(np.sum(w * (y - mu) ** 2) / ((k - 1) * w.sum())); q = tdist.ppf(0.975, k - 1)
    hk_lo, hk_hi, hk_p = mu - q * se_hk, mu + q * se_hk, 2 * tdist.sf(abs(mu / se_hk), k - 1)
    w0 = 1 / v; m0 = np.sum(w0 * y) / w0.sum(); Q = np.sum(w0 * (y - m0) ** 2)
    t2d = max(0.0, (Q - (k - 1)) / (w0.sum() - np.sum(w0 ** 2) / w0.sum())); wd = 1 / (v + t2d)
    md = np.sum(wd * y) / wd.sum(); sd = np.sqrt(1 / wd.sum())
    zq = norm.ppf(0.975); dl_lo, dl_hi, dl_p = md - zq * sd, md + zq * sd, 2 * norm.sf(abs(md / sd))
    hk_wider = (hk_hi - hk_lo) >= (dl_hi - dl_lo)
    return dict(est=mu, lo=hk_lo if hk_wider else dl_lo, hi=hk_hi if hk_wider else dl_hi, p=max(hk_p, dl_p), k=k)


# ------------------------------------------------------------------ data (copied from E5, E10, E21)
TR = {}

print("=" * 90, "\nIMmotion150 (E5)")
C = pd.read_csv(f"{P5}/v4_immotion150_clinical_api.csv", index_col=0)
M = pd.read_csv(f"{P5}/ici/rcc_iatlas_immotion150_2018.csv", index_col=0)
keep = [s for s in M.columns if s in C.index]; C, M = C.loc[keep], M[keep]
Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
S = scores(Z)
t = pd.to_numeric(C["PFS_MONTHS"], errors="coerce")
e = pd.to_numeric(C["PFS_STATUS"].astype(str).str.extract(r"^\s*([01])\s*:")[0], errors="coerce").fillna(0).astype(int)
atezo = C["ICI_RX"].astype(str) == "Atezolizumab"
suni = C["NON_ICI_RX"].astype(str).str.contains("Sunitinib", na=False)
sel = ((atezo | suni) & t.notna() & (t > 0)).values
TR["IMM150"] = dict(kind="cox", Z=Z.loc[:, sel], X=M.loc[:, sel], S=S.loc[sel], arm=atezo.values[sel].astype(int),
                    t=t.values[sel].astype(float), e=e.values[sel].astype(int), firth=True)

for key in ["ISPY2_pembro", "ISPY2_durva"]:
    print("=" * 90, f"\n{key} (E10)")
    X = pd.read_csv(f"{V10}/{key}.expr.tsv.gz", sep="\t", index_col=0); X.index.name = None
    X.columns = X.columns.astype(str)
    C = pd.read_csv(f"{V10}/{key}.clin.tsv", sep="\t"); C["sample"] = C["sample"].astype(str)
    C = C.set_index("sample").loc[X.columns]
    Z = ((X.T - X.mean(1)) / X.std(1).replace(0, 1)).T
    S = scores(Z)
    TR[key] = dict(kind="logit", Z=Z, X=X, S=S, arm=C.arm.values.astype(float), y=C.y.values.astype(float),
                   cov=C.hr.values.astype(float), firth=True)

print("=" * 90, "\nCheckMate 025 (E21)")
C = pd.read_csv(f"{V21}/Braun2020_CheckMate/S1_clinical.tsv", sep="\t")
X = pd.read_csv(f"{V21}/Braun2020_CheckMate/S4A_RNA_expression.tsv.gz", sep="\t", index_col=0)
X, _ = collapse(X)
c = C[C.Cohort == "CM-025"]; c = c[c.RNA_ID.isin(X.columns)].set_index("RNA_ID")
Z = rowz(X[c.index]); S = scores(Z)
t = pd.to_numeric(c.PFS, errors="coerce").values; e = (pd.to_numeric(c.PFS_CNSR, errors="coerce").values == 1).astype(int)
sel = np.isfinite(t) & (t > 0)
TR["CM025"] = dict(kind="cox", Z=Z.loc[:, sel], X=X[c.index].loc[:, sel], S=S.loc[sel], arm=(c.Arm == "NIVOLUMAB").astype(int).values[sel],
                   t=t[sel], e=e[sel], firth=True)

print("=" * 90, "\nJAVELIN Renal 101 (E21)")
C = pd.read_csv(f"{V21}/JAVELIN_Renal101/JR101_S11_clinical.tsv", sep="\t")
X = pd.read_csv(f"{V21}/JAVELIN_Renal101/JR101_S13_expression_log2TPM.tsv.gz", sep="\t", index_col=0)
X, _ = collapse(X)
c = C.set_index("ID").loc[X.columns]
Z = rowz(X); S = scores(Z)
t = pd.to_numeric(c.PFS_P, errors="coerce").values; e = (pd.to_numeric(c.PFS_P_CNSR, errors="coerce").values == 0).astype(int)
sel = np.isfinite(t) & (t > 0)
TR["JR101"] = dict(kind="cox", Z=Z.loc[:, sel], X=X.loc[:, sel], C=c.loc[sel], S=S.loc[sel], arm=(c.TRT01P == "Avelumab+Axitinib").astype(int).values[sel],
                   t=t[sel], e=e[sel], firth=True)

print("=" * 90, "\nNeoTRIP (E21)")
Pn = pd.read_csv(f"{V21}/NeoTRIP_GSE319641/github_Early_GEP_NeoTRIP/data/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt", sep="\t")
X = pd.read_csv(f"{V21}/NeoTRIP_GSE319641/GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz", sep="\t")
X = X.drop(columns=["HGNC_id", "entrez_gene_id"]).set_index("HGNC_approved_symbol")
c = Pn[Pn.Timepoint == "Baseline"].set_index("Sample_ID")
Xb = np.log2(X[c.index].clip(lower=0) + 1); Xb, _ = collapse(Xb)
Z = rowz(Xb); S = scores(Z)
bat = pd.get_dummies(c.Batch_correlation)[["Batch 1", "Batch 2"]].astype(float).values
TR["NeoTRIP"] = dict(kind="logit", Z=Z, X=Xb, S=S, arm=(c.Arm == "CT/A").astype(int).values.astype(float),
                     y=c["pCR.num"].astype(float).values, cov=bat, firth=True)

for k in ALL6:
    d = TR[k]; n1 = int(d["arm"].sum()); n0 = len(d["arm"]) - n1
    print(f"  {LABEL[k]:28s} checkpoint {n1:4d}  control {n0:4d}  genes {d['Z'].shape[0]}")


# ------------------------------------------------------------------ step 0: reproduce published values
def ok_mask(d, k):
    return np.isfinite(d["S"][k].values)


def sub(d, m):
    out = dict(d); out["arm"] = d["arm"][m]
    for f in ["t", "e", "y"]:
        if f in d: out[f] = d[f][m]
    if d.get("cov") is not None: out["cov"] = d["cov"][m]
    return out


E5 = pd.read_csv(f"{O}/E5_trials.csv"); E10 = pd.read_csv(f"{O}/E10_ispy2_interactions.csv")
E21 = pd.read_csv(f"{O}/E21_interactions.csv"); A10 = pd.read_csv(f"{O}/E10_ispy2_within_arm.csv")
A21 = pd.read_csv(f"{O}/E21_within_arm.csv")


def published(k, sc):
    """(interaction y, checkpoint-arm y, control-arm y) as published, benefit-oriented."""
    if k == "IMM150":
        r = E5[(E5.trial == "IMmotion150: atezolizumab (+/- bev) vs sunitinib") & (E5.score == sc)].iloc[0]
        return -np.log(r.int_HR), -np.log(r.ici_HR), -np.log(r.ctrl_HR)
    if k.startswith("ISPY2"):
        r = E10[(E10.cohort == k) & (E10.score == sc)].iloc[0]
        a = A10[(A10.cohort == k) & (A10.score == sc) & (A10.model == "HRadj")].set_index("arm").logOR
        return r.HRadj_logOR, a["checkpoint"], a["control"]
    ep, mod = ("pCR", "logistic + batch") if k == "NeoTRIP" else ("PFS", "Cox")
    r = E21[(E21.key == k) & (E21.endpoint == ep) & (E21.model == mod) & (E21.score == sc)].iloc[0]
    a = A21[(A21.key == k) & (A21.endpoint == ep) & (A21.model == ("logistic + batch" if k == "NeoTRIP" else "unadjusted"))
            & (A21.score == sc)].set_index("arm").log
    sgn = 1 if k == "NeoTRIP" else -1
    return r.y, sgn * a["checkpoint"], sgn * a["control"]


rep, within_rows = [], []
for k in ALL6:
    d = TR[k]
    for sc in sorted(set(CHECK + Q1)):
        m = ok_mask(d, sc); dd = sub(d, m); s = d["S"][sc].values[m]
        yi, si = m_int(dd, s); y1, s1 = m_within(dd, s, 1); y0, s0 = m_within(dd, s, 0)
        pi, p1, p0 = published(k, sc)
        rep.append(dict(trial=LABEL[k], score=sc, int_e25=yi, int_pub=pi, chk_e25=y1, chk_pub=p1, ctl_e25=y0, ctl_pub=p0,
                        maxdiff=max(abs(yi - pi), abs(y1 - p1), abs(y0 - p0))))
        if sc in Q1:
            for a, lab, yy, ss in [(1, "checkpoint", y1, s1), (0, "control", y0, s0)]:
                within_rows.append(dict(trial=LABEL[k], key=k, cancer="breast" if k in BREAST else "renal", score=sc,
                                        arm=lab, y=yy, se=ss, p=2 * norm.sf(abs(yy / ss))))
REP = pd.DataFrame(rep); REP.to_csv(f"{O}/E28_reproduction_check.csv", index=False)
print(f"\nstep 0: {len(REP)} score x comparison values re-derived; largest difference from the published "
      f"interaction / within-arm estimates {REP.maxdiff.max():.2e}")
assert REP.maxdiff.max() < 1e-6, REP.sort_values("maxdiff").tail()
