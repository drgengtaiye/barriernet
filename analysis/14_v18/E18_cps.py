"""E18 - does NetImmune add to MSI/EBV status and PD-L1 CPS in PRJEB25780? (pre-registered: E17_E22_prereg.md)

CPS (22C3 combined positive score) was transcribed from Kim et al. 2018 Table 2
(data_download/validation/kim2018_cps/); 40 of the 44 analysed patients have it.
Primary: Firth logistic response ~ MSI/EBV + log(CPS + 1) + NetImmune, penalised likelihood-ratio
test for NetImmune (restricted fit keeps the full design in the penalty, as logistf and E7c).
Secondary: CPS >= 1 in place of log(CPS + 1). Scores are the log2-scale scores of E7a, re-standardised
within the 40 patients. Firth fitting is copied from E7c_firth_combo.py.
"""
import numpy as np, pandas as pd
from scipy.stats import chi2, spearmanr
from sklearn.metrics import roc_auc_score

G = "/groups/adv2105_gp/yichen/Yi/com/gastric/2"
O = f"{G}/analysis/results"
rng = np.random.default_rng(2026)


def firth(X, y, fixed=None, tol=1e-9, it=200):
    n, k = X.shape
    fixed = np.zeros(k, bool) if fixed is None else np.asarray(fixed)
    free = ~fixed
    b = np.zeros(k)
    def pll(b):
        p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
        sign, logdet = np.linalg.slogdet(X.T @ (X * W[:, None]))
        return np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)) + 0.5 * logdet
    for _ in range(it):
        p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
        I = X.T @ (X * W[:, None]); Ii = np.linalg.inv(I)
        h = np.einsum("ij,jk,ik->i", X * np.sqrt(W)[:, None], Ii, X * np.sqrt(W)[:, None])
        U = X.T @ (y - p + h * (0.5 - p))
        step = np.zeros(k)
        step[free] = np.linalg.solve(I[np.ix_(free, free)], U[free])
        f0, t = pll(b), 1.0
        while t > 1e-6 and pll(b + t * step) < f0 - 1e-12:
            t /= 2
        b = b + t * step
        if np.max(np.abs(t * step)) < tol:
            break
    p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    return b, cov, pll(b)


def plrt(X, y, j):
    b, cov, ll1 = firth(X, y)
    m = np.zeros(X.shape[1], bool); m[j] = True
    _, _, ll0 = firth(X, y, fixed=m)
    stat = 2 * (ll1 - ll0)
    return b[j], np.sqrt(cov[j, j]), stat, chi2.sf(stat, 1)


def boot_auc(y, s, B=2000):
    out = []
    i1, i0 = np.where(y == 1)[0], np.where(y == 0)[0]
    for _ in range(B):
        idx = np.r_[rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]
        out.append(roc_auc_score(y[idx], s[idx]))
    return np.percentile(out, [2.5, 97.5])


d = pd.read_csv(f"{O}/E7_scores_log2.csv")
c = pd.read_csv(f"{G}/data_download/validation/kim2018_cps/kim2018_cps_clinical.csv")
m = d.merge(c[["seq_id", "patient_id", "cps", "pdl1_cps_ge1", "responder", "msi_h", "ebv_pos"]],
            left_on="SEQ_ID", right_on="seq_id", how="left")
assert len(m) == 44 and (m.Responder == m.responder).all()
assert ((m.MSI_EBV == 1) == (m.msi_h | m.ebv_pos)).all()
m = m[m.cps.notna()].copy()
print(f"patients with CPS: {len(m)}; responders {int(m.Responder.sum())}; "
      f"CPS>=1 {int((m.cps >= 1).sum())}; MSI/EBV {int(m.MSI_EBV.sum())}")

y = m.Responder.values.astype(float)
z = lambda v: (v - v.mean()) / v.std()
net = z(m.NetImmune.values)
lcps = np.log1p(m.cps.values)
rows = []

# correlations
for name, v in [("CPS", m.cps.values), ("log(CPS+1)", lcps)]:
    r, p = spearmanr(m.NetImmune, v)
    rows.append(dict(analysis="Spearman with NetImmune", term=name, estimate=r, p=p, n=len(m)))

# AUCs in the same 40 patients
for name, v in [("NetImmune", net), ("log(CPS+1)", lcps), ("CPS>=1", (m.cps.values >= 1).astype(float)),
                ("MSI/EBV", m.MSI_EBV.values.astype(float))]:
    a = roc_auc_score(y, v); lo, hi = boot_auc(y, v)
    rows.append(dict(analysis="AUC (n = 40)", term=name, estimate=a, lo=lo, hi=hi, n=len(m)))

# Firth models
one = np.ones(len(y))
for label, cpsvar in [("primary: log(CPS+1)", z(lcps)), ("secondary: CPS>=1", (m.cps.values >= 1).astype(float))]:
    X = np.c_[one, m.MSI_EBV.values.astype(float), cpsvar, net]
    b, se, stat, p = plrt(X, y, 3)
    rows.append(dict(analysis=f"Firth: response ~ MSI/EBV + CPS + NetImmune ({label})", term="NetImmune per SD",
                     estimate=b, lo=b - 1.96 * se, hi=b + 1.96 * se, OR=np.exp(b), chi2=stat, p=p, n=len(m)))
    bc, sec, statc, pc = plrt(X, y, 2)
    rows.append(dict(analysis=f"Firth: response ~ MSI/EBV + CPS + NetImmune ({label})", term="CPS term",
                     estimate=bc, lo=bc - 1.96 * sec, hi=bc + 1.96 * sec, OR=np.exp(bc), chi2=statc, p=pc, n=len(m)))
# the earlier model (MSI/EBV + NetImmune) in the same 40, for comparison
X = np.c_[one, m.MSI_EBV.values.astype(float), net]
b, se, stat, p = plrt(X, y, 2)
rows.append(dict(analysis="Firth: response ~ MSI/EBV + NetImmune (same 40 patients)", term="NetImmune per SD",
                 estimate=b, lo=b - 1.96 * se, hi=b + 1.96 * se, OR=np.exp(b), chi2=stat, p=p, n=len(m)))

# descriptive: responders by CPS and MSI/EBV
tab = pd.crosstab([m.MSI_EBV.rename("MSI/EBV"), (m.cps >= 1).rename("CPS>=1")], m.Responder.rename("responder"))
print(tab)
res = pd.DataFrame(rows)
res.to_csv(f"{O}/E18_cps.csv", index=False)
tab.to_csv(f"{O}/E18_cps_crosstab.csv")
m[["SEQ_ID", "patient_id", "cps", "MSI_EBV", "Responder", "NetImmune"]].to_csv(f"{O}/E18_cps_patients.csv", index=False)
pd.set_option("display.width", 200)
print(res.round(4).to_string(index=False))
