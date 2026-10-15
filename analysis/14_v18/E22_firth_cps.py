"""E22 - PRJEB60680: does NetImmune add to MSI/EBV status and PD-L1 CPS? (pre-registered: E17_E22_prereg.md, E22)

Same models as E18 (PRJEB25780), with the Firth fitting and penalised likelihood-ratio test copied unchanged from
E18_cps.py / E7c_firth_combo.py (restricted fit keeps the full design in the penalty, as logistf).
Scores: results/E22_scores_PRJEB60680.csv (primary scores, z within library batch; written by E22_gastric_icb.R),
re-standardised within the analysed patients.
MSI/EBV = MSI-H or EBV-positive (Wright 2026 Table S1; 3 patients, all responders -> separation, hence Firth).
CPS = PD-L1 22C3 CPS at baseline (all 46 patients). Response = CR/PR vs SD/PD.
  primary   : response ~ MSI/EBV + z(log(CPS + 1)) + NetImmune      (penalised LRT for NetImmune)
  secondary : CPS >= 1 in place of log(CPS + 1); MSI/EBV + NetImmune; NetImmune alone
  subgroup  : MSS/EBV-negative patients (pre-registered sub-hypothesis NetImmune log-OR > 0):
              response ~ NetImmune and response ~ z(log(CPS + 1)) + NetImmune
  sensitivity (not pre-registered): MSI/EBV also counting E02 (MSI-L by PCR, TCGA subtype MSI in Table S1).
"""
import numpy as np, pandas as pd
from scipy.stats import chi2, spearmanr
from sklearn.metrics import roc_auc_score

G = "/groups/adv2105_gp/yichen/Yi/com/gastric/2"
import os
O = os.environ.get("E22_OUT", f"{G}/analysis/results")
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


# toy check carried over from E7c: Firth log-OR on a 2x2 with a zero cell = Haldane-corrected log-OR
x = np.r_[np.ones(7), np.zeros(13)]; yy = np.r_[np.ones(7), np.ones(3), np.zeros(10)]
bt, _, _ = firth(np.c_[np.ones_like(x), x], yy)
assert abs(bt[1] - np.log(7.5 * 10.5 / (0.5 * 3.5))) < 1e-6

d = pd.read_csv(f"{O}/E22_scores_PRJEB60680.csv")
z = lambda v: (v - v.mean()) / v.std()
MAIN = ["CAF_ECM", "TAM_SUP", "APM", "IFN", "TLS_B", "T_CELL", "Barrier", "Immune", "NetImmune"]
rows = []
print(f"PRJEB60680 (QC-pass): n = {len(d)}, responders {int(d.y.sum())}, MSI/EBV {int(d.msi_ebv.sum())} "
      f"(responders {int(d.y[d.msi_ebv == 1].sum())}), CPS >= 1 {int((d.cps >= 1).sum())}")


def block(m, tag, msi_col="msi_ebv"):
    y = m.y.values.astype(float); one = np.ones(len(m)); lcps = np.log1p(m.cps.values)
    net = z(m.NetImmune.values)
    for name, v in [("CPS", m.cps.values), ("log(CPS+1)", lcps)]:
        r, p = spearmanr(m.NetImmune, v)
        rows.append(dict(set=tag, analysis="Spearman with NetImmune", term=name, estimate=r, p=p, n=len(m)))
    terms = [("NetImmune", net), ("log(CPS+1)", lcps), ("CPS>=1", (m.cps.values >= 1).astype(float))]
    if m[msi_col].nunique() > 1:
        terms.append(("MSI/EBV", m[msi_col].values.astype(float)))
    for name, v in terms:
        a = roc_auc_score(y, v); lo, hi = boot_auc(y, v)
        rows.append(dict(set=tag, analysis=f"AUC (n = {len(m)})", term=name, estimate=a, lo=lo, hi=hi, n=len(m)))
    has_msi = m[msi_col].nunique() > 1
    base = [m[msi_col].values.astype(float)] if has_msi else []
    lab_m = "MSI/EBV + " if has_msi else ""
    for label, cpsvar in [("primary: log(CPS+1)", z(lcps)), ("secondary: CPS>=1", (m.cps.values >= 1).astype(float))]:
        for k in MAIN:
            sv = z(m[k].values)
            X = np.c_[one, *base, cpsvar, sv]
            b, se, stat, p = plrt(X, y, X.shape[1] - 1)
            rows.append(dict(set=tag, analysis=f"Firth: response ~ {lab_m}CPS + score ({label})", term=f"{k} per SD",
                             estimate=b, lo=b - 1.96 * se, hi=b + 1.96 * se, OR=np.exp(b), chi2=stat, p=p, n=len(m)))
            if k == "NetImmune":
                bc, sec, statc, pc = plrt(X, y, X.shape[1] - 2)
                rows.append(dict(set=tag, analysis=f"Firth: response ~ {lab_m}CPS + score ({label})", term="CPS term (NetImmune model)",
                                 estimate=bc, lo=bc - 1.96 * sec, hi=bc + 1.96 * sec, OR=np.exp(bc), chi2=statc, p=pc, n=len(m)))
                if has_msi:
                    bm, sem, statm, pm = plrt(X, y, 1)
                    rows.append(dict(set=tag, analysis=f"Firth: response ~ {lab_m}CPS + score ({label})", term="MSI/EBV term (NetImmune model)",
                                     estimate=bm, lo=bm - 1.96 * sem, hi=bm + 1.96 * sem, OR=np.exp(bm), chi2=statm, p=pm, n=len(m)))
    if has_msi:
        X = np.c_[one, *base, net]; b, se, stat, p = plrt(X, y, 2)
        rows.append(dict(set=tag, analysis="Firth: response ~ MSI/EBV + NetImmune", term="NetImmune per SD",
                         estimate=b, lo=b - 1.96 * se, hi=b + 1.96 * se, OR=np.exp(b), chi2=stat, p=p, n=len(m)))
    X = np.c_[one, net]; b, se, stat, p = plrt(X, y, 1)
    rows.append(dict(set=tag, analysis="Firth: response ~ NetImmune", term="NetImmune per SD",
                     estimate=b, lo=b - 1.96 * se, hi=b + 1.96 * se, OR=np.exp(b), chi2=stat, p=p, n=len(m)))


block(d, f"all QC-pass baseline (n = {len(d)})")
block(d[d.msi_ebv == 0].copy(), "MSS/EBV-negative (pre-registered subgroup)")
block(d, "sensitivity: MSI/EBV incl. E02 (TCGA-MSI)", msi_col="msi_ebv_tcga")
block(d[d.msi_ebv_tcga == 0].copy(), "sensitivity: MSS/EBV-negative excl. E02", msi_col="msi_ebv_tcga")
# sensitivity: primary model with library batch as a further covariate (scores are already z-scored within batch)
y = d.y.values.astype(float); one = np.ones(len(d))
X = np.c_[one, d.msi_ebv.values.astype(float), d.nextera.values.astype(float), z(np.log1p(d.cps.values)), z(d.NetImmune.values)]
b, se, stat, p = plrt(X, y, 4)
rows.append(dict(set="sensitivity: + library batch covariate", analysis="Firth: response ~ MSI/EBV + batch + CPS + score (primary: log(CPS+1))",
                 term="NetImmune per SD", estimate=b, lo=b - 1.96 * se, hi=b + 1.96 * se, OR=np.exp(b), chi2=stat, p=p, n=len(d)))
m = d[d.msi_ebv == 0]; y = m.y.values.astype(float)
X = np.c_[np.ones(len(m)), m.nextera.values.astype(float), z(np.log1p(m.cps.values)), z(m.NetImmune.values)]
b, se, stat, p = plrt(X, y, 3)
rows.append(dict(set="sensitivity: + library batch covariate, MSS/EBV-negative", analysis="Firth: response ~ batch + CPS + score (primary: log(CPS+1))",
                 term="NetImmune per SD", estimate=b, lo=b - 1.96 * se, hi=b + 1.96 * se, OR=np.exp(b), chi2=stat, p=p, n=len(m)))
tab = pd.crosstab([d.msi_ebv.rename("MSI/EBV"), (d.cps >= 1).rename("CPS>=1")], d.y.rename("responder"))
print(tab)
res = pd.DataFrame(rows)
res.to_csv(f"{O}/E22_firth_cps.csv", index=False)
tab.to_csv(f"{O}/E22_firth_cps_crosstab.csv")
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
sel = res[res.term.isin(["NetImmune per SD", "CPS term (NetImmune model)", "MSI/EBV term (NetImmune model)",
                         "NetImmune", "log(CPS+1)", "CPS>=1", "MSI/EBV", "CPS"])]
print(sel.round(4).to_string(index=False))
