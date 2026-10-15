"""E7c - Table 5 re-run on log2 scores with Firth-penalised logistic regression.

MSI/EBV status nearly separates responders in PRJEB25780, so ordinary maximum likelihood is
unstable. Firth (1993) regression maximises l(b) + 0.5 log|X'WX|; logistf/brglm2 are not
installed, so it is implemented here (Newton-Raphson on the modified score
U*(b) = X'(y - p + h(1/2 - p)), h = leverages) and checked on a 2x2 table, where the Firth
estimate of the log odds ratio equals the Haldane-corrected log((a+.5)(d+.5)/((b+.5)(c+.5))).
The penalised likelihood-ratio test for a coefficient follows logistf: the restricted fit
keeps the full design in the penalty and fixes that coefficient at 0.
Cross-validation mirrors revision/scripts/A3_cv_combo.py (apparent, LOOCV, 200 x stratified
5-fold, pooled out-of-fold AUC per repetition), with Firth fits in place of unpenalised ML.
"""
import numpy as np, pandas as pd
from scipy.stats import chi2, norm
from sklearn.model_selection import StratifiedKFold, LeaveOneOut
from sklearn.metrics import roc_auc_score

O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"


def firth(X, y, fixed=None, tol=1e-9, it=200):
    """X includes the intercept column. fixed: boolean mask of coefficients held at 0."""
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
        # step halving on the penalised log-likelihood
        f0, t = pll(b), 1.0
        while t > 1e-6 and pll(b + t * step) < f0 - 1e-12:
            t /= 2
        b = b + t * step
        if np.max(np.abs(t * step)) < tol:
            break
    p = 1 / (1 + np.exp(-X @ b)); W = p * (1 - p)
    cov = np.linalg.inv(X.T @ (X * W[:, None]))
    return b, cov, pll(b)


def predict(b, X):
    return 1 / (1 + np.exp(-X @ b))


# ---- check on a 2x2 table: Firth log-OR = Haldane-corrected log-OR
a_, b_, c_, d_ = 7, 0, 3, 10                     # exposed: 7 events, 0 non; unexposed: 3, 10
x = np.r_[np.ones(a_ + b_), np.zeros(c_ + d_)]
yy = np.r_[np.ones(a_), np.zeros(b_), np.ones(c_), np.zeros(d_)]
bt, _, _ = firth(np.c_[np.ones_like(x), x], yy)
haldane = np.log((a_ + .5) * (d_ + .5) / ((b_ + .5) * (c_ + .5)))
print(f"toy 2x2 with a zero cell: Firth log-OR {bt[1]:.6f} vs Haldane {haldane:.6f}")
assert abs(bt[1] - haldane) < 1e-6

# ---- data
d = pd.read_csv(f"{O}/E7_scores_log2.csv", index_col=0)
y = d.Responder.values.astype(float)
me = d.MSI_EBV.values.astype(float)
net = d.NetImmune.values; net = (net - net.mean()) / net.std()
MODELS = {"MSI/EBV only": np.c_[me], "NetImmune only": np.c_[net], "NetImmune + MSI/EBV": np.c_[net, me]}
add1 = lambda X: np.c_[np.ones(len(X)), X]

def fit_pred(Xtr, ytr, Xte):
    mu, sd = Xtr.mean(0), Xtr.std(0); sd[sd == 0] = 1           # standardise within fold, as A3
    b, _, _ = firth(add1((Xtr - mu) / sd), ytr)
    return predict(b, add1((Xte - mu) / sd))

rows = []
for name, X in MODELS.items():
    app = roc_auc_score(y, fit_pred(X, y, X))
    p = np.zeros(len(y))
    for tr, te in LeaveOneOut().split(X):
        p[te] = fit_pred(X[tr], y[tr], X[te])
    loo = roc_auc_score(y, p)
    cv = []
    for i in range(200):
        p = np.zeros(len(y))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=i).split(X, y):
            p[te] = fit_pred(X[tr], y[tr], X[te])
        cv.append(roc_auc_score(y, p))
    cv = np.array(cv)
    rows.append(dict(Model=name, Apparent=app, LOOCV=loo, CV5_mean=cv.mean(), CV5_lo=np.percentile(cv, 2.5),
                     CV5_hi=np.percentile(cv, 97.5), Optimism=app - cv.mean()))
r = pd.DataFrame(rows)
r.to_csv(f"{O}/E7_combo_firth_log2.csv", index=False)
print(r.round(3).to_string(index=False))

# ---- penalised LRT for NetImmune beyond MSI/EBV (and ML LRT for reference)
Xf = add1(np.c_[me, net])
bf, cf, lf = firth(Xf, y)
b0, _, l0 = firth(Xf, y, fixed=[False, False, True])
plr = 2 * (lf - l0); p_plr = chi2.sf(plr, 1)
se = np.sqrt(cf[2, 2]); z = bf[2] / se
import statsmodels.api as sm
m1 = sm.Logit(y, Xf).fit(disp=0); m0 = sm.Logit(y, Xf[:, :2]).fit(disp=0)
lr_ml = 2 * (m1.llf - m0.llf)
txt = (f"Firth: NetImmune per SD log-OR {bf[2]:.3f} (OR {np.exp(bf[2]):.2f}, Wald 95% CI "
       f"{np.exp(bf[2] - 1.96 * se):.2f}-{np.exp(bf[2] + 1.96 * se):.2f}); MSI/EBV log-OR {bf[1]:.3f} (OR {np.exp(bf[1]):.1f})\n"
       f"penalised LRT NetImmune beyond MSI/EBV: chi2 = {plr:.3f}, df 1, P = {p_plr:.4f}\n"
       f"(unpenalised ML LRT on log2 scores for reference: chi2 = {lr_ml:.3f}, P = {chi2.sf(lr_ml, 1):.4f}; "
       f"ML MSI/EBV coefficient {m1.params[1]:.2f})\n")
open(f"{O}/E7_lrt_firth_log2.txt", "w").write(txt)
print(txt)
