"""E44b addendum (post hoc, descriptive): leave-one-out, cohort weights, heterogeneity, main-text-rule counts and the power
of the gastric trials. Protocol: analysis/E44b_prereg.md (hash checked below). Reads E44 outputs only; nothing refitted.
Writes results/E44b_checks.csv and results/E44b.log."""
import hashlib, sys, datetime, pathlib
import numpy as np, pandas as pd
from scipy.stats import norm, chi2

A = pathlib.Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"); R = A / "results"
h = open(A / "E44b_prereg.sha256").read().split()[0]
assert hashlib.sha256(open(A / "E44b_prereg.md", "rb").read()).hexdigest() == h, "E44b protocol hash mismatch"
log = open(R / "E44b.log", "w")
def out(*a):
    s = " ".join(str(x) for x in a); print(s); log.write(s + "\n")
out("E44b protocol sha256 verified:", h[:12])

PC = pd.read_csv(R / "E44_B_per_cohort.csv"); PO = pd.read_csv(R / "E44_B_pooled.csv")
E32 = pd.read_csv(R / "E32_primary.csv"); E32 = E32[E32.set == "six (primary)"]; LIB = set(E32[E32.kind == "library"].signature)
COH = ["PRJEB25780", "PRJEB60680", "PERFECT"]
rows = []


def fe(y, se):
    w = 1 / se ** 2; est = (w * y).sum() / w.sum(); s = w.sum() ** -0.5
    return est, s, est - 1.96 * s, est + 1.96 * s, 2 * norm.sf(abs(est / s)), ((y - est) ** 2 * w).sum()


# 1. weights, shares, leave-one-out, heterogeneity
for nm in ["CD8A", "barriernet 0.3", "Immune", "TMEscore_z"]:
    d = PC[PC.name == nm].set_index("cohort").loc[COH]
    y, se = d.logOR.values, d.se.values; w = 1 / se ** 2
    est, s, lo, hi, p, Q = fe(y, se); I2 = max(0.0, (Q - 2) / Q)
    out(f"\n{nm}: pooled {est:+.3f} ({lo:+.3f} to {hi:+.3f}) P = {p:.3g}; Q = {Q:.2f} (P = {chi2.sf(Q, 2):.2f}), I2 = {100 * I2:.0f}%")
    rows.append(dict(item="heterogeneity", score=nm, cohort="all three", value=I2, est=est, lo=lo, hi=hi, p=p, Q=Q))
    for i, c in enumerate(COH):
        keep = [j for j in range(3) if j != i]
        e2, s2, lo2, hi2, p2, _ = fe(y[keep], se[keep])
        out(f"  {c}: weight {100 * w[i] / w.sum():.1f}%, share of estimate {100 * w[i] * y[i] / (w * y).sum():.0f}%; "
            f"without it {e2:+.3f} ({lo2:+.3f} to {hi2:+.3f}) P = {p2:.3g}")
        rows.append(dict(item="weight", score=nm, cohort=c, value=w[i] / w.sum()))
        rows.append(dict(item="share of pooled estimate", score=nm, cohort=c, value=w[i] * y[i] / (w * y).sum()))
        rows.append(dict(item="leave one out", score=nm, cohort=c, est=e2, lo=lo2, hi=hi2, p=p2))

# 2. significant under the main-text rule
for st in PO.set.unique():
    q = PO[PO.set == st]
    s52 = q[q.name.isin(LIB) & (q.p < 0.05)]; c52 = q[q.name.isin(LIB) & (q.cons_p < 0.05)]
    oth = q[~q.name.isin(LIB) & (q.cons_p < 0.05)]
    out(f"\n{st}: of 52, {len(s52)} at P < 0.05 by fixed effect, {len(c52)} by the main-text rule "
        f"({', '.join(f'{r.name} {r.cons_est:+.2f} P = {r.cons_p:.3f}' for r in c52.itertuples())}); "
        f"other scores/genes by the main-text rule: {', '.join(f'{r.name} P = {r.cons_p:.3f}' for r in oth.itertuples()) or 'none'}")
    rows.append(dict(item="signatures P<0.05, fixed effect", score="52 signatures", cohort=st, value=len(s52)))
    rows.append(dict(item="signatures P<0.05, main-text rule", score="; ".join(c52.name), cohort=st, value=len(c52)))
    rows.append(dict(item="other scores P<0.05, main-text rule", score="; ".join(oth.name), cohort=st, value=len(oth)))

# 3. power of the gastric trials for CD8A's six-comparison interaction
T = pd.read_csv(R / "E44_A_gastric_trials.csv").set_index("row")
d0 = float(E32[E32.signature == "CD8A"].est.iloc[0]); z = norm.ppf(0.975)
out(f"\npower for a per-SD interaction of {d0:+.4f} (CD8A, six comparisons):")
for k in ["CheckMate 649", "KEYNOTE-061", "gastric two, fixed effect"]:
    se = T.loc[k, "mdd80"] / 2.8016; pw = norm.cdf(d0 / se - z) + norm.cdf(-d0 / se - z)
    out(f"  {k}: SE {se:.3f}, power {100 * pw:.0f}%")
    rows.append(dict(item="power for CD8A interaction", score=k, cohort="", value=pw))

# 4. file dates
G = pathlib.Path("/groups/adv2105_gp/yichen/Yi/com/gastric")
utc = lambda p: datetime.datetime.fromtimestamp(p.stat().st_mtime, datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
out("\nfile dates:")
out(f"  barriernet_0.2.0.zip: {utc(G / 'barriernet_0.2.0.zip')}")
hits = [p for pat in ("*PRJEB25780*", "*Kim2018*") for p in G.rglob(pat) if "final_v" not in str(p)]
for p in sorted(hits, key=lambda p: p.stat().st_mtime)[:4]:
    out(f"  {p.relative_to(G)}: {utc(p)}")
rows.append(dict(item="file date", score="barriernet_0.2.0.zip", cohort=utc(G / "barriernet_0.2.0.zip")))
first = min(hits, key=lambda p: p.stat().st_mtime)
rows.append(dict(item="file date", score=f"earliest PRJEB25780/Kim2018 file ({first.name})", cohort=utc(first)))
pd.DataFrame(rows).to_csv(R / "E44b_checks.csv", index=False)
log.close()
