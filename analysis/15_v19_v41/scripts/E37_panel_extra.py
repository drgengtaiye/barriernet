"""E37 addition beyond the protocol (written after the E37 test results were seen): barriernet 0.3 against each of the
27 single genes of the two panels by the paired bootstrap, over the six comparisons and the three breast comparisons.

The protocol compared barriernet 0.3 only with CD8A and NetImmune 0.2.1; the manuscript's framework compares every
score with the whole single-gene panel, and in breast CD8A is weak (E33/E35), so this is reported as a declared
post-protocol addition whatever it shows. Setup and resamples exactly as E37_test.py (seeds [3032, comparison]).
"""
import hashlib, sys
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_hw = hashlib.sha256(open(f"{A}/results/E37_frozen_weights.tsv", "rb").read()).hexdigest()
assert _hw == open(f"{A}/E37_frozen.sha256").read().split()[0], "frozen weights hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E37x_step0_reproduction").replace("E34_coverage_check.csv", "E37x_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
NP = int(sys.argv[1]) if len(sys.argv) > 1 else 4
FW = pd.read_csv(f"{O}/E37_frozen_weights.tsv", sep="\t"); WG = dict(zip(FW.gene, FW.weight))
PANEL = T_PANEL + I_PANEL
SC = {}
for k in ALL6:
    Z = TR[k]["Z"]; nc = set(nonconst(Z))
    pres = [g for g in WG if g in nc and np.isfinite(Z.loc[g].to_numpy(float)).all()]
    SC[("BN03", k)] = wscore(Z, pres, np.array([WG[g] for g in pres]))
    for g in PANEL:
        if g in nc and np.isfinite(Z.loc[g].to_numpy(float)).all() and Z.loc[g].std() > 0:
            SC[(g, k)] = Z.loc[g].to_numpy(float)
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
with mp.get_context("fork").Pool(NP) as pool:
    res = pool.map(run, tasks, chunksize=2)
R = dict(zip(tasks, res))
prev = pd.read_csv(f"{O}/E37_test_per_comparison.csv")
chk = max(abs(R[("BN03", k)][0] - prev[(prev.comparison == k) & (prev.score == "BN03")].y.iloc[0]) for k in ALL6)
print(f"barriernet 0.3 per-comparison interactions reproduce E37_test.py: max |diff| {chk:.1e}"); assert chk < 1e-10
e35 = pd.read_csv(f"{O}/E35_by_gene.csv")


def paired(a, b, keys):
    ks = [k for k in keys if (a, k) in R and (b, k) in R]
    dd, vv = [], []
    for k in ks:
        ya, _, ba = R[(a, k)]; yb, _, bb = R[(b, k)]; db = ba - bb; ok = np.isfinite(db)
        dd.append(ya - yb); vv.append(np.var(db[ok], ddof=1))
    dd, vv = np.array(dd), np.array(vv); w = 1 / vv; D_ = np.sum(w * dd) / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(k=len(ks), diff=D_, se=se, lo=D_ - 1.96 * se, hi=D_ + 1.96 * se, p=2 * norm.sf(abs(D_ / se)))

rows = []
for lab, keys, kmin in [("six (primary)", ALL6, 5), ("breast three", BREAST, 3)]:
    for g in PANEL:
        r = paired("BN03", g, keys, )
        if r["k"] >= kmin:
            q = pool_rule([R[(g, k)][0] for k in keys if (g, k) in R], [R[(g, k)][1] for k in keys if (g, k) in R])
            rows.append(dict(set=lab, gene=g, panel="T cell" if g in T_PANEL else "IFN-γ/checkpoint", gene_est=q["est"], **r))
P = pd.DataFrame(rows); P.to_csv(f"{O}/E37_panel_extra.csv", index=False)
for lab in ("six (primary)", "breast three"):
    s = P[P.set == lab]
    print(f"\n{lab}: barriernet 0.3 significantly above {int(((s['diff'] > 0) & (s.p < 0.05)).sum())} and below "
          f"{int(((s['diff'] < 0) & (s.p < 0.05)).sum())} of {len(s)} genes")
    print(s.sort_values("gene_est", ascending=False)[["gene", "gene_est", "diff", "lo", "hi", "p"]].round(3).to_string(index=False))
