"""E37 step 4 - test the frozen barriernet 0.3 in the randomized comparisons.

Protocol: analysis/E37_prereg.md (hash checked); frozen weights results/E37_frozen_weights.tsv (hash in E37_frozen.sha256,
logged before this script was written). Setup: E34_trials.py up to its fitting stage (E32 setup + NeoPembrOv), executed
unchanged except that its step-0 and coverage tables are written under E37 names. Bootstrap resamples as E32/E33
(seeds [3032, comparison]) and E34/E35 (NeoPembrOv, [3034, 7]); random reference seeds [2037, comparison, size].
Step 0: the NetImmune 0.2.1 paired difference from CD8A over the six must reproduce E32 (|diff| < 1e-8).
"""
import hashlib, zlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h37 = hashlib.sha256(open(f"{A}/E37_prereg.md", "rb").read()).hexdigest()
assert _h37 == open(f"{A}/E37_prereg.sha256").read().split()[0], "E37 protocol hash mismatch"
_hw = hashlib.sha256(open(f"{A}/results/E37_frozen_weights.tsv", "rb").read()).hexdigest()
assert _hw == open(f"{A}/E37_frozen.sha256").read().split()[0], "frozen weights hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E37_step0_reproduction").replace("E34_coverage_check.csv", "E37_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
print("E37 protocol sha256 verified:", _h37[:12], "| frozen weights sha256 verified:", _hw[:12])

FW = pd.read_csv(f"{O}/E37_frozen_weights.tsv", sep="\t")
WG = dict(zip(FW.gene, FW.weight))
K7 = ALL6 + ["NeoPembrOV"]


def bn03(k):
    Z = TR[k]["Z"]; nc = set(nonconst(Z))
    pres = [g for g in WG if g in nc and np.isfinite(Z.loc[g].to_numpy(float)).all()]
    if len(pres) < 3 or len(pres) < 0.5 * len(WG):
        return None, pres
    return wscore(Z, pres, np.array([WG[g] for g in pres])), pres


SC, PRES = {}, {}
for k in K7:
    s, p = bn03(k); PRES[k] = p
    if s is not None:
        SC[("BN03", k)] = s
    SC[("CD8A", k)] = OBS[("CD8A", k)][0] if k in ALL6 else OB[("gene", "CD8A")]
    SC[("NetImmune", k)] = OBS[("NetImmune", k)][0] if k in ALL6 else OB[("sig", "NetImmune")]
print("barriernet 0.3 genes present:", {k: len(PRES[k]) for k in K7})

BIDX = {}
for k in ALL6:
    rng_ = np.random.default_rng([3032, KIDX[k]]); a = TR[k]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
    BIDX[k] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]
rng_ = np.random.default_rng([3034, 7]); a = TR["NeoPembrOV"]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
BIDX["NeoPembrOV"] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]


def run(task):
    name, k = task
    sc = SC[(name, k)]; d = TR[k]
    y, se = m_int(d, sc); bs = []
    for idx in BIDX[k]:
        try:
            bs.append(m_int(take(d, idx), sc[idx])[0])
        except Exception:
            bs.append(np.nan)
    return y, se, np.array(bs)

tasks = [t for t in SC]
with mp.get_context("fork").Pool(NPROC) as pool:
    res = pool.map(run, tasks, chunksize=1)
R = dict(zip(tasks, res))


def paired(a, b, keys):
    dd, vv = [], []
    for k in keys:
        ya, _, ba = R[(a, k)]; yb, _, bb = R[(b, k)]; db = ba - bb; ok = np.isfinite(db)
        dd.append(ya - yb); vv.append(np.var(db[ok], ddof=1))
    dd, vv = np.array(dd), np.array(vv); w = 1 / vv; D_ = np.sum(w * dd) / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(diff=D_, se=se, lo=D_ - 1.96 * se, hi=D_ + 1.96 * se, p=2 * norm.sf(abs(D_ / se)), k=len(keys))

s0 = paired("NetImmune", "CD8A", ALL6)
e32 = pd.read_csv(f"{O}/E32_primary.csv").set_index("signature").loc["NetImmune"]
print(f"step 0: NetImmune - CD8A {s0['diff']:+.6f} (E32 {e32.diff_vs_CD8A:+.6f})"); assert abs(s0["diff"] - e32.diff_vs_CD8A) < 1e-8

# ------------------------------------------------------------------ random reference (Q4)
UNIV = {}
for k in K7:
    Z = TR[k]["Z"]; g = [x for x in nonconst(Z) if np.isfinite(Z.loc[x].to_numpy(float)).all()]
    UNIV[k] = Z.loc[g].to_numpy(float)
B_REF = 1000


def ref_task(args):
    k, b0, b1 = args
    w = np.array([WG[g] for g in PRES[k]]); n = len(w); U = UNIV[k]; d = TR[k]
    rng = np.random.default_rng([2037, K7.index(k), n])
    out = []
    for b in range(b1):
        idx = rng.choice(U.shape[0], n, replace=False); ww = w[rng.permutation(n)]
        if b >= b0:
            try:
                out.append(m_int(d, U[idx].T @ ww / np.abs(ww).sum()))
            except Exception:
                out.append((np.nan, np.nan))
    return out

rt = [(k, b0, b0 + 250) for k in K7 if ("BN03", k) in SC for b0 in range(0, B_REF, 250)]
with mp.get_context("fork").Pool(NPROC) as pool:
    rr = pool.map(ref_task, rt, chunksize=1)
REF = {}
for (k, b0, b1), r in zip(rt, rr):
    REF.setdefault(k, [None] * (B_REF // 250))[b0 // 250] = r
REF = {k: np.array(sum(v, [])) for k, v in REF.items()}


def pooled(name, keys):
    ys = [R[(name, k)][0] for k in keys]; ses = [R[(name, k)][1] for k in keys]
    return pool_rule(ys, ses)


def refcmp(keys):
    D = np.stack([REF[k] for k in keys], 1)
    est_b = np.array([pool_rule(D[b, :, 0], D[b, :, 1])["est"] for b in range(B_REF)])
    obs = pooled("BN03", keys)["est"]
    return dict(ref_median=np.nanmedian(est_b), ref_p95=np.nanpercentile(est_b, 95), emp_p=(1 + np.sum(est_b >= obs)) / (B_REF + 1),
                exceeds=obs > np.nanpercentile(est_b, 95))

SETS = [("six (primary)", ALL6), ("renal three", RENAL), ("breast three", BREAST), ("NeoPembrOv", ["NeoPembrOV"]),
        ("seven (with NeoPembrOv)", K7), ("four chemotherapy-backbone", BREAST + ["NeoPembrOV"])]
ROWS = []
for lab, keys in SETS:
    for nm in ("BN03", "NetImmune", "CD8A"):
        if len(keys) == 1:
            y, se, _ = R[(nm, keys[0])]; q = dict(est=y, lo=y - 1.96 * se, hi=y + 1.96 * se, p=2 * norm.sf(abs(y / se)), k=1)
        else:
            q = pooled(nm, keys)
        row = dict(set=lab, score={"BN03": "barriernet 0.3", "NetImmune": "NetImmune 0.2.1", "CD8A": "CD8A"}[nm], **q)
        if nm == "BN03":
            for ref in ("CD8A", "NetImmune"):
                pp = paired("BN03", ref, keys)
                row.update({f"diff_vs_{ref}": pp["diff"], f"lo_vs_{ref}": pp["lo"], f"hi_vs_{ref}": pp["hi"], f"p_vs_{ref}": pp["p"]})
            if len(keys) > 1:
                row.update(refcmp(keys))
        ROWS.append(row)
RS = pd.DataFrame(ROWS); RS.to_csv(f"{O}/E37_test_pooled.csv", index=False)
PC = pd.DataFrame([dict(comparison=k, score=nm, y=R[(nm, k)][0], se=R[(nm, k)][1], genes=len(PRES[k]) if nm == "BN03" else None)
                   for k in K7 for nm in ("BN03", "NetImmune", "CD8A")])
PC.to_csv(f"{O}/E37_test_per_comparison.csv", index=False)

six = RS[(RS.set == "six (primary)") & (RS.score == "barriernet 0.3")].iloc[0]
EXP = pd.DataFrame([
    dict(q="Q1 (primary)", expectation="pooled six-comparison interaction > 0 with P < 0.05",
         result=f"{six.est:+.3f} ({six.lo:+.3f} to {six.hi:+.3f}); P = {six.p:.4f}", held=bool(six.est > 0 and six.p < 0.05)),
    dict(q="Q2", expectation="paired difference from CD8A not significant",
         result=f"{six.diff_vs_CD8A:+.3f} ({six.lo_vs_CD8A:+.3f} to {six.hi_vs_CD8A:+.3f}); P = {six.p_vs_CD8A:.3f}", held=bool(six.p_vs_CD8A >= 0.05)),
    dict(q="Q3", expectation="paired difference from NetImmune 0.2.1 not significant",
         result=f"{six.diff_vs_NetImmune:+.3f} ({six.lo_vs_NetImmune:+.3f} to {six.hi_vs_NetImmune:+.3f}); P = {six.p_vs_NetImmune:.3f}",
         held=bool(six.p_vs_NetImmune >= 0.05)),
    dict(q="Q4", expectation="above the 95th percentile of 1,000 size-matched random gene sets",
         result=f"reference median {six.ref_median:+.3f}, 95th {six.ref_p95:+.3f}; empirical P = {six.emp_p:.4f}", held=bool(six.exceeds)),
    dict(q="Q5", expectation="descriptive", result="see E37_test_pooled.csv", held=None)])
EXP.to_csv(f"{O}/E37_test_expectations.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print("\n== expectations\n", EXP.to_string(index=False))
print("\n== pooled\n", RS.round(4).to_string(index=False))
print("\n== per comparison\n", PC.round(3).to_string(index=False))
