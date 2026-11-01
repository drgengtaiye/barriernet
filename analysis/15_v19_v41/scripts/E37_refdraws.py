"""E37 figure support: the random-reference draws of barriernet 0.3, recomputed exactly as E37_test.py (same seeds
[2037, comparison, size]) and saved for plotting (E37_test.py kept only their summaries). Checks that the six-comparison
summary equals E37_test_pooled.csv. Writes results/E37_reference_draws.csv (draw x set: pooled estimate)."""

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



P = pd.read_csv(f"{O}/E37_test_pooled.csv"); out = {}
SETS = [("six (primary)", ALL6), ("renal three", RENAL), ("breast three", BREAST), ("seven (with NeoPembrOv)", K7),
        ("four chemotherapy-backbone", BREAST + ["NeoPembrOV"])]
for lab, keys in SETS:
    D = np.stack([REF[k] for k in keys], 1)
    out[lab] = np.array([pool_rule(D[b, :, 0], D[b, :, 1])["est"] for b in range(B_REF)])
    r = P[(P.set == lab) & (P.score == "barriernet 0.3")].iloc[0]
    assert abs(np.nanmedian(out[lab]) - r.ref_median) < 1e-12 and abs(np.nanpercentile(out[lab], 95) - r.ref_p95) < 1e-12, lab
pd.DataFrame(out).to_csv(f"{O}/E37_reference_draws.csv", index=False)
print("reference draws reproduce E37_test_pooled.csv for", list(out))
