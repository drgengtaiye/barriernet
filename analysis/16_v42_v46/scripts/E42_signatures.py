"""E42 - barriernet 0.3 against each of the 52 library signatures by the paired bootstrap (six randomized comparisons).

Protocol: analysis/E42_prereg.md (hash checked); frozen weights results/E37_frozen_weights.tsv (hash checked). Setup:
E34_trials.py up to its fitting stage, executed unchanged as in E37_test.py except that its step-0 and coverage tables
are written under E42 names. Bootstrap resamples as E32/E37 (seeds [3032, comparison]; NeoPembrOv [3034, 7]).
Step 0: the unoriented paired differences of NetImmune 0.2.1 and the 52 signatures from CD8A reproduce E32_primary.csv,
and barriernet 0.3 - CD8A reproduces E37_test_pooled.csv (|diff| < 1e-8).
Orientations: conservative (sign of the E32 six-comparison estimate; primary), declared direction, unoriented.
CD8A against each signature is computed the same way for reference.
"""
import hashlib, sys
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h42 = hashlib.sha256(open(f"{A}/E42_prereg.md", "rb").read()).hexdigest()
assert _h42 == open(f"{A}/E42_prereg.sha256").read().split()[0], "E42 protocol hash mismatch"
_hw = hashlib.sha256(open(f"{A}/results/E37_frozen_weights.tsv", "rb").read()).hexdigest()
assert _hw == open(f"{A}/E37_frozen.sha256").read().split()[0], "frozen weights hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E42_step0_reproduction").replace("E34_coverage_check.csv", "E42_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
print("E42 protocol sha256 verified:", _h42[:12], "| frozen weights sha256 verified:", _hw[:12])

FW = pd.read_csv(f"{O}/E37_frozen_weights.tsv", sep="\t"); WG = dict(zip(FW.gene, FW.weight))
K7 = ALL6 + ["NeoPembrOV"]
SC = {}
for k in K7:
    Z = TR[k]["Z"]; nc = set(nonconst(Z))
    pres = [g for g in WG if g in nc and np.isfinite(Z.loc[g].to_numpy(float)).all()]
    if len(pres) >= 3 and len(pres) >= 0.5 * len(WG):
        SC[("BN03", k)] = wscore(Z, pres, np.array([WG[g] for g in pres]))
    for s in list(LIB) + ["CD8A", "NetImmune"]:
        if k in ALL6:
            if (s, k) in OBS:
                SC[(s, k)] = OBS[(s, k)][0]
        else:
            key = ("gene", "CD8A") if s == "CD8A" else ("sig", s)
            if key in OB:
                SC[(s, k)] = OB[key]
print("barriernet 0.3 scored in:", [k for k in K7 if ("BN03", k) in SC])
print("signatures scored per comparison:", {k: sum(1 for s in LIB if (s, k) in SC) for k in K7})

BIDX = {}
for k in ALL6:
    rng_ = np.random.default_rng([3032, KIDX[k]]); a = TR[k]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
    BIDX[k] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]
rng_ = np.random.default_rng([3034, 7]); a = TR["NeoPembrOV"]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
BIDX["NeoPembrOV"] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]


def run(task):
    name, k = task; sc = SC[(name, k)]; d = TR[k]; y, se = m_int(d, sc); bs = []
    for idx in BIDX[k]:
        try:
            bs.append(m_int(take(d, idx), sc[idx])[0])
        except Exception:
            bs.append(np.nan)
    return y, se, np.array(bs)

NP = int(sys.argv[1]) if len(sys.argv) > 1 else NPROC
tasks = list(SC)
with mp.get_context("fork").Pool(NP) as pool:
    res = pool.map(run, tasks, chunksize=1)
R = dict(zip(tasks, res))


def paired(a, b, keys, ob=1.0, kmin=4):
    """a - ob * b, fixed-effect inverse variance of the bootstrap variances over the comparisons scoring both."""
    ks = [k for k in keys if (a, k) in R and (b, k) in R]
    if len(ks) < kmin:
        return None
    dd, vv = [], []
    for k in ks:
        ya, _, ba = R[(a, k)]; yb, _, bb = R[(b, k)]; db = ba - ob * bb; ok = np.isfinite(db)
        dd.append(ya - ob * yb); vv.append(np.var(db[ok], ddof=1))
    dd, vv = np.array(dd), np.array(vv); w = 1 / vv; D_ = np.sum(w * dd) / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(k=len(ks), diff=D_, se=se, lo=D_ - 1.96 * se, hi=D_ + 1.96 * se, p=2 * norm.sf(abs(D_ / se)), mdd=2.80 * se)

# ------------------------------------------------------------------ step 0
E32P = pd.read_csv(f"{O}/E32_primary.csv"); E32P = E32P[E32P.set == "six (primary)"].set_index("signature")
gap = max(abs(paired(s, "CD8A", ALL6)["diff"] - E32P.loc[s, "diff_vs_CD8A"]) for s in list(LIB) + ["NetImmune"])
e37 = pd.read_csv(f"{O}/E37_test_pooled.csv"); e37 = e37[(e37.set == "six (primary)") & (e37.score == "barriernet 0.3")].iloc[0]
g37 = abs(paired("BN03", "CD8A", ALL6)["diff"] - e37.diff_vs_CD8A)
print(f"step 0: 52 signatures + NetImmune vs CD8A reproduce E32 (max |diff| {gap:.1e}); barriernet 0.3 - CD8A reproduces "
      f"E37 ({g37:.1e})"); assert gap < 1e-8 and g37 < 1e-8

# ------------------------------------------------------------------ paired differences
ORI = {"conservative": {s: (1.0 if E32P.loc[s, "est"] >= 0 else -1.0) for s in LIB},
       "declared": {s: float(LIB[s]["dir"]) for s in LIB}, "unoriented": {s: 1.0 for s in LIB}}
rows = []
for lab, keys in [("six (primary)", ALL6), ("seven (with NeoPembrOv)", K7)]:
    for ori, o in ORI.items():
        if lab.startswith("seven") and ori != "conservative":
            continue
        for s in LIB:
            r = paired("BN03", s, keys, o[s]); c = paired("CD8A", s, keys, o[s])
            if r is None:
                continue
            rows.append(dict(set=lab, orientation=ori, signature=s, family=LIB[s]["family"], declared_dir=LIB[s]["dir"],
                             sign=o[s], sig_est_E32=o[s] * E32P.loc[s, "est"], **r,
                             **({f"cd8a_{x}": v for x, v in c.items() if x != "k"} if c else {})))
P = pd.DataFrame(rows)
P["above"] = (P["diff"] > 0) & (P.p < 0.05); P["below"] = (P["diff"] < 0) & (P.p < 0.05)
P["within010"] = (P.lo > -0.10) & (P.hi < 0.10)
P["cd8a_above"] = (P.cd8a_diff > 0) & (P.cd8a_p < 0.05); P["cd8a_below"] = (P.cd8a_diff < 0) & (P.cd8a_p < 0.05)
P["cd8a_within010"] = (P.cd8a_lo > -0.10) & (P.cd8a_hi < 0.10)
P.to_csv(f"{O}/E42_paired.csv", index=False)
CNT = (P.groupby(["set", "orientation"], sort=False)
        .agg(n=("signature", "size"), above=("above", "sum"), below=("below", "sum"), within010=("within010", "sum"),
             median_mdd=("mdd", "median"), cd8a_above=("cd8a_above", "sum"), cd8a_below=("cd8a_below", "sum"),
             cd8a_within010=("cd8a_within010", "sum"), cd8a_median_mdd=("cd8a_mdd", "median")).reset_index())
fam = (P[P.set == "six (primary)"].groupby(["orientation", "family"]).agg(n=("signature", "size"), above=("above", "sum"),
       below=("below", "sum"), cd8a_above=("cd8a_above", "sum")).reset_index())
CNT.to_csv(f"{O}/E42_counts.csv", index=False); fam.to_csv(f"{O}/E42_counts_by_family.csv", index=False)

bn = pool_rule([R[("BN03", k)][0] for k in ALL6], [R[("BN03", k)][1] for k in ALL6])
cons = P[(P.set == "six (primary)") & (P.orientation == "conservative")]
rank = 1 + int((cons.sig_est_E32 > bn["est"]).sum())
cd8 = paired("BN03", "CD8A", ALL6); nim = paired("BN03", "NetImmune", ALL6)
c6 = CNT[(CNT.set == "six (primary)")].set_index("orientation")
EXP = pd.DataFrame([
    dict(q="A1 (primary)", expectation="conservative orientation: significantly above >= 5 of the 52, below none",
         result=f"above {int(c6.loc['conservative', 'above'])}, below {int(c6.loc['conservative', 'below'])} of "
                f"{int(c6.loc['conservative', 'n'])}",
         held=bool(c6.loc["conservative", "above"] >= 5 and c6.loc["conservative", "below"] == 0)),
    dict(q="A2", expectation="declared direction: significantly above >= 10 of the 52, below none",
         result=f"above {int(c6.loc['declared', 'above'])}, below {int(c6.loc['declared', 'below'])} of "
                f"{int(c6.loc['declared', 'n'])}",
         held=bool(c6.loc["declared", "above"] >= 10 and c6.loc["declared", "below"] == 0))])
EXP.to_csv(f"{O}/E42_expectations.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print(f"\nbarriernet 0.3 six-comparison pooled {bn['est']:+.3f} ({bn['lo']:+.3f} to {bn['hi']:+.3f}); rank {rank} of "
      f"{len(cons) + 1} with every signature in its conservative orientation")
print(f"barriernet 0.3 - CD8A {cd8['diff']:+.3f} ({cd8['lo']:+.3f} to {cd8['hi']:+.3f}), MDD {cd8['mdd']:.3f}; "
      f"- NetImmune {nim['diff']:+.3f} ({nim['lo']:+.3f} to {nim['hi']:+.3f})")
print("\n== expectations\n", EXP.to_string(index=False))
print("\n== counts\n", CNT.round(3).to_string(index=False))
print("\n== by family (six)\n", fam.to_string(index=False))
print("\n== six, conservative orientation\n", cons.sort_values("diff")[["signature", "family", "sign", "sig_est_E32", "k", "diff",
      "lo", "hi", "p", "mdd", "cd8a_diff", "cd8a_p"]].round(3).to_string(index=False))
