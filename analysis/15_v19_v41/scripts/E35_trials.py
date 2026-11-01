"""E35 - single genes as the reference in the breast comparisons (and with NeoPembrOV).

Protocol: analysis/E35_prereg.md (sha256 in E35_prereg.sha256, logged 2026-10-27T03:16:39Z, before this script ran).
Setup: E34_trials.py up to its fitting stage (E32 setup plus NeoPembrOV), executed unchanged except that its step-0
and coverage tables are written under E35 names. Paired differences as E33 Q3 with the E32/E33 seeds (NeoPembrOV: E34).
"""
import hashlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h35 = hashlib.sha256(open(f"{A}/E35_prereg.md", "rb").read()).hexdigest()
assert _h35 == open(f"{A}/E35_prereg.sha256").read().split()[0], "E35 protocol hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E35_step0_reproduction").replace("E34_coverage_check.csv", "E35_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
print("E35 protocol sha256 verified:", _h35[:12])

PANEL = T_PANEL + I_PANEL
KEYS4 = BREAST + ["NeoPembrOV"]


def score(kind, name, k):
    if k == "NeoPembrOV":
        return OB.get((kind, name))
    if kind == "sig":
        return OBS[(name, k)][0] if (name, k) in OBS else None
    Z = TR[k]["Z"]
    return Z.loc[name].to_numpy(float) if name in Z.index and np.isfinite(Z.loc[name].to_numpy(float)).all() and Z.loc[name].std() > 0 else None


ITEMS = [("sig", s) for s in LIB] + [("gene", g) for g in PANEL]
BIDX = {}
for k in BREAST:
    rng_ = np.random.default_rng([3032, KIDX[k]]); a = TR[k]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
    BIDX[k] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]
rng_ = np.random.default_rng([3034, 7]); a = TR["NeoPembrOV"]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
BIDX["NeoPembrOV"] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]


def run(task):
    kind, name, k = task
    sc = score(kind, name, k); d = TR[k]
    if sc is None:
        return None
    try:
        y, se = m_int(d, sc)
    except Exception:
        y, se = np.nan, np.nan
    bs = []
    for idx in BIDX[k]:
        try:
            bs.append(m_int(take(d, idx), sc[idx])[0])
        except Exception:
            bs.append(np.nan)
    return y, se, np.array(bs)

tasks = [(kd, nm, k) for kd, nm in ITEMS for k in KEYS4]
with mp.get_context("fork").Pool(NPROC) as pool:
    res = pool.map(run, tasks, chunksize=2)
R = {t: r for t, r in zip(tasks, res) if r is not None}
print(f"fits with bootstrap: {len(R)} of {len(tasks)}")


def paired(sig, gene, keys, kmin):
    ks = [k for k in keys if ("sig", sig, k) in R and ("gene", gene, k) in R]
    if len(ks) < kmin:
        return None
    dd, vv = [], []
    for k in ks:
        ys, _, bs = R[("sig", sig, k)]; yg, _, bg = R[("gene", gene, k)]
        db = bs - bg; ok = np.isfinite(db); dd.append(ys - yg); vv.append(np.var(db[ok], ddof=1))
    dd, vv = np.array(dd), np.array(vv); w = 1 / vv
    D_ = np.sum(w * dd) / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(k=len(ks), diff=D_, se=se, lo=D_ - 1.96 * se, hi=D_ + 1.96 * se, p=2 * norm.sf(abs(D_ / se)))


rows = []
for setname, keys, kmin in [("breast three", BREAST, 3), ("four chemotherapy-backbone", KEYS4, 4)]:
    for g in PANEL:
        for s in LIB:
            r = paired(s, g, keys, kmin)
            if r:
                rows.append(dict(set=setname, signature=s, reference=g, **r))
PP = pd.DataFrame(rows); PP.to_csv(f"{O}/E35_paired.csv", index=False)
E33 = pd.read_csv(f"{O}/E33_paired.csv"); e = E33[(E33.set == "breast three") & (E33.reference == "CD8A")].set_index("signature")
mine = PP[(PP.set == "breast three") & (PP.reference == "CD8A")].set_index("signature")
dmax = float(np.nanmax(np.abs(mine["diff"] - e.loc[mine.index, "diff"])))
print(f"step 1: breast differences from CD8A reproduce E33 (max |diff| {dmax:.2e})"); assert dmax < 1e-8

# per-gene summary, with each gene's pooled breast interaction
G = []
for setname, keys in [("breast three", BREAST), ("four chemotherapy-backbone", KEYS4)]:
    for g in PANEL:
        ys = [(R[("gene", g, k)][0], R[("gene", g, k)][1]) for k in keys if ("gene", g, k) in R]
        if len(ys) < len(keys):
            continue
        q = pool_rule([y for y, _ in ys], [s for _, s in ys])
        sub = PP[(PP.set == setname) & (PP.reference == g)]
        G.append(dict(set=setname, gene=g, panel="T cell" if g in T_PANEL else "IFN-γ/checkpoint", est=q["est"], lo=q["lo"],
                      hi=q["hi"], p=q["p"], n_sig=len(sub), exceeded_by=int(((sub["diff"] > 0) & (sub.p < 0.05)).sum()),
                      exceeding=", ".join(sub[(sub["diff"] > 0) & (sub.p < 0.05)].signature)))
G = pd.DataFrame(G); G.to_csv(f"{O}/E35_by_gene.csv", index=False)
b = G[G.set == "breast three"].set_index("gene"); f4 = G[G.set == "four chemotherapy-backbone"].set_index("gene")
q1 = {g: int(b.loc[g, "exceeded_by"]) for g in ("GZMB", "IDO1", "CXCL10")}
q2 = int((b.loc[[g for g in T_PANEL if g in b.index], "exceeded_by"] > 2).sum())
q3 = b.exceeded_by.idxmin(), int(b.exceeded_by.min())
q4 = {g: int(f4.loc[g, "exceeded_by"]) for g in ("GZMB", "IDO1", "CXCL10")}
EXP = pd.DataFrame([
    dict(q="Q1 (primary)", expectation="breast: at most 2 of 52 above each of GZMB, IDO1, CXCL10", result=str(q1), held=max(q1.values()) <= 2),
    dict(q="Q2", expectation="breast: at least 5 of 13 T-cell genes exceeded by > 2 signatures", result=f"{q2}/13", held=q2 >= 5),
    dict(q="Q3", expectation="descriptive", result=f"fewest: {q3[0]} ({q3[1]})", held=None),
    dict(q="Q4", expectation="four chemotherapy-backbone: at most 2 of 52 above each of GZMB, IDO1, CXCL10", result=str(q4), held=max(q4.values()) <= 2)])
EXP.to_csv(f"{O}/E35_expectations.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 160)
print("\n== expectations\n", EXP.to_string(index=False))
print("\n== breast, by gene (sorted by interaction)\n", b.sort_values("est", ascending=False)[["panel", "est", "p", "exceeded_by", "exceeding"]].round(3).to_string())
print("\n== four with NeoPembrOV, by gene\n", f4.sort_values("est", ascending=False)[["est", "p", "exceeded_by"]].round(3).to_string())
