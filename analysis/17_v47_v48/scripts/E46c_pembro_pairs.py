"""E46c - addendum to E46: near-identical samples inside ISPY2_pembro. Protocol: analysis/E46c_addendum.md (hash checked)."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "4"
import hashlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h = hashlib.sha256(open(f"{A}/E46c_addendum.md", "rb").read()).hexdigest(); assert _h == open(f"{A}/E46c_addendum.sha256").read().split()[0]
_s = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0].replace("E34_step0_reproduction", "E46c_step0_reproduction").replace("E34_coverage_check.csv", "E46c_coverage_check.csv")
exec(compile(_s, f"{A}/E34_trials.py", "exec"))
d = TR["ISPY2_pembro"]; X = d["X"]; X = X[~X.index.duplicated()]; ids = X.columns.astype(str); arm = np.asarray(d["arm"])
M = np.nan_to_num(X.to_numpy(float)); M = M - M.mean(1, keepdims=True); Z_ = M / np.maximum(np.sqrt((M ** 2).sum(0, keepdims=True)), 1e-12); C = Z_.T @ Z_
mean_prof = M.mean(1); iu = np.triu_indices(C.shape[0], 1); order = np.argsort(-C[iu]); rows = []
adj = globals().get("ADJ_P", None)
print("pairwise r: median", round(float(np.median(C[iu])), 3), "99.9th pct", round(float(np.percentile(C[iu], 99.9)), 3), "max", round(float(C[iu].max()), 3))
nearest = (C - 2 * np.eye(C.shape[0])).max(1); print("best partner r per sample: median", round(float(np.median(nearest)), 3), "n > 0.8:", int((nearest > 0.8).sum()), "of", len(nearest))
for o in order:
    a_, b_ = iu[0][o], iu[1][o]
    if C[a_, b_] <= 0.8:
        break
    rows.append(dict(sample_a=ids[a_], sample_b=ids[b_], r=C[a_, b_], arm_a=int(arm[a_]), arm_b=int(arm[b_]), proc_a=None if adj is None else str(adj[a_]), proc_b=None if adj is None else str(adj[b_]), sd_a=float(M[:, a_].std()), sd_b=float(M[:, b_].std())))
P = pd.DataFrame(rows); P.to_csv(f"{O}/E46c_pembro_pairs.csv", index=False); print(P.round(3).to_string(index=False)); print("median sample SD of gene-centered profile:", round(float(np.median(M.std(0))), 3))
if len(P):
    drop = set(P.sample_b); keep = np.array([i not in drop for i in ids]); dq = take(d, np.where(keep)[0]); Z = d["Z"]; SENS = []; print("removed", sorted(drop))
    for nm in ["CD8A", "R3", "Tcell-inflamed_Ayers"]:
        ys, ss = [], []
        for k in ALL6:
            Zk = TR[k]["Z"]
            sc = (wscore(Zk, ["CD8A", "PRF1", "GZMA"], np.ones(3)) if nm == "R3" else (Zk.loc["CD8A"].to_numpy(float) if nm == "CD8A" else OBS[(nm, k)][0]))
            y, se = m_int(TR[k], sc) if k != "ISPY2_pembro" else m_int(dq, sc[keep]); ys.append(y); ss.append(se)
        r = pool_rule(ys, ss); SENS.append(dict(score=nm, removed=len(drop), est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"]))
    SENS = pd.DataFrame(SENS); SENS.to_csv(f"{O}/E46c_sensitivity.csv", index=False); print(SENS.round(3).to_string(index=False))
