"""E46b - addendum to E46 part D: negative control, gene-centered profiles and positive control for the I-SPY2 flag. Reads no outcome."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "4"
import hashlib, sys
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h = hashlib.sha256(open(f"{A}/E46b_addendum.md", "rb").read()).hexdigest(); assert _h == open(f"{A}/E46b_addendum.sha256").read().split()[0]
_s = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0].replace("E34_step0_reproduction", "E46b_step0_reproduction").replace("E34_coverage_check.csv", "E46b_coverage_check.csv")
exec(compile(_s, f"{A}/E34_trials.py", "exec"))
L = open(f"{O}/E46b_ispy2.log", "w")
def say(*a):
    t = " ".join(str(x) for x in a); print(t); L.write(t + "\n"); L.flush()
F = pd.read_csv(f"{O}/E46_D_ispy2_flags.csv", dtype={"sample": str, "partner": str})
Xd, Xp = TR["ISPY2_durva"]["X"], TR["ISPY2_pembro"]["X"]; Xd = Xd[~Xd.index.duplicated()]; Xp = Xp[~Xp.index.duplicated()]; g = Xd.index.intersection(Xp.index)
ad, ap = np.asarray(TR["ISPY2_durva"]["arm"]), np.asarray(TR["ISPY2_pembro"]["arm"]); idd, idp = Xd.columns.astype(str), Xp.columns.astype(str)
say("shared genes", len(g), "| durva controls", int((ad == 0).sum()), "treated", int((ad == 1).sum()), "| pembro controls", int((ap == 0).sum()), "treated", int((ap == 1).sum()))
fa = F[(F.set.str.startswith("all")) & (F.direction == "durva to pembro")]; ctrl = set(idd[ad == 0]); tr = set(idd[ad == 1])
fc, ft = fa[fa["sample"].isin(ctrl)], fa[fa["sample"].isin(tr)]
ROW = [dict(check="1 flag rate, controls (against all pembro)", value=len(fc), of=len(ctrl)), dict(check="1 flag rate, treated (cannot be shared)", value=len(ft), of=len(tr)), dict(check="1 distinct partners of flagged", value=fa.partner.nunique(), of=len(fa)),
       dict(check="1 most frequent partner, flags", value=int(fa.partner.value_counts().iloc[0]), of=len(fa))]
def z(M):
    M = np.nan_to_num(M); M = M - M.mean(1, keepdims=True); return M / np.maximum(np.sqrt((M ** 2).sum(0, keepdims=True)), 1e-12)
Za, Zb = z(Xd.loc[g].to_numpy(float)), z(Xp.loc[g].to_numpy(float)); C = Za.T @ Zb
for lab, ia in [("controls", ad == 0), ("treated", ad == 1)]:
    best = C[np.ix_(ia, ap == 0)].max(1); bestall = C[ia].max(1)
    ROW += [dict(check=f"2 best gene-centered r, durva {lab} vs pembro controls: median", value=float(np.median(best))), dict(check=f"2 best gene-centered r, durva {lab} vs pembro controls: max", value=float(best.max())), dict(check=f"2 durva {lab} above 0.8", value=int((best > 0.8).sum()), of=int(ia.sum()))]
CC = C[np.ix_(ad == 0, ap == 0)]; a_ = CC.argmax(1); b_ = CC.argmax(0); mut = sum(b_[a_[i]] == i for i in range(CC.shape[0])); ROW.append(dict(check="1 mutual best pairs among durva controls", value=int(mut), of=CC.shape[0]))
# positive control
rng = np.random.default_rng(4646); Xq = Xp.loc[g].to_numpy(float); Xq = np.nan_to_num(Xq); Xn = np.sign(Xq) * np.abs(Xq) ** 1.1 + rng.normal(0, 0.3, Xq.shape); Zn = z(Xn); Cp = Zn.T @ Zb; dia = np.diag(Cp).copy(); np.fill_diagonal(Cp, -2)
ROW += [dict(check="3 positive control: r with own copy, median", value=float(np.median(dia))), dict(check="3 positive control: r with own copy, min", value=float(dia.min())), dict(check="3 positive control: best r with another sample, max", value=float(Cp.max())), dict(check="3 copies above 0.8", value=int((dia > 0.8).sum()), of=len(dia)), dict(check="3 copies whose own original is the best match", value=int((dia > Cp.max(0)).sum()), of=len(dia))]
R = pd.DataFrame(ROW); R.to_csv(f"{O}/E46b_ispy2.csv", index=False); say(R.round(3).to_string(index=False))
