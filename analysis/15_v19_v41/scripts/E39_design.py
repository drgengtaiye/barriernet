"""E39 - design analysis: information and events a signature-versus-CD8A comparison needs.

Protocol: analysis/E39_prereg.md (hash checked). Setup: E34_trials.py up to its fitting stage (E32 setup + NeoPembrOv),
executed unchanged except that its step-0 and coverage tables are written under E39 names, as in E37_test.py.
Bootstrap resamples as E32/E33 (seeds [3032, comparison]); barriernet 0.3 from the frozen weights (hash checked).
Step 0: the pooled paired differences of the 52 signatures from CD8A over the six reproduce E33_paired.csv (< 1e-8).
"""
import hashlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h39 = hashlib.sha256(open(f"{A}/E39_prereg.md", "rb").read()).hexdigest()
assert _h39 == open(f"{A}/E39_prereg.sha256").read().split()[0], "E39 protocol hash mismatch"
_hw = hashlib.sha256(open(f"{A}/results/E37_frozen_weights.tsv", "rb").read()).hexdigest()
assert _hw == open(f"{A}/E37_frozen.sha256").read().split()[0], "frozen weights hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E39_step0_reproduction").replace("E34_coverage_check.csv", "E39_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))
from scipy.stats import spearmanr
print("E39 protocol sha256 verified:", _h39[:12], "| frozen weights sha256 verified:", _hw[:12])

FW = pd.read_csv(f"{O}/E37_frozen_weights.tsv", sep="\t"); WG = dict(zip(FW.gene, FW.weight))


def bn03(k):
    Z = TR[k]["Z"]; nc = set(nonconst(Z))
    pres = [g for g in WG if g in nc and np.isfinite(Z.loc[g].to_numpy(float)).all()]
    if len(pres) < 3 or len(pres) < 0.5 * len(WG):
        return None
    return wscore(Z, pres, np.array([WG[g] for g in pres]))


SC = {}
for k in ALL6:
    SC[("CD8A", k)] = TR[k]["Z"].loc["CD8A"].to_numpy(float)
    for s in LIB:
        if (s, k) in OBS:
            SC[(s, k)] = OBS[(s, k)][0]
    SC[("NetImmune", k)] = OBS[("NetImmune", k)][0]
    b = bn03(k)
    if b is not None:
        SC[("barriernet 0.3", k)] = b

BIDX = {}
for k in ALL6:
    rng_ = np.random.default_rng([3032, KIDX[k]]); a = TR[k]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
    BIDX[k] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]


def run(task):
    sc = SC[task]; d = TR[task[1]]
    try:
        y = m_int(d, sc)[0]
    except Exception:
        y = np.nan
    bs = []
    for idx in BIDX[task[1]]:
        try:
            bs.append(m_int(take(d, idx), sc[idx])[0])
        except Exception:
            bs.append(np.nan)
    return y, np.array(bs)

tasks = list(SC)
with mp.get_context("fork").Pool(NPROC) as pool:
    res = pool.map(run, tasks, chunksize=1)
R = dict(zip(tasks, res))


def within_rho(x, c, arm):
    xr, cr = x.copy(), c.copy()
    for a_ in (0, 1):
        m = arm == a_; xr[m] -= x[m].mean(); cr[m] -= c[m].mean()
    return float(np.corrcoef(xr, cr)[0, 1])


rows = []
for (s, k) in tasks:
    if s == "CD8A":
        continue
    d = TR[k]; arm = d["arm"].astype(int)
    db = R[(s, k)][1] - R[("CD8A", k)][1]; ok = np.isfinite(db)
    v = float(np.var(db[ok], ddof=1)); rho = within_rho(SC[(s, k)], SC[("CD8A", k)], arm)
    if d["kind"] == "cox":
        D = int(d["e"].sum()); p = float(arm.mean()); info0 = D * p * (1 - p); n_a = pi_a = None
    else:
        y = d["y"]; info0 = 1 / sum(1 / ((arm == a_).sum() * y[arm == a_].mean() * (1 - y[arm == a_].mean())) for a_ in (0, 1))
        D = None; p = float(arm.mean())
    vt = 2 * (1 - rho) / info0
    rows.append(dict(score=s, comparison=k, kind=d["kind"], diff=R[(s, k)][0] - R[("CD8A", k)][0], v_boot=v, rho=rho,
                     events=D, p_checkpoint=p, info_unit=info0, v_theory=vt, ratio=v / vt, n_boot_ok=int(ok.sum())))
PC = pd.DataFrame(rows); PC.to_csv(f"{O}/E39_per_comparison.csv", index=False)

# pooled (fixed-effect inverse variance, as E32/E33), scored in >= 4 of the six
Z_ = 1.959964 + 0.841621
pooled = []
for s, g in PC.groupby("score"):
    if len(g) < 4:
        continue
    w = 1 / g.v_boot; D_ = float(np.sum(w * g["diff"]) / w.sum()); se = float(np.sqrt(1 / w.sum()))
    pooled.append(dict(score=s, k=len(g), diff=D_, se=se, info=float(w.sum()), mdd80=Z_ * se, rho_mean=float(g.rho.mean()),
                       rho_median=float(g.rho.median()), ratio_median=float(g.ratio.median())))
PO = pd.DataFrame(pooled).set_index("score"); PO.to_csv(f"{O}/E39_pooled.csv")

E33 = pd.read_csv(f"{O}/E33_paired.csv"); E33 = E33[(E33.set == "six (primary)") & (E33.reference == "CD8A")].set_index("signature")
lib = [s for s in LIB if s in PO.index]
assert set(lib) == set(E33.index), "library coverage differs from E33"
d0 = float(np.max(np.abs(PO.loc[lib, "diff"] - E33.loc[lib, "diff"]))); s0 = float(np.max(np.abs(PO.loc[lib, "se"] - E33.loc[lib, "se"])))
print(f"step 0: 52 signatures vs CD8A, max |diff - E33| {d0:.2e}, max |se - E33| {s0:.2e}"); assert d0 < 1e-8 and s0 < 1e-8

# ------------------------------------------------------------------ outputs
L = PC[PC.score.isin(LIB)]
print("\nratio bootstrap / theory, library signatures, by comparison:")
print(L.groupby("comparison").ratio.describe()[["count", "25%", "50%", "75%"]].round(3).to_string())
print("all score x comparison pairs (library):", L.ratio.describe()[["count", "25%", "50%", "75%"]].round(3).to_dict())
sp = {k: spearmanr(g.v_boot, g.v_theory)[0] for k, g in L.groupby("comparison")}
print("Spearman(v_boot, v_theory) by comparison:", {k: round(v, 3) for k, v in sp.items()})
x1 = float(L.ratio.median()); x2 = float(np.median(list(sp.values())))

mdd = PO.loc[lib, "mdd80"]; rho_med = float(PO.loc[lib, "rho_mean"].median())
print(f"\ncurrent MDD (80%, two-sided 0.05): library median {mdd.median():.3f} (range {mdd.min():.3f}-{mdd.max():.3f}); "
      f"barriernet 0.3 {PO.loc['barriernet 0.3', 'mdd80']:.3f}; NetImmune {PO.loc['NetImmune', 'mdd80']:.3f}")
print(f"median correlation with CD8A (mean over comparisons): library {rho_med:.3f}; barriernet 0.3 "
      f"{PO.loc['barriernet 0.3', 'rho_mean']:.3f}; NetImmune {PO.loc['NetImmune', 'rho_mean']:.3f}")
DELTAS = [0.05, 0.10, 0.15, 0.20]; RHOS = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
info_med = float(PO.loc[lib, "info"].median())
req = []
for dl in DELTAS:
    I_req = Z_ ** 2 / dl ** 2
    r = dict(delta=dl, info_required=I_req)
    for rh in RHOS:
        r[f"events_rho{rh}"] = 8 * (1 - rh) * I_req
    r["events_rho_median_library"] = 8 * (1 - rho_med) * I_req
    r["events_rho_barriernet03"] = 8 * (1 - PO.loc["barriernet 0.3", "rho_mean"]) * I_req
    r["x_current_info_median_library"] = I_req / info_med
    r["x_current_info_barriernet03"] = I_req / PO.loc["barriernet 0.3", "info"]
    req.append(r)
REQ = pd.DataFrame(req); REQ.to_csv(f"{O}/E39_required.csv", index=False)
print("\nrequired information / events (1:1 survival comparison):"); print(REQ.round(1).to_string(index=False))
ev6 = {k: int(TR[k]["e"].sum()) for k in ALL6 if TR[k]["kind"] == "cox"}
print("events in the current survival comparisons:", ev6, "| pCR comparisons n:", {k: len(TR[k]["y"]) for k in ALL6 if TR[k]["kind"] != "cox"})

EXP = pd.DataFrame([
    dict(q="X1", expectation="median ratio of bootstrap to theoretical variance in [0.67, 1.5]", result=f"{x1:.3f}",
         held="yes" if 0.67 <= x1 <= 1.5 else "no"),
    dict(q="X2", expectation="median over comparisons of Spearman(v_boot, v_theory) >= 0.8", result=f"{x2:.3f}",
         held="yes" if x2 >= 0.8 else "no")])
EXP.to_csv(f"{O}/E39_expectations.csv", index=False); print(EXP.to_string(index=False))
