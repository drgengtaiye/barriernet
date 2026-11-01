"""E33 - single-gene reference for CD8A, common-comparison and breast-only checks (fifth review).

Protocol: analysis/E33_prereg.md (sha256 in E33_prereg.sha256, logged 2026-10-24T17:57:57Z, before this script ran).
Setup (data, I-SPY2 MammaPrint covariate, library, observed scores) is E32_trials.py up to its Q1 stage, executed
unchanged except that its step-0 reproduction and coverage tables are written under E33 names. Bootstrap resamples
use the E32 seeds, so the paired differences from CD8A reproduce E32 (checked below).
"""
import hashlib, json
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h33 = hashlib.sha256(open(f"{A}/E33_prereg.md", "rb").read()).hexdigest()
assert _h33 == open(f"{A}/E33_prereg.sha256").read().split()[0], "E33 protocol hash mismatch"
_s32 = open(f"{A}/E32_trials.py").read().split("# ------------------------------------------------------------------ Q1:")[0]
_s32 = _s32.replace("E32_step0_reproduction", "E33_step0_reproduction").replace("E32_coverage.csv", "E33_coverage_check.csv")
exec(compile(_s32, f"{A}/E32_trials.py", "exec"))
print("E33 protocol sha256 verified:", _h33[:12])

T_PANEL = ["CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "CD2", "CD247", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7"]
I_PANEL = ["IFNG", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCR5", "STAT1", "IDO1", "CD274", "PDCD1", "LAG3", "TIGIT",
           "CTLA4", "HAVCR2"]
PANEL = T_PANEL + I_PANEL
FIVE = [k for k in ALL6 if k != "IMM150"]

# ------------------------------------------------------------------ single-gene scores
GOBS, MEAS = {}, {}
for k in ALL6:
    Z = TR[k]["Z"]; nc = nonconst(Z)
    ok = [g for g in nc if np.isfinite(Z.loc[g].to_numpy(float)).all()]
    MEAS[k] = set(ok)
ALLG = sorted(set.intersection(*MEAS.values()))
GENES = sorted(set(ALLG) | {g for g in PANEL if sum(g in MEAS[k] for k in ALL6) >= 4})
print(f"genes measured in all six comparisons: {len(ALLG)}; panel genes missing somewhere:",
      {g: [k for k in ALL6 if g not in MEAS[k]] for g in PANEL if any(g not in MEAS[k] for k in ALL6)})


def gscore(g, k):
    return TR[k]["Z"].loc[g].to_numpy(float)


def fit_any(task):
    """task = (kind, name, k, variant); kind 'gene' or 'sig' (library / extra scores of E32)."""
    kind, s, k, variant = task
    d = CM_OS if variant == "CM025_OS" else TR[k]
    if variant == "hr_only" and k.startswith("ISPY2"):
        d = dict(d); d["cov"] = d["hr"]
    sc = gscore(s, k) if kind == "gene" else OBS[(s, k)][0]
    try:
        if variant == "angio_adj":
            return m_adj(d, sc, zmean(d["Z"], ANGIO))
        return m_int(d, sc)
    except Exception:
        return (np.nan, np.nan)


tasks = [("gene", g, k, "primary") for g in GENES for k in ALL6 if g in MEAS[k]]
tasks += [("gene", g, k, "hr_only") for g in PANEL for k in ("ISPY2_pembro", "ISPY2_durva") if g in MEAS[k]]
tasks += [("gene", g, "CM025", "CM025_OS") for g in PANEL if g in MEAS["CM025"]]
tasks += [("gene", g, k, "angio_adj") for g in PANEL for k in RENAL if g in MEAS[k]]
tasks += [("sig", s, k, "primary") for (s, k) in OBS if k in ALL6 and (s in LIB or s == "CD8A")]
with mp.get_context("fork").Pool(NPROC) as pool:
    res = pool.map(fit_any, tasks, chunksize=8)
PER = pd.DataFrame([dict(kind=kd, name=s, comparison=k, variant=v, y=r[0], se=r[1]) for (kd, s, k, v), r in zip(tasks, res)])
PER["p"] = 2 * norm.sf(np.abs(PER.y / PER.se))
PER.to_csv(f"{O}/E33_per_comparison.csv", index=False)
print(f"fits: {len(tasks)}; failed {int(PER.y.isna().sum())}")

# reproduction: library and CD8A per-comparison values equal E32
P32 = pd.read_csv(f"{O}/E32_per_comparison.csv"); P32 = P32[P32.variant == "primary"].set_index(["signature", "comparison"]).y
chk = PER[(PER.kind == "sig")].set_index(["name", "comparison"]).y
dmax = float(np.nanmax(np.abs(chk.values - P32.reindex(chk.index).values)))
gchk = PER[(PER.kind == "gene") & (PER.name == "CD8A") & (PER.variant == "primary")].set_index("comparison").y
dmax_g = float(np.nanmax(np.abs(gchk.values - P32.reindex([("CD8A", k) for k in gchk.index]).values)))
print(f"step 1 reproduction vs E32: library max |diff| {dmax:.2e}; CD8A as a gene {dmax_g:.2e}")
assert dmax < 1e-8 and dmax_g < 1e-8


def pooled(kind, s, keys, vmap=None, kmin=4):
    vmap = vmap or {}
    parts = []
    for k in keys:
        v = vmap.get(k, "primary"); kk = "CM025" if v == "CM025_OS" else k
        parts.append(PER[(PER.kind == kind) & (PER.name == s) & (PER.comparison == kk) & (PER.variant == v)])
    rows = pd.concat(parts); rows = rows[np.isfinite(rows.y) & np.isfinite(rows.se)]
    if len(rows) < kmin:
        return None
    r = pool_rule(rows.y.values, rows.se.values)
    return dict(est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"], k=len(rows), k_pos=int((rows.y > 0).sum()))


SETS = {"six (primary)": (ALL6, {}, 4), "five without IMmotion150": (FIVE, {}, 4), "breast three": (BREAST, {}, 3),
        "six, CheckMate 025 overall survival": (ALL6, {"CM025": "CM025_OS"}, 4),
        "six, I-SPY2 hormone receptor only": (ALL6, {"ISPY2_pembro": "hr_only", "ISPY2_durva": "hr_only"}, 4),
        "six, renal adjusted for angiogenesis": (ALL6, {k: "angio_adj" for k in RENAL}, 4)}
rows = []
for g in GENES:
    r = pooled("gene", g, ALL6)
    if r:
        rows.append(dict(set="six (primary)", kind="gene", name=g, all_six=g in ALLG,
                         panel="T cell" if g in T_PANEL else "IFN-γ/checkpoint" if g in I_PANEL else "", **r))
for lab, (keys, vmap, kmin) in SETS.items():
    for g in PANEL:
        if lab == "six (primary)":
            continue
        r = pooled("gene", g, keys, vmap, kmin)
        if r:
            rows.append(dict(set=lab, kind="gene", name=g, all_six=g in ALLG,
                             panel="T cell" if g in T_PANEL else "IFN-γ/checkpoint", **r))
    if lab in ("five without IMmotion150", "breast three"):
        for s in list(LIB) + ["CD8A"]:
            r = pooled("sig", s, keys, vmap, kmin)
            if r:
                rows.append(dict(set=lab, kind="library" if s in LIB else "extra", name=s, all_six=None, panel="", **r))
GP = pd.DataFrame(rows); GP.to_csv(f"{O}/E33_pooled.csv", index=False)
six = GP[(GP.set == "six (primary)") & (GP.kind == "gene")].set_index("name")

# ------------------------------------------------------------------ Q1 / Q2
tp = six.loc[[g for g in T_PANEL if g in six.index]].sort_values("est", ascending=False)
rank_cd8a = int(list(tp.index).index("CD8A")) + 1
n_sig_t = int(((tp.p < 0.05) & (tp.est > 0)).sum())
dist = six[six.all_six].est.values
pct = lambda x: float(np.mean(dist < x))
pct_cd8a, pct_tmed = pct(six.loc["CD8A", "est"]), pct(np.median(tp.est))
print(f"Q1: CD8A rank {rank_cd8a}/{len(tp)} among T-cell genes; {n_sig_t}/{len(tp)} at P < 0.05 favoring the checkpoint arm")
print(f"Q2: {len(dist)} genes in all six; CD8A percentile {pct_cd8a:.3f}; T-cell median percentile {pct_tmed:.3f}")

# ------------------------------------------------------------------ Q3 / Q5 / Q6: paired differences (E32 bootstrap)
BIDX = {}
for k in ALL6:
    rng_ = np.random.default_rng([3032, KIDX[k]]); a = TR[k]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
    BIDX[k] = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]


def boot_any(args):
    kind, s, k = args
    sc = gscore(s, k) if kind == "gene" else OBS[(s, k)][0]; d = TR[k]; out = []
    for idx in BIDX[k]:
        try:
            out.append(m_int(take(d, idx), sc[idx])[0])
        except Exception:
            out.append(np.nan)
    return out

btasks = [("sig", s, k) for (s, k) in OBS if k in ALL6 and s in LIB]
btasks += [("gene", g, k) for g in PANEL for k in ALL6 if g in MEAS[k]]
with mp.get_context("fork").Pool(NPROC) as pool:
    bres = pool.map(boot_any, btasks, chunksize=1)
BOOT = {(kd, s, k): np.array(r) for (kd, s, k), r in zip(btasks, bres)}
Y = PER[PER.variant == "primary"].set_index(["kind", "name", "comparison"]).y


def paired(s, ref, keys, kmin):
    ks = [k for k in keys if ("sig", s, k) in BOOT and ("gene", ref, k) in BOOT]
    if len(ks) < kmin:
        return None
    dd, vv = [], []
    for k in ks:
        db = BOOT[("sig", s, k)] - BOOT[("gene", ref, k)]; ok = np.isfinite(db)
        dd.append(Y[("sig", s, k)] - Y[("gene", ref, k)]); vv.append(np.var(db[ok], ddof=1))
    dd, vv = np.array(dd), np.array(vv); w = 1 / vv
    D = np.sum(w * dd) / w.sum(); se = np.sqrt(1 / w.sum())
    return dict(k=len(ks), diff=D, se=se, lo=D - 1.96 * se, hi=D + 1.96 * se, p=2 * norm.sf(abs(D / se)))


PAIR = []
for lab, keys, kmin in [("six (primary)", ALL6, 4), ("five without IMmotion150", FIVE, 4), ("breast three", BREAST, 3)]:
    refs = PANEL if lab == "six (primary)" else ["CD8A"]
    for ref in refs:
        for s in LIB:
            r = paired(s, ref, keys, kmin)
            if r:
                PAIR.append(dict(set=lab, signature=s, reference=ref, **r))
PAIR = pd.DataFrame(PAIR); PAIR.to_csv(f"{O}/E33_paired.csv", index=False)
E32P = pd.read_csv(f"{O}/E32_primary.csv").set_index("signature")
c6 = PAIR[(PAIR.set == "six (primary)") & (PAIR.reference == "CD8A")].set_index("signature")
dmax_q3 = float(np.nanmax(np.abs(c6["diff"] - E32P.loc[c6.index, "diff_vs_CD8A"])))
print(f"Q3 reproduction of E32 paired differences from CD8A: max |diff| {dmax_q3:.2e}")
assert dmax_q3 < 1e-8
q3 = (PAIR[PAIR.set == "six (primary)"].assign(exceeds=lambda x: (x["diff"] > 0) & (x.p < 0.05))
      .groupby("reference").agg(n=("signature", "size"), exceeds=("exceeds", "sum"),
                                below=("p", lambda p: 0)).drop(columns="below"))
below = PAIR[(PAIR.set == "six (primary)") & (PAIR["diff"] < 0) & (PAIR.p < 0.05)].groupby("reference").size()
q3["below"] = below.reindex(q3.index).fillna(0).astype(int)
q3 = q3.join(six[["est", "lo", "hi", "p"]], how="left")
q3.to_csv(f"{O}/E33_Q3_by_reference.csv")
print("Q3 by reference gene:\n", q3.sort_values("est", ascending=False).round(3).to_string())
n_ok_t = int((q3.loc[[g for g in T_PANEL if g in q3.index], "exceeds"] <= 2).sum())

# ------------------------------------------------------------------ Q4 / Q5 / Q6
L32 = E32P[E32P.kind == "library"]
se32 = (L32.cd8a_hi - L32.cd8a_lo) / (2 * 1.96)
q4 = dict(max_upper=float(L32.cd8a_hi.max()), max_upper_sig=str(L32.cd8a_hi.idxmax()),
          median_upper=float(L32.cd8a_hi.median()), median_se=float(se32.median()), mdd80=float(2.8 * se32.median()),
          cd8a_est=float(six.loc["CD8A", "est"]))
print("Q4:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in q4.items()})
c5 = PAIR[(PAIR.set == "five without IMmotion150") & (PAIR.reference == "CD8A")]
q5n = int(((c5["diff"] > 0) & (c5.p < 0.05)).sum())
br = GP[(GP.set == "breast three") & (GP.kind == "library")]
q6n = int(((br.p < 0.05) & (br.est > 0)).sum())
cb = PAIR[(PAIR.set == "breast three") & (PAIR.reference == "CD8A")]
q6x = int(((cb["diff"] > 0) & (cb.p < 0.05)).sum())
print(f"Q5: {q5n}/{len(c5)} exceed CD8A over the five; Q6: breast {q6n}/{len(br)} favor at P < 0.05; "
      f"{q6x}/{len(cb)} exceed CD8A there; CD8A breast:",
      GP[(GP.set == "breast three") & (GP.name == "CD8A") & (GP.kind == "gene")][["est", "lo", "hi", "p"]].round(3).to_dict("records"))

# ------------------------------------------------------------------ Q7: pooled reference draws per signature (for audited directions)
REFD = np.load(f"{O}/E32_reference_draws.npz")
estb = {}
for s in LIB:
    ks = [k for k in ALL6 if (s, k) in OBS]
    if len(ks) < 4:
        continue
    D = np.stack([REFD[f"{k}|{len(OBS[(s, k)][1])}|{zlib.crc32(repr(tuple(np.round(np.sort(OBS[(s, k)][2]), 12))).encode())}"]
                  for k in ks], 1)
    estb[s] = np.array([pool_rule(D[b, :, 0], D[b, :, 1])["est"] for b in range(B_REF)])
np.savez_compressed(f"{O}/E33_pooled_reference_estimates.npz", **estb)
chk7 = []
for s, eb in estb.items():
    dr = LIB[s]["dir"]; obs = E32P.loc[s, "est"]
    chk7.append(abs((1 + np.sum(dr * eb >= dr * obs)) / (B_REF + 1) - E32P.loc[s, "emp_p_declared"]))
print(f"Q7 setup: pooled reference estimates for {len(estb)} signatures; max |emp P - E32| {max(chk7):.2e}")
assert max(chk7) < 1e-12

EXP = pd.DataFrame([
    dict(q="Q1 (primary)", expectation="CD8A not ranked first among the T-cell genes, and >= 7 of 13 at P < 0.05 favoring the checkpoint arm",
         result=f"CD8A rank {rank_cd8a}/{len(tp)}; {n_sig_t}/{len(tp)} at P < 0.05", held=(rank_cd8a > 1) and (n_sig_t >= 7)),
    dict(q="Q2", expectation="CD8A and the T-cell panel median above the 90th percentile of all genes measured in all six",
         result=f"{len(dist)} genes; CD8A {pct_cd8a:.1%}; T-cell median {pct_tmed:.1%}", held=(pct_cd8a > 0.9) and (pct_tmed > 0.9)),
    dict(q="Q3", expectation="for >= 7 of the 13 T-cell genes, at most 2 of 52 signatures exceed the gene at P < 0.05",
         result=f"{n_ok_t}/{len(tp)} genes", held=n_ok_t >= 7),
    dict(q="Q4", expectation="descriptive", result=json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in q4.items()}), held=None),
    dict(q="Q5", expectation="at most 2 of 52 exceed CD8A over the five comparisons without IMmotion150", result=f"{q5n}/{len(c5)}", held=q5n <= 2),
    dict(q="Q6", expectation="fewer than 10 of 52 signatures at P < 0.05 favoring the checkpoint arm in the three breast comparisons",
         result=f"{q6n}/{len(br)}; exceeding CD8A {q6x}/{len(cb)}", held=q6n < 10),
])
EXP.to_csv(f"{O}/E33_expectations.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 140)
print("\n== expectations\n", EXP.to_string(index=False))
print("\n== T-cell and IFN panels, six comparisons\n",
      six.loc[[g for g in PANEL if g in six.index], ["panel", "k", "est", "lo", "hi", "p", "k_pos"]].sort_values("est", ascending=False).round(3).to_string())
