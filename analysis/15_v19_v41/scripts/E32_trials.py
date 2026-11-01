"""E32 - calibrating signature x treatment interaction in the six randomized checkpoint comparisons.

Protocol: analysis/E32_prereg.md (sha256 in E32_prereg.sha256, logged 2026-10-22T18:01:38Z, before this script ran).
Data, models and pooling are E28_base.py (E25/E21 code), executed with its step-0 reproduction gate (every published
interaction and within-arm value re-derived, < 1e-6) and its output redirected to results/E32_step0_reproduction.csv.

Library: the companion's 52 signatures (SignatureSets IO gene lists and weights, family from icb/v5_signature_clusters,
declared direction from icb/sig_methods_clean), plus ESTIMATE, TMEscore (z-mean), CD8A and the barriernet modules
and axes (published scores). Common engine: genes z-scored within each comparison across arms; score = sum(w z) /
sum(|w|) over the genes present (published weights; weight 1 where none is given), NaN-aware per sample; scored when
>= 3 genes and >= 50% of the list are present (non-constant); pooled when scored in >= 4 of the 6 comparisons.

Model change for E32 (protocol): I-SPY2 covariates are hormone receptor and MammaPrint class (0 = MP1, 1 = MP2).

Stages, all in this run: Q1 per-comparison and pooled interactions; Q2 size-matched random reference (1,000 draws per
signature and comparison from that comparison's complete non-constant genes, the signature's weights in random order,
draw b pooled across comparisons; seed 2032); Q3 paired difference from CD8A (300 within-arm patient bootstrap
resamples per comparison, common to all scores; inverse-variance pooled); Q4 interaction adjusted for the ESTIMATE
immune score and its interaction with arm; Q5/Q6 I-SPY2 MammaPrint and processing ("adjustment method") checks;
sensitivity sets (five comparisons without IMmotion150; + PRINCE and DUTRENEO; I-SPY2 hormone receptor only;
CheckMate 025 overall survival; renal trials adjusted for the angiogenesis score).
The published-engine sensitivity is E32_native.py; the locked gastric immune/stromal null (Q7) is E32_locked_bionull.R.
"""
import hashlib, zlib, json
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h = hashlib.sha256(open(f"{A}/E32_prereg.md", "rb").read()).hexdigest()
assert _h == open(f"{A}/E32_prereg.sha256").read().split()[0], "protocol hash mismatch"
_src = open(f"{A}/E28_base.py").read()
_src = _src.replace('f"{O}/E28_reproduction_check.csv"', 'f"{O}/E32_step0_reproduction.csv"')
assert "E32_step0_reproduction" in _src
exec(compile(_src, f"{A}/E28_base.py", "exec"))
from scipy.stats import binomtest
print("E32 protocol sha256 verified:", _h[:12])

G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
B_REF, B_BOOT, NPROC = 1000, 300, 96
KIDX = {k: i for i, k in enumerate(ALL6 + ["PRINCE", "DUTRENEO"])}

# ------------------------------------------------------------------ I-SPY2: MammaPrint and processing
def soft(acc):
    out = []
    for r in open(f"{G}/validation/{acc}_gsm.soft", encoding="utf-8", errors="ignore").read().split("^SAMPLE = ")[1:]:
        d = {"gsm": r.split("\n", 1)[0].strip()}
        for l in r.splitlines():
            if l.startswith("!Sample_characteristics_ch1 = ") and ": " in l:
                k_, v = l.split(" = ", 1)[1].split(": ", 1); d[k_] = v
        out.append(d)
    return pd.DataFrame(out)

S194 = soft("GSE194040").drop_duplicates("patient id").set_index("patient id")
ids = TR["ISPY2_pembro"]["Z"].columns.astype(str)
mp_p = pd.to_numeric(S194.loc[ids, "mp"], errors="coerce").values
ADJ_P = S194.loc[ids, "adjustment method"].values
BIO = pd.read_csv(f"{G}/validation/raw/GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz")
BIO["ResearchID"] = BIO.ResearchID.astype(str)
ids_d = TR["ISPY2_durva"]["Z"].columns.astype(str)
mp_d = pd.to_numeric(BIO.drop_duplicates("ResearchID").set_index("ResearchID").loc[ids_d, "MP.1.2..class..0..MP1..1..MP2."],
                     errors="coerce").values
for k, mp_ in [("ISPY2_pembro", mp_p), ("ISPY2_durva", mp_d)]:
    assert np.isfinite(mp_).all(), f"{k}: MammaPrint missing for {int((~np.isfinite(mp_)).sum())} patients"
    d = TR[k]; d["hr"] = d["cov"].copy(); d["mp"] = mp_
    d["cov"] = np.column_stack([d["hr"], mp_])                      # E32 primary: hormone receptor + MammaPrint
    print(f"{k}: MP2 {int(mp_.sum())}/{len(mp_)}; arm x MP2 = {int(mp_[d['arm'] == 1].sum())}/{int((d['arm'] == 1).sum())} "
          f"vs {int(mp_[d['arm'] == 0].sum())}/{int((d['arm'] == 0).sum())}")
print("I-SPY2 pembrolizumab, arm x adjustment method:", pd.crosstab(TR["ISPY2_pembro"]["arm"], ADJ_P).to_dict())

# ------------------------------------------------------------------ library
ALIAS_ = {"WARS1": "WARS", "HLA_DPB1": "HLA-DPB1"}                # gene aliases as in cmp_scores.py
SIGDIR = f"{G}/icb/SignatureSets-main/data-raw/IO"
FAM = pd.read_csv(f"{G}/icb/v5_signature_clusters.csv").set_index("signature").family
METH = pd.read_csv(f"{G}/icb/sig_methods_clean.csv").set_index("signature")
LIB, wnote = {}, []
for s in FAM.index:
    d = pd.read_csv(f"{SIGDIR}/{s}.csv", encoding="latin1")
    g = d.gene_name.astype(str).map(lambda x: ALIAS_.get(x, x))
    w = pd.to_numeric(d["weight"], errors="coerce") if "weight" in d else pd.Series(np.nan, index=d.index)
    if w.isna().all():
        w[:] = 1.0
    elif w.isna().any():
        wnote.append(f"{s}: {int(w.isna().sum())} of {len(w)} weights missing, set to 1"); w = w.fillna(1.0)
    t_ = pd.DataFrame({"g": g, "w": w.values})
    t_ = t_[(t_.g != "") & (t_.g != "nan")].drop_duplicates("g")
    LIB[s] = dict(genes=t_.g.tolist(), w=t_.w.values.astype(float), family=FAM[s],
                  dir=1 if METH.loc[s, "association"] == "sensitive" else -1, engine=METH.loc[s, "method"], kind="library")
print(f"library: {len(LIB)} signatures; families {pd.Series({s: v['family'] for s, v in LIB.items()}).value_counts().to_dict()}; "
      f"declared sensitive {sum(v['dir'] == 1 for v in LIB.values())}")
for x in wnote:
    print("  weights:", x)
CG = pd.read_csv(f"{O}/E23_comparator_genes.csv")
cg = lambda n: list(dict.fromkeys(ALIAS_.get(x, x) for x in CG.gene[CG.set == n]))
# published scores (observed) with a group structure for the random reference: list of (genes, coefficient)
EXTRA = {"ESTIMATE": ([cg("ESTIMATE_immune"), cg("ESTIMATE_stromal")], [1, -1], 1),
         "TMEscore_z": ([cg("TMEscore_A"), cg("TMEscore_B")], [1, -1], 1),
         "CD8A": ([["CD8A"]], [1], 1)}
for m_, G_ in MOD.items():
    EXTRA[m_] = ([G_], [1], -1 if m_ in ("CAF_ECM", "TAM_SUP") else 1)
EXTRA["Barrier"] = ([MOD["CAF_ECM"], MOD["TAM_SUP"]], [0.5, 0.5], -1)
EXTRA["Immune"] = ([MOD[m] for m in ("APM", "IFN", "TLS_B", "T_CELL")], [0.25] * 4, 1)
EXTRA["NetImmune"] = ([MOD[m] for m in ("APM", "IFN", "TLS_B", "T_CELL", "CAF_ECM", "TAM_SUP")], [0.25] * 4 + [-0.5] * 2, 1)
ESTIMMUNE = cg("ESTIMATE_immune"); ANGIO = ["VEGFA", "KDR", "ESM1", "PECAM1", "ANGPTL4", "CD34"]


def nonconst(Z):
    sd = Z.std(1, skipna=True)
    return Z.index[(sd > 0) & sd.notna()]


def wscore(Z, genes, w):
    """sum(w z) / sum(|w|) over the genes present, NaN-aware per sample."""
    M = Z.loc[genes].to_numpy(float); W = np.where(np.isfinite(M), 1.0, 0.0) * np.abs(w)[:, None]
    return np.nansum(M * w[:, None], 0) / W.sum(0)


def zmean(Z, genes):
    g = [x for x in genes if x in Z.index]
    return Z.loc[g].mean(0).to_numpy(float) if len(g) >= 3 else None


# ------------------------------------------------------------------ sensitivity comparisons: PRINCE, DUTRENEO, CM025 OS
C = pd.read_csv(f"{P5}/v2_prince_clinical_api.csv", index_col=0)
M = pd.read_csv(f"{P5}/ici/paad_iatlas_prince_2022.csv", index_col=0)
keep = [s for s in M.columns if s in C.index]; C, M = C.loc[keep], M[keep]
Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T; S = scores(Z)
t = pd.to_numeric(C["PFS_MONTHS"], errors="coerce")
e = pd.to_numeric(C["PFS_STATUS"].astype(str).str.extract(r"^\s*([01])\s*:")[0], errors="coerce").fillna(0).astype(int)
nivo = C["ICI_RX"].astype(str).str.contains("Nivolumab", case=False, na=False)
sel = (t.notna() & (t > 0)).values
TR["PRINCE"] = dict(kind="cox", Z=Z.loc[:, sel], X=M.loc[:, sel], S=S.loc[sel], arm=nivo.values[sel].astype(int),
                    t=t.values[sel].astype(float), e=e.values[sel].astype(int), firth=True)
Md = pd.read_csv(f"{V21}/DUTRENEO_GSE328930/GSE328930_metadata.tsv.gz", sep="\t", index_col=0)
Xd = pd.read_csv(f"{V21}/DUTRENEO_GSE328930/GSE328930_counts.tsv.gz", sep="\t", index_col=0)
cc = Md[Md.Group.isin(["HOT-ICI", "HOT-NAC"])]; lib = Xd.sum(0)
Xd = np.log2(Xd[cc.index] / lib[cc.index] * 1e6 + 1); Xd, _ = collapse(Xd)
Z = rowz(Xd); S = scores(Z)
TR["DUTRENEO"] = dict(kind="logit", Z=Z, X=Xd, S=S, arm=(cc.Group == "HOT-ICI").astype(int).values.astype(float),
                      y=(cc.Pathological_response == "Pathological complete response").astype(float).values, cov=None, firth=True)
C25 = pd.read_csv(f"{V21}/Braun2020_CheckMate/S1_clinical.tsv", sep="\t")
C25 = C25[C25.Cohort == "CM-025"].set_index("RNA_ID").loc[TR["CM025"]["Z"].columns]
CM_OS = dict(TR["CM025"]); CM_OS["t"] = pd.to_numeric(C25.OS, errors="coerce").values
CM_OS["e"] = (pd.to_numeric(C25.OS_CNSR, errors="coerce").values == 1).astype(int)
assert np.isfinite(CM_OS["t"]).all() and (CM_OS["t"] > 0).all()
KEYS = ALL6 + ["PRINCE", "DUTRENEO"]
for k in ["PRINCE", "DUTRENEO"]:
    d = TR[k]; print(f"  {k:12s} checkpoint {int(d['arm'].sum())} control {int(len(d['arm']) - d['arm'].sum())} genes {d['Z'].shape[0]}")


def take(d, idx):
    out = dict(d); out["arm"] = d["arm"][idx]
    for f in ["t", "e", "y"]:
        if f in d: out[f] = d[f][idx]
    if d.get("cov") is not None: out["cov"] = d["cov"][idx]
    return out


def m_adj(d, s, x):
    """Interaction of s with arm, adjusted for x and its interaction with arm (benefit-oriented y, se)."""
    s, x = std(s), std(x); a = d["arm"]
    if d["kind"] == "cox":
        b, se = cox({"s": s, "x": x, "arm": a, "sx": s * a, "xx": x * a}, d["t"], d["e"], ["sx"])[0]; return -b, se
    X = np.column_stack([np.ones_like(s), s, x, a, s * a, x * a] + ([d["cov"]] if d["cov"] is not None else []))
    return lfit(X, d["y"], [4], d["firth"])[0]


# ------------------------------------------------------------------ observed scores
COV, OBS = [], {}                      # OBS[(name, k)] = score vector (finite for all patients of k)
for k in KEYS:
    Z = TR[k]["Z"]; nc = set(nonconst(Z))
    for s, v in LIB.items():
        pres = [g for g in v["genes"] if g in nc]
        ok = len(pres) >= 3 and len(pres) >= 0.5 * len(v["genes"])
        COV.append(dict(signature=s, comparison=k, genes=len(v["genes"]), present=len(pres), scored=ok))
        if ok:
            w = v["w"][[v["genes"].index(g) for g in pres]]
            sc = wscore(Z, pres, w)
            if np.isfinite(sc).all():
                OBS[(s, k)] = (sc, pres, w)
    for s, (groups, coef, _) in EXTRA.items():
        sc = TR[k]["S"][s].to_numpy(float) if s in TR[k]["S"] else None
        if sc is not None and np.isfinite(sc).all():
            pres_g = [[g for g in grp if g in nc] for grp in groups]
            w = np.concatenate([np.full(len(pg), c / len(pg)) for pg, c in zip(pres_g, coef) if len(pg)])
            OBS[(s, k)] = (sc, sum(pres_g, []), w)
COVDF = pd.DataFrame(COV); COVDF.to_csv(f"{O}/E32_coverage.csv", index=False)
NAMES = list(LIB) + list(EXTRA)
meta = lambda s: (LIB[s]["family"], LIB[s]["dir"], "library") if s in LIB else ("barriernet / comparator", EXTRA[s][2], "extra")


# ------------------------------------------------------------------ Q1: observed interactions (primary model) + sensitivities
def fit_obs(task):
    s, k, variant = task
    d = CM_OS if variant == "CM025_OS" else TR[k]
    if variant == "hr_only" and k.startswith("ISPY2"):
        d = dict(d); d["cov"] = d["hr"]
    sc = OBS[(s, k)][0]
    if variant == "estimate_adj":
        x = zmean(d["Z"], ESTIMMUNE); return m_adj(d, sc, x)
    if variant == "angio_adj":
        x = zmean(d["Z"], ANGIO); return m_adj(d, sc, x)
    return m_int(d, sc)

tasks = [(s, k, "primary") for (s, k) in OBS]
tasks += [(s, k, "hr_only") for (s, k) in OBS if k.startswith("ISPY2")]
tasks += [(s, "CM025", "CM025_OS") for (s, k) in OBS if k == "CM025"]
tasks += [(s, k, "estimate_adj") for (s, k) in OBS if k in ALL6]
tasks += [(s, k, "angio_adj") for (s, k) in OBS if k in RENAL]
with mp.get_context("fork").Pool(NPROC) as pool:
    res = pool.map(fit_obs, tasks, chunksize=4)
PER = pd.DataFrame([dict(signature=s, comparison=k, variant=v, y=r[0], se=r[1]) for (s, k, v), r in zip(tasks, res)])
PER["p"] = 2 * norm.sf(np.abs(PER.y / PER.se))
PER.to_csv(f"{O}/E32_per_comparison.csv", index=False)
print(f"Q1: {len(tasks)} interaction fits")


def pooled_row(s, rows):
    rows = rows[np.isfinite(rows.y) & np.isfinite(rows.se)]
    r = pool_rule(rows.y.values, rows.se.values); kpos = int((rows.y > 0).sum()); k = len(rows)
    return dict(est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"], k=k, k_pos=kpos,
                sign_p=binomtest(kpos, k, 0.5).pvalue if k else np.nan, comparisons=" ".join(rows.comparison))


def pool_set(variant_map, keys, label):
    out = []
    for s in NAMES:
        parts = []
        for k in keys:
            v = variant_map.get(k, "primary")
            kk = "CM025" if v == "CM025_OS" else k
            r = PER[(PER.signature == s) & (PER.comparison == kk) & (PER.variant == v)]
            parts.append(r)
        rows = pd.concat(parts) if parts else pd.DataFrame()
        if len(rows) >= 4 or (meta(s)[2] == "extra" and len(rows) >= 2):
            fam, dr, kind = meta(s)
            out.append(dict(set=label, signature=s, family=fam, declared_dir=dr, kind=kind, **pooled_row(s, rows)))
    return out

POOL = []
POOL += pool_set({}, ALL6, "six (primary)")
POOL += pool_set({}, [k for k in ALL6 if k != "IMM150"], "five without IMmotion150")
POOL += pool_set({}, KEYS, "six + PRINCE + DUTRENEO")
POOL += pool_set({"ISPY2_pembro": "hr_only", "ISPY2_durva": "hr_only"}, ALL6, "six, I-SPY2 hormone receptor only")
POOL += pool_set({"CM025": "CM025_OS"}, ALL6, "six, CheckMate 025 overall survival")
POOL += pool_set({k: "angio_adj" for k in RENAL}, ALL6, "six, renal adjusted for angiogenesis")
POOL += pool_set({k: "estimate_adj" for k in ALL6}, ALL6, "six, adjusted for ESTIMATE immune (Q4)")
POOL = pd.DataFrame(POOL)

# ------------------------------------------------------------------ Q2: size-matched random reference
UNIV = {}
for k in ALL6:
    Z = TR[k]["Z"]; g = [x for x in nonconst(Z) if np.isfinite(Z.loc[x].to_numpy(float)).all()]
    UNIV[k] = Z.loc[g].to_numpy(float)
REFKEYS = {}
for (s, k), (_, pres, w) in OBS.items():
    if k in ALL6 and s != "CD8A":
        wk = tuple(np.round(np.sort(w), 12)); REFKEYS.setdefault((k, len(pres), wk), []).append(s)
print(f"Q2: {len(REFKEYS)} distinct (comparison, size, weights) reference keys x {B_REF} draws")


def ref_task(args):
    k, n, wk, b0, b1 = args
    w = np.array(wk); U = UNIV[k]; d = TR[k]
    rng = np.random.default_rng([2032, KIDX[k], n, zlib.crc32(repr(wk).encode())])
    out = []
    for b in range(b1):
        idx = rng.choice(U.shape[0], n, replace=False); ww = w[rng.permutation(n)]
        if b >= b0:
            try:
                out.append(m_int(d, U[idx].T @ ww))
            except Exception:                                            # non-converging fit: recorded as missing
                out.append((np.nan, np.nan))
    return out

rtasks = [(k, n, wk, b0, b0 + 250) for (k, n, wk) in REFKEYS for b0 in range(0, B_REF, 250)]
with mp.get_context("fork").Pool(NPROC) as pool:
    rres = pool.map(ref_task, rtasks, chunksize=1)
REF = {}
for (k, n, wk, b0, b1), r in zip(rtasks, rres):
    REF.setdefault((k, n, wk), [None] * (B_REF // 250))[b0 // 250] = r
REF = {key: np.array(sum(v, [])) for key, v in REF.items()}            # (B, 2): y, se
_nf = sum(int((~np.isfinite(v[:, 0])).sum()) for v in REF.values())
print(f"Q2: failed reference fits {_nf} of {len(REF) * B_REF} (missing; pooled over the remaining comparisons of that draw)")
np.savez_compressed(f"{O}/E32_reference_draws.npz", **{f"{k}|{n}|{zlib.crc32(repr(wk).encode())}": v for (k, n, wk), v in REF.items()})

EMP = []
for s in NAMES:
    if s == "CD8A":
        continue
    ks = [k for k in ALL6 if (s, k) in OBS]
    if len(ks) < (4 if s in LIB else 2):
        continue
    draws = []
    for k in ks:
        _, pres, w = OBS[(s, k)]; draws.append(REF[(k, len(pres), tuple(np.round(np.sort(w), 12)))])
    D = np.stack(draws, 1)                                               # (B, k, 2)
    est_b = np.array([pool_rule(D[b, :, 0], D[b, :, 1])["est"] for b in range(B_REF)])
    obs = POOL[(POOL.set == "six (primary)") & (POOL.signature == s)].iloc[0]
    dr = meta(s)[1]
    emp = (1 + np.sum(dr * est_b >= dr * obs.est)) / (B_REF + 1)
    nominal_b = np.array([pool_rule(D[b, :, 0], D[b, :, 1])["p"] for b in range(0, B_REF, 4)])
    EMP.append(dict(signature=s, ref_k=len(ks), emp_p_declared=emp, ref_median=np.nanmedian(est_b),
                    ref_p95_declared=np.nanpercentile(dr * est_b, 95) * dr, ref_frac_nominal_p05=np.mean(nominal_b < 0.05)))
EMP = pd.DataFrame(EMP)

# ------------------------------------------------------------------ Q3: paired difference from CD8A (bootstrap)
BIDX = {}
for k in ALL6:
    rng = np.random.default_rng([3032, KIDX[k]]); a = TR[k]["arm"]; i1, i0 = np.where(a == 1)[0], np.where(a == 0)[0]
    BIDX[k] = [np.concatenate([rng.choice(i1, len(i1)), rng.choice(i0, len(i0))]) for _ in range(B_BOOT)]


def boot_task(args):
    s, k = args
    sc = OBS[(s, k)][0]; d = TR[k]; out = []
    for idx in BIDX[k]:
        try:
            out.append(m_int(take(d, idx), sc[idx])[0])
        except Exception:
            out.append(np.nan)
    return out

btasks = [(s, k) for (s, k) in OBS if k in ALL6]
with mp.get_context("fork").Pool(NPROC) as pool:
    bres = pool.map(boot_task, btasks, chunksize=1)
BOOT = {t_: np.array(r) for t_, r in zip(btasks, bres)}
Q3 = []
for s in NAMES:
    if s == "CD8A":
        continue
    ks = [k for k in ALL6 if (s, k) in OBS and ("CD8A", k) in OBS]
    if len(ks) < (4 if s in LIB else 2):
        continue
    dd, vv = [], []
    for k in ks:
        y1 = PER[(PER.signature == s) & (PER.comparison == k) & (PER.variant == "primary")].y.iloc[0]
        y0 = PER[(PER.signature == "CD8A") & (PER.comparison == k) & (PER.variant == "primary")].y.iloc[0]
        diff_b = BOOT[(s, k)] - BOOT[("CD8A", k)]; ok = np.isfinite(diff_b)
        dd.append(y1 - y0); vv.append(np.var(diff_b[ok], ddof=1))
    dd, vv = np.array(dd), np.array(vv); wgt = 1 / vv
    D_ = np.sum(wgt * dd) / wgt.sum(); se_ = np.sqrt(1 / wgt.sum())
    Q3.append(dict(signature=s, k=len(ks), diff_vs_CD8A=D_, lo=D_ - 1.96 * se_, hi=D_ + 1.96 * se_, p=2 * norm.sf(abs(D_ / se_))))
Q3 = pd.DataFrame(Q3)

# ------------------------------------------------------------------ Q5 / Q6: I-SPY2
Q56 = []
for sname in ["NetImmune", "Immune", "CAF_ECM", "ESTIMATE", "CD8A"]:
    for lab, vmap in [("HR + MammaPrint (E32 primary)", {}), ("hormone receptor only (published)", {"ISPY2_pembro": "hr_only", "ISPY2_durva": "hr_only"})]:
        r = pool_set(vmap, ALL6, lab); r = [x for x in r if x["signature"] == sname][0]
        Q56.append(dict(question="Q5", score=sname, model=lab, est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"], k_pos=r["k_pos"], k=r["k"]))
d = TR["ISPY2_pembro"]; lin = np.array([str(x).startswith("Linear") for x in ADJ_P]); keep = lin | (d["arm"] == 1)
print(f"Q6: pembrolizumab {int((d['arm'] == 1).sum())} (Linear {int((lin & (d['arm'] == 1)).sum())}); controls kept {int((keep & (d['arm'] == 0)).sum())}")
dk = take(d, np.where(keep)[0])
Xk = d["X"].loc[:, keep]; Zk = ((Xk.T - Xk.mean(1)) / Xk.std(1).replace(0, 1)).T; Sk = scores(Zk)
for sname in ["NetImmune", "Immune", "CAF_ECM", "ESTIMATE", "CD8A"]:
    for lab, sc in [("all controls", d["S"][sname].to_numpy(float)),
                    ("Linear-adjusted controls, published z", d["S"][sname].to_numpy(float)[keep]),
                    ("Linear-adjusted controls, z within subset", Sk[sname].to_numpy(float))]:
        dd_ = d if lab == "all controls" else dk
        y, se = m_int(dd_, sc)
        Q56.append(dict(question="Q6", score=sname, model=lab, est=y, lo=y - 1.96 * se, hi=y + 1.96 * se, p=2 * norm.sf(abs(y / se))))
    # six-comparison pool with the subset in place of all pembrolizumab controls (z within subset)
    rows = PER[(PER.signature == sname) & (PER.variant == "primary") & (PER.comparison.isin(ALL6)) & (PER.comparison != "ISPY2_pembro")]
    y, se = m_int(dk, Sk[sname].to_numpy(float))
    rows = pd.concat([rows, pd.DataFrame([dict(signature=sname, comparison="ISPY2_pembro (Linear controls)", y=y, se=se)])])
    r = pooled_row(sname, rows)
    Q56.append(dict(question="Q6", score=sname, model="six-comparison pool, pembrolizumab with Linear-adjusted controls",
                    est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"], k_pos=r["k_pos"], k=r["k"]))
Q56 = pd.DataFrame(Q56)

# ------------------------------------------------------------------ outputs and expectations
PRIM = POOL[POOL.set == "six (primary)"].merge(EMP, on="signature", how="left").merge(
    Q3.rename(columns={"k": "k_cd8a", "lo": "cd8a_lo", "hi": "cd8a_hi", "p": "cd8a_p"}), on="signature", how="left")
adj = POOL[POOL.set == "six, adjusted for ESTIMATE immune (Q4)"][["signature", "est", "p"]].rename(columns={"est": "est_estimate_adj", "p": "p_estimate_adj"})
PRIM = PRIM.merge(adj, on="signature", how="left")
POOL.to_csv(f"{O}/E32_pooled_all_sets.csv", index=False); PRIM.to_csv(f"{O}/E32_primary.csv", index=False)
Q56.to_csv(f"{O}/E32_ispy2.csv", index=False)
L = PRIM[PRIM.kind == "library"]; NLIB = len(LIB)
inflam = L[L.family == "Inflamed/IFN"]
q1a = (inflam.est > 0).mean(); sig1 = L[L.p < 0.05]
q1b = len(sig1) / NLIB; q1c = (sig1.est > 0).all()
q2n = int((L.emp_p_declared < 0.05).sum())
q3n = int(((L.diff_vs_CD8A > 0) & (L.cd8a_p < 0.05)).sum())
keepq4 = int((sig1.p_estimate_adj < 0.05).sum())
net5 = Q56[(Q56.question == "Q5") & (Q56.score == "NetImmune") & (Q56.model.str.startswith("HR + Mamma"))].iloc[0]
q6 = Q56[(Q56.question == "Q6") & (Q56.model == "Linear-adjusted controls, z within subset")].set_index("score")
EXP = pd.DataFrame([
    dict(q="Q1a", expectation=">= 70% of Inflamed/IFN favor the checkpoint arm", result=f"{(inflam.est > 0).sum()}/{len(inflam)} ({q1a:.0%}) of those pooled", held=q1a >= 0.7),
    dict(q="Q1b", expectation="10-40% of the 52 reach P < 0.05, all favoring the checkpoint arm",
         result=f"{len(sig1)}/{NLIB} ({q1b:.0%}); pooled {len(L)}; all favoring: {q1c}", held=(0.10 <= q1b <= 0.40) and bool(q1c)),
    dict(q="Q2 (primary)", expectation="fewer than 20% of the 52 exceed their size-matched random reference (declared direction)",
         result=f"{q2n}/{NLIB} ({q2n / NLIB:.0%}); evaluable {int(L.emp_p_declared.notna().sum())}", held=q2n < 0.2 * NLIB),
    dict(q="Q3", expectation="no signature exceeds CD8A at P < 0.05 (at most 2 of 52)", result=f"{q3n}/{NLIB}", held=q3n <= 2),
    dict(q="Q4", expectation="fewer than one-third of the Q1-significant keep P < 0.05 after ESTIMATE-immune adjustment",
         result=f"{keepq4}/{len(sig1)}", held=(keepq4 < len(sig1) / 3) if len(sig1) else None),
    dict(q="Q5", expectation="NetImmune six-comparison pool keeps its sign with I-SPY2 adjusted for HR + MammaPrint",
         result=f"{net5.est:+.3f} ({net5.lo:.2f} to {net5.hi:.2f}); P = {net5.p:.3g}", held=net5.est > 0),
    dict(q="Q6", expectation="NetImmune and Immune keep their sign in I-SPY2 pembrolizumab with Linear-adjusted controls only",
         result=f"NetImmune {q6.loc['NetImmune', 'est']:+.3f} (P = {q6.loc['NetImmune', 'p']:.3g}); Immune {q6.loc['Immune', 'est']:+.3f} (P = {q6.loc['Immune', 'p']:.3g})",
         held=(q6.loc["NetImmune", "est"] > 0) and (q6.loc["Immune", "est"] > 0)),
])
EXP.to_csv(f"{O}/E32_expectations.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 120)
print("\n== expectations\n", EXP.to_string(index=False))
cols = ["signature", "family", "declared_dir", "k", "est", "lo", "hi", "p", "k_pos", "sign_p", "emp_p_declared", "ref_frac_nominal_p05",
        "diff_vs_CD8A", "cd8a_p", "est_estimate_adj", "p_estimate_adj"]
print("\n== six-comparison results (sorted by P)\n", PRIM.sort_values("p")[cols].round(4).to_string(index=False))
print("\n== I-SPY2 (Q5/Q6)\n", Q56.round(4).to_string(index=False))
