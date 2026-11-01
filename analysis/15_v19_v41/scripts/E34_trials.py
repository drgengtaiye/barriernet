"""E34 - NeoPembrOV (GSE227666) as a seventh randomized comparison (sensitivity).

Protocol: analysis/E34_prereg.md (sha256 in E34_prereg.sha256, logged 2026-10-25T16:36:59Z, before this script ran).
Setup is E32_trials.py up to its Q1 stage (as in E33). The six comparisons' per-comparison values are read from the
E32/E33 result files; paired differences from CD8A over seven (or four) comparisons are inverse-variance combinations
of the six- (or three-breast-) comparison pooled difference with NeoPembrOV's, which equals pooling all of them at once.
"""
import hashlib
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h34 = hashlib.sha256(open(f"{A}/E34_prereg.md", "rb").read()).hexdigest()
assert _h34 == open(f"{A}/E34_prereg.sha256").read().split()[0], "E34 protocol hash mismatch"
_s32 = open(f"{A}/E32_trials.py").read().split("# ------------------------------------------------------------------ Q1:")[0]
_s32 = _s32.replace("E32_step0_reproduction", "E34_step0_reproduction").replace("E32_coverage.csv", "E34_coverage_check.csv")
exec(compile(_s32, f"{A}/E32_trials.py", "exec"))
print("E34 protocol sha256 verified:", _h34[:12])
import openpyxl

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/validation/NeoPembrOV_GSE227666"
T_PANEL = ["CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "CD2", "CD247", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7"]
I_PANEL = ["IFNG", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCR5", "STAT1", "IDO1", "CD274", "PDCD1", "LAG3", "TIGIT",
           "CTLA4", "HAVCR2"]

# ------------------------------------------------------------------ NeoPembrOV data
import gzip
titles = [l.rstrip("\n").split("\t")[1:] for l in gzip.open(f"{D}/GSE227666_series_matrix.txt.gz", "rt") if l.startswith("!Sample_title")][0]
meta = pd.DataFrame({"title": [t.strip('"') for t in titles]})
meta["rid"] = meta.title.str.extract(r"\((R\d+)\)"); meta["arm"] = meta.title.str.extract(r", (NACT\+P|NACT), ")
meta["time"] = meta.title.str.extract(r", (Pre|Post) ")
link = pd.read_csv(f"{D}/id_match_final.csv")
link["margin"] = link.r - link.second
flag = set(link.loc[link.margin < 0.10, "geo"])          # 2 links (chemotherapy arm); same set as at the original run
wb = openpyxl.load_workbook(f"{D}/nc2024_source_data.xlsx", read_only=True)
h, *rs = list(wb["Figures 4 a,b,c,d"].iter_rows(values_only=True)); cl = pd.DataFrame(rs, columns=h).set_index("sample_name")
h, *rs = list(wb["Table 1"].iter_rows(values_only=True)); st = pd.DataFrame(rs, columns=h).set_index("sample_name")
pre = meta[meta.time == "Pre"].merge(link[["src", "geo"]], left_on="rid", right_on="geo", how="left")
assert pre.src.notna().all() and pre.src.is_unique, "unlinked or duplicated pre-treatment samples"
assert (pre.arm.values == cl.loc[pre.src, "Treatment_arm"].values).all(), "arm mismatch between GEO and source data"
Cn = pd.read_csv(f"{D}/GSE227666_neopembrov_counts.txt.gz", sep="\t", index_col=0); Cn = Cn.groupby(level=0).sum()
Cn = Cn[pre.rid]
Xn = np.log2(Cn / Cn.sum() * 1e6 + 0.5); Xn = Xn[(Cn > 0).sum(1) >= 3]
Zn = rowz(Xn); Sn = scores(Zn)
t_ = pd.to_numeric(cl.loc[pre.src, "del_pfs2"]).values.astype(float); e_ = pd.to_numeric(cl.loc[pre.src, "pfs2"]).values.astype(int)
arm_ = (pre.arm == "NACT+P").astype(int).values
npr = st.reindex(pre.src)["Patient_status"].map({"NPr": 1.0, "Pr": 0.0}).values
TR["NeoPembrOV"] = dict(kind="cox", Z=Zn, X=Xn, S=Sn.set_index(Zn.columns) if not Sn.index.equals(Zn.columns) else Sn, arm=arm_,
                        t=t_, e=e_, firth=True)
okr = np.isfinite(npr)
NPR = dict(kind="logit", Z=Zn.loc[:, okr], arm=arm_[okr].astype(float), y=npr[okr], cov=None, firth=True)
print(f"NeoPembrOV: {len(arm_)} pre-treatment ({arm_.sum()} NACT+P, {len(arm_) - arm_.sum()} NACT); PFS events {e_.sum()}; "
      f"24-month status for {okr.sum()}; genes {Zn.shape[0]}; flagged links {len(flag)}")

# ------------------------------------------------------------------ scores in NeoPembrOV
k = "NeoPembrOV"; nc = set(nonconst(Zn)); OB = {}
for s, v in LIB.items():
    pres = [g for g in v["genes"] if g in nc]
    if len(pres) >= 3 and len(pres) >= 0.5 * len(v["genes"]):
        w = v["w"][[v["genes"].index(g) for g in pres]]; sc = wscore(Zn, pres, w)
        if np.isfinite(sc).all():
            OB[("sig", s)] = sc
for s in EXTRA:
    if s in Sn and np.isfinite(Sn[s].to_numpy(float)).all():
        OB[("sig", s)] = Sn[s].to_numpy(float)
for g in T_PANEL + I_PANEL:
    if g in nc:
        OB[("gene", g)] = Zn.loc[g].to_numpy(float)
print(f"scored: {sum(1 for kk in OB if kk[0] == 'sig' and kk[1] in LIB)}/52 signatures; extras "
      f"{[s for s in EXTRA if ('sig', s) in OB]}; panel genes {sum(1 for kk in OB if kk[0] == 'gene')}/27")


def fit(task):
    key, variant = task
    sc = OB[key]; d = TR[k]
    try:
        if variant == "npr":
            return m_int(NPR, sc[okr])
        if variant == "dropflag":
            keep = np.where(~pre.rid.isin(flag).values)[0]
            return m_int(take(d, keep), sc[keep])
        return m_int(d, sc)
    except Exception:
        return (np.nan, np.nan)


tasks = [(key, v) for key in OB for v in ("primary", "npr", "dropflag")]
with mp.get_context("fork").Pool(NPROC) as pool:
    res = pool.map(fit, tasks, chunksize=4)
PN = pd.DataFrame([dict(kind=kk[0], name=kk[1], comparison=k, variant=v, y=r[0], se=r[1]) for (kk, v), r in zip(tasks, res)])
PN["p"] = 2 * norm.sf(np.abs(PN.y / PN.se)); PN.to_csv(f"{O}/E34_neopembrov_per_score.csv", index=False)

# ------------------------------------------------------------------ paired difference from CD8A in NeoPembrOV
rng_ = np.random.default_rng([3034, 7]); i1, i0 = np.where(arm_ == 1)[0], np.where(arm_ == 0)[0]
BI = [np.concatenate([rng_.choice(i1, len(i1)), rng_.choice(i0, len(i0))]) for _ in range(B_BOOT)]


def boot(key):
    sc = OB[key]; out = []
    for idx in BI:
        try:
            out.append(m_int(take(TR[k], idx), sc[idx])[0])
        except Exception:
            out.append(np.nan)
    return np.array(out)

bkeys = [kk for kk in OB if (kk[0] == "sig" and kk[1] in LIB) or kk == ("gene", "CD8A")]
with mp.get_context("fork").Pool(NPROC) as pool:
    bres = pool.map(boot, bkeys, chunksize=1)
BT = dict(zip(bkeys, bres))
Y = PN[PN.variant == "primary"].set_index(["kind", "name"]).y
DN = {}
for kk in bkeys:
    if kk[0] == "sig":
        db = BT[kk] - BT[("gene", "CD8A")]; ok = np.isfinite(db)
        DN[kk[1]] = (Y[kk] - Y[("gene", "CD8A")], np.sqrt(np.var(db[ok], ddof=1)))

# ------------------------------------------------------------------ pools
P32 = pd.read_csv(f"{O}/E32_per_comparison.csv"); P32 = P32[P32.variant == "primary"]
P33 = pd.read_csv(f"{O}/E33_per_comparison.csv"); P33 = P33[(P33.variant == "primary") & (P33.kind == "gene")]
PR = pd.read_csv(f"{O}/E33_paired.csv")


def per(name, kind, keys):
    src = P33[(P33.name == name)] if kind == "gene" else P32[P32.signature == name]
    r = src[src.comparison.isin(keys)][["y", "se"]]
    q = PN[(PN.kind == ("gene" if kind == "gene" else "sig")) & (PN.name == name) & (PN.variant == "primary")][["y", "se"]]
    return pd.concat([r, q])


ROWS = []
for setname, keys, kmin in [("seven (six + NeoPembrOV)", ALL6, 5), ("four chemotherapy-backbone (breast + NeoPembrOV)", BREAST, 4)]:
    for name, kind in [(s, "sig") for s in LIB] + [(s, "sig") for s in ("ESTIMATE", "TMEscore_z", "NetImmune", "Immune", "CAF_ECM")] + \
                      [(g, "gene") for g in T_PANEL + I_PANEL]:
        rr = per(name, kind, keys); rr = rr[np.isfinite(rr.y) & np.isfinite(rr.se)]
        if len(rr) < kmin:
            continue
        q = pool_rule(rr.y.values, rr.se.values)
        row = dict(set=setname, name=name, kind="library" if name in LIB else ("gene" if kind == "gene" else "extra"),
                   k=len(rr), est=q["est"], lo=q["lo"], hi=q["hi"], p=q["p"], k_pos=int((rr.y > 0).sum()))
        if name in LIB and name in DN:
            base = PR[(PR.set == ("six (primary)" if keys is ALL6 else "breast three")) & (PR.reference == "CD8A") & (PR.signature == name)]
            if len(base):
                b = base.iloc[0]; d7, s7 = DN[name]
                w1, w2 = 1 / b.se ** 2, 1 / s7 ** 2; Dd = (w1 * b["diff"] + w2 * d7) / (w1 + w2); se = np.sqrt(1 / (w1 + w2))
                row.update(diff_vs_CD8A=Dd, cd8a_lo=Dd - 1.96 * se, cd8a_hi=Dd + 1.96 * se, cd8a_p=2 * norm.sf(abs(Dd / se)))
        ROWS.append(row)
PO = pd.DataFrame(ROWS); PO.to_csv(f"{O}/E34_pooled.csv", index=False)

# ------------------------------------------------------------------ expectations
L7 = PO[(PO.set.str.startswith("seven")) & (PO.kind == "library")]; L4 = PO[(PO.set.str.startswith("four")) & (PO.kind == "library")]
f7, a7 = int(((L7.p < 0.05) & (L7.est > 0)).sum()), int(((L7.p < 0.05) & (L7.est < 0)).sum())
ab7 = int(((L7.diff_vs_CD8A > 0) & (L7.cd8a_p < 0.05)).sum())
f4, a4 = int(((L4.p < 0.05) & (L4.est > 0)).sum()), int(((L4.p < 0.05) & (L4.est < 0)).sum())
ab4 = int(((L4.diff_vs_CD8A > 0) & (L4.cd8a_p < 0.05)).sum())
EXP = pd.DataFrame([
    dict(q="Q1", expectation="seven-comparison counts within ±3 of 17 favoring / 1 against", result=f"{f7} favoring / {a7} against",
         held=abs(f7 - 17) <= 3 and abs(a7 - 1) <= 3),
    dict(q="Q2", expectation="at most 2 of 52 significantly above CD8A over seven", result=f"{ab7}/52", held=ab7 <= 2),
    dict(q="Q3", expectation="fewer than 10 favoring in the four chemotherapy-backbone comparisons",
         result=f"{f4} favoring / {a4} against; above CD8A {ab4}/52", held=f4 < 10),
    dict(q="Q4", expectation="descriptive", result="see E34_neopembrov_per_score.csv", held=None)])
EXP.to_csv(f"{O}/E34_expectations.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120)
print("\n== expectations\n", EXP.to_string(index=False))
q4 = PN[(PN.name.isin(["CD8A", "NetImmune", "Immune", "CAF_ECM", "ESTIMATE"] + T_PANEL)) & (PN.variant == "primary")]
print("\n== NeoPembrOV alone (primary, PFS)\n", q4[["kind", "name", "y", "se", "p"]].round(3).to_string(index=False))
nl = PN[(PN.kind == "sig") & PN.name.isin(LIB) & (PN.variant == "primary")]
print(f"\nNeoPembrOV library: {int(((nl.p < 0.05) & (nl.y > 0)).sum())} favoring, {int(((nl.p < 0.05) & (nl.y < 0)).sum())} against at P < 0.05 of {len(nl)}")
for nm in ["CD8A", "NetImmune", "CAF_ECM"]:
    for v in ("npr", "dropflag"):
        r = PN[(PN.name == nm) & (PN.variant == v)]
        if len(r):
            print(f"  {nm} {v}: {r.y.iloc[0]:+.3f} (SE {r.se.iloc[0]:.3f}, P = {r.p.iloc[0]:.3g})")
c7 = PO[(PO.name == "CD8A") & (PO.kind == "gene")]
print("\nCD8A pools:\n", c7[["set", "k", "est", "lo", "hi", "p", "k_pos"]].round(3).to_string(index=False))
