"""E26 - an independent gastric Visium atlas (GSE270678; 72 sections, 20 donors), pre-registered in
E26_prereg.md (2026-10-19) before any expression value, coordinate or image was read.

Primary  P1  SPP1_arm Ring1 - Far < 0 in more donors than not (one-sided sign test)
         P2  CAF_ECM Ring1 - Far not consistently > 0 (one-sided sign test for > 0, P >= 0.05)
         P3  BarrierIndex not associated with stage (Spearman over sections, P >= 0.05)
Second.  S1/S2  section T-cell infiltration of tumor ~ fibroblast number + continuity excess + tumor
         fraction, donor random intercept

Rings, peaks and the BarrierIndex are barriernet 0.2.1 (assign_rings / ring_summary / barrier_index) with
their defaults; module scores are the bulk rule (mean per-gene z over spots, >= 3 genes) applied per section.
"""
import warnings, functools, gzip, glob, os, re, sys; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
import scanpy as sc
from scipy.stats import binomtest, spearmanr
import statsmodels.formula.api as smf
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/repo_barriernet/barriernet")
from barriernet.spatial import assign_rings, ring_summary, barrier_index, build_hex_graph, infer_spatial_keys

G = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/geo/GSE270678"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
MIN_UMI, MIN_SPOTS, NPERM, MIN_BOUNDARY = 500, 300, 200, 30
MOD = dict(
    CAF_ECM=["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
    TAM_SUP=["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
    APM=["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    IFN=["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
    TLS_B=["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
    T_CELL=["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"],
    C1Q_arm=["C1QA", "C1QB", "C1QC"],
    SPP1_arm=["SPP1", "TREM2", "MARCO", "APOE", "LGALS3"])
EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4", "CDH1"]
RING_SCORES = ["SPP1_arm", "CAF_ECM", "TAM_SUP", "C1Q_arm", "IFN", "T_CELL"]

# ------------------------------------------------------------------ metadata (series matrix)
with gzip.open(f"{G}/series_matrix.txt.gz", "rt", errors="ignore") as f:
    lines = [l.rstrip("\n").split("\t") for l in f if l.startswith("!Sample_")]
acc = [v.strip('"') for v in next(l for l in lines if l[0] == "!Sample_geo_accession")[1:]]
META = pd.DataFrame({"gsm": acc}).set_index("gsm")
META["title"] = [v.strip('"') for v in next(l for l in lines if l[0] == "!Sample_title")[1:]]
for l in lines:
    if l[0] == "!Sample_characteristics_ch1":
        vals = [v.strip('"') for v in l[1:]]
        key = vals[0].split(":")[0].strip()
        META[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
META = META.rename(columns={"subject id": "donor", "Stage": "stage", "subtype": "subtype", "sampling location": "location"})
print(f"metadata: {len(META)} sections, {META.donor.nunique()} donors")


# ------------------------------------------------------------------ one section
def read_section(gsm):
    fs = glob.glob(f"{G}/{gsm}_*")
    h5 = [f for f in fs if f.endswith("filtered_feature_bc_matrix.h5")][0]
    pos = [f for f in fs if "tissue_positions" in f][0]
    ad = sc.read_10x_h5(h5); ad.var_names_make_unique()
    with gzip.open(pos, "rt") as f:
        first = f.readline()
    hdr = 0 if "barcode" in first else None
    P = pd.read_csv(pos, header=hdr)
    P.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    P = P.set_index("barcode")
    keep = [b for b in ad.obs_names if b in P.index and P.loc[b, "in_tissue"] == 1]
    ad = ad[keep].copy()
    ad.obs = ad.obs.join(P[["array_row", "array_col", "pxl_row", "pxl_col"]])
    ad.obs["umi"] = np.asarray(ad.X.sum(1)).ravel()
    ad = ad[ad.obs.umi >= MIN_UMI].copy()
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    return ad


def module_scores(ad):
    need = [g for g in dict.fromkeys(sum(MOD.values(), []) + EPI) if g in ad.var_names]
    sub = ad[:, need]
    X = pd.DataFrame(sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X, index=ad.obs_names, columns=need)
    absS, zS, used = {}, {}, {}
    for k, gs in list(MOD.items()) + [("EPI", EPI)]:
        g = [x for x in gs if x in X.columns]
        sd = X[g].std(0)
        gz = [x for x in g if sd[x] > 0]
        absS[k] = X[g].mean(1) if len(g) >= 3 else pd.Series(np.nan, index=X.index)
        Z = (X[gz] - X[gz].mean(0)) / X[gz].std(0)
        zS[k] = Z.mean(1) if len(gz) >= 3 else pd.Series(np.nan, index=X.index)
        used[k] = len(gz)
    return pd.DataFrame(absS), pd.DataFrame(zS), used


SEC = {}
for gsm in META.index:
    if not glob.glob(f"{G}/{gsm}_*filtered_feature_bc_matrix.h5"):
        print(f"  {gsm}: files missing, skipped"); continue
    ad = read_section(gsm)
    if ad.n_obs < MIN_SPOTS:
        print(f"  {gsm}: {ad.n_obs} spots after QC, excluded"); continue
    A, Zs, used = module_scores(ad)
    SEC[gsm] = dict(ad=ad, A=A, Z=Zs, used=used)
print(f"sections analyzed: {len(SEC)} of {len(META)}; spots {sum(v['ad'].n_obs for v in SEC.values()):,}")

# global thresholds for the secondary analysis (pooled over all spots of all sections)
allA = pd.concat([v["A"] for v in SEC.values()])
EPI_THR, CAF_THR = allA.EPI.median(), allA.CAF_ECM.quantile(0.75)
print(f"tumor-spot threshold (epithelial, median) {EPI_THR:.3f}; fibroblast-high threshold (CAF_ECM, 75th pct) {CAF_THR:.3f}")

rng = np.random.default_rng(2026)
rows = []
for gsm, v in SEC.items():
    ad, A, Zs = v["ad"], v["A"], v["Z"]
    obs = ad.obs.copy()
    obs["TLS_score"] = Zs.TLS_B.values
    obs["score_CAF_ECM"] = Zs.CAF_ECM.values
    obs["score_TAM_SUP"] = Zs.TAM_SUP.values
    xy = obs[["pxl_col", "pxl_row"]].to_numpy(float)
    R = assign_rings(obs, xy, "TLS_score")                    # package defaults
    keys = infer_spatial_keys(obs)
    rs = ring_summary(obs, R, keys)
    r = dict(gsm=gsm, donor=META.loc[gsm, "donor"], stage=META.loc[gsm, "stage"], subtype=META.loc[gsm, "subtype"],
             location=META.loc[gsm, "location"], spots=ad.n_obs, tls_peaks=int(R.tls_peak.sum()),
             ring1=int((R.ring == "Ring1").sum()), far=int((R.ring == "Far").sum()),
             tam_key=keys.tam, ecm_key=keys.ecm,
             ECM_Ring1_minus_Far=rs.ECM_Ring1_minus_Far.iloc[0], TAM_Ring1_minus_Far=rs.TAM_Ring1_minus_Far.iloc[0])
    for s in RING_SCORES:
        z = Zs[s].values
        r[f"{s}_R1mF"] = np.nanmean(z[R.ring.values == "Ring1"]) - np.nanmean(z[R.ring.values == "Far"])
    # secondary: number versus arrangement
    nbrs, keep = build_hex_graph(obs)
    tumor = (A.EPI.values >= EPI_THR); fib = (A.CAF_ECM.values >= CAF_THR)
    nontumor = ~tumor
    bnd = np.array([tumor[i] and any(nontumor[j] for j in nbrs[i]) for i in range(len(obs))])
    def cont(fmask):
        b = np.where(bnd)[0]
        return np.mean([any(fmask[j] for j in nbrs[i]) for i in b]) if len(b) else np.nan
    obs_c = cont(fib)
    nt_idx = np.where(nontumor)[0]; nfib_nt = int((fib & nontumor).sum())
    perm = []
    for _ in range(NPERM):
        f2 = np.zeros(len(obs), bool); f2[rng.choice(nt_idx, nfib_nt, replace=False)] = True
        perm.append(cont(f2))
    r.update(boundary=int(bnd.sum()), fib_fraction=fib.mean(), tumor_fraction=tumor.mean(),
             continuity=obs_c, continuity_expected=np.mean(perm), continuity_excess=obs_c - np.mean(perm),
             tcell_tumor=np.nanmean(A.T_CELL.values[tumor]) if tumor.any() else np.nan,
             **{f"genes_{k}": v["used"][k] for k in ["CAF_ECM", "TAM_SUP", "TLS_B", "SPP1_arm", "T_CELL", "EPI"]})
    rows.append(r)
T = pd.DataFrame(rows)
T["BarrierIndex"] = barrier_index(T.set_index("gsm")).values
T.to_csv(f"{O}/E26_sections.csv", index=False)

# ------------------------------------------------------------------ primary tests (donor level)
D = T.groupby("donor")[[f"{s}_R1mF" for s in RING_SCORES]].mean()
D.to_csv(f"{O}/E26_donor.csv")
tests = []
for s in RING_SCORES:
    x = D[f"{s}_R1mF"].dropna(); neg, pos = int((x < 0).sum()), int((x > 0).sum())
    tests.append(dict(test=f"{s} Ring1 - Far (donors)", n=len(x), negative=neg, positive=pos,
                      median=x.median(), p_less=binomtest(neg, neg + pos, 0.5, alternative="greater").pvalue,
                      p_greater=binomtest(pos, neg + pos, 0.5, alternative="greater").pvalue))
st = T.stage.map({"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4})
rho, p = spearmanr(T.BarrierIndex, st, nan_policy="omit")
tests.append(dict(test="BarrierIndex vs stage (sections, Spearman)", n=int(st.notna().sum()), median=rho, p_less=p, p_greater=np.nan))
TT = pd.DataFrame(tests); TT.to_csv(f"{O}/E26_tests.csv", index=False)
print(TT.round(4).to_string(index=False))
p1 = TT.iloc[0]; p2 = TT.iloc[1]
print(f"\nP1 SPP1_arm Ring1 < Far: {p1.negative}/{p1.negative + p1.positive} donors, one-sided P = {p1.p_less:.4f} -> "
      f"{'held' if p1.p_less < 0.05 else 'not held'}")
print(f"P2 CAF_ECM Ring1 > Far: {p2.positive}/{p2.negative + p2.positive} donors, one-sided P = {p2.p_greater:.4f} -> "
      f"{'held' if p2.p_greater >= 0.05 else 'not held'}")
print(f"P3 BarrierIndex vs stage: rho {rho:.2f}, P = {p:.3f} -> {'held' if p >= 0.05 else 'not held'}")

# ------------------------------------------------------------------ secondary: number versus arrangement
M = T[T.boundary >= MIN_BOUNDARY].copy()
for c in ["fib_fraction", "continuity_excess", "tumor_fraction"]:
    M[f"z_{c}"] = (M[c] - M[c].mean()) / M[c].std()
fit = smf.mixedlm("tcell_tumor ~ z_fib_fraction + z_continuity_excess + z_tumor_fraction", M, groups=M["donor"]).fit(reml=True)
with open(f"{O}/E26_mixedlm.txt", "w") as f:
    f.write(str(fit.summary()))
print(f"\nsecondary model: {len(M)} sections, {M.donor.nunique()} donors")
print(fit.summary().tables[1])
b = fit.params; pv = fit.pvalues
print(f"S1 fibroblast number: {b['z_fib_fraction']:+.4f}, P = {pv['z_fib_fraction']:.4f} -> "
      f"{'held' if (b['z_fib_fraction'] < 0 and pv['z_fib_fraction'] < 0.05) else 'not held'}")
print(f"S2 continuity excess: {b['z_continuity_excess']:+.4f}, P = {pv['z_continuity_excess']:.4f} -> "
      f"{'held' if not (b['z_continuity_excess'] < 0 and pv['z_continuity_excess'] < 0.05) else 'not held'}")
print("\nE26 done")
