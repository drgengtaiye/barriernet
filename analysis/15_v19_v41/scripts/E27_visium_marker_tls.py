"""E27 - TLS-centered Visium rings recomputed with a marker-based TLS definition, in the published gastric
sections (GSE251950, 10) and the independent atlas (GSE270678, 69 sections passing E26 QC, 20 donors).

Why: E26b showed that the TLS_niche_score behind the published rings (Tangram B cell + plasma cell + DC
fractions) is unrelated to lymphoid markers at spot level (MS4A1, CD79A, CXCL13, CR2, CCL19, CCL21 all |r| <
0.1), so the published rings were not centered on lymphoid aggregates. Here TLS peaks are called on the
barriernet TLS_B module (CXCL13, LTB, MS4A1, CD79A, CD74, MZB1), validated against CR2, CCL19 and CCL21,
which are not in it.

Identical in both atlases: module score = mean per-gene z over the section's in-tissue spots (>= 3 genes);
rings = barriernet 0.2.1 assign_rings defaults (kNN-9 local maxima above the 95th percentile, Far >= 3 hops);
BarrierIndex = z(CAF_ECM Ring1 - Far) + z(SPP1_arm Ring1 - Far), standardized within the atlas (population SD),
the published definition with the SPP1+ macrophage module in place of the Tangram-era state.
Tests: published sections, two-sided paired Wilcoxon over the 10 sections (as published); GSE270678, donor
means (20 donors), two-sided Wilcoxon and sign test. Sensitivity (published sections): lattice peaks, and Far
>= 4 hops with tissue-edge spots excluded.
"""
import warnings, glob, gzip, sys; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
from scipy.stats import wilcoxon, binomtest, spearmanr
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/repo_barriernet/barriernet")
from barriernet.spatial import assign_rings

O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
V = "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/results/ring_barrier_state_v2"
G = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/geo/GSE270678"
MIN_UMI, MIN_SPOTS = 500, 300                                     # as E26
MOD = dict(
    CAF_ECM=["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
    TAM_SUP=["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
    APM=["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    IFN=["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
    TLS_B=["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
    T_CELL=["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"],
    C1Q_arm=["C1QA", "C1QB", "C1QC"],
    SPP1_arm=["SPP1", "TREM2", "MARCO", "APOE", "LGALS3"])
SCORES = ["SPP1_arm", "CAF_ECM", "T_CELL", "TAM_SUP", "C1Q_arm", "IFN", "APM", "TLS_B"]
STATES = ["state_SPP1_TAM", "state_CAF_ECM", "state_IL10_STAT3"]         # stored Tangram-era states, published sections
VALID = ["CR2", "CCL19", "CCL21", "MS4A1", "CD79A", "CXCL13"]
ZONES = ["Hotspot", "Ring1", "Ring2", "Far"]
DEFS = [("knn9", "far3", dict()), ("lattice", "far3", dict(peak_method="lattice")),
        ("knn9", "far4_noedge", dict(far_min_dist=4, far_exclude_edge=True))]


def zscores(X):
    out = {}
    for k, gs in MOD.items():
        g = [x for x in gs if x in X.columns and X[x].std() > 0]
        out[k] = ((X[g] - X[g].mean()) / X[g].std()).mean(1) if len(g) >= 3 else pd.Series(np.nan, index=X.index)
    return pd.DataFrame(out)


def contrasts(o, ring, cols):
    r = {}
    for c in cols:
        v = o[c].to_numpy(float)
        far = np.nanmean(v[ring == "Far"])
        for z in ["Ring1", "Ring2"]:
            r[f"{c}_{z}_minus_Far"] = np.nanmean(v[ring == z]) - far
    return r


def tls_validation(o, X, atlas, section):
    r = dict(atlas=atlas, section=section)
    for g in VALID:
        if g in X.columns and X[g].std() > 0:
            r[f"TLS_B|{g}"] = np.corrcoef(o.TLS_B, X[g])[0, 1]
            if "TLS_niche_score" in o:
                r[f"TLS_niche|{g}"] = np.corrcoef(o.TLS_niche_score.astype(float), X[g])[0, 1]
    return r


SEC, ZM, VAL, SPOTS = [], [], [], []
need = list(dict.fromkeys(sum(MOD.values(), []) + VALID))

# ------------------------------------------------------------------ published sections (GSE251950)
clin = pd.read_csv(f"{V}/BarrierIndex_by_slice_with_clinical.csv")
GC = dict(zip(clin.slice, clin.GC_raw)); MET = dict(zip(clin.slice, clin.is_metastasis))
for f in sorted(glob.glob(f"{V}/*.with_state_and_ring.h5ad")):
    s = f.split("/")[-1].split(".")[0]
    a = ad.read_h5ad(f)
    a = a[a.obs.in_tissue.astype(int) == 1].copy()
    g = [x for x in need if x in a.var_names]
    X = pd.DataFrame(a[:, g].X.toarray() if hasattr(a.X, "toarray") else a[:, g].X, index=a.obs_names, columns=g)
    o = a.obs.copy(); o = pd.concat([o, zscores(X)], axis=1)
    xy = np.asarray(a.obsm["spatial"], float)
    VAL.append(tls_validation(o, X, "GSE251950", GC[s]))
    for pk, fr, kw in DEFS:
        R = assign_rings(o, xy, "TLS_B", **kw)
        ring = R.ring.values
        r = dict(atlas="GSE251950", section=s, unit=GC[s], GC_raw=GC[s], is_metastasis=int(MET[s]), peaks=pk, far=fr,
                 spots=a.n_obs, n_peaks=int(R.tls_peak.sum()), n_Ring1=int((ring == "Ring1").sum()),
                 n_Ring2=int((ring == "Ring2").sum()), n_Far=int((ring == "Far").sum()))
        r.update(contrasts(o, ring, SCORES + STATES))
        SEC.append(r)
        if pk == "knn9" and fr == "far3":
            for c in SCORES + STATES:
                for z in ZONES:
                    ZM.append(dict(atlas="GSE251950", unit=GC[s], score=c, zone=z,
                                   mean=float(np.nanmean(o[c].to_numpy(float)[ring == z])), n=int((ring == z).sum())))
            sp = pd.DataFrame({"slice": s, "GC": GC[s], "x": xy[:, 0], "y": -xy[:, 1], "ring": ring})
            for c in ["SPP1_arm", "CAF_ECM", "T_CELL", "TLS_B"]:
                sp[c] = o[c].to_numpy(float)
            SPOTS.append(sp)
print("published sections:", sum(1 for r in SEC if r["atlas"] == "GSE251950" and r["peaks"] == "knn9" and r["far"] == "far3"))

# ------------------------------------------------------------------ independent atlas (GSE270678), as E26
with gzip.open(f"{G}/series_matrix.txt.gz", "rt", errors="ignore") as fh:
    lines = [l.rstrip("\n").split("\t") for l in fh if l.startswith("!Sample_")]
acc = [v.strip('"') for v in next(l for l in lines if l[0] == "!Sample_geo_accession")[1:]]
META = pd.DataFrame({"gsm": acc}).set_index("gsm")
for l in lines:
    if l[0] == "!Sample_characteristics_ch1":
        vals = [v.strip('"') for v in l[1:]]
        META[vals[0].split(":")[0].strip()] = [v.split(":", 1)[1].strip() for v in vals]
META = META.rename(columns={"subject id": "donor", "Stage": "stage", "sampling location": "location"})
for gsm in META.index:
    fs = glob.glob(f"{G}/{gsm}_*")
    h5 = [x for x in fs if x.endswith("filtered_feature_bc_matrix.h5")][0]
    pos = [x for x in fs if "tissue_positions" in x][0]
    a = sc.read_10x_h5(h5); a.var_names_make_unique()
    with gzip.open(pos, "rt") as fh:
        hdr = 0 if "barcode" in fh.readline() else None
    P = pd.read_csv(pos, header=hdr); P.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
    P = P.set_index("barcode")
    a = a[[b for b in a.obs_names if b in P.index and P.loc[b, "in_tissue"] == 1]].copy()
    a.obs = a.obs.join(P[["array_row", "array_col", "pxl_row", "pxl_col"]])
    a = a[np.asarray(a.X.sum(1)).ravel() >= MIN_UMI].copy()
    if a.n_obs < MIN_SPOTS:
        continue
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    g = [x for x in need if x in a.var_names]
    X = pd.DataFrame(a[:, g].X.toarray(), index=a.obs_names, columns=g)
    o = pd.concat([a.obs.copy(), zscores(X)], axis=1)
    xy = o[["pxl_col", "pxl_row"]].to_numpy(float)
    VAL.append(tls_validation(o, X, "GSE270678", gsm))
    R = assign_rings(o, xy, "TLS_B"); ring = R.ring.values
    r = dict(atlas="GSE270678", section=gsm, unit=META.loc[gsm, "donor"], donor=META.loc[gsm, "donor"],
             stage=META.loc[gsm, "stage"], subtype=META.loc[gsm, "subtype"], location=META.loc[gsm, "location"],
             is_metastasis=0, peaks="knn9", far="far3", spots=a.n_obs, n_peaks=int(R.tls_peak.sum()),
             n_Ring1=int((ring == "Ring1").sum()), n_Ring2=int((ring == "Ring2").sum()), n_Far=int((ring == "Far").sum()))
    r.update(contrasts(o, ring, SCORES)); SEC.append(r)
    for c in SCORES:
        for z in ZONES:
            ZM.append(dict(atlas="GSE270678", unit=META.loc[gsm, "donor"], section=gsm, score=c, zone=z,
                           mean=float(np.nanmean(o[c].to_numpy(float)[ring == z])), n=int((ring == z).sum())))
T = pd.DataFrame(SEC)
print("independent atlas sections:", int((T.atlas == "GSE270678").sum()), "donors:", T[T.atlas == "GSE270678"].donor.nunique())

# BarrierIndex within each atlas (primary ring definition)
z0 = lambda v: (v - v.mean()) / v.std(ddof=0)
T["z_ecm"] = np.nan; T["z_tam"] = np.nan; T["BarrierIndex"] = np.nan
for atl in ["GSE251950", "GSE270678"]:
    k = (T.atlas == atl) & (T.peaks == "knn9") & (T.far == "far3")
    T.loc[k, "z_ecm"] = z0(T.loc[k, "CAF_ECM_Ring1_minus_Far"]); T.loc[k, "z_tam"] = z0(T.loc[k, "SPP1_arm_Ring1_minus_Far"])
    T.loc[k, "BarrierIndex"] = T.loc[k, "z_ecm"] + T.loc[k, "z_tam"]
T.to_csv(f"{O}/E27_rings_sections.csv", index=False)

# zone means: donor means for the independent atlas
Z = pd.DataFrame(ZM)
Z = pd.concat([Z[Z.atlas == "GSE251950"],
               Z[Z.atlas == "GSE270678"].groupby(["atlas", "unit", "score", "zone"], as_index=False)[["mean", "n"]].mean()])
Z.to_csv(f"{O}/E27_zone_means.csv", index=False)
pd.concat(SPOTS).to_parquet(f"{O}/E27_spots_GSE251950.parquet")
VV = pd.DataFrame(VAL); VV.to_csv(f"{O}/E27_tls_validation.csv", index=False)

# ------------------------------------------------------------------ tests
rows = []
for (atl, pk, fr), g in T.groupby(["atlas", "peaks", "far"]):
    if atl == "GSE270678":
        g = g.groupby("donor")[[c for c in g.columns if c.endswith("_minus_Far")]].mean()
    for c in [c for c in g.columns if c.endswith("_minus_Far")]:
        v = g[c].dropna().to_numpy()
        if len(v) < 5:
            continue
        neg, pos = int((v < 0).sum()), int((v > 0).sum())
        rows.append(dict(atlas=atl, peaks=pk, far=fr, score=c.rsplit("_Ring", 1)[0], contrast="Ring" + c.rsplit("_Ring", 1)[1].replace("_minus_Far", "-Far"),
                         n=len(v), n_below=neg, n_above=pos, median=float(np.median(v)),
                         P=float(wilcoxon(v).pvalue), P_sign=float(binomtest(neg, neg + pos, 0.5).pvalue)))
TS = pd.DataFrame(rows); TS.to_csv(f"{O}/E27_rings_summary.csv", index=False)
show = TS[(TS.peaks == "knn9") & (TS.far == "far3") & (TS.contrast == "Ring1-Far") & TS.score.isin(["T_CELL", "SPP1_arm", "CAF_ECM", "TAM_SUP", "state_SPP1_TAM", "state_IL10_STAT3", "state_CAF_ECM"])]
print(show[["atlas", "score", "n", "n_below", "n_above", "median", "P", "P_sign"]].round(4).to_string(index=False))

# BarrierIndex against clinical variables
p = T[(T.atlas == "GSE251950") & (T.peaks == "knn9") & (T.far == "far3")].merge(
    clin[["slice", "Age", "pT", "pN"]], left_on="section", right_on="slice")
p = p[p.is_metastasis == 0]
p["pT_o"] = p.pT.str.extract(r"T(\d)").astype(float); p["pN_o"] = p.pN.str.extract(r"N(\d)").astype(float)
bi = []
for v in ["Age", "pT_o", "pN_o"]:
    r_, p_ = spearmanr(p.BarrierIndex, p[v]); bi.append(dict(atlas="GSE251950", variable=v, n=len(p), rho=r_, P=p_))
q = T[T.atlas == "GSE270678"].groupby("donor").agg(BarrierIndex=("BarrierIndex", "mean"), stage=("stage", "first"))
r_, p_ = spearmanr(q.BarrierIndex, q.stage.map({"Stage I": 1, "Stage II": 2, "Stage III": 3}))
bi.append(dict(atlas="GSE270678", variable="stage (donor level)", n=len(q), rho=r_, P=p_))
BI = pd.DataFrame(bi); BI.to_csv(f"{O}/E27_barrierindex_clinical.csv", index=False)
print(); print(BI.round(3).to_string(index=False))
print("\nTLS definition against lymphoid markers (median spot-level r):")
print(VV.groupby("atlas").median(numeric_only=True).round(2).T.to_string())
print("\nE27 done")
