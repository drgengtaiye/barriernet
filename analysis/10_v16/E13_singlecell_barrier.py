"""E13 - Single-cell-resolution tests of the physical-barrier idea (public IMC; not gastric).

HYPOTHESES (fixed before the first run)
  H1 (physical barrier): at fixed stromal abundance, a more continuous stromal/fibroblast
     layer at the tumour-stroma interface is accompanied by fewer CD8 T cells inside
     tumour nests (negative continuity coefficient in
     infiltration ~ continuity + stromal density + total CD8 density).
  H0: CD8 infiltration tracks immune abundance and stromal composition, not continuity
     (continuity coefficient ~ 0 once stromal and CD8 density are in the model).
  Secondary (replication of the Visium SPP1+ TAM / TLS result): macrophage density is
     lower next to B/T-cell aggregates than further away.
  Tertiary (breast only): stromal continuity or stromal density relates to
     disease-specific survival in METABRIC.

DATA (public; neither dataset is gastric)
  PDAC IMC   spatial_public/pdac_cells.csv, exported from the authors' Seurat object
             (Zenodo record 10246315). 34 ROIs from 9 patients.
  Breast IMC public_imc/SingleCells.csv, METABRIC imaging mass cytometry
             (Zenodo record 5850952), linked by metabric_id to METABRIC clinical data
             (cBioPortal brca_metabric).
  Not used   codex_crc_B004.csv is healthy small bowel and colon (tissue SB/CL, donor
             B004), so it has no tumour boundary. Basel/Zurich IMC metaclusters do not
             separate CD8 from other T cells, so METABRIC IMC is the breast cohort.
             The Xenium breast section is a single sample and was not processed.

DESIGN
  Cells are grouped as tumour/epithelium (T), fibroblast stroma (S), endothelium (V),
  other non-immune (O) and immune (CD8 T, other T, B, other lymphocytes, macrophages,
  granulocytes/other myeloid).
  - The "scaffold" is the set of non-immune cells. Tumour regions, the boundary and the
    continuity metric are computed on the scaffold only, so where the CD8 cells sit
    cannot feed back into the barrier measure.
  - Tumour region: a scaffold cell whose 10 nearest scaffold cells (self included) are
    >= 50% tumour.
  - Graph: Delaunay triangulation of the scaffold, edges <= 40 um.
  - Boundary edges join a tumour-region cell to a non-tumour-region cell. The outer end
    of each boundary edge is a "rind" cell.
  - Continuity: the fraction of boundary edges whose rind cell is a fibroblast lying in
    a run of >= 5 rind fibroblasts connected along the rind. This is the perimeter
    covered by a contiguous fibroblast layer.
  - Abundance-preserving null: fibroblast labels are permuted among the non-tumour-region
    scaffold cells (100 draws). This keeps the fibroblast count and the tumour geometry
    fixed and moves only where the fibroblasts sit. Excess continuity = observed - null
    mean.
  - Area: tissue is rasterised at 5 um (pixels within 15 um of a cell); each pixel takes
    the region of its nearest scaffold cell.
  - Infiltration: log2 of [CD8 per mm2 in tumour regions] / [CD8 per mm2 elsewhere],
    with +0.5 cell in each. Stromal abundance is fibroblasts per mm2 of tissue (log10);
    CD8 abundance is CD8 per mm2 of tissue (log10). All are densities, not fractions
    of a total.
  - Model: infiltration ~ z(continuity) + z(log stromal density) + z(log CD8 density).
      PDAC: ROI level, with a patient random intercept (MixedLM) and a patient cluster
            bootstrap; patient-pooled OLS as a sensitivity analysis.
      METABRIC: counts and areas pooled per patient, OLS with HC3 intervals.
  - Positive control of the whole pipeline: CD8-sized sets of random walkers start in
    the non-tumour region of each image and take 100 steps on the scaffold graph.
    Crossings between a rind fibroblast and the tumour region are blocked with
    probability p in {0, 0.5, 1}. The same model is then fitted to the simulated
    infiltration (5 replicates). With p > 0 the continuity coefficient should turn
    negative; p = 0 shows what tissue geometry alone produces.
  - Aggregates: lymphocytes (CD8, other T, B, plasma/CD38) with >= 8 lymphocyte
    neighbours within 30 um, linked at <= 20 um, in components of >= 25 cells that
    contain >= 3 B cells.
      Near = non-member cells or pixels <= 50 um from an aggregate cell.
      Far  = > 150 um.
      Densities are pooled per patient and summarised as log2(near/far). Fibroblast
      and all-scaffold densities are comparators for displacement. Stroma-region-only
      contrasts are a sensitivity analysis.
Outputs: analysis/results/E13_*
"""
import warnings, functools, json, time, sys
warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from scipy.spatial import cKDTree, Delaunay
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy import stats

R = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
SP = "/groups/adv2105_gp/yichen/Yi/com/gastric/spatial_public"
PI = "/groups/adv2105_gp/yichen/Yi/com/gastric/public_imc"

LMAX, KREG, RUN, NPERM, PIX, TISR = 40.0, 10, 5, 100, 5.0, 15.0
RUNS = (3, 5, 10)          # run-length sensitivity; 5 is the primary definition
PBLOCK, NREP, NSTEP = (0.0, 0.5, 1.0), 5, 100
NEAR, FAR = 50.0, 150.0

# class codes
T, S, V, O, CD8, TC, B, LO, M, G, IO = range(11)
SCAF = (T, S, V, O)
LYMPH = (CD8, TC, B, LO)


# ------------------------------------------------------------------ loaders
def load_pdac():
    d = pd.read_csv(f"{SP}/pdac_cells.csv")
    d = d[d.type != "Artifact"]
    m = {"Tumor/Ductal Cells": T, "Proliferating Tumor/Ductal Cells": T,
         "Stromal/ECM-1": S, "Stromal/ECM-2": S, "Endothelial/Vascular Cells": V,
         "Alpha Cells": O, "Beta Cells": O, "Delta Cells": O,
         "CD8+ T Cells": CD8, "CD4+ T Cells": TC, "B Cells": B,
         "CD68+ Macrophages": M, "Proliferating Macrophages": M,
         "CD68+ CD44+ Macrophages": M, "CD11b+ Granulocytes": G}
    return pd.DataFrame({"dataset": "PDAC", "image": d.roi.values,
                         "patient": d.roi.str.split("-").str[0].values,
                         "x": d.x.values, "y": d.y.values,
                         "cls": d.type.map(m).astype(int).values, "cd163": np.nan})


def load_metabric():
    cols = ["ImageNumber", "metabric_id", "cellPhenotype", "is_epithelial", "is_tumour",
            "CD163", "Location_Center_X", "Location_Center_Y"]
    d = pd.read_csv(f"{PI}/SingleCells.csv", usecols=cols)
    ph = {"CD8^{+} T cells": CD8, "CD4^{+} T cells": TC, "CD4^{+} T cells & APCs": TC,
          "T_{Reg} & T_{Ex}": TC, "B cells": B, "CD38^{+} lymphocytes": LO,
          "Macrophages": M, "Macrophages & granulocytes": G, "Granulocytes": G,
          "CD15^{+}": G, "MHC^{hi}CD15^{+}": G, "MHC I & II^{hi}": IO, "CD57^{+}": IO,
          "MHC I^{hi}CD57^{+}": IO, "Fibroblasts": S, "Myofibroblasts": S,
          "Myofibroblasts PDPN^{+}": S, "Fibroblasts FSP1^{+}": S, "Endothelial": V}
    c = d.cellPhenotype.map(ph)
    c = c.where(c.notna(), np.where(d.is_epithelial == 1, T, O)).astype(int)
    tum = d.groupby("ImageNumber").is_tumour.transform("mean") > 0.5
    out = pd.DataFrame({"dataset": "METABRIC", "image": "MB_img" + d.ImageNumber.astype(str),
                        "patient": d.metabric_id.values, "x": d.Location_Center_X.values,
                        "y": d.Location_Center_Y.values, "cls": c.values,
                        "cd163": d.CD163.values})
    return out[tum.values & d.metabric_id.notna().values].reset_index(drop=True)


# ------------------------------------------------------------------ geometry
def delaunay_edges(xy, lmax=LMAX):
    s = Delaunay(xy).simplices
    e = np.vstack([s[:, [0, 1]], s[:, [1, 2]], s[:, [0, 2]]])
    e.sort(1); e = np.unique(e, axis=0)
    L = np.linalg.norm(xy[e[:, 0]] - xy[e[:, 1]], axis=1)
    return e[L <= lmax]


def scaffold_regions(xy, cls):
    """Scaffold indices, tumour-region flag, Delaunay edges and boundary edges."""
    sc = np.where(np.isin(cls, SCAF))[0]
    X, c = xy[sc], cls[sc]
    kd = cKDTree(X)
    nn = kd.query(X, k=min(KREG, len(X)))[1]
    inT = (c[nn] == T).mean(1) >= 0.5
    E = delaunay_edges(X)
    a, b = E.T
    cr = inT[a] != inT[b]
    bt = np.where(inT[a], a, b)[cr]
    bo = np.where(inT[a], b, a)[cr]
    return sc, X, c, kd, inT, E, bt, bo


class Continuity:
    """Perimeter fraction covered by a contiguous rind fibroblast layer, fast for permutations."""
    def __init__(self, n, E, bo):
        self.bo = bo
        self.rind = np.unique(bo)
        pos = -np.ones(n, int); pos[self.rind] = np.arange(len(self.rind))
        a, b = pos[E[:, 0]], pos[E[:, 1]]
        k = (a >= 0) & (b >= 0)
        self.ra, self.rb, self.pos, self.nr = a[k], b[k], pos, len(self.rind)

    def __call__(self, isS):
        s = isS[self.rind]
        m = s[self.ra] & s[self.rb]
        g = coo_matrix((np.ones(m.sum()), (self.ra[m], self.rb[m])), shape=(self.nr, self.nr))
        lab = connected_components(g, directed=False)[1]
        size = np.bincount(lab)[lab]
        fr = {}
        for L in RUNS:
            gbL = (s & (size >= L))[self.pos[self.bo]]
            fr[L] = gbL.mean()
            if L == RUN: gb = gbL
        return fr, s[self.pos[self.bo]].mean(), gb


def tissue_grid(xy):
    lo, hi = xy.min(0), xy.max(0)
    gx = np.arange(lo[0], hi[0] + PIX, PIX); gy = np.arange(lo[1], hi[1] + PIX, PIX)
    G = np.stack(np.meshgrid(gx, gy), -1).reshape(-1, 2)
    d, j = cKDTree(xy).query(G)
    return G[d <= TISR], j[d <= TISR]


def walk(indptr, indices, blk_edge, inT, start_pool, n, p, rng):
    """n random walkers on the scaffold graph; blocked crossings are refused with prob p."""
    pos = rng.choice(start_pool, n, replace=True)
    deg = np.diff(indptr)
    for _ in range(NSTEP):
        dg = deg[pos]; ok = dg > 0
        j = indptr[pos] + (rng.random(n) * np.maximum(dg, 1)).astype(int)
        j = np.where(ok, j, 0)
        nxt = indices[j]
        stop = ~ok | (blk_edge[j] & (rng.random(n) < p))
        pos = np.where(stop, pos, nxt)
    return inT[pos]


def lir(n_in, a_in, n_out, a_out):
    return np.log2((n_in + 0.5) / a_in) - np.log2((n_out + 0.5) / a_out)


# ------------------------------------------------------------------ per image
def image_metrics(args):
    key, xy, cls, cd163, cd163_thr, seed, keep_cells = args
    rng = np.random.default_rng(seed)
    out = {"image": key}
    if (cls == T).sum() < 50 or np.isin(cls, SCAF).sum() < 100:
        return out, None
    sc, X, c, kd, inT, E, bt, bo = scaffold_regions(xy, cls)
    n = len(X); isS = c == S
    # areas
    G, jn = tissue_grid(xy)
    js = kd.query(G)[1]; regT = inT[js]
    px = PIX * PIX / 1e6
    A_T, A_O = regT.sum() * px, (~regT).sum() * px
    # every cell's region = region of its nearest scaffold cell
    reg_all = inT[kd.query(xy)[1]]
    cd8 = cls == CD8
    out.update(n_cells=len(xy), n_T=int((cls == T).sum()), n_S=int((cls == S).sum()),
               n_CD8=int(cd8.sum()), n_M=int((cls == M).sum()), n_B=int((cls == B).sum()),
               A_tissue=A_T + A_O, A_T=A_T, A_O=A_O,
               cd8_in=int((cd8 & reg_all).sum()), cd8_out=int((cd8 & ~reg_all).sum()),
               n_bedge=len(bo), n_bcell=len(np.unique(bt)))
    ok = len(bo) >= 20 and A_T >= 0.01 and A_O >= 0.01
    out["ok_barrier"] = ok
    cont_ex = None
    if ok:
        C = Continuity(n, E, bo)
        fr, cover, gb = C(isS)
        nonT = np.where(~inT)[0]; sv = isS[nonT]
        nc = {L: np.empty(NPERM) for L in RUNS}; nv = np.empty(NPERM)
        for p in range(NPERM):
            s2 = isS.copy(); s2[nonT] = rng.permutation(sv)
            f2, nv[p], _ = C(s2)
            for L in RUNS: nc[L][p] = f2[L]
        out.update(cont=fr[RUN], cover=cover, cont_null=nc[RUN].mean(),
                   cont_null_sd=nc[RUN].std(), cover_null=nv.mean(), good_edges=int(gb.sum()),
                   S_frac_nonT=float(sv.mean()))
        for L in RUNS:
            out[f"cont_r{L}"] = fr[L]; out[f"cont_null_r{L}"] = nc[L].mean()
        # positive control: random walkers, crossings rind-fibroblast <-> tumour blocked
        a, b = E.T
        A = csr_matrix((np.ones(2 * len(E)), (np.r_[a, b], np.r_[b, a])), shape=(n, n))
        A.sort_indices()
        src = np.repeat(np.arange(n), np.diff(A.indptr)); dst = A.indices
        outer = np.where(inT[src], dst, src)
        blk = (inT[src] != inT[dst]) & isS[outer]
        nw = max(int(cd8.sum()), 1)
        for p in PBLOCK:
            for r in range(NREP):
                w = walk(A.indptr, A.indices, blk, inT, nonT, nw, p, rng)
                out[f"sim_in_p{p}_r{r}"] = int(w.sum()); out[f"sim_out_p{p}_r{r}"] = int((~w).sum())
        if keep_cells:
            rs = np.zeros(n, bool); rs[C.rind] = True
            goodrind = np.zeros(n, bool); goodrind[C.bo[gb]] = True
            cont_ex = dict(sc=sc, inT=inT, rind=rs, goodrind=goodrind, E=E)
    # macrophage density near vs far from B/T aggregates
    lym = np.isin(cls, LYMPH)
    agg_members = np.zeros(len(xy), bool)
    if lym.sum() >= 25:
        li = np.where(lym)[0]; kl = cKDTree(xy[li])
        cnt = kl.query_ball_point(xy[li], r=30.0, return_length=True) - 1
        dn = li[cnt >= 8]
        if len(dn) >= 25:
            pr = cKDTree(xy[dn]).query_pairs(20.0, output_type="ndarray")
            g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(dn), len(dn)))
            lab = connected_components(g, directed=False)[1]
            for k in np.unique(lab):
                mem = dn[lab == k]
                if len(mem) >= 25 and (cls[mem] == B).sum() >= 3:
                    agg_members[mem] = True
    out["n_agg_cells"] = int(agg_members.sum())
    if agg_members.any():
        ka = cKDTree(xy[agg_members])
        dc = ka.query(xy)[0]; dp = ka.query(G)[0]
        pin = ~agg_members[jn]              # pixel's nearest cell is not an aggregate member
        near_c = (~agg_members) & (dc <= NEAR); far_c = dc > FAR
        near_p = pin & (dp <= NEAR); far_p = dp > FAR
        m163 = (cls == M) & (cd163 > cd163_thr) if cd163_thr == cd163_thr else np.zeros(len(xy), bool)
        for tag, rc, rp in [("all", np.ones(len(xy), bool), np.ones(len(G), bool)),
                            ("str", ~reg_all, ~regT)]:
            out[f"agg_Anear_{tag}"] = (near_p & rp).sum() * px
            out[f"agg_Afar_{tag}"] = (far_p & rp).sum() * px
            for nm, msk in [("M", cls == M), ("M163", m163), ("S", cls == S),
                            ("scaf", np.isin(cls, SCAF)), ("CD8", cd8)]:
                out[f"agg_{nm}_near_{tag}"] = int((msk & near_c & rc).sum())
                out[f"agg_{nm}_far_{tag}"] = int((msk & far_c & rc).sum())
    return out, cont_ex


# ------------------------------------------------------------------ models
def zs(v):
    v = np.asarray(v, float); return (v - v.mean()) / v.std(ddof=1)


def ols(y, X, names, hc="HC3"):
    import statsmodels.api as sm
    Xm = sm.add_constant(np.column_stack(X))
    f = sm.OLS(y, Xm).fit(cov_type=hc)
    ci = f.conf_int()
    return [dict(term=nm, coef=f.params[i + 1], lo=ci[i + 1][0], hi=ci[i + 1][1],
                 p=f.pvalues[i + 1]) for i, nm in enumerate(names)]


def cluster_boot(df, ycol, xcols, B=2000, seed=1):
    rng = np.random.default_rng(seed)
    pats = df.patient.unique(); co = []
    for _ in range(B):
        pk = rng.choice(pats, len(pats), replace=True)
        d = pd.concat([df[df.patient == q] for q in pk])
        X = np.column_stack([np.ones(len(d))] + [d[c].values for c in xcols])
        try:
            co.append(np.linalg.lstsq(X, d[ycol].values, rcond=None)[0][1:])
        except Exception:
            pass
    co = np.array(co)
    return np.percentile(co, [2.5, 97.5], axis=0), co


def pool_patient(df):
    """Sum counts/areas per patient, then recompute the metrics."""
    sums = ["n_cells", "n_T", "n_S", "n_CD8", "A_tissue", "A_T", "A_O", "cd8_in", "cd8_out",
            "n_bedge", "good_edges"] + [c for c in df.columns if c.startswith("sim_")]
    g = df.copy()
    ew = ["cont", "cont_null", "cover"] + [f"cont{x}_r{L}" for L in RUNS for x in ("", "_null")]
    for c in ew:                       # edge-weighted pooling of per-image fractions
        g["w_" + c] = g[c] * g.n_bedge
    P = g.groupby("patient")[sums + ["w_" + c for c in ew]].sum()
    P["n_img"] = g.groupby("patient").size()
    for c in ew:
        P[c] = P["w_" + c] / P.n_bedge
    P["cont_exc"] = P.cont - P.cont_null
    for L in RUNS:
        P[f"cont_exc_r{L}"] = P[f"cont_r{L}"] - P[f"cont_null_r{L}"]
    P["lir"] = lir(P.cd8_in, P.A_T, P.cd8_out, P.A_O)
    P["logS"] = np.log10((P.n_S + 0.5) / P.A_tissue)
    P["logCD8"] = np.log10((P.n_CD8 + 0.5) / P.A_tissue)
    return P.reset_index()


def add_image_vars(df):
    df = df.copy()
    df["cont_exc"] = df.cont - df.cont_null
    for L in RUNS:
        df[f"cont_exc_r{L}"] = df[f"cont_r{L}"] - df[f"cont_null_r{L}"]
    df["lir"] = lir(df.cd8_in, df.A_T, df.cd8_out, df.A_O)
    df["logS"] = np.log10((df.n_S + 0.5) / df.A_tissue)
    df["logCD8"] = np.log10((df.n_CD8 + 0.5) / df.A_tissue)
    return df


SPECS = {  # name -> predictors (z-scored); outcome = log2 infiltration ratio
    "continuity": ["cont", "logS", "logCD8"],
    "excess continuity": ["cont_exc", "logS", "logCD8"],
    "decomposed": ["cont_exc", "cont_null", "logS", "logCD8"],
    "perimeter cover": ["cover", "logS", "logCD8"],
    "continuity, run>=3": ["cont_r3", "logS", "logCD8"],
    "continuity, run>=10": ["cont_r10", "logS", "logCD8"],
    "excess, run>=3": ["cont_exc_r3", "logS", "logCD8"],
    "excess, run>=10": ["cont_exc_r10", "logS", "logCD8"],
}


def fit_h1(D, level, dataset, spec, ycol="lir", label="observed", boot=True):
    """ycol ~ z(predictors of SPECS[spec]); ROI level = patient random intercept."""
    d = D.copy(); preds = SPECS[spec]
    for c in preds:
        d["z_" + c] = zs(d[c])
    xc = ["z_" + c for c in preds]
    rows = []
    if level == "roi":
        import statsmodels.formula.api as smf
        f = smf.mixedlm(f"{ycol} ~ " + " + ".join(xc), d, groups=d.patient).fit(reml=True)
        ci = f.conf_int()
        blo = bhi = [np.nan] * len(xc)
        if boot:
            (blo, bhi), _ = cluster_boot(d, ycol, xc)
        for i, c in enumerate(xc):
            rows.append(dict(dataset=dataset, level="ROI, patient random intercept", outcome=ycol,
                             label=label, spec=spec, term=c[2:], coef=f.params[c],
                             lo=ci.loc[c, 0], hi=ci.loc[c, 1], p=f.pvalues[c],
                             boot_lo=blo[i], boot_hi=bhi[i], n=len(d), n_pat=d.patient.nunique()))
    else:
        for r in ols(d[ycol].values, [d[c].values for c in xc], [c[2:] for c in xc]):
            rows.append(dict(dataset=dataset, level="patient (pooled), OLS HC3", outcome=ycol,
                             label=label, spec=spec, n=len(d), n_pat=len(d), **r))
    return rows


# ------------------------------------------------------------------ main
def run_dataset(df, name, pool, cd163_thr, examples_wanted=()):
    t0 = time.time()
    imgs = df.image.unique()
    tasks = []
    for i, k in enumerate(imgs):
        g = df[df.image == k]
        tasks.append((k, g[["x", "y"]].to_numpy(float), g.cls.to_numpy(),
                      g.cd163.to_numpy(float), cd163_thr, [{"PDAC": 1, "METABRIC": 2}[name], i],
                      k in examples_wanted))
    res = list(pool.imap(image_metrics, tasks, chunksize=4))
    rows = [r[0] for r in res]
    ex = {r[0]["image"]: r[1] for r in res if r[1] is not None}
    M = pd.DataFrame(rows).merge(df.groupby("image").patient.first().reset_index(), on="image")
    M.insert(0, "dataset", name)
    print(f"  {name}: {len(imgs)} images, {time.time() - t0:.0f} s")
    return M, ex


def macro_summary(M, name, tag, cellset=("M", "M163", "S", "scaf", "CD8")):
    need = [f"agg_Anear_{tag}", f"agg_Afar_{tag}"]
    if need[0] not in M: return [], None
    d = M[M[need[0]].notna()]
    cols = need + [f"agg_{c}_{w}_{tag}" for c in cellset for w in ("near", "far")]
    P = d.groupby("patient")[cols].sum()
    P = P[(P[need[0]] >= 0.005) & (P[need[1]] >= 0.02)]
    out = []
    for c in cellset:
        if P[f"agg_{c}_near_{tag}"].sum() + P[f"agg_{c}_far_{tag}"].sum() == 0: continue
        lr = np.log2((P[f"agg_{c}_near_{tag}"] + 0.5) / P[need[0]]) - \
             np.log2((P[f"agg_{c}_far_{tag}"] + 0.5) / P[need[1]])
        P[f"lr_{c}"] = lr
    if "lr_M" in P and "lr_scaf" in P:
        P["lr_M_adj"] = P.lr_M - P.lr_scaf
    if "lr_M163" in P and "lr_scaf" in P and P[f"agg_M163_near_{tag}"].sum() > 0:
        P["lr_M163_adj"] = P.lr_M163 - P.lr_scaf
        P["lr_M163_share"] = P.lr_M163 - P.lr_M   # CD163-high share of macrophages, near vs far
    rng = np.random.default_rng(3)
    for c in [x for x in P.columns if x.startswith("lr_")]:
        v = P[c].dropna().values
        if len(v) < 5: continue
        bs = [np.median(rng.choice(v, len(v))) for _ in range(2000)]
        w = stats.wilcoxon(v)
        out.append(dict(dataset=name, region=tag, contrast=c[3:], n_pat=len(v),
                        median_log2=np.median(v), lo=np.percentile(bs, 2.5),
                        hi=np.percentile(bs, 97.5), n_lower_near=int((v < 0).sum()),
                        p_wilcoxon=w.pvalue))
    return out, P.reset_index().assign(dataset=name, region=tag)


def survival(P):
    from statsmodels.duration.hazard_regression import PHReg
    raw = json.load(open(f"{PI}/metabric_clinical_raw.json"))
    cl = pd.DataFrame(raw).pivot_table(index="patientId", columns="clinicalAttributeId",
                                       values="value", aggfunc="first")
    cl = cl[["AGE_AT_DIAGNOSIS", "NPI", "ER_IHC", "OS_MONTHS", "VITAL_STATUS"]].dropna()
    cl["time"] = cl.OS_MONTHS.astype(float); cl["age"] = cl.AGE_AT_DIAGNOSIS.astype(float)
    cl["npi"] = cl.NPI.astype(float); cl["er"] = (cl.ER_IHC.str.startswith("Pos")).astype(float)
    cl["dss"] = (cl.VITAL_STATUS == "Died of Disease").astype(int)
    cl["death"] = (cl.VITAL_STATUS != "Living").astype(int)
    d = P.merge(cl[["time", "age", "npi", "er", "dss", "death"]], left_on="patient",
                right_index=True)
    d = d[d.time > 0]
    rows = []
    for c in ("cont", "cont_exc", "logS", "lir", "logCD8", "cover"):
        d["z_" + c] = zs(d[c])
    specs = [(f"{c} (+age, NPI, ER)", ["z_" + c, "age", "npi", "er"])
             for c in ("cont", "cont_exc", "logS", "lir", "logCD8")]
    specs += [("joint: cont + logS + logCD8 (+age, NPI, ER)",
               ["z_cont", "z_logS", "z_logCD8", "age", "npi", "er"])]
    for ev in ("dss", "death"):
        for lab, xs in specs:
            f = PHReg(d.time.values, d[xs].values, status=d[ev].values).fit()
            ci = f.conf_int()
            for i, x in enumerate(xs):
                if not x.startswith("z_"): continue
                rows.append(dict(endpoint=ev, model=lab, term=x[2:], HR=np.exp(f.params[i]),
                                 lo=np.exp(ci[i][0]), hi=np.exp(ci[i][1]), p=f.pvalues[i],
                                 n=len(d), events=int(d[ev].sum())))
    return pd.DataFrame(rows), d


def main():
    from multiprocessing import Pool
    t0 = time.time()
    pdac = load_pdac(); print(f"PDAC cells {len(pdac):,}  ROIs {pdac.image.nunique()}  "
                              f"patients {pdac.patient.nunique()}")
    mb = load_metabric()
    thr = float(np.nanmedian(mb.cd163[mb.cls == M]))
    print(f"METABRIC cells {len(mb):,}  images {mb.image.nunique()}  patients "
          f"{mb.patient.nunique()}  (CD163 threshold, median among macrophages {thr:.3f})")
    with Pool(8) as pool:
        Mp, exP = run_dataset(pdac, "PDAC", pool, np.nan, examples_wanted=set(pdac.image.unique()))
        Mb, exB = run_dataset(mb, "METABRIC", pool, thr)
    IM = pd.concat([Mp, Mb], ignore_index=True)
    IM.to_csv(f"{R}/E13_image_metrics.csv", index=False)

    summ = {"hypotheses": {"H1": "continuous interface fibroblast layer -> fewer CD8 in tumour at fixed stromal density",
                           "H0": "CD8 infiltration tracks CD8 and stromal density, not continuity"},
            "data": {}}
    for nm, d, raw in [("PDAC", Mp, pdac), ("METABRIC", Mb, mb)]:
        summ["data"][nm] = dict(cells=int(len(raw)), images=int(raw.image.nunique()),
                                patients=int(raw.patient.nunique()),
                                images_barrier=int(d.ok_barrier.fillna(False).sum()))

    models, pos = [], []
    # PDAC, ROI level (mixed model + cluster bootstrap) and patient level
    dP = add_image_vars(Mp[Mp.ok_barrier.fillna(False) & (Mp.n_CD8 >= 10)])
    PP = pool_patient(Mp[Mp.ok_barrier.fillna(False)])
    PB = pool_patient(Mb[Mb.ok_barrier.fillna(False)])
    PB = PB[PB.n_CD8 >= 10]
    for sp in SPECS:
        models += fit_h1(dP, "roi", "PDAC", sp)
        models += fit_h1(PB, "patient", "METABRIC", sp)
    models += fit_h1(PP, "patient", "PDAC", "continuity")
    summ["data"]["PDAC"]["rois_model"] = len(dP)
    summ["data"]["METABRIC"]["patients_model"] = len(PB)
    PP["dataset"] = "PDAC"; PB["dataset"] = "METABRIC"
    pd.concat([PP, PB]).to_csv(f"{R}/E13_patient_metrics.csv", index=False)
    dP.to_csv(f"{R}/E13_pdac_roi_model_data.csv", index=False)

    # positive control: the same models on simulated infiltration
    for nm, D, lev in [("PDAC", dP, "roi"), ("METABRIC", PB, "patient")]:
        for p in PBLOCK:
            for r in range(NREP):
                d = D.copy()
                d["simlir"] = lir(d[f"sim_in_p{p}_r{r}"], d.A_T, d[f"sim_out_p{p}_r{r}"], d.A_O)
                for sp in ("continuity", "excess continuity", "decomposed"):
                    for rr in fit_h1(d, lev, nm, sp, ycol="simlir", boot=False):
                        pos.append(dict(dataset=nm, p_block=p, rep=r, spec=sp, term=rr["term"],
                                        coef=rr["coef"], lo=rr["lo"], hi=rr["hi"], p=rr["p"],
                                        mean_simlir=d.simlir.mean()))
    MO = pd.DataFrame(models); MO.to_csv(f"{R}/E13_models.csv", index=False)
    PO = pd.DataFrame(pos); PO.to_csv(f"{R}/E13_poscontrol.csv", index=False)

    # macrophages next to B/T aggregates
    macro, mp = [], []
    for nm, D in [("PDAC", Mp), ("METABRIC", Mb)]:
        for tag in ("all", "str"):
            o, P = macro_summary(D, nm, tag)
            macro += o
            if P is not None: mp.append(P)
    MA = pd.DataFrame(macro); MA.to_csv(f"{R}/E13_macro_aggregates.csv", index=False)
    pd.concat(mp).to_csv(f"{R}/E13_macro_aggregates_patients.csv", index=False)

    # survival (METABRIC)
    SV, dS = survival(PB)
    SV.to_csv(f"{R}/E13_survival.csv", index=False)

    # descriptive correlations
    cor = []
    for nm, D in [("PDAC ROI", dP), ("METABRIC patient", PB)]:
        for a, b in [("cont", "logS"), ("cont_exc", "logS"), ("cont", "lir"), ("cont_exc", "lir"),
                     ("logS", "lir"), ("logCD8", "lir"), ("cont", "cont_null")]:
            r, p = stats.spearmanr(D[a], D[b]); cor.append(dict(set=nm, a=a, b=b, rho=r, p=p))
    pd.DataFrame(cor).to_csv(f"{R}/E13_correlations.csv", index=False)

    # example maps for the figure: PDAC ROIs (all kept) - store cell tables of three
    keep = []
    for k, e in exP.items():
        g = pdac[pdac.image == k].reset_index(drop=True)
        g["scaf"] = False; g["inT"] = False; g["rind"] = False; g["goodrind"] = False
        g.loc[e["sc"], "scaf"] = True; g.loc[e["sc"], "inT"] = e["inT"]
        g.loc[e["sc"], "rind"] = e["rind"]; g.loc[e["sc"], "goodrind"] = e["goodrind"]
        keep.append(g)
    pd.concat(keep).to_csv(f"{R}/E13_pdac_cells_annotated.csv.gz", index=False)

    summ["models"] = MO.to_dict("records"); summ["poscontrol_mean"] = \
        PO.groupby(["dataset", "p_block", "spec", "term"])[["coef", "lo", "hi", "p"]].mean().reset_index().to_dict("records")
    summ["macro"] = MA.to_dict("records"); summ["survival"] = SV.to_dict("records")
    summ["survival_n"] = dict(n=len(dS), dss=int(dS.dss.sum()), death=int(dS.death.sum()))
    json.dump(summ, open(f"{R}/E13_summary.json", "w"), indent=1, default=float)

    pd.set_option("display.width", 200)
    print("\n== H1 models =="); print(MO.drop(columns=["outcome", "label"]).round(3).to_string())
    print("\n== positive control (mean over replicates) ==")
    q = PO[PO.term.isin(["cont", "cont_exc", "cont_null"])]
    print(q.groupby(["dataset", "spec", "term", "p_block"])[["coef", "lo", "hi", "p", "mean_simlir"]].mean().round(3).to_string())
    print("\n== macrophage near/far ==")
    print(MA.round(3).to_string())
    print("\n== survival ==")
    print(SV.round(3).to_string())
    print("\n== correlations =="); print(pd.DataFrame(cor).round(3).to_string())
    print(json.dumps(summ["data"], indent=1), summ["survival_n"])
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
