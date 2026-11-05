"""v19 Fig. 6 - immune scores favor the checkpoint arm in randomized comparisons.

a  design of the six primary and two sensitivity comparisons (Table 3; n from
   manuscript/make_v18_tables.trial_effect)
b  benefit-oriented score x arm interaction per SD, every score and comparison
   (E21_meta_inputs.csv, set "primary: 6 comparisons" and "+ PRINCE + DUTRENEO (8)"),
   pooled rows from E21_meta.csv (the estimates of Table 4)
c  within-arm NetImmune and CAF_ECM effects in each primary comparison (E5_trials.csv,
   E10_ispy2_within_arm.csv, E21_within_arm.csv): what each interaction is made of
d  pCR rate by NetImmune tertile in each arm of NeoTRIP and I-SPY2 pembrolizumab (descriptive;
   per-patient scores from E21_scores.csv and, for I-SPY2, rescored exactly as in
   analysis/E10_ispy2_trials.py; tertiles over both arms)
e  progression-free survival by NetImmune tertile in each arm of JAVELIN Renal 101 (descriptive
   Kaplan-Meier; per-patient scores from E21_scores.csv, PFS from JR101_S11_clinical.tsv)
f  robustness of the pooled interaction for NetImmune and CAF_ECM (E21_meta.csv)
g  structure-matched random-gene null for the three added trials (E21_null_summary.csv)
"""
import sys
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/cns")
sys.path.insert(0, "/groups/adv2105_gp/yichen/Yi/com/gastric/2/manuscript")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, Polygon
from matplotlib.colors import TwoSlopeNorm, to_rgb
from scipy import stats
import cns_style as S
from cns_style import INK, SLATE, GREY, LIGHT, PAPER, IMMUNE, STROMA, MACRO, fmt_p
import make_v18_tables as T

R = S.RES
W, H = 180, 230
fig = plt.figure(figsize=(W * S.MM, H * S.MM))
cv = fig.add_axes([0, 0, 1, 1]); cv.set_xlim(0, W); cv.set_ylim(H, 0); cv.axis("off")
CMAP = S.DIV.reversed()


def txt(x, y, s, size=5.0, color=INK, ha="left", va="center", **kw):
    return cv.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, **kw)


from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
PT = 25.4 / 72                    # mm per point


def _ink_top(s, size):
    """Height of the inked part of a bold string above its baseline, mm."""
    return TextPath((0, 0), s, size=size, prop=FontProperties(family=S.fs.SANS, weight="bold")).get_extents().y1 * PT


def heading(xl, xt, ytop, letter, title):
    """Panel letter (8 pt) and bold title (6.5 pt) with their inked tops level at ytop (mm from the top)."""
    W_, H_ = fig.get_size_inches() / S.MM
    fig.text(xl / W_, 1 - (ytop + _ink_top(letter, 8)) / H_, letter, fontsize=8, fontweight="bold",
             va="baseline", ha="left", color=INK)
    cv.text(xt, ytop + _ink_top(title, 6.5), title, fontsize=6.5, fontweight="bold", va="baseline", ha="left",
            color=INK)


def fmtp(p):
    """fmt_p, keeping the trailing zero of two significant figures (0.060, not 0.06)."""
    return f"{p:#.2g}" if 0.001 <= p < 0.1 else fmt_p(p, lead="")


def sgn(v, nd=2):
    return f"{v:+.{nd}f}".replace("-", "−")


def ink_on(c):
    r, g, b = to_rgb(c)
    return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.45 else INK


# =========================================================================================
# data for d (cached)
# =========================================================================================
TERT = S.OUT / "data/fig6_tertiles.csv"
if not TERT.exists():
    sc = pd.read_csv(R / "E21_scores.csv")
    nt = sc[sc.trial == "NeoTRIP"].set_index("sample")
    P = pd.read_csv(S.ROOT / "data_download/validation/NeoTRIP_GSE319641/github_Early_GEP_NeoTRIP/data/"
                    "NeoTRIP_baseline_D1C2_pheno_and_clinical.txt", sep="\t")
    P = P[P.Timepoint == "Baseline"].set_index("Sample_ID")
    nt = nt.join(P[["pCR.num", "Arm"]])
    assert ((nt.Arm == "CT/A").astype(int) == nt.arm).all()
    d1 = pd.DataFrame(dict(trial="NeoTRIP", arm=nt.arm.values, y=nt["pCR.num"].values, s=nt.NetImmune.values))
    V = "/groups/adv2105_gp/yichen/Yi/com/gastric/validation/cohorts"
    X = pd.read_csv(f"{V}/ISPY2_pembro.expr.tsv.gz", sep="\t", index_col=0); X.columns = X.columns.astype(str)
    C = pd.read_csv(f"{V}/ISPY2_pembro.clin.tsv", sep="\t"); C["sample"] = C["sample"].astype(str)
    C = C.set_index("sample").loc[X.columns]
    Z = ((X.T - X.mean(1)) / X.std(1).replace(0, 1)).T
    MOD = dict(
        CAF_ECM=["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
        TAM_SUP=["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
        APM=["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
        IFN=["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
        TLS_B=["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
        T_CELL=["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"])
    Sm = pd.DataFrame({k: Z.loc[[g for g in G if g in Z.index]].mean(0) for k, G in MOD.items()})
    net = Sm[["APM", "IFN", "TLS_B", "T_CELL"]].mean(1) - Sm[["CAF_ECM", "TAM_SUP"]].mean(1)
    d2 = pd.DataFrame(dict(trial="I-SPY2 pembrolizumab", arm=C.arm.values, y=C.y.values, s=net.values))
    d = pd.concat([d1, d2], ignore_index=True)
    d["tertile"] = d.groupby("trial").s.transform(lambda v: pd.qcut(v, 3, labels=False))
    d.to_csv(TERT, index=False)
tert = pd.read_csv(TERT)

# =========================================================================================
# a + b  design and interaction matrix
# =========================================================================================
inp = pd.read_csv(R / "E21_meta_inputs.csv")
meta = pd.read_csv(R / "E21_meta.csv")
PRIM = inp[inp.set == "primary: 6 comparisons"]
SENS = inp[inp.set == "+ PRINCE + DUTRENEO (8)"]
COLS = [("CAF_ECM", "CAF_ECM", STROMA), ("Barrier", "Barrier", STROMA), ("TAM_SUP", "TAM_SUP", MACRO),
        ("C1Q_arm", "C1Q subset", MACRO), ("SPP1_arm", "SPP1 subset", MACRO), ("CXCL9_SPP1", "CXCL9:SPP1", MACRO),
        ("APM", "APM", IMMUNE), ("IFN", "IFN", IMMUNE), ("TLS_B", "TLS_B", IMMUNE), ("T_CELL", "T_CELL", IMMUNE),
        ("Immune", "Immune", IMMUNE), ("NetImmune", "NetImmune", IMMUNE),
        ("ESTIMATE", "ESTIMATE", SLATE), ("TMEscore_PC", "TMEscore PC1", SLATE),
        ("TMEscore_z", "TMEscore z", SLATE), ("CD8A", "CD8A", SLATE)]   # E23 published comparators
CW = 5.2                             # E23: 16 columns instead of 12
XM = W - 2.5 - len(COLS) * CW        # matrix left edge
RH = 5.4
norm = TwoSlopeNorm(vcenter=0, vmin=-1.0, vmax=1.0)
# (row label, trial key in E21_meta_inputs, Table-3 index, cancer, checkpoint arm, control, endpoint)
TR = [("Renal cell carcinoma", None),
      ("IMmotion150", "IMmotion150", 0, "Atezolizumab ± bevacizumab", "sunitinib", "PFS"),
      ("CheckMate 025", "CheckMate 025", 3, "Nivolumab", "everolimus", "PFS"),
      ("JAVELIN Renal 101", "JAVELIN Renal 101", 4, "Avelumab + axitinib", "sunitinib", "PFS"),
      ("Breast cancer", None),
      ("I-SPY2", "I-SPY2 pembrolizumab", 1, "Pembrolizumab + paclitaxel", "paclitaxel", "pCR"),
      ("I-SPY2", "I-SPY2 durvalumab-olaparib", 2, "Durvalumab + olaparib + paclitaxel", "paclitaxel", "pCR"),
      ("NeoTRIP", "NeoTRIP", 5, "Atezolizumab + chemotherapy", "chemotherapy", "pCR"),
      ("pool", "primary: 6 comparisons", "Pooled, six primary"),
      ("Sensitivity comparisons", None),
      ("PRINCE", "PRINCE", 6, "Nivolumab-containing arms", "sotigalimab; all with chemotherapy", "PFS"),
      ("DUTRENEO", "DUTRENEO", 7, "Durvalumab + tremelimumab", "cisplatin-based chemotherapy", "pCR"),
      ("pool", "+ PRINCE + DUTRENEO (8)", "Pooled, all eight")]

XT, XARM, XEP, XN = 4.0, 31.0, 75.5, 85.0   # E23: n moved left of the wider matrix
heading(1.0, 5.0, 1.7, "a", "Randomized comparisons")

heading(XM - 5.0, XM, 1.7, "b", "Score × arm interaction")

YH = 20.0
for i, (_, lab, col) in enumerate(COLS):
    xc = XM + (i + 0.5) * CW
    txt(xc, YH, lab, size=5.4, ha="center", va="bottom", rotation=90)
    cv.add_patch(Rectangle((XM + i * CW + 0.3, YH + 0.6), CW - 0.6, 0.9, color=col, lw=0))
for x, s_, ha in [(XT, "trial", "left"), (XARM, "checkpoint arm / control", "left"), (XEP, "endpoint", "left"),
                  (XN + 6.5, "n", "right")]:
    txt(x, YH + 0.8, s_, size=5.0, color=GREY, va="bottom", ha=ha)

y = YH + 3.4
cells = []
for row in TR:
    if row[1] is None:                                   # section header
        if y > YH + 4:
            y += 1.2
        txt(2.0, y + 1.6, row[0], size=5.4, fontweight="bold")
        cv.plot([2.0, W - 2.5], [y + 3.3, y + 3.3], color=LIGHT, lw=0.4)
        y += 3.8
        continue
    if row[0] == "pool":
        _, setname, lab = row
        y += 0.4
        cv.plot([2.0, W - 2.5], [y, y], color=INK, lw=0.4)
        y += 0.6
        yc = y + 2.5
        txt(XT, yc, lab, size=5.2, fontweight="bold")
        m = meta[meta.set == setname].set_index("score")
        for i, (k, _, _) in enumerate(COLS):
            r = m.loc[k]
            ring = r.p < 0.05
            fc = CMAP(norm(np.clip(r.est, -1, 1)))
            cv.add_patch(Rectangle((XM + i * CW + 0.3, yc - 2.2), CW - 0.6, 4.4, fc=fc,
                                   ec=INK if ring else "none", lw=0.6, zorder=2))
            tc = ink_on(fc)
            txt(XM + (i + 0.5) * CW, yc - 0.85, sgn(r.est), size=5.0, color=tc, ha="center",
                fontweight="bold" if ring else "normal", zorder=3)
            txt(XM + (i + 0.5) * CW, yc + 1.05, f"{int(r.k_benefit)}/{int(r.k)}", size=5.0, color=tc, ha="center",
                zorder=3)
        y += 5.4
        continue
    lab, key, ti, arm1, arm0, ep = row
    src = T.TRIALS[ti]
    kind, e, lo, hi, p, n1, n0 = T.trial_effect(src[4], src[5], "NetImmune")
    yc = y + RH / 2
    txt(XT, yc - 0.95, lab, size=5.2, fontweight="bold" if ti < 6 else "normal")
    txt(XT, yc + 1.15, src[1], size=5.0, color=SLATE)
    txt(XARM, yc - 0.95, arm1, size=5.0)
    txt(XARM, yc + 1.15, "vs " + arm0, size=5.0, color=SLATE)
    txt(XEP, yc, ep, size=5.0)
    txt(XN + 6.5, yc, f"{n1} / {n0}", size=5.0, ha="right")
    src_df = PRIM if ti < 6 else SENS
    for i, (k, _, _) in enumerate(COLS):
        r = src_df[(src_df.trial == key) & (src_df.score == k)].iloc[0]
        xc = XM + (i + 0.5) * CW
        rad = float(np.clip(0.32 / r.se, 0.75, 1.75))
        cv.add_patch(Circle((xc, yc), rad, fc=CMAP(norm(np.clip(r.y, -1, 1))), ec="none", zorder=3))
        if r.p_trial < 0.05:
            cv.add_patch(Circle((xc, yc), rad + 0.35, fc="none", ec=INK, lw=0.5, zorder=4))
        cells.append(dict(trial=key, score=k, y=r.y, se=r.se, p=r.p_trial))
    y += RH
YAB = y

# legend strip for b
YL = YAB + 2.6
n_ = 120; x0, w0 = XM, 32
for j in range(n_):
    v = -1 + 2 * (j + 0.5) / n_
    cv.add_patch(Rectangle((x0 + w0 * j / n_, YL), w0 / n_ + 0.02, 1.5, color=CMAP(norm(v)), lw=0))
for t, l in [(-1, "−1"), (0, "0"), (1, "+1")]:
    xx = x0 + w0 * (t + 1) / 2
    cv.plot([xx, xx], [YL + 1.5, YL + 2.1], color=INK, lw=0.4)
    txt(xx, YL + 3.1, l, size=5.0, ha="center")
txt(x0, YL - 1.3, "favors control", size=5.0, color=SLATE)
txt(x0 + w0, YL - 1.3, "favors checkpoint", size=5.0, color=SLATE, ha="right")
xs = x0 + w0 + 6
txt(xs - 0.6, YL + 0.8, "SE", size=5.0, color=SLATE, ha="right")
for j, se in enumerate([0.2, 0.4]):
    rad = float(np.clip(0.32 / se, 0.75, 1.75))
    cv.add_patch(Circle((xs + 1.6 + j * 5.5, YL + 0.8), rad, fc=GREY, ec="none"))
    txt(xs + 1.6 + j * 5.5, YL + 3.3, f"{se:g}", size=5.0, ha="center")
cv.add_patch(Circle((xs + 13.5, YL + 0.8), 1.0, fc=GREY, ec="none"))
cv.add_patch(Circle((xs + 13.5, YL + 0.8), 1.35, fc="none", ec=INK, lw=0.5))
txt(xs + 15.4, YL + 0.8, "P < 0.05", size=5.0, color=SLATE)



# =========================================================================================
# row geometry
# =========================================================================================
Y2 = YL + 11.0                    # title line, row 2 (c, d, e)
Y2A = Y2 + 8.5; h2 = 38.0; B2 = Y2A + h2
Y3 = B2 + 13.5                    # title line, row 3 (f, g)
Y3A = Y3 + 8.5; h3 = H - 11.5 - Y3A

# =========================================================================================
# c  within-arm effects: what each interaction is made of
# =========================================================================================
heading(1.0, 5.0, Y2 - 0.85, "c", "Within-arm effects")

e5 = pd.read_csv(R / "E5_trials.csv")
e5 = e5[e5.trial == "IMmotion150: atezolizumab (+/- bev) vs sunitinib"].set_index("score")
e10 = pd.read_csv(R / "E10_ispy2_within_arm.csv"); e10 = e10[e10.model == "HRadj"]
e21 = pd.read_csv(R / "E21_within_arm.csv")


def within(trial, k, arm):
    """Benefit-oriented within-arm effect per SD (-log HR for PFS, log OR for pCR) and its 95% CI."""
    if trial == "IMmotion150":
        r = e5.loc[k]; pre = "ici" if arm == "checkpoint" else "ctrl"
        return -np.log(r[f"{pre}_HR"]), -np.log(r[f"{pre}_hi"]), -np.log(r[f"{pre}_lo"]), int(r.n_ici if pre == "ici" else r.n_ctrl)
    if trial.startswith("I-SPY2"):
        key = "ISPY2_pembro" if "pembro" in trial else "ISPY2_durva"
        r = e10[(e10.cohort == key) & (e10.score == k) & (e10.arm == arm)].iloc[0]
        return r.logOR, r.logOR - 1.96 * r.se, r.logOR + 1.96 * r.se, int(r.n)
    if trial == "NeoTRIP":
        r = e21[(e21.trial == trial) & (e21.model == "logistic + batch") & (e21.score == k) & (e21.arm == arm)].iloc[0]
        return r.log, r.log - 1.96 * r.se, r.log + 1.96 * r.se, int(r.n)
    r = e21[(e21.trial == trial) & (e21.endpoint == "PFS") & (e21.score == k) & (e21.arm == arm)].iloc[0]
    return -r.log, -r.log - 1.96 * r.se, -r.log + 1.96 * r.se, int(r.n)


CT_ = [("IMmotion150", "IMmotion150", "PFS"), ("CheckMate 025", "CheckMate 025", "PFS"),
       ("JAVELIN Renal 101", "JAVELIN Renal 101", "PFS"), ("I-SPY2 pembrolizumab", "I-SPY2 pembro.", "pCR"),
       ("I-SPY2 durvalumab-olaparib", "I-SPY2 durva.–olap.", "pCR"), ("NeoTRIP", "NeoTRIP", "pCR")]
rows_c = []
for j, (k, col, lab) in enumerate([("NetImmune", IMMUNE, "NetImmune"), ("CAF_ECM", STROMA, "CAF_ECM")]):
    ax = S.ax_mm(fig, 26.0 + j * 25.0, Y2A, 21.0, h2)
    for i, (tk, tl, ep) in enumerate(CT_):
        yy = -i - (0.4 if i >= 3 else 0)
        v1 = within(tk, k, "checkpoint"); v0 = within(tk, k, "control")
        ax.plot([v0[0], v1[0]], [yy, yy], color=LIGHT, lw=1.6, solid_capstyle="butt", zorder=1)
        for v, dy, filled in [(v1, 0.17, True), (v0, -0.17, False)]:
            ax.plot([v[1], v[2]], [yy + dy, yy + dy], color=col if filled else GREY, lw=0.6, solid_capstyle="butt",
                    zorder=2)
            ax.scatter([v[0]], [yy + dy], s=11, marker="o", color=col if filled else "white",
                       ec=col if filled else GREY, lw=0.6, zorder=3)
        rows_c.append(dict(trial=tk, score=k, endpoint=ep, chk=v1[0], chk_lo=v1[1], chk_hi=v1[2], n_chk=v1[3],
                           ctrl=v0[0], ctrl_lo=v0[1], ctrl_hi=v0[2], n_ctrl=v0[3]))
        if j == 0:
            ax.text(-0.05, yy, tl, transform=ax.get_yaxis_transform(), fontsize=5.0, ha="right", va="center")
    ax.axhline(-2.7, color=LIGHT, lw=0.4)
    if j == 1:
        ax.text(1.72, -2.7 + 0.22, "renal", fontsize=5.0, color=SLATE, ha="right", va="bottom")
        ax.text(1.72, -2.7 - 0.22, "breast", fontsize=5.0, color=SLATE, ha="right", va="top")
    S.rule(ax, 0)
    ax.set_ylim(-5.4 - 0.55, 0.55); ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(-1.75, 1.75); ax.set_xticks([-1, 0, 1]); ax.set_xticklabels(["−1", "0", "+1"])
    ax.text(0.0, 1.03, lab, transform=ax.transAxes, fontsize=5.4, fontweight="bold", va="bottom")
txt(26.0 + 23.0, B2 + 6.4, "Better outcome per SD", size=6.0, ha="center")
txt(26.0 + 23.0, B2 + 9.2, "−log HR (PFS) or log OR (pCR)", size=5.0, color=SLATE, ha="center")
# key: filled = checkpoint arm, open = control arm
kx = 3.0; ky = Y2A + h2 + 3.0
for j_, (lab_, filled) in enumerate([("checkpoint arm", True), ("control arm", False)]):
    cv.scatter([kx + 1.0], [ky + j_ * 2.6], s=11, marker="o", color=INK if filled else "white", ec=INK if filled else GREY,
               lw=0.6, zorder=5)
    txt(kx + 2.4, ky + j_ * 2.6, lab_, size=5.0)

# =========================================================================================
# e  PFS by NetImmune tertile in each arm of JAVELIN Renal 101 (descriptive Kaplan-Meier)
# =========================================================================================
KMF = S.OUT / "data/fig6_javelin_pfs.csv"
if not KMF.exists():
    sc = pd.read_csv(R / "E21_scores.csv"); jv = sc[sc.trial == "JAVELIN Renal 101"].set_index("sample")
    Cj = pd.read_csv(S.ROOT / "data_download/validation/JAVELIN_Renal101/JR101_S11_clinical.tsv", sep="\t").set_index("ID")
    jv = jv.join(Cj[["PFS_P", "PFS_P_CNSR", "TRT01P"]])
    assert ((jv.TRT01P == "Avelumab+Axitinib").astype(int) == jv.arm).all()
    kd = pd.DataFrame(dict(arm=jv.arm.values, s=jv.NetImmune.values, time=jv.PFS_P.values,
                           event=(jv.PFS_P_CNSR.values == 0).astype(int)))
    kd["tertile"] = pd.qcut(kd.s, 3, labels=False)
    kd.to_csv(KMF, index=False)
kd = pd.read_csv(KMF)


def km(t, e):
    o = np.argsort(t); t, e = t[o], e[o]
    ut = np.unique(t[e == 1]); S_, tt, ss = 1.0, [0.0], [1.0]
    for u in ut:
        n = (t >= u).sum(); d = ((t == u) & (e == 1)).sum()
        S_ *= 1 - d / n; tt.append(u); ss.append(S_)
    return np.array(tt), np.array(ss)


def km_median(t, e):
    tt, ss = km(t, e)
    below = np.where(ss <= 0.5)[0]
    return tt[below[0]] if len(below) else np.nan


XE6 = 130
heading(XE6 - 4.0, XE6, Y2 - 0.85, "e", "Renal: PFS by NetImmune tertile")

TCOL = ["#9CC3E6", "#3F8FD8", "#0A3D6E"]
km_rows = []
for j, (a, lab) in enumerate([(1, "Avelumab + axitinib"), (0, "Sunitinib")]):
    ax = S.ax_mm(fig, XE6 + 6.5 + j * 22.5, Y2A, 18.0, h2)
    for q in range(3):
        g = kd[(kd.arm == a) & (kd.tertile == q)]
        tt, ss = km(g.time.values, g.event.values)
        tcut = np.sort(g.time.values)[::-1][9]          # last time with at least 10 patients at risk
        keep = tt <= tcut
        ax.step(np.r_[tt[keep], tcut], np.r_[ss[keep], ss[keep][-1]], where="post", color=TCOL[q], lw=0.9)
        km_rows.append(dict(arm="avelumab + axitinib" if a else "sunitinib", tertile=["low", "mid", "high"][q],
                            n=len(g), events=int(g.event.sum()), median_pfs_months=km_median(g.time.values, g.event.values),
                            curve_stops_months=tcut))
        yk = 0.19 - q * 0.08
        ax.plot([1.5, 4.5], [yk, yk], color=TCOL[q], lw=0.9)
        ax.text(5.3, yk, f"{['low', 'mid', 'high'][q]} ({len(g)})", fontsize=5.0, va="center")
    ax.set_xlim(0, 24); ax.set_ylim(0, 1.02)
    ax.set_xticks([0, 6, 12, 18, 24]); ax.set_xticklabels(["0", "", "12", "", "24"])
    ax.set_yticks([0, 0.5, 1]); ax.set_yticklabels(["0", "0.5", "1"] if j == 0 else [])
    ax.text(0.0, 1.03, lab, transform=ax.transAxes, fontsize=5.4, fontweight="bold", va="bottom")
    if j == 0:
        ax.set_ylabel("Progression-free")
txt(XE6 + 6.5 + 20.25, B2 + 6.4, "Months", size=6.0, ha="center")

# =========================================================================================
# f  robustness of the pooled interaction (row 3)
# =========================================================================================
SETS = [("primary: 6 comparisons", "All six (primary)", True),
        ("primary without IMmotion150", "− IMmotion150", False),
        ("primary without CheckMate 025", "− CheckMate 025", False),
        ("primary without JAVELIN Renal 101", "− JAVELIN Renal 101", False),
        ("primary without I-SPY2 pembrolizumab", "− I-SPY2 pembrolizumab", False),
        ("primary without I-SPY2 durvalumab-olaparib", "− I-SPY2 durvalumab–olaparib", False),
        ("primary without NeoTRIP", "− NeoTRIP", False),
        ("renal only (IMmotion150, CheckMate 025, JAVELIN)", "Renal only (3)", True),
        ("breast only (I-SPY2 x 2, NeoTRIP)", "Breast only (3)", False),
        ("+ PRINCE + DUTRENEO (8)", "+ PRINCE, DUTRENEO (8)", True),
        ("primary, CheckMate 025 OS for PFS", "CheckMate 025 with OS", False),
        ("primary, NeoTRIP without batch", "NeoTRIP without batch term", False)]
heading(1.0, 5.0, Y3 - 0.85, "f", "Robustness of the pooled interaction")

axs = []
for j, (k, col, lab) in enumerate([("NetImmune", IMMUNE, "NetImmune"), ("CAF_ECM", STROMA, "CAF_ECM")]):
    ax = S.ax_mm(fig, 32 + j * 27, Y3A, 19, h3)
    axs.append(ax)
    m = meta[meta.score == k].set_index("set")
    for i, (sk, sl, sep) in enumerate(SETS):
        r = m.loc[sk]
        yy = -i
        sig = r.p < 0.05
        ax.plot([r.lo, r.hi], [yy, yy], color=col, lw=0.8, solid_capstyle="butt")
        ax.scatter([r.est], [yy], s=13 if i == 0 else 9, marker="D" if i == 0 else "o",
                   color=col if sig else "white", ec=col, lw=0.7, zorder=3)
        ax.text(1.02, yy, fmtp(r.p), transform=ax.get_yaxis_transform(), fontsize=5.0,
                color=INK if sig else SLATE, fontweight="bold" if sig else "normal", va="center", ha="left")
        if sep and i > 0:
            ax.axhline(yy + 0.5, color=LIGHT, lw=0.4)
    S.rule(ax, 0)
    ax.set_ylim(-len(SETS) + 0.4, 0.6); ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-1.15, 1.3) if k == "NetImmune" else ax.set_xlim(-1.15, 1.3)
    ax.set_xticks([-1, 0, 1]); ax.set_xticklabels(["−1", "0", "+1"])
    ax.text(0.0, 1.0, lab, transform=ax.transAxes, fontsize=5.6, fontweight="bold", va="bottom")
    if j == 0:
        for i, (sk, sl, sep) in enumerate(SETS):
            ax.text(-0.06, -i, sl, transform=ax.get_yaxis_transform(), fontsize=5.0, ha="right", va="center",
                    fontweight="bold" if i == 0 else "normal")
txt(32 + 23, H - 5.0, "Pooled interaction per SD", size=6.0, ha="center")
txt(32 + 19 + 0.6, Y3A - 1.3, "P", size=5.0, color=GREY)
txt(59 + 19 + 0.6, Y3A - 1.3, "P", size=5.0, color=GREY)

# =========================================================================================
# d  pCR by NetImmune tertile in each arm (descriptive)
# =========================================================================================
XD = 83
heading(XD - 4.0, XD, Y2 - 0.85, "d", "Breast: pCR by tertile")



def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


tab = []
axd = []
for j, (t, lab) in enumerate([("NeoTRIP", "NeoTRIP"), ("I-SPY2 pembrolizumab", "I-SPY2 pembro.")]):
    ax = S.ax_mm(fig, XD + 7 + j * 19, Y2A, 15.5, h2)
    axd.append(ax)
    d = tert[tert.trial == t]
    for a, col, mk, dx in [(1, IMMUNE, "o", 0.09), (0, GREY, "s", -0.09)]:
        xs_, ps = [], []
        for q in range(3):
            g = d[(d.arm == a) & (d.tertile == q)]
            k_, n_ = int(g.y.sum()), len(g)
            lo, hi = wilson(k_, n_)
            ax.plot([q + dx, q + dx], [lo * 100, hi * 100], color=col, lw=0.6, solid_capstyle="butt")
            xs_.append(q + dx); ps.append(100 * k_ / n_)
            tab.append(dict(trial=t, arm="checkpoint" if a else "control", tertile=["low", "mid", "high"][q],
                            n=n_, pcr=k_, rate=k_ / n_, wilson_lo=lo, wilson_hi=hi))
        ax.plot(xs_, ps, color=col, lw=0.9, zorder=3)
        ax.scatter(xs_, ps, s=11, marker=mk, color=col, ec="white", lw=0.4, zorder=4)
    ax.set_xlim(-0.45, 2.45); ax.set_ylim(0, 100)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["low", "mid", "high"])
    ax.set_yticks([0, 25, 50, 75, 100]); ax.set_yticklabels(["0", "", "50", "", "100"] if j == 0 else [])
    ax.text(0.0, 1.03, lab, transform=ax.transAxes, fontsize=5.4, fontweight="bold", va="bottom")
    if j == 0:
        ax.set_ylabel("pCR (%)")
txt(XD + 7 + 17, B2 + 6.4, "NetImmune tertile", size=6.0, ha="center")
for j_, (lab_, col_, mk_) in enumerate([("checkpoint arm", IMMUNE, "o"), ("control", GREY, "s")]):
    ky = 95 - j_ * 8
    axd[0].plot([-0.3, 0.05], [ky, ky], color=col_, lw=0.9)
    axd[0].scatter([-0.125], [ky], s=9, marker=mk_, color=col_, ec="white", lw=0.4, zorder=5)
    axd[0].text(0.15, ky, lab_, fontsize=5.0, va="center")

# =========================================================================================
# g  structure-matched null, three added trials (row 3)
# =========================================================================================
XE = 100
heading(XE - 4.0, XE, Y3 - 0.85, "g", "Random-gene reference")

nul = pd.read_csv(R / "E21_null_summary.csv")
XEA, WEA = 128, 38
axe = S.ax_mm(fig, XEA, Y3A, WEA, h3)
rows_e = []
ORD = [("CM025", "CheckMate 025"), ("JR101", "JAVELIN Renal 101"), ("NeoTRIP", "NeoTRIP"), ("combined", "Pooled (3)")]
yy = 0
for k, col in [("NetImmune", IMMUNE), ("CAF_ECM", STROMA)]:
    axe.text(-0.02, yy + 0.15, k, transform=axe.get_yaxis_transform(), fontsize=5.4, fontweight="bold",
             ha="right", va="center")
    yy -= 1
    for c, lab in ORD:
        r = nul[(nul.score == k) & (nul.contrast == c)].iloc[0]
        axe.barh(yy, r.null_p95_abs, height=0.62, color=LIGHT, lw=0)
        axe.plot([r.null_median_abs] * 2, [yy - 0.31, yy + 0.31], color=GREY, lw=0.6)
        sig = r.p_emp < 0.05
        axe.scatter([abs(r.obs_y)], [yy], s=12, marker="D", color=col, ec="white", lw=0.4, zorder=3)
        axe.text(-0.04, yy, lab, transform=axe.get_yaxis_transform(), fontsize=5.0,
                 ha="right", va="center")
        axe.text(1.03, yy, fmtp(r.p_emp), transform=axe.get_yaxis_transform(), fontsize=5.0,
                 color=INK if sig else SLATE, fontweight="bold" if sig else "normal", ha="left", va="center")
        rows_e.append(dict(score=k, contrast=c, obs_abs=abs(r.obs_y), null_median=r.null_median_abs,
                           null_p95=r.null_p95_abs, p_emp=r.p_emp))
        yy -= 1
    yy -= 0.5
axe.set_ylim(yy + 0.5, 0.7); axe.set_yticks([]); axe.spines["left"].set_visible(False)
axe.set_xlim(0, 0.7); axe.set_xticks([0, 0.2, 0.4, 0.6]); axe.set_xticklabels(["0", "0.2", "0.4", "0.6"])
txt(XEA + WEA / 2, H - 5.0, "Absolute interaction per SD", size=6.0, ha="center")
txt(XEA + WEA + 0.6, Y3A - 1.3, "P", size=5.0, color=GREY)

S.save(fig, "Fig6")
pd.DataFrame(cells).to_csv(S.OUT / "data/fig6b_cells.csv", index=False)
pd.DataFrame(tab).to_csv(S.OUT / "data/fig6d_tertiles_summary.csv", index=False)
pd.DataFrame(rows_e).to_csv(S.OUT / "data/fig6g_null.csv", index=False)
pd.DataFrame(rows_c).to_csv(S.OUT / "data/fig6c_within_arm.csv", index=False)
pd.DataFrame(km_rows).to_csv(S.OUT / "data/fig6e_km_summary.csv", index=False)
print("saved; matrix ends", YAB)
