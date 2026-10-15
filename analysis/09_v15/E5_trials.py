"""E5 (trials) - score x arm interaction in randomised trials, outside gastric cancer.

Primary: IMmotion150 (renal cell carcinoma; atezolizumab +/- bevacizumab vs sunitinib),
iAtlas expression and cBioPortal clinical data as used in pancan/v4_arm_split.py. Contrast
as published: both atezolizumab arms pooled vs sunitinib; atezolizumab monotherapy vs
sunitinib as a secondary contrast (the only one without anti-VEGF; smallest).
Secondary: PRINCE (pancreatic cancer; nivolumab-containing arms vs sotigalimab arm, all with
gemcitabine/nab-paclitaxel) - not a checkpoint-vs-no-immunotherapy contrast, so reported
for completeness only.

Model: Cox PFS ~ s + arm + s:arm, s = score standardised to SD 1 over the analysed patients;
interaction HR = exp(coef of s:arm), < 1 meaning the score predicts benefit specifically in
the checkpoint arm. Within-arm HR per SD from separate Cox fits.

Pre-stated expectation (written before running), from the gastric results:
  TAM_SUP, IFN, Immune   - interaction HR < 1 (benefit direction), IFN/Immune the strongest
  CAF_ECM                - little interaction (HR near 1); adverse (HR > 1) or null in both arms
  Barrier                - attenuated interaction, halves pulling in opposite directions
  NetImmune              - interaction HR < 1, no stronger than Immune
Peer-review addition (defined before any outcome for them was computed):
  C1Q_arm  (C1QA, C1QB, C1QC)                 - interaction HR < 1
  SPP1_arm (SPP1, TREM2, MARCO, APOE, LGALS3) - no benefit-direction interaction, or the opposite
  CXCL9_SPP1 polarity = z(CXCL9) - z(SPP1)    - interaction HR < 1 (reference)
"""
import warnings, functools; warnings.filterwarnings("ignore")
print = functools.partial(print, flush=True)
import numpy as np, pandas as pd
from lifelines import CoxPHFitter

P = "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"
O = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
MOD = dict(
    CAF_ECM=["COL1A1", "COL1A2", "LUM", "DCN", "POSTN", "TAGLN", "ACTA2", "THY1", "FAP", "FN1", "COL3A1"],
    TAM_SUP=["SPP1", "APOE", "TREM2", "C1QA", "C1QB", "C1QC", "LGALS3", "MARCO", "IL10", "SOCS3", "STAT3"],
    APM=["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "PSMB8", "PSMB9"],
    IFN=["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG"],
    TLS_B=["CXCL13", "LTB", "MS4A1", "CD79A", "CD74", "MZB1"],
    T_CELL=["CD3D", "CD3E", "TRAC", "TRBC1", "CD8A", "GZMB", "NKG7"],
    C1Q_arm=["C1QA", "C1QB", "C1QC"],
    SPP1_arm=["SPP1", "TREM2", "MARCO", "APOE", "LGALS3"])
ORDER = ["CAF_ECM", "TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1", "IFN", "APM", "TLS_B", "T_CELL",
         "Barrier", "Immune", "NetImmune"]


def scores(Z):
    S = {}
    for k, G in MOD.items():
        g = [x for x in G if x in Z.index]
        S[k] = Z.loc[g].mean(0) if len(g) >= 3 else pd.Series(np.nan, index=Z.columns)
    S = pd.DataFrame(S)
    S["Barrier"] = S[["CAF_ECM", "TAM_SUP"]].mean(1)
    S["Immune"] = S[["APM", "IFN", "TLS_B", "T_CELL"]].mean(1)
    S["NetImmune"] = S.Immune - S.Barrier
    S["CXCL9_SPP1"] = Z.loc["CXCL9"] - Z.loc["SPP1"]
    return S


def pfs(C):
    t = pd.to_numeric(C["PFS_MONTHS"], errors="coerce")
    e = pd.to_numeric(C["PFS_STATUS"].astype(str).str.extract(r"^\s*([01])\s*:")[0], errors="coerce").fillna(0).astype(int)
    return t, e


def fit(s, arm, t, e):
    s = (s - s.mean()) / s.std()
    df = pd.DataFrame({"T": t.values, "E": e.values, "s": s.values, "arm": arm.values})
    df["sx"] = df.s * df.arm
    f = CoxPHFitter().fit(df, "T", "E")
    r = f.summary.loc["sx"]
    out = dict(int_HR=np.exp(r["coef"]), int_lo=np.exp(r["coef lower 95%"]), int_hi=np.exp(r["coef upper 95%"]), int_p=r["p"])
    for a, lab in [(1, "ici"), (0, "ctrl")]:
        d = df[df.arm == a][["T", "E", "s"]]
        g = CoxPHFitter().fit(d, "T", "E").summary.loc["s"]
        out.update({f"{lab}_HR": np.exp(g["coef"]), f"{lab}_lo": np.exp(g["coef lower 95%"]),
                    f"{lab}_hi": np.exp(g["coef upper 95%"]), f"{lab}_p": g["p"]})
    return out


rows = []
# ------------------------------------------------------------------ IMmotion150
C = pd.read_csv(f"{P}/v4_immotion150_clinical_api.csv", index_col=0)
M = pd.read_csv(f"{P}/ici/rcc_iatlas_immotion150_2018.csv", index_col=0)
keep = [s for s in M.columns if s in C.index]; C, M = C.loc[keep], M[keep]
Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
S = scores(Z); t, e = pfs(C)
atezo = C["ICI_RX"].astype(str) == "Atezolizumab"
bev = C["NON_ICI_RX"].astype(str).str.contains("Bevacizumab", na=False)
suni = C["NON_ICI_RX"].astype(str).str.contains("Sunitinib", na=False)
for label, (tr, ct) in {"IMmotion150: atezolizumab (+/- bev) vs sunitinib": (atezo, suni),
                        "IMmotion150: atezolizumab mono vs sunitinib": (atezo & ~bev, suni)}.items():
    sel = (tr | ct) & t.notna() & (t > 0)
    arm = tr[sel].astype(int)
    print(f"{label}: checkpoint n={arm.sum()}, control n={(1 - arm).sum()}, events={int(e[sel].sum())}")
    for k in ORDER:
        rows.append(dict(trial=label, score=k, n_ici=int(arm.sum()), n_ctrl=int((1 - arm).sum()),
                         events=int(e[sel].sum()), **fit(S.loc[sel[sel].index, k], arm, t[sel], e[sel])))
# ------------------------------------------------------------------ PRINCE (secondary)
C = pd.read_csv(f"{P}/v2_prince_clinical_api.csv", index_col=0)
M = pd.read_csv(f"{P}/ici/paad_iatlas_prince_2022.csv", index_col=0)
keep = [s for s in M.columns if s in C.index]; C, M = C.loc[keep], M[keep]
Z = ((M.T - M.mean(1)) / M.std(1).replace(0, 1)).T
S = scores(Z); t, e = pfs(C)
nivo = C["ICI_RX"].astype(str).str.contains("Nivolumab", case=False, na=False)
sel = t.notna() & (t > 0); arm = nivo[sel].astype(int)
label = "PRINCE: nivolumab-containing vs sotigalimab (all chemo)"
print(f"{label}: checkpoint n={arm.sum()}, control n={(1 - arm).sum()}, events={int(e[sel].sum())}")
for k in ORDER:
    rows.append(dict(trial=label, score=k, n_ici=int(arm.sum()), n_ctrl=int((1 - arm).sum()),
                     events=int(e[sel].sum()), **fit(S.loc[sel[sel].index, k], arm, t[sel], e[sel])))

R = pd.DataFrame(rows)
R.to_csv(f"{O}/E5_trials.csv", index=False)
pd.set_option("display.width", 200)
for tr, g in R.groupby("trial", sort=False):
    print(f"\n=== {tr} ===")
    print(g[["score", "int_HR", "int_lo", "int_hi", "int_p", "ici_HR", "ici_p", "ctrl_HR", "ctrl_p"]].round(3).to_string(index=False))
