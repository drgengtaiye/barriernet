"""E34 linking - NeoPembrOv patient identifiers (Le Saux et al. 2024 source data) to GEO GSE227666 library identifiers,
by expression values only (no outcome column is used).

Run before E34_prereg.md was written (the interactive version of this code produced id_match_expression_only.csv and
id_match_final.csv on 2026-10-25); this script reproduces both files and checks that the assignments are unchanged.
- Pembrolizumab arm (NACT+P): 40 genes reported per patient in the source data (log2 TPM), z-scored across patients in
  each dataset, Pearson correlation, one-to-one assignment (Hungarian).
- Chemotherapy arm (NACT): the 10 genes reported per patient for both arms (TIGIT, FOXP3, HAVCR2, ENTPD1, TREM2, SMARCD3,
  CXCL13, EPCAM, PTPRC, JCHAIN), rank-transformed, Pearson correlation of ranks, one-to-one assignment. The same 10 genes
  are applied to the pembrolizumab arm as a check against the 40-gene links.
"""
import gzip, sys
import numpy as np, pandas as pd, openpyxl
from scipy.optimize import linear_sum_assignment

D = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/validation/NeoPembrOV_GSE227666"
titles = [l.rstrip("\n").split("\t")[1:] for l in gzip.open(f"{D}/GSE227666_series_matrix.txt.gz", "rt") if l.startswith("!Sample_title")][0]
meta = pd.DataFrame({"title": [t.strip('"') for t in titles]})
meta["rid"] = meta.title.str.extract(r"\((R\d+)\)"); meta["arm"] = meta.title.str.extract(r", (NACT\+P|NACT), ")
meta["time"] = meta.title.str.extract(r", (Pre|Post) ")
pre = meta[meta.time == "Pre"]
C = pd.read_csv(f"{D}/GSE227666_neopembrov_counts.txt.gz", sep="\t", index_col=0); C = C.groupby(level=0).sum()
cpm = C / C.sum() * 1e6
wb = openpyxl.load_workbook(f"{D}/nc2024_source_data.xlsx", read_only=True)


def sheet(n):
    h, *rs = list(wb[n].iter_rows(values_only=True)); return pd.DataFrame(rs, columns=h)


arm_src = sheet("Figures 4 a,b,c,d").set_index("sample_name")["Treatment_arm"].to_dict()     # identifiers and arm only


def assign(M, ss, gs, arm, method):
    r_, c_ = linear_sum_assignment(-np.nan_to_num(M, nan=-1))
    out = []
    for i, j in zip(r_, c_):
        row = np.sort(np.nan_to_num(M[i], nan=-1))[::-1]
        out.append(dict(arm=arm, src=ss[i], geo=gs[j], r=M[i, j], second=row[1], method=method))
    return pd.DataFrame(out)


# ---- pembrolizumab arm, 40 genes
s6 = sheet("Supp Figures 6 b,d"); s6["log2_tpm"] = pd.to_numeric(s6.log2_tpm, errors="coerce")
S40 = s6.pivot_table(index="SUBJID", columns="symbol", values="log2_tpm")
s5 = sheet("Supp Figures 5 a,c,d").set_index("sample_name")
for g in ("TIGIT", "FOXP3", "HAVCR2"):
    S40[g] = pd.to_numeric(s5[f"{g}_log2_tpm"], errors="coerce").reindex(S40.index)
g40 = [g for g in S40.columns if g in cpm.index]
gsP = pre[pre.arm == "NACT+P"].rid.tolist(); ssP = [s for s in S40.index if arm_src.get(s) == "NACT+P"]
A = np.log2(cpm.loc[g40, gsP].T + 1); B = S40.loc[ssP, g40]
za, zb = (A - A.mean()) / A.std(), (B - B.mean()) / B.std()
M = np.array([[np.corrcoef(zb.loc[s][zb.loc[s].notna() & za.loc[g].notna()], za.loc[g][zb.loc[s].notna() & za.loc[g].notna()])[0, 1]
               for g in gsP] for s in ssP])
P40 = assign(M, ssP, gsP, "NACT+P", f"{len(g40)} genes")

# ---- 10 genes reported for both arms (ranks)
L = lambda x: np.log2(x + 0.5)
F = {}
for g in ("TIGIT", "FOXP3", "HAVCR2"):
    F[g] = pd.to_numeric(s5[f"{g}_log2_tpm"], errors="coerce")
f5 = sheet("Figures 5 a,b,c").set_index("sample_name"); F["ENTPD1"] = pd.to_numeric(f5["ENTPD1_tpm"], errors="coerce")
f4 = sheet("Figures 4 a,b,c,d").set_index("sample_name")
for g in ("TREM2", "SMARCD3"):
    F[g] = pd.to_numeric(f4[g], errors="coerce")
f1 = sheet("Figure 1 a,b"); f1 = f1[f1.time == "Pre"].copy(); f1["p"] = f1.sample_name.str.replace("_Pre", "")
for g, lab in [("CXCL13", "CXCL13"), ("EPCAM", "EPCAM"), ("PTPRC", "PTPRC (CD45)"), ("JCHAIN", "JCHAIN")]:
    F[g] = pd.to_numeric(f1[f1.immune_pop == lab].set_index("p")["score"], errors="coerce")
G10 = list(F)


def match10(arm):
    gs = pre[pre.arm == arm].rid.tolist(); ss = [s for s, a in arm_src.items() if a == arm]
    Aa = pd.DataFrame({g: L(cpm.loc[g, gs]).values for g in G10}, index=gs)
    Bb = pd.DataFrame({g: F[g].reindex(ss).values for g in G10}, index=ss)
    ra, rb = Aa.rank() / len(Aa), Bb.rank() / Bb.notna().sum()
    M = np.array([[np.corrcoef(rb.loc[s][rb.loc[s].notna()], ra.loc[g][rb.loc[s].notna()])[0, 1] for g in gs] for s in ss])
    return assign(M, ss, gs, arm, "10 genes (ranks)")


C10P, C10N = match10("NACT+P"), match10("NACT")
agree = (C10P.set_index("src").geo == P40.set_index("src").geo).sum()
best = sum(1 for _, r in C10P.iterrows() if np.isclose(r.r, max(r.r, r.second)))
print(f"pembrolizumab arm: 40-gene links r {P40.r.min():.2f}-{P40.r.max():.2f}, second-best max {P40.second.max():.2f}, "
      f"margin min {(P40.r - P40.second).min():.2f}; 10-gene check reproduces {agree}/{len(P40)} links")
print(f"chemotherapy arm (10 genes): r {C10N.r.min():.2f}-{C10N.r.max():.2f}, second-best max {C10N.second.max():.2f}, "
      f"margin min {(C10N.r - C10N.second).min():.3f}; links with margin < 0.10: {int(((C10N.r - C10N.second) < 0.10).sum())}")
final = pd.concat([P40, C10N], ignore_index=True)
old = pd.read_csv(f"{D}/id_match_final.csv").set_index("src").geo
same = (final.set_index("src").geo == old.reindex(final.src).values).all()
print("assignments identical to the files used by E34:", bool(same)); assert same
final.to_csv(f"{D}/id_match_final.csv", index=False)
pd.concat([P40, C10P.assign(method="10 genes (ranks), check")]).to_csv(f"{D}/id_match_checks.csv", index=False)
