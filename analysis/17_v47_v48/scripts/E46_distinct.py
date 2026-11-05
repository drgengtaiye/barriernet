"""E46 part E - every count of signatures in the main text, as 52 rows and as 51 distinct signatures (TIS_Damotte removed; it has the
same 18 genes as Tcell-inflamed_Ayers). No new model. Protocol: analysis/E46_prereg.md (hash checked). Writes results/E46_E_distinct.csv."""
import hashlib, numpy as np, pandas as pd
from scipy.stats import spearmanr
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"; O = f"{A}/results"
_h = hashlib.sha256(open(f"{A}/E46_prereg.md", "rb").read()).hexdigest(); assert _h == open(f"{A}/E46_prereg.sha256").read().split()[0]
DUP, AY = "TIS_Damotte", "Tcell-inflamed_Ayers"; ROW = []


def add(item, full, dist, source):
    ROW.append(dict(item=item, rows_52=full, distinct_51=dist, source=source))


def both(item, df, name, f, source):
    add(item, f(df), f(df[df[name] != DUP]), source)


P = pd.read_csv(f"{O}/E32_primary.csv"); L = P[P.kind == "library"].copy(); assert len(L) == 52
a, b = L.set_index("signature").loc[DUP], L.set_index("signature").loc[AY]; assert abs(a.est - b.est) < 1e-10 and a.family == b.family and a.declared_dir == b.declared_dir
print("TIS and Ayers: est", round(a.est, 4), "family", a.family, "| emp P", a.emp_p_declared, b.emp_p_declared, "| adj p", a.p_estimate_adj, b.p_estimate_adj)
s = "E32_primary.csv"
both("signatures", L, "signature", len, s)
both("sensitivity label", L, "signature", lambda d: int((d.declared_dir > 0).sum()), s); both("resistance label", L, "signature", lambda d: int((d.declared_dir < 0).sum()), s)
for fam in L.family.unique():
    both(f"family {fam}", L, "signature", lambda d, fam=fam: int((d.family == fam).sum()), s)
both("favor P<0.05", L, "signature", lambda d: int(((d.p < 0.05) & (d.est > 0)).sum()), s); both("against P<0.05", L, "signature", lambda d: int(((d.p < 0.05) & (d.est < 0)).sum()), s)
both("significant any direction", L, "signature", lambda d: int((d.p < 0.05).sum()), s)
both("scored in all six", L, "signature", lambda d: int((d.k == 6).sum()), s); both("favor P<0.05 among all six", L, "signature", lambda d: int(((d.k == 6) & (d.p < 0.05) & (d.est > 0)).sum()), s)
both("not all six", L, "signature", lambda d: int((d.k < 6).sum()), s); both("favor P<0.05 among not all six", L, "signature", lambda d: int(((d.k < 6) & (d.p < 0.05) & (d.est > 0)).sum()), s)
INF = L.family.isin(["Inflamed/IFN"])
both("inflamed est>0", L[INF], "signature", lambda d: int((d.est > 0).sum()), s)
both("exceed random ref", L, "signature", lambda d: int((d.emp_p_declared < 0.05).sum()), s)
for fam in L.family.unique():
    both(f"exceed random ref {fam}", L, "signature", lambda d, fam=fam: int(((d.emp_p_declared < 0.05) & (d.family == fam)).sum()), s)
both("exceed random ref, not nominal", L, "signature", lambda d: int(((d.emp_p_declared < 0.05) & (d.p >= 0.05)).sum()), s)
both("median of ref medians", L, "signature", lambda d: round(float(d.ref_median.median()), 3), s)
both("positive diff from CD8A", L, "signature", lambda d: int((d.diff_vs_CD8A > 0).sum()), s); both("above CD8A P<0.05", L, "signature", lambda d: int(((d.diff_vs_CD8A > 0) & (d.cd8a_p < 0.05)).sum()), s)
both("below CD8A P<0.05", L, "signature", lambda d: int(((d.diff_vs_CD8A < 0) & (d.cd8a_p < 0.05)).sum()), s); both("below CD8A P<0.05, resistance label", L, "signature", lambda d: int(((d.diff_vs_CD8A < 0) & (d.cd8a_p < 0.05) & (d.declared_dir < 0)).sum()), s)
both("median paired SE proxy max hi", L, "signature", lambda d: round(float(d.cd8a_hi.max()), 3), s); both("median cd8a_hi", L, "signature", lambda d: round(float(d.cd8a_hi.median()), 3), s)
sg = L[L.p < 0.05]; both("significant, ESTIMATE-adjusted P<0.05", sg, "signature", lambda d: int((d.p_estimate_adj < 0.05).sum()), s)
# coexpression
C = pd.read_csv(f"{O}/E45_B_coexpression.csv"); s = "E45_B_coexpression.csv"
both("exceed coexpression ref", C, "signature", lambda d: int((d.emp_p_coex < 0.05).sum()), s)
for fam in C.family.unique():
    both(f"exceed coexpression ref {fam}", C, "signature", lambda d, fam=fam: int(((d.emp_p_coex < 0.05) & (d.family == fam)).sum()), s)
both("lost by coexpression", C, "signature", lambda d: int(((d.emp_p_size < 0.05) & (d.emp_p_coex >= 0.05)).sum()), s); both("gained by coexpression", C, "signature", lambda d: int(((d.emp_p_size >= 0.05) & (d.emp_p_coex < 0.05)).sum()), s)
both("median coexpression ref median", C, "signature", lambda d: round(float(d.ref_median_coex.median()), 3), s)
# single-arm concordance
E = pd.read_csv(f"{O}/E40_concordance.csv"); s = "E40_concordance.csv"
both("Spearman single-arm vs interaction", E, "signature", lambda d: round(float(spearmanr(d.sa_est, d.int_est).statistic), 3), s)
both("interaction significant", E, "signature", lambda d: int((d.int_p < 0.05).sum()), s); both("interaction significant and single-arm significant", E, "signature", lambda d: int(((d.int_p < 0.05) & (d.sa_p < 0.05)).sum()), s)
both("single-arm significant, interaction not", E, "signature", lambda d: int(((d.int_p >= 0.05) & (d.sa_p < 0.05)).sum()), s); both("checkpoint arm P<0.05", E, "signature", lambda d: int((d.icb_p < 0.05).sum()), s)
both("checkpoint arm P<0.05, interaction not", E, "signature", lambda d: int(((d.icb_p < 0.05) & (d.int_p >= 0.05)).sum()), s); both("control arm P<0.05", E, "signature", lambda d: int((d.ctl_p < 0.05).sum()), s)
# sensitivity sets
AS = pd.read_csv(f"{O}/E32_pooled_all_sets.csv"); AS = AS[AS.kind == "library"]; s = "E32_pooled_all_sets.csv"
for st, q in AS.groupby("set", sort=False):
    both(f"{st}: favor", q, "signature", lambda d: int(((d.p < 0.05) & (d.est > 0)).sum()), s); both(f"{st}: against", q, "signature", lambda d: int(((d.p < 0.05) & (d.est < 0)).sum()), s)
NV = pd.read_csv(f"{O}/E32_native.csv"); s = "E32_native.csv (TIS and Ayers may differ under published engines)"
print("native engines: ", NV[NV.signature.isin([DUP, AY])][["signature", "est", "p"]].round(4).to_dict("records"))
both("published engines: favor", NV, "signature", lambda d: int(((d.p < 0.05) & (d.est > 0)).sum()), s); both("published engines: against", NV, "signature", lambda d: int(((d.p < 0.05) & (d.est < 0)).sum()), s)
# panel genes
PQ = pd.read_csv(f"{O}/E33_paired.csv"); print("E33_paired columns", PQ.columns.tolist()); s = "E33_paired.csv"
Q6 = PQ[PQ.set == "six (primary)"]; refcol = [c for c in PQ.columns if c in ("gene", "reference", "name")][0]; assert Q6.signature.nunique() == 52
lib = set(L.signature); Q6 = Q6[Q6.signature.isin(lib)]; TP = ["CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "CD2", "CD247", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "NKG7"]
for nm, q in [("rows_52", Q6), ("distinct_51", Q6[Q6.signature != DUP])]:
    ab = q[(q["diff"] > 0) & (q.p < 0.05)].groupby(refcol).size().reindex(q[refcol].unique(), fill_value=0); t = ab.reindex(TP, fill_value=0)
    print(nm, "T-cell genes signatures above:", t.to_dict()); print(nm, "all panel signatures above:", ab.sort_values(ascending=False).to_dict())
    ROW.append(dict(item=f"T-cell genes exceeded by at most 2 ({nm})", rows_52=int((t <= 2).sum()) if nm == "rows_52" else None, distinct_51=int((t <= 2).sum()) if nm != "rows_52" else None, source=s))
    ROW.append(dict(item=f"T-cell genes exceeded by none ({nm})", rows_52=int((t == 0).sum()) if nm == "rows_52" else None, distinct_51=int((t == 0).sum()) if nm != "rows_52" else None, source=s))
B3 = PQ[(PQ.set == "breast three") & PQ.signature.isin(lib)] if (PQ.set == "breast three").any() else None
GP = pd.read_csv(f"{O}/E33_pooled.csv"); BL = GP[(GP.set == "breast three") & (GP.kind == "library")]; s = "E33_pooled.csv"
both("breast: favor", BL, "name", lambda d: int(((d.p < 0.05) & (d.est > 0)).sum()), s); both("breast: est>0", BL, "name", lambda d: int((d.est > 0).sum()), s)
print("Ayers/TIS in breast: ", BL[BL.name.isin([DUP, AY])][["name", "est", "p"]].round(3).to_dict("records"))
BX = pd.read_csv(f"{O}/E33_breast_extras.csv"); print("E33_breast_extras columns", BX.columns.tolist())
# references
PA = pd.read_csv(f"{O}/E45_A_paired.csv"); s = "E45_A_paired.csv"
for st in PA.set.unique():
    q = PA[(PA.set == st) & (PA.reference == "CXCL10") & (~PA.score.isin(["NetImmune", "barriernet 0.3"]))]; assert len(q) == 52, len(q)
    both(f"{st}: above CXCL10", q, "score", lambda d: int(d.above.sum()), s); both(f"{st}: below CXCL10", q, "score", lambda d: int(d.below.sum()), s)
    q = PA[(PA.set == st) & (PA.reference == "CD8A") & (~PA.score.isin(["NetImmune", "barriernet 0.3"]))]
    both(f"{st}: above CD8A", q, "score", lambda d: int(d.above.sum()), s); both(f"{st}: below CD8A", q, "score", lambda d: int(d.below.sum()), s)
# gastric
GB = pd.read_csv(f"{O}/E44_B_pooled.csv"); s = "E44_B_pooled.csv"; G3 = GB[(GB.set.str.startswith("three")) & (GB.kind == "sig") & GB.name.isin(lib)]; assert len(G3) == 52
both("gastric three: fixed P<0.05", G3, "name", lambda d: int((d.p < 0.05).sum()), s); both("gastric three: rule P<0.05", G3, "name", lambda d: int((d.cons_p < 0.05).sum()), s)
both("gastric three: significant fixed and rule", G3, "name", lambda d: int(((d.cons_p < 0.05) & (d.p < 0.05)).sum()), s)
L2 = L.set_index("signature")
for nm, q in [("rows_52", G3), ("distinct_51", G3[G3.name != DUP])]:
    print(nm, "gastric vs randomized Spearman", round(float(spearmanr(q.est, L2.loc[q.name, "est"]).statistic), 3))
GPp = pd.read_csv(f"{O}/E44_B_paired.csv"); print("E44_B_paired columns", GPp.columns.tolist()); s = "E44_B_paired.csv"
gq = GPp[GPp.set.str.startswith("three") & (GPp.name if "name" in GPp else GPp.score).isin(lib)] if "set" in GPp else None
if gq is not None:
    nmc = "name" if "name" in gq else "score"; dc = "diff" if "diff" in gq else "diff_vs_CD8A"; pc = "p" if "p" in gq else "diff_p"
    both("gastric three: above CD8A", gq, nmc, lambda d: int(((d[dc] > 0) & (d[pc] < 0.05)).sum()), s); both("gastric three: below CD8A", gq, nmc, lambda d: int(((d[dc] < 0) & (d[pc] < 0.05)).sum()), s)
GA = pd.read_csv(f"{O}/E44_B_adjusted.csv"); GA = GA[(GA.kind == "sig") & GA.name.isin(lib)]; s = "E44_B_adjusted.csv"; assert len(GA) == 52
both("two gastric: P<0.05 unadjusted", GA, "name", lambda d: int((d.p < 0.05).sum()), s); both("two gastric: P<0.05 adjusted", GA, "name", lambda d: int((d.p_msiebv < 0.05).sum()), s); both("two gastric: adjusted est>0", GA, "name", lambda d: int((d.est_msiebv > 0).sum()), s)
q = GA[GA.p < 0.05]; both("two gastric: median kept fraction among significant", q, "name", lambda d: round(float((d.est_msiebv / d.est).median()), 3), s)
print("TIS gastric: ", GA[GA.name.isin([DUP, AY])].round(3).to_dict("records"))
R = pd.DataFrame(ROW); R.to_csv(f"{O}/E46_E_distinct.csv", index=False); pd.set_option("display.width", 200, "display.max_rows", 200, "display.max_colwidth", 70); print(R[["item", "rows_52", "distinct_51"]].to_string(index=False))
