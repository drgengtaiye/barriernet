"""E32 sensitivity - the 52 companion signatures scored with their published engines (GSVA, ssGSEA, weighted mean), as in
the companion's score_matrices.R (gsvaParam / ssgseaParam, minSize 3; weighted mean of row-z with published weights,
missing weights 0), in each of the six comparisons; interaction with the E32 primary model (I-SPY2 adjusted for
hormone receptor and MammaPrint) and the manuscript pooling rule. Protocol: E32_prereg.md (sensitivity, no expectation).
Expression passed to R = each comparison's log matrix as used by E28_base (same patients); the R step is
E32_native_scores.R. Output: results/E32_native.csv.
"""
import hashlib, subprocess, os
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h = hashlib.sha256(open(f"{A}/E32_prereg.md", "rb").read()).hexdigest()
assert _h == open(f"{A}/E32_prereg.sha256").read().split()[0], "protocol hash mismatch"
_src = open(f"{A}/E28_base.py").read().replace('f"{O}/E28_reproduction_check.csv"', 'f"{O}/E32_native_step0.csv"')
exec(compile(_src, f"{A}/E28_base.py", "exec"))
G = "/groups/adv2105_gp/yichen/Yi/com/gastric"
T = f"{O}/E32_native_tmp"; os.makedirs(T, exist_ok=True)


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
BIO = pd.read_csv(f"{G}/validation/raw/GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz"); BIO["ResearchID"] = BIO.ResearchID.astype(str)
mp = {"ISPY2_pembro": pd.to_numeric(S194.loc[TR["ISPY2_pembro"]["Z"].columns.astype(str), "mp"], errors="coerce").values,
      "ISPY2_durva": pd.to_numeric(BIO.drop_duplicates("ResearchID").set_index("ResearchID").loc[
          TR["ISPY2_durva"]["Z"].columns.astype(str), "MP.1.2..class..0..MP1..1..MP2."], errors="coerce").values}
for k, v in mp.items():
    TR[k]["cov"] = np.column_stack([TR[k]["cov"], v])

for k in ALL6:
    X = TR[k]["X"].loc[:, TR[k]["Z"].columns]
    X = X[np.isfinite(X.to_numpy(float)).all(1)]                       # GSVA needs complete rows
    X.to_csv(f"{T}/{k}.tsv.gz", sep="\t", compression="gzip")
r = subprocess.run(["Rscript", f"{A}/E32_native_scores.R", T] + ALL6, capture_output=True, text=True)
print(r.stdout[-3000:], r.stderr[-2000:]); assert r.returncode == 0

rows = []
for k in ALL6:
    S = pd.read_csv(f"{T}/{k}.scores.csv", index_col=0)               # samples x signatures
    S.index = S.index.astype(str); S = S.loc[TR[k]["Z"].columns.astype(str)]
    for s in S.columns:
        sc = S[s].to_numpy(float)
        if np.isfinite(sc).all() and np.std(sc) > 0:
            y, se = m_int(TR[k], sc)
            rows.append(dict(signature=s, comparison=k, y=y, se=se))
PER = pd.DataFrame(rows)
FAM = pd.read_csv(f"{G}/icb/v5_signature_clusters.csv").set_index("signature").family
out = []
for s, g in PER.groupby("signature"):
    g = g[np.isfinite(g.y) & np.isfinite(g.se)]
    if len(g) >= 4:
        p = pool_rule(g.y.values, g.se.values)
        out.append(dict(signature=s, family=FAM.get(s, ""), k=len(g), est=p["est"], lo=p["lo"], hi=p["hi"], p=p["p"],
                        k_pos=int((g.y > 0).sum())))
OUT = pd.DataFrame(out).sort_values("p")
PER.to_csv(f"{O}/E32_native_per_comparison.csv", index=False); OUT.to_csv(f"{O}/E32_native.csv", index=False)
pd.set_option("display.width", 200)
print(OUT.round(4).to_string(index=False))
print(f"\npublished engines: {len(OUT)} pooled; P < 0.05: {int((OUT.p < 0.05).sum())}, all favoring: {bool((OUT[OUT.p < 0.05].est > 0).all())}")
