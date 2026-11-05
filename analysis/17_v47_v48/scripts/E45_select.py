"""E45 step 2 - choose the reference genes in the 20 single-arm cohorts and freeze them. Protocol: analysis/E45_prereg.md (hash
checked). Setup: E34_trials.py up to its fitting stage, for lfit, std and pool_rule only; no randomized result is read or fitted
here apart from the list of candidate genes. Step 0: Firth log-odds of CD8A per cohort from this export equal those from the E40
export. Writes results/E45_single_arm_gene_pooled.csv and results/E45_reference_genes.tsv."""
import hashlib, sys
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"
_h45 = hashlib.sha256(open(f"{A}/E45_prereg.md", "rb").read()).hexdigest()
assert _h45 == open(f"{A}/E45_prereg.sha256").read().split()[0], "E45 protocol hash mismatch"
_s34 = open(f"{A}/E34_trials.py").read().split("\n\ndef fit(task):")[0]
_s34 = _s34.replace("E34_step0_reproduction", "E45_step0_reproduction").replace("E34_coverage_check.csv", "E45_coverage_check.csv")
exec(compile(_s34, f"{A}/E34_trials.py", "exec"))


class Tee:
    def __init__(self, *f): self.f = f
    def write(self, s): [x.write(s) for x in self.f]
    def flush(self): [x.flush() for x in self.f]


sys.stdout = Tee(sys.__stdout__, open(f"{O}/E45_select.log", "w"))
print("E45 protocol sha256 verified:", _h45[:12])
SA = pd.read_csv(f"{O}/E45_single_arm_genes.csv.gz")
fl = pd.read_csv(f"{O}/E37_twin_flags.csv")
drop = SA.set_index(["cohort", "sample"]).index.isin(list(zip(fl.development, fl["sample"]))); assert drop.sum() == 3
SA = SA[~drop]; print(f"single-arm: {SA.cohort.nunique()} cohorts, {len(SA)} patients, {int(SA.response.sum())} responders")
assert SA.cohort.nunique() == 20 and len(SA) == 1195 and int(SA.response.sum()) == 368
P33 = pd.read_csv(f"{O}/E33_pooled.csv")
C = P33[(P33.set == "six (primary)") & (P33.kind == "gene") & (P33.all_six.astype(str) == "True")].name.tolist(); assert len(C) == 2082
PANEL = T_PANEL + I_PANEL; assert len(PANEL) == 27
COH = list(SA.cohort.unique())
Y = {c: SA[SA.cohort == c].response.to_numpy(float) for c in COH}
GZ = {c: SA[SA.cohort == c] for c in COH}


def one(g):
    out = []
    for c in COH:
        x = GZ[c][g].to_numpy(float)
        if not np.isfinite(x).all() or np.std(x) == 0:
            continue
        try:
            b, se = lfit(np.column_stack([np.ones(len(x)), std(x)]), Y[c], [1], True)[0]
        except Exception:
            b, se = np.nan, np.nan
        out.append((c, b, se))
    return g, out


with mp.get_context("fork").Pool(NPROC) as pool:
    res = pool.map(one, C, chunksize=8)
rows = []
for g, out in res:
    ys = np.array([o[1] for o in out]); ses = np.array([o[2] for o in out]); ok = np.isfinite(ys) & np.isfinite(ses)
    r = pool_rule(ys[ok], ses[ok]) if ok.sum() >= 2 else dict(est=np.nan, lo=np.nan, hi=np.nan, p=np.nan, k=int(ok.sum()))
    rows.append(dict(gene=g, k=int(ok.sum()), k_pos=int((ys[ok] > 0).sum()), est=r["est"], lo=r["lo"], hi=r["hi"], p=r["p"], panel=g in PANEL))
T = pd.DataFrame(rows); T["eligible"] = T.k >= 16
T["rank"] = T.est.where(T.eligible).rank(ascending=False); T = T.sort_values(["eligible", "est"], ascending=[False, False])
T.to_csv(f"{O}/E45_single_arm_gene_pooled.csv", index=False)
# step 0 against the E40 export, CD8A per cohort
S40 = pd.read_csv(f"{O}/E40_single_arm_genes.csv.gz"); S40 = S40[~S40.set_index(["cohort", "sample"]).index.isin(list(zip(fl.development, fl["sample"])))]
dmax = max(abs(lfit(np.column_stack([np.ones((S40.cohort == c).sum()), std(S40[S40.cohort == c]["CD8A"].to_numpy(float))]), S40[S40.cohort == c].response.to_numpy(float), [1], True)[0][0] -
               lfit(np.column_stack([np.ones(len(Y[c])), std(GZ[c]["CD8A"].to_numpy(float))]), Y[c], [1], True)[0][0]) for c in COH)
print(f"step 0: CD8A log-odds per cohort, this export against E40, max |diff| {dmax:.1e}"); assert dmax < 1e-8
E = T[T.eligible]; print(f"eligible candidates: {len(E)} of {len(T)}; measured in 20 cohorts: {(T.k == 20).sum()}")
print(E.head(15)[["gene", "k", "k_pos", "est", "lo", "hi", "p", "panel"]].round(3).to_string(index=False))
G1 = E.iloc[0].gene; G2 = E[E.panel].iloc[0].gene
print("panel genes:"); print(E[E.panel][["gene", "k", "est", "lo", "hi", "p", "rank"]].round(3).to_string(index=False))
print("G1:", G1, "| G2:", G2, "| rank of CD8A:", int(T.set_index("gene").loc["CD8A", "rank"]), "of", len(E))
pd.DataFrame([dict(reference="G1", gene=G1, rule="highest pooled log-odds among eligible candidates"),
              dict(reference="G2", gene=G2, rule="highest pooled log-odds among the 27 panel genes")]).to_csv(f"{O}/E45_reference_genes.tsv", sep="\t", index=False)
