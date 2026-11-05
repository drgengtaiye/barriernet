"""E38 - reconcile the two independent extraction passes, then count and cross-tabulate (protocol E38_prereg.md).

Inputs: E38_pass1.csv, E38_pass2.csv (independent AI-assistant passes, blind to each other and to the results).
Disagreements are resolved in RESOLVE below from the quoted evidence (and the paper where needed); each resolution
carries its reason. Writes results/E38_adjudicated.csv, E38_agreement.csv, E38_crosstab.csv, E38_expectations.csv.
"""
import hashlib, sys
import numpy as np, pandas as pd
A = "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"; O = f"{A}/results"
_h = hashlib.sha256(open(f"{A}/E38_prereg.md", "rb").read()).hexdigest()
assert _h == open(f"{A}/E38_prereg.sha256").read().split()[0], "E38 protocol hash mismatch"
print("E38 protocol sha256 verified:", _h[:12])
ITEMS = ["D1", "D2", "D3", "D4"]
P1 = pd.read_csv(f"{A}/E38_pass1.csv", dtype=str, keep_default_na=False).set_index("signature")
P2 = pd.read_csv(f"{A}/E38_pass2.csv", dtype=str, keep_default_na=False).set_index("signature")
SRC = pd.read_csv(f"{A}/E38_sources.csv", dtype=str).set_index("signature")
assert set(P1.index) == set(P2.index) == set(SRC.index) and len(SRC) == 52
norm_ = lambda x: "NA" if str(x).strip().upper() in ("", "NAN", "N/A", "NA", "NONE") else str(x).strip().upper()
for P in (P1, P2):
    for c in ITEMS:
        P[c] = P[c].map(norm_)


def kappa(a, b):
    cats = sorted(set(a) | set(b)); n = len(a)
    po = np.mean(np.array(a) == np.array(b))
    pe = sum((np.sum(np.array(a) == c) / n) * (np.sum(np.array(b) == c) / n) for c in cats)
    return (po - pe) / (1 - pe) if pe < 1 else 1.0


AG = []
for c in ITEMS:
    a, b = P1.loc[SRC.index, c].tolist(), P2.loc[SRC.index, c].tolist()
    AG.append(dict(item=c, agree=int(np.sum(np.array(a) == np.array(b))), n=52, pct=100 * np.mean(np.array(a) == np.array(b)),
                   kappa=kappa(a, b)))
AG = pd.DataFrame(AG); AG.to_csv(f"{O}/E38_agreement.csv", index=False); print(AG.round(3).to_string(index=False))

# (signature, item) -> (value, reason); filled after reading the disagreements
RESOLVE = {
    ("ADO_Sidders", "D2"): ("YES", "pass 1 read the full text (journal site): CD274 mRNA hazard ratio shown beside the signature; pass 2 had the abstract only"),
    ("APM_Thompson", "D2"): ("NO", "single-gene results are cited from the authors' earlier study, not reported in this paper"),
    ("CRMA_Shukla", "D2"): ("NO", "only the signature's own component genes are reported individually; a decomposition is not a comparison with a single-gene marker"),
    ("CRMA_Shukla", "D4"): ("NO", "adjusted for neoantigen load, which is neither tumor mutational burden nor PD-L1 immunohistochemistry as the rubric specifies"),
    ("ImmunoScore_Hao", "D2"): ("YES", "pass 1 read the full text (journal site): predictive power stated as comparable with expression of the immunotherapy target genes; pass 2 had the abstract only"),
    ("MHC-II_Liu", "D2"): ("NO", "only component HLA class II genes reported individually (decomposition, as for CRMA_Shukla)"),
    ("MHC-I_Liu", "D2"): ("NO", "only component MHC-I genes reported individually (decomposition, as for CRMA_Shukla)"),
    ("NonResponse_Chen", "D1"): ("D", "the full text reports a test in ICB-treated bulk cohorts (no significant association), so checkpoint-outcome data exist"),
    ("NonResponse_Chen", "D2"): ("NO", "follows D1 = D; no single-gene comparison reported"),
    ("NonResponse_Chen", "D3"): ("NO", "follows D1 = D; no random gene sets"),
    ("NonResponse_Chen", "D4"): ("NO", "follows D1 = D; no PD-L1 or TMB comparison"),
    ("Response_Chen", "D1"): ("D", "as NonResponse_Chen (same paper and test)"),
    ("Response_Chen", "D2"): ("NO", "follows D1 = D"),
    ("Response_Chen", "D3"): ("NO", "follows D1 = D"),
    ("Response_Chen", "D4"): ("NO", "follows D1 = D"),
}

rows, unresolved = [], []
for s in SRC.index:
    r = dict(signature=s, PMID=SRC.loc[s, "PMID"], first_author_year=SRC.loc[s, "first_author_year"])
    for c in ITEMS:
        v1, v2 = P1.loc[s, c], P2.loc[s, c]
        if v1 == v2:
            r[c], r[f"{c}_how"] = v1, "both passes"
        elif (s, c) in RESOLVE:
            r[c], r[f"{c}_how"] = RESOLVE[(s, c)][0], f"adjudicated ({v1} / {v2}): {RESOLVE[(s, c)][1]}"
        else:
            r[c], r[f"{c}_how"] = None, f"UNRESOLVED ({v1} / {v2})"; unresolved.append((s, c, v1, v2))
    r["quote_pass1"], r["quote_pass2"] = P1.loc[s, "quote_D1"], P2.loc[s, "quote_D1"]
    r["note_pass1"], r["note_pass2"] = P1.loc[s, "note"], P2.loc[s, "note"]
    r["source_pass1"], r["source_pass2"] = P1.loc[s, "source"], P2.loc[s, "source"]
    r["abstract_only_both"] = (P1.loc[s, "source"] == "abstract_only_no_fulltext") and (P2.loc[s, "source"] == "abstract_only_no_fulltext")
    rows.append(r)
if unresolved:
    print(f"\n{len(unresolved)} unresolved disagreements:")
    for s, c, v1, v2 in unresolved:
        print(f"  {s} {c}: pass1 {v1} | pass2 {v2}")
        print(f"     1: {P1.loc[s, 'quote_D1']} || {P1.loc[s, 'note']}")
        print(f"     2: {P2.loc[s, 'quote_D1']} || {P2.loc[s, 'note']}")
    sys.exit(0)
ADJ = pd.DataFrame(rows); ADJ.to_csv(f"{O}/E38_adjudicated.csv", index=False)

# ------------------------------------------------------------------ counts and cross-tabulation
print("\nD1:", ADJ.D1.value_counts().reindex(list("ABCDE"), fill_value=0).to_dict())
for c in ["D2", "D3", "D4"]:
    print(c, ADJ[c].value_counts().to_dict())
E32 = pd.read_csv(f"{O}/E32_primary.csv").set_index("signature")
ADJ = ADJ.set_index("signature")
ADJ["favors_checkpoint_p05"] = (E32.loc[ADJ.index, "est"] > 0) & (E32.loc[ADJ.index, "p"] < 0.05)
ADJ["exceeds_random"] = E32.loc[ADJ.index, "emp_p_declared"] < 0.05
ADJ["above_CD8A_p05"] = (E32.loc[ADJ.index, "diff_vs_CD8A"] > 0) & (E32.loc[ADJ.index, "cd8a_p"] < 0.05)
ADJ["group"] = ADJ.D1.map({"A": "randomized (A/B)", "B": "randomized (A/B)", "C": "treated vs untreated (C)",
                           "D": "checkpoint-treated only (D)", "E": "no checkpoint data (E)"})
CT = ADJ.groupby("group")[["favors_checkpoint_p05", "exceeds_random", "above_CD8A_p05"]].agg(["sum", "count"])
CT.to_csv(f"{O}/E38_crosstab.csv"); print(CT.to_string())
nA = int((ADJ.D1 == "A").sum()); nDE = int(ADJ.D1.isin(["D", "E"]).sum())
nD2 = int((ADJ.D2 == "YES").sum()); nD3 = int((ADJ.D3 == "YES").sum())
EXP = pd.DataFrame([
    dict(q="L1", expectation="at most 3 of 52 with D1 = A", result=nA, held="yes" if nA <= 3 else "no"),
    dict(q="L2", expectation="at least 38 of 52 with D1 = D or E", result=nDE, held="yes" if nDE >= 38 else "no"),
    dict(q="L3", expectation="at most 12 of 52 with D2 = yes", result=nD2, held="yes" if nD2 <= 12 else "no"),
    dict(q="L4", expectation="at most 3 of 52 with D3 = yes", result=nD3, held="yes" if nD3 <= 3 else "no")])
EXP.to_csv(f"{O}/E38_expectations.csv", index=False); print(EXP.to_string(index=False))
ADJ.reset_index().to_csv(f"{O}/E38_adjudicated.csv", index=False)
