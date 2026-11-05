"""E41 adjudication (protocol E41_prereg.md): reconcile the two blind extraction passes into results/E41_adjudicated.csv.

Rows come from pass 2 for the comparisons both passes or either pass judged eligible with numbers; every value was
checked against the verbatim quotes of both passes, and the decisions below (RESOLVE) record each disagreement and how
it was settled. Writes results/E41_adjudicated.csv, E41_agreement.csv and E41_decisions.csv."""
import re
from pathlib import Path
import pandas as pd

H = Path(__file__).resolve().parent
A = pd.read_csv(H / "E41_pass1/extraction.csv", dtype=str).fillna("")
B = pd.read_csv(H / "E41_pass2/extraction.csv", dtype=str).fillna("")

# eligible, estimable comparisons (independent of the six patient-level comparisons): pass-2 row name -> tumor group
KEEP = {"ALICE": "breast", "CheckMate 459": "liver", "CheckMate 649": "gastric", "CLEAR (KEYNOTE-581)": "renal",
        "COMBI-i (part 3)": "melanoma", "ICON": "breast", "IMbassador250": "prostate", "IMbrave150": "liver",
        "IMpower133": "small-cell lung", "IMspire150": "melanoma", "IMvigor010": "urothelial", "KEYNOTE-061": "gastric",
        "KEYNOTE-119": "breast", "KEYNOTE-426": "renal", "OAK": "non-small-cell lung", "ORIENT-11": "non-small-cell lung",
        "POPLAR": "non-small-cell lung", "RATIONALE-309": "nasopharyngeal"}
RESOLVE = [
    ("IMvigor010", "signature, PD-L1 and TMB HRs",
     "pass 1: OS adjusted for baseline prognostic factors (ED Fig. 6a, values mapped by text-layer order); pass 2: "
     "univariable OS (ED Fig. 2b, printed values)", "pass 2 (univariable; most trials report unadjusted or stratified HRs)"),
    ("CheckMate 649", "PD-L1 HRs", "pass 1: 0.79 (0.60-1.05) / 0.83 (0.58-1.18), WES-evaluable; pass 2: 0.77 (0.59-1.01) "
     "/ 0.91 (0.65-1.27), RNA-seq-evaluable (ED Fig. 6, quoted)", "pass 2 (same population as the signature)"),
    ("CheckMate 649", "primary signature", "10-gene inflammatory (both passes); the 12-gene chemokine signature has more "
     "genes", "10-gene inflammatory: an IFN-gamma signature (Ayers), the class named first in the protocol's rule; the "
     "chemokine signature is recorded as secondary"),
    ("CheckMate 459", "low-group HR and comparators", "pass 1: abstract only (high group); pass 2: Fig. 2A, 3A, 1A printed "
     "values, quoted", "pass 2"),
    ("CLEAR (KEYNOTE-581)", "signature HRs", "pass 1: numbers in a figure it could not retrieve; pass 2: Fig. 3B printed "
     "values, quoted", "pass 2"),
    ("IMbassador250", "low-group HR and comparators", "pass 1: text only (high group); pass 2: Fig. 5a printed values, "
     "quoted", "pass 2"),
    ("IMbrave150", "eligibility and HRs", "pass 1: Zhu 2022 not accessible; an immune-signature-score paper pooling a "
     "non-randomized arm recorded and excluded; pass 2: Zhu 2022 Fig. 2e printed labels (Teff, median, adjusted HRs)",
     "pass 2 (randomized IMbrave150 patients only; HRs adjusted as published)"),
    ("IMpower133", "eligibility and HRs", "pass 1: not accessible; pass 2: Fig. 6D read from point and whisker positions "
     "(no printed values, about +/-0.02)", "pass 2, flagged from_figure = estimated; dropped in the sensitivity "
     "analysis without figure-read numbers"),
    ("ORIENT-11", "eligibility and split", "pass 1: not accessible; pass 2: Fig. 2 printed labels for high, medium and "
     "low immune-infiltration clusters", "pass 2; signature class eligible as an immune-infiltration score (doubt "
     "recorded); high vs low with the medium cluster set aside, as for tertile splits (frac_high 40/171, frac_low "
     "53/171 from numbers at risk)"),
    ("POPLAR", "eligibility and HRs", "pass 1: not accessible; pass 2: Fig. 4A printed values, quoted", "pass 2"),
    ("OAK", "signature class", "AP/T score (XCR1, CXCL9, CXCL10, CXCL11)", "eligible: three of four genes are "
     "IFN-gamma-induced chemokines (doubt recorded)"),
    ("RATIONALE-309", "signature class", "GEP-defined hot vs cold tumor-microenvironment class", "eligible as an "
     "immune-infiltration class (doubt recorded)"),
    ("AtezoTRIBE", "eligibility", "both passes: abstract only, HRs 0.39 and 0.83 without 95% CIs (interaction P 0.066)",
     "not eligible: no 95% CI (protocol)"),
    ("IMmotion151; IMpower150", "eligibility", "both passes: only the signature-high HR found in accessible text",
     "not eligible: low-group HR not reported in the accessible text"),
    ("CheckMate 214; CASPIAN; JAVELIN Bladder 100; CheckMate 274; CheckMate 743; KEYNOTE-045", "eligibility",
     "within-arm associations and/or interaction P only (pass 1 unclear for three not accessible)",
     "not eligible (protocol: no derived interactions)"),
    ("RATIONALE-303; spartalizumab NPC phase 2; IMvigor pooled urothelial trials; GeparNuevo", "eligibility",
     "pass 1 unclear (not accessible) or no; pass 2 no", "not eligible: no signature-stratified treatment HRs found"),
]
rows = []
for name, grp in KEEP.items():
    r = B[(B.trial == name) & (B.eligible == "yes")]
    assert len(r) == 1, name
    rows.append(r.iloc[0].to_dict() | {"tumor_group": grp})
D = pd.DataFrame(rows)
D["frac_low"] = ""
D.loc[D.trial == "ORIENT-11", ["frac_high", "frac_low"]] = [f"{40 / 171:.3f}", f"{53 / 171:.3f}"]
D.loc[D.trial == "CheckMate 459", ["frac_high", "frac_low"]] = [f"{156 / 469:.3f}", f"{157 / 469:.3f}"]
D.loc[D.trial == "CheckMate 649", ["frac_high", "frac_low"]] = [f"{134 / 402:.3f}", f"{134 / 402:.3f}"]
D.loc[D.trial == "IMpower133", "from_figure"] = "estimated"
D["eligible"] = "yes"; D["in_ipd_six"] = "no"; D["int_scale"] = ""
D.loc[D.trial == "CheckMate 649", ["int_test_reported", "int_p"]] = ["yes", "0.045"]
# sensitivity 4: the published subgroup result of one of the six (IMmotion150, atezolizumab + bevacizumab vs sunitinib)
m = B[B.trial.str.startswith("IMmotion150")].iloc[0].to_dict()
m |= {"eligible": "yes", "in_ipd_six": "yes", "tumor_group": "renal", "frac_low": "", "int_scale": "",
      "checkpoint_arm": "atezolizumab + bevacizumab", "from_figure": "yes"}
D = pd.concat([D, pd.DataFrame([m])], ignore_index=True)
for c in ["hr_high", "lo_high", "hi_high", "hr_low", "lo_low", "hi_low"]:
    assert D[c].str.len().gt(0).all(), c
out = H / "results"
D.to_csv(out / "E41_adjudicated.csv", index=False)
pd.DataFrame(RESOLVE, columns=["comparison", "item", "passes", "decision"]).to_csv(out / "E41_decisions.csv", index=False)

# agreement on eligibility over every comparison that either pass rated eligible or unclear
def k(t):
    return re.sub(r"[^a-z0-9]", "", re.sub(r"\(.*?\)", "", t.lower()).replace("keynote-581", ""))
def verdict(df):
    v = {}
    for _, x in df.iterrows():
        e = x.eligible.lower(); e = "yes" if e.startswith("yes") else "unclear" if e.startswith("unclear") else "no"
        kk = k(x.trial); v[kk] = "yes" if "yes" in (v.get(kk), e) else (e if v.get(kk) in (None, "no") else v[kk])
    return v
va, vb = verdict(A), verdict(B)
cand = sorted(x for x in set(va) | set(vb) if "yes" in (va.get(x), vb.get(x)) or "unclear" in (va.get(x), vb.get(x)))
G = pd.DataFrame({"comparison": cand, "pass1": [va.get(x, "not assessed") for x in cand],
                  "pass2": [vb.get(x, "not assessed") for x in cand]})
G["agree"] = G.pass1 == G.pass2
G.to_csv(out / "E41_agreement.csv", index=False)
both = [x for x in cand if va.get(x) == "yes" and vb.get(x) == "yes"]
print(f"candidates {len(G)}; eligibility agreed {G.agree.sum()}; pass1 unclear {(G.pass1 == 'unclear').sum()}; "
      f"both eligible {len(both)}")
print(D[["trial", "endpoint", "frac_high", "frac_low", "hr_high", "lo_high", "hi_high", "hr_low", "lo_low", "hi_low",
         "from_figure"]].to_string())
