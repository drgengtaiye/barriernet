"""Assemble lean_supp_src.md (Supplementary Information for the lean manuscript) from v17_src.md:
Supplementary Methods (the full v17 Methods, without version-history remarks), the five v17 main
tables that move to the supplement, a guide to the workbook, and legends S1-S20 in the lean order.
Citations stay as [@key]; build_lean.py numbers them separately for the supplement.

One-time generator: the v17 legend input below was a scratch export of the v17 online document and no
longer exists (session scratch space was cleared on 2026-10-03). lean_supp_src.md is now the maintained
source and is edited directly; do not re-run this script."""
import re, pathlib
HERE = pathlib.Path(__file__).resolve().parent
v17 = (HERE / "v17_src.md").read_text()
supp_tab = pathlib.Path("/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-com-gastric-2/"
                        "727be15e-8d67-407e-a66a-a54d26cc1ebb/scratchpad/doc17/supp.md").read_text()


def between(s, a, b):
    i = s.index(a); j = s.index(b, i)
    return s[i:j].strip()


# ---- Supplementary Methods: v17 Methods minus the version-history remarks ----------------------
meth = between(v17, "### Study design and data", "### Software and reproducibility")
for a, b in [
    (" (Version 15 of this manuscript transformed it a second time; the corrected results are given here, and Supplementary Table S2 shows both.)", ""),
    (" (Version 15 reported DerSimonian–Laird intervals for the first 11 cohorts; the more conservative intervals are given here.)", ""),
    ("(544 patients, 149 responders; Supplementary Table S5)", "(544 patients, 149 responders; Supplementary Table S10)"),
    ("(Supplementary Table S9)", "(Supplementary Table S5)"),
]:
    assert a in meth, a[:60]
    meth = meth.replace(a, b)
meth = meth.replace("\n### ", "\n### ").replace("### Analysis hierarchy", "### Analysis hierarchy (as in the main text)")
meth = re.sub(r"(?m)^### ", "### ", meth)

# ---- fixes from VERIFY_lean ----------------------------------------------------------------------
for a, b in [
    ("TIDE[@jiang2018] requires an online service and could not be computed. This proxy was used in earlier drafts and is shown only for completeness.",
     "TIDE[@jiang2018] requires an online service and could not be computed; this proxy is shown only for completeness."),
]:
    assert a in meth, a[:60]
    meth = meth.replace(a, b)
i = meth.index("The replication cohorts, the I-SPY2 comparisons, the TCGA analysis and the IMC analysis were added after the results above were known.")
j = meth.index("before it was run.", i) + len("before it was run.")
meth = meth[:i] + ("The locked gastric cohorts, replication cohorts, PERFECT, I-SPY2, TCGA, METABRIC and IMC analyses were added later; "
                   "for each, the expected direction of every primary score was written into the analysis code before it was run.") + meth[j:]

# ---- the five v17 main tables that become supplementary tables ----------------------------------
tabs = between(v17, "**Table 2. Head-to-head", "**Table 7.")
tabs = (tabs.replace("**Table 2. Head-to-head", "**Supplementary Table S2. Head-to-head")
            .replace("**Table 3. Incremental", "**Supplementary Table S3. Incremental")
            .replace("**Table 4. Fitted", "**Supplementary Table S4. Fitted")
            .replace("**Table 5. Prediction", "**Supplementary Table S8. Prediction")
            .replace("**Table 6. Combination", "**Supplementary Table S9. Combination"))
old_b = "(b) A 16-gene proxy used in earlier drafts in place of TIDE, which could not be computed offline."
assert old_b in tabs
tabs = tabs.replace(old_b, "(b) A 16-gene CAF/TAM exclusion proxy, shown in place of TIDE, which could not be computed offline.")
# order: S2, S3, S4, then S8, S9 (S5-S7 are workbook sheets)
i8 = tabs.index("**Supplementary Table S8.")
tabs = tabs[:i8] + "Supplementary Tables S5–S7 are sheets of the accompanying workbook (see below).\n\n" + tabs[i8:]

# ---- legends --------------------------------------------------------------------------------------
old_leg = {}
for m in re.finditer(r"\*\*Supplementary Fig\. S(\d+)\.(.*?)(?=\n+!\[|\Z)", supp_tab, flags=re.S):
    old_leg[int(m.group(1))] = m.group(2).strip()
main_leg = {}
for m in re.finditer(r"\*\*Figure (\d)\.(.*?)(?=\n\n!\[|\n\n---|\Z)", v17[v17.index("## Figures"):], flags=re.S):
    main_leg[int(m.group(1))] = m.group(2).strip()
src = {1: ("main", 2), 2: ("old", 1), 3: ("old", 2), 4: ("old", 3), 5: ("old", 4), 6: ("old", 9), 7: ("main", 5),
       8: ("old", 10), 9: ("old", 16), 10: ("old", 13), 11: ("old", 14), 12: ("old", 17), 13: ("main", 3),
       14: ("old", 5), 15: ("old", 6), 16: ("old", 7), 17: ("old", 8), 18: ("old", 15), 19: ("old", 11), 20: ("old", 12)}
legs = []
for new, (kind, k) in src.items():
    t = main_leg[k] if kind == "main" else old_leg[k]
    t = t.replace("(as in version 15; DerSimonian–Laird intervals)", "(DerSimonian–Laird intervals; the text and Table 2 give the wider Hartung–Knapp intervals)")
    if new == 6:
        t = t.rstrip() + " Pooled diamonds are DerSimonian–Laird; the text and Table 1 give the wider Hartung–Knapp intervals."
    t = t.replace(" (deposited scale)", "")
    legs.append(f"![](figures/lean/supp/FigS{new}.png){{width=100%}}\n\n**Supplementary Fig. S{new}.{(' ' + t) if not t.startswith('**') else t}")
legs = "\n\n".join(legs)
legs = legs.replace("**Supplementary Fig. S1. **", "**Supplementary Fig. S1. ")

guide = """### Supplementary Tables in the workbook

`Supplementary_Tables_lean.xlsx` holds the following tables as sheets; a README sheet lists the source script for each.

- **Supplementary Table S1.** Head-to-head differences in Harrell's C by cohort, and C of each score by cohort.
- **Supplementary Table S5.** Locked gastric cohorts: cohort table, overlap check, module coverage, per-cohort and pooled hazard ratios, random-gene null and ΔC.
- **Supplementary Table S6.** Prognostic association of each module in the gastric benchmark cohorts: per-cohort HRs and pooled estimates, DerSimonian–Laird and Hartung–Knapp.
- **Supplementary Table S7.** PRJEB25780: structure-matched random-axis nulls for the AUCs, and the full-transcriptome classifier by penalty.
- **Supplementary Table S10.** Checkpoint cohorts: audit (included and excluded cohorts, platform, module coverage), duplicate checks, per-cohort effects, pooled effects for every set and leave-one-out ranges.
- **Supplementary Table S11.** PERFECT per-score effects and pooled checkpoint estimates with chemotherapy–checkpoint cohorts; METABRIC hazard ratios.
- **Supplementary Table S12.** TCGA per-type and pooled hazard ratios for every score and endpoint, the random-gene null and TCGA-BRCA.
- **Supplementary Table S13.** I-SPY2 interactions (unadjusted, adjusted, stratified), within-arm effects, random-gene null and pooled cross-trial estimates.
- **Supplementary Table S14.** Single-cell imaging metrics, model estimates, positive control, macrophage–aggregate contrasts and METABRIC survival.
- **Supplementary Table S15.** GDSC2 associations for all 11 compounds.
"""

out = f"""# Supplementary Information

**barriernet: marker-level dissociation of prognostic and predictive information in the tumor microenvironment**

Contents: Supplementary Methods; Supplementary Tables S2–S4, S8 and S9 (the other supplementary tables are sheets of the accompanying workbook); Supplementary Figures S1–S20; Supplementary references.

---

## Supplementary Methods

{meth}

---

## Supplementary Tables

{tabs}

{guide}

---

## Supplementary Figures

{legs}

---

## Supplementary references

{{{{REFERENCES}}}}
"""
(HERE / "lean_supp_src.md").write_text(out)
print("written", len(out.split()), "words; legends", len(src))
