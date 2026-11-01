"""Writes results/E30_dedup_summary.md from the E30 Q1/Q2 result files (tables are generated, not transcribed)."""
from pathlib import Path
import pandas as pd

A = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
DUP = pd.read_csv(A / "E30_dedup_duplicates.csv")
GATE = pd.read_csv(A / "E30_dedup_gate.csv")
C = pd.read_csv(A / "E30_dedup_compare.csv")
Q2 = pd.read_csv(A / "E30_estimators_gastric.csv")


def md(df):
    df = df.fillna("")
    h = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns) + "\n"
    return h + "".join("| " + " | ".join(str(v).replace("|", "/") for v in r) + " |\n" for r in df.values)


out = []
w = out.append
w("# E30 Q1 (duplicates) and Q2 (estimators, gastric) - summary\n")
w("Protocol `analysis/E30_prereg.md`, sha256 8e76e799...cdd67 (verified with `sha256sum -c` before the run and again inside "
  "`E30_dedup.R`; not edited). Run 2026-10-21 on c0107 (R 4.5.2, metafor 5.0.1). Published scripts and results were not "
  "modified; every re-run is a copy (`analysis/E30_dedup_*.R`) whose only analytic change is the GSE26253 sample filter.\n")

w("## 1. Duplicates selected by the rule\n")
w("Rule implemented in `analysis/E30_dedup.R`: (i) mutual best match, Pearson on all 11,520 shared genes after E20 `dup_pair()` "
  "preprocessing (GSE26253 loaded as in E15, 432 arrays; GSE62254 = benchmark matrix B1, 283 arrays); (ii) same sex and age within "
  "1 year (GSE26253 age/sex/OS from Oh et al. 2018 Suppl. Data 2, sheet SMC, linked by the unique key RFS/recurrence/stage, 402 of 432 "
  "arrays; GSE62254 from sheet ACRG); (iii) time compatibility with 0.5-month tolerance as written in the protocol.\n")
n_recip = len(DUP); n_ii = int(DUP.rule_ii.sum())
w(f"Mutual best matches: {n_recip}; with sex and age agreement: {n_ii}; meeting all three conditions: {int(DUP.duplicate.sum())}. "
  f"**The rule selects exactly the 8 expected arrays; the 3-array sensitivity set (matching event times) is also exactly as expected.** "
  "Four mutual best matches have no unique SMC link (two SMC rows share their key); for none of them does either candidate SMC row agree "
  "with the GSE62254 partner on sex and age, so the missing link does not hide a duplicate.\n")
t = DUP[DUP.rule_ii].copy()
t["set"] = t.apply(lambda r: "primary + sens3" if r.set_sens3 else ("primary" if r.set_primary else "not duplicate (time)"), axis=1)
t["age26"] = t["age26"].map(lambda x: f"{x:.0f}")
for col in ["max_r", "row_z", "rfs26", "dfs62", "os26", "os62"]:
    t[col] = t[col].map(lambda x: f"{x:.2f}")
w(md(t[["sample", "best_match", "max_r", "row_z", "sex26", "age26", "age62", "rec26", "rfs26", "dfs62", "os26", "os62", "set"]]))
w("\n(rfs26/os26: GSE26253 recurrence-free time and overall survival, months; dfs62/os62: GSE62254 disease-free time and overall survival. "
  "GSM644782 and GSM644878 fail (iii) because GSE26253 records recurrence at 18.4 / 25.9 months against GSE62254 DFS 0.5 / 58.8; "
  "GSM644832 is censored at 94.2 months in GSE26253 against GSE62254 DFS 91.9.)\n")
w("Effect on cohort counts: GSE26253 432 patients (177 recurrences) -> **424 (174)**; sensitivity 429 (174). Locked first three "
  "679 (299) -> 671 (296); locked seven 1,055 (501) -> 1,047 (498); all ten 2,152 (999) -> 2,144 (996).\n")

w("## 2. Analyses that use GSE26253 and re-runs\n")
w("Found by grepping every script in `analysis/` and every script that reads E15/E20/E29 outputs: **E15_gastric_validation.R** "
  "(locked first three), **E20_extended_replication.R** (locked seven, all ten, sensitivity sets, RFS secondary set, structure-matched "
  "null, ΔC; GSE26253 also enters its duplicate screen as a comparator), **E29_dissociation.R** (survival contrasts, ten gastric) and "
  "**E29b_posthoc_adjusted.R**. Not affected (benchmark cohorts or checkpoint data only): E28_gastric.R, E19_robustness.R, E4, E1, E2, "
  "E8, E8b, E5, E9, E23, E11, E12 and the trial scripts. Figure and workbook scripts (`figures/cns/make_fig3/5/8.py`, "
  "`supp/make_supp_gastric.py`, `manuscript/make_v30_workbook.py`) only read these results files.\n")
w("Copies: `E30_dedup_E15.R`, `E30_dedup_E20.R`, `E30_dedup_E29.R`, `E30_dedup_E29b.R` (`E30_MODE=off|primary|sens3`). The filter is "
  "one line in E15 and E20 (`cl <- cl[!cl$sample %in% DROP, ]`); E20 reads the locked-first-three estimates and ΔC from the E15 copy "
  "of the same mode, E29/E29b read the scores of the same mode; outputs go to `results/E30_dedup_runs/<mode>/` (RDS intermediates to "
  "the session scratchpad). Seeds unchanged (E15/E20 2026, E29 2029). Random-gene nulls were re-run with the same draws: the gene pools "
  "are unchanged (11,520 for E15, 8,343 for E20) and the gene sets are drawn first, so every random axis is identical; only GSE26253's "
  "scores (and the pools that contain it) differ. The bootstrap streams after GSE26253 shift because it has fewer patients (E15 ΔC for "
  "GSE26253; E29 bootstrap SEs of the four added cohorts, changes ≤ 0.002). Runtime: E15 ~3 min, E20 ~6 min, E29 ~3 min, E29b ~1 min "
  "(three modes in parallel); no null needed approval.\n")

w("## 3. Step 0 gate (filter off)\n")
g = GATE.groupby("script").agg(files=("file", "size"), numeric_cells=("n_numeric", "sum"),
                               max_abs_diff=("max_abs_diff", "max"), text_mismatches=("text_mismatch", "sum"),
                               all_passed=("passed", "all")).reset_index()
w(md(g))
w("\nEvery published results file (including null draws and bootstrap output) is reproduced exactly (max |difference| 0; "
  "tolerance 1e-6). In addition, E20's internal check (E15 cohorts re-scored reproduce E15_per_cohort.csv, 5e-16) and E29's own "
  "step 0 (4.9e-15) passed in all three modes. Details: `results/E30_dedup_gate.csv`. No branch was stopped.\n")

w("## 4. Changed reported numbers (main text, SI, tables, figure text)\n")
T = C[C.section == "text"].copy()
ch = T[T.reported_precision_change == "yes"][["quantity", "reported_old", "reported_dedup", "reported_sens3", "where_reported"]]
ch.columns = ["quantity", "published", "duplicate-free (8)", "sensitivity (3)", "where reported"]
w(md(ch))
un = T[T.reported_precision_change == "no"]
w("\nUnchanged at reported precision: " + "; ".join(f"{q} ({o})" for q, o in zip(un.quantity, un.reported_old)) + ".\n")

w("## 5. Pooled estimates and conclusions\n")
W = C[C.section == "workbook"]
pts = W[W.column.isin(["HR", "est"]) & W.file.str.contains("pooled") & ((W.abs_change > 0.01) | (W.abs_change_sens3 > 0.01))]
w("Pooled point estimates that move by more than 0.01 on their reported (HR) scale:\n")
w(md(pts[["quantity", "published_value", "dedup_value", "sens3_value", "abs_change", "where_reported"]].round(4)))
w("\nThe locked-first-three NetImmune pool (unadjusted) moves by 0.0097 (just under 0.01), but its conservative interval and P cross the "
  "threshold. P values crossing 0.05 (duplicate-free or sensitivity):\n")
pc = W[W.column.isin(["p", "P"])].copy()
pc = pc[((pc.published_value < 0.05) != (pc.dedup_value < 0.05)) | ((pc.published_value < 0.05) != (pc.sens3_value < 0.05))]
w(md(pc[["quantity", "published_value", "dedup_value", "sens3_value", "where_reported"]].round(4)))
w("\nWorkbook: changed cells by results file / sheet (all changed cells are listed in `E30_dedup_compare.csv`, section 'workbook'):\n")
wb = W.groupby("file").agg(cells_changed=("quantity", "size"), changed_at_sheet_precision=("reported_precision_change", lambda x: (x == "yes").sum()),
                           sheet=("where_reported", lambda x: sorted(set(s.split(" (")[0] for s in x))[0])).reset_index()
w(md(wb))
w("\nUnchanged sheets among the gastric ones: S7c (coverage), S7g (extension cohorts), S7h (extension duplicate exclusions - the "
  "16 + 4 exclusions are identical in all modes), S7i (added-cohort estimates), S8, S2a/S2b (benchmark).\n")

w("## 6. Protocol expectation\n")
w("Expectation: *No reported pooled estimate changes by more than 0.01 on its reported scale; no conclusion changes.* "
  "**Not held.**\n")
w("- *Estimates:* every pooled estimate quoted in the main-text paragraphs moves by less than 0.01 in its point estimate "
  "(largest: locked first three NetImmune 0.7344 -> 0.7247, Δ 0.0097). But pooled HRs printed on Fig. 3d and in the Supplementary "
  "Tables move by more than 0.01: locked first three CAF_ECM 1.258 -> 1.269 (Fig. 3d label 1.26 (0.74–2.14), P = 0.21 -> 1.27 "
  "(0.75–2.15), P = 0.19; S7e, S7j) and the secondary RFS set (GSE26253 + three added cohorts with RFS, adjusted model; both extension-duplicate "
  "versions: NetImmune 0.764 -> 0.752 and 0.765 -> 0.753, CAF_ECM 1.228 -> 1.244 and 1.217 -> 1.232, Barrier 1.176 -> 1.187; S7j only). Interval bounds move by up to 0.044 "
  "(locked-three CAF_ECM adjusted upper bound), and the locked-seven I² falls from 7% to 1% (MS [120]).\n")
w("- *Conclusion:* the locked-first-three NetImmune pool becomes significant under the conservative rule: 0.73 (0.53–1.02; P = 0.057) "
  "-> **0.72 (0.53–1.00; P = 0.049)**; the upper bound is 0.9964 (sensitivity, 3 arrays: 0.72 (0.53–0.99; P = 0.048)). MS [120] states "
  "that this estimate 'was not significant under the conservative rule', and Fig. 3d would draw its diamond filled instead of open. "
  "The P = 0.057 quoted as the reason for the extension (MS [14], [154], [272], [281], [285]; SI [14], [163], [288]) is the value on "
  "which that decision was taken, so it remains historically accurate, but under the decision rule the primary locked-three value "
  "would now be P = 0.049. How to word this is an authors' decision. Under REML-HK the published value was already 0.047 (now 0.041).\n")
w("- Unchanged conclusions: locked seven NetImmune 0.70 (0.62–0.80 -> 0.62–0.79), P 4.3 × 10⁻⁴ -> 2.9 × 10⁻⁴, beyond every random axis; "
  "all ten 0.72; CAF_ECM adverse in six of seven; ten of ten HRs below 1, nine significant; E29 survival contrast −0.08 (−0.22 to 0.06; "
  "P = 0.22) and every E29 statement unchanged at reported precision; E20 expectations X1–X3 held as before.\n")
w("- Other P values crossing 0.05 occur only in workbook rows (GSE26253 per-cohort TLS_B and CD8A unadjusted, CAF_ECM and ESTIMATE "
  "adjusted; a no-exclusion RFS-set IFN row; GSE26253 ΔC vs CD8A in the sensitivity set only); none is cited in the text.\n")
w("- Text statements affected by the decision rule (not numbers): per the protocol, the Methods (MS [41]-[42], Duplicate screening) "
  "add that the 99.9th-percentile rule cannot detect cross-platform duplicates and describe the E30 rule; the Supplementary Methods "
  "statement in SI [13] ('Expression matching found no duplicate ... between GSE26253 and GSE62254') is withdrawn. MS [42] 'none of the "
  "flagged pairs matched on age, sex and survival' and MS [120] 'none of which matched a patient of another cohort' refer to the "
  "extension screen and remain true for those arrays.\n")

w("## 7. Figures with plotted values moving by more than 0.01 (not re-rendered)\n")
F = C[C.section == "figure"]
fm = F[F.abs_change > 0.01][["quantity", "published_value", "dedup_value", "abs_change"]].round(4)
w(md(fm))
w("\n- **Fig. 3d**: GSE26253 rows (NetImmune 0.75 (0.64–0.87) -> 0.73 (0.63–0.86); CAF_ECM 1.23 (1.04–1.47) -> 1.26 (1.05–1.50)), "
  "the 'Pooled, first three' diamonds (NetImmune upper bound 1.02 -> 1.00 and fill open -> filled because P < 0.05; CAF_ECM 1.26 -> 1.27) "
  "and the printed labels ('432/177' -> '424/174'; pooled labels for first three, locked 7 and all 10).\n"
  "- **Fig. 3j**: GSE26253 tiles CAF_ECM, IFN, TLS_B, Immune, NetImmune (> 0.01), APM and T_CELL (< 0.01; printed value changes for APM); "
  "pooled-row printed values and P change at reported precision but move < 0.01.\n"
  "- **Fig. 3e / 3f, SI Fig. S10b**: plotted moves < 0.01 (locked-7 |log HR| 0.347 -> 0.355; GSE26253 ΔC vs CD8A 0.044 -> 0.050 "
  "(Δ 0.006); locked-7 ΔC vs CD8A +0.062 -> +0.063), but printed labels change.\n"
  "- **Fig. 5c, Fig. 8a** (all-ten pooled HRs): largest move 0.004 on the −log HR scale; no re-render needed; Fig. 8a column label "
  "'10 cohorts; 2,152 patients' and Fig. 1b '2,152 patients' -> 2,144 (text only). **Fig. 8b** row text '0.70 (0.62–0.80)' -> "
  "'0.70 (0.62–0.79)'.\n"
  "- SI Fig. S10a (added cohorts) is unchanged.\n")

w("## 8. Q2 (gastric): estimators, with and without the duplicates\n")
q = Q2[Q2.reported_in_main_text.notna() & (Q2.duplicates != "without_duplicates_3_sens")].copy()
q["value"] = q.apply(lambda r: f"{r.HR:.2f} ({r.HR_lo:.2f}–{r.HR_hi:.2f}); P = {r.p:.2g}", axis=1)
pv = q.pivot_table(index=["set", "score", "model", "duplicates"], columns="estimator", values="value", aggfunc="first").reset_index()
pv = pv[["set", "score", "model", "duplicates", "DL", "REML", "DL-HK", "REML-HK", "manuscript rule"]]
w(md(pv))
w("\nAll k ≤ 7 gastric pools reported in the main text (benchmark three, locked first three, extension four, locked seven); "
  "full table with τ², I², the HK t critical value (4.30 for k = 3) and the 3-array sensitivity in `results/E30_estimators_gastric.csv` "
  "(1,230 rows: 14 scores, unadjusted/adjusted). The manuscript-rule rows reproduce the published E20_pooled.csv values (max |Δ| 5e-15). "
  "Benchmark and extension pools contain no GSE26253 and are identical with and without the duplicates.\n")

w("## Files\n")
for f in ["analysis/E30_dedup.R", "analysis/E30_dedup_E15.R", "analysis/E30_dedup_E20.R", "analysis/E30_dedup_E29.R",
          "analysis/E30_dedup_E29b.R", "analysis/E30_dedup_gate.py", "analysis/E30_dedup_compare.py",
          "analysis/E30_dedup_estimators.R", "analysis/E30_dedup_summary.py",
          "results/E30_dedup.log", "results/E30_dedup_duplicates.csv", "results/E30_dedup_gate.csv",
          "results/E30_dedup_compare.csv", "results/E30_estimators_gastric.csv", "results/E30_dedup_summary.md",
          "results/E30_dedup_runs/{off,primary,sens3}/ (full re-run outputs and logs)"]:
    w(f"- `{f}`")
(A / "E30_dedup_summary.md").write_text("\n".join(out))
print("\n".join(out)[:200])
