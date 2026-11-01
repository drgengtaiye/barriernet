"""E30 Q1: published against duplicate-free values for every quantity of the analyses that use GSE26253.

Inputs (read only): results/E15_*, E20_*, E29_*, E29b_* (published); results/E30_dedup_runs/{primary,sens3}/ (re-runs
with the 8 / 3 GSE26253 duplicates removed; the 'off' re-runs reproduce the published files, see E30_dedup_gate.csv);
the v31 manuscript and SI text (paragraph-indexed .txt) and the v31 Supplementary Tables workbook.

Output: results/E30_dedup_compare.csv
  section 'text'     - every number in the main text, SI, Tables and figure text that comes from these analyses,
                       formatted as reported (old -> new), with paragraph [n] / table [Tn] indices;
  section 'workbook' - every numeric cell of the affected results files that changes (primary or sens3), with the
                       workbook sheet that carries it and whether it changes at the sheet's precision;
  section 'figure'   - plotted quantities (scale as plotted) and their change.
"""
import re, math
from pathlib import Path
import numpy as np, pandas as pd, openpyxl

A = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results")
RUN = {m: A / "E30_dedup_runs" / m for m in ("off", "primary", "sens3")}
REV = Path("/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-com-gastric-2/"
           "727be15e-8d67-407e-a66a-a54d26cc1ebb/scratchpad/rev")
WB = Path("/groups/adv2105_gp/yichen/Yi/com/gastric/2/final_v31/BarrierNet_Supplementary_Tables.xlsx")
PRIM = "pre-registered exclusions (primary)"
SCRIPT = {"E15": "E15_gastric_validation.R -> E30_dedup_E15.R", "E20": "E20_extended_replication.R -> E30_dedup_E20.R",
          "E29": "E29_dissociation.R -> E30_dedup_E29.R", "E29b": "E29b_posthoc_adjusted.R -> E30_dedup_E29b.R"}
SHEETS = {"E15_cohort_table.csv": ["S7a_locked_cohorts"], "E15_overlap_check.csv": ["S7b_locked_overlap"],
          "E15_module_coverage.csv": ["S7c_locked_coverage"], "E15_per_cohort.csv": ["S7d_locked_per_cohort"],
          "E15_pooled.csv": ["S7e_locked_pooled"], "E15_null.csv": ["S7f_locked_null_dC"], "E15_deltaC.csv": ["S7f_locked_null_dC"],
          "E15_cindex.csv": [], "E20_cohort_table.csv": ["S7g_ext_cohorts"], "E20_duplicate_exclusions.csv": ["S7h_ext_duplicates"],
          "E20_per_cohort.csv": ["S7i_ext_per_cohort"], "E20_estimates_used_for_pooling.csv": [], "E20_pooled.csv": ["S7j_ext_pooled"],
          "E20_null.csv": ["S7k_ext_null"], "E20_deltaC.csv": ["S7l_ext_deltaC"], "E20_cindex.csv": [],
          "E20_expectations.csv": ["S7m_ext_expectations"], "E29_per_cohort.csv": ["S7n_direct_per_cohort", "S13f_direct_per_cohort"],
          "E29_pooled.csv": ["S7o_direct_pooled", "S13g_direct_pooled"], "E29_Q3_crossover.csv": ["S7o_direct_pooled", "S13g_direct_pooled"],
          "E29b_posthoc_adjusted_per_cohort.csv": ["S7n_direct_per_cohort", "S13f_direct_per_cohort"],
          "E29b_posthoc_adjusted_pooled.csv": ["S7o_direct_pooled", "S13g_direct_pooled"],
          "E20_duplicate_check.csv": [], "E20_duplicate_flagged.csv": [], "E20_module_coverage.csv": []}


def rd(f, m=None):
    return pd.read_csv((RUN[m] if m else A) / f)


# ------------------------------------------------------------------ formatting as in the manuscript
SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def fp(p):
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return "NA"
    if p < 0.001:
        e = math.floor(math.log10(p)); m = p / 10 ** e
        if round(m, 1) >= 10: m, e = m / 10, e + 1
        return f"{m:.1f} × 10{str(e).translate(SUP)}"
    return f"{p:.{2 - math.floor(math.log10(p)) - 1}f}"      # two significant digits


f2 = lambda x: f"{x:.2f}"
f3 = lambda x: f"{x:+.3f}"
s2 = lambda x: f"{x:+.2f}".replace("+0.00", "0.00")
ci2 = lambda h, l, u: f"{h:.2f} ({l:.2f}–{u:.2f})"
cnt = lambda x: f"{int(round(x)):,}"

# ------------------------------------------------------------------ data access per mode
def D(mode):
    g = lambda f: rd(f, mode) if mode != "pub" else rd(f)
    d = {k: g(f"{k}.csv") for k in ["E15_per_cohort", "E15_pooled", "E15_null", "E15_deltaC", "E15_cohort_table",
                                    "E20_pooled", "E20_null", "E20_deltaC", "E20_cohort_table", "E20_estimates_used_for_pooling",
                                    "E29_pooled", "E29_per_cohort", "E29b_posthoc_adjusted_pooled"]}
    return d


def pooled20(d, st, sc, model="unadjusted", rule=PRIM):
    p = d["E20_pooled"]; r = p[(p.set == st) & (p.score == sc) & (p.model == model) & (p.duplicate_rule == rule)]
    assert len(r) == 1, (st, sc, model, rule, len(r)); return r.iloc[0]


def pc15(d, coh, sc, model="unadjusted"):
    p = d["E15_per_cohort"]; r = p[(p.cohort == coh) & (p.score == sc) & (p.model == model)]
    assert len(r) == 1; return r.iloc[0]


def est20(d, coh, sc, model="unadjusted"):
    e = d["E20_estimates_used_for_pooling"]
    r = e[(e.cohort == coh) & (e.score == sc) & (e.model == model) & (e.duplicate_rule == PRIM)
          & (e.group != "locked (added), RFS") & (e.group != "sensitivity")]
    assert len(r) == 1, (coh, sc, model, len(r)); return r.iloc[0]


def e29(d, st, imm="Immune", meth="a"):
    p = d["E29_pooled"]; r = p[(p.set == st) & (p.immune == imm) & (p.method == meth)]
    assert len(r) == 1; return r.iloc[0]


def counts(d):
    c15 = d["E15_cohort_table"]; c15 = c15[c15.status == "analysed"].set_index("cohort")
    c20 = d["E20_cohort_table"].set_index("cohort"); new = ["GSE13861", "GSE26901", "GSE26899", "GSE29272"]
    l3n, l3e = c15.n.sum(), c15.events.sum()
    exn, exe = c20.loc[new, "n_OS"].sum(), c20.loc[new, "deaths"].sum()
    bmn, bme = 383 + 431 + 283, 156 + 207 + 135          # benchmark (B1/E8; unchanged)
    return dict(g26n=c15.loc["GSE26253", "n"], g26e=c15.loc["GSE26253", "events"], l3n=l3n, l3e=l3e,
                l7n=l3n + exn, l7e=l3e + exe, a10n=l3n + exn + bmn, a10e=l3e + exe + bme)


# ------------------------------------------------------------------ section 'text': numbers as reported
SENS_SETS = ["locked 7 (primary)|no exclusions (sensitivity)", "locked 7, OS only (sensitivity)|" + PRIM,
             "locked 7 + GSE57303 (sensitivity)|" + PRIM]          # the three locked-seven sensitivity versions (0.68-0.72)
TEN = ["TCGA", "GSE84437", "GSE62254", "GSE15459", "GSE34942", "GSE26253", "GSE13861", "GSE26901", "GSE26899", "GSE29272"]
L7 = TEN[3:]


def text_items(d):
    """(quantity, script, raw value, reported string, where). One list per mode, same order."""
    c = counts(d); it = []
    add = lambda q, s, v, rep, w: it.append((q, s, v, rep, w))
    add("GSE26253 patients (n)", "E15", c["g26n"], cnt(c["g26n"]), "MS [37]; [T1] GSE26253; SI [13]; S7a_locked_cohorts; Fig. 1b (cohort icon 432, if drawn)")
    add("GSE26253 recurrences (events)", "E15", c["g26e"], cnt(c["g26e"]), "MS [37]; [T1] GSE26253; SI [13]; S7a_locked_cohorts")
    add("Locked first three: patients", "E15", c["l3n"], cnt(c["l3n"]), "MS [36] (679 patients); [T1] Pooled, first three")
    add("Locked first three: events", "E15", c["l3e"], cnt(c["l3e"]), "[T1] Pooled, first three")
    add("Locked seven: patients", "E15+E20", c["l7n"], cnt(c["l7n"]), "MS [14] abstract (1,055 patients); [T1] Pooled, all seven locked")
    add("Locked seven: events", "E15+E20", c["l7e"], cnt(c["l7e"]), "[T1] Pooled, all seven locked")
    add("All ten gastric: patients", "E15+E20", c["a10n"], cnt(c["a10n"]), "[T1] All ten cohorts; MS [281] Fig. 1b legend (2,152 patients); Fig. 1b; Fig. 8a column label (10 cohorts; 2,152 patients)")
    add("All ten gastric: events", "E15+E20", c["a10e"], cnt(c["a10e"]), "[T1] All ten cohorts")
    for mdl, lab, w in [("unadjusted", "univariable", "[T1] GSE26253; MS [120] (0.65, 0.92 and 0.75)"), ("adjusted", "adjusted (stage)", "[T1] GSE26253")]:
        r = pc15(d, "GSE26253", "NetImmune", mdl)
        add(f"GSE26253 NetImmune HR per SD, {lab}", "E15", r.HR, ci2(r.HR, r.lo, r.hi), w)
        add(f"GSE26253 NetImmune P, {lab}", "E15", r.p, fp(r.p), "[T1] GSE26253")
    n = d["E15_null"].set_index("level")
    add("GSE26253 NetImmune null rank (all genes)", "E15", n.loc["GSE26253", "rank"], f"{int(n.loc['GSE26253','rank'])} (P {'≤ 0.001' if n.loc['GSE26253','rank']==1 else '= '+fp(n.loc['GSE26253','P'])})", "[T1] GSE26253")
    p = d["E15_pooled"]; p = p[p.set == "new cohorts"].set_index(["score", "model"])
    for mdl, w in [("unadjusted", "MS [120] 0.73 (0.53–1.02; P = 0.057); [T1] Pooled, first three; P = 0.057 also in MS [14], [154], [272], [281], [285]; SI [14], [163], [288]"),
                   ("adjusted", "[T1] Pooled, first three")]:
        r = p.loc[("NetImmune", mdl)]
        add(f"Locked first three pooled NetImmune HR, {mdl}", "E15", r.HR, ci2(r.HR, r.lo, r.hi), w)
        add(f"Locked first three pooled NetImmune P, {mdl}", "E15", r.p, fp(r.p), w)
        add(f"Locked first three pooled NetImmune: conservative CI excludes 1, {mdl}", "E15", float(r.hi < 1), "yes" if r.hi < 1 else "no",
            "MS [120] 'was not significant under the conservative rule'" if mdl == "unadjusted" else "[T1]")
    add("Locked first three pooled NetImmune null rank", "E15", n.loc["pooled (inverse variance)", "rank"],
        f"{int(n.loc['pooled (inverse variance)','rank'])} (P ≤ 0.001)" if n.loc['pooled (inverse variance)','rank'] == 1 else str(n.loc['pooled (inverse variance)','rank']), "[T1] Pooled, first three")
    for mdl in ["unadjusted", "adjusted"]:
        r = pooled20(d, "locked 7 (primary)", "NetImmune", mdl)
        w = "MS [15] abstract (0.70, 0.62–0.80); MS [120]; MS [148] (pooled 0.70); [T1] Pooled, all seven locked; Fig. 8b row text" if mdl == "unadjusted" else "[T1] Pooled, all seven locked"
        add(f"Locked seven pooled NetImmune HR, {mdl}", "E20", r.HR, ci2(r.HR, r.lo, r.hi), w)
        add(f"Locked seven pooled NetImmune P, {mdl}", "E20", r.p, fp(r.p), "MS [120]; [T1]" if mdl == "unadjusted" else "[T1]")
    r = pooled20(d, "locked 7 (primary)", "NetImmune")
    add("Locked seven pooled NetImmune I²", "E20", r.I2, f"{r.I2:.0f}%", "MS [120] (I² = 7%)")
    nn = d["E20_null"].set_index("level")
    add("Locked seven pooled NetImmune null rank", "E20", nn.loc["locked 7, pooled", "rank"],
        "1 (P ≤ 0.001)" if nn.loc["locked 7, pooled", "rank"] == 1 else str(int(nn.loc["locked 7, pooled", "rank"])), "MS [120] (P ≤ 0.001); [T1]; MS [148]; Fig. 3e")
    r = pooled20(d, "locked 7 (primary)", "CAF_ECM")
    add("Locked seven pooled CAF_ECM HR", "E20", r.HR, ci2(r.HR, r.lo, r.hi), "MS [120] (pooled 1.34, 1.08–1.65)")
    add("Locked seven pooled CAF_ECM P", "E20", r.p, fp(r.p), "S7m_ext_expectations (P = 0.0151)")
    e = d["E20_estimates_used_for_pooling"]
    caf7 = sum(est20(d, k, "CAF_ECM").beta > 0 for k in L7)
    add("Locked seven: cohorts with CAF_ECM HR > 1", "E20", caf7, f"{caf7} of 7", "MS [120], [123] (six of the seven)")
    net10 = sum((est20(d, k, "NetImmune").beta < 0) for k in TEN)
    pz = lambda r: math.erfc(abs(r.beta / r.se) / math.sqrt(2))
    sig10 = sum((est20(d, k, "NetImmune").beta < 0) and pz(est20(d, k, "NetImmune")) < 0.05 for k in TEN)
    add("Ten cohorts: NetImmune HR < 1 (significant)", "E15+E20", sig10, f"{net10} of 10 ({sig10} significant)", "MS [15] abstract (all ten; significant in nine)")
    sig7 = sum((est20(d, k, "NetImmune").beta < 0) and pz(est20(d, k, "NetImmune")) < 0.05 for k in L7)
    add("Locked seven: NetImmune HR < 1 significant", "E15+E20", sig7, f"{sig7} of 7", "MS [148] (significant in six)")
    for sc, w in [("NetImmune", "MS [120] 0.72 (0.67–0.79); [T1] All ten cohorts; Fig. 5c; Fig. 8a"),
                  ("CAF_ECM", "MS [15] abstract (HR 1.30, 1.17–1.45); MS [120]; MS [123]; Fig. 5c; Fig. 8a"),
                  ("Immune", "MS [122] Immune axis 0.84, 0.78–0.91; Fig. 3j; Fig. 5c; Fig. 8a")]:
        r = pooled20(d, "all 10 gastric", sc)
        add(f"All ten pooled {sc} HR", "E20", r.HR, ci2(r.HR, r.lo, r.hi), w)
    r = pooled20(d, "all 10 gastric", "NetImmune")
    add("All ten pooled NetImmune P", "E20", r.p, fp(r.p), "[T1] All ten cohorts (8.8 × 10⁻⁶)")
    mods = [pooled20(d, "all 10 gastric", m).HR for m in ["APM", "IFN", "TLS_B", "T_CELL"]]
    add("All ten pooled immune modules HR range", "E20", max(mods) - min(mods), f"{min(mods):.2f}–{max(mods):.2f}", "MS [122] (HR per SD 0.86–0.89)")
    sens = []
    for k in SENS_SETS:
        st, rule = k.split("|"); sens.append(pooled20(d, st, "NetImmune", rule=rule).HR)
    add("Locked-seven NetImmune sensitivity sets HR range (no exclusions / OS only / + GSE57303)", "E20", max(sens) - min(sens), f"{min(sens):.2f}–{max(sens):.2f}",
        "MS [120] (estimates were similar (0.68–0.72) without the duplicate exclusions, with OS only, or with GSE57303)")
    dc = d["E20_deltaC"].set_index(["cohort", "comparator"])
    for cmp in ["CD8A", "ESTIMATE"]:
        r = dc.loc[("locked 7, pooled (inverse variance)", cmp)]
        add(f"Locked seven pooled ΔC NetImmune − {cmp}", "E20 (+E15 ΔC)", r.dC,
            f"{r.dC:+.3f} ({r.lo:.3f}–{r.hi:.3f})" if r.lo > 0 else f"{r.dC:+.3f} ({r.lo:.3f} to {r.hi:.3f})", "MS [120]; Fig. 3f; SI Fig. S10b")
    r = e29(d, "Q2 survival, 10 gastric")
    add("E29 survival contrast (Immune vs CAF_ECM), ten gastric", "E29", r.est, f"{s2(r.est)} ({s2(r.lo)} to {s2(r.hi)}); P = {fp(r.p)}", "MS [122]; SI [188]; S7o")
    add("E29 survival contrast P", "E29", r.p, fp(r.p), "MS [122]; SI [188]")
    add("E29 survival within-cohort median Spearman rho (Immune, CAF_ECM)", "E29", r.rho_median, f"{r.rho_median:.2f}", "SI [188] (0.01)")
    r = e29(d, "S5 survival, 10 gastric, adjusted")
    add("E29 survival contrast, adjusted for age and stage", "E29", r.est, f"{s2(r.est)}; P = {fp(r.p)}", "MS [122] (+0.02); SI [188] (+0.02, P = 0.66)")
    r = e29(d, "Q2 survival, 10 gastric", meth="b")
    add("E29 survival contrast, marginal with bootstrap covariance", "E29", r.est, f"{s2(r.est)}; P = {fp(r.p)}", "SI [188] (−0.10, P = 0.18)")
    mm = [e29(d, "Q2 survival, 10 gastric", imm=m) for m in ["APM", "IFN", "TLS_B", "T_CELL"]]
    add("E29 survival contrast, single immune modules (range)", "E29", max(x.est for x in mm) - min(x.est for x in mm),
        f"{s2(max(x.est for x in mm))} to {s2(min(x.est for x in mm))}; significant: " + ",".join(x.immune for x in mm if x.p < 0.05),
        "MS [122] (only TLS_B significantly weaker); SI [188] (−0.10 to −0.12, significant only for TLS_B)")
    # ---- numbers printed on figures (figures/cns/Fig3_numbers.md, Fig8_numbers.md, supp/FigS_locked_numbers.md)
    add("Fig. 3d label: GSE26253 n/events", "E15", c["g26n"], f"{int(c['g26n'])}/{int(c['g26e'])}", "Fig. 3d (432/177)")
    for sc in ["NetImmune", "CAF_ECM"]:
        for st, lab in [("locked (first 3), E15 reference", "Pooled, first three"), ("locked 7 (primary)", "Pooled, locked 7"),
                        ("all 10 gastric", "Pooled, all 10 gastric")]:
            r = pooled20(d, st, sc)
            add(f"Fig. 3d label: {sc} {lab}", "E20", r.HR, f"{ci2(r.HR, r.lo, r.hi)}, P = {fp(r.p)}", f"Fig. 3d ({lab})")
        r = pc15(d, "GSE26253", sc)
        add(f"Fig. 3d label: {sc} GSE26253", "E15", r.HR, ci2(r.HR, r.lo, r.hi), "Fig. 3d (GSE26253 row)")
    add("Fig. 3e label: locked 7 observed |log HR|", "E20", abs(nn.loc["locked 7, pooled", "observed_logHR"]),
        f"{abs(nn.loc['locked 7, pooled', 'observed_logHR']):.3f}", "Fig. 3e (Locked 7, pooled)")
    for cmp in ["CD8A", "ESTIMATE"]:
        r = dc.loc[("locked 7, pooled (inverse variance)", cmp)]
        add(f"Fig. 3f / S10b label: locked 7 ΔC vs {cmp}", "E20", r.dC, f"{r.dC:+.3f} ({r.lo:+.3f} to {r.hi:+.3f}), P {fp(r.p)}", "Fig. 3f; SI Fig. S10b")
    for sc in ["CAF_ECM", "TAM_SUP", "APM", "IFN", "TLS_B", "T_CELL", "Immune", "NetImmune"]:
        r = est20(d, "GSE26253", sc)
        add(f"Fig. 3j tile: GSE26253 {sc} HR", "E15", math.exp(r.beta), f"{math.exp(r.beta):.2f}", "Fig. 3j (GSE26253 row)")
        r = pooled20(d, "all 10 gastric", sc)
        add(f"Fig. 3j pooled row: {sc}", "E20", r.HR, f"{ci2(r.HR, r.lo, r.hi)}, P {fp(r.p)}", "Fig. 3j (pooled row); Fig. 5c (x position); Fig. 8a (−log HR)")
    r = pooled20(d, "all 10 gastric", "CXCL9_SPP1")
    add("Fig. 5c: CXCL9_SPP1 all-nine pooled HR", "E20", r.HR, ci2(r.HR, r.lo, r.hi), "Fig. 5c (CXCL9:SPP1, nine cohorts)")
    r = pooled20(d, "all 10 gastric", "Barrier")
    add("Fig. 5c / 8a: Barrier all-ten pooled HR", "E20", r.HR, ci2(r.HR, r.lo, r.hi), "Fig. 5c; Fig. 8a")
    b = d["E29b_posthoc_adjusted_pooled"]; b = b[b.data == "survival"].set_index("term")
    for t in ["Immune", "CAF_ECM"]:
        r = b.loc[t]
        add(f"E29b mutually adjusted {t} HR, ten gastric (post hoc)", "E29b", math.exp(r.est),
            f"{math.exp(r.est):.2f} ({math.exp(r.lo):.2f}–{math.exp(r.hi):.2f}); P = {fp(r.p)}", "MS [122]; SI [188]")
    return it


# ------------------------------------------------------------------ section 'figure': plotted values, plotted scale
def figure_items(d):
    it = []; add = lambda q, v, w: it.append((q, v, w))
    for sc in ["NetImmune", "CAF_ECM"]:
        r = pc15(d, "GSE26253", sc); add(f"Fig. 3d GSE26253 {sc} HR", r.HR, "Fig. 3d (HR axis)")
        add(f"Fig. 3d GSE26253 {sc} HR lower", r.lo, "Fig. 3d"); add(f"Fig. 3d GSE26253 {sc} HR upper", r.hi, "Fig. 3d")
        for st, lab in [("locked (first 3), E15 reference", "Pooled, first three"), ("locked 7 (primary)", "Pooled, locked 7")]:
            r = pooled20(d, st, sc)
            for col in ["HR", "lo", "hi"]:
                add(f"Fig. 3d {lab} {sc} {col}", r[col], "Fig. 3d (diamond)")
    n = d["E20_null"].set_index("level")
    add("Fig. 3e Locked 7 pooled |log HR| (observed)", abs(n.loc["locked 7, pooled", "observed_logHR"]), "Fig. 3e (|log HR| axis)")
    add("Fig. 3e Locked 7 null 95th percentile", n.loc["locked 7, pooled", "null_p95_abs"], "Fig. 3e")
    d15 = d["E15_deltaC"]; d20 = d["E20_deltaC"].set_index(["cohort", "comparator"])
    for cmp in ["ESTIMATE", "CD8A"]:
        r = d15[(d15.cohort == "GSE26253") & (d15.comparator == cmp)].iloc[0]
        add(f"Fig. 3f GSE26253 ΔC vs {cmp} (gray dot)", r.dC, "Fig. 3f (ΔC axis)")
        r = d20.loc[("locked 7, pooled (inverse variance)", cmp)]
        for col in ["dC", "lo", "hi"]:
            add(f"Fig. 3f / Fig. S10b locked 7 pooled ΔC vs {cmp} {col}", r[col], "Fig. 3f; SI Fig. S10b")
    J = ["CAF_ECM", "TAM_SUP", "APM", "IFN", "TLS_B", "T_CELL", "Immune", "NetImmune"]          # Fig. 3j rows
    F5 = J + ["Barrier", "CXCL9_SPP1"]                                                           # Fig. 5c points
    F8 = J + ["Barrier"]                                                                         # Fig. 8a rows
    for sc in J:
        r = est20(d, "GSE26253", sc); add(f"Fig. 3j GSE26253 {sc} HR", math.exp(r.beta), "Fig. 3j (HR per SD, tile)")
    for sc in F5:
        r = pooled20(d, "all 10 gastric", sc)
        add(f"Fig. 3j / 5c all-ten pooled {sc} HR", r.HR, ("Fig. 3j (pooled row); " if sc in J else "") + "Fig. 5c (x = pooled HR)")
        if sc in ("CAF_ECM", "Immune"):
            add(f"Fig. 5c all-ten pooled {sc} HR lower", r.lo, "Fig. 5c (crosshair)")
            add(f"Fig. 5c all-ten pooled {sc} HR upper", r.hi, "Fig. 5c (crosshair)")
        if sc in F8:
            add(f"Fig. 8a all-ten {sc} −log HR", -math.log(r.HR), "Fig. 8a (−log HR per SD)")
            add(f"Fig. 8a all-ten {sc} CI end −log(upper)", -math.log(r.hi), "Fig. 8a")
            add(f"Fig. 8a all-ten {sc} CI end −log(lower)", -math.log(r.lo), "Fig. 8a")
    r = pooled20(d, "locked 7 (primary)", "NetImmune")
    add("Fig. 8b row text: locked 7 NetImmune HR", r.HR, "Fig. 8b (text 0.70 (0.62–0.80))")
    return it


# ------------------------------------------------------------------ workbook: precision per sheet column
wb = openpyxl.load_workbook(WB, read_only=True, data_only=True)
WBV = {}                                   # sheet -> {column header -> [numeric values]}
for name in set(sum(SHEETS.values(), [])):
    ws = wb[name]; rows_ = list(ws.iter_rows(values_only=True)); hdr = rows_[0]
    cols_ = {}
    for r in rows_[1:]:
        for j, v in enumerate(r):
            if j >= len(hdr) or hdr[j] is None:
                continue
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                cols_.setdefault(hdr[j], []).append(float(v))
            elif isinstance(v, str):
                try: cols_.setdefault(hdr[j], []).append(float(v))
                except ValueError: pass
    WBV[name] = cols_


def decimals(x):
    s = repr(float(x))
    if "e" in s:
        mant, e = s.split("e"); dm = len(mant.split(".")[1]) if "." in mant else 0
        return max(dm - int(e), 0)
    return len(s.split(".")[1].rstrip("0")) if "." in s else 0


COLPREC = {sh: {c: min(max(decimals(x) for x in v), 15) for c, v in cols_.items() if v} for sh, cols_ in WBV.items()}


def wb_lookup(f, col, v):
    """sheet holding the published value v in the column of the same name, and that column's displayed precision
    (the largest number of decimals among its cells; 12 or more = full precision)."""
    for sh in SHEETS.get(f, []):
        dw = COLPREC[sh].get(col)
        if dw is None:
            continue
        for w in WBV[sh][col]:
            ok = abs(w - v) <= 1e-12 * max(1, abs(v)) if dw >= 12 else abs(w - v) <= 0.5 * 10 ** -dw + 1e-12
            if ok:
                return [(sh, col, dw)]
    return []


rows = []
# --- text
pub_t, pri_t, se3_t = text_items(D("pub")), text_items(D("primary")), text_items(D("sens3"))
for a, b, c in zip(pub_t, pri_t, se3_t):
    q, s, v0, r0, w = a
    rows.append(dict(section="text", quantity=q, analysis_script=" + ".join(SCRIPT[x] for x in re.findall(r"E\d+b?", s)),
                     file="", row_key="", column="", published_value=v0, dedup_value=b[2], sens3_value=c[2],
                     abs_change=abs(b[2] - v0), abs_change_sens3=abs(c[2] - v0), reported_old=r0, reported_dedup=b[3], reported_sens3=c[3],
                     reported_precision_change="yes" if b[3] != r0 else "no",
                     reported_precision_change_sens3="yes" if c[3] != r0 else "no", where_reported=w))
# --- figures (plotted scale; > 0.01 flagged)
for a, b, c in zip(figure_items(D("pub")), figure_items(D("primary")), figure_items(D("sens3"))):
    q, v0, w = a
    rows.append(dict(section="figure", quantity=q, analysis_script="(plotted)", file="", row_key="", column="",
                     published_value=v0, dedup_value=b[1], sens3_value=c[1], abs_change=abs(b[1] - v0), abs_change_sens3=abs(c[1] - v0),
                     reported_old=f"{v0:.3f}", reported_dedup=f"{b[1]:.3f}", reported_sens3=f"{c[1]:.3f}",
                     reported_precision_change="yes" if abs(b[1] - v0) > 0.01 else "no",
                     reported_precision_change_sens3="yes" if abs(c[1] - v0) > 0.01 else "no",
                     where_reported=w + "; moves > 0.01 on plotted scale" if abs(b[1] - v0) > 0.01 else w))
# --- workbook / results files: every numeric cell that changes
for f in SHEETS:
    pub = rd(f); pri = rd(f, "primary"); se3 = rd(f, "sens3")
    tag = f.split("_")[0]
    if pub.shape != pri.shape or pub.shape != se3.shape:
        rows.append(dict(section="workbook", quantity=f"{f}: shape changed {pub.shape} -> {pri.shape}/{se3.shape}", analysis_script=SCRIPT[tag], file=f))
        continue
    keycols = [c for c in pub.columns if not pd.api.types.is_numeric_dtype(pub[c])]
    for c in keycols:
        diff = pub[c].fillna("").astype(str) != pri[c].fillna("").astype(str)
        for i in np.where(diff)[0]:
            rows.append(dict(section="workbook", quantity=f"{f} row {i + 1} {c} (text)", analysis_script=SCRIPT[tag], file=f,
                             row_key=" | ".join(f"{k}={pub.at[i, k]}" for k in keycols[:5] if k != c), column=c,
                             reported_old=str(pub.at[i, c]), reported_dedup=str(pri.at[i, c]), reported_sens3=str(se3.at[i, c]),
                             reported_precision_change="yes" if str(pub.at[i, c]) != str(pri.at[i, c]) else "no",
                             reported_precision_change_sens3="yes" if str(pub.at[i, c]) != str(se3.at[i, c]) else "no",
                             where_reported="; ".join(SHEETS[f]) or "results file only"))
    for c in [c for c in pub.columns if c not in keycols]:
        a, b, s3 = pub[c].astype(float).values, pri[c].astype(float).values, se3[c].astype(float).values
        for i in range(len(a)):
            if (np.isnan(a[i]) and np.isnan(b[i]) and np.isnan(s3[i])):
                continue
            ch = abs(b[i] - a[i]); ch3 = abs(s3[i] - a[i])
            if not (ch > 1e-12 or ch3 > 1e-12):
                continue
            hits = wb_lookup(f, c, a[i]) if not np.isnan(a[i]) else []
            if hits:
                sh, col, dw = hits[0]
                prec = "yes" if round(b[i], dw) != round(a[i], dw) else "no"
                prec3 = "yes" if round(s3[i], dw) != round(a[i], dw) else "no"
                where = f"{sh} (column {col}; {dw if dw < 12 else 'full'} decimals)"
            else:
                prec = prec3 = "n/a (not in workbook)"; where = "results file only" if not SHEETS[f] else f"{'; '.join(SHEETS[f])} (value not located)"
            rows.append(dict(section="workbook", quantity=f"{f}: " + " | ".join(f"{k}={pub.at[i, k]}" for k in keycols[:5]) + f" | {c}",
                             analysis_script=SCRIPT[tag], file=f, row_key=" | ".join(f"{k}={pub.at[i, k]}" for k in keycols), column=c,
                             published_value=a[i], dedup_value=b[i], sens3_value=s3[i], abs_change=ch, abs_change_sens3=ch3,
                             reported_precision_change=prec, reported_precision_change_sens3=prec3, where_reported=where))
C = pd.DataFrame(rows)
cols = ["section", "quantity", "analysis_script", "published_value", "dedup_value", "sens3_value", "abs_change", "reported_precision_change",
        "where_reported", "abs_change_sens3", "reported_precision_change_sens3", "reported_old", "reported_dedup", "reported_sens3",
        "file", "row_key", "column"]
C = C[cols]
C.to_csv(A / "E30_dedup_compare.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70); pd.set_option("display.max_rows", 500)
t = C[C.section == "text"]
print(t[["quantity", "reported_old", "reported_dedup", "reported_sens3", "reported_precision_change"]].to_string(index=False))
f_ = C[(C.section == "figure") & (C.abs_change > 0.01)]
print("\nfigure quantities moving > 0.01 (primary):"); print(f_[["quantity", "published_value", "dedup_value", "abs_change"]].to_string(index=False))
w = C[C.section == "workbook"]
print("\nworkbook/results cells changed:", len(w), "; at sheet precision (primary):", (w.reported_precision_change == "yes").sum())
print(w.groupby("file").agg(n=("quantity", "size"), max_change=("abs_change", "max"),
                            prec_yes=("reported_precision_change", lambda x: (x == "yes").sum())).to_string())
