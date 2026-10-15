"""Markdown for the v18 Table 1 rows of the locked extension and for Tables 3 and 4 (randomised comparisons),
generated from the analysis outputs so that no number is retyped. Used by apply_v18_edits.py."""
import math, pathlib
import pandas as pd

R = pathlib.Path(__file__).resolve().parent.parent / "analysis" / "results"
SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"


def fp(p):
    if p < 0.001:
        e = math.floor(math.log10(p)); m = p / 10 ** e
        if round(m, 1) >= 10: m, e = m / 10, e + 1
        sup = "".join(SUP[int(c)] if c.isdigit() else SUP[10] for c in str(e))
        return f"{m:.1f} × 10{sup}"
    if p < 0.01: return f"{p:.4f}" if p < 0.00995 else "0.010"
    if p < 0.1: return f"{p:.3f}"
    return f"{p:.2f}"


def f2(x):
    s = f"{x:.2f}"
    return "−" + s[1:] if s.startswith("-") else s


def ci(lo, hi, signed=False):
    a, b = f2(lo), f2(hi)
    if signed or lo < 0 or hi < 0:
        return f"{a} to {b}"
    return f"{a}–{b}"


def sgn(x):
    s = f2(x)
    return s if s.startswith("−") else "+" + s


def table1_rows():
    pc = pd.read_csv(R / "E20_per_cohort.csv")
    pc = pc[(pc.score == "NetImmune") & (pc.duplicate_rule == "pre-registered exclusions (primary)") & (pc.endpoint == "OS")]
    nl = pd.read_csv(R / "E20_null.csv").set_index("level")
    pl = pd.read_csv(R / "E20_pooled.csv")
    pl = pl[(pl.score == "NetImmune") & (pl.duplicate_rule == "pre-registered exclusions (primary)")]
    rows = ["| **Locked cohorts, pre-registered extension** | | | | | | | |"]
    for c in ["GSE13861", "GSE26899", "GSE26901", "GSE29272"]:
        u = pc[(pc.cohort == c) & (pc.model == "unadjusted")].iloc[0]
        a = pc[(pc.cohort == c) & (pc.model == "adjusted")].iloc[0]
        k = nl.loc[c]
        rows.append(f"| {c} | {int(u.n)} ({int(u.events)}) | {f2(u.HR)} ({ci(u.lo, u.hi)}) | {fp(u.p)} | "
                    f"{f2(a.HR)} ({ci(a.lo, a.hi)}) | {fp(a.p)} | {int(k['rank'])} (*P* = {k.P:.3f}) | — |")
    def prow(label, setname, n, ev, rank):
        u = pl[(pl.set == setname) & (pl.model == "unadjusted")].iloc[0]
        a = pl[(pl.set == setname) & (pl.model == "adjusted")]
        adj = f"{f2(a.iloc[0].HR)} ({ci(a.iloc[0].lo, a.iloc[0].hi)}) | {fp(a.iloc[0].p)}" if len(a) else "— | —"
        return f"| {label} | {n} | {f2(u.HR)} ({ci(u.lo, u.hi)}) | {fp(u.p)} | {adj} | {rank} | — |"
    rows.append(prow("Pooled, added four", "locked (added 4)", "376 (202)", None, "1 (*P* ≤ 0.001)"))
    rows.append(prow("**Pooled, all seven locked**", "locked 7 (primary)", "1,055 (501)", None, "1 (*P* ≤ 0.001)"))
    rows.append(prow("**All ten cohorts, pooled**", "all 10 gastric", "2,152 (999)", None, "—"))
    return "\n".join(rows)


TRIALS = [  # (label, cancer, comparison, endpoint, source, key)
    ("IMmotion150", "Renal", "Atezolizumab ± bevacizumab against sunitinib", "PFS", "E5", "IMmotion150: atezolizumab (+/- bev) vs sunitinib"),
    ("I-SPY2", "Breast (HER2−)", "Paclitaxel + pembrolizumab against paclitaxel", "pCR", "E10", "pembro"),
    ("I-SPY2", "Breast (HER2−)", "Paclitaxel + durvalumab + olaparib against paclitaxel", "pCR", "E10", "durva"),
    ("CheckMate 025", "Renal", "Nivolumab against everolimus", "PFS", "E21", "CheckMate 025"),
    ("JAVELIN Renal 101", "Renal", "Avelumab + axitinib against sunitinib", "PFS", "E21", "JAVELIN Renal 101"),
    ("NeoTRIP", "Breast (TNBC)", "Chemotherapy + atezolizumab against chemotherapy", "pCR", "E21", "NeoTRIP"),
    ("PRINCE", "Pancreas", "Nivolumab-containing arms against sotigalimab, all with chemotherapy", "PFS", "E5", "PRINCE: nivolumab-containing vs sotigalimab (all chemo)"),
    ("DUTRENEO", "Bladder", "Durvalumab + tremelimumab against cisplatin-based chemotherapy", "pCR", "E21", "DUTRENEO"),
]


def trial_effect(src, key, score):
    if src == "E5":
        t = pd.read_csv(R / "E5_trials.csv"); r = t[(t.trial == key) & (t.score == score)].iloc[0]
        return "HR", r.int_HR, r.int_lo, r.int_hi, r.int_p, int(r.n_ici), int(r.n_ctrl)
    if src == "E10":
        t = pd.read_csv(R / "E10_ispy2_interactions.csv")
        r = t[t.cohort.str.contains(key, case=False) & (t.score == score)].iloc[0]
        return "OR", r.HRadj_OR, r.HRadj_lo, r.HRadj_hi, r.HRadj_p, int(r.n_ici), int(r.n_ctrl)   # hormone-receptor adjusted
    t = pd.read_csv(R / "E21_interactions.csv")
    r = t[(t.trial == key) & (t.score == score) & (t.role.isin(["primary", "sensitivity"]))]
    r = r[r.endpoint.isin(["PFS", "pCR"])].iloc[0]
    kind = "HR" if r.model == "Cox" else "OR"
    return kind, r.int_est, r.int_lo, r.int_hi, r.int_p, int(r.n_ici), int(r.n_ctrl)


def table3():
    hdr = ("| Trial | Cancer | Checkpoint arm against control | Endpoint | n | NetImmune interaction (95% CI); *P* | CAF_ECM interaction (95% CI); *P* |\n"
           "|---|---|---|---|---|---|---|")
    rows = [hdr, "| **Primary** | | | | | | |"]
    for i, (lab, can, comp, ep, src, key) in enumerate(TRIALS):
        if i == 6: rows.append("| **Sensitivity** | | | | | | |")
        k, e, lo, hi, p, n1, n0 = trial_effect(src, key, "NetImmune")
        k2, e2, lo2, hi2, p2, _, _ = trial_effect(src, key, "CAF_ECM")
        rows.append(f"| {lab} | {can} | {comp} | {ep} | {n1} / {n0} | {k} {f2(e)} ({ci(lo, hi)}); {fp(p)} | {k2} {f2(e2)} ({ci(lo2, hi2)}); {fp(p2)} |")
    return "\n".join(rows)


def table4():
    m = pd.read_csv(R / "E21_meta.csv")
    a = m[m.set == "primary: 6 comparisons"].set_index("score")
    b = m[m.set == "+ PRINCE + DUTRENEO (8)"].set_index("score")
    names = [("CAF_ECM", "CAF_ECM"), ("TAM_SUP", "TAM_SUP"), ("Barrier", "Barrier"), ("IFN", "IFN"), ("Immune", "Immune"),
             ("NetImmune", "NetImmune"), ("C1Q_arm", "C1Q arm"), ("SPP1_arm", "SPP1 arm"), ("CXCL9_SPP1", "CXCL9:SPP1 polarity")]
    rows = ["| Score | Pooled interaction, six primary comparisons (95% CI) | *P* | Favouring checkpoint arm | Stouffer *P* | With PRINCE and DUTRENEO (95% CI); *P* |",
            "|---|---|---|---|---|---|"]
    for k, lab in names:
        r, s = a.loc[k], b.loc[k]
        rows.append(f"| {lab} | {sgn(r.est)} ({ci(r.lo, r.hi, True)}) | {fp(r.p)} | {int(r.k_benefit)}/{int(r.k)} | {fp(r.stouffer_p)} | "
                    f"{sgn(s.est)} ({ci(s.lo, s.hi, True)}); {fp(s.p)} |")
    return "\n".join(rows)


if __name__ == "__main__":
    print(table1_rows()); print(); print(table3()); print(); print(table4())
