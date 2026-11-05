"""Writes the v48 figure scripts from the previous ones: TIS_Damotte removed so the counts are of 51 distinct signatures (protocol
analysis/E46_prereg.md part E), "profile" replaced by "gene list" in Fig. 2, the simple T-cell mean and the MSS/EBV-negative patients added."""
from pathlib import Path
H = Path(__file__).resolve().parent; C = H.parent; DUP = '"TIS_Damotte"'


def patch(src, dst, subs):
    s = (C / src).read_text()
    for a, b in subs:
        assert s.count(a) >= 1, (src, a)
        s = s.replace(a, b)
    (H / dst).write_text(s); print("wrote", dst)


patch("v43/make_fig2.py", "make_fig2.py", [
    ("T-cell-inflamed expression profile", "T-cell-inflamed gene list"),
    ('f"profile\\n{PR.est:+.2f}"', 'f"gene list\\n{PR.est:+.2f}"'),
    ('a3.set_ylabel("Profile − gene")', 'a3.set_ylabel("Difference per SD")'),
    ('"median signature', '"typical signature'),
    ('a3.set_xlim(-1, len(PQ)); a3.set_xticks([])', 'R3_ = pd.read_csv(R / "E46_A_paired.csv"); R3_ = R3_[(R3_.set == "six (primary)") & (R3_.reference == "R3") & (R3_.score == SIG)].iloc[0]\n'
     'a3.plot([len(PQ) + 1] * 2, [R3_.lo, R3_.hi], color=INK, lw=0.6); a3.plot(len(PQ) + 1, R3_["diff"], "D", ms=2.8, mfc="white", mec=INK, mew=0.7, zorder=3)\n'
     'a3.annotate("T-cell mean", xy=(len(PQ) + 1, R3_.lo), xytext=(len(PQ) - 3.5, -0.125), fontsize=6, color=INK, ha="center", va="center", arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4))\n'
     'rec("a3 minus T-cell mean", f"{R3_[\'diff\']:+.3f} ({R3_.lo:.3f} to {R3_.hi:.3f}), P = {R3_.p:.3g}", "E46_A_paired.csv")\n'
     'a3.set_xlim(-1, len(PQ) + 2.2); a3.set_xticks([])'),
    ('CC = pd.read_csv(R / "E40_concordance.csv", index_col=0)\naxH', f'CC = pd.read_csv(R / "E40_concordance.csv", index_col=0).drop({DUP})\naxH')])

patch("v47/make_fig3.py", "make_fig3.py", [
    ('L = P[P.kind == "library"].copy()', f'L = P[(P.kind == "library") & (P.signature != {DUP})].copy(); assert len(L) == 51'),
    ('PER = PER[PER.variant == "primary"]', f'PER = PER[(PER.variant == "primary") & (PER.signature != {DUP})]'),
    ('"Of 52 published signatures"', '"Of 51 distinct signatures"'),
    ('ROWS = [("favor the checkpoint arm,\\nP < 0.05", nfav, FAV), ("favor the control arm,\\nP < 0.05", nag, UNFAV),\n        ("beyond random genes\\nof their size", nbey, INK), ("above CD8A, P < 0.05", nabove, IMMUNE),\n        ("positive difference\\nfrom CD8A, n.s.", npt, "#9CC3E6"), ("below CD8A, P < 0.05", nbelow, GREY)]',
     'XC = pd.read_csv(R / "E45_B_coexpression.csv"); XC = XC[XC.signature != "TIS_Damotte"]; ncox = int((XC.emp_p_coex < 0.05).sum()); assert ncox == 15\n'
     'R3C = pd.read_csv(R / "E46_A_counts.csv"); R3C = R3C[(R3C.set == "six (primary)") & (R3C.reference == "R3")].iloc[0]\n'
     'ROWS = [("favor checkpoint arm", nfav, FAV), ("favor control arm", nag, UNFAV), ("beyond size-matched", nbey, INK), ("beyond coexpression", ncox, INK),\n        ("above CD8A", nabove, IMMUNE), ("above T-cell mean", int(R3C.above), IMMUNE), ("below T-cell mean", int(R3C.below), GREY), ("below CD8A", nbelow, GREY)]'),
    ('f"{v}/52"', 'f"{v}/51"'),
    ('fontweight="bold" if lab.startswith("above CD8A, P") else "normal"', 'fontweight="bold" if lab.startswith("above ") else "normal"'),
    ('axB.set_xlabel("Signatures")', 'axB.set_xlabel("Signatures (P < 0.05)")'),
    ('of 52 beyond their reference', 'of 51 beyond their reference'),
    ('f"{nb}/52"', 'f"{nb}/51"'),
    ('median upper limit {Q4[\'median_upper\']:.3f}', 'median upper limit {L.cd8a_hi.median():.3f}'),
    ('27 scored', '26 scored')])

patch("v43/make_fig4.py", "make_fig4.py", [
    ('Q3 = pd.read_csv(R / "E33_Q3_by_reference.csv").set_index("reference")',
     'Q3 = pd.read_csv(R / "E33_Q3_by_reference.csv").set_index("reference")\n'
     '_L0 = pd.read_csv(R / "E32_primary.csv"); _LIB = set(_L0[_L0.kind == "library"].signature) - {"TIS_Damotte"}; assert len(_LIB) == 51\n'
     '_PQ = pd.read_csv(R / "E33_paired.csv"); _PQ = _PQ[(_PQ.set == "six (primary)") & _PQ.signature.isin(_LIB)]\n'
     'Q3["exceeds"] = _PQ[(_PQ["diff"] > 0) & (_PQ.p < 0.05)].groupby("reference").size().reindex(Q3.index, fill_value=0)'),
    ('L = pd.read_csv(R / "E32_primary.csv"); L = L[L.kind == "library"]', 'L = pd.read_csv(R / "E32_primary.csv"); L = L[(L.kind == "library") & (L.signature != "TIS_Damotte")]'),
    ('"52 signatures"', '"51 signatures"'), ('"more than 2 of 52"', '"more than 2 of 51"'), ('of 52)', 'of 51)'),
    ('GB = pd.read_csv(R / "E35_by_gene.csv");', 'GB = pd.read_csv(R / "E35_by_gene.csv"); _E35 = pd.read_csv(R / "E35_paired.csv"); _E35 = _E35[_E35.signature.isin(_LIB) & (_E35["diff"] > 0) & (_E35.p < 0.05)].groupby(["set", "reference"]).size()\n'
     'GB["exceeded_by"] = [int(_E35.get((s_, g_), 0)) for s_, g_ in zip(GB.set, GB.gene)];')])

patch("v44/make_fig5.py", "make_fig5.py", [
    ('L = P[P.kind == "library"].copy()', f'L = P[(P.kind == "library") & (P.signature != {DUP})].copy(); assert len(L) == 51'),
    ('A = pd.read_csv(R / "E32_pooled_all_sets.csv"); NV = pd.read_csv(R / "E32_native.csv")', f'A = pd.read_csv(R / "E32_pooled_all_sets.csv"); A = A[A.signature != {DUP}]; NV = pd.read_csv(R / "E32_native.csv"); NV = NV[NV.signature != {DUP}]'),
    ('G = pd.read_csv(R / "E33_pooled.csv"); PR = pd.read_csv(R / "E33_paired.csv")', f'G = pd.read_csv(R / "E33_pooled.csv"); G = G[G.name != {DUP}]; PR = pd.read_csv(R / "E33_paired.csv"); PR = PR[PR.signature != {DUP}]'),
    ('P34 = pd.read_csv(R / "E34_pooled.csv")', f'P34 = pd.read_csv(R / "E34_pooled.csv"); P34 = P34[P34.name != {DUP}]'),
    ('(of 52)', '(of 51)'), ('of 52 signatures', 'of 51 signatures'), ('rec("e adjusted/unadjusted (18 significant)"', 'rec("e adjusted/unadjusted (significant)"')])

patch("v47/make_fig6.py", "make_fig6.py", [
    ('E20 = E20[(E20.set == "all 10 gastric")', 'E7 = E20[(E20.set == "locked 7 (primary)") & (E20.model == "unadjusted") & (E20.duplicate_rule == "pre-registered exclusions (primary)")].set_index("score")\nE20 = E20[(E20.set == "all 10 gastric")'),
    ('    r = E20.loc[k]; ci_line(axA, i, r.HR, r.lo, r.hi, XP, col, mk="D", filled=r.p < 0.05)\n',
     '    r7 = E7.loc[k]; ci_line(axA, i - 0.17, r7.HR, r7.lo, r7.hi, XP, col, mk="D", filled=r7.p < 0.05)\n    r = E20.loc[k]; ci_line(axA, i + 0.17, r.HR, r.lo, r.hi, XP, col, mk="o", filled=r.p < 0.05)\n    rec(f"a {k} seven", f"HR {r7.HR:.3f} ({r7.lo:.3f}-{r7.hi:.3f}), P = {r7.p:.3g}", "E30_dedup_runs/primary/E20_pooled.csv")\n'),
    ('axA.text(XP[1] * 1.03, i, f"{r.HR:.2f}", fontsize=FN, color=SLATE, va="center", ha="left")', 'axA.text(XP[1] * 1.03, i - 0.17, f"{r7.HR:.2f}", fontsize=FN, color=SLATE, va="center", ha="left"); axA.text(XP[1] * 1.03, i + 0.17, f"{r.HR:.2f}", fontsize=FN, color=SLATE, va="center", ha="left")'),
    ('axA.set_xlabel("HR per SD, ten gastric cohorts (2,144 patients)")', 'axA.set_xlabel("HR per SD"); axA.text(0.615, 4.0, "diamond, seven\\ncircle, all ten", fontsize=FN, color=SLATE, va="center", ha="left", linespacing=1.05)'),
    ('XP = (0.6, 1.6)', 'XP = (0.6, 1.75)'),
    ('("fe2", "Two gastric"), ("adj2", "Two gastric, MSI-H/EBV-adjusted")]', '("fe2", "Two gastric"), ("adj2", "Two gastric, MSI-H/EBV-adjusted"), ("neg2", "Two gastric, MSS/EBV-negative (73; 28)")]\nC46 = pd.read_csv(R / "E46_C_pooled.csv")'),
    ('        elif k == "adj2":\n', '        elif k == "neg2":\n            r = C46[C46.name == nm].iloc[0]; e, lo, hi, p = r.est, r.lo, r.hi, r.p\n        elif k == "adj2":\n'),
    ('LIB = E[E.kind == "library"]', 'LIB = E[(E.kind == "library") & (E.index != "TIS_Damotte")]; assert len(LIB) == 51'),
    ('assert (n_un, n_adj) == (17, 1)', 'assert (n_un, n_adj) == (16, 1)'), ('{n_un} of 52 unadjusted', '{n_un} of 51 unadjusted'),
    ('assert len(q) == 52 and', 'assert len(q) == 51 and'), ('assert abs(rho.statistic - 0.760) < 0.001', 'assert abs(rho.statistic - 0.758) < 0.001'), ('(52 signatures)', '(51 signatures)')])

patch("v47/make_figS36.py", "make_figS36.py", [
    ('fig = plt.figure(figsize=(180 * MM, 128 * MM))', 'fig = plt.figure(figsize=(180 * MM, 186 * MM))'),
    ('P = P[(P.set == "six (primary)") & (P.family != "composite")]', 'P = P[(P.set == "six (primary)") & (P.family != "composite") & (P.score != "TIS_Damotte")]'),
    ('assert len(j) == 52', 'assert len(j) == 51'), ('assert (n_up, n_dn) == (8, 2)', 'assert (n_up, n_dn) == (7, 2)'),
    ('X = pd.read_csv(R / "E45_B_coexpression.csv")', 'X = pd.read_csv(R / "E45_B_coexpression.csv"); X = X[X.signature != "TIS_Damotte"]; assert len(X) == 51'),
    ('assert (nA, nB) == (26, 16)', 'assert (nA, nB) == (25, 15)'), ('{nA} of 52, size only', '{nA} of 51, size only'),
    ('K = S.ax_mm(fig, 0, 122, 180, 5)', '''# ---- e  steadiness across tumor type
SD = pd.read_csv(R / "E46_A_steadiness.csv")
axE = S.ax_mm(fig, 16, 130, 62, 42); title(axE, "Renal against breast comparisons", -0.2); LE = (-0.2, 0.6)
axE.plot(LE, LE, color=LIGHT, lw=0.5, zorder=0)
g_ = SD[SD.group == "panel gene"]; s_ = SD[SD.group == "inflamed signature"]
axE.plot(g_.renal, g_.breast, "o", ms=2.6, mfc="white", mec=OTHER, mew=0.6, zorder=2, label="panel gene"); axE.plot(s_.renal, s_.breast, "o", ms=2.6, mfc=IMMUNE, mec="white", mew=0.3, zorder=3)
for nm, lab, dx, dy in [("R3", "three-gene mean", 0.09, 0.05), ("R13", "13-gene mean", -0.155, -0.27), ("CD8A", "CD8A", 0.09, 0.0), ("CXCL10", "CXCL10", 0.10, 0.06)]:
    q_ = SD[SD.score == nm].iloc[0]; axE.plot(q_.renal, q_.breast, "D", ms=3.0, mfc="white" if nm in ("CD8A", "CXCL10") else INK, mec=INK, mew=0.7, zorder=4)
    axE.annotate(lab, (q_.renal, q_.breast), (q_.renal + dx, q_.breast + dy), fontsize=FN, color=INK, va="center", ha="left", arrowprops=dict(arrowstyle="-", color=GREY, lw=0.4, shrinkA=0, shrinkB=2.2))
    rec(f"e {nm}", f"renal {q_.renal:+.3f}; breast {q_.breast:+.3f}; D {q_.D:.3f}", "E46_A_steadiness.csv")
md_g, md_s = g_.D.median(), s_.D.median(); rec("e median D", f"panel genes {md_g:.3f}; inflamed signatures {md_s:.3f}", "E46_A_steadiness.csv")
axE.text(0.98, 0.03, f"median |renal − breast|:\\n{md_g:.2f}, 25 panel genes (open)\\n{md_s:.2f}, 19 inflamed signatures", transform=axE.transAxes, fontsize=FSS, color=SLATE, va="bottom", ha="right", linespacing=1.1)
axE.set_xlim(*LE); axE.set_ylim(*LE); ticks(axE, [-0.2, 0, 0.2, 0.4, 0.6]); ticks(axE, [-0.2, 0, 0.2, 0.4], "y"); S.rule(axE, 0); S.hrule(axE, 0); axE.set_xlabel("Pooled interaction, three renal"); axE.set_ylabel("Pooled interaction, three breast")
S.letter_at(fig, 2, 124, "e")
# ---- f  choosing a gene in advance
BP = pd.read_csv(R / "E46_B_pooled.csv").iloc[0]; BH = pd.read_csv(R / "E46_B_held_out.csv"); BH = BH[BH.candidates.str.startswith("25")]; RP_ = pd.read_csv(R / "E46_A_reference_pooled.csv"); RP_ = RP_[RP_.set == "six (primary)"].set_index("score")
E32_ = pd.read_csv(R / "E32_primary.csv").set_index("signature").loc["Tcell-inflamed_Ayers"]
RW = [("Best panel gene, apparent (CXCL13)", BP.apparent_best, None, None, "D"), ("Best gene chosen on the other five", BP.est, BP.lo, BP.hi, "D"), ("CD8A, PRF1, GZMA mean", RP_.loc["R3", "est"], RP_.loc["R3", "lo"], RP_.loc["R3", "hi"], "o"),
      ("13 T-cell gene mean", RP_.loc["R13", "est"], RP_.loc["R13", "lo"], RP_.loc["R13", "hi"], "o"), ("CD8A", RP_.loc["CD8A", "est"], RP_.loc["CD8A", "lo"], RP_.loc["CD8A", "hi"], "o"), ("T-cell-inflamed gene list", E32_.est, E32_.lo, E32_.hi, "o")]
axF = S.ax_mm(fig, 134, 130, 42, 42); title(axF, "Choosing a gene in advance", -1.05)
for i, (lab, e_, lo_, hi_, mk) in enumerate(RW):
    if lo_ is not None:
        axF.plot([lo_, hi_], [i, i], color=INK, lw=0.6)
    axF.plot(e_, i, mk, ms=2.8, mfc=INK if lo_ is not None else "white", mec=INK, mew=0.7, zorder=3); axF.text(-0.13, i, lab, fontsize=FN, color=INK, ha="right", va="center")
    rec(f"f {lab}", f"{e_:+.3f}" + (f" ({lo_:.3f} to {hi_:.3f})" if lo_ is not None else ""), "E46_B_pooled.csv; E46_A_reference_pooled.csv; E32_primary.csv")
axF.axhline(1.5, color=LIGHT, lw=0.5); axF.set_ylim(len(RW) - 0.5, -0.5); axF.set_yticks([]); axF.spines["left"].set_visible(False); axF.set_xlim(-0.1, 0.5); ticks(axF, [0, 0.2, 0.4]); S.rule(axF, 0); axF.set_xlabel("Interaction per SD, six comparisons")
rec("f chosen genes", "; ".join(f"{a_}: {b_}" for a_, b_ in zip(BH.held_out, BH.gene)), "E46_B_held_out.csv"); S.letter_at(fig, 96, 124, "f")
K = S.ax_mm(fig, 0, 180, 180, 5)'''),
    ('S.letter_at(fig, 96, 5, "b")', 'S.letter_at(fig, 96, 5, "b")')])
