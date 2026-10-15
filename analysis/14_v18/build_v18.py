"""Build MANUSCRIPT_v18 and SUPPLEMENTARY_v18 (.md, .docx) from v18_src.md and v18_supp_src.md.
Source order: apply_v18_edits.py -> apply_v18_fixes.py (-> E22 edits) -> build_v18.py -> make_v18_workbook.py -> build_final.py v18.

Compared with build_lean.py:
- supplementary items are cited as {SF:key} / {ST:key} and numbered by first citation in the main
  text, then in the supplement; consecutive numbers collapse to ranges (S3–S5);
- figures are referenced as ![](FIG:key); the files are copied to figures/v18/ under their numbers;
- supplementary legends and supplementary tables are reordered by number;
- references are numbered separately in each document, as before (texts in refs.py plus R18 below).
Writes figures/v18/numbering.json, which make_v18_workbook.py reads.
"""
import json, re, shutil, subprocess, pathlib
from refs import R, HERE, OUT

R18 = {
    "motzer2015": "Motzer RJ, Escudier B, McDermott DF, George S, Hammers HJ, Srinivas S, et al. Nivolumab versus everolimus in advanced renal-cell carcinoma. N Engl J Med. 2015;373:1803–1813.",
    "motzer2019": "Motzer RJ, Penkov K, Haanen J, Rini B, Albiges L, Campbell MT, et al. Avelumab plus axitinib versus sunitinib for advanced renal-cell carcinoma. N Engl J Med. 2019;380:1103–1115.",
    "motzer2020": "Motzer RJ, Robbins PB, Powles T, Albiges L, Haanen JB, Larkin J, et al. Avelumab plus axitinib versus sunitinib in advanced renal cell carcinoma: biomarker analysis of the phase 3 JAVELIN Renal 101 trial. Nat Med. 2020;26:1733–1741.",
    "gianni2022": "Gianni L, Huang CS, Egle D, Bermejo B, Zamagni C, Thill M, et al. Pathologic complete response (pCR) to neoadjuvant treatment with or without atezolizumab in triple-negative, early high-risk and locally advanced breast cancer: NeoTRIP Michelangelo randomized study. Ann Oncol. 2022;33:534–543.",
    "dugo2026": "Dugo M, Huang CS, Egle D, Bermejo B, Seitz RS, Nielsen TJ, et al. Longitudinal transcriptomic profiling identifies predictors of response to neoadjuvant chemoimmunotherapy in triple-negative breast cancer: results from the NeoTRIPaPDL1 trial. Ann Oncol. 2026;37:1278–1292.",
    "grande2026": "Grande E, Sibai M, Grases D, Andrada E, Reig O, Escobosa M, et al. Spatial architecture contributes to failure of bulk biomarker-guided neoadjuvant immunotherapy selection in bladder cancer: the DUTRENEO study. Cell Rep Med. 2026;7:102878.",
    "oh2018": "Oh SC, Sohn BH, Cheong JH, Kim SB, Lee JE, Park KC, et al. Clinical and genomic landscape of gastric cancer with a mesenchymal phenotype. Nat Commun. 2018;9:1777.",
    "cho2011": "Cho JY, Lim JY, Cheong JH, Park YY, Yoon SL, Kim SM, et al. Gene expression signature-based prognostic risk score in gastric cancer. Clin Cancer Res. 2011;17:1850–1857.",
    "wang2013": "Wang G, Hu N, Yang HH, Wang L, Su H, Wang C, et al. Comparison of global gene expression of gastric cardia and noncardia cancers from a high-risk population in China. PLoS One. 2013;8:e63826.",
    "li2023": "Li G, Huo D, Guo N, Li Y, Ma H, Liu L, et al. Integrating multiple machine learning algorithms for prognostic prediction of gastric cancer based on immune-related lncRNAs. Front Genet. 2023;14:1106724.",
    "sun2022": "Sun K, Xu R, Ma F, Yang N, Li Y, Sun X, et al. scRNA-seq of gastric tumor shows complex intercellular interaction with an alternative T cell exhaustion trajectory. Nat Commun. 2022;13:4943.",
    "an2024": "An M, Mehta A, Min BH, Heo YJ, Wright SJ, Parikh M, et al. Early immune remodeling steers clinical response to first-line chemoimmunotherapy in advanced gastric cancer. Cancer Discov. 2024;14:766–785.",
    "wright2026": "Wright SJ, Kang S, An M, Heo YJ, Parikh M, Bi L, et al. Differential T cell clonal dynamics underlie outcomes to frontline chemoimmunotherapy in advanced gastric cancer. Cell Rep Med. 2026;7:102910.",
    "mamdani2021": "Mamdani H, Schneider B, Perkins SM, Burney HN, Kasi PM, Abushahin LI, et al. A phase II trial of adjuvant durvalumab following trimodality therapy for locally advanced esophageal and gastroesophageal junction adenocarcinoma: a Big Ten Cancer Research Consortium study. Front Oncol. 2021;11:736620.",
    "chen2018fastp": "Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34:i884–i890.",
    "patro2017": "Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017;14:417–419.",
    "kang2022": "Kang B, Camps J, Fan B, Jiang H, Ibrahim MM, Hu X, et al. Parallel single-cell and bulk transcriptome analyses reveal key features of the gastric tumor microenvironment. Genome Biol. 2022;23:265.",
}
REFS = {**R, **R18}

F = OUT / "figures"
LEAN = {1: "sc_atlas", 2: "sc_programmes", 3: "sc_composition", 4: "sc_fractions", 5: "sc_share",
        6: "meta_forest", 7: "fitted_models", 8: "upgrades", 9: "locked", 10: "icb11", 11: "tcga",
        12: "perfect_metabric", 13: "visium", 14: "visium_sens", 15: "barrierindex", 16: "moran",
        17: "transport", 18: "imc", 19: "proteomics", 20: "drugs"}
FIGS = {k: F / f"lean/supp/FigS{n}" for n, k in LEAN.items()}          # stem; .png and .pdf both copied
FIGS.update({
    "locked": F / "Fig_E20_extended_replication",       # supersedes the three-cohort locked figure
    "sc_indep": F / "Fig_E17_attribution",
    "robust": F / "Fig_E19_robustness",
    "trials_detail": F / "Fig_E10_trials",               # the three-comparison figure, kept for the within-arm panels
    "gastric_icb": F / "Fig_E22_gastric_icb",
    "main1": F / "v18/src/Fig1", "main2": F / "lean/Fig2", "main3": F / "lean/Fig3", "main4": F / "lean/Fig4",
    "main5": F / "Fig_E21_trials",
})
TOKEN = re.compile(r"\{(SF|ST):([a-z0-9_]+)\}")


def fmt_cite(ns):
    ns = sorted(set(ns)); out = []; i = 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1: j += 1
        out.append(str(ns[i]) if j == i else (f"{ns[i]},{ns[j]}" if j == i + 1 else f"{ns[i]}–{ns[j]}"))
        i = j + 1
    return "^" + ",".join(out) + "^"


def number_items(main_body, supp_text):
    num = {"SF": {}, "ST": {}}
    for text in (main_body, supp_text):
        for m in TOKEN.finditer(text):
            d = num[m.group(1)]
            if m.group(2) not in d: d[m.group(2)] = len(d) + 1
    return num


def resolve(text, num):
    text = TOKEN.sub(lambda m: f"S{num[m.group(1)][m.group(2)]}", text)
    def collapse(m):                                    # "S3, S4, S5" -> "S3–S5"; "S3, S4" stays
        ns = [int(x) for x in re.findall(r"S(\d+)", m.group(0))]
        if ns != sorted(ns) or len(ns) < 3: return m.group(0)
        out, i = [], 0
        while i < len(ns):
            j = i
            while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1: j += 1
            if j == i + 1: out += [f"S{ns[i]}", f"S{ns[j]}"]
            else: out.append(f"S{ns[i]}" if j == i else f"S{ns[i]}–S{ns[j]}")
            i = j + 1
        return ", ".join(out) if " and " not in m.group(0) else ", ".join(out[:-1]) + (" and " + out[-1] if len(out) > 1 else out[0])
    return re.sub(r"S\d+(?:(?:, | and )S\d+){2,}", collapse, text)


def reorder_blocks(section, kind, num):
    """Split a section into blocks that each start with a FIG/legend or table caption, and sort them by number."""
    if kind == "SF":
        parts = re.split(r"(?m)^(?=!\[\]\(FIG:)", section)
        key = lambda b: num["SF"][re.search(r"FIG:([a-z0-9_]+)", b).group(1)]
    else:
        parts = re.split(r"(?m)^(?=\*\*Supplementary Table \{ST:)", section)
        key = lambda b: num["ST"][re.search(r"\{ST:([a-z0-9_]+)\}", b).group(1)]
    head = [p for p in parts if not (("FIG:" in p[:12]) if kind == "SF" else p.startswith("**Supplementary Table {ST:"))]
    blocks = [re.sub(r"(?m)^---\s*$\n?", "", p).rstrip() + "\n\n" for p in parts if p not in head]   # section rules move to the end
    return "".join(head) + "".join(sorted(blocks, key=key)) + "---\n\n"


def sort_guide(text, num):
    """Sort the bullets of '### Supplementary Tables in the workbook' (with their indented sub-lines) by number."""
    i = text.index("### Supplementary Tables in the workbook")
    j = text.index("\n---", i)
    sec = text[i:j]
    lines = sec.split("\n")
    k = next(n for n, l in enumerate(lines) if l.startswith("- **Supplementary Table {ST:"))
    items, cur = [], None
    for l in lines[k:]:
        if l.startswith("- **Supplementary Table {ST:"):
            cur = [l]; items.append(cur)
        elif cur is not None and (l.startswith("  ") or l == "") and l.strip():
            cur.append(l)
    items.sort(key=lambda it: num["ST"][re.search(r"\{ST:([a-z0-9_]+)\}", it[0]).group(1)])
    new = "\n".join(lines[:k]) + "\n" + "\n".join("\n".join(it) for it in items) + "\n"
    return text[:i] + new + text[j:]


def place_figures(text, num, sub):
    def rep(m):
        k = m.group(1)
        src = FIGS[k]
        if k.startswith("main"):
            dst = F / "v18" / f"Fig{k[4:]}"
        else:
            dst = F / "v18" / "supp" / f"FigS{num['SF'][k]}"
        dst.parent.mkdir(parents=True, exist_ok=True)
        for ext in (".png", ".pdf"):
            shutil.copy(src.with_suffix(ext), dst.with_suffix(ext))
        return f"![]({dst.relative_to(OUT)}.png){{width=100%}}"
    return re.sub(r"!\[\]\(FIG:([a-z0-9_]+)\)", rep, text)


def cite(body):
    c = re.compile(r"\[(@[^\]]+)\]")
    order = []
    for m in c.finditer(body):
        for k in re.findall(r"@(\w+)", m.group(1)):
            if k not in order: order.append(k)
    missing = [k for k in order if k not in REFS]; assert not missing, missing
    n = {k: i + 1 for i, k in enumerate(order)}
    return c.sub(lambda m: fmt_cite([n[k] for k in re.findall(r"@(\w+)", m.group(1))]), body), order


def build():
    main = (HERE / "v18_src.md").read_text()
    supp = (HERE / "v18_supp_src.md").read_text()
    mbody = main.split("{{REFERENCES}}")[0]
    num = number_items(mbody, supp)
    # every legend/table must be cited, every citation must have a legend/table
    legend_keys = set(re.findall(r"FIG:([a-z0-9_]+)", supp))
    table_keys = set(re.findall(r"\*\*Supplementary Table \{ST:([a-z0-9_]+)\}", supp)) | \
                 set(re.findall(r"- \*\*Supplementary Table \{ST:([a-z0-9_]+)\}", supp))
    assert set(num["SF"]) == legend_keys, (set(num["SF"]) ^ legend_keys)
    assert set(num["ST"]) <= table_keys | set(), set(num["ST"]) - table_keys
    i, j = supp.index("## Supplementary Figures"), supp.index("## Supplementary references")
    supp = supp[:i] + reorder_blocks(supp[i:j], "SF", num) + supp[j:]
    i, j = supp.index("## Supplementary Tables"), supp.index("### Supplementary Tables in the workbook")
    supp = supp[:i] + reorder_blocks(supp[i:j], "ST", num).rstrip("-\n") + "\n\n" + supp[j:]
    supp = sort_guide(supp, num)
    out_files = {}
    for name, text in [("MANUSCRIPT_v18", main), ("SUPPLEMENTARY_v18", supp)]:
        text = place_figures(text, num, name)
        text = resolve(text, num)
        body, tail = text.split("{{REFERENCES}}")
        body, order = cite(body)
        text = body + "\n\n".join(f"{k + 1}\\. {REFS[r]}" for k, r in enumerate(order)) + tail
        assert not TOKEN.search(text) and "FIG:" not in text and "[@" not in text
        (OUT / f"{name}.md").write_text(text)
        subprocess.run(["pandoc", f"{name}.md", "-o", f"{name}.docx"], check=True, cwd=OUT)
        out_files[name] = (text, order)
    (F / "v18").mkdir(exist_ok=True)
    (F / "v18" / "numbering.json").write_text(json.dumps(num, indent=1))
    m, o = out_files["MANUSCRIPT_v18"]
    words = lambda t: len(t.split())
    abstract = m.split("## Abstract")[1].split("**Keywords")[0]
    body = m.split("## Introduction")[1].split("## Declarations")[0]
    print(f"main: {len(o)} references | abstract {words(abstract) - 4} words (excl. 4 labels) | "
          f"Introduction–Conclusion {words(body)} words | supp figs {len(num['SF'])} | supp tables {len(num['ST'])}")
    print("supplement:", len(out_files["SUPPLEMENTARY_v18"][1]), "references")
    print("numbering:", {k: v for k, v in num["SF"].items()})


if __name__ == "__main__":
    build()
