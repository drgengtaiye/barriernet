"""Build MANUSCRIPT_v14.{md,docx} from v14_src.md.

Citations in the source are written [@key] or [@key1; @key2]; they are numbered by first
appearance, rendered as superscripts (ranges collapsed), and the reference list is inserted
at {{REFERENCES}}. Entries shared with MANUSCRIPT_GM.md are taken verbatim from its list.
"""
import re, subprocess, pathlib
HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE.parent
gm = {}
for line in open(OUT / "analysis/results/gm_refs_subset.txt"):
    n, txt = line.split("\\. ", 1); gm[int(n)] = txt.strip()
gm[73] = gm[73].replace("De Carvalho DD, De Carvalho DD.", "De Carvalho DD.")  # duplicated author in GM list
R = {
 "bray2024": "Bray F, Laversanne M, Sung H, Ferlay J, Siegel RL, Soerjomataram I, Jemal A. Global cancer statistics 2022: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. CA Cancer J Clin. 2024;74:229–263.",
 "janjigian2021": "Janjigian YY, Shitara K, Moehler M, Garrido M, Salman P, Shen L, Wyrwicz L, Yamaguchi K, Skoczylas T, Campos Bragagnoli A, et al. First-line nivolumab plus chemotherapy versus chemotherapy alone for advanced gastric, gastro-oesophageal junction, and oesophageal adenocarcinoma (CheckMate 649): a randomised, open-label, phase 3 trial. Lancet. 2021;398:27–40.",
 "kim2018": gm[16],
 "sautes2019": "Sautès-Fridman C, Petitprez F, Calderaro J, Fridman WH. Tertiary lymphoid structures in the era of cancer immunotherapy. Nat Rev Cancer. 2019;19:307–325.",
 "cabrita2020": gm[67], "helmink2020": gm[61], "mariathasan2018": gm[13],
 "chen2017": "Chen DS, Mellman I. Elements of cancer immunity and the cancer–immune set point. Nature. 2017;541:321–330.",
 "chakravarthy2018": gm[73], "lee2025": gm[40],
 "palla2022": "Palla G, Spitzer H, Klein M, Fischer D, Schaar AC, Kuemmerle LB, Rybakov S, Ibarra IL, Holmberg O, Virshup I, Lotfollahi M, Richter S, Theis FJ. Squidpy: a scalable framework for spatial omics analysis. Nat Methods. 2022;19:171–178.",
 "singhal2024": "Singhal V, Chou N, Lee J, Yue Y, Liu J, Chock WK, Lin L, Chang YC, Teo EML, Aow J, Lee HK, Chen KH, Prabhakar S. BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis. Nat Genet. 2024;56:431–441.",
 "varrone2024": "Varrone M, Tavernari D, Santamaria-Martínez A, Walsh LA, Ciriello G. CellCharter reveals spatial cell niches associated with tissue remodeling and cell plasticity. Nat Genet. 2024;56:74–84.",
 "chen2020lda": "Chen Z, Soifer I, Hilton H, Keren L, Jojic V. Modeling multiplexed images with Spatial-LDA reveals novel tissue microenvironments. J Comput Biol. 2020;27:1204–1218.",
 "venet2011": gm[6], "yoshihara2013": gm[11], "kumar2022": gm[12], "tcga2014": gm[30],
 "yoon2020": gm[31], "cristescu2015": gm[32],
 "mun2019": "Mun DG, Bhin J, Kim S, Kim H, Jung JH, Jung Y, et al. Proteogenomic characterization of human early-onset gastric cancer. Cancer Cell. 2019;35:111–124.e10.",
 "yang2013": "Yang W, Soares J, Greninger P, Edelman EJ, Lightfoot H, Forbes S, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. Nucleic Acids Res. 2013;41:D955–D961.",
 "iorio2016": "Iorio F, Knijnenburg TA, Vis DJ, Bignell GR, Menden MP, Schubert M, et al. A landscape of pharmacogenomic interactions in cancer. Cell. 2016;166:740–754.",
 "wolf2018": "Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol. 2018;19:15.",
 "biancalani2021": "Biancalani T, Scalia G, Buffoni L, Avasthi R, Lu Z, Sanger A, et al. Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. Nat Methods. 2021;18:1352–1362.",
 "moran1950": "Moran PAP. Notes on continuous stochastic phenomena. Biometrika. 1950;37:17–23.",
 "hanzelmann2013": gm[7],
 "dersimonian1986": "DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7:177–188.",
 "viechtbauer2010": "Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1–48.",
 "harrell1982": gm[91],
 "zou2005": "Zou H, Hastie T. Regularization and variable selection via the elastic net. J R Stat Soc B. 2005;67:301–320.",
 "simon2011": "Simon N, Friedman J, Hastie T, Tibshirani R. Regularization paths for Cox's proportional hazards model via coordinate descent. J Stat Softw. 2011;39(5):1–13.",
 "ayers2017": gm[1],
 "charoentong2017": "Charoentong P, Finotello F, Angelova M, Mayer C, Efremova M, Rieder D, Hackl H, Trajanoski Z. Pan-cancer immunogenomic analyses reveal genotype-immunophenotype relationships and predictors of response to checkpoint blockade. Cell Rep. 2017;18:248–262.",
 "jiang2018": gm[62],
 "delong1988": "DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics. 1988;44:837–845.",
 "bh1995": gm[34],
 "pertpy": "pertpy: perturbation analysis in Python. https://github.com/scverse/pertpy.",
 "barriernet": "barriernet: marker-level Barrier/NetImmune scoring for spatial, single-cell, bulk and protein data (Python package, version 0.2.0). https://github.com/drgengtaiye/barriernet.",
 "kalluri2016": "Kalluri R. The biology and function of fibroblasts in cancer. Nat Rev Cancer. 2016;16:582–598.",
 "companion": "[Companion manuscript] Immune gene-expression signatures carry treatment-specific information: a matched-null benchmark of 52 published signatures across checkpoint cohorts, randomized trials and 20 cancer types. Submitted.",
 "mantovani2017": "Mantovani A, Marchesi F, Malesci A, Laghi L, Allavena P. Tumour-associated macrophages as treatment targets in oncology. Nat Rev Clin Oncol. 2017;14:399–416.",
 "qi2022": "Qi J, Sun H, Zhang Y, Wang Z, Xun Z, Li Z, et al. Single-cell and spatial analysis reveal interaction of FAP⁺ fibroblasts and SPP1⁺ macrophages in colorectal cancer. Nat Commun. 2022;13:1742.",
 "ballman2015": gm[9],
 "tauriello2018": "Tauriello DVF, Palomo-Ponce S, Stork D, Berenguer-Llergo A, Badia-Ramentol J, Iglesias M, et al. TGFβ drives immune evasion in genetically reconstituted colon cancer metastasis. Nature. 2018;554:538–543.",
}

def fmt(ns):
    ns = sorted(set(ns)); out = []; i = 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1: j += 1
        out.append(str(ns[i]) if j == i else (f"{ns[i]},{ns[j]}" if j == i + 1 else f"{ns[i]}–{ns[j]}"))
        i = j + 1
    return "^" + ",".join(out) + "^"

src = (HERE / "v14_src.md").read_text()
body, tail = src.split("{{REFERENCES}}")
cite = re.compile(r"\[(@[^\]]+)\]")
order = []
for m in cite.finditer(body):
    for k in re.findall(r"@(\w+)", m.group(1)):
        if k not in order: order.append(k)
missing = [k for k in order if k not in R]; assert not missing, missing
unused = [k for k in R if k not in order]; assert not unused, unused
num = {k: i + 1 for i, k in enumerate(order)}
body = cite.sub(lambda m: fmt([num[k] for k in re.findall(r"@(\w+)", m.group(1))]), body)
out = body + "\n\n".join(f"{num[k]}\\. {R[k]}" for k in order) + tail
(OUT / "MANUSCRIPT_v14.md").write_text(out)
subprocess.run(["pandoc", "MANUSCRIPT_v14.md", "-o", "MANUSCRIPT_v14.docx"], check=True, cwd=OUT)
abstract = out.split("## Abstract")[1].split("**Keywords")[0]
main = out.split("## Introduction")[1].split("## Declarations")[0]
print(f"{len(order)} references | abstract {len(abstract.split())} words | "
      f"Introduction-Conclusion {len(main.split())} words (incl. tables)")
