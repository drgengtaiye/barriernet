"""Build MANUSCRIPT_v17.{md,docx} from v17_src.md.

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
for line in open(OUT / "analysis/results/gm_refs_v15.txt"):
    n, txt = line.split("\\. ", 1); gm[int(n)] = txt.strip()
gm[73] = gm[73].replace("De Carvalho DD, De Carvalho DD.", "De Carvalho DD.")  # duplicated author in GM list
gm[36] = gm[36].replace("Hollmann TJ, Hollmann TJ,", "Hollmann TJ,")  # duplicated author in GM list
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
 "barriernet": "barriernet: marker-level Barrier/NetImmune scoring for spatial, single-cell, bulk and protein data (Python package, version 0.2.1). https://github.com/drgengtaiye/barriernet.",
 "kalluri2016": "Kalluri R. The biology and function of fibroblasts in cancer. Nat Rev Cancer. 2016;16:582–598.",
 "companion": "[Companion manuscript] Immune gene-expression signatures carry treatment-specific information: a matched-null benchmark of 52 published signatures across checkpoint cohorts, randomized trials and 20 cancer types. Submitted.",
 "mantovani2017": "Mantovani A, Marchesi F, Malesci A, Laghi L, Allavena P. Tumour-associated macrophages as treatment targets in oncology. Nat Rev Clin Oncol. 2017;14:399–416.",
 "qi2022": "Qi J, Sun H, Zhang Y, Wang Z, Xun Z, Li Z, et al. Single-cell and spatial analysis reveal interaction of FAP⁺ fibroblasts and SPP1⁺ macrophages in colorectal cancer. Nat Commun. 2022;13:1742.",
 "ballman2015": gm[9],
 "tauriello2018": "Tauriello DVF, Palomo-Ponce S, Stork D, Berenguer-Llergo A, Badia-Ramentol J, Iglesias M, et al. TGFβ drives immune evasion in genetically reconstituted colon cancer metastasis. Nature. 2018;554:538–543.",
 "zeng2019": "Zeng D, Li M, Zhou R, Zhang J, Sun H, Shi M, et al. Tumor microenvironment characterization in gastric cancer identifies prognostic and immunotherapeutically relevant gene signatures. Cancer Immunol Res. 2019;7:737–750.",
 "beck2013": "Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, Culhane AC, et al. Significance analysis of prognostic signatures. PLoS Comput Biol. 2013;9:e1002875.",
 "yuan2014": "Yuan Y, Van Allen EM, Omberg L, Wagle N, Amin-Mansour A, Sokolov A, et al. Assessing the clinical utility of cancer genomic and proteomic data across tumor types. Nat Biotechnol. 2014;32:644–652.",
 "bagaev2021": "Bagaev A, Kotlov N, Nomie K, Svekolkin V, Gafurov A, Isaeva O, et al. Conserved pan-cancer microenvironment subtypes predict response to immunotherapy. Cancer Cell. 2021;39:845–865.e7.",
 "thorsson2018": "Thorsson V, Gibbs DL, Brown SD, Wolf D, Bortone DS, Ou Yang TH, et al. The immune landscape of cancer. Immunity. 2018;48:812–830.e14.",
 "kieffer2020": "Kieffer Y, Hocine HR, Gentric G, Pelon F, Bernard C, Bourachot B, et al. Single-cell analysis reveals fibroblast clusters linked to immunotherapy resistance in cancer. Cancer Discov. 2020;10:1330–1351.",
 "zhang2020": "Zhang L, Li Z, Skrzypczynska KM, Fang Q, Zhang W, O'Brien SA, et al. Single-cell analyses inform mechanisms of myeloid-targeted therapies in colon cancer. Cell. 2020;181:442–459.e29.",
 "bill2023": "Bill R, Wirapati P, Messemaker M, Roh W, Zitti B, Duval F, et al. CXCL9:SPP1 macrophage polarity identifies a network of cellular programs that control human cancers. Science. 2023;381:515–524.",
 "hartung2001": "Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. Stat Med. 2001;20:3875–3889.",
 "firth1993": "Firth D. Bias reduction of maximum likelihood estimates. Biometrika. 1993;80:27–38.",
 "barbie2009": gm[8], "mcdermott2018": gm[35], "bareche2022": gm[4], "gide2019": gm[18], "hugo2016": gm[19],
 "liu2019": gm[15], "riaz2017": gm[23], "vanallen2015": gm[24], "braun2020": gm[17], "miao2018": gm[22],
 "jung2019": gm[20], "limagne2022": gm[21], "padron2022": gm[36],
 # v16 additions (replication cohorts, I-SPY2, TCGA-CDR, IMC). The GM list pairs several of these
 # author lists with the wrong titles, so they are written out here.
 "nathanson2017": "Nathanson T, Ahuja A, Rubinsteyn A, Aksoy BA, Hellmann MD, Miao D, et al. Somatic mutations and neoepitope homology in melanomas treated with CTLA-4 blockade. Cancer Immunol Res. 2017;5:84–91.",
 "snyder2017": "Snyder A, Nathanson T, Funt SA, Ahuja A, Buros Novik J, Hellmann MD, et al. Contribution of systemic and somatic factors to clinical response and resistance to PD-L1 blockade in urothelial cancer: an exploratory multi-omic analysis. PLoS Med. 2017;14:e1002309.",
 "rose2021": "Rose TL, Weir WH, Mayhew GM, Shibata Y, Eulitt P, Uronis JM, et al. Fibroblast growth factor receptor 3 alterations and response to immune checkpoint inhibition in metastatic urothelial cancer: a real world experience. Br J Cancer. 2021;125:1251–1260.",
 "eddy2020": "Eddy JA, Thorsson V, Lamb AE, Gibbs DL, Heimann C, Yu JX, et al. CRI iAtlas: an interactive portal for immuno-oncology research. F1000Res. 2020;9:1028.",
 "foy2022a": "Foy JP, Karabajakian A, Ortiz-Cuaran S, Boussageon M, Michon L, Bouaoud J, et al. Immunologically active phenotype by gene expression profiling is associated with clinical benefit from PD-1/PD-L1 inhibitors in real-world head and neck and lung cancer patients. Eur J Cancer. 2022;174:287–298.",
 "foy2022b": "Foy JP, Karabajakian A, Ortiz-Cuaran S, Boussageon M, Michon L, Bouaoud J, et al. Datasets for gene expression profiles of head and neck squamous cell carcinoma and lung cancer treated or not by PD1/PD-L1 inhibitors. Data Brief. 2022;44:108556.",
 "vathiotis2022": "Vathiotis IA, Salichos L, Martinez-Morilla S, Gavrielatou N, Aung TN, Shafi S, et al. Baseline gene expression profiling determines long-term benefit to programmed cell death protein 1 axis blockade. NPJ Precis Oncol. 2022;6:92.",
 "prat2017": "Prat A, Navarro A, Paré L, Reguart N, Galván P, Pascual T, et al. Immune-related gene expression profiling after PD-1 blockade in non-small cell lung carcinoma, head and neck squamous cell carcinoma, and melanoma. Cancer Res. 2017;77:3540–3550.",
 "hu2023": "Hu J, Zhang L, Xia H, Yan Y, Zhu X, Sun F, et al. Tumor microenvironment remodeling after neoadjuvant immunotherapy in non-small cell lung cancer revealed by single-cell RNA sequencing. Genome Med. 2023;15:14.",
 "thibaudin2023": "Thibaudin M, Fumet JD, Chibaudel B, Bennouna J, Borg C, Martin-Babau J, et al. First-line durvalumab and tremelimumab with chemotherapy in RAS-mutated metastatic colorectal cancer: a phase 1b/2 trial. Nat Med. 2023;29:2087–2098.",
 "vandenende2021": "van den Ende T, de Clercq NC, van Berge Henegouwen MI, Gisbertz SS, Geijsen ED, Verhoeven RHA, et al. Neoadjuvant chemoradiotherapy combined with atezolizumab for resectable esophageal adenocarcinoma: a single-arm phase II feasibility trial (PERFECT). Clin Cancer Res. 2021;27:3351–3359.",
 "wolf2022": "Wolf DM, Yau C, Wulfkuhle J, Brown-Swigart L, Gallagher RI, Lee PRE, et al. Redefining breast cancer subtypes to guide treatment prioritization and maximize response: predictive biomarkers across 10 cancer therapies. Cancer Cell. 2022;40:609–623.e6.",
 "nanda2020": "Nanda R, Liu MC, Yau C, Shatsky R, Pusztai L, Wallace A, et al. Effect of pembrolizumab plus neoadjuvant chemotherapy on pathologic complete response in women with early-stage breast cancer: an analysis of the ongoing phase 2 adaptively randomized I-SPY2 trial. JAMA Oncol. 2020;6:676–684.",
 "pusztai2021": "Pusztai L, Yau C, Wolf DM, Han HS, Du L, Wallace AM, et al. Durvalumab with olaparib and paclitaxel for high-risk HER2-negative stage II/III breast cancer: results from the adaptively randomized I-SPY2 trial. Cancer Cell. 2021;39:989–998.e5.",
 "liu2018cdr": "Liu J, Lichtenberg T, Hoadley KA, Poisson LM, Lazar AJ, Cherniack AD, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. Cell. 2018;173:400–416.e11.",
 "danenberg2022": "Danenberg E, Bardwell H, Zanotelli VRT, Provenzano E, Chin SF, Rueda OM, et al. Breast tumor microenvironment structures are associated with genomic features and clinical outcome. Nat Genet. 2022;54:660–669.",
 "curtis2012": "Curtis C, Shah SP, Chin SF, Turashvili G, Rueda OM, Dunning MJ, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. Nature. 2012;486:346–352.",
 "cerami2012": "Cerami E, Gao J, Dogrusoz U, Gross BE, Sumer SO, Aksoy BA, et al. The cBio cancer genomics portal: an open platform for exploring multidimensional cancer genomics data. Cancer Discov. 2012;2:401–404.",
 "pdacimc": "Sussman JH, Kim N, Kemp SB, Traum D, Katsuda T, Kahn BM, et al. Multiplexed imaging mass cytometry analysis characterizes the vascular niche in pancreatic cancer. Cancer Res. 2024;84:2364–2376.",
 "anders2022": "Anders CK, Woodcock MG, Van Swearingen AED, Moore DT, Sambade MJ, Laurie S, et al. Evaluating the efficacy of a priming dose of cyclophosphamide prior to pembrolizumab to treat metastatic triple negative breast cancer. J Immunother Cancer. 2022;10:e003427.",
 "ooi2009": "Ooi CH, Ivanova T, Wu J, Lee M, Tan IB, Tao J, et al. Oncogenic pathway combinations predict clinical prognosis in gastric cancer. PLoS Genet. 2009;5:e1000676.",
 "lei2013": "Lei Z, Tan IB, Das K, Deng N, Zouridis H, Pattison S, et al. Identification of molecular subtypes of gastric cancer with different responses to PI3-kinase inhibitors and 5-fluorouracil. Gastroenterology. 2013;145:554–565.",
 "lee2014": "Lee J, Sohn I, Do IG, Kim KM, Park SH, Park JO, et al. Nanostring-based multigene assay to predict recurrence for gastric cancer patients after surgery. PLoS One. 2014;9:e90133.",
 "shiuan2021": "Shiuan E, Reddy A, Dudzinski SO, Lim AR, Sugiura A, Hongo R, et al. Clinical features and multiplatform molecular analysis assist in understanding patient response to anti-PD-1/PD-L1 in renal cell carcinoma. Cancers (Basel). 2021;13:1475.",
}

def fmt(ns):
    ns = sorted(set(ns)); out = []; i = 0
    while i < len(ns):
        j = i
        while j + 1 < len(ns) and ns[j + 1] == ns[j] + 1: j += 1
        out.append(str(ns[i]) if j == i else (f"{ns[i]},{ns[j]}" if j == i + 1 else f"{ns[i]}–{ns[j]}"))
        i = j + 1
    return "^" + ",".join(out) + "^"

src = (HERE / "v17_src.md").read_text()
body, tail = src.split("{{REFERENCES}}")
cite = re.compile(r"\[(@[^\]]+)\]")
order = []
for m in cite.finditer(body):
    for k in re.findall(r"@(\w+)", m.group(1)):
        if k not in order: order.append(k)
missing = [k for k in order if k not in R]; assert not missing, missing
unused = [k for k in R if k not in order]; print("not cited (dropped):", unused)
num = {k: i + 1 for i, k in enumerate(order)}
body = cite.sub(lambda m: fmt([num[k] for k in re.findall(r"@(\w+)", m.group(1))]), body)
out = body + "\n\n".join(f"{num[k]}\\. {R[k]}" for k in order) + tail
(OUT / "MANUSCRIPT_v17.md").write_text(out)
subprocess.run(["pandoc", "MANUSCRIPT_v17.md", "-o", "MANUSCRIPT_v17.docx"], check=True, cwd=OUT)
abstract = out.split("## Abstract")[1].split("**Keywords")[0]
main = out.split("## Introduction")[1].split("## Declarations")[0]
print(f"{len(order)} references | abstract {len(abstract.split())} words | "
      f"Introduction-Conclusion {len(main.split())} words (incl. tables)")
