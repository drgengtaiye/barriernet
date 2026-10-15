## E12a - audit of candidate checkpoint cohorts for the expanded cross-cancer analysis (E12).
## Nothing here relates a score to response. It records, for every candidate: cancer type,
## treatment, platform, expression scale, labelled n / responders, module gene coverage, and
## sample-level duplicate checks against the 11 E5 cohorts. Writes results/E12_cohort_audit.csv,
## results/E12_duplicate_check.csv and a cached list of standardised cohorts
## (results/E12_cohorts.rds) that E12_icb_expanded.R reads.
suppressMessages({library(MultiAssayExperiment); library(SummarizedExperiment); library(matrixStats)
                  library(org.Hs.eg.db)})
ICB <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
VAL <- "/groups/adv2105_gp/yichen/Yi/com/gastric/validation/cohorts"
IAT <- "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan/ici"
KIM <- "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"

MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"),
  C1Q_arm = c("C1QA","C1QB","C1QC"),
  SPP1_arm = c("SPP1","TREM2","MARCO","APOE","LGALS3"))

## symbol-level log-scale matrix (genes x samples) and labels; duplicated symbols -> highest mean
collapse <- function(M, sym) {
  ok <- !is.na(sym) & sym != ""; M <- M[ok, , drop = FALSE]; sym <- sym[ok]
  o <- order(-rowMeans(M, na.rm = TRUE)); M <- M[o, , drop = FALSE]; sym <- sym[o]
  M <- M[!duplicated(sym), , drop = FALSE]; rownames(M) <- sym[!duplicated(sym)]; M
}
ens2sym <- function(ids) {
  e <- sub("\\..*$", "", ids)
  unname(suppressMessages(mapIds(org.Hs.eg.db, keys = e, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")))
}
from_mae <- function(nm) {
  x <- suppressWarnings(readRDS(file.path(ICB, paste0(nm, ".rds")))); ex <- experiments(x)
  pick <- c("expr_gene_tpm", "expr")[c("expr_gene_tpm", "expr") %in% names(ex)][1]
  se <- ex[[pick]]; M <- as.matrix(assay(se)); rd <- as.data.frame(rowData(se))
  sym <- if ("gene_name" %in% colnames(rd)) rd$gene_name else ens2sym(rownames(M))
  miss <- is.na(sym); if (any(miss)) sym[miss] <- ens2sym(rownames(M)[miss])
  cd <- as.data.frame(colData(x)); r <- as.character(cd$response)
  y <- setNames(ifelse(r == "R", 1, ifelse(r == "NR", 0, NA)), rownames(cd))
  list(M = collapse(M, sym), y = y, assay = pick, cd = cd)
}
from_tsv <- function(nm) {
  M <- read.delim(gzfile(file.path(VAL, paste0(nm, ".expr.tsv.gz"))), check.names = FALSE, row.names = 1)
  cl <- read.delim(file.path(VAL, paste0(nm, ".clin.tsv")), colClasses = c("character", "numeric"))
  list(M = collapse(as.matrix(M), rownames(M)), y = setNames(cl$y, cl$sample), assay = "expr.tsv")
}
from_iatlas <- function(nm) {
  M <- read.csv(file.path(IAT, paste0(nm, ".csv")), check.names = FALSE, row.names = 1)
  cl <- read.csv(file.path(IAT, paste0(nm, "_clinical.csv")), check.names = FALSE)
  if ("SAMPLE_TREATMENT" %in% names(cl)) cl <- cl[cl$SAMPLE_TREATMENT == "Pre", ]
  r <- if ("RESPONSE" %in% names(cl)) cl$RESPONSE else rep(NA, nrow(cl)); y <- ifelse(r %in% c("Complete Response", "Partial Response"), 1,
                                ifelse(r %in% c("Stable Disease", "Progressive Disease"), 0, NA))
  list(M = collapse(as.matrix(M), rownames(M)), y = setNames(y, cl$sampleId), assay = "iAtlas csv")
}
from_kim <- function() {
  b <- read.delim(file.path(KIM, "benchmark_scores.tsv"), row.names = 1)
  tpm <- read.delim(file.path(KIM, "mGC_Pembrolizumab_Kim2018.TPM"), row.names = 1, check.names = FALSE)[, rownames(b)]
  ## the deposited ".TPM" file is already log2(TPM + 1) (column sums ~5.5e4, max 16.8), so it is used
  ## as is; E7a applies a second log2(x + 1), which E12 reports as a sensitivity
  list(M = collapse(as.matrix(tpm), rownames(tpm)), y = setNames(b$Responder, rownames(b)),
       assay = "log2(TPM+1) as deposited (manuscript PRJEB25780 labels)")
}

## candidate list: source, cancer, treatment, platform, group (whole transcriptome / panel), role
CAND <- read.csv(text = '
id,src,cancer,treatment,platform,tech,role
ICB_VanDenEnde,mae,oesophageal adenocarcinoma,neoadjuvant atezolizumab + chemoradiotherapy (PERFECT),RNA-seq (GEO GSE165252 via compendium),WT,new
ICB_Roh,mae,melanoma,anti-CTLA-4 then anti-PD-1 (response to anti-CTLA-4),targeted (721-gene extract),panel,new
ICB_Nathanson,mae,melanoma,ipilimumab (anti-CTLA-4),RNA-seq,WT,new
ICB_Snyder,mae,urothelial,atezolizumab (anti-PD-L1),RNA-seq,WT,new
ICB_Hwang,mae,NSCLC,anti-PD-1,targeted (384-gene extract),panel,new
ICB_Jerby_Arnon,mae,melanoma,pembrolizumab (anti-PD-1),targeted (587-gene extract),panel,new
ICB_Limagne2,mae,NSCLC,anti-PD-1/PD-L1,RNA-seq,WT,new
ICB_Shiuan,mae,renal cell,anti-PD-1,RNA-seq,WT,new
ICB_Padron,mae,pancreatic,chemo + anti-PD-1 +/- CD40 agonist (PRINCE),RNA-seq,WT,excluded
GSE176307,tsv,urothelial,anti-PD-1/PD-L1,RNA-seq (Salmon TPM),WT,new
GSE161537,tsv,NSCLC,anti-PD-1/PD-L1 (2nd line+),HTG EdgeSeq OBP (~2.5k genes),panel,new
GSE159067,tsv,HNSCC,anti-PD-1/PD-L1,HTG EdgeSeq OBP (~2.5k genes),panel,new
GSE303300,tsv,cervical SCC,checkpoint window before chemoradiation,HTG EdgeSeq OBP (~2.5k genes),panel,new
GSE215868,tsv,melanoma,anti-PD-1 axis,NanoString IO360 (750 genes),panel,new
GSE93157,tsv,NSCLC / HNSCC / melanoma,anti-PD-1,NanoString PanCancer Immune (770 genes),panel,new
GSE207422,tsv,NSCLC,neoadjuvant anti-PD-1 + chemotherapy (MPR),RNA-seq,WT,chemoICB
GSE235919,tsv,colorectal (MSS),durvalumab + tremelimumab + FOLFOX (MEDITREME),RNA-seq,WT,chemoICB
brca_iatlas_anders_2022,iatlas,breast (metastatic),ICB + cyclophosphamide (iAtlas NEOICI_RX),RNA-seq (iAtlas 2475-gene extract),WT,new
gbm_iatlas_prins_2019,iatlas,glioblastoma,neoadjuvant/adjuvant pembrolizumab,RNA-seq (iAtlas extract),WT,new
mel_dfci_2019,iatlas,melanoma,anti-CTLA-4 / anti-PD-1,RNA-seq (cBioPortal extract),WT,new
PRJEB25780,kim,gastric,pembrolizumab (Kim 2018),RNA-seq,WT,gastric_discovery
', stringsAsFactors = FALSE, strip.white = TRUE)

load_any <- function(r) switch(r$src, mae = from_mae(r$id), tsv = from_tsv(r$id), iatlas = from_iatlas(r$id), kim = from_kim())

cov_of <- function(M) sapply(MOD, function(g) sum(g %in% rownames(M)))
CO <- list(); rows <- list()
for (i in seq_len(nrow(CAND))) {
  r <- CAND[i, ]; d <- load_any(r)
  y <- d$y[!is.na(d$y)]; keep <- intersect(colnames(d$M), names(y))
  M <- d$M[, keep, drop = FALSE]; y <- y[keep]
  M <- M[rowSums(is.finite(M)) == ncol(M), , drop = FALSE]; M <- M[rowSds(M) > 0, , drop = FALSE]
  cv <- cov_of(M); rng <- range(M)
  rows[[i]] <- data.frame(r[, c("id", "cancer", "treatment", "platform", "tech", "role")], genes = nrow(M),
    n_labelled = length(y), n_resp = sum(y), n_nonresp = sum(1 - y), min = round(rng[1], 2), max = round(rng[2], 2),
    t(cv), modules_ok = paste(names(cv)[cv >= 3], collapse = " "), stringsAsFactors = FALSE)
  CO[[r$id]] <- list(M = M, y = y, meta = r)
  cat(sprintf("%-26s genes %6d  n %3d  R %3d  range %6.1f..%5.1f  cov %s\n", r$id, nrow(M), length(y), sum(y),
              rng[1], rng[2], paste(cv, collapse = "/")))
}
aud <- do.call(rbind, rows)

## ---- duplicate checks against the 11 E5 cohorts (and among candidates)
## (a) sample/patient identifiers; (b) expression: for each candidate sample, the best Spearman
## correlation with any sample of an earlier cohort of the same cancer type, over shared genes,
## against the best correlation with other samples of its own cohort (a re-profiled tumour would
## correlate with its twin as closely as replicate libraries, ~0.97+).
E5 <- c("ICB_Braun","ICB_Gide","ICB_Hugo","ICB_Jung","ICB_Limagne1","ICB_Liu","ICB_Mariathasan","ICB_Miao1",
        "ICB_Puch","ICB_Riaz","ICB_Van_Allen")
OLD <- lapply(setNames(E5, E5), function(nm) { d <- from_mae(nm); list(M = d$M, cd = d$cd) })
old_cancer <- c(ICB_Braun = "renal", ICB_Gide = "melanoma", ICB_Hugo = "melanoma", ICB_Jung = "NSCLC",
                ICB_Limagne1 = "NSCLC", ICB_Liu = "melanoma", ICB_Mariathasan = "urothelial", ICB_Miao1 = "renal",
                ICB_Puch = "melanoma", ICB_Riaz = "melanoma", ICB_Van_Allen = "melanoma")
match_old <- list(ICB_Nathanson = c("ICB_Van_Allen", "ICB_Liu", "ICB_Riaz", "ICB_Gide", "ICB_Hugo", "ICB_Puch"),
                  ICB_Jerby_Arnon = c("ICB_Gide", "ICB_Liu", "ICB_Riaz", "ICB_Hugo", "ICB_Puch", "ICB_Van_Allen"),
                  ICB_Roh = c("ICB_Gide", "ICB_Liu", "ICB_Van_Allen", "ICB_Riaz"),
                  GSE215868 = c("ICB_Gide", "ICB_Liu", "ICB_Riaz", "ICB_Hugo", "ICB_Puch"),
                  ICB_Snyder = c("ICB_Mariathasan"), GSE176307 = c("ICB_Mariathasan"),
                  ICB_Limagne2 = c("ICB_Limagne1", "ICB_Jung"), GSE161537 = c("ICB_Jung", "ICB_Limagne1"),
                  ICB_Shiuan = c("ICB_Braun", "ICB_Miao1"), GSE93157 = c("ICB_Jung", "ICB_Limagne1", "ICB_Gide", "ICB_Hugo"),
                  brca_iatlas_anders_2022 = character(0))
## candidate-vs-candidate pairs of one cancer type
OLD$GSE215868 <- list(M = CO$GSE215868$M, cd = data.frame(patientid = character(0)))
OLD$GSE176307 <- list(M = CO$GSE176307$M, cd = data.frame(patientid = character(0)))
OLD$GSE161537 <- list(M = CO$GSE161537$M, cd = data.frame(patientid = character(0)))
match_old$GSE93157 <- c(match_old$GSE93157, "GSE215868", "GSE161537")
match_old$ICB_Snyder <- c(match_old$ICB_Snyder, "GSE176307")
match_old$ICB_Limagne2 <- c(match_old$ICB_Limagne2, "GSE161537")
## best and second-best cross-cohort correlation per candidate sample: a re-profiled tumour has one
## partner that stands clear of the rest (gap best - second >> the ~0.005 seen between unrelated
## cohorts of one cancer type)
cross_cor <- function(A, B, same = FALSE) {
  g <- intersect(rownames(A), rownames(B)); if (length(g) < 200) return(NULL)
  C <- cor(apply(A[g, , drop = FALSE], 2, rank), apply(B[g, , drop = FALSE], 2, rank))
  if (same) diag(C) <- NA
  t(apply(C, 1, function(v) { v <- unname(sort(v, decreasing = TRUE)); c(best = v[1], gap = v[1] - v[2]) }))
}
normid <- function(v) sub("^(patient|pt|p)", "", sub("_t(_[a-z])?$", "", tolower(as.character(v))))
dup <- list()
for (nm in names(match_old)) {
  A <- CO[[nm]]$M
  own <- cross_cor(A, A, same = TRUE)
  for (o in match_old[[nm]]) {
    B <- OLD[[o]]$M; g <- length(intersect(rownames(A), rownames(B)))
    cc <- cross_cor(A, B)
    idA <- normid(colnames(A)); idB <- normid(c(colnames(B), as.character(OLD[[o]]$cd$patientid)))
    dup[[length(dup) + 1]] <- data.frame(candidate = nm, against = o, shared_genes = g,
      shared_ids = length(intersect(idA, idB)),
      median_best_r_other = if (is.null(cc)) NA else round(median(cc[, "best"]), 3),
      median_gap_other = if (is.null(cc)) NA else round(median(cc[, "gap"]), 3),
      n_gap_gt_0.03 = if (is.null(cc)) NA else sum(cc[, "gap"] > 0.03),
      median_best_r_own = round(median(own[, "best"]), 3), median_gap_own = round(median(own[, "gap"]), 3))
  }
}
## mel_dfci_2019 sample names (Sample1...) against the iAtlas Liu extract (Liu_Sample1...)
liu_iat <- colnames(read.csv(file.path(IAT, "mel_iatlas_liu_2019.csv"), nrows = 1, check.names = FALSE))[-1]
dfci <- colnames(read.csv(file.path(IAT, "mel_dfci_2019.csv"), nrows = 1, check.names = FALSE))[-1]
dup[[length(dup) + 1]] <- data.frame(candidate = "mel_dfci_2019", against = "mel_iatlas_liu_2019 (= ICB_Liu) ids",
  shared_genes = NA, shared_ids = length(intersect(dfci, sub("^Liu_", "", liu_iat))), median_best_r_other = NA,
  median_gap_other = NA, n_gap_gt_0.03 = NA, median_best_r_own = NA, median_gap_own = NA)
dupt <- do.call(rbind, dup); print(dupt)
write.csv(aud, file.path(O, "E12_cohort_audit.csv"), row.names = FALSE)
write.csv(dupt, file.path(O, "E12_duplicate_check.csv"), row.names = FALSE)
saveRDS(CO, file.path(O, "E12_cohorts.rds"))
