## E45 step 1 - export the 20 single-arm checkpoint cohorts (barriernet 0.3 development set) for the 2,082 genes measured
## in all six randomized comparisons (E40_export.R with a wider gene list). Protocol: analysis/E45_prereg.md (hash checked below). Loading and z-scoring exactly as
## E37_export.R (E5 load_cohort for the compendium cohorts, E12 z-scoring for the replication cohorts).
## Step 0: the 48 barriernet gene z-scores must equal results/E37_dev_genes.csv (< 1e-10).
## Writes results/E45_single_arm_genes.csv.gz (cohort, set, sample, response, one column per gene; NA = not measured).
suppressMessages({library(matrixStats); library(MultiAssayExperiment); library(org.Hs.eg.db); library(AnnotationDbi)})
A <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"; O <- file.path(A, "results")
h <- strsplit(readLines(file.path(A, "E45_prereg.sha256"))[1], " +")[[1]][1]
h2 <- strsplit(system(sprintf("sha256sum %s", file.path(A, "E45_prereg.md")), intern = TRUE), " +")[[1]][1]
if (h != h2) stop("E45 protocol hash mismatch"); cat("E45 protocol sha256 verified:", substr(h, 1, 12), "\n")

G0 <- "/groups/adv2105_gp/yichen/Yi/com/gastric"; D <- file.path(G0, "icb")
## library genes, read as E32_trials.py (aliases as cmp_scores.py)
ALIAS <- c(WARS1 = "WARS", HLA_DPB1 = "HLA-DPB1"); RALIAS <- c(WARS = "WARS1", "HLA-DPB1" = "HLA_DPB1")
FAM <- read.csv(file.path(D, "v5_signature_clusters.csv"), stringsAsFactors = FALSE)
LIBG <- unique(unlist(lapply(FAM$signature, function(s) {
  d <- read.csv(file.path(D, "SignatureSets-main/data-raw/IO", paste0(s, ".csv")), fileEncoding = "latin1", stringsAsFactors = FALSE)
  g <- as.character(d$gene_name); g <- g[!is.na(g) & g != "" & g != "nan"]; ifelse(g %in% names(ALIAS), ALIAS[g], g) })))
G48 <- c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1",
         "SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3",
         "B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9", "STAT1","IRF1","CXCL9","CXCL10","IFNG",
         "CXCL13","LTB","MS4A1","CD79A","CD74","MZB1", "CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7")
P33 <- read.csv(file.path(O, "E33_pooled.csv"), stringsAsFactors = FALSE)
G2082 <- P33$name[P33$set == "six (primary)" & P33$kind == "gene" & P33$all_six %in% c(TRUE, "True")]; stopifnot(length(G2082) == 2082)
GX <- unique(c(G48, G2082)); cat("candidate genes:", length(G2082), "; exported genes:", length(GX), "\n")
cand <- function(g) unique(c(g, if (g %in% names(RALIAS)) RALIAS[[g]] else character(0)))
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, keys = unique(unlist(lapply(GX, cand))), keytype = "SYMBOL",
                                   column = "ENSEMBL", multiVals = "list"))

E5C <- c("ICB_Gide", "ICB_Hugo", "ICB_Jung", "ICB_Limagne1", "ICB_Liu", "ICB_Mariathasan", "ICB_Miao1", "ICB_Puch",
         "ICB_Riaz", "ICB_Van_Allen")
load_cohort <- function(nm) {                                          # E5 load_cohort, unchanged
  x <- suppressWarnings(readRDS(file.path(D, paste0(nm, ".rds"))))
  ex <- experiments(x); an <- names(ex)
  pick <- c("expr_gene_tpm", "expr")[c("expr_gene_tpm", "expr") %in% an][1]
  M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  r <- as.character(cd$response); y <- ifelse(r == "R", 1, ifelse(r == "NR", 0, NA)); names(y) <- rownames(cd)
  keep <- intersect(colnames(M), names(y)[!is.na(y)]); M <- M[, keep, drop = FALSE]; y <- y[keep]
  rownames(M) <- sub("\\..*$", "", rownames(M)); M <- M[!duplicated(rownames(M)), , drop = FALSE]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  s <- rowSds(M); list(Z = (M - rowMeans(M)) / s, M = M, y = y)
}
gz_ens <- function(Z) {
  G <- t(sapply(GX, function(g) {
    e <- intersect(unlist(lapply(cand(g), function(x) sym2ens[[x]])), rownames(Z))
    if (length(e)) Z[e[1], ] else rep(NA_real_, ncol(Z)) }))
  colnames(G) <- colnames(Z); G
}
CO <- readRDS(file.path(O, "E12_cohorts.rds"))
tab <- read.csv(file.path(O, "E12_cohort_table.csv"), stringsAsFactors = FALSE)
NEW <- tab$id[tab$status == "included"]; stopifnot(length(NEW) == 10)
gz_sym <- function(M) {                                                # E12 scores_sym z-scoring, unchanged
  M <- M[rowSds(M) > 0, , drop = FALSE]; Z <- (M - rowMeans(M)) / rowSds(M)
  G <- t(sapply(GX, function(g) { x <- intersect(cand(g), rownames(Z)); if (length(x)) Z[x[1], ] else rep(NA_real_, ncol(Z)) }))
  colnames(G) <- colnames(Z); G
}
DEV <- list()
for (nm in E5C) { d <- load_cohort(nm); DEV[[nm]] <- list(G = gz_ens(d$Z), y = d$y, set = "compendium") }
for (nm in NEW) { d <- CO[[nm]]; DEV[[nm]] <- list(G = gz_sym(d$M), y = d$y, set = "replication") }

rows <- do.call(rbind, lapply(names(DEV), function(nm) {
  G <- DEV[[nm]]$G; data.frame(cohort = nm, set = DEV[[nm]]$set, sample = colnames(G), response = unname(DEV[[nm]]$y[colnames(G)]),
                                t(G), check.names = FALSE) }))
## step 0 against the E37 export (48 genes)
E37 <- read.csv(file.path(O, "E37_dev_genes.csv"), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(E37$cohort, rows$cohort), identical(E37$sample, rows$sample), all(E37$response == rows$response))
a <- as.matrix(E37[, G48]); b <- as.matrix(rows[, G48]); stopifnot(all(is.na(a) == is.na(b)))  # element-wise (identical() also compared row names)
d0 <- max(abs(a - b), na.rm = TRUE); cat("step 0: 48 genes vs E37_dev_genes.csv, max |diff|", signif(d0, 3), "\n"); stopifnot(d0 < 1e-10)
write.csv(rows, gzfile(file.path(O, "E45_single_arm_genes.csv.gz")), row.names = FALSE)
cat("cohorts", length(DEV), "; samples", nrow(rows), "; responders", sum(rows$response), "; genes measured per cohort:\n")
print(sapply(DEV, function(d) sum(rowSums(is.finite(d$G)) == ncol(d$G))))
