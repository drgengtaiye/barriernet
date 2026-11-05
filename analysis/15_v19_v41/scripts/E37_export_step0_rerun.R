## Re-run (2026-11-01) of the step 0 of E37_export.R only, because its table results/E37_step0_reproduction.csv was
## later overwritten by E37_test.py, which writes a file of the same name. Same code, output renamed; nothing else written.
## E37 step 1 - export the development data (20 single-arm checkpoint cohorts) for barriernet 0.3.
## Protocol: analysis/E37_prereg.md (hash checked below). Loading and z-scoring exactly as E5_icb_split.R (compendium
## cohorts) and E12_icb_expanded.R (replication cohorts); ICB_Braun is excluded (CheckMate 025 patients).
## Step 0: the NetImmune and Immune log-odds per cohort are re-derived and must equal the E5/E12 values (< 1e-8).
## Writes results/E37_dev_genes.csv (cohort, sample, response, 48 gene z-scores; NA = not measured) and the full log
## matrices of the kidney and breast cohorts for the twin check (results/E37_twin_<cohort>.csv.gz).
suppressMessages({library(matrixStats); library(MultiAssayExperiment); library(org.Hs.eg.db); library(AnnotationDbi)})
A <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"; O <- file.path(A, "results")
h <- strsplit(readLines(file.path(A, "E37_prereg.sha256"))[1], " +")[[1]][1]
h2 <- system(sprintf("sha256sum %s", file.path(A, "E37_prereg.md")), intern = TRUE); h2 <- strsplit(h2, " +")[[1]][1]
if (h != h2) stop("E37 protocol hash mismatch"); cat("E37 protocol sha256 verified:", substr(h, 1, 12), "\n")

D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"))
G48 <- unlist(MOD, use.names = FALSE); stopifnot(length(G48) == 48, !anyDuplicated(G48))
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, keys = G48, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "list"))

## ---- compendium cohorts (E5 load_cohort, unchanged)
E5C <- c("ICB_Gide", "ICB_Hugo", "ICB_Jung", "ICB_Limagne1", "ICB_Liu", "ICB_Mariathasan", "ICB_Miao1", "ICB_Puch",
         "ICB_Riaz", "ICB_Van_Allen")                                  # E5's 11 minus ICB_Braun
load_cohort <- function(nm) {
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
gz_ens <- function(Z) {                                                # 48 x n, symbol rows, NA where not measured
  G <- t(sapply(G48, function(g) { e <- intersect(sym2ens[[g]], rownames(Z)); if (length(e)) Z[e[1], ] else rep(NA_real_, ncol(Z)) }))
  colnames(G) <- colnames(Z); G
}
## ---- replication cohorts (E12 scores_sym z-scoring, unchanged)
CO <- readRDS(file.path(O, "E12_cohorts.rds"))
tab <- read.csv(file.path(O, "E12_cohort_table.csv"), stringsAsFactors = FALSE)
NEW <- tab$id[tab$status == "included"]; stopifnot(length(NEW) == 10)
gz_sym <- function(M) {
  M <- M[rowSds(M) > 0, , drop = FALSE]; Z <- (M - rowMeans(M)) / rowSds(M)
  G <- t(sapply(G48, function(g) if (g %in% rownames(Z)) Z[g, ] else rep(NA_real_, ncol(Z))))
  colnames(G) <- colnames(Z); G
}

DEV <- list()
for (nm in E5C) { d <- load_cohort(nm); DEV[[nm]] <- list(G = gz_ens(d$Z), y = d$y, set = "compendium") }
for (nm in NEW) { d <- CO[[nm]]; DEV[[nm]] <- list(G = gz_sym(d$M), y = d$y, set = "replication") }

## ---- step 0: NetImmune and Immune per cohort as published (E5 / E12 per-cohort tables)
modsc <- function(G) {
  S <- sapply(names(MOD), function(k) { g <- MOD[[k]]; g <- g[rowSums(is.finite(G[g, , drop = FALSE])) == ncol(G)]
    if (length(g) < 3) rep(NA_real_, ncol(G)) else colMeans(G[g, , drop = FALSE]) })
  S <- as.data.frame(S); S$Barrier <- rowMeans(S[, c("CAF_ECM", "TAM_SUP")])
  S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")]); S$NetImmune <- S$Immune - S$Barrier; S
}
logit_eff <- function(y, s) { s <- (s - mean(s)) / sd(s); f <- suppressWarnings(glm(y ~ s, family = binomial()))
  c(b = unname(coef(f)[2]), se = unname(sqrt(diag(vcov(f)))[2])) }
PUB <- read.csv(file.path(O, "E12_icb_per_cohort.csv"), stringsAsFactors = FALSE)
chk <- list()
for (nm in names(DEV)) {
  S <- modsc(DEV[[nm]]$G)
  for (k in c("Immune", "NetImmune", "T_CELL", "CAF_ECM")) {
    p <- PUB[PUB$cohort == nm & PUB$score == k, ]
    if (all(is.na(S[[k]]))) { chk[[length(chk) + 1]] <- data.frame(cohort = nm, score = k, mine = NA, pub = if (nrow(p)) p$logOR else NA, diff = NA); next }
    e <- logit_eff(DEV[[nm]]$y, S[[k]])
    chk[[length(chk) + 1]] <- data.frame(cohort = nm, score = k, mine = e[["b"]], pub = if (nrow(p)) p$logOR else NA,
                                         diff = if (nrow(p)) abs(e[["b"]] - p$logOR) else NA)
  }
}
chk <- do.call(rbind, chk); write.csv(chk, file.path(O, "E37_export_step0_rerun.csv"), row.names = FALSE)
cat("step 0: max |difference| from the published per-cohort log-odds:", signif(max(chk$diff, na.rm = TRUE), 3),
    "over", sum(is.finite(chk$diff)), "values; not scorable as published:", sum(is.na(chk$mine) & is.na(chk$pub)), "\n")
stopifnot(max(chk$diff, na.rm = TRUE) < 1e-8, all(is.na(chk$mine) == is.na(chk$pub)))

