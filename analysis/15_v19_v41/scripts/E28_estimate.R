## E28 Q1 (E28_prereg.md). Official ESTIMATE (estimate::estimateScore, platform "affymetrix" as in
## B2b_estimate_fix.R) in the gastric benchmark cohorts and the six randomized comparisons, run three ways:
## as is; genes with any missing value removed; missing values set to the gene's minimum. estimateScore ranks
## each sample with rank(), whose default na.last = TRUE puts missing values at the top of the ranking.
## Usage: Rscript E28_estimate.R <workdir with the E28_dump.py matrices>
suppressMessages(library(estimate))
args <- commandArgs(TRUE); W <- args[1]
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
IMM4 <- c("PTPRC", "CD2", "CD3E", "CD8A"); STR4 <- c("COL1A2", "DCN", "LUM", "PDGFRB")
E <- new.env(); data("SI_geneset", package = "estimate", envir = E); sig <- E$SI_geneset
gi <- as.character(unlist(sig["ImmuneSignature", -1])); gs <- as.character(unlist(sig["StromalSignature", -1]))
stopifnot(length(gi) == 141, length(gs) == 141)

official <- function(m, tag) {
  fin <- file.path(W, paste0(tag, "_in.txt")); fc <- file.path(W, paste0(tag, "_c.gct")); fs <- file.path(W, paste0(tag, "_e.gct"))
  write.table(cbind(GeneSymbol = rownames(m), as.data.frame(m, check.names = FALSE)), fin, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(capture.output(filterCommonGenes(fin, fc, id = "GeneSymbol")))
  invisible(capture.output(estimateScore(fc, fs, platform = "affymetrix")))
  g <- read.delim(fs, skip = 2, row.names = 1, check.names = FALSE)
  g <- g[, setdiff(colnames(g), "Description"), drop = FALSE]
  stopifnot(identical(make.names(colnames(g)), make.names(colnames(m))))
  data.frame(imm = as.numeric(g["ImmuneScore", ]), str = as.numeric(g["StromalScore", ]))
}
zmean <- function(m, genes) { g <- intersect(genes, rownames(m)); z <- t(scale(t(m[g, , drop = FALSE]))); colMeans(z, na.rm = TRUE) }
canon <- function(m, genes) zmean(m, intersect(genes, rownames(m)))

run <- function(m, tag, cohort) {
  m <- m[!duplicated(rownames(m)), , drop = FALSE]
  na_g <- rowSums(is.na(m)) > 0
  vars <- list(as_is = m, drop_na = m[!na_g, , drop = FALSE],
               na_min = { x <- m; mn <- apply(x, 1, min, na.rm = TRUE); i <- which(is.na(x), arr.ind = TRUE); x[i] <- mn[i[, 1]]; x })
  ci <- canon(m, IMM4); cs <- canon(m, STR4); zi <- zmean(m, gi); zs <- zmean(m, gs)
  out <- list(); sc <- data.frame(sample = colnames(m))
  for (v in names(vars)) {
    e <- official(vars[[v]], paste0(tag, "_", v))
    sc[[paste0("imm_", v)]] <- e$imm; sc[[paste0("str_", v)]] <- e$str
    out[[v]] <- data.frame(cohort = cohort, variant = v,
      n_samples = ncol(m), genes_with_missing = sum(na_g), missing_values = sum(is.na(m)),
      immune_genes_with_missing = sum(na_g[intersect(gi, rownames(m))]), stromal_genes_with_missing = sum(na_g[intersect(gs, rownames(m))]),
      immune_genes_present = length(intersect(gi, rownames(vars[[v]]))), stromal_genes_present = length(intersect(gs, rownames(vars[[v]]))),
      r_official_immune_canonical = cor(e$imm, ci, use = "pair"), r_official_stromal_canonical = cor(e$str, cs, use = "pair"),
      r_official_vs_zmean_immune = cor(e$imm, zi, use = "pair"), r_official_vs_zmean_stromal = cor(e$str, zs, use = "pair"),
      r_official_vs_zmean_contrast = cor(scale(e$imm)[, 1] - scale(e$str)[, 1], zi - zs, use = "pair"))
  }
  sc$zmean_imm <- zi; sc$zmean_str <- zs
  write.csv(sc, file.path(W, paste0(tag, "_official_scores.csv")), row.names = FALSE)
  out <- do.call(rbind, out)
  out$r_zmean_immune_canonical <- cor(zi, ci, use = "pair"); out$r_zmean_stromal_canonical <- cor(zs, cs, use = "pair")
  out
}

res <- list()
for (t in c("TCGA", "GSE84437", "GSE62254")) {
  e <- readRDS(file.path(R, paste0("B1_", t, ".rds")))$expr
  res[[t]] <- run(e, t, t)
}
for (k in c("IMM150", "ISPY2_pembro", "ISPY2_durva", "CM025", "JR101", "NeoTRIP")) {
  m <- as.matrix(read.delim(file.path(W, paste0(k, ".tsv")), row.names = 1, check.names = FALSE))
  res[[k]] <- run(m, k, k)
}
R1 <- do.call(rbind, res); rownames(R1) <- NULL
write.csv(R1, file.path(O, "E28_Q1_estimate_audit.csv"), row.names = FALSE)
options(width = 250)
print(format(R1[, c("cohort", "variant", "genes_with_missing", "immune_genes_with_missing", "immune_genes_present",
                    "r_official_immune_canonical", "r_zmean_immune_canonical", "r_official_stromal_canonical",
                    "r_official_vs_zmean_contrast")], digits = 3), row.names = FALSE)
