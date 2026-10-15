## E21 step 1 - check without outcomes (pre-registered, E17_E22_prereg.md, Deviations 2026-10-05).
## For the CheckMate 009/010/025 nivolumab patients present both in the Braun et al. 2020
## supplementary matrix (S4A; log2 UQ-TPM, ComBat, +21) and in the compendium object used for
## the first 11 checkpoint cohorts (icb/ICB_Braun.rds, assay "expr"), score both matrices and
## correlate the scores per patient (Spearman). No outcome value is read or used.
##
## Compendium scoring is copied verbatim from E5_icb_split.R (load_cohort / scores_of: Ensembl
## IDs without version, first Ensembl match per symbol via org.Hs.eg.db, non-finite and
## zero-SD genes dropped, row-z, module mean of >= 3 genes). S4A scoring follows the manuscript
## by symbol (duplicate symbols collapsed by mean, row-z, >= 3 genes).
##
## Comparisons (z always computed over the patients named, in each source separately):
##   A  all overlapping nivolumab patients (n = 181), z within that set in both sources  [primary]
##   B  as used in the 11-cohort analysis: compendium z over its response-labelled patients
##      (response R/NR present; the label's value is not used), S4A z over the same patients
##   C  CM-025 nivolumab patients: S4A z over the whole CM-025 randomised set (both arms, as in
##      the E21 interaction models) vs the compendium as used in B (patients in both)
## Transform diagnostics (scale only): per-gene SD ratio, per-gene Pearson and Spearman between
## the two sources across patients, curvature (quadratic coefficient of compendium-z on S4A-z),
## and the across-gene agreement of levels within each sample.
suppressMessages({library(MultiAssayExperiment); library(SummarizedExperiment); library(matrixStats)
                  library(org.Hs.eg.db)})
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/icb"
V <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/validation/Braun2020_CheckMate"
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
SC <- c("CAF_ECM", "TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1", "IFN", "APM", "TLS_B", "T_CELL",
        "Barrier", "Immune", "NetImmune")
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, keys = unique(unlist(MOD)), keytype = "SYMBOL",
                                   column = "ENSEMBL", multiVals = "list"))
sym2ens <- c(sym2ens, suppressMessages(mapIds(org.Hs.eg.db, keys = "CXCL9", keytype = "SYMBOL",
                                              column = "ENSEMBL", multiVals = "list")))

## ---- compendium: E5_icb_split.R load_cohort(), with the sample set made selectable
x <- suppressWarnings(readRDS(file.path(D, "ICB_Braun.rds")))
ex <- experiments(x); M0 <- as.matrix(assay(ex[["expr"]])); cd <- as.data.frame(colData(x))
cd <- cd[match(colnames(M0), cd$patientid), ]
comp_Z <- function(keep) {
  M <- M0[, keep, drop = FALSE]
  rownames(M) <- sub("\\..*$", "", rownames(M)); M <- M[!duplicated(rownames(M)), , drop = FALSE]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  s <- rowSds(M); (M - rowMeans(M)) / s
}
## E5_icb_split.R scores_of(), unchanged
scores_of <- function(Z) {
  used <- list()
  m <- sapply(names(MOD), function(k) {
    ids <- unlist(lapply(MOD[[k]], function(g) { e <- intersect(sym2ens[[g]], rownames(Z)); if (length(e)) e[1] else NULL }))
    used[[k]] <<- length(ids)
    if (length(ids) < 3) return(rep(NA_real_, ncol(Z)))
    colMeans(Z[ids, , drop = FALSE])
  })
  S <- as.data.frame(m)
  S$Barrier <- rowMeans(S[, c("CAF_ECM", "TAM_SUP")])
  S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")])
  S$NetImmune <- S$Immune - S$Barrier
  gi <- function(g) { e <- intersect(sym2ens[[g]], rownames(Z)); if (length(e)) e[1] else NA }
  S$CXCL9_SPP1 <- if (!is.na(gi("CXCL9")) && !is.na(gi("SPP1"))) Z[gi("CXCL9"), ] - Z[gi("SPP1"), ] else NA_real_
  list(S = S, used = unlist(used))
}

## ---- S4A (symbol level; duplicate symbols collapsed by mean; already log2, not re-logged)
S4 <- read.delim(gzfile(file.path(V, "S4A_RNA_expression.tsv.gz")), check.names = FALSE)
g4 <- as.character(S4[[1]]); S4 <- as.matrix(S4[, -1]); S4 <- rowsum(S4, g4); S4 <- S4 / as.vector(table(g4)[rownames(S4)])
stopifnot(!anyDuplicated(rownames(S4)))
cat(sprintf("S4A: %d symbols x %d samples; range %.2f-%.2f, median %.2f (log2 already)\n",
            nrow(S4), ncol(S4), min(S4), max(S4), median(S4)))
s4_scores <- function(ids) {
  M <- S4[, ids, drop = FALSE]; sdv <- rowSds(M); sdv[sdv == 0] <- 1
  Z <- (M - rowMeans(M)) / sdv
  used <- c(); S <- sapply(names(MOD), function(k) { g <- intersect(MOD[[k]], rownames(Z)); used[k] <<- length(g)
    if (length(g) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[g, , drop = FALSE]) })
  S <- as.data.frame(S); rownames(S) <- ids
  S$Barrier <- rowMeans(S[, c("CAF_ECM", "TAM_SUP")]); S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")])
  S$NetImmune <- S$Immune - S$Barrier; S$CXCL9_SPP1 <- Z["CXCL9", ] - Z["SPP1", ]
  list(S = S, used = used)
}

## ---- patient matching: compendium RNA_ID = S4A column; cross-checked against S1 SUBJID
S1 <- read.delim(file.path(V, "S1_clinical.tsv"), check.names = FALSE)
ov <- !is.na(cd$RNA_ID) & cd$RNA_ID %in% colnames(S4)
cat(sprintf("compendium RNA samples: %d; matched to S4A by RNA_ID: %d\n", ncol(M0), sum(ov)))
subj <- S1$SUBJID[match(cd$RNA_ID, S1$RNA_ID)]
norm_id <- function(z) gsub("[-_]", "", z)
cat(sprintf("  of these, S1 SUBJID agrees with compendium patientid (ignoring - and _): %d; arm in S1: %s\n",
            sum(norm_id(subj) == norm_id(cd$patientid), na.rm = TRUE),
            paste(names(table(S1$Arm[match(cd$RNA_ID[ov], S1$RNA_ID)])), table(S1$Arm[match(cd$RNA_ID[ov], S1$RNA_ID)]), collapse = " ")))
cm025_all <- S1$RNA_ID[S1$Cohort == "CM-025" & !is.na(S1$RNA_ID) & S1$RNA_ID %in% colnames(S4)]
labelled <- as.character(cd$response) %in% c("R", "NR")

compare <- function(lab, comp_keep, s4_ids, on) {
  cs <- scores_of(comp_Z(cd$patientid[comp_keep]))
  rownames(cs$S) <- cd$RNA_ID[comp_keep]
  ss <- s4_scores(s4_ids)
  do.call(rbind, lapply(SC, function(k) {
    a <- cs$S[on, k]; b <- ss$S[on, k]
    data.frame(comparison = lab, score = k, n = length(on),
               genes_compendium = if (k %in% names(cs$used)) cs$used[[k]] else NA,
               genes_S4A = if (k %in% names(ss$used)) ss$used[[k]] else NA,
               spearman = cor(a, b, method = "spearman"), pearson = cor(a, b))
  }))
}
A_ids <- cd$RNA_ID[ov]
B_keep <- ov & labelled; B_ids <- cd$RNA_ID[B_keep]
C_on <- intersect(cd$RNA_ID[B_keep], cm025_all)
RES <- rbind(
  compare("A: all overlapping nivolumab patients, z within that set", ov, A_ids, A_ids),
  compare("B: as used in the 11-cohort analysis (response-labelled), S4A z over the same patients", B_keep, B_ids, B_ids),
  compare("C: CM-025 nivolumab; S4A z over the CM-025 randomised set (both arms) vs compendium as used", B_keep, cm025_all, C_on))
RES$flag_below_0.95 <- RES$spearman < 0.95
write.csv(RES, file.path(O, "E21_braun_consistency.csv"), row.names = FALSE)
options(width = 200)
cat(sprintf("\nA: n = %d | B: n = %d | C: n = %d (CM-025 set for S4A z: %d)\n", length(A_ids), length(B_ids), length(C_on), length(cm025_all)))
for (cmp in unique(RES$comparison)) {
  cat("\n==", cmp, "\n"); r <- RES[RES$comparison == cmp, -1]; r$spearman <- round(r$spearman, 3); r$pearson <- round(r$pearson, 3)
  print(r, row.names = FALSE)
}

## ---- transform diagnostics (scale only, all 181 overlapping patients)
Mc <- M0[, ov]; colnames(Mc) <- cd$RNA_ID[ov]
rd <- as.data.frame(rowData(ex[["expr"]])); sym <- rd$gene_name
keep <- !is.na(sym) & !(sym %in% sym[duplicated(sym)]) & sym %in% rownames(S4)
Mc <- Mc[keep, ]; rownames(Mc) <- sym[keep]; Ms <- S4[rownames(Mc), A_ids]
ok <- rowSds(Mc) > 0 & rowSds(Ms) > 0; Mc <- Mc[ok, ]; Ms <- Ms[ok, ]
curv <- function(a, b) { za <- (a - mean(a)) / sd(a); zb <- (b - mean(b)) / sd(b); unname(coef(lm(za ~ zb + I(zb^2)))[3]) }
modg <- unique(unlist(MOD)); modg <- modg[modg %in% rownames(Mc)]
TG <- do.call(rbind, lapply(modg, function(g) data.frame(gene = g, sd_compendium = sd(Mc[g, ]), sd_S4A = sd(Ms[g, ]),
  mean_compendium = mean(Mc[g, ]), mean_S4A = mean(Ms[g, ]),
  pearson = cor(Mc[g, ], Ms[g, ]), spearman = cor(Mc[g, ], Ms[g, ], method = "spearman"), curvature = curv(Mc[g, ], Ms[g, ]))))
TG$sd_ratio_S4A_over_compendium <- TG$sd_S4A / TG$sd_compendium
set.seed(2026); smp <- sample(rownames(Mc), 3000)
allg <- data.frame(sd_ratio = rowSds(Ms) / rowSds(Mc),
                   pearson = sapply(seq_len(nrow(Mc)), function(i) cor(Mc[i, ], Ms[i, ])))
allg$spearman <- NA; allg[smp, "spearman"] <- sapply(smp, function(g) cor(Mc[g, ], Ms[g, ], method = "spearman"))
lev <- sapply(seq_len(ncol(Mc)), function(j) cor(Mc[, j], Ms[, j], method = "spearman"))
write.csv(TG, file.path(O, "E21_braun_transform_genes.csv"), row.names = FALSE)
SUM <- data.frame(
  quantity = c("genes compared (unique symbols, non-constant in both)", "median per-gene SD ratio, S4A / compendium (all genes)",
               "median per-gene SD ratio, module genes", "median per-gene Pearson across patients (all genes)",
               "median per-gene Spearman across patients (3,000 random genes)", "median per-gene Pearson, module genes",
               "median per-gene Spearman, module genes", "range per-gene Pearson, module genes",
               "median |curvature| (quadratic coef., z units), module genes", "max |curvature|, module genes",
               "median within-sample Spearman of gene levels across genes (compendium vs S4A)",
               "compendium ACTB mean (S4A mean)", "compendium GAPDH mean (S4A mean)"),
  value = c(nrow(Mc), median(allg$sd_ratio), median(TG$sd_ratio_S4A_over_compendium), median(allg$pearson),
            median(allg$spearman, na.rm = TRUE), median(TG$pearson), median(TG$spearman),
            sprintf("%.3f to %.3f", min(TG$pearson), max(TG$pearson)),
            median(abs(TG$curvature)), max(abs(TG$curvature)), median(lev),
            sprintf("%.2f (%.2f)", mean(Mc["ACTB", ]), mean(Ms["ACTB", ])), sprintf("%.2f (%.2f)", mean(Mc["GAPDH", ]), mean(Ms["GAPDH", ]))))
write.csv(SUM, file.path(O, "E21_braun_transform_summary.csv"), row.names = FALSE)
cat("\n== transform diagnostics (181 overlapping patients)\n"); print(SUM, row.names = FALSE)
cat("\nmodule genes, lowest per-gene Pearson:\n"); print(head(TG[order(TG$pearson), c("gene", "pearson", "spearman", "curvature", "sd_ratio_S4A_over_compendium")], 8), row.names = FALSE, digits = 3)

## ---- is the disagreement a shared per-sample component? (comparison A; scale only)
## If the two sources differ mainly by a sample-level offset common to all genes (normalisation /
## batch handling), the S4A - compendium discrepancy of Immune and of Barrier move together and
## cancel in NetImmune = Immune - Barrier.
csA <- scores_of(comp_Z(cd$patientid[ov]))$S; ssA <- s4_scores(A_ids)$S
zs <- function(v) (v - mean(v)) / sd(v)
dI <- zs(ssA$Immune) - zs(csA$Immune); dB <- zs(ssA$Barrier) - zs(csA$Barrier)
ZcA <- comp_Z(cd$patientid[ov]); M4 <- S4[, A_ids]; sd4 <- rowSds(M4); sd4[sd4 == 0] <- 1; Z4 <- (M4 - rowMeans(M4)) / sd4
gm_diff <- colMeans(Z4[intersect(rownames(Z4), unique(sym[keep])), ]) - colMeans(ZcA)
SH <- data.frame(quantity = c("cor(Immune discrepancy, Barrier discrepancy)",
                              "cor(Immune discrepancy, all-gene mean-z discrepancy)",
                              "cor(Barrier discrepancy, all-gene mean-z discrepancy)"),
                 value = c(cor(dI, dB), cor(dI, gm_diff), cor(dB, gm_diff)))
write.csv(SH, file.path(O, "E21_braun_shared_component.csv"), row.names = FALSE)
cat("\n== shared per-sample component (comparison A)\n"); print(SH, row.names = FALSE, digits = 3)
