## E21 step 1, POST HOC follow-up (not pre-registered; run because E21_braun_check.R found module-level
## score correlations < 0.95 between the compendium Braun matrix and the authors' S4A matrix).
## How much would the 11-cohort Braun response estimate (E5_icb_split.R) change if Braun were scored
## from S4A instead of the compendium object? Same 142 response-labelled nivolumab patients, same
## response labels (compendium colData), same model (logistic, log-odds of response per SD) and same
## DerSimonian-Laird pooling of the 11 cohorts; only the Braun matrix changes. Uses no outcome that was
## not already analysed in E5 (no new trial outcome). Earlier result files are not modified.
suppressMessages({library(MultiAssayExperiment); library(SummarizedExperiment); library(matrixStats); library(metafor)})
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
SC <- c("CAF_ECM", "TAM_SUP", "IFN", "APM", "TLS_B", "T_CELL", "Barrier", "Immune", "NetImmune",
        "C1Q_arm", "SPP1_arm", "CXCL9_SPP1")
x <- suppressWarnings(readRDS(file.path(D, "ICB_Braun.rds")))
M0 <- as.matrix(assay(experiments(x)[["expr"]])); cd <- as.data.frame(colData(x))
cd <- cd[match(colnames(M0), cd$patientid), ]
r <- as.character(cd$response); keep <- r %in% c("R", "NR")
y <- ifelse(r[keep] == "R", 1, 0); ids <- cd$RNA_ID[keep]
S4 <- read.delim(gzfile(file.path(V, "S4A_RNA_expression.tsv.gz")), check.names = FALSE)
g4 <- as.character(S4[[1]]); S4 <- as.matrix(S4[, -1]); S4 <- rowsum(S4, g4); S4 <- S4 / as.vector(table(g4)[rownames(S4)])
M <- S4[, ids]; M <- M[rowSds(M) > 0, ]; Z <- (M - rowMeans(M)) / rowSds(M)
S <- as.data.frame(sapply(names(MOD), function(k) colMeans(Z[intersect(MOD[[k]], rownames(Z)), , drop = FALSE])))
S$Barrier <- rowMeans(S[, c("CAF_ECM", "TAM_SUP")]); S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")])
S$NetImmune <- S$Immune - S$Barrier; S$CXCL9_SPP1 <- Z["CXCL9", ] - Z["SPP1", ]
logit_eff <- function(y, s) {             # E5_icb_split.R, unchanged
  s <- (s - mean(s)) / sd(s)
  f <- suppressWarnings(glm(y ~ s, family = binomial()))
  c(b = unname(coef(f)[2]), se = unname(sqrt(diag(vcov(f)))[2]))
}
per <- read.csv(file.path(O, "E5_icb_per_cohort.csv"))
cat(sprintf("Braun response-labelled patients: %d (responders %d); E5 used n = %d, responders %d\n",
            length(y), sum(y), per$n[per$cohort == "ICB_Braun"][1], per$n_resp[per$cohort == "ICB_Braun"][1]))
out <- do.call(rbind, lapply(SC, function(k) {
  e <- logit_eff(y, S[[k]]); old <- per[per$cohort == "ICB_Braun" & per$score == k, ]
  p_old <- per[per$score == k, ]; p_new <- p_old; p_new$logOR[p_new$cohort == "ICB_Braun"] <- e[["b"]]
  p_new$se[p_new$cohort == "ICB_Braun"] <- e[["se"]]
  m_old <- rma(yi = p_old$logOR, sei = p_old$se, method = "DL"); m_new <- rma(yi = p_new$logOR, sei = p_new$se, method = "DL")
  data.frame(score = k, braun_logOR_compendium = old$logOR, braun_se_compendium = old$se,
             braun_logOR_S4A = e[["b"]], braun_se_S4A = e[["se"]],
             pooled11_DL_compendium = m_old$b[1], pooled11_lo_compendium = m_old$ci.lb, pooled11_hi_compendium = m_old$ci.ub,
             pooled11_DL_S4A = m_new$b[1], pooled11_lo_S4A = m_new$ci.lb, pooled11_hi_S4A = m_new$ci.ub)
}))
write.csv(out, file.path(O, "E21_braun_posthoc_11cohort.csv"), row.names = FALSE)
options(width = 220); o <- out; o[, -1] <- round(o[, -1], 3); print(o, row.names = FALSE)
