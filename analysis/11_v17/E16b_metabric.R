## E16b - METABRIC bulk prognosis (Curtis 2012; cBioPortal brca_metabric bundle downloaded
## 2026-09-24). E11 could only use TCGA-BRCA as a stand-in because the local METABRIC file held
## 3 of the 48 module genes.
##
## Expected directions, written before any METABRIC outcome was computed (from the gastric
## cohorts, TCGA pan-cancer E11 and the TCGA-BRCA stand-in):
##   CAF_ECM    HR > 1 (adverse), small
##   NetImmune  HR < 1 (protective)
##   Immune, IFN, TLS_B, T_CELL  HR < 1, stronger in ER- than ER+ tumours
##   TAM_SUP, C1Q_arm, SPP1_arm  near 1
##   CXCL9_SPP1 HR < 1
##
## Expression: Illumina HT-12 v3, "log2 intensity" per the bundle's meta file (checked below);
## probes are already collapsed to Hugo symbols; remaining duplicate symbols collapsed to the row
## with the highest mean. Scoring as in the manuscript (row-z within cohort; module = mean z of
## >= 3 genes). Endpoints: overall survival (OS_MONTHS / OS_STATUS) and disease-specific survival
## (VITAL_STATUS "Died of Disease" = event; other deaths censored). Cox HR per SD, unadjusted and
## adjusted for age, Nottingham Prognostic Index and ER status; also within ER+ and ER- strata
## (adjusted for age and NPI).
suppressMessages({library(data.table); library(survival); library(matrixStats)})
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
D <- file.path(O, "E16_metabric_tmp/brca_metabric")

MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"),
  C1Q_arm = c("C1QA","C1QB","C1QC"),
  SPP1_arm = c("SPP1","TREM2","MARCO","APOE","LGALS3"))
SC <- c("CAF_ECM", "TAM_SUP", "APM", "IFN", "TLS_B", "T_CELL", "Barrier", "Immune", "NetImmune",
        "C1Q_arm", "SPP1_arm", "CXCL9_SPP1")

x <- fread(file.path(D, "data_mrna_illumina_microarray.txt"), data.table = FALSE)
x <- x[!is.na(x$Hugo_Symbol) & x$Hugo_Symbol != "", ]
M <- as.matrix(x[, -(1:2)]); rownames(M) <- x$Hugo_Symbol
M <- M[rowSums(is.na(M)) == 0, ]
o <- order(-rowMeans(M)); M <- M[o, ]; M <- M[!duplicated(rownames(M)), ]
cat("expression:", nrow(M), "genes x", ncol(M), "samples; range", signif(range(M), 4), "median", signif(median(M), 4),
    "(log2 intensity if roughly 5-15)\n")

pat <- fread(file.path(D, "data_clinical_patient.txt"), skip = 4, data.table = FALSE)
smp <- fread(file.path(D, "data_clinical_sample.txt"), skip = 4, data.table = FALSE)
cl <- merge(pat, smp[, c("PATIENT_ID", "SAMPLE_ID", "ER_STATUS", "TUMOR_STAGE", "TUMOR_SIZE", "GRADE")], by = "PATIENT_ID")
cl <- cl[cl$SAMPLE_ID %in% colnames(M), ]
cl$os_t <- as.numeric(cl$OS_MONTHS); cl$os_e <- as.integer(grepl("^1", cl$OS_STATUS))
cl$dss_e <- as.integer(cl$VITAL_STATUS == "Died of Disease")
cl$age <- as.numeric(cl$AGE_AT_DIAGNOSIS); cl$npi <- as.numeric(cl$NPI)
cl$er <- factor(ifelse(cl$ER_STATUS == "Positive", "pos", ifelse(cl$ER_STATUS == "Negative", "neg", NA)))
cat("clinical matched:", nrow(cl), "| OS events", sum(cl$os_e[!is.na(cl$os_t)]), "| DSS events", sum(cl$dss_e[!is.na(cl$os_t)]),
    "| ER", paste(names(table(cl$er)), table(cl$er), collapse = " "), "\n")

M <- M[, cl$SAMPLE_ID]
M <- M[rowSds(M) > 0, ]; Z <- (M - rowMeans(M)) / rowSds(M)
used <- sapply(MOD, function(g) sum(g %in% rownames(Z)))
cat("module genes present:\n"); print(used)
cat("missing:", setdiff(unlist(MOD), rownames(Z)), "\n")
S <- as.data.frame(sapply(names(MOD), function(k) { g <- intersect(MOD[[k]], rownames(Z))
  if (length(g) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[g, , drop = FALSE]) }))
S$Barrier <- rowMeans(S[, c("CAF_ECM", "TAM_SUP")]); S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")])
S$NetImmune <- S$Immune - S$Barrier
S$CXCL9_SPP1 <- Z["CXCL9", ] - Z["SPP1", ]
d <- cbind(cl, S)

cox1 <- function(dd, k, endpoint, adj) {
  ev <- if (endpoint == "OS") dd$os_e else dd$dss_e
  s <- (dd[[k]] - mean(dd[[k]], na.rm = TRUE)) / sd(dd[[k]], na.rm = TRUE)
  df <- data.frame(t = dd$os_t, e = ev, s = s, age = dd$age, npi = dd$npi, er = dd$er)
  fml <- switch(adj, none = "Surv(t, e) ~ s", full = "Surv(t, e) ~ s + age + npi + er", age_npi = "Surv(t, e) ~ s + age + npi")
  df <- df[complete.cases(df[, all.vars(as.formula(fml))[-(1:2)], drop = FALSE]) & !is.na(df$t) & df$t > 0, ]
  f <- coxph(as.formula(fml), data = df)
  b <- coef(f)[["s"]]; se <- sqrt(vcov(f)["s", "s"])
  data.frame(score = k, endpoint = endpoint, model = adj, n = nrow(df), events = sum(df$e),
             HR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))
}
rows <- list()
for (ep in c("OS", "DSS")) for (k in SC) {
  rows[[length(rows) + 1]] <- cbind(stratum = "all", cox1(d, k, ep, "none"))
  rows[[length(rows) + 1]] <- cbind(stratum = "all", cox1(d, k, ep, "full"))
  for (e in c("pos", "neg")) rows[[length(rows) + 1]] <- cbind(stratum = paste0("ER", e), cox1(d[d$er %in% e, ], k, ep, "age_npi"))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(O, "E16_metabric_cox.csv"), row.names = FALSE)
cat("\n== METABRIC, HR per SD\n"); print(res, digits = 3, row.names = FALSE)

## comparison with the TCGA-BRCA stand-in (E11: age, stage and inferred ER)
tb <- read.csv(file.path(O, "E11_brca_substitute.csv"), stringsAsFactors = FALSE)
cmp <- merge(res[res$stratum == "all" & res$model == "full", c("score", "endpoint", "HR", "lo", "hi", "p")],
             tb[, c("score", "endpoint", "HR", "lo", "hi", "p")], by = c("score", "endpoint"), suffixes = c("_METABRIC", "_TCGA_BRCA"))
write.csv(cmp, file.path(O, "E16_metabric_vs_tcga_brca.csv"), row.names = FALSE)
cat("\n== METABRIC (adjusted) vs TCGA-BRCA stand-in\n"); print(cmp, digits = 3, row.names = FALSE)
cat("\nsign agreement:", sum(sign(log(cmp$HR_METABRIC)) == sign(log(cmp$HR_TCGA_BRCA))), "of", nrow(cmp), "\n")
