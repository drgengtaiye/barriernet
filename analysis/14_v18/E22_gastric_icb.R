## E22 - gastric and gastro-oesophageal checkpoint cohorts (pre-registered: analysis/E17_E22_prereg.md, section E22
## and the Deviations entries "2026-10-05: E22 data ..."; further deviations in analysis/results/E22_README.md, written
## before this script was run on the full cohort). Run after the PRJEB60680 quantification QC
## (data_download/validation/bulk_PRJEB60680_Samsung_NCT04249739/quant/aggregate_qc.py).
##
## Cohorts
##   PRJEB60680 (Samsung NCT04249739; An 2024 Cancer Discov; clinical data Wright 2026 Cell Rep Med Table S1):
##     5-FU/platinum x 1 cycle, then pembrolizumab added; baseline (pre-treatment) biopsies; fastp + salmon 2.0.0 on
##     Ensembl 110 cDNA; gene-symbol TPM re-normalised over nuclear protein-coding + IG/TR genes; log2(TPM + 1) once.
##     Two library chemistries (TruSeq-type HiSeq 2x101; Nextera-type NovaSeq/NextSeq 2x151/2x149).
##     Response = best overall response CR/PR vs SD/PD.
##   GSE183924 (BTCRC-ESO14-012, adjuvant durvalumab, EAC/GEJ, 37 patients): edgeR FPKM, linear -> log2(FPKM + 1)
##     once. Endpoint relapse-free survival (reported separately, never pooled with response cohorts).
## Expected (prereg):
##   response cohorts: Immune, IFN, NetImmune log-OR > 0; CAF_ECM <= 0; pooled with PRJEB25780 (+ PERFECT as
##   sensitivity); MSS/EBV-negative: NetImmune log-OR > 0.  GSE183924: Immune, IFN, NetImmune HR < 1; CAF_ECM HR >= 1.
## Scoring as E12/E16 (module = mean z of >= 3 listed genes present; Barrier = mean(CAF_ECM, TAM_SUP); Immune =
## mean(APM, IFN, TLS_B, T_CELL); NetImmune = Immune - Barrier; each score standardised to SD 1 over the analysed
## patients before modelling). PRJEB60680 primary (deviation, see README): gene z computed WITHIN library batch; the
## literal pre-registered version (all-transcript TPM, z over the whole cohort) is reported as a sensitivity analysis.
## Logistic log-odds of response per SD (E12 logit_eff); AUC with DeLong CI (pROC, as E7b); Cox (Efron) per SD.
## Structure-matched random-gene null for the AUCs of NetImmune, Immune and IFN as E3 (random modules of the same sizes
## from the 14,364-gene bulk pool, same z as the observed scores, 2,000 draws, random axes credited with their better
## direction). Pooling: the manuscript's conservative rule (E21_meta.R pool(): REML estimate; the wider of the REML
## Hartung-Knapp and DerSimonian-Laird intervals; the larger P). PRJEB25780 inputs are the E12 per-cohort rows
## (deposited log2 scale); PERFECT inputs are E16 primary_log2_vst rows; both read unchanged.
## Firth / MSI-EBV / CPS models: E22_firth_cps.py (E18 code).
suppressMessages({library(matrixStats); library(metafor); library(survival); library(pROC); library(readxl)})
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2"
RES <- file.path(G, "analysis/results")                # inputs from earlier experiments (E12, E16, bulk pool)
O <- Sys.getenv("E22_OUT", RES)                          # outputs (override only for test runs)
PJ <- file.path(G, "data_download/validation/bulk_PRJEB60680_Samsung_NCT04249739")
GS <- file.path(G, "data_download/validation/GSE183924_adjDurvalumab_EAC")
QCD <- Sys.getenv("E22_QC_OUT", file.path(PJ, "quant"))   # override only for test runs
set.seed(2026)

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
MAIN <- c("CAF_ECM", "TAM_SUP", "APM", "IFN", "TLS_B", "T_CELL", "Barrier", "Immune", "NetImmune")

## ---- helpers copied from E12_icb_expanded.R / E16a_perfect.R (scores_sym unchanged; Z returned for the null)
scores_from_Z <- function(Z) {
  used <- sapply(MOD, function(g) sum(g %in% rownames(Z)))
  S <- as.data.frame(sapply(names(MOD), function(k) {
    g <- intersect(MOD[[k]], rownames(Z)); if (length(g) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[g, , drop = FALSE]) }))
  S$Barrier <- rowMeans(S[, c("CAF_ECM", "TAM_SUP")])
  S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")])
  S$NetImmune <- S$Immune - S$Barrier
  S$CXCL9_SPP1 <- if (all(c("CXCL9", "SPP1") %in% rownames(Z))) Z["CXCL9", ] - Z["SPP1", ] else NA_real_
  rownames(S) <- colnames(Z)
  list(S = S, used = used, Z = Z)
}
scores_sym <- function(M) {                               # E12: row-z over the whole cohort
  M <- M[rowSds(M) > 0, , drop = FALSE]; scores_from_Z((M - rowMeans(M)) / rowSds(M))
}
scores_batch <- function(M, batch) {                      # deviation: row-z within library batch
  M <- M[rowSds(M) > 0, , drop = FALSE]
  Z <- M; Z[] <- NA_real_
  for (b in unique(batch)) {
    i <- which(batch == b); X <- M[, i, drop = FALSE]; sdv <- rowSds(X)
    Zb <- (X - rowMeans(X)) / sdv; Zb[sdv == 0, ] <- 0   # gene constant within a batch carries no information there
    Z[, i] <- Zb
  }
  scores_from_Z(Z)
}
logit_eff <- function(y, s) {
  s <- (s - mean(s)) / sd(s)
  f <- suppressWarnings(glm(y ~ s, family = binomial()))
  c(b = unname(coef(f)[2]), se = unname(sqrt(diag(vcov(f)))[2]))
}
auc <- function(y, v) { r <- rank(v); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
auc_ci <- function(y, v) {
  r <- suppressMessages(roc(y, v, levels = c(0, 1), direction = "<", quiet = TRUE))
  ci <- as.numeric(ci.auc(r, method = "delong")); c(auc = ci[2], lo = ci[1], hi = ci[3])
}
cox_eff <- function(t, e, s, adj = NULL) {
  s <- (s - mean(s)) / sd(s)
  f <- if (is.null(adj)) coxph(Surv(t, e) ~ s, ties = "efron") else coxph(Surv(t, e) ~ s + adj, ties = "efron")
  cf <- summary(f)$coefficients
  c(logHR = cf["s", "coef"], se = cf["s", "se(coef)"], p = cf["s", "Pr(>|z|)"], C = unname(f$concordance["concordance"]))
}
## the manuscript's conservative pooling rule (E21_meta.R pool())
pool <- function(d) {
  note <- ""
  re <- withCallingHandlers(rma(yi = d$y, sei = d$se, method = "REML", test = "knha"),
                            warning = function(w) { note <<- gsub("\n", " ", conditionMessage(w)); invokeRestart("muffleWarning") })
  dl <- rma(yi = d$y, sei = d$se, method = "DL")
  fe <- rma(yi = d$y, sei = d$se, method = "FE")
  hk_wider <- (re$ci.ub - re$ci.lb) >= (dl$ci.ub - dl$ci.lb)
  data.frame(k = nrow(d), k_pos = sum(d$y > 0), est = re$b[1], lo = if (hk_wider) re$ci.lb else dl$ci.lb,
             hi = if (hk_wider) re$ci.ub else dl$ci.ub, p = max(re$pval, dl$pval), interval = if (hk_wider) "HK" else "DL",
             hk_lo = re$ci.lb, hk_hi = re$ci.ub, hk_p = re$pval, tau2_reml = re$tau2, I2 = re$I2,
             dl_est = dl$b[1], dl_lo = dl$ci.lb, dl_hi = dl$ci.ub, dl_p = dl$pval,
             fe_est = fe$b[1], fe_lo = fe$ci.lb, fe_hi = fe$ci.ub, fe_p = fe$pval, reml_warning = note)
}
f2 <- function(x) sprintf("%+.2f", x)
rd <- function(f) as.matrix(read.delim(gzfile(file.path(QCD, f)), row.names = 1, check.names = FALSE))

## =============================================================================================================
## 1. PRJEB60680
## =============================================================================================================
qc <- read.delim(file.path(QCD, "PRJEB60680_baseline_qc.tsv"), stringsAsFactors = FALSE)
qc$QC_pass <- as.character(qc$QC_pass) %in% c("True", "TRUE")       # pandas writes True/False
cl_all <- as.data.frame(read_excel(file.path(PJ, "Wright2026_CellRepMed_TableS1_clinical_47pts.xlsx")))
names(cl_all) <- c("id", "scr", "sex", "age", "BOR", "fast", "pfs_ev", "os_ev", "pfs_d", "os_d", "HER2", "EBV", "PDL1",
                   "cps_B", "cps_F1", "cps_F2", "MSI", "TCGA", "scRNA")
cl_all$y <- as.integer(cl_all$BOR %in% c("CR", "PR"))
cl_all$msi_ebv <- as.integer(cl_all$MSI == "MSI-H" | cl_all$EBV == "positive")
cl_all$msi_ebv_tcga <- as.integer(cl_all$msi_ebv == 1 | cl_all$TCGA %in% c(3, 4))   # sensitivity: + E02 (MSI-L, TCGA MSI)
cl_all$cps <- as.numeric(cl_all$cps_B)
cl_all$her2 <- as.integer(cl_all$HER2 == "positive")
cl_all$batch <- qc$lib_batch[match(cl_all$id, qc$sample)]
cl_all$nextera <- as.integer(cl_all$batch == "Nextera-type")
keep <- qc$sample[qc$QC_pass]
excl <- qc[!qc$QC_pass, c("sample", "exclusion_reason")]
cat("PRJEB60680: QC-pass baseline samples", length(keep), "of", nrow(qc), "\n"); print(excl, row.names = FALSE)
cl <- cl_all[match(keep, cl_all$id), ]; stopifnot(identical(cl$id, keep), all(cl$BOR %in% c("CR", "PR", "SD", "PD")), !anyNA(cl$batch))

TPMNC <- rd("PRJEB60680_baseline_gene_TPM_nuclear_coding.tsv.gz")         # primary (re-normalised)
TPMALL <- rd("PRJEB60680_baseline_gene_TPM.tsv.gz")                        # literal (all transcripts)
TPMUNC <- rd("PRJEB60680_baseline_gene_TPM_untrimmed_nuclear_coding.tsv.gz")
E36NC <- rd("PRJEB60680_E36_single_mate_gene_TPM_nuclear_coding.tsv.gz")
sc_chk <- data.frame(cohort = "PRJEB60680",
  matrix = c("gene TPM, all transcripts (salmon, summed by symbol)", "gene TPM, nuclear protein-coding + IG/TR (primary)", "log2(TPM + 1) used (primary)"),
  min = c(min(TPMALL[, keep]), min(TPMNC[, keep]), min(log2(TPMNC[, keep] + 1))),
  median = c(median(TPMALL[, keep]), median(TPMNC[, keep]), median(log2(TPMNC[, keep] + 1))),
  max = c(max(TPMALL[, keep]), max(TPMNC[, keep]), max(log2(TPMNC[, keep] + 1))),
  colsum_min = c(min(colSums(TPMALL[, keep])), min(colSums(TPMNC[, keep])), NA),
  colsum_max = c(max(colSums(TPMALL[, keep])), max(colSums(TPMNC[, keep])), NA),
  verdict = c("linear (column sums 1e6)", "linear (re-normalised to 1e6)", "log2 applied once"))

X1 <- log2(TPMNC[, keep] + 1)
s1 <- scores_batch(X1, cl$batch)                          # PRIMARY
S1 <- s1$S
cat(sprintf("response CR/PR %d vs SD/PD %d; MSI-H/EBV+ %d (all responders? %s); MSS/EBV-negative %d (%d responders)\n",
            sum(cl$y), sum(1 - cl$y), sum(cl$msi_ebv), all(cl$y[cl$msi_ebv == 1] == 1),
            sum(cl$msi_ebv == 0), sum(cl$y[cl$msi_ebv == 0])))
cat("library batch x best response:\n"); print(table(cl$batch, cl$BOR))
cov1 <- data.frame(cohort = "PRJEB60680", module = names(s1$used), genes_listed = sapply(MOD, length), genes_used = s1$used,
                   missing = sapply(MOD, function(g) paste(setdiff(g, rownames(s1$Z)), collapse = ";")))

## sensitivity score sets
S_lit <- scores_sym(log2(TPMALL[, keep] + 1))$S                                  # (a) literally as pre-registered
S_nc_cohortz <- scores_sym(X1)$S                                                  # (b) nuclear coding, cohort-wide z (+ batch term)
S_untr <- scores_batch(log2(TPMUNC[, keep] + 1), cl$batch)$S                     # (d) untrimmed reads
ids_e36 <- c(keep, "E36")                                                         # (e) + E36 (R2-only library)
cl_e36 <- cl_all[match(ids_e36, cl_all$id), ]; cl_e36$batch[cl_e36$id == "E36"] <- "Nextera-type"
S_e36 <- scores_batch(log2(cbind(TPMNC[, keep], E36 = E36NC[rownames(TPMNC), "E36_R2only"]) + 1), cl_e36$batch)$S
excl_other <- setdiff(excl$sample, "E36")                                        # (f) + samples excluded by the QC rules
if (length(excl_other)) {
  ids_x <- c(keep, excl_other); cl_x <- cl_all[match(ids_x, cl_all$id), ]
  S_x <- scores_batch(log2(TPMNC[, ids_x] + 1), cl_x$batch)$S
}
AG <- data.frame(score = SC,
  rho_literal = sapply(SC, function(k) cor(S1[keep, k], S_lit[keep, k], method = "spearman")),
  rho_cohortz = sapply(SC, function(k) cor(S1[keep, k], S_nc_cohortz[keep, k], method = "spearman")),
  rho_untrimmed = sapply(SC, function(k) cor(S1[keep, k], S_untr[keep, k], method = "spearman")),
  batch_R2_literal = sapply(SC, function(k) summary(lm(S_lit[keep, k] ~ cl$batch))$r.squared),
  batch_R2_cohortz = sapply(SC, function(k) summary(lm(S_nc_cohortz[keep, k] ~ cl$batch))$r.squared))
write.csv(AG, file.path(O, "E22_PRJEB60680_score_agreement.csv"), row.names = FALSE)
cat("\nscore agreement with the primary scores (Spearman) and share of score variance explained by batch:\n"); print(AG, digits = 3, row.names = FALSE)

out1 <- cbind(cl[, c("id", "BOR", "y", "pfs_d", "pfs_ev", "os_d", "os_ev", "MSI", "EBV", "TCGA", "msi_ebv", "msi_ebv_tcga",
                     "cps", "her2", "batch", "nextera", "sex", "age")], S1[cl$id, SC])
write.csv(out1, file.path(O, "E22_scores_PRJEB60680.csv"), row.names = FALSE)

resp_rows <- function(d, S, label, adj = NULL) {
  if (sum(d$y) < 3 || sum(1 - d$y) < 3 || (!is.null(adj) && length(unique(d[[adj]])) < 2)) {
    cat("skipped (too few patients per class or constant covariate):", label, "\n"); return(NULL) }
  do.call(rbind, lapply(SC, function(k) {
    v <- S[d$id, k]
    if (is.null(adj)) { e <- logit_eff(d$y, v); b <- e[["b"]]; se <- e[["se"]] } else {
      s <- (v - mean(v)) / sd(v); f <- suppressWarnings(glm(d$y ~ s + d[[adj]], family = binomial()))
      b <- unname(coef(f)["s"]); se <- unname(sqrt(diag(vcov(f)))["s"]) }
    a <- auc_ci(d$y, v)
    data.frame(cohort = "PRJEB60680", analysis = label, score = k, n = nrow(d), n_resp = sum(d$y),
               genes = if (k %in% names(s1$used)) s1$used[[k]] else NA, logOR = b, se = se, lo = b - 1.96 * se,
               hi = b + 1.96 * se, p = 2 * pnorm(-abs(b / se)), AUC = auc(d$y, v), AUC_lo = a[["lo"]], AUC_hi = a[["hi"]])
  }))
}
PRIM <- "primary: all QC-pass baseline, z within library batch"
MSSL <- "MSS/EBV-negative (pre-registered subgroup)"
R1 <- rbind(
  resp_rows(cl, S1, PRIM),
  resp_rows(cl[cl$msi_ebv == 0, ], S1, MSSL),
  resp_rows(cl, S_lit, "sensitivity (a): as pre-registered literally (all-transcript TPM, z over cohort, no batch term)"),
  resp_rows(cl[cl$msi_ebv == 0, ], S_lit, "sensitivity (a): as pre-registered literally, MSS/EBV-negative"),
  resp_rows(cl, S_nc_cohortz, "sensitivity (b): nuclear-coding TPM, z over cohort, batch covariate", adj = "nextera"),
  resp_rows(cl, S1, "sensitivity: primary scores + batch covariate", adj = "nextera"),
  resp_rows(cl[cl$batch == "TruSeq-type", ], S1, "sensitivity (c): TruSeq-type batch only (HiSeq 2x101)"),
  resp_rows(cl[cl$batch == "Nextera-type", ], S1, "sensitivity (c): Nextera-type batch only (NovaSeq/NextSeq)"),
  resp_rows(cl, S_untr, "sensitivity (d): untrimmed reads"),
  resp_rows(cl_e36, S_e36, "sensitivity (e): + E36 (R2-only single-end library)"),
  if (length(excl_other)) resp_rows(cl_x, S_x, paste0("sensitivity (f): + QC-excluded ", paste(excl_other, collapse = ","))),
  resp_rows(cl[cl$msi_ebv_tcga == 0, ], S1, "sensitivity: MSS/EBV-negative also excluding E02 (TCGA subtype MSI)"),
  resp_rows(cl[cl$her2 == 0, ], S1, "sensitivity: HER2-negative only (no trastuzumab)"))
write.csv(R1, file.path(O, "E22_PRJEB60680_response.csv"), row.names = FALSE)
cat("\n== PRJEB60680 response (CR/PR vs SD/PD), log-OR per SD\n")
for (a in unique(R1$analysis)) { cat("--", a, "\n"); q <- R1[R1$analysis == a & R1$score %in% MAIN, ]
  for (i in seq_len(nrow(q))) with(q[i, ], cat(sprintf("  %-10s n=%d/%d  %s (%s to %s) P=%.3g  AUC %.3f (%.3f-%.3f)\n",
    score, n_resp, n, f2(logOR), f2(lo), f2(hi), p, AUC, AUC_lo, AUC_hi))) }

## ---- survival in PRJEB60680 (PFS, OS; secondary, not pre-registered for this cohort)
C1 <- do.call(rbind, lapply(list(c("PFS", "pfs_d", "pfs_ev"), c("OS", "os_d", "os_ev")), function(ep)
  do.call(rbind, lapply(SC, function(k) {
    e <- cox_eff(cl[[ep[2]]], cl[[ep[3]]], S1[cl$id, k]); el <- cox_eff(cl[[ep[2]]], cl[[ep[3]]], S_lit[cl$id, k])
    data.frame(cohort = "PRJEB60680", endpoint = ep[1], score = k, n = nrow(cl), events = sum(cl[[ep[3]]]),
               HR = exp(e[["logHR"]]), lo = exp(e[["logHR"]] - 1.96 * e[["se"]]), hi = exp(e[["logHR"]] + 1.96 * e[["se"]]),
               logHR = e[["logHR"]], se = e[["se"]], p = e[["p"]], C = e[["C"]],
               HR_literal = exp(el[["logHR"]]), p_literal = el[["p"]]) }))))
write.csv(C1, file.path(O, "E22_PRJEB60680_cox.csv"), row.names = FALSE)
cat("\n== PRJEB60680 Cox per SD (primary scores; literal-version HR in brackets)\n")
for (i in which(C1$score %in% MAIN)) with(C1[i, ], cat(sprintf("  %-3s %-10s HR %.2f (%.2f-%.2f) P=%.3g  events %d/%d  [literal HR %.2f P=%.3g]\n",
  endpoint, score, HR, lo, hi, p, events, n, HR_literal, p_literal)))

## ---- structure-matched random-gene null (E3 design) for NetImmune / Immune / IFN AUC, same z as the primary scores
NRAND <- 2000
pool_all <- scan(file.path(RES, "bulk_gene_pool.txt"), what = "", quiet = TRUE)
Zs <- s1$Z[, cl$id]
pool_g <- intersect(pool_all, rownames(Zs))
Zp <- Zs[pool_g, ]
si <- sapply(MOD[c("APM", "IFN", "TLS_B", "T_CELL")], length); sb <- sapply(MOD[c("CAF_ECM", "TAM_SUP")], length)
axis_score <- function(idx, sizes) { k <- 0; out <- 0
  for (s in sizes) { out <- out + colMeans(Zp[idx[(k + 1):(k + s)], , drop = FALSE]); k <- k + s }; out / length(sizes) }
y1 <- cl$y
obs <- c(NetImmune = auc(y1, S1[cl$id, "NetImmune"]), Immune = auc(y1, S1[cl$id, "Immune"]), IFN = auc(y1, S1[cl$id, "IFN"]))
draws <- t(replicate(NRAND, {
  g <- sample(nrow(Zp), sum(si) + sum(sb))
  imm <- axis_score(g[1:sum(si)], si); bar <- axis_score(g[(sum(si) + 1):(sum(si) + sum(sb))], sb)
  ifn <- colMeans(Zp[sample(nrow(Zp), 5), ])
  a <- c(auc(y1, imm - bar), auc(y1, imm), auc(y1, ifn)); pmax(a, 1 - a) }))
colnames(draws) <- names(obs)
write.csv(draws, file.path(O, "E22_PRJEB60680_null_draws.csv"), row.names = FALSE)
N1 <- data.frame(cohort = "PRJEB60680", score = names(obs), AUC_obs = obs, null_median = apply(draws, 2, median),
                 null_p95 = apply(draws, 2, quantile, 0.95), p_emp = (colSums(sweep(draws, 2, obs, ">=")) + 1) / (NRAND + 1),
                 pool_genes = length(pool_g), draws = NRAND)
write.csv(N1, file.path(O, "E22_PRJEB60680_null.csv"), row.names = FALSE)
cat("\n== structure-matched null (PRJEB60680; pool", length(pool_g), "genes)\n"); print(N1[, -1], digits = 3, row.names = FALSE)

## =============================================================================================================
## 2. pooling with PRJEB25780 (primary) and + PERFECT (sensitivity)
## =============================================================================================================
e12 <- read.csv(file.path(RES, "E12_icb_per_cohort.csv"), stringsAsFactors = FALSE)
pj25 <- e12[e12$cohort == "PRJEB25780", c("cohort", "score", "n", "n_resp", "logOR", "se", "AUC")]
e16 <- read.csv(file.path(RES, "E16_perfect_per_score.csv"), stringsAsFactors = FALSE)
pf <- e16[e16$scale == "primary_log2_vst", c("cohort", "score", "n", "n_resp", "logOR", "se", "AUC")]
pj60 <- R1[R1$analysis == PRIM, c("cohort", "score", "n", "n_resp", "logOR", "se", "AUC")]
pj60l <- R1[startsWith(R1$analysis, "sensitivity (a): as pre-registered literally ("), c("cohort", "score", "n", "n_resp", "logOR", "se", "AUC")]
pj60l$cohort <- "PRJEB60680 (literal)"
IN <- rbind(pj60, pj25, pf, pj60l); names(IN)[names(IN) == "logOR"] <- "y"
IN$p_cohort <- 2 * pnorm(-abs(IN$y / IN$se))
write.csv(IN, file.path(O, "E22_meta_inputs.csv"), row.names = FALSE)
SETS <- list("PRJEB60680 + PRJEB25780 (primary)" = c("PRJEB60680", "PRJEB25780"),
             "PRJEB60680 + PRJEB25780 + PERFECT (sensitivity)" = c("PRJEB60680", "PRJEB25780", "PERFECT_GSE165252"),
             "PRJEB60680 literal + PRJEB25780 (sensitivity)" = c("PRJEB60680 (literal)", "PRJEB25780"))
P <- do.call(rbind, lapply(names(SETS), function(s) do.call(rbind, lapply(SC, function(k)
  cbind(set = s, score = k, pool(IN[IN$score == k & IN$cohort %in% SETS[[s]], ]))))))
write.csv(P, file.path(O, "E22_pooled.csv"), row.names = FALSE)
cat("\n== pooled log-OR per SD (conservative rule: REML estimate, wider of HK / DL interval, larger P)\n")
for (s in names(SETS)) { cat("--", s, "\n"); q <- P[P$set == s & P$score %in% MAIN, ]
  for (i in seq_len(nrow(q))) with(q[i, ], cat(sprintf("  %-10s k=%d pos %d  %s (%s to %s; %s) P=%.3g | DL %s (%s to %s) P=%.3g | FE P=%.3g  I2=%.0f%%\n",
    score, k, k_pos, f2(est), f2(lo), f2(hi), interval, p, f2(dl_est), f2(dl_lo), f2(dl_hi), dl_p, fe_p, I2))) }

## =============================================================================================================
## 3. GSE183924 (adjuvant durvalumab; RFS)
## =============================================================================================================
fp <- read.delim(gzfile(file.path(GS, "GSE183924_FPKM.txt.gz")), check.names = FALSE, stringsAsFactors = FALSE)
names(fp)[1] <- "gene"; stopifnot(!any(duplicated(fp$gene)))
FP <- as.matrix(fp[, -1]); rownames(FP) <- fp$gene
gc <- as.data.frame(read_excel(file.path(GS, "GSE183924_Clinical_Data_with_CIBERSORT_analysis_for_GEO.xlsx")))
gc <- gc[, c("Patient ID", "RNA-Seq ID", "Relapse", "RFS Censor", "RFS (months)")]
names(gc) <- c("patient", "rna", "relapse", "rfs_ev", "rfs_m")
stopifnot(all(gc$rna %in% colnames(FP)), all(gc$rfs_ev[gc$relapse == "Yes"] == 1))
cat("\nGSE183924: FPKM columns without a clinical row (excluded):", setdiff(colnames(FP), gc$rna), "\n")
FP <- FP[, gc$rna]
sc_chk <- rbind(sc_chk, data.frame(cohort = "GSE183924", matrix = c("FPKM (edgeR, deposited)", "log2(FPKM + 1) used"),
  min = c(min(FP), min(log2(FP + 1))), median = c(median(FP), median(log2(FP + 1))), max = c(max(FP), max(log2(FP + 1))),
  colsum_min = c(min(colSums(FP)), NA), colsum_max = c(max(colSums(FP)), NA),
  verdict = c("linear (max 1.8e5, median < 1)", "log2 applied once")))
write.csv(sc_chk, file.path(O, "E22_scale_checks.csv"), row.names = FALSE)
s2 <- scores_sym(log2(FP + 1))
cov2 <- data.frame(cohort = "GSE183924", module = names(s2$used), genes_listed = sapply(MOD, length), genes_used = s2$used,
                   missing = sapply(MOD, function(g) paste(setdiff(g, rownames(s2$Z)), collapse = ";")))
write.csv(rbind(cov1, cov2), file.path(O, "E22_gene_coverage.csv"), row.names = FALSE)
S2 <- s2$S
write.csv(cbind(gc, S2[gc$rna, SC]), file.path(O, "E22_scores_GSE183924.csv"), row.names = FALSE)
C2 <- do.call(rbind, lapply(SC, function(k) {
  e <- cox_eff(gc$rfs_m, gc$rfs_ev, S2[gc$rna, k])
  data.frame(cohort = "GSE183924", endpoint = "RFS", score = k, n = nrow(gc), events = sum(gc$rfs_ev),
             genes = if (k %in% names(s2$used)) s2$used[[k]] else NA,
             HR = exp(e[["logHR"]]), lo = exp(e[["logHR"]] - 1.96 * e[["se"]]), hi = exp(e[["logHR"]] + 1.96 * e[["se"]]),
             logHR = e[["logHR"]], se = e[["se"]], p = e[["p"]], C = e[["C"]]) }))
write.csv(C2, file.path(O, "E22_GSE183924_cox.csv"), row.names = FALSE)
cat("\n== GSE183924 Cox RFS per SD (", nrow(gc), "patients,", sum(gc$rfs_ev), "events)\n")
for (i in which(C2$score %in% MAIN)) with(C2[i, ], cat(sprintf("  %-10s HR %.2f (%.2f-%.2f) P=%.3g\n", score, HR, lo, hi, p)))

## =============================================================================================================
## 4. expectations
## =============================================================================================================
r <- function(a, k) R1[R1$analysis == a & R1$score == k, ]
pp <- function(s, k) P[P$set == s & P$score == k, ]
c2 <- function(k) C2[C2$score == k, ]
EX <- list()
add <- function(item, expected, est, lo, hi, p, held) {
  hr <- grepl("HR", expected)
  EX[[length(EX) + 1]] <<- data.frame(item = item, expected = expected, estimate = est, lo = lo, hi = hi, p = p,
    direction_held = held, ci_excludes_null = if (hr) (lo > 1 | hi < 1) else (lo > 0 | hi < 0))
}
LIT <- "sensitivity (a): as pre-registered literally (all-transcript TPM, z over cohort, no batch term)"
for (a in c(PRIM, LIT)) {
  tag <- if (a == PRIM) "PRJEB60680" else "PRJEB60680 [literal]"
  for (k in c("Immune", "IFN", "NetImmune")) { q <- r(a, k); add(paste(tag, k, "log-OR"), "> 0", q$logOR, q$lo, q$hi, q$p, q$logOR > 0) }
  q <- r(a, "CAF_ECM"); add(paste(tag, "CAF_ECM log-OR"), "<= 0", q$logOR, q$lo, q$hi, q$p, q$logOR <= 0)
}
q <- r(MSSL, "NetImmune"); add("PRJEB60680 MSS/EBV-negative NetImmune log-OR", "> 0", q$logOR, q$lo, q$hi, q$p, q$logOR > 0)
q <- r("sensitivity (a): as pre-registered literally, MSS/EBV-negative", "NetImmune")
add("PRJEB60680 [literal] MSS/EBV-negative NetImmune log-OR", "> 0", q$logOR, q$lo, q$hi, q$p, q$logOR > 0)
for (s in names(SETS)) {
  for (k in c("Immune", "IFN", "NetImmune")) { q <- pp(s, k); add(paste(s, k, "pooled log-OR"), "> 0", q$est, q$lo, q$hi, q$p, q$est > 0) }
  q <- pp(s, "CAF_ECM"); add(paste(s, "CAF_ECM pooled log-OR"), "<= 0", q$est, q$lo, q$hi, q$p, q$est <= 0) }
for (k in c("Immune", "IFN", "NetImmune")) { q <- c2(k); add(paste("GSE183924", k, "RFS"), "HR < 1", q$HR, q$lo, q$hi, q$p, q$HR < 1) }
q <- c2("CAF_ECM"); add("GSE183924 CAF_ECM RFS", "HR >= 1", q$HR, q$lo, q$hi, q$p, q$HR >= 1)
EX <- do.call(rbind, EX)
write.csv(EX, file.path(O, "E22_expectations.csv"), row.names = FALSE)
cat("\n== expectations\n"); print(EX, digits = 3, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(O, "E22_sessionInfo.txt"))
