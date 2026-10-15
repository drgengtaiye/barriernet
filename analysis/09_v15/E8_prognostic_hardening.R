## E8 - hardening the prognostic statistics (peer-review M2, M3, M7).
##
## 1. Hartung-Knapp random-effects pooling next to DerSimonian-Laird for every pooled HR.
## 2. Likelihood-ratio tests for incremental value over age + stage (+ pN in GSE84437).
## 3. TMEscore (Zeng et al. 2019) as a gastric-specific comparator: TMEscoreA_CIR and
##    TMEscoreB_CIR gene lists from icb/SignatureSets-main; (i) PC1 version - PC1 of the
##    z-scored A genes minus PC1 of the B genes, each oriented to correlate positively with the
##    mean of its genes and scaled to unit variance, as in the published score; (ii) z-mean
##    version - mean z(A) - mean z(B), the same engine as NetImmune. High = favourable.
## 4. Adjustment for ACRG molecular subtype and Lauren type (GSE62254) and for Lauren-type
##    histology (TCGA-STAD, from GDC primary_diagnosis; TCGA molecular subtypes are not on disk).
## 5. Biology-matched null: NetImmune-shaped random axes drawn from immune/stromal signature
##    genes (ESTIMATE immune + stromal, SignatureSets IO + TME minus tumour-intrinsic sets),
##    excluding the 48 barriernet genes, each draw matched to the mean-expression decile of the
##    gene it replaces. Same gene set in all three cohorts; per-cohort |log HR| rank (sign per
##    cohort) and inverse-variance pooled |log HR| (sign chosen once).
## 6. NetImmune versus CAF_ECM: HK-pooled HRs and nested LR tests in both directions.
suppressMessages({library(survival); library(metafor); library(matrixStats); library(readxl)})
set.seed(2026)
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
SIG <- file.path(G, "icb/SignatureSets-main/data-raw")
TAGS <- c("TCGA", "GSE84437", "GSE62254")
NB <- 2000; NRAND <- 1000

res <- readRDS(file.path(R, "B2_scores.rds"))
B1 <- lapply(TAGS, function(t) readRDS(file.path(R, paste0("B1_", t, ".rds")))); names(B1) <- TAGS
zrow <- function(m) { s <- rowSds(m, na.rm = TRUE); s[s == 0 | is.na(s)] <- 1
                      (m - rowMeans(m, na.rm = TRUE)) / s }

## ---------------- TMEscore ----------------
rd <- function(f) { d <- read.csv(f); unique(na.omit(d$gene_name)) }
TA <- rd(file.path(SIG, "TME/TMEscoreA_CIR.csv")); TB <- rd(file.path(SIG, "TME/TMEscoreB_CIR.csv"))
ALIAS <- c(WARS1 = "WARS", HLA_DPB1 = "HLA-DPB1")
cov <- list()
for (t in TAGS) {
  m <- B1[[t]]$expr; rn <- rownames(m)
  fix <- function(g) { g2 <- ifelse(g %in% rn, g, ifelse(!is.na(ALIAS[g]) & ALIAS[g] %in% rn, ALIAS[g], NA)); unique(na.omit(g2)) }
  a <- fix(TA); b <- fix(TB); cov[[t]] <- c(A = length(a), A_total = length(TA), B = length(b), B_total = length(TB))
  Z <- zrow(m[c(a, b), , drop = FALSE])
  pc1 <- function(g) { x <- t(Z[g, , drop = FALSE]); x[is.na(x)] <- 0
                       p <- prcomp(x, center = TRUE, scale. = FALSE)$x[, 1]
                       mz <- colMeans(Z[g, , drop = FALSE], na.rm = TRUE)
                       if (cor(p, mz) < 0) p <- -p; as.numeric(scale(p)) }
  d <- res[[t]]; idx <- match(d$sample, colnames(Z))
  tpc <- pc1(a) - pc1(b); tzm <- colMeans(Z[a, , drop = FALSE], na.rm = TRUE) - colMeans(Z[b, , drop = FALSE], na.rm = TRUE)
  res[[t]]$TMEscore_PC <- tpc[idx]; res[[t]]$TMEscore_z <- tzm[idx]
}
cat("TMEscore gene coverage:\n"); print(do.call(rbind, cov))

## ---------------- scores, direction fixed a priori (risk orientation) ----------------
SC <- c(NetImmune = "NetImmune", Barrier = "Barrier", Immune = "Immune", CD8A = "CD8A",
        ESTIMATE = "ESTsig_diff", NetImmune_ssGSEA = "NetImmune_ss", TMEscore_PC = "TMEscore_PC",
        TMEscore_z = "TMEscore_z", CAF_ECM = "CAF_ECM", TAM_SUP = "TAM_SUP", APM = "APM",
        IFN = "IFN", TLS_B = "TLS_B", T_CELL = "T_CELL")
SIGN <- c(NetImmune = -1, Barrier = 1, Immune = -1, CD8A = -1, ESTIMATE = -1, NetImmune_ssGSEA = -1,
          TMEscore_PC = -1, TMEscore_z = -1, CAF_ECM = 1, TAM_SUP = 1, APM = -1, IFN = -1,
          TLS_B = -1, T_CELL = -1)

## ---------------- 1. HK vs DL pooling ----------------
uni <- do.call(rbind, lapply(TAGS, function(t) {
  d <- res[[t]]
  do.call(rbind, lapply(names(SC), function(k) {
    f <- coxph(Surv(time, event) ~ scale(d[[SC[k]]])[, 1], data = d)
    data.frame(cohort = t, score = k, beta = unname(coef(f)), se = sqrt(vcov(f)[1, 1]))
  }))
}))
meta <- do.call(rbind, lapply(names(SC), function(k) {
  u <- uni[uni$score == k, ]
  dl <- rma(yi = u$beta, sei = u$se, method = "DL")
  hk <- rma(yi = u$beta, sei = u$se, method = "DL", test = "knha")
  data.frame(score = k, HR = exp(dl$b[1]), DL_lo = exp(dl$ci.lb), DL_hi = exp(dl$ci.ub), DL_p = dl$pval,
             HK_lo = exp(hk$ci.lb), HK_hi = exp(hk$ci.ub), HK_p = hk$pval, I2 = dl$I2, tau2 = dl$tau2)
}))
write.csv(uni, file.path(O, "E8_univariate.csv"), row.names = FALSE)
write.csv(meta, file.path(O, "E8_meta_DL_HK.csv"), row.names = FALSE)
cat("\n=== 1. pooled HR per SD: DerSimonian-Laird vs Hartung-Knapp ===\n"); print(format(meta, digits = 3), row.names = FALSE)

## ---------------- 2 & 6. nested likelihood-ratio tests ----------------
CLIN <- list(TCGA = c("age", "stage_ord"), GSE84437 = c("age", "stage_ord", "n_ord"), GSE62254 = c("age", "stage_ord"))
res$GSE84437$n_ord <- as.integer(gsub("N", "", res$GSE84437$n_st))
for (t in TAGS) for (k in names(SC)) res[[t]][[paste0(k, "_z")]] <- as.numeric(scale(res[[t]][[SC[k]]]))
lr <- function(d, base, add) {
  f0 <- if (length(base)) coxph(as.formula(paste("Surv(time, event) ~", paste(base, collapse = " + "))), data = d)
        else coxph(Surv(time, event) ~ 1, data = d)
  f1 <- coxph(as.formula(paste("Surv(time, event) ~", paste(c(base, add), collapse = " + "))), data = d)
  a <- anova(f0, f1); c(chisq = a$Chisq[2], p = a[["Pr(>|Chi|)"]][2])
}
lrt <- do.call(rbind, lapply(TAGS, function(t) {
  cv <- CLIN[[t]]; d <- res[[t]]
  d <- d[complete.cases(d[, c("time", "event", cv, paste0(names(SC), "_z"))]), ]
  rows <- list()
  for (k in c("NetImmune", "ESTIMATE", "CAF_ECM", "TMEscore_PC", "TMEscore_z")) {
    r <- lr(d, cv, paste0(k, "_z")); rows[[length(rows) + 1]] <- data.frame(cohort = t, n = nrow(d), base = "age+stage", added = k, chisq = r["chisq"], p = r["p"])
  }
  d2 <- res[[t]]; d2 <- d2[complete.cases(d2[, c("time", "event", "NetImmune_z", "CAF_ECM_z")]), ]
  for (pr in list(c("CAF_ECM_z", "NetImmune_z"), c("NetImmune_z", "CAF_ECM_z"))) {
    r <- lr(d2, pr[1], pr[2]); rows[[length(rows) + 1]] <- data.frame(cohort = t, n = nrow(d2), base = sub("_z", "", pr[1]), added = sub("_z", "", pr[2]), chisq = r["chisq"], p = r["p"])
    r <- lr(d, c(cv, pr[1]), pr[2]); rows[[length(rows) + 1]] <- data.frame(cohort = t, n = nrow(d), base = paste0("age+stage+", sub("_z", "", pr[1])), added = sub("_z", "", pr[2]), chisq = r["chisq"], p = r["p"])
  }
  do.call(rbind, rows)
}))
rownames(lrt) <- NULL
write.csv(lrt, file.path(O, "E8_lrt.csv"), row.names = FALSE)
cat("\n=== 2/6. nested likelihood-ratio tests ===\n"); print(format(lrt, digits = 3), row.names = FALSE)

## ---------------- 3. head-to-head C: NetImmune vs TMEscore ----------------
hc <- function(tt, ee, r) concordance(Surv(tt, ee) ~ r, reverse = TRUE)$concordance
h2h <- list(); cidx <- list()
for (t in TAGS) {
  d <- res[[t]]; d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$event), ]
  M <- sapply(c("NetImmune", "TMEscore_PC", "TMEscore_z", "CAF_ECM", "ESTIMATE"), function(k) SIGN[k] * d[[SC[k]]])
  C0 <- apply(M, 2, function(r) hc(d$time, d$event, r))
  B <- t(replicate(NB, { i <- sample(nrow(d), replace = TRUE); apply(M[i, ], 2, function(r) hc(d$time[i], d$event[i], r)) }))
  cidx[[t]] <- data.frame(cohort = t, score = names(C0), C = C0)
  for (b in c("TMEscore_PC", "TMEscore_z")) {
    dd <- B[, "NetImmune"] - B[, b]
    h2h[[length(h2h) + 1]] <- data.frame(cohort = t, comparator = b, dC = C0["NetImmune"] - C0[b], se = sd(dd),
                                         lo = quantile(dd, 0.025), hi = quantile(dd, 0.975), p = 2 * min(mean(dd <= 0), mean(dd >= 0)))
  }
}
h2h <- do.call(rbind, h2h); cidx <- do.call(rbind, cidx)
h2h <- rbind(h2h, do.call(rbind, lapply(c("TMEscore_PC", "TMEscore_z"), function(b) {
  s <- h2h[h2h$comparator == b, ]; w <- 1 / s$se^2; e <- sum(w * s$dC) / sum(w); se <- sqrt(1 / sum(w))
  data.frame(cohort = "pooled", comparator = b, dC = e, se = se, lo = e - 1.96 * se, hi = e + 1.96 * se, p = 2 * pnorm(-abs(e / se)))
})))
rownames(h2h) <- NULL; rownames(cidx) <- NULL
write.csv(h2h, file.path(O, "E8_h2h_TMEscore.csv"), row.names = FALSE)
write.csv(cidx, file.path(O, "E8_cindex_TMEscore.csv"), row.names = FALSE)
cat("\n=== 3. C-index and NetImmune minus TMEscore ===\n"); print(format(cidx, digits = 3), row.names = FALSE); print(format(h2h, digits = 3), row.names = FALSE)

## ---------------- 4. subtype and Lauren adjustment ----------------
adj <- list()
## GSE62254: GEO patient number <-> ACRG 'Tumor ID'
ph <- read.csv(file.path(G, "GC_TME_Reprogrammer/data_GEO/GSE62254/GSE62254_pheno.csv"))
ph <- ph[grepl("^T", ph$title), c("geo_accession", "patient.ch1")]
x <- as.data.frame(read_excel(file.path(G, "GC_TME_Reprogrammer/data/41591_2015_BFnm3850_MOESM34_ESM.xls")))
names(x) <- gsub("\n", " ", names(x))
mcol <- grep("^Mol. Subtype", names(x), value = TRUE)
x$subtype <- factor(c("MSS/TP53-", "MSS/TP53+", "MSI", "EMT")[x[[mcol]] + 1], levels = c("MSS/TP53-", "MSS/TP53+", "MSI", "EMT"))
x$lauren <- ifelse(x$Lauren == "intestinal", "intestinal", ifelse(x$Lauren == "diffuse", "diffuse", "mixed/other"))
ph$Tumor.ID <- as.integer(ph$patient.ch1)
d <- res$GSE62254
d <- merge(d, ph, by.x = "sample", by.y = "geo_accession", all.x = TRUE)
d <- merge(d, x[, c("Tumor ID", "subtype", "lauren", "OS (months)")], by.x = "Tumor.ID", by.y = "Tumor ID", all.x = TRUE)
match_rate <- mean(!is.na(d$subtype)); os_agree <- mean(abs(d$time - d$`OS (months)`) < 0.05, na.rm = TRUE)
cat(sprintf("\nGSE62254: matched %d/%d samples to ACRG table (%.1f%%); OS agrees within 0.05 months in %.1f%% of matches\n",
            sum(!is.na(d$subtype)), nrow(d), 100 * match_rate, 100 * os_agree))
print(table(d$subtype, useNA = "ifany")); print(table(d$lauren, useNA = "ifany"))
fitHR <- function(d, rhs, lab, t) {
  f <- coxph(as.formula(paste("Surv(time, event) ~ NetImmune_z +", rhs)), data = d); s <- summary(f)
  data.frame(cohort = t, model = lab, n = s$n, HR = s$conf.int[1, 1], lo = s$conf.int[1, 3], hi = s$conf.int[1, 4], p = s$coefficients[1, 5])
}
adj[[1]] <- fitHR(d, "age + stage_ord", "age + stage", "GSE62254")
adj[[2]] <- fitHR(d, "age + stage_ord + subtype", "age + stage + ACRG subtype", "GSE62254")
adj[[3]] <- fitHR(d, "age + stage_ord + subtype + lauren", "age + stage + ACRG subtype + Lauren", "GSE62254")
adj[[4]] <- fitHR(d, "age + stage_ord + lauren", "age + stage + Lauren", "GSE62254")
f_emt <- coxph(Surv(time, event) ~ NetImmune_z + age + stage_ord + subtype, data = d)
cat("\nGSE62254, NetImmune + age + stage + ACRG subtype:\n"); print(round(summary(f_emt)$conf.int, 3))
## TCGA: Lauren-type histology from GDC clinical
gc <- read.csv(file.path(G, "GC_TME_Reprogrammer/data/raw/TCGA_STAD/TCGA_STAD_clinical.csv"))
gc$hist <- with(gc, ifelse(primary_diagnosis %in% c("Adenocarcinoma, intestinal type", "Tubular adenocarcinoma", "Papillary adenocarcinoma, NOS"), "intestinal",
                    ifelse(primary_diagnosis %in% c("Carcinoma, diffuse type", "Signet ring cell carcinoma"), "diffuse", "other/NOS")))
dt <- res$TCGA; dt$hist <- gc$hist[match(dt$patient, gc$submitter_id)]
cat(sprintf("\nTCGA: histology matched %d/%d\n", sum(!is.na(dt$hist)), nrow(dt))); print(table(dt$hist, useNA = "ifany"))
dt$hist <- factor(dt$hist, levels = c("intestinal", "diffuse", "other/NOS"))
adj[[5]] <- fitHR(dt, "age + stage_ord", "age + stage", "TCGA")
adj[[6]] <- fitHR(dt, "age + stage_ord + hist", "age + stage + Lauren-type histology", "TCGA")
adj <- do.call(rbind, adj); rownames(adj) <- NULL
write.csv(adj, file.path(O, "E8_subtype_adjusted.csv"), row.names = FALSE)
cat("\n=== 4. NetImmune HR adjusted for subtype / Lauren ===\n"); print(format(adj, digits = 3), row.names = FALSE)

## ---------------- 5. biology-matched null ----------------
data(SI_geneset, package = "estimate", envir = environment())
est <- unlist(SI_geneset[c("ImmuneSignature", "StromalSignature"), -1]); est <- as.character(est[est != "" & !is.na(est)])
EXCL <- c("Base_excision_repair", "CellCycle_Reg", "Cell_cycle", "DDR", "DNA_replication", "FGFR3_related", "GPAGs", "PPAGs",
          "Histones", "Homologous_recombination", "Mismatch_Repair", "Nucleotide_excision_repair", "SW480_cancer_cells_Bindea_et_al",
          "Normal_mucosa_Bindea_et_al", "WNT_target", "CIN25_Carter", "CIN70_Carter", "proliferation_Pabla", "KDM5A_Wang",
          "PTEN-MITF_Cabrita", "WNT-beta-catenin_Takeuchi", "Blood_Friedlander", "Tcell-exclusion_Arnon")
files <- c(list.files(file.path(SIG, "IO"), full.names = TRUE), list.files(file.path(SIG, "TME"), full.names = TRUE))
files <- files[!sub("\\.csv$", "", basename(files)) %in% EXCL]
sig_genes <- unique(unlist(lapply(files, rd)))
MOD <- list(APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"), IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
            TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"), T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"),
            CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
            TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"))
pool_all <- readLines(file.path(O, "bulk_gene_pool.txt"))
bio <- setdiff(intersect(unique(c(est, sig_genes)), pool_all), unlist(MOD))
cat(sprintf("\nbiology pool: %d ESTIMATE + %d signature genes -> %d in the 14,364-gene pool after removing the 48 module genes\n",
            length(unique(est)), length(sig_genes), length(bio)))
## mean-expression decile, averaged over cohorts, among the 14,364 pool genes
dec <- rowMeans(sapply(TAGS, function(t) { mu <- rowMeans(B1[[t]]$expr[pool_all, ], na.rm = TRUE); ceiling(10 * rank(mu) / length(mu)) }))
dec <- round(dec); names(dec) <- pool_all
prep <- lapply(TAGS, function(t) {
  s <- res[[t]]; z <- zrow(B1[[t]]$expr[c(bio, intersect(unlist(MOD), pool_all)), ])
  idx <- match(s$sample, colnames(z)); ok <- !is.na(idx) & !is.na(s$time) & s$time > 0 & !is.na(s$event)
  list(z = z[, idx[ok]], time = s$time[ok], event = s$event[ok], net = s$NetImmune[ok])
}); names(prep) <- TAGS
fit <- function(p, x) { f <- coxph(Surv(p$time, p$event) ~ scale(x)[, 1]); c(b = unname(coef(f)), se = sqrt(vcov(f)[1, 1])) }
axis_of <- function(z, mods) { s <- sapply(mods, function(g) colMeans(z[g, , drop = FALSE])); rowMeans(s[, 1:4]) - rowMeans(s[, 5:6]) }
obs <- sapply(TAGS, function(t) fit(prep[[t]], prep[[t]]$net))
mod_in_pool <- lapply(MOD, function(g) intersect(g, pool_all))
bio_by_dec <- split(bio, dec[bio])
draw <- function() {
  used <- character(0)
  lapply(mod_in_pool, function(g) sapply(g, function(x) {
    cand <- setdiff(bio_by_dec[[as.character(dec[x])]], used)
    if (!length(cand)) cand <- setdiff(bio, used)
    y <- sample(cand, 1); used <<- c(used, y); y }))
}
rnd <- array(NA, c(NRAND, 3, 2), dimnames = list(NULL, TAGS, c("b", "se")))
for (k in seq_len(NRAND)) { mods <- draw(); for (t in TAGS) rnd[k, t, ] <- fit(prep[[t]], axis_of(prep[[t]]$z, mods)) }
fe <- function(b, se) { w <- 1 / se^2; sum(w * b) / sum(w) }
obs_pool <- fe(obs["b", ], obs["se", ]); rnd_pool <- apply(rnd, 1, function(m) fe(m[, "b"], m[, "se"]))
bnull <- rbind(
  do.call(rbind, lapply(TAGS, function(t) data.frame(cohort = t, HR_obs = exp(obs["b", t]),
    rank = sum(abs(rnd[, t, "b"]) >= abs(obs["b", t])) + 1, p_emp = (sum(abs(rnd[, t, "b"]) >= abs(obs["b", t])) + 1) / (NRAND + 1),
    null_absHR_p95 = exp(quantile(abs(rnd[, t, "b"]), 0.95))))),
  data.frame(cohort = "pooled", HR_obs = exp(obs_pool), rank = sum(abs(rnd_pool) >= abs(obs_pool)) + 1,
             p_emp = (sum(abs(rnd_pool) >= abs(obs_pool)) + 1) / (NRAND + 1), null_absHR_p95 = exp(quantile(abs(rnd_pool), 0.95))))
rownames(bnull) <- NULL
same <- apply(sign(rnd[, , "b"]), 1, function(s) length(unique(s)) == 1)
write.csv(bnull, file.path(O, "E8_biology_null.csv"), row.names = FALSE)
saveRDS(list(obs = obs, rnd = rnd, rnd_pool = rnd_pool, bio = bio), file.path(O, "E8_biology_null_draws.rds"))
cat("\n=== 5. biology-matched null ===\n"); print(format(bnull, digits = 3), row.names = FALSE)
cat(sprintf("random axes with the same sign in all three cohorts: %.1f%%; median pooled |HR| %.3f\n", 100 * mean(same), exp(median(abs(rnd_pool)))))
