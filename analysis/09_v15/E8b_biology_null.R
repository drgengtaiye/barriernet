## E8b - biology-matched null, split out of E8 (see E8_prognostic_hardening.R header, item 5).
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

rd <- function(f) { d <- read.csv(f); unique(na.omit(d$gene_name)) }
## ---------------- 5. biology-matched null ----------------
data(SI_geneset, package = "estimate", envir = environment())
est <- unlist(SI_geneset[c("ImmuneSignature", "StromalSignature"), -1]); est <- as.character(est[est != "" & !is.na(est)])
EXCL <- c("Base_excision_repair", "CellCycle_Reg", "Cell_cycle", "DDR", "DNA_replication", "FGFR3_related", "GPAGs", "PPAGs",
          "Histones", "Homologous_recombination", "Mismatch_Repair", "Nucleotide_excision_repair", "SW480_cancer_cells_Bindea_et_al",
          "Normal_mucosa_Bindea_et_al", "WNT_target", "CIN25_Carter", "CIN70_Carter", "proliferation_Pabla", "KDM5A_Wang",
          "PTEN-MITF_Cabrita", "WNT-beta-catenin_Takeuchi", "Blood_Friedlander", "Tcell-exclusion_Arnon")
files <- c(list.files(file.path(SIG, "IO"), pattern = "\\.csv$", full.names = TRUE, recursive = TRUE),
           list.files(file.path(SIG, "TME"), pattern = "\\.csv$", full.names = TRUE, recursive = TRUE))
files <- files[!sub("\\.csv$", "", basename(files)) %in% EXCL & !basename(dirname(files)) %in% EXCL]
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
