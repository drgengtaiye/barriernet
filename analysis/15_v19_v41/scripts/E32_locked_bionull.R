## E32 Q7 - NetImmune in the seven locked gastric cohorts against structure-matched random axes drawn from immune and
## stromal genes (third/fourth external review, point 2). Protocol: analysis/E32_prereg.md (sha256 in E32_prereg.sha256).
## Method of E8b_biology_null.R (benchmark cohorts), applied to the locked cohorts: the biology pool is E8b's (ESTIMATE
## immune + stromal and SignatureSets IO + TME genes minus tumour-intrinsic sets, minus the 48 module genes) restricted
## to genes measured in all seven cohorts; each module gene is replaced by a pool gene of the same mean-expression decile
## (deciles among the genes measured in all seven, averaged over cohorts); the same draw in every cohort; Cox per cohort;
## inverse-variance pooled log HR, credited with its better direction (|log HR|). Observed NetImmune = the current
## (duplicate-free) patient-level scores. Expression loading as in E31_dedup_all_pairs.R (E15/E20 preparation).
suppressMessages({library(matrixStats); library(readxl); library(survival)})
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
D <- file.path(G, "2/data_download/geo")
V <- file.path(G, "2/data_download/validation/e20_outcomes")
R <- file.path(G, "revision/results")
O <- file.path(G, "2/analysis/results")
sink(file.path(O, "E32_locked_bionull.log"), split = TRUE)
cat("E32 Q7 run", format(Sys.time()), "\n")
h <- system2("sha256sum", file.path(G, "2/analysis/E32_prereg.md"), stdout = TRUE)
stopifnot(sub(" .*", "", h) == sub(" .*", "", readLines(file.path(G, "2/analysis/E32_prereg.sha256"))[1]))
cat("protocol hash verified\n")
set.seed(2032); NRAND <- 1000

zrow <- function(m) { s <- rowSds(m, na.rm = TRUE); s[s == 0 | is.na(s)] <- 1
                      (m - rowMeans(m, na.rm = TRUE)) / s }
collapse_max <- function(mat, sym) {
  keep <- !is.na(sym) & sym != ""
  mat <- mat[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(sym, -rowMeans(mat, na.rm = TRUE)); mat <- mat[o, , drop = FALSE]; sym <- sym[o]
  mat <- mat[!duplicated(sym), , drop = FALSE]; rownames(mat) <- sym[!duplicated(sym)]; mat
}
read_matrix <- function(g) {
  L <- readLines(gzfile(file.path(D, g, paste0(g, "_series_matrix.txt.gz"))))
  b <- grep("^!series_matrix_table_begin", L); e <- grep("^!series_matrix_table_end", L)
  tab <- read.delim(text = L[(b + 1):(e - 1)], check.names = FALSE, quote = "\"")
  m <- as.matrix(tab[, -1]); rownames(m) <- tab[[1]]; storage.mode(m) <- "double"
  list(m = m, meta = L[1:b])
}
meta_field <- function(meta, key) {
  l <- meta[grepl(paste0("^!", key, "\t"), meta)]
  lapply(l, function(x) gsub('"', "", strsplit(x, "\t")[[1]][-1]))
}
prep <- function(g, symmap, keep = NULL) {
  x <- read_matrix(g); m <- x$m
  rng <- range(m, na.rm = TRUE)
  if (rng[2] > 50) m <- log2(pmax(m, 0) + 1)
  if (!is.null(keep)) m <- m[, intersect(keep, colnames(m)), drop = FALSE]
  m <- m[rowSums(is.na(m)) == 0, , drop = FALSE]
  m <- collapse_max(m, unname(symmap[rownames(m)]))
  list(m = m, meta = x$meta)
}
f570 <- read.csv(file.path(G, "GC_TME_Reprogrammer/data_GEO/GSE62254/GSE62254_fdata.csv"), stringsAsFactors = FALSE)
sym570 <- setNames(sub(" ///.*", "", f570$Gene.Symbol), f570$ID)
a8432 <- read.delim(file.path(D, "platforms/GPL8432_id_symbol.tsv"), stringsAsFactors = FALSE)
sym8432 <- setNames(a8432$Symbol, a8432$ID)
a96 <- read.delim(file.path(D, "platforms/GPL96_id_symbol.tsv"), stringsAsFactors = FALSE, check.names = FALSE, quote = "")
sym96 <- setNames(sub(" ///.*", "", a96[["Gene Symbol"]]), a96$ID)
a6884 <- read.delim(file.path(D, "platforms/GPL6884_id_symbol.tsv"), stringsAsFactors = FALSE, quote = "")
sym6884 <- setNames(a6884$Symbol, a6884$ID)
a6947 <- read.delim(file.path(D, "platforms/GPL6947_id_symbol.tsv"), stringsAsFactors = FALSE, quote = "")
sym6947 <- setNames(a6947$Symbol, a6947$ID)

M <- list(); CL <- list()
## benchmark
BTC <- readRDS(file.path(R, "B1_TCGA.rds")); B84 <- readRDS(file.path(R, "B1_GSE84437.rds")); B62 <- readRDS(file.path(R, "B1_GSE62254.rds"))
M[["TCGA-STAD"]] <- BTC$expr; M[["GSE84437"]] <- B84$expr; M[["GSE62254"]] <- B62$expr
## locked 3 (E15)
out_sg <- function(g) {
  o <- as.data.frame(read_excel(file.path(D, g, paste0(g, "_outcome.xls"))))
  names(o) <- c("gsm","id","cel","subtype","age","gender","lauren","stage","os_m","dead")
  o[!is.na(o$gsm) & grepl("^GSM", o$gsm), ]
}
for (g in c("GSE15459", "GSE34942")) {
  p <- prep(g, sym570); o <- out_sg(g)
  o <- o[o$gsm %in% colnames(p$m) & !is.na(as.numeric(o$os_m)) & as.numeric(o$os_m) > 0 & !is.na(o$dead), ]
  M[[g]] <- p$m[, o$gsm]
  CL[[g]] <- data.frame(sample = o$gsm, age = as.numeric(o$age), sex = toupper(substr(o$gender, 1, 1)),
                        time_m = as.numeric(o$os_m), event = as.integer(o$dead), endpoint = "OS",
                        stage = as.character(o$stage), stringsAsFactors = FALSE)
}
p <- prep("GSE26253", sym8432)
st <- meta_field(p$meta, "Sample_characteristics_ch1")
key <- function(k) { i <- which(sapply(st, function(v) any(grepl(k, v, fixed = TRUE)))); sub(".*: ", "", st[[i[1]]]) }
gsm <- meta_field(p$meta, "Sample_geo_accession")[[1]]
c26 <- data.frame(sample = gsm, age = NA_real_, sex = NA_character_,
                  time_m = as.numeric(key("recurrence free survival time")),
                  event = as.integer(key("status (0=non-recurrence")), endpoint = "RFS",
                  stage = key("pathological stage"), stringsAsFactors = FALSE)
c26 <- c26[!is.na(c26$time_m) & c26$time_m > 0 & !is.na(c26$event), ]
M[["GSE26253"]] <- p$m[, c26$sample]
## GSE26253 age/sex: Oh 2018 Nat Commun Suppl Data 2 sheet SMC (same 432 patients); keyed by unique (RFS.m, recurrence, stage roman)
smc <- as.data.frame(read_excel(file.path(V, "oh2018_natcommun/41467_2018_4179_MOESM5_ESM.xlsx"), sheet = "SMC"))
rom <- c("1"="I","1A"="I","1B"="I","2"="II","3A"="III","3B"="III","4"="IV")
k26 <- paste(sprintf("%.1f", c26$time_m), c26$event, rom[c26$stage])
ksm <- paste(sprintf("%.1f", as.numeric(smc$RFS.m)), smc[["Recurrence (1=yes, 0=no)"]], smc$AJCC.stage)
uniq <- k26 %in% names(which(table(k26) == 1)) & k26 %in% names(which(table(ksm) == 1))
j <- match(k26, ksm)
c26$age[uniq] <- as.numeric(smc$Age[j[uniq]]); c26$sex[uniq] <- toupper(substr(smc$Sex[j[uniq]], 1, 1))
c26$os_m_smc <- NA; c26$os_m_smc[uniq] <- as.numeric(smc$OS.m[j[uniq]])
c26$death_smc <- NA; c26$death_smc[uniq] <- as.integer(smc[["Death (1=yes, 0=no)"]][j[uniq]])
cat("GSE26253 -> SMC sheet unique key matches:", sum(uniq), "of", nrow(c26), "\n")
CL[["GSE26253"]] <- c26
## added 4 (E20)
NEW <- c("GSE13861", "GSE26901", "GSE26899", "GSE29272")
SYM <- list(GSE13861 = sym6884, GSE26901 = sym6947, GSE26899 = sym6947, GSE29272 = sym96)
for (g in NEW) {
  oc <- read.csv(file.path(V, paste0(g, "_outcome.csv")), stringsAsFactors = FALSE)
  M[[g]] <- prep(g, SYM[[g]], keep = oc$gsm)$m
  tm <- as.numeric(oc$os_time); if (oc$time_unit[1] == "days") tm <- tm / 30.4375
  CL[[g]] <- data.frame(sample = oc$gsm, age = as.numeric(oc$age), sex = toupper(substr(oc$sex, 1, 1)),
                        time_m = tm, event = as.integer(oc$os_event), endpoint = "OS",
                        stage = as.character(if ("stage_ajcc" %in% names(oc)) oc$stage_ajcc else oc$stage), stringsAsFactors = FALSE)
}
## benchmark clinical (descriptive only)
tc <- read.csv(file.path(G, "GC_TME_Reprogrammer/data/raw/TCGA_STAD/TCGA_STAD_clinical.csv"), stringsAsFactors = FALSE)
CL[["TCGA-STAD"]] <- data.frame(sample = BTC$clin$sample, age = BTC$clin$age,
                                sex = toupper(substr(tc$gender[match(BTC$clin$patient, tc$submitter_id)], 1, 1)),
                                time_m = BTC$clin$time / 30.4375, event = BTC$clin$event, endpoint = "OS",
                                stage = BTC$clin$stage, stringsAsFactors = FALSE)
p84 <- read.csv(file.path(G, "GC_TME_Reprogrammer/data_GEO/GSE84437/GSE84437_pheno.csv"), stringsAsFactors = FALSE)
CL[["GSE84437"]] <- data.frame(sample = p84$geo_accession, age = as.numeric(p84[["age.ch1"]]),
                               sex = toupper(substr(p84[["Sex.ch1"]], 1, 1)),
                               time_m = as.numeric(p84[["duration.overall.survival.ch1"]]), event = as.integer(p84[["death.ch1"]]),
                               endpoint = "OS", stage = paste(p84[["ptstage.ch1"]], p84[["pnstage.ch1"]]), stringsAsFactors = FALSE)
acrg <- as.data.frame(read_excel(file.path(V, "oh2018_natcommun/41467_2018_4179_MOESM5_ESM.xlsx"), sheet = "ACRG"))
CL[["GSE62254"]] <- data.frame(sample = acrg$GEO_ID, age = as.numeric(acrg$age), sex = toupper(substr(acrg$sex, 1, 1)),
                               time_m = as.numeric(acrg$OS.m), event = as.integer(acrg$Death), endpoint = "OS",
                               stage = as.character(acrg$pStage), dfs_m = as.numeric(acrg$DFS.m), recur = acrg$Recur,
                               stringsAsFactors = FALSE)
CLX <- do.call(rbind, lapply(names(CL), function(g) { d <- CL[[g]][, c("sample","age","sex","time_m","event","endpoint","stage")]; d$cohort <- g; d }))
## GSE62254 DFS and GSE26253 SMC OS for the cross-endpoint comparison
DFS62 <- setNames(CL$GSE62254$dfs_m, CL$GSE62254$sample); OS26 <- setNames(c26$os_m_smc, c26$sample)



LOCK <- c("GSE15459", "GSE34942", "GSE26253", "GSE13861", "GSE26901", "GSE26899", "GSE29272")
SC <- file.path(O, "E30_dedup_runs/primary/scores")
S15 <- readRDS(file.path(SC, "E15_scores.rds")); S20 <- readRDS(file.path(SC, "E20_scores.rds"))$OS
SURV <- c(S15[c("GSE15459", "GSE34942", "GSE26253")], S20[c("GSE13861", "GSE26901", "GSE26899", "GSE29272")])
MOD <- list(APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"), IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
            TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"), T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"),
            CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
            TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"))
L7 <- Reduce(intersect, lapply(LOCK, function(g) rownames(M[[g]])[rowSds(M[[g]]) > 0]))
bio_all <- readRDS(file.path(O, "E8_biology_null_draws.rds"))$bio
bio <- setdiff(intersect(bio_all, L7), unlist(MOD))
cat(sprintf("genes measured in all seven locked cohorts: %d; E8b biology pool %d -> %d here\n", length(L7), length(bio_all), length(bio)))
dec <- rowMeans(sapply(LOCK, function(g) { mu <- rowMeans(M[[g]][L7, ], na.rm = TRUE); ceiling(10 * rank(mu) / length(mu)) }))
dec <- round(dec); names(dec) <- L7
mod_in <- lapply(MOD, function(g) intersect(g, L7))
cat("module genes measured in all seven:", paste(names(mod_in), sapply(mod_in, length), sep = " ", collapse = "; "), "\n")
prep <- lapply(LOCK, function(g) {
  s <- SURV[[g]]; idx <- match(s$sample, colnames(M[[g]]))
  stopifnot(!anyNA(idx))
  z <- zrow(M[[g]][c(bio, unlist(mod_in)), idx])
  list(z = z, time = s$time, event = s$event, net = s$NetImmune)
}); names(prep) <- LOCK
for (g in LOCK) cat(sprintf("  %-9s n = %d, events = %d\n", g, length(prep[[g]]$time), sum(prep[[g]]$event)))
fit <- function(p, x) { f <- coxph(Surv(p$time, p$event) ~ scale(x)[, 1]); c(b = unname(coef(f)), se = sqrt(vcov(f)[1, 1])) }
axis_of <- function(z, mods) { s <- sapply(mods, function(g) colMeans(z[g, , drop = FALSE])); rowMeans(s[, 1:4]) - rowMeans(s[, 5:6]) }
obs <- sapply(LOCK, function(g) fit(prep[[g]], prep[[g]]$net))
## check: observed per-cohort NetImmune log HR equals the current published value
pub <- read.csv(file.path(O, "E30_dedup_runs/primary/E20_estimates_used_for_pooling.csv"))
pub <- pub[pub$model == "unadjusted" & grepl("primary", pub$duplicate_rule) & pub$score == "NetImmune" &
           pub$group %in% c("locked (first 3)", "locked (added)"), ]
cat("max |observed - published| NetImmune log HR:", max(abs(obs["b", LOCK] - pub$beta[match(LOCK, pub$cohort)])), "\n")
stopifnot(max(abs(obs["b", LOCK] - pub$beta[match(LOCK, pub$cohort)])) < 1e-6)
## check: the module axis rebuilt from the loaded expression (module genes measured in all seven) tracks the published score
cat("Spearman rho, rebuilt axis vs published NetImmune:", paste(LOCK, round(sapply(LOCK, function(g)
    cor(axis_of(prep[[g]]$z, mod_in), prep[[g]]$net, method = "spearman")), 3), collapse = "; "), "\n")
bio_by_dec <- split(bio, dec[bio])
draw <- function() {
  used <- character(0)
  lapply(mod_in, function(g) sapply(g, function(x) {
    cand <- setdiff(bio_by_dec[[as.character(dec[x])]], used)
    if (!length(cand)) cand <- setdiff(bio, used)
    y <- sample(cand, 1); used <<- c(used, y); y }))
}
rnd <- array(NA, c(NRAND, length(LOCK), 2), dimnames = list(NULL, LOCK, c("b", "se")))
for (k in seq_len(NRAND)) { mods <- draw(); for (g in LOCK) rnd[k, g, ] <- fit(prep[[g]], axis_of(prep[[g]]$z, mods)) }
fe <- function(b, se) { w <- 1 / se^2; sum(w * b) / sum(w) }
obs_pool <- fe(obs["b", ], obs["se", ]); rnd_pool <- apply(rnd, 1, function(m) fe(m[, "b"], m[, "se"]))
out <- rbind(
  do.call(rbind, lapply(LOCK, function(g) data.frame(cohort = g, HR_obs = exp(obs["b", g]),
    rank = sum(abs(rnd[, g, "b"]) >= abs(obs["b", g])) + 1, p_emp = (sum(abs(rnd[, g, "b"]) >= abs(obs["b", g])) + 1) / (NRAND + 1),
    null_absHR_p95 = exp(quantile(abs(rnd[, g, "b"]), 0.95))))),
  data.frame(cohort = "pooled, seven locked", HR_obs = exp(obs_pool), rank = sum(abs(rnd_pool) >= abs(obs_pool)) + 1,
             p_emp = (sum(abs(rnd_pool) >= abs(obs_pool)) + 1) / (NRAND + 1), null_absHR_p95 = exp(quantile(abs(rnd_pool), 0.95))))
rownames(out) <- NULL
print(out, digits = 4)
cat("EXPECTATION (Q7: pooled empirical P <= 0.05):", if (out$p_emp[nrow(out)] <= 0.05) "HELD" else "NOT HELD", "\n")
write.csv(out, file.path(O, "E32_locked_bionull.csv"), row.names = FALSE)
saveRDS(list(obs = obs, rnd = rnd, rnd_pool = rnd_pool, bio = bio, L7 = L7), file.path(O, "E32_locked_bionull_draws.rds"))
sink()
