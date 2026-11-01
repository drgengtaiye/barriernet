## E31 - the E30 duplicate rule applied to every pair of the ten gastric survival cohorts.
## Protocol: analysis/E31_prereg.md (sha256 in E31_prereg.sha256). Expression and clinical data only; no outcome model.
## Cohort loading as in E15/E20 (benchmark: revision/results/B1_*.rds; first three locked: E15 prep; extension: E20 prep),
## taken from the read-only check of the second review (scratchpad rev/agentB/dup_screen_all10.R), unchanged except
## that GSE26253 also takes death status from the SMC sheet. dup_pair() is E20's, verbatim, plus a row z column.
suppressMessages({library(matrixStats); library(readxl)})
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
D <- file.path(G, "2/data_download/geo")
V <- file.path(G, "2/data_download/validation/e20_outcomes")
R <- file.path(G, "revision/results")
O <- file.path(G, "2/analysis/results")
sink(file.path(O, "E31.log"), split = TRUE)
cat("E31 run", format(Sys.time()), "\n")
h <- system2("sha256sum", file.path(G, "2/analysis/E31_prereg.md"), stdout = TRUE)
stopifnot(sub(" .*", "", h) == sub(" .*", "", readLines(file.path(G, "2/analysis/E31_prereg.sha256"))[1]))
cat("protocol hash verified\n")

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


## ---- dup_pair (E20 verbatim + row z) ----
dup_pair <- function(A, B, a, b, top = NULL) {
  g <- intersect(rownames(A), rownames(B))
  g <- g[rowSums(is.na(A[g, , drop = FALSE])) == 0 & rowSums(is.na(B[g, , drop = FALSE])) == 0]
  g <- g[rowSds(A[g, , drop = FALSE]) > 0 & rowSds(B[g, , drop = FALSE]) > 0]
  if (!is.null(top))
    g <- g[order(rank(-rowVars(A[g, ])) + rank(-rowVars(B[g, ])))][1:min(top, length(g))]
  ZA <- zrow(A[g, ]); ZB <- zrow(B[g, ])
  W <- cor(ZA); w <- W[upper.tri(W)]; thr <- unname(quantile(w, 0.999))
  WB <- cor(ZB); thrB <- unname(quantile(WB[upper.tri(WB)], 0.999))
  C <- cor(ZA, ZB)
  best <- apply(C, 1, which.max); r1 <- apply(C, 1, max)
  r2 <- apply(C, 1, function(r) sort(r, decreasing = TRUE)[2])
  recip <- sapply(seq_along(best), function(i) which.max(C[, best[i]]) == i)
  zbest <- (r1 - rowMeans(C)) / rowSds(C)           # extra: how far the best match stands out of its row
  smp <- data.frame(cohort = a, sample = colnames(A), comparator = b, genes = length(g),
                    max_r = r1, second_r = r2, gap = r1 - r2, best_match = colnames(B)[best], reciprocal = recip,
                    best_match_n_partners_over_thr = sapply(best, function(j) sum(C[, j] > thr)),
                    within_p999 = thr, flagged = r1 > thr, z_best = zbest, stringsAsFactors = FALSE)
  list(smp = smp, pair = data.frame(cohort = a, comparator = b, genes = length(g),
         n_samples = ncol(A), n_comparator = ncol(B), within_p999 = thr, within_max = max(w),
         comparator_within_p999 = thrB, cross_median = median(C), cross_p999 = unname(quantile(C, 0.999)),
         cross_max = max(C), n_flagged = sum(r1 > thr), flagged = paste(colnames(A)[r1 > thr], collapse = " "),
         max_gap = max(r1 - r2), max_z_best = max(zbest), n_recip_gap_gt_0.10 = sum(recip & (r1 - r2) > 0.10),
         stringsAsFactors = FALSE))
}

COH <- c("TCGA-STAD", "GSE84437", "GSE62254", "GSE15459", "GSE34942", "GSE26253", "GSE13861", "GSE26901", "GSE26899", "GSE29272")
for (g in COH) cat(sprintf("%-10s %6d genes x %4d arrays\n", g, nrow(M[[g]]), ncol(M[[g]])))
## overall survival for rule iii: GSE26253 from the SMC sheet (time and death); every other cohort its OS as used
OSX <- CLX[, c("sample", "age", "sex", "time_m", "event", "endpoint", "cohort")]
i26 <- OSX$cohort == "GSE26253"
OSX$time_m[i26] <- c26$os_m_smc[match(OSX$sample[i26], c26$sample)]
OSX$event[i26] <- c26$death_smc[match(OSX$sample[i26], c26$sample)]
OSX$endpoint[i26] <- "OS (SMC sheet)"
TOL <- 0.5
no_conflict <- function(t1, e1, t2, e2) {
  ifelse(is.na(t1) | is.na(t2) | is.na(e1) | is.na(e2), FALSE,
  ifelse(e1 == 1 & e2 == 1, abs(t1 - t2) <= TOL,
  ifelse(e1 == 1 & e2 == 0, t2 <= t1 + TOL,
  ifelse(e1 == 0 & e2 == 1, t1 <= t2 + TOL, TRUE))))
}
PAIRS <- list(); CAND <- list()
for (ia in seq_along(COH)) for (ib in seq_along(COH)) {
  if (ia == ib) next
  a <- COH[ia]; b <- COH[ib]
  x <- dup_pair(M[[a]], M[[b]], a, b)
  S <- x$smp; m <- S[S$reciprocal, ]
  m$age <- OSX$age[match(m$sample, OSX$sample)]; m$m_age <- OSX$age[match(m$best_match, OSX$sample)]
  m$sex <- OSX$sex[match(m$sample, OSX$sample)]; m$m_sex <- OSX$sex[match(m$best_match, OSX$sample)]
  m$os <- OSX$time_m[match(m$sample, OSX$sample)]; m$m_os <- OSX$time_m[match(m$best_match, OSX$sample)]
  m$death <- OSX$event[match(m$sample, OSX$sample)]; m$m_death <- OSX$event[match(m$best_match, OSX$sample)]
  m$rule_ii <- !is.na(m$age) & !is.na(m$m_age) & !is.na(m$sex) & !is.na(m$m_sex) & m$sex == m$m_sex & abs(m$age - m$m_age) <= 1
  m$rule_iii <- no_conflict(m$os, m$death, m$m_os, m$m_death)
  m$duplicate <- m$rule_ii & m$rule_iii
  PAIRS[[length(PAIRS) + 1]] <- cbind(x$pair, n_mutual_best = nrow(m), n_rule_ii = sum(m$rule_ii),
                                      n_duplicate = sum(m$duplicate),
                                      duplicates = paste(paste(m$sample, m$best_match, sep = "=")[m$duplicate], collapse = " "),
                                      unordered_primary = ia < ib)
  if (ia < ib && any(m$rule_ii)) CAND[[length(CAND) + 1]] <- m[m$rule_ii, ]
  cat(sprintf("%-10s vs %-10s genes %5d | E20 flagged %3d | mutual best %3d | +sex/age %2d | +OS (duplicate) %d\n",
              a, b, x$pair$genes, x$pair$n_flagged, nrow(m), sum(m$rule_ii), sum(m$duplicate)))
}
P <- do.call(rbind, PAIRS); C <- do.call(rbind, CAND)
write.csv(P, file.path(O, "E31_pairs.csv"), row.names = FALSE)
write.csv(C, file.path(O, "E31_candidates.csv"), row.names = FALSE)
U <- P[P$unordered_primary, ]
other <- U[!(U$cohort == "GSE62254" & U$comparator == "GSE26253") & !(U$cohort == "GSE26253" & U$comparator == "GSE62254"), ]
cat("\nunordered pairs screened:", nrow(U), "\n")
cat("mutual best matches with sex and age agreement, outside GSE26253-GSE62254: max per pair",
    max(other$n_rule_ii), "; total", sum(other$n_rule_ii), "\n")
cat("duplicates (rules i-iii) outside GSE26253-GSE62254:", sum(other$n_duplicate), "\n")
chk <- U[(U$cohort == "GSE62254" & U$comparator == "GSE26253") | (U$cohort == "GSE26253" & U$comparator == "GSE62254"), ]
E30 <- c("GSM644726", "GSM644801", "GSM644802", "GSM644593", "GSM644733", "GSM644810", "GSM644871", "GSM644940")
found <- unlist(regmatches(chk$duplicates, gregexpr("GSM644[0-9]{3}", chk$duplicates)))
cat("GSE26253-GSE62254 check (OS-only rule iii):", chk$n_duplicate, "pairs;", "E30 eight recovered:", all(E30 %in% found),
    "; extra:", paste(setdiff(found, E30), collapse = " "), "\n")
cat("EXPECTATION (no other pair meets i-iii):", if (sum(other$n_duplicate) == 0) "HELD" else "NOT HELD", "\n")
options(width = 250)
cat("\nall mutual best matches with sex and age agreement:\n")
print(transform(C[, c("cohort", "sample", "comparator", "best_match", "max_r", "gap", "z_best", "sex", "m_sex", "age", "m_age",
                      "os", "death", "m_os", "m_death", "rule_iii", "duplicate")],
                max_r = round(max_r, 3), gap = round(gap, 3), z_best = round(z_best, 2), os = round(os, 1), m_os = round(m_os, 1)),
      row.names = FALSE)
sink()
