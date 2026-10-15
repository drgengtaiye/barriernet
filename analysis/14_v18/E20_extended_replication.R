## E20 - EXTENDED LOCKED gastric prognostic replication (pre-registered in E17_E22_prereg.md, section E20,
## and the dated "2026-10-04: E20 data" entry under Deviations; that file is not edited here).
##
## Status (from the pre-registration): this extension was decided after the three-cohort locked
## replication (E15) gave P = 0.057 and is reported as such.
##
## PRE-SPECIFIED EXPECTATIONS (E20):
##   X1  NetImmune HR per SD < 1 in each added cohort (OS primary; RFS only where it is the only outcome).
##   X2  CAF_ECM  HR per SD > 1 in each added cohort.
##   X3  PRIMARY: NetImmune pooled across all 7 locked cohorts (GSE15459, GSE34942, GSE26253 from E15 +
##       GSE13861, GSE26901, GSE29272, GSE26899): HR < 1.
##   Also reported: all gastric cohorts pooled (+ TCGA-STAD, GSE84437, GSE62254); structure-matched
##   random-gene null as in E15; Delta C of NetImmune vs ESTIMATE and CD8A.
##   Sensitivity: + GSE57303 (third-party OS re-curation); OS-only sets.
##
## LOCKED METHODS - copied from E15_gastric_validation.R without change:
##   - scale check per cohort: max > 50 -> linear, log2(x+1) applied once; otherwise used as provided;
##   - probes -> gene symbol by highest mean expression (GPL570: GSE62254 annotation, first symbol of
##     '///'; GPL96: GEO platform table 'Gene Symbol', first symbol of '///'; GPL6884/GPL6947: GEO
##     platform table 'Symbol' (GPL6947 table identical to the GSE84437 fdata used in B1));
##   - genes row-z within cohort; module = mean z of its genes present (>= 3); Barrier = mean(CAF_ECM,
##     TAM_SUP); Immune = mean(APM, IFN, TLS_B, T_CELL); NetImmune = Immune - Barrier;
##   - Cox HR per SD; adjusted = + stage (ordinal I-IV) + age where recorded (> 20 non-missing);
##   - random-effects pooling (metafor), wider of HK / DL interval, larger P;
##   - structure-matched random-gene null (1,000 axes, 8/5/6/7 minus 11/11 genes, same gene set in every
##     cohort, all genes measured in every cohort of the set as pool, better direction credited once on the
##     inverse-variance pooled log HR);
##   - Harrell's C with direction fixed; paired bootstrap (2,000) Delta C, pooled by inverse variance.
## Previously analysed cohorts are NOT re-estimated for pooling: GSE15459/GSE34942/GSE26253 come from
## results/E15_per_cohort.csv and the benchmark cohorts from results/E8_univariate.csv (+ E5 for
## CXCL9_SPP1), exactly as E15 did. (E15's three cohorts are re-loaded only for the duplicate check and the
## random-gene null, and the script checks that their observed estimates reproduce E15_per_cohort.csv.)
##
## STEP ORDER: the expression duplicate check (step 1) runs before any outcome model. Run with
## E20_STAGE=dup to stop after it.
suppressMessages({library(survival); library(metafor); library(matrixStats); library(readxl)})
set.seed(2026)
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
D <- file.path(G, "2/data_download/geo")
V <- file.path(G, "2/data_download/validation/e20_outcomes")
O <- file.path(G, "2/analysis/results")
R <- file.path(G, "revision/results")
STAGE <- Sys.getenv("E20_STAGE", "all")
LOG <- file.path(O, if (STAGE == "dup") "E20_log_dupcheck.txt" else "E20_log.txt"); sink(LOG, split = TRUE)
cat("E20 run", format(Sys.time()), "stage:", STAGE, "\n")

## ---------------- locked definitions and helpers (verbatim from E15) ----------------
MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"))
SCORES <- c("NetImmune","CAF_ECM","TAM_SUP","Barrier","Immune","APM","IFN","TLS_B","T_CELL",
            "CXCL9_SPP1","ESTIMATE","CD8A")

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
  meta <- L[1:b]
  list(m = m, meta = meta)
}
meta_field <- function(meta, key) {
  l <- meta[grepl(paste0("^!", key, "\t"), meta)]
  lapply(l, function(x) gsub('"', "", strsplit(x, "\t")[[1]][-1]))
}
## E15's prep(), with one addition for series that also deposit non-tumour arrays (normal mucosa, GIST):
## `keep` restricts the matrix to the cohort's gastric-cancer arrays before probes are collapsed, as the
## benchmark pipeline did for TCGA (primary tumours selected before collapsing). For E15's cohorts keep = NULL.
SCALE <- list()
prep <- function(g, symmap, keep = NULL) {
  x <- read_matrix(g); m <- x$m
  rng <- range(m, na.rm = TRUE)
  scale_note <- if (rng[2] > 50) "linear (max > 50): log2(x+1) applied once" else "log2 as deposited: used as provided"
  if (rng[2] > 50) m <- log2(pmax(m, 0) + 1)
  n_all <- ncol(m)
  if (!is.null(keep)) m <- m[, intersect(keep, colnames(m)), drop = FALSE]
  m <- m[rowSums(is.na(m)) == 0, , drop = FALSE]
  q <- quantile(m, c(0.01, 0.5, 0.99))
  m <- collapse_max(m, unname(symmap[rownames(m)]))
  cat(sprintf("%-9s stored range %.2f..%.2f (1%%/50%%/99%% %.2f/%.2f/%.2f) -> %s; %d genes x %d of %d arrays\n",
              g, rng[1], rng[2], q[1], q[2], q[3], scale_note, nrow(m), ncol(m), n_all))
  SCALE[[g]] <<- data.frame(cohort = g, min = rng[1], q01 = q[1], median = q[2], q99 = q[3], max = rng[2],
                            scale = scale_note, arrays_deposited = n_all, arrays_used = ncol(m), genes = nrow(m))
  list(m = m, meta = x$meta, scale = scale_note, range = rng)
}

## ---------------- annotations ----------------
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

## ---------------- E15 locked cohorts (loaded exactly as in E15) ----------------
cat("== loading E15 locked cohorts (duplicate check and null only)\n")
X15 <- list()
out_sg <- function(g) {
  o <- as.data.frame(read_excel(file.path(D, g, paste0(g, "_outcome.xls"))))
  names(o) <- c("gsm","id","cel","subtype","age","gender","lauren","stage","os_m","dead")
  o <- o[!is.na(o$gsm) & grepl("^GSM", o$gsm), ]
  data.frame(sample = o$gsm, pid = as.character(o$id), time = as.numeric(o$os_m),
             event = as.integer(o$dead), age = as.numeric(o$age),
             stage_ord = suppressWarnings(as.integer(o$stage)), endpoint = "OS",
             stringsAsFactors = FALSE)
}
for (g in c("GSE15459", "GSE34942")) {
  p <- prep(g, sym570); cl <- out_sg(g)
  cl <- cl[cl$sample %in% colnames(p$m) & !is.na(cl$time) & cl$time > 0 & !is.na(cl$event), ]
  X15[[g]] <- list(m = p$m[, cl$sample], clin = cl, platform = "GPL570 (U133 Plus 2.0)", scale = p$scale)
}
p <- prep("GSE26253", sym8432)
st <- meta_field(p$meta, "Sample_characteristics_ch1")
key <- function(k) { i <- which(sapply(st, function(v) any(grepl(k, v, fixed = TRUE)))); sub(".*: ", "", st[[i[1]]]) }
gsm <- meta_field(p$meta, "Sample_geo_accession")[[1]]
stg <- key("pathological stage")
cl <- data.frame(sample = gsm, pid = gsm, time = as.numeric(key("recurrence free survival time")),
                 event = as.integer(key("status (0=non-recurrence")), age = NA_real_,
                 stage_ord = c("1"=1, "1A"=1, "1B"=1, "2"=2, "3A"=3, "3B"=3, "4"=4)[stg],
                 endpoint = "RFS", stringsAsFactors = FALSE)
cl <- cl[!is.na(cl$time) & cl$time > 0 & !is.na(cl$event), ]
X15[["GSE26253"]] <- list(m = p$m[, cl$sample], clin = cl, platform = "GPL8432 (DASL, FFPE)", scale = p$scale)

## ---------------- added cohorts (E20) ----------------
cat("\n== loading added cohorts\n")
NEW <- c("GSE13861", "GSE26901", "GSE26899", "GSE29272")
SENS <- "GSE57303"
PLAT <- c(GSE13861 = "GPL6884 (Illumina WG-6 v3)", GSE26901 = "GPL6947 (Illumina HT-12 v3)",
          GSE26899 = "GPL6947 (Illumina HT-12 v3)", GSE29272 = "GPL96 (U133A)", GSE57303 = "GPL570 (U133 Plus 2.0)")
SYM <- list(GSE13861 = sym6884, GSE26901 = sym6947, GSE26899 = sym6947, GSE29272 = sym96, GSE57303 = sym570)
OUTF <- c(GSE13861 = "GSE13861_outcome.csv", GSE26901 = "GSE26901_outcome.csv", GSE26899 = "GSE26899_outcome.csv",
          GSE29272 = "GSE29272_outcome.csv", GSE57303 = "GSE57303_outcome_secondary.csv")
RAW <- list(); OC <- list()
for (g in c(NEW, SENS)) {
  oc <- read.csv(file.path(V, OUTF[[g]]), stringsAsFactors = FALSE)
  OC[[g]] <- oc
  keep <- if (g == "GSE57303") NULL else oc$gsm        # GSE57303: all 70 deposited arrays are tumours
  RAW[[g]] <- prep(g, SYM[[g]], keep = keep)
}
## clinical tables: OS (primary) and RFS (secondary, where recorded); time > 0 rule as in E15
mk_clin <- function(g, endpoint) {
  oc <- OC[[g]]
  tt <- if (endpoint == "OS") oc$os_time else oc$rfs_time
  ee <- if (endpoint == "OS") oc$os_event else oc$rfs_event
  stage <- if ("stage_num" %in% names(oc)) oc$stage_num else if ("stage" %in% names(oc)) oc$stage else NA
  data.frame(sample = oc$gsm, pid = if ("array_id" %in% names(oc)) oc$array_id else oc$gsm,
             time = as.numeric(tt), event = as.integer(ee),
             age = if ("age" %in% names(oc)) as.numeric(oc$age) else NA_real_,
             stage_ord = suppressWarnings(as.integer(stage)), endpoint = endpoint, stringsAsFactors = FALSE)
}

## =====================================================================================
## STEP 1. Expression duplicate check (before any outcome model)
## =====================================================================================
cat("\n== STEP 1: duplicate check\n")
B84 <- readRDS(file.path(R, "B1_GSE84437.rds")); B62 <- readRDS(file.path(R, "B1_GSE62254.rds"))
BTC <- readRDS(file.path(R, "B1_TCGA.rds"))
COMP <- list(GSE84437 = B84$expr, GSE62254 = B62$expr, GSE26253 = X15$GSE26253$m,
             GSE15459 = X15$GSE15459$m, GSE34942 = X15$GSE34942$m, `TCGA-STAD` = BTC$expr)
## new-cohort arrays checked: every gastric-cancer array in the outcome table (with or without outcome)
NEWM <- lapply(c(NEW, SENS), function(g) RAW[[g]]$m); names(NEWM) <- c(NEW, SENS)
dup_pair <- function(A, B, a, b, top = NULL) {
  g <- intersect(rownames(A), rownames(B))
  g <- g[rowSums(is.na(A[g, , drop = FALSE])) == 0 & rowSums(is.na(B[g, , drop = FALSE])) == 0]  # GSE84437: 35 genes with NA
  g <- g[rowSds(A[g, , drop = FALSE]) > 0 & rowSds(B[g, , drop = FALSE]) > 0]
  if (!is.null(top))   # most variable in both cohorts (sum of within-cohort variance ranks; platform-free)
    g <- g[order(rank(-rowVars(A[g, ])) + rank(-rowVars(B[g, ])))][1:min(top, length(g))]
  ## genes are z-scored within each cohort (as in E15's overlap check), so platform and batch offsets cancel
  ZA <- zrow(A[g, ]); ZB <- zrow(B[g, ])
  W <- cor(ZA); w <- W[upper.tri(W)]; thr <- unname(quantile(w, 0.999))   # within-cohort pairs, same genes
  WB <- cor(ZB); thrB <- unname(quantile(WB[upper.tri(WB)], 0.999))
  C <- cor(ZA, ZB)                                                          # Pearson
  best <- apply(C, 1, which.max); r1 <- apply(C, 1, max)
  r2 <- apply(C, 1, function(r) sort(r, decreasing = TRUE)[2])
  recip <- sapply(seq_along(best), function(i) which.max(C[, best[i]]) == i)  # a is also b's best match
  smp <- data.frame(cohort = a, sample = colnames(A), comparator = b, genes = length(g),
                    max_r = r1, second_r = r2, gap = r1 - r2, best_match = colnames(B)[best], reciprocal = recip,
                    best_match_n_partners_over_thr = sapply(best, function(j) sum(C[, j] > thr)),
                    within_p999 = thr, flagged = r1 > thr, stringsAsFactors = FALSE)
  list(smp = smp, pair = data.frame(cohort = a, comparator = b, genes = length(g),
         n_samples = ncol(A), n_comparator = ncol(B), within_p999 = thr, within_max = max(w),
         comparator_within_p999 = thrB, cross_median = median(C), cross_p999 = unname(quantile(C, 0.999)),
         cross_max = max(C), n_flagged = sum(r1 > thr), flagged = paste(colnames(A)[r1 > thr], collapse = " "),
         max_gap = max(r1 - r2), stringsAsFactors = FALSE))
}
## comparators used for exclusion: the five pre-registered used cohorts, TCGA-STAD, and the sister series
## GSE26901 <-> GSE26899. Other pairs of added cohorts are computed as diagnostics only.
EXPAIR <- function(a, b) b %in% names(COMP) || (a == "GSE26901" && b == "GSE26899") || (a == "GSE26899" && b == "GSE26901")
DS <- list(); DP <- list(); DS2 <- list(); DP2 <- list()
for (a in names(NEWM)) {
  targets <- c(COMP, NEWM[setdiff(names(NEWM), a)])
  for (b in names(targets)) {
    x <- dup_pair(NEWM[[a]], targets[[b]], a, b); DS[[paste(a, b)]] <- x$smp; DP[[paste(a, b)]] <- x$pair
    y <- dup_pair(NEWM[[a]], targets[[b]], a, b, top = 2000)
    DS2[[paste(a, b)]] <- y$smp; DP2[[paste(a, b)]] <- y$pair
  }
}
DS <- do.call(rbind, DS); DP <- do.call(rbind, DP); DS2 <- do.call(rbind, DS2); DP2 <- do.call(rbind, DP2)
DS$used_for_exclusion <- mapply(EXPAIR, DS$cohort, DS$comparator); DP$used_for_exclusion <- mapply(EXPAIR, DP$cohort, DP$comparator)
DS2$used_for_exclusion <- FALSE; DP2$used_for_exclusion <- FALSE
DP$genes_set <- "all shared genes (pre-registered rule)"; DP2$genes_set <- "2,000 most variable shared genes (diagnostic)"
DUP <- rbind(DP, DP2)
DUP$comparator_type <- ifelse(DUP$comparator %in% names(COMP), "used cohort", "added/sensitivity cohort")
write.csv(DUP, file.path(O, "E20_duplicate_check.csv"), row.names = FALSE)
DSall <- rbind(cbind(DS, genes_set = "all shared"), cbind(DS2, genes_set = "top 2000 variable"))
write.csv(DSall, file.path(O, "E20_duplicate_check_samples.csv"), row.names = FALSE)
cat("pre-registered rule (all shared genes, z-scored within cohort; flag if max cross-cohort Pearson r > within-cohort 99.9th pct):\n")
print(transform(DP[, c("cohort","comparator","used_for_exclusion","genes","n_samples","within_p999","cross_p999","cross_max","n_flagged","max_gap")],
                within_p999 = round(within_p999, 3), cross_p999 = round(cross_p999, 3),
                cross_max = round(cross_max, 3), max_gap = round(max_gap, 3)), row.names = FALSE)
cat("\ndiagnostic (2,000 most variable shared genes):\n")
print(transform(DP2[, c("cohort","comparator","genes","within_p999","cross_max","n_flagged","max_gap")],
                within_p999 = round(within_p999, 3), cross_max = round(cross_max, 3), max_gap = round(max_gap, 3)),
      row.names = FALSE)

## clinical concordance of flagged pairs (age, sex where both cohorts record them)
p84 <- read.csv(file.path(G, "GC_TME_Reprogrammer/data_GEO/GSE84437/GSE84437_pheno.csv"), stringsAsFactors = FALSE)
CLX <- rbind(
  data.frame(sample = p84$geo_accession, age = as.numeric(p84[["age.ch1"]]), sex = tolower(p84[["Sex.ch1"]])),
  do.call(rbind, lapply(c("GSE15459", "GSE34942"), function(g) { o <- as.data.frame(read_excel(file.path(D, g, paste0(g, "_outcome.xls"))))
    names(o) <- c("gsm","id","cel","subtype","age","gender","lauren","stage","os_m","dead"); o <- o[!is.na(o$gsm) & grepl("^GSM", o$gsm), ]
    data.frame(sample = o$gsm, age = as.numeric(o$age), sex = tolower(substr(o$gender, 1, 1))) })),
  data.frame(sample = B62$clin$sample, age = B62$clin$age, sex = NA), data.frame(sample = BTC$clin$sample, age = BTC$clin$age, sex = NA),
  do.call(rbind, lapply(c(NEW, SENS), function(g) data.frame(sample = OC[[g]]$gsm,
    age = if ("age" %in% names(OC[[g]])) as.numeric(OC[[g]]$age) else NA,
    sex = if ("sex" %in% names(OC[[g]])) tolower(substr(OC[[g]]$sex, 1, 1)) else NA))))
CLX$sex <- substr(CLX$sex, 1, 1)
FL <- DS[DS$flagged, ]
FL$n_comparators_flagged <- sapply(FL$sample, function(s) sum(DS$flagged[DS$sample == s]))
FL$age <- CLX$age[match(FL$sample, CLX$sample)]; FL$sex <- CLX$sex[match(FL$sample, CLX$sample)]
FL$match_age <- CLX$age[match(FL$best_match, CLX$sample)]; FL$match_sex <- CLX$sex[match(FL$best_match, CLX$sample)]
FL$age_within_1y <- abs(FL$age - FL$match_age) <= 1
FL$same_sex <- ifelse(is.na(FL$sex) | is.na(FL$match_sex), NA, FL$sex == FL$match_sex)
write.csv(FL, file.path(O, "E20_duplicate_flagged.csv"), row.names = FALSE)
cat("\nflagged sample-comparator pairs (all pairs; used_for_exclusion marks the pre-registered/requested comparators):\n")
print(transform(FL[, c("cohort","sample","comparator","used_for_exclusion","best_match","max_r","within_p999","gap","reciprocal",
                       "best_match_n_partners_over_thr","n_comparators_flagged","age","match_age","sex","match_sex")],
                max_r = round(max_r, 3), within_p999 = round(within_p999, 3), gap = round(gap, 3)), row.names = FALSE)
cat("best matches shared by several flagged samples (hub arrays):\n"); print(sort(table(FL$best_match), decreasing = TRUE)[1:6])

## exclusions (pre-registered rule, applied literally): every added-cohort sample flagged against an exclusion comparator
FE <- FL[FL$used_for_exclusion, ]
EXCL <- do.call(rbind, lapply(split(FE, FE$sample), function(d) data.frame(cohort = d$cohort[1], sample = d$sample[1],
          reason = paste(sprintf("%s %s r=%.3f (thr %.3f)", d$comparator, d$best_match, d$max_r, d$within_p999), collapse = "; "))))
if (is.null(EXCL)) EXCL <- data.frame(cohort = character(), sample = character(), reason = character())
cat("\nsamples excluded under the pre-registered rule:", nrow(EXCL), "\n"); print(table(EXCL$cohort))
write.csv(EXCL, file.path(O, "E20_duplicate_exclusions.csv"), row.names = FALSE)
if (STAGE == "dup") { sink(); quit(save = "no") }

## =====================================================================================
## STEP 2. Scores and per-cohort Cox models
## =====================================================================================
build <- function(g, endpoint, excl = EXCL$sample) {
  cl <- mk_clin(g, endpoint)
  cl <- cl[cl$sample %in% colnames(RAW[[g]]$m) & !is.na(cl$time) & cl$time > 0 & !is.na(cl$event), ]
  cl <- cl[!cl$sample %in% excl, ]
  list(m = RAW[[g]]$m[, cl$sample], clin = cl, platform = PLAT[[g]], scale = RAW[[g]]$scale)
}
XO <- lapply(c(NEW, SENS), build, endpoint = "OS"); names(XO) <- c(NEW, SENS)
XR <- lapply(c("GSE13861", "GSE26901", "GSE26899"), build, endpoint = "RFS"); names(XR) <- c("GSE13861", "GSE26901", "GSE26899")
## sensitivity fixed after the duplicate check and before any outcome model: no duplicate exclusions (the flagged
## arrays match tumours from several countries at once and none is age/sex-concordant with its match)
XOn <- lapply(c(NEW, SENS), build, endpoint = "OS", excl = character()); names(XOn) <- c(NEW, SENS)
XRn <- lapply(c("GSE13861", "GSE26901", "GSE26899"), build, endpoint = "RFS", excl = character()); names(XRn) <- names(XR)

data(SI_geneset, package = "estimate", envir = environment())
imm <- as.character(unlist(SI_geneset["ImmuneSignature", -1])); imm <- imm[imm != "" & !is.na(imm)]
str_ <- as.character(unlist(SI_geneset["StromalSignature", -1])); str_ <- str_[str_ != "" & !is.na(str_)]
COV <- list()
score <- function(x, tag) {
  Z <- zrow(x$m); cl <- x$clin
  for (k in names(MOD)) {
    g <- intersect(MOD[[k]], rownames(Z))
    COV[[paste(tag, k)]] <<- data.frame(cohort = tag, module = k, n_genes = length(g), of = length(MOD[[k]]),
                                         missing = paste(setdiff(MOD[[k]], g), collapse = " "))
    cl[[k]] <- if (length(g) >= 3) colMeans(Z[g, , drop = FALSE]) else NA
  }
  cl$Barrier <- rowMeans(cl[, c("CAF_ECM", "TAM_SUP")])
  cl$Immune <- rowMeans(cl[, c("APM", "IFN", "TLS_B", "T_CELL")])
  cl$NetImmune <- cl$Immune - cl$Barrier
  cl$CXCL9_SPP1 <- if (all(c("CXCL9", "SPP1") %in% rownames(Z))) Z["CXCL9", ] - Z["SPP1", ] else NA
  cl$ESTIMATE <- colMeans(Z[intersect(imm, rownames(Z)), ]) - colMeans(Z[intersect(str_, rownames(Z)), ])
  cl$CD8A <- if ("CD8A" %in% rownames(Z)) Z["CD8A", ] else NA
  cl$cohort <- tag; cl
}
SO <- Map(score, XO, names(XO))
covtab <- do.call(rbind, COV); COV <- list()
SR <- Map(score, XR, names(XR)); COV <- list()
S15 <- Map(score, X15, names(X15)); COV <- list()
SOn <- Map(score, XOn, names(XOn)); SRn <- Map(score, XRn, names(XRn)); COV <- list()
write.csv(covtab, file.path(O, "E20_module_coverage.csv"), row.names = FALSE)
cat("\n== module coverage (added cohorts)\n"); print(covtab[covtab$n_genes < covtab$of, ], row.names = FALSE)

cox1 <- function(d, s, adj) {
  x <- d[[s]]; if (all(is.na(x))) return(NULL)
  d$x <- as.numeric(scale(x))
  cov <- c()
  if (adj) { if (sum(!is.na(d$stage_ord)) > 20) cov <- c(cov, "stage_ord")
             if (sum(!is.na(d$age)) > 20) cov <- c(cov, "age")
             if (!length(cov)) return(NULL) }
  f <- as.formula(paste("Surv(time, event) ~ x", if (length(cov)) paste("+", paste(cov, collapse = "+")) else ""))
  dd <- d[complete.cases(d[, c("time", "event", "x", cov)]), ]
  fit <- coxph(f, data = dd); co <- summary(fit)$coefficients["x", ]
  data.frame(cohort = d$cohort[1], score = s, model = if (adj) "adjusted" else "unadjusted",
             covariates = paste(cov, collapse = "+"), n = nrow(dd), events = sum(dd$event),
             endpoint = d$endpoint[1], beta = co[["coef"]], se = co[["se(coef)"]],
             HR = exp(co[["coef"]]), lo = exp(co[["coef"]] - 1.96 * co[["se(coef)"]]),
             hi = exp(co[["coef"]] + 1.96 * co[["se(coef)"]]), p = co[["Pr(>|z|)"]])
}
fit_all <- function(S) do.call(rbind, lapply(S, function(d) do.call(rbind, c(lapply(SCORES, cox1, d = d, adj = FALSE),
                                                                             lapply(SCORES, cox1, d = d, adj = TRUE)))))
PCO <- fit_all(SO); PCR <- fit_all(SR)
PCO$role <- ifelse(PCO$cohort == SENS, "sensitivity (secondary OS source)", "locked (added)")
PCR$role <- "locked (added), secondary endpoint"
## sensitivity for the RFS anomalies (RFS > OS in GSE13861 YG24T, GSE26899 KG028C/KG114C): drop them
anom <- c(OC$GSE13861$gsm[OC$GSE13861$rfs_time > OC$GSE13861$os_time], OC$GSE26899$gsm[which(OC$GSE26899$rfs_time > OC$GSE26899$os_time)])
SRa <- lapply(SR, function(d) d[!d$sample %in% anom, ])
PCRa <- do.call(rbind, lapply(SRa, function(d) do.call(rbind, lapply(c("NetImmune", "CAF_ECM"), cox1, d = d, adj = FALSE))))
PCRa$role <- "RFS, RFS>OS records dropped (sensitivity)"
PCOn <- fit_all(SOn); PCRn <- fit_all(SRn)
PCOn$role <- ifelse(PCOn$cohort == SENS, "sensitivity (secondary OS source)", "locked (added)")
PCRn$role <- "locked (added), secondary endpoint"
PC <- rbind(cbind(rbind(PCO, PCR, PCRa), duplicate_rule = "pre-registered exclusions (primary)"),
            cbind(rbind(PCOn, PCRn), duplicate_rule = "no exclusions (sensitivity)"))
write.csv(PC, file.path(O, "E20_per_cohort.csv"), row.names = FALSE)
cat("\n== per-cohort HR per SD (OS, unadjusted and adjusted)\n")
shw <- function(a) print(transform(a[, c("cohort","score","model","endpoint","n","events","HR","lo","hi","p")],
                                   HR = round(HR, 2), lo = round(lo, 2), hi = round(hi, 2), p = signif(p, 2)), row.names = FALSE)
shw(PCO[PCO$score %in% c("NetImmune", "CAF_ECM"), ])
cat("\n== per-cohort HR per SD (RFS, secondary)\n"); shw(PCR[PCR$score %in% c("NetImmune", "CAF_ECM"), ])
cat("\n== RFS without RFS>OS records:", paste(anom, collapse = " "), "\n"); shw(PCRa)
cat("\n== OS, no duplicate exclusions (sensitivity)\n"); shw(PCOn[PCOn$score %in% c("NetImmune", "CAF_ECM") & PCOn$model == "unadjusted", ])

## reproducibility check: E15 cohorts re-scored here must give E15's estimates
P15 <- read.csv(file.path(O, "E15_per_cohort.csv"), stringsAsFactors = FALSE)
chk <- fit_all(S15); mm <- merge(chk, P15, by = c("cohort", "score", "model"))
cat("\nE15 cohorts reproduce E15_per_cohort.csv: max |beta difference| =", signif(max(abs(mm$beta.x - mm$beta.y)), 3),
    "over", nrow(mm), "estimates\n")
stopifnot(max(abs(mm$beta.x - mm$beta.y)) < 1e-8)

## =====================================================================================
## STEP 3. Pooling (E15 pool(), unchanged)
## =====================================================================================
pool <- function(b, se) {
  dl <- rma(yi = b, sei = se, method = "DL"); hk <- rma(yi = b, sei = se, method = "DL", test = "knha")
  wide_hk <- (hk$ci.ub - hk$ci.lb) >= (dl$ci.ub - dl$ci.lb)
  data.frame(k = length(b), HR = exp(dl$b[1]),
             lo = exp(if (wide_hk) hk$ci.lb else dl$ci.lb), hi = exp(if (wide_hk) hk$ci.ub else dl$ci.ub),
             p = max(dl$pval, hk$pval), interval = if (wide_hk) "HK" else "DL",
             DL_lo = exp(dl$ci.lb), DL_hi = exp(dl$ci.ub), DL_p = dl$pval,
             HK_lo = exp(hk$ci.lb), HK_hi = exp(hk$ci.ub), HK_p = hk$pval, I2 = dl$I2,
             n_adverse = sum(b > 0))
}
orig <- read.csv(file.path(O, "E8_univariate.csv"), stringsAsFactors = FALSE)
orig <- rbind(orig, setNames(read.csv(file.path(O, "E5_gastric_arms_per_cohort.csv"))[, c("cohort","score","b","se")],
                             c("cohort","score","beta","se")))
orig <- orig[!duplicated(orig[, c("cohort", "score")]), ]
orig$endpoint <- "OS"
## per-cohort estimate table for pooling: benchmark (E8), locked first 3 (E15), added 4 + GSE57303 (here)
L3 <- c("GSE15459", "GSE34942", "GSE26253"); BM <- c("TCGA", "GSE84437", "GSE62254")
SETS <- list(
  "locked 7 (primary)"              = c(L3, NEW),
  "locked (added 4)"                = NEW,
  "locked (first 3), E15 reference" = L3,
  "all 10 gastric"                  = c(BM, L3, NEW),
  "benchmark 3 (reference)"         = BM,
  "locked 7 + GSE57303 (sensitivity)" = c(L3, NEW, SENS),
  "all 10 + GSE57303 (sensitivity)" = c(BM, L3, NEW, SENS),
  "locked 7, OS only (sensitivity)" = c("GSE15459", "GSE34942", NEW),
  "all 10, OS only (sensitivity)"   = c(BM, "GSE15459", "GSE34942", NEW),
  "locked (added 4) + GSE57303 (sensitivity)" = c(NEW, SENS))
make_est <- function(PCO, PCR) rbind(
  data.frame(cohort = orig$cohort, score = orig$score, model = "unadjusted", endpoint = "OS", beta = orig$beta,
             se = orig$se, group = "benchmark", stringsAsFactors = FALSE),
  data.frame(P15[, c("cohort", "score", "model", "endpoint", "beta", "se")], group = "locked (first 3)"),
  data.frame(PCO[PCO$cohort %in% NEW, c("cohort", "score", "model", "endpoint", "beta", "se")], group = "locked (added)"),
  data.frame(PCO[PCO$cohort == SENS, c("cohort", "score", "model", "endpoint", "beta", "se")], group = "sensitivity"),
  data.frame(PCR[, c("cohort", "score", "model", "endpoint", "beta", "se")], group = "locked (added), RFS"))
pool_sets <- function(EST) {
  PL <- list()
  for (s in SCORES) for (nm in names(SETS)) for (m in c("unadjusted", "adjusted")) {
    a <- EST[EST$score == s & EST$model == m & EST$cohort %in% SETS[[nm]] & EST$group != "locked (added), RFS", ]
    if (m == "adjusted" && any(!SETS[[nm]] %in% a$cohort)) next     # adjusted only where every cohort has it
    if (nrow(a) >= 2) PL[[length(PL) + 1]] <- cbind(set = nm, score = s, model = m,
                                                     cohorts = paste(a$cohort, collapse = " "), pool(a$beta, a$se))
  }
  ## secondary: RFS set (GSE26253 + the three added cohorts with RFS)
  for (s in SCORES) for (m in c("unadjusted", "adjusted")) {
    a <- rbind(EST[EST$score == s & EST$model == m & EST$cohort == "GSE26253", ],
               EST[EST$score == s & EST$model == m & EST$group == "locked (added), RFS", ])
    if (nrow(a) >= 2) PL[[length(PL) + 1]] <- cbind(set = "locked, RFS (secondary)", score = s, model = m,
                                                     cohorts = paste(a$cohort, collapse = " "), pool(a$beta, a$se))
  }
  do.call(rbind, PL)
}
EST <- make_est(PCO, PCR); ESTn <- make_est(PCOn, PCRn)
write.csv(rbind(cbind(EST, duplicate_rule = "pre-registered exclusions (primary)"),
                cbind(ESTn[ESTn$group %in% c("locked (added)", "sensitivity", "locked (added), RFS"), ], duplicate_rule = "no exclusions (sensitivity)")),
          file.path(O, "E20_estimates_used_for_pooling.csv"), row.names = FALSE)
PL <- pool_sets(EST); PLn <- pool_sets(ESTn)
write.csv(rbind(cbind(PL, duplicate_rule = "pre-registered exclusions (primary)"),
                cbind(PLn[!PLn$set %in% c("locked (first 3), E15 reference", "benchmark 3 (reference)"), ], duplicate_rule = "no exclusions (sensitivity)")),
          file.path(O, "E20_pooled.csv"), row.names = FALSE)
cat("\n== pooled (conservative interval)\n")
print(transform(PL[PL$score %in% c("NetImmune", "CAF_ECM"), c("set","score","model","k","HR","lo","hi","p","interval","I2","n_adverse")],
                HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 3), I2 = round(I2)), row.names = FALSE)
cat("\n== pooled NetImmune / CAF_ECM, no duplicate exclusions (sensitivity)\n")
print(transform(PLn[PLn$score %in% c("NetImmune", "CAF_ECM") & PLn$model == "unadjusted" & !grepl("reference", PLn$set), c("set","score","k","HR","lo","hi","p","interval","I2")],
                HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 3), I2 = round(I2)), row.names = FALSE)
cat("\n== pooled, all scores, locked 7 (primary) and all 10, unadjusted\n")
print(transform(PL[PL$set %in% c("locked 7 (primary)", "all 10 gastric", "locked (added 4)") & PL$model == "unadjusted",
                   c("set","score","k","HR","lo","hi","p","I2")],
                HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 3), I2 = round(I2)), row.names = FALSE)

## =====================================================================================
## STEP 4. Structure-matched random-gene null (NetImmune, unadjusted, as E15)
## =====================================================================================
cat("\n== STEP 4: random-gene null\n")
LK <- c(S15, SO[NEW]); XL <- c(X15, XO[NEW])                     # the 7 locked cohorts (added: OS)
Zs <- lapply(XL, function(x) zrow(x$m))
gpool <- Reduce(intersect, lapply(Zs, rownames))
gpool <- gpool[sapply(gpool, function(g) all(sapply(Zs, function(z) sd(z[g, ]) > 0)))]
cat("null gene pool (measured in all 7 locked cohorts):", length(gpool), "\n")
Z57 <- zrow(XO[[SENS]]$m); do57 <- all(gpool %in% rownames(Z57)) && all(rowSds(Z57[intersect(gpool, rownames(Z57)), ]) > 0)
cat("GSE57303 carries the whole pool:", do57, "\n")
if (do57) { Zs[[SENS]] <- Z57; LK[[SENS]] <- SO[[SENS]] }
SIZES <- c(8, 5, 6, 7, 11, 11)
obsb <- sapply(names(LK), function(k) { r <- EST[EST$cohort == k & EST$score == "NetImmune" & EST$model == "unadjusted" &
                                               EST$group != "locked (added), RFS", ]; c(r$beta, r$se) })
fitb <- function(d, x) { co <- summary(coxph(Surv(d$time, d$event) ~ as.numeric(scale(x))))$coefficients; co[1, c(1, 3)] }
iv <- function(b, se) { w <- 1 / se^2; sum(w * b) / sum(w) }
POOLS <- list("locked 7, pooled" = c(L3, NEW), "locked (added 4), pooled" = NEW)
if (do57) POOLS[["locked 7 + GSE57303, pooled (sensitivity)"]] <- c(L3, NEW, SENS)
NR <- 1000; K <- names(LK)
draws <- matrix(NA, NR, length(K) + length(POOLS), dimnames = list(NULL, c(K, names(POOLS))))
for (r in 1:NR) {
  gs <- split(sample(gpool, sum(SIZES)), rep(1:6, SIZES))
  bs <- ses <- setNames(numeric(length(K)), K)
  for (k in K) {
    z <- Zs[[k]]; mods <- sapply(gs, function(g) colMeans(z[g, , drop = FALSE]))
    ax <- rowMeans(mods[, 1:4]) - rowMeans(mods[, 5:6])
    ax <- ax[match(LK[[k]]$sample, colnames(z))]
    fb <- fitb(LK[[k]], ax); bs[k] <- fb[1]; ses[k] <- fb[2]
  }
  draws[r, ] <- c(abs(bs), sapply(POOLS, function(p) abs(iv(bs[p], ses[p]))))
}
obs_lev <- c(obsb[1, K], sapply(POOLS, function(p) iv(obsb[1, p], obsb[2, p])))
nullp <- data.frame(level = colnames(draws), observed_logHR = obs_lev,
                    null_median_abs = apply(draws, 2, median), null_p95_abs = apply(draws, 2, quantile, 0.95),
                    rank = sapply(seq_along(obs_lev), function(j) 1 + sum(draws[, j] >= abs(obs_lev[j]))))
nullp$P <- nullp$rank / (NR + 1)
nullp$role <- ifelse(nullp$level %in% L3, "locked (first 3): byproduct with the 7-cohort pool; E15_null.csv is the reported value",
              ifelse(nullp$level %in% NEW, "locked (added)", ifelse(nullp$level == SENS, "sensitivity", "pooled")))
nullp$pool_genes <- length(gpool)
write.csv(nullp, file.path(O, "E20_null.csv"), row.names = FALSE)
saveRDS(draws, file.path(O, "E20_null_draws.rds"))
print(nullp[, 1:6], row.names = FALSE)

## =====================================================================================
## STEP 5. Harrell's C and paired bootstrap Delta C (added cohorts, OS)
## =====================================================================================
cat("\n== STEP 5: Harrell's C (direction fixed: higher NetImmune/ESTIMATE/CD8A = lower risk)\n")
cidx <- function(t, e, risk) concordance(Surv(t, e) ~ risk, reverse = TRUE)$concordance
CC <- list(); DC <- list()
for (k in c(NEW, SENS)) {
  d <- SO[[k]]; d <- d[complete.cases(d[, c("time", "event", "NetImmune", "ESTIMATE", "CD8A")]), ]
  cs <- sapply(c("NetImmune", "ESTIMATE", "CD8A"), function(s) cidx(d$time, d$event, -d[[s]]))
  CC[[k]] <- data.frame(cohort = k, score = names(cs), C = cs)
  B <- replicate(2000, { i <- sample(nrow(d), replace = TRUE); dd <- d[i, ]
    cn <- cidx(dd$time, dd$event, -dd$NetImmune)
    c(cn - cidx(dd$time, dd$event, -dd$ESTIMATE), cn - cidx(dd$time, dd$event, -dd$CD8A)) })
  for (j in 1:2) DC[[paste(k, j)]] <- data.frame(cohort = k, comparator = c("ESTIMATE", "CD8A")[j],
        dC = cs[1] - cs[j + 1], se = sd(B[j, ]), lo = quantile(B[j, ], 0.025), hi = quantile(B[j, ], 0.975))
}
CC <- do.call(rbind, CC); DC <- do.call(rbind, DC)
D15 <- read.csv(file.path(O, "E15_deltaC.csv"), stringsAsFactors = FALSE); D15 <- D15[D15$cohort %in% L3, ]
ivp <- function(a, lab, cmp) { w <- 1 / a$se^2; est <- sum(w * a$dC) / sum(w); se <- sqrt(1 / sum(w))
  data.frame(cohort = lab, comparator = cmp, dC = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se) }
for (cmp in c("ESTIMATE", "CD8A")) {
  a <- DC[DC$comparator == cmp & DC$cohort %in% NEW, ]
  DC <- rbind(DC, ivp(a, "locked (added 4), pooled (inverse variance)", cmp),
              ivp(rbind(a, D15[D15$comparator == cmp, names(a)]), "locked 7, pooled (inverse variance)", cmp))
}
DC$p <- 2 * pnorm(-abs(DC$dC / DC$se))
DC$role <- ifelse(DC$cohort == SENS, "sensitivity", ifelse(DC$cohort %in% NEW, "locked (added)", "pooled"))
write.csv(CC, file.path(O, "E20_cindex.csv"), row.names = FALSE)
write.csv(DC, file.path(O, "E20_deltaC.csv"), row.names = FALSE)
print(CC, row.names = FALSE)
print(transform(DC, dC = round(dC, 4), lo = round(lo, 4), hi = round(hi, 4), p = signif(p, 2)), row.names = FALSE)

## =====================================================================================
## cohort table and expectation checks
## =====================================================================================
CT <- do.call(rbind, lapply(c(NEW, SENS), function(k) {
  s <- SO[[k]]; r <- if (k %in% names(SR)) SR[[k]] else NULL
  data.frame(cohort = k, platform = PLAT[[k]], scale = XO[[k]]$scale, arrays_used = ncol(RAW[[k]]$m),
             duplicates_excluded = sum(EXCL$cohort == k), n_OS = nrow(s), deaths = sum(s$event),
             n_RFS = if (is.null(r)) NA else nrow(r), RFS_events = if (is.null(r)) NA else sum(r$event),
             n_OS_no_exclusion = nrow(SOn[[k]]), deaths_no_exclusion = sum(SOn[[k]]$event),
             stage = sum(!is.na(s$stage_ord)), age = sum(!is.na(s$age)),
             role = if (k == SENS) "sensitivity (secondary OS source)" else "locked (added)",
             outcome_source = OC[[k]]$source[1]) }))
CT <- merge(CT, do.call(rbind, SCALE)[, c("cohort", "min", "median", "max", "arrays_deposited", "genes")], by = "cohort", sort = FALSE)
write.csv(CT, file.path(O, "E20_cohort_table.csv"), row.names = FALSE)
cat("\n== cohort table\n"); print(CT[, setdiff(names(CT), "outcome_source")], row.names = FALSE)

fmt_pool <- function(r) sprintf("HR %.3f (%.3f-%.3f, %s), P = %.3g, I2 = %.0f%%", r$HR, r$lo, r$hi, r$interval, r$p, r$I2)
expect <- function(PCO, PL, tag) {
  g <- function(set, s, m = "unadjusted") PL[PL$set == set & PL$score == s & PL$model == m, ]
  u <- PCO[PCO$cohort %in% NEW & PCO$model == "unadjusted", ]
  per <- function(sc) paste(sapply(NEW, function(k) { r <- u[u$cohort == k & u$score == sc, ]
           sprintf("%s %.2f (%.2f-%.2f)", k, r$HR, r$lo, r$hi) }), collapse = "; ")
  pr <- g("locked 7 (primary)", "NetImmune")
  rbind(
    data.frame(analysis = tag, expectation = "X1 NetImmune HR < 1 in each added cohort (OS, unadjusted)", result = per("NetImmune"),
               held = sprintf("%d/4 (%d/4 with P < 0.05)", sum(u$beta[u$score == "NetImmune"] < 0), sum(u$beta[u$score == "NetImmune"] < 0 & u$p[u$score == "NetImmune"] < 0.05))),
    data.frame(analysis = tag, expectation = "X2 CAF_ECM HR > 1 in each added cohort (OS, unadjusted)", result = per("CAF_ECM"),
               held = sprintf("%d/4 (%d/4 with P < 0.05)", sum(u$beta[u$score == "CAF_ECM"] > 0), sum(u$beta[u$score == "CAF_ECM"] > 0 & u$p[u$score == "CAF_ECM"] < 0.05))),
    data.frame(analysis = tag, expectation = "X3 PRIMARY pooled NetImmune, locked 7 (unadjusted): HR < 1", result = fmt_pool(pr),
               held = if (pr$HR < 1) (if (pr$hi < 1) "yes; conservative CI excludes 1" else "direction only; conservative CI includes 1") else "no"),
    data.frame(analysis = tag, expectation = "Reported: pooled NetImmune, locked 7, adjusted", result = fmt_pool(g("locked 7 (primary)", "NetImmune", "adjusted")), held = NA),
    data.frame(analysis = tag, expectation = "Reported: pooled NetImmune, all 10 gastric", result = fmt_pool(g("all 10 gastric", "NetImmune")), held = NA),
    data.frame(analysis = tag, expectation = "Reported: pooled CAF_ECM, locked 7", result = fmt_pool(g("locked 7 (primary)", "CAF_ECM")), held = NA),
    data.frame(analysis = tag, expectation = "Reported: pooled CAF_ECM, all 10 gastric", result = fmt_pool(g("all 10 gastric", "CAF_ECM")), held = NA))
}
EXP <- rbind(expect(PCO, PL, "primary (pre-registered duplicate exclusions)"),
             data.frame(analysis = "primary (pre-registered duplicate exclusions)", expectation = "Reported: random-gene null, locked 7 pooled / added 4 pooled",
                        result = sprintf("P = %.4f / %.4f", nullp$P[nullp$level == "locked 7, pooled"], nullp$P[nullp$level == "locked (added 4), pooled"]), held = NA),
             expect(PCOn, PLn, "sensitivity (no duplicate exclusions)"))
write.csv(EXP, file.path(O, "E20_expectations.csv"), row.names = FALSE)
cat("\n== expectations\n"); print(EXP, row.names = FALSE)
saveRDS(list(OS = SO, RFS = SR), file.path(O, "E20_scores.rds"))
sink()
