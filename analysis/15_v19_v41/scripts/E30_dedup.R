## E30 Q1 - duplicate rule between GSE26253 (SMC, Illumina WG-DASL, FFPE) and GSE62254 (ACRG, Affymetrix U133 Plus 2).
## Protocol: analysis/E30_prereg.md (sha256 in E30_prereg.sha256), section Q1. Expression and clinical data only;
## no outcome model is fitted here.
##
## Rule (verbatim from the protocol): an array of GSE26253 is a duplicate of an array of GSE62254 if
##   (i)   the two are each other's best cross-cohort expression match,
##   (ii)  sex agrees and age differs by at most one year, and
##   (iii) their times do not contradict each other: where GSE26253 records recurrence, GSE62254 disease-free time
##         and overall survival are each within 0.5 months of GSE26253 recurrence-free time and overall survival;
##         where GSE26253 is censored, GSE62254 disease-free time and overall survival are each at least the
##         corresponding GSE26253 time minus 0.5 months.
## Duplicates are removed from GSE26253. Sensitivity: only the arrays with matching event times (rule met with
## GSE26253 recurrence recorded).
##
## Implementation:
##   - GSE26253 loaded exactly as in E15_gastric_validation.R (series matrix, GPL8432 Symbol, max-mean collapse,
##     samples with RFS time > 0 and status recorded: 432 arrays);
##   - GSE62254 = the benchmark analysis matrix (revision/results/B1_GSE62254.rds, 283 arrays), as used by E15/E20;
##   - (i) dup_pair() copied verbatim from E20_extended_replication.R: all shared genes with no NA and SD > 0 in both
##     cohorts, genes z-scored within each cohort, Pearson; 'reciprocal' = each other's best match;
##   - (ii) GSE26253 age/sex (and OS) are not in GEO: taken from Oh et al. 2018 Nat Commun Suppl. Data 2, sheet SMC
##     (the same 432 patients), linked by the unique key (RFS months to 0.1, recurrence status, AJCC stage group) as in
##     the read-only QC check (rev/agentB/dup_screen_all10.R). GSE62254 age/sex/OS/DFS: same file, sheet ACRG (GEO_ID).
##   - (iii) as written; unrounded times.
suppressMessages({library(matrixStats); library(readxl)})
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
D <- file.path(G, "2/data_download/geo")
V <- file.path(G, "2/data_download/validation/e20_outcomes")
R <- file.path(G, "revision/results")
O <- file.path(G, "2/analysis/results")
sink(file.path(O, "E30_dedup.log"), split = TRUE)
cat("E30 Q1 duplicate rule run", format(Sys.time()), "\nprotocol sha256:", readLines(file.path(G, "2/analysis/E30_prereg.sha256"))[1], "\n")
h <- system2("sha256sum", file.path(G, "2/analysis/E30_prereg.md"), stdout = TRUE)
stopifnot(sub(" .*", "", h) == sub(" .*", "", readLines(file.path(G, "2/analysis/E30_prereg.sha256"))[1]))
cat("protocol hash verified\n")

## ---- E15 helpers (verbatim) ----
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
a8432 <- read.delim(file.path(D, "platforms/GPL8432_id_symbol.tsv"), stringsAsFactors = FALSE)
sym8432 <- setNames(a8432$Symbol, a8432$ID)
prep <- function(g, symmap) {
  x <- read_matrix(g); m <- x$m
  rng <- range(m, na.rm = TRUE)
  if (rng[2] > 50) m <- log2(pmax(m, 0) + 1)
  m <- m[rowSums(is.na(m)) == 0, , drop = FALSE]
  m <- collapse_max(m, unname(symmap[rownames(m)]))
  list(m = m, meta = x$meta)
}
p <- suppressWarnings(prep("GSE26253", sym8432))
st <- meta_field(p$meta, "Sample_characteristics_ch1")
key <- function(k) { i <- which(sapply(st, function(v) any(grepl(k, v, fixed = TRUE)))); sub(".*: ", "", st[[i[1]]]) }
gsm <- meta_field(p$meta, "Sample_geo_accession")[[1]]
stg <- key("pathological stage")
cl <- data.frame(sample = gsm, time = as.numeric(key("recurrence free survival time")),
                 event = as.integer(key("status (0=non-recurrence")), stage = stg, stringsAsFactors = FALSE)
cl <- cl[!is.na(cl$time) & cl$time > 0 & !is.na(cl$event), ]
A <- p$m[, cl$sample]
B62 <- readRDS(file.path(R, "B1_GSE62254.rds")); B <- B62$expr
cat(sprintf("GSE26253: %d genes x %d arrays; GSE62254 (B1): %d genes x %d arrays\n", nrow(A), ncol(A), nrow(B), ncol(B)))

## ---- (i) dup_pair() verbatim from E20_extended_replication.R ----
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
         max_gap = max(r1 - r2), stringsAsFactors = FALSE), C = C)
}
x <- dup_pair(A, B, "GSE26253", "GSE62254")
cat(sprintf("shared genes %d; within-GSE26253 99.9th pct %.3f; cross-cohort max r %.3f; flagged by the E20 rule: %d\n",
            x$pair$genes, x$pair$within_p999, x$pair$cross_max, x$pair$n_flagged))
S <- x$smp
## z of the best match within its row (descriptive)
S$row_z <- (S$max_r - rowMeans(x$C)) / rowSds(x$C)
M <- S[S$reciprocal, ]
cat("mutual best matches (rule i):", nrow(M), "\n")

## ---- (ii) clinical data ----
f <- file.path(V, "oh2018_natcommun/41467_2018_4179_MOESM5_ESM.xlsx")
smc <- as.data.frame(read_excel(f, sheet = "SMC"))
acrg <- as.data.frame(read_excel(f, sheet = "ACRG"))
rom <- c("1"="I","1A"="I","1B"="I","2"="II","3A"="III","3B"="III","4"="IV")
k26 <- paste(sprintf("%.1f", cl$time), cl$event, rom[cl$stage])
ksm <- paste(sprintf("%.1f", as.numeric(smc$RFS.m)), smc[["Recurrence (1=yes, 0=no)"]], smc$AJCC.stage)
uniq <- k26 %in% names(which(table(k26) == 1)) & k26 %in% names(which(table(ksm) == 1))
j <- match(k26, ksm)
cl$smc_id <- NA; cl$smc_id[uniq] <- smc$Patients_ID[j[uniq]]
cl$age <- NA_real_; cl$age[uniq] <- as.numeric(smc$Age[j[uniq]])
cl$sex <- NA_character_; cl$sex[uniq] <- toupper(substr(smc$Sex[j[uniq]], 1, 1))
cl$os <- NA_real_; cl$os[uniq] <- as.numeric(smc$OS.m[j[uniq]])
cat(sprintf("GSE26253 arrays linked to the SMC sheet by a unique key: %d of %d\n", sum(uniq), nrow(cl)))
M$age26 <- cl$age[match(M$sample, cl$sample)]; M$sex26 <- cl$sex[match(M$sample, cl$sample)]
M$smc_id <- cl$smc_id[match(M$sample, cl$sample)]
M$rec26 <- cl$event[match(M$sample, cl$sample)]; M$rfs26 <- cl$time[match(M$sample, cl$sample)]
M$os26 <- cl$os[match(M$sample, cl$sample)]; M$stage26 <- cl$stage[match(M$sample, cl$sample)]
ia <- match(M$best_match, acrg$GEO_ID)
M$age62 <- as.numeric(acrg$age[ia]); M$sex62 <- toupper(substr(acrg$sex[ia], 1, 1))
M$dfs62 <- as.numeric(acrg$DFS.m[ia]); M$os62 <- as.numeric(acrg$OS.m[ia])
M$recur62 <- acrg$Recur[ia]; M$death62 <- acrg$Death[ia]; M$stage62 <- acrg$pStage[ia]
M$rule_i <- TRUE
M$rule_ii <- !is.na(M$age26) & !is.na(M$age62) & !is.na(M$sex26) & !is.na(M$sex62) &
             M$sex26 == M$sex62 & abs(M$age26 - M$age62) <= 1
TOL <- 0.5
M$rule_iii <- with(M, ifelse(is.na(rec26) | is.na(rfs26) | is.na(os26) | is.na(dfs62) | is.na(os62), FALSE,
               ifelse(rec26 == 1, abs(dfs62 - rfs26) <= TOL & abs(os62 - os26) <= TOL,
                                  dfs62 >= rfs26 - TOL & os62 >= os26 - TOL)))
M$duplicate <- M$rule_i & M$rule_ii & M$rule_iii
M$matching_event_times <- M$duplicate & M$rec26 == 1
M <- M[order(-M$duplicate, -M$rule_ii, -M$row_z), ]

## unresolved clinical linkage among mutual best matches: list every SMC row sharing the array's key and test
## whether ANY of them would satisfy (ii) and (iii) with the GSE62254 partner (diagnostic only; rule needs a link)
unl <- M[is.na(M$age26), ]
if (nrow(unl)) {
  cat("\nmutual best matches without a unique SMC link (rule ii cannot be evaluated):", nrow(unl), "\n")
  for (i in seq_len(nrow(unl))) {
    kk <- k26[match(unl$sample[i], cl$sample)]; cand <- smc[ksm == kk, ]
    ok <- if (nrow(cand)) with(cand, toupper(substr(Sex, 1, 1)) == unl$sex62[i] & abs(as.numeric(Age) - unl$age62[i]) <= 1) else logical()
    cat(sprintf("  %s (key %s; %d SMC rows with this key, of which %d agree with %s on sex and age +/-1: %s)\n",
                unl$sample[i], kk, nrow(cand), sum(ok), unl$best_match[i],
                if (any(ok)) paste(cand$Patients_ID[ok], collapse = " ") else "-"))
  }
}

dups <- M$sample[M$duplicate]; sens3 <- M$sample[M$matching_event_times]
EXPECT <- c("GSM644726", "GSM644801", "GSM644802", "GSM644593", "GSM644733", "GSM644810", "GSM644871", "GSM644940")
EXPECT3 <- c("GSM644726", "GSM644801", "GSM644802")
cat("\nmutual best matches agreeing on sex and age +/- 1 year (i + ii):", sum(M$rule_ii), "\n")
cat("duplicates (i + ii + iii):", length(dups), ":", paste(sort(dups), collapse = " "), "\n")
cat("matching event times (sensitivity):", length(sens3), ":", paste(sort(sens3), collapse = " "), "\n")
cat("rule selects the expected 8:", setequal(dups, EXPECT), "; expected 3:", setequal(sens3, EXPECT3), "\n")
if (!setequal(dups, EXPECT)) cat("  DIFFERENCE: selected not expected:", setdiff(dups, EXPECT), "; expected not selected:", setdiff(EXPECT, dups), "\n")
options(width = 250)
cols <- c("sample", "best_match", "max_r", "gap", "row_z", "smc_id", "sex26", "sex62", "age26", "age62", "stage26", "stage62",
          "rec26", "rfs26", "dfs62", "os26", "os62", "death62", "rule_ii", "rule_iii", "duplicate", "matching_event_times")
cat("\nmutual best matches with sex and age agreement:\n")
print(transform(M[M$rule_ii, cols], max_r = round(max_r, 3), gap = round(gap, 3), row_z = round(row_z, 2),
                rfs26 = round(rfs26, 2), dfs62 = round(dfs62, 2), os26 = round(os26, 2), os62 = round(os62, 2)), row.names = FALSE)
out <- M[, c("cohort", "comparator", "genes", cols)]
out$set_primary <- out$duplicate; out$set_sens3 <- out$matching_event_times
write.csv(out, file.path(O, "E30_dedup_duplicates.csv"), row.names = FALSE)
cat("\nwritten: results/E30_dedup_duplicates.csv (all", nrow(out), "mutual best matches with rule columns)\n")
sink()
