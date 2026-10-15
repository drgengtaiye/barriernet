## E16a - PERFECT oesophageal adenocarcinoma checkpoint cohort (GSE165252; van den Ende 2021):
## neoadjuvant atezolizumab + chemoradiotherapy. The compendium copy (ICB_VanDenEnde) kept only
## 3 of the 48 module genes (E12), so the deposited GEO matrices are scored here.
##
## Expected directions, written before any PERFECT outcome was computed (from E5/E12):
##   IFN, Immune, APM, TLS_B, T_CELL, NetImmune  > 0 (higher in responders)
##   CAF_ECM                                     ~ 0 (no consistent association with response)
##   TAM_SUP, C1Q_arm, SPP1_arm, CXCL9_SPP1      leaning > 0
## PERFECT gives chemoradiotherapy with the checkpoint antibody, so under the E12 entry rules it
## belongs with the chemo-ICB sensitivity cohorts (GSE207422, GSE235919), not the primary sets.
##
## Samples: pre-treatment ("baseline") endoscopic biopsies, one per patient, labelled
## responder / non-responder in the series matrix (the compendium's R/NR labels on the same GSMs
## are compared as a check). Expression: GSE165252_vst_PERFECT (values on a linear, count-like
## scale; log2 taken, primary) and GSE165252_norm.cnt_PERFECT (log2(x + 1), sensitivity).
## Ensembl IDs (version stripped) -> HGNC symbol with org.Hs.eg.db; duplicate symbols collapsed
## to the row with the highest mean. Scoring and effect exactly as E12 (row-z within cohort,
## module = mean z of >= 3 genes, logistic log-odds of response per SD).
## Pooling: the E12 per-cohort table (results/E12_icb_per_cohort.csv, unchanged) plus PERFECT;
## DerSimonian-Laird and Hartung-Knapp, reporting the wider interval and the larger P.
suppressMessages({library(matrixStats); library(metafor); library(org.Hs.eg.db); library(MultiAssayExperiment)})
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/data_download/geo/GSE165252"

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

## ---- sample annotation from the series matrix
sm <- readLines(gzfile(file.path(G, "GSE165252_series_matrix.txt.gz")))
field <- function(tag, which = 1) {
  l <- grep(paste0("^!", tag, "\t"), sm, value = TRUE)[which]
  gsub('"', "", strsplit(l, "\t")[[1]][-1])
}
ch <- grep("^!Sample_characteristics_ch1\t", sm, value = TRUE)
resp_line <- ch[grepl("response:", ch)]
ann <- data.frame(gsm = field("Sample_geo_accession"), title = field("Sample_title"),
                  id = field("Sample_description"),
                  response = sub("^response: ", "", gsub('"', "", strsplit(resp_line, "\t")[[1]][-1])),
                  stringsAsFactors = FALSE)
ann$patient <- sub("^Patient tumor sample_([0-9]+)_.*$", "\\1", ann$title)
ann$time <- sub("^Patient tumor sample_[0-9]+_", "", ann$title)
write.csv(ann, file.path(O, "E16_perfect_samples_all.csv"), row.names = FALSE)
cat("GEO samples:", nrow(ann), "\n"); print(table(ann$time, ann$response, useNA = "ifany"))

base <- ann[ann$time == "baseline" & ann$response %in% c("responder", "non-responder"), ]
stopifnot(!any(duplicated(base$patient)))
base$y <- as.integer(base$response == "responder")

## check against the compendium labels on the same GSMs
cmp <- as.data.frame(colData(readRDS("/groups/adv2105_gp/yichen/Yi/com/gastric/icb/ICB_VanDenEnde.rds")))
cmp$gsm <- rownames(cmp)
chk <- merge(base[, c("gsm", "response")], cmp[, c("gsm", "response")], by = "gsm", all = TRUE, suffixes = c("_geo", "_compendium"))
cat("GEO baseline labelled:", nrow(base), " compendium labelled:", sum(!is.na(cmp$response)),
    " overlap:", sum(!is.na(chk$response_geo) & !is.na(chk$response_compendium)),
    " agree:", sum((chk$response_geo == "responder") == (chk$response_compendium == "R"), na.rm = TRUE), "\n")
write.csv(chk, file.path(O, "E16_perfect_label_check.csv"), row.names = FALSE)

## ---- expression
read_mat <- function(f) {
  x <- read.delim(gzfile(file.path(G, f)), row.names = 1, check.names = FALSE)
  as.matrix(x)
}
to_symbol <- function(M) {
  ens <- sub("\\..*$", "", rownames(M))
  keep <- !grepl("_PAR_Y", rownames(M))
  M <- M[keep, , drop = FALSE]; ens <- ens[keep]
  sym <- suppressMessages(mapIds(org.Hs.eg.db, keys = unique(ens), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
  s <- unname(sym[ens]); ok <- !is.na(s)
  M <- M[ok, , drop = FALSE]; s <- s[ok]
  o <- order(-rowMeans(M)); M <- M[o, , drop = FALSE]; s <- s[o]
  M <- M[!duplicated(s), , drop = FALSE]; rownames(M) <- s[!duplicated(s)]
  M
}
vst <- read_mat("GSE165252_vst_PERFECT.txt.gz")
cnt <- read_mat("GSE165252_norm.cnt_PERFECT.txt.gz")
cat("vst range:", signif(range(vst), 4), " median:", signif(median(vst), 4),
    "| norm.cnt range:", signif(range(cnt), 4), " median:", signif(median(cnt), 4), "\n")
cat("column sums (first 3) vst:", signif(colSums(vst)[1:3], 4), " norm.cnt:", signif(colSums(cnt)[1:3], 4), "\n")
stopifnot(all(base$id %in% colnames(vst)), all(base$id %in% colnames(cnt)))
SCALES <- list(primary_log2_vst = log2(to_symbol(vst[, base$id])),
               sens_log2_normcnt = log2(to_symbol(cnt[, base$id]) + 1))

scores_sym <- function(M) {
  M <- M[rowSds(M) > 0, , drop = FALSE]; Z <- (M - rowMeans(M)) / rowSds(M)
  used <- sapply(MOD, function(g) sum(g %in% rownames(Z)))
  S <- as.data.frame(sapply(names(MOD), function(k) {
    g <- intersect(MOD[[k]], rownames(Z)); if (length(g) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[g, , drop = FALSE]) }))
  S$Barrier <- rowMeans(S[, c("CAF_ECM", "TAM_SUP")])
  S$Immune <- rowMeans(S[, c("APM", "IFN", "TLS_B", "T_CELL")])
  S$NetImmune <- S$Immune - S$Barrier
  S$CXCL9_SPP1 <- if (all(c("CXCL9", "SPP1") %in% rownames(Z))) Z["CXCL9", ] - Z["SPP1", ] else NA_real_
  list(S = S, used = used)
}
logit_eff <- function(y, s) {
  s <- (s - mean(s)) / sd(s)
  f <- suppressWarnings(glm(y ~ s, family = binomial()))
  c(b = unname(coef(f)[2]), se = unname(sqrt(diag(vcov(f)))[2]))
}
auc <- function(y, v) { r <- rank(v); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

rows <- list()
for (sc_name in names(SCALES)) {
  sc <- scores_sym(SCALES[[sc_name]])
  if (sc_name == "primary_log2_vst") { cat("module genes present:\n"); print(sc$used) }
  for (k in SC) {
    v <- sc$S[[k]]; if (all(is.na(v))) next
    e <- logit_eff(base$y, v)
    rows[[length(rows) + 1]] <- data.frame(scale = sc_name, cohort = "PERFECT_GSE165252", tumour = "oesophageal adenocarcinoma",
      n = nrow(base), n_resp = sum(base$y), score = k, genes = if (k %in% names(sc$used)) sc$used[[k]] else NA,
      logOR = e[["b"]], se = e[["se"]], lo = e[["b"]] - 1.96 * e[["se"]], hi = e[["b"]] + 1.96 * e[["se"]],
      p = 2 * pnorm(-abs(e[["b"]] / e[["se"]])), AUC = auc(base$y, v), tech = "WT", set = "chemoICB")
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(O, "E16_perfect_per_score.csv"), row.names = FALSE)
cat("\n== PERFECT, log-odds of response per SD\n")
print(res[, c("scale", "score", "genes", "logOR", "lo", "hi", "p", "AUC")], digits = 3, row.names = FALSE)

## ---- pooling with the E12 per-cohort table
per <- read.csv(file.path(O, "E12_icb_per_cohort.csv"), stringsAsFactors = FALSE)
per <- per[per$cohort != "PRJEB25780 (E7a double log)", ]
pf <- res[res$scale == "primary_log2_vst", names(per)]
all_per <- rbind(per, pf)
E5C <- unique(per$cohort[per$set == "E5"]); NEW <- unique(per$cohort[per$set == "new"])
CHEMO <- unique(per$cohort[per$set == "chemoICB"])
pool_one <- function(p, set) {
  if (nrow(p) < 2) return(NULL)
  m <- rma(yi = p$logOR, sei = p$se, method = "DL")
  h <- rma(yi = p$logOR, sei = p$se, method = "DL", test = "knha")
  hk_wider <- (h$ci.ub - h$ci.lb) >= (m$ci.ub - m$ci.lb)
  data.frame(set = set, score = p$score[1], k = nrow(p), n = sum(p$n), n_resp = sum(p$n_resp),
             logOR = m$b[1], lo_DL = m$ci.lb, hi_DL = m$ci.ub, p_DL = m$pval,
             lo_HK = h$ci.lb, hi_HK = h$ci.ub, p_HK = h$pval,
             lo = if (hk_wider) h$ci.lb else m$ci.lb, hi = if (hk_wider) h$ci.ub else m$ci.ub, p = max(m$pval, h$pval),
             I2 = m$I2, n_pos = sum(p$logOR > 0))
}
SETS <- list(
  "all + chemo-ICB (E12)" = c(E5C, NEW, CHEMO),
  "all + chemo-ICB + PERFECT" = c(E5C, NEW, CHEMO, "PERFECT_GSE165252"),
  "chemo-ICB only (GSE207422, GSE235919, PERFECT)" = c(CHEMO, "PERFECT_GSE165252"),
  "gastro-oesophageal (PRJEB25780 + PERFECT)" = c("PRJEB25780", "PERFECT_GSE165252"))
pooled <- do.call(rbind, lapply(names(SETS), function(s) do.call(rbind, lapply(SC, function(k)
  pool_one(all_per[all_per$score == k & all_per$cohort %in% SETS[[s]], ], s)))))
write.csv(pooled, file.path(O, "E16_perfect_pooled.csv"), row.names = FALSE)
f <- function(x) sprintf("%+.2f", x)
for (s in names(SETS)) {
  cat("\n==", s, "\n"); q <- pooled[pooled$set == s, ]
  for (i in seq_len(nrow(q))) with(q[i, ], cat(sprintf("  %-11s k=%2d  %s (%s to %s; P=%.3g)  I2=%2.0f%%  pos %d/%d\n",
    score, k, f(logOR), f(lo), f(hi), p, I2, n_pos, k)))
}
