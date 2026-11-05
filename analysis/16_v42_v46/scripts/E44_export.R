## E44 step 1 - export the three single-arm gastroesophageal checkpoint cohorts (z-scored, all genes) with response and
## MSI-H/EBV status. Protocol: analysis/E44_prereg.md (hash checked below). Loading and z-scoring exactly as before:
## PRJEB25780 as E12 (deposited log2, z over the cohort), PRJEB60680 as E22 primary (log2(TPM + 1) over nuclear coding,
## z within library batch, constant-within-batch genes 0), PERFECT as E16a primary (log2 vst, z over the cohort).
## Step 0: the Immune axis log-odds per SD from the exported matrices equals E12 / E16 / E22 (< 1e-6).
## Writes results/E44_Z_<cohort>.csv.gz (genes x samples) and results/E44_labels.csv.
suppressMessages({library(matrixStats); library(readxl); library(org.Hs.eg.db); library(AnnotationDbi)})
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2"; A <- file.path(G, "analysis"); O <- file.path(A, "results")
h <- strsplit(readLines(file.path(A, "E44_prereg.sha256"))[1], " +")[[1]][1]
h2 <- strsplit(system(sprintf("sha256sum %s", file.path(A, "E44_prereg.md")), intern = TRUE), " +")[[1]][1]
if (h != h2) stop("E44 protocol hash mismatch"); cat("E44 protocol sha256 verified:", substr(h, 1, 12), "\n")
sink(file.path(O, "E44_export.log"), split = TRUE)

MOD <- list(APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
            IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
            TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
            T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"))
immune <- function(Z) rowMeans(sapply(MOD, function(g) { g <- intersect(g, rownames(Z))
  if (length(g) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[g, , drop = FALSE]) }))
logit_eff <- function(y, s) {                                     # E12 / E16 / E22 logit_eff, unchanged
  s <- (s - mean(s)) / sd(s); f <- suppressWarnings(glm(y ~ s, family = binomial()))
  c(b = unname(coef(f)[2]), se = unname(sqrt(diag(vcov(f)))[2]))
}
zcohort <- function(M) { M <- M[rowSds(M) > 0, , drop = FALSE]; (M - rowMeans(M)) / rowSds(M) }
zbatch <- function(M, batch) {                                    # E22 scores_batch, unchanged
  M <- M[rowSds(M) > 0, , drop = FALSE]; Z <- M; Z[] <- NA_real_
  for (b in unique(batch)) {
    i <- which(batch == b); X <- M[, i, drop = FALSE]; sdv <- rowSds(X)
    Zb <- (X - rowMeans(X)) / sdv; Zb[sdv == 0, ] <- 0; Z[, i] <- Zb
  }
  Z
}
OUT <- list(); LAB <- list()

## ---- 1. PRJEB25780 (E12)
d <- readRDS(file.path(O, "E12_cohorts.rds"))$PRJEB25780
Z1 <- zcohort(d$M); y1 <- unname(d$y[colnames(Z1)])
k <- read.csv(file.path(G, "data_download/validation/kim2018_cps/kim2018_cps_clinical.csv"), stringsAsFactors = FALSE)
k <- k[match(colnames(Z1), k$seq_id), ]; stopifnot(!anyNA(k$seq_id), all(k$responder == y1))
me1 <- as.integer(k$msi_h %in% c("True", TRUE) | k$ebv_pos %in% c("True", TRUE))
OUT$PRJEB25780 <- Z1
LAB$PRJEB25780 <- data.frame(cohort = "PRJEB25780", sample = colnames(Z1), response = y1, msi_ebv = me1)

## ---- 2. PRJEB60680 (E22 primary)
PJ <- file.path(G, "data_download/validation/bulk_PRJEB60680_Samsung_NCT04249739"); QCD <- file.path(PJ, "quant")
qc <- read.delim(file.path(QCD, "PRJEB60680_baseline_qc.tsv"), stringsAsFactors = FALSE)
qc$QC_pass <- as.character(qc$QC_pass) %in% c("True", "TRUE")
cl <- as.data.frame(read_excel(file.path(PJ, "Wright2026_CellRepMed_TableS1_clinical_47pts.xlsx")))
names(cl) <- c("id", "scr", "sex", "age", "BOR", "fast", "pfs_ev", "os_ev", "pfs_d", "os_d", "HER2", "EBV", "PDL1",
               "cps_B", "cps_F1", "cps_F2", "MSI", "TCGA", "scRNA")
keep <- qc$sample[qc$QC_pass]; cl <- cl[match(keep, cl$id), ]; stopifnot(identical(cl$id, keep))
batch <- qc$lib_batch[match(keep, qc$sample)]
TPMNC <- as.matrix(read.delim(gzfile(file.path(QCD, "PRJEB60680_baseline_gene_TPM_nuclear_coding.tsv.gz")), row.names = 1,
                              check.names = FALSE))
Z2 <- zbatch(log2(TPMNC[, keep] + 1), batch)
y2 <- as.integer(cl$BOR %in% c("CR", "PR"))
OUT$PRJEB60680 <- Z2
LAB$PRJEB60680 <- data.frame(cohort = "PRJEB60680", sample = keep, response = y2,
                             msi_ebv = as.integer(cl$MSI == "MSI-H" | cl$EBV == "positive"))

## ---- 3. PERFECT (E16a primary)
GP <- file.path(G, "data_download/geo/GSE165252")
sm <- readLines(gzfile(file.path(GP, "GSE165252_series_matrix.txt.gz")))
field <- function(tag) gsub('"', "", strsplit(grep(paste0("^!", tag, "\t"), sm, value = TRUE)[1], "\t")[[1]][-1])
ch <- grep("^!Sample_characteristics_ch1\t", sm, value = TRUE); rl <- ch[grepl("response:", ch)]
ann <- data.frame(gsm = field("Sample_geo_accession"), title = field("Sample_title"), id = field("Sample_description"),
                  response = sub("^response: ", "", gsub('"', "", strsplit(rl, "\t")[[1]][-1])), stringsAsFactors = FALSE)
ann$patient <- sub("^Patient tumor sample_([0-9]+)_.*$", "\\1", ann$title)
ann$time <- sub("^Patient tumor sample_[0-9]+_", "", ann$title)
base <- ann[ann$time == "baseline" & ann$response %in% c("responder", "non-responder"), ]
stopifnot(!any(duplicated(base$patient))); base$y <- as.integer(base$response == "responder")
to_symbol <- function(M) {                                        # E16a to_symbol, unchanged
  ens <- sub("\\..*$", "", rownames(M)); keep <- !grepl("_PAR_Y", rownames(M))
  M <- M[keep, , drop = FALSE]; ens <- ens[keep]
  sym <- suppressMessages(mapIds(org.Hs.eg.db, keys = unique(ens), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
  s <- unname(sym[ens]); ok <- !is.na(s); M <- M[ok, , drop = FALSE]; s <- s[ok]
  o <- order(-rowMeans(M)); M <- M[o, , drop = FALSE]; s <- s[o]
  M <- M[!duplicated(s), , drop = FALSE]; rownames(M) <- s[!duplicated(s)]; M
}
vst <- as.matrix(read.delim(gzfile(file.path(GP, "GSE165252_vst_PERFECT.txt.gz")), row.names = 1, check.names = FALSE))
Z3 <- zcohort(log2(to_symbol(vst[, base$id])))
OUT$PERFECT <- Z3
LAB$PERFECT <- data.frame(cohort = "PERFECT", sample = base$id, response = base$y, msi_ebv = NA_integer_)

## ---- step 0
ref <- c(PRJEB25780 = 1.320745, PRJEB60680 = 0.397973, PERFECT = 0.400282)
r12 <- read.csv(file.path(O, "E12_icb_per_cohort.csv")); r12 <- r12[r12$cohort == "PRJEB25780" & r12$score == "Immune", ]
r16 <- read.csv(file.path(O, "E16_perfect_per_score.csv")); r16 <- r16[r16$scale == "primary_log2_vst" & r16$score == "Immune", ]
r22 <- read.csv(file.path(O, "E22_PRJEB60680_response.csv")); r22 <- r22[r22$score == "Immune" & grepl("^primary", r22$analysis), ]
ref[] <- c(r12$logOR, r22$logOR, r16$logOR)
s0 <- do.call(rbind, lapply(names(OUT), function(nm) {
  e <- logit_eff(LAB[[nm]]$response, immune(OUT[[nm]]))
  data.frame(cohort = nm, n = ncol(OUT[[nm]]), responders = sum(LAB[[nm]]$response), genes = nrow(OUT[[nm]]),
             immune_logOR = e[["b"]], reference = ref[[nm]], diff = e[["b"]] - ref[[nm]]) }))
print(s0, digits = 7); stopifnot(all(abs(s0$diff) < 1e-6)); cat("step 0 passed (max |diff|", signif(max(abs(s0$diff)), 3), ")\n")
for (nm in names(OUT)) {
  stopifnot(identical(colnames(OUT[[nm]]), LAB[[nm]]$sample))
  write.csv(data.frame(gene = rownames(OUT[[nm]]), OUT[[nm]], check.names = FALSE),
            gzfile(file.path(O, sprintf("E44_Z_%s.csv.gz", nm))), row.names = FALSE)
}
L <- do.call(rbind, LAB); write.csv(L, file.path(O, "E44_labels.csv"), row.names = FALSE)
print(aggregate(cbind(n = 1, responders = response, msi_ebv = msi_ebv) ~ cohort, data = transform(L, msi_ebv = ifelse(is.na(msi_ebv), 0, msi_ebv)), FUN = sum))
sink()
