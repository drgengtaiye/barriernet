## E32 sensitivity: the 52 companion signatures scored with their published engines, following the companion's
## scripts/score_matrices.R (load_sig, engine(), gsvaParam / ssgseaParam with minSize = 3, weighted mean of row-z with
## published weights and missing weights 0). Usage: Rscript E32_native_scores.R <dir> <comparison keys...>
## Reads <dir>/<key>.tsv.gz (genes x samples, log scale), writes <dir>/<key>.scores.csv (samples x signatures).
suppressMessages({library(GSVA); library(matrixStats)})
args <- commandArgs(TRUE); T <- args[1]; KEYS <- args[-1]
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
D <- file.path(G, "icb")
meta <- read.csv(file.path(D, "sig_methods_clean.csv"), stringsAsFactors = FALSE)
lib52 <- read.csv(file.path(D, "v5_signature_clusters.csv"), stringsAsFactors = FALSE)$signature
ALIAS <- c(WARS1 = "WARS", HLA_DPB1 = "HLA-DPB1")
load_sig <- function(nm) {
  d <- read.csv(file.path(D, "SignatureSets-main/data-raw/IO", paste0(nm, ".csv")), stringsAsFactors = FALSE, fileEncoding = "latin1")
  g <- d$gene_name; g <- ifelse(g %in% names(ALIAS), ALIAS[g], g)
  w <- if ("weight" %in% names(d)) suppressWarnings(as.numeric(d$weight)) else rep(NA_real_, nrow(d))
  ok <- !is.na(g) & g != ""; list(genes = g[ok], weights = w[ok])
}
engine <- function(m) if (grepl("Weighted", m)) "wmean" else if (grepl("ssGSEA", m)) "ssgsea" else "gsva"
for (k in KEYS) {
  M <- as.matrix(read.delim(gzfile(file.path(T, paste0(k, ".tsv.gz"))), row.names = 1, check.names = FALSE))
  M <- M[rowSds(M) > 0, ]
  Z <- (M - rowMeans(M)) / rowSds(M)
  obs <- list(); gs <- list(); ss <- list()
  for (s_ in lib52) {
    S <- load_sig(s_); g <- intersect(unique(S$genes), rownames(M)); if (length(g) < 3) next
    e <- engine(meta$method[meta$signature == s_])
    if (e == "wmean") {
      w <- S$weights[match(g, S$genes)]; if (all(is.na(w))) w <- rep(1, length(g)); w[is.na(w)] <- 0
      if (sum(abs(w)) == 0) next
      obs[[s_]] <- as.numeric(crossprod(Z[g, , drop = FALSE], w)) / sum(abs(w))
    } else if (e == "ssgsea") ss[[s_]] <- g else gs[[s_]] <- g
  }
  if (length(gs)) { rr <- gsva(gsvaParam(M, gs, minSize = 3), verbose = FALSE); for (n in rownames(rr)) obs[[n]] <- as.numeric(rr[n, ]) }
  if (length(ss)) { rr <- gsva(ssgseaParam(M, ss, minSize = 3), verbose = FALSE); for (n in rownames(rr)) obs[[n]] <- as.numeric(rr[n, ]) }
  out <- as.data.frame(obs, check.names = FALSE); rownames(out) <- colnames(M)
  write.csv(out, file.path(T, paste0(k, ".scores.csv")))
  cat(k, ":", ncol(out), "signatures scored on", nrow(M), "genes\n")
}
