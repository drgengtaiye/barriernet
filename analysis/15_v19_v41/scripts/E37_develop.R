## E37 step 3 - develop barriernet 0.3 on the 20 single-arm checkpoint cohorts (no randomized data), select, freeze.
## Protocol: analysis/E37_prereg.md (hash checked). Input: results/E37_dev_genes.csv (E37_export.R) minus the samples
## flagged by E37_twins.py. Candidates C0-C4, leave-one-cohort-out (LOCO) selection, refit on all 20, frozen weights in
## results/E37_frozen_weights.tsv (hash logged by the caller in E37_frozen.sha256 before E37_test.py runs).
suppressMessages({library(glmnet); library(metafor)})
A <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"; O <- file.path(A, "results")
h <- strsplit(readLines(file.path(A, "E37_prereg.sha256"))[1], " +")[[1]][1]
h2 <- strsplit(system(sprintf("sha256sum %s", file.path(A, "E37_prereg.md")), intern = TRUE), " +")[[1]][1]
if (h != h2) stop("E37 protocol hash mismatch"); cat("E37 protocol sha256 verified:", substr(h, 1, 12), "\n")

MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"))
G48 <- unlist(MOD, use.names = FALSE)

d <- read.csv(file.path(O, "E37_dev_genes.csv"), check.names = FALSE, stringsAsFactors = FALSE)
fl <- read.csv(file.path(O, "E37_twin_flags.csv"), stringsAsFactors = FALSE)
drop <- paste(d$cohort, d$sample) %in% paste(fl$development, fl$sample)
cat("removed as possible twins:", sum(drop), "\n"); d <- d[!drop, ]
COH <- unique(d$cohort); stopifnot(length(COH) == 20)
Graw <- as.matrix(d[, G48]); G <- Graw; G[!is.finite(G)] <- 0            # gene not measured -> 0
y <- d$response; coh <- factor(d$cohort, levels = COH)

## module scores: 0.2.1 rule (NA below three genes) for C0/C1; 0 below three genes for C2
modsc <- function(Gr, zero) {
  S <- sapply(names(MOD), function(k) {
    out <- rep(NA_real_, nrow(Gr))
    for (cc in COH) { i <- which(d$cohort == cc); g <- MOD[[k]][colSums(is.finite(Gr[i, MOD[[k]], drop = FALSE])) == length(i)]
      out[i] <- if (length(g) >= 3) rowMeans(Gr[i, g, drop = FALSE]) else if (zero) 0 else NA_real_ }
    out })
  as.data.frame(S)
}
S0 <- modsc(Graw, FALSE); S0$Immune <- rowMeans(S0[, c("APM", "IFN", "TLS_B", "T_CELL")])
S0$NetImmune <- S0$Immune - rowMeans(S0[, c("CAF_ECM", "TAM_SUP")])
M6 <- as.matrix(modsc(Graw, TRUE))

logit_eff <- function(y, s) { if (sd(s) == 0) return(c(b = NA, se = NA)); s <- (s - mean(s)) / sd(s)
  f <- suppressWarnings(glm(y ~ s, family = binomial())); c(b = unname(coef(f)[2]), se = unname(sqrt(diag(vcov(f)))[2])) }
dl_pool <- function(b, se) { ok <- is.finite(b) & is.finite(se); r <- rma(yi = b[ok], sei = se[ok], method = "DL")
  c(est = as.numeric(r$b), lo = r$ci.lb, hi = r$ci.ub, p = r$pval, k = sum(ok)) }

## design: features + unpenalized cohort indicators (first cohort of the training set is the reference)
fitpen <- function(X, yy, cc, alpha, foldid = NULL) {
  cc <- droplevels(cc); D <- model.matrix(~ cc)[, -1, drop = FALSE]
  XX <- cbind(X, D); pf <- c(rep(1, ncol(X)), rep(0, ncol(D)))
  fid <- as.integer(cc)
  cv <- cv.glmnet(XX, yy, family = "binomial", alpha = alpha, penalty.factor = pf, standardize = FALSE,
                  foldid = fid, type.measure = "deviance", nfolds = max(fid))
  b <- as.numeric(coef(cv, s = "lambda.min"))[-1][seq_len(ncol(X))]; names(b) <- colnames(X)
  list(b = b, lambda = cv$lambda.min)
}
set.seed(2037)
CAND <- list(C2 = list(X = M6, alpha = 0), C3 = list(X = G, alpha = 0.5), C4 = list(X = G, alpha = 0))

## ---- LOCO
loco <- list()
for (cc in COH) {
  te <- d$cohort == cc; tr <- !te
  add <- function(name, s) { e <- logit_eff(y[te], s)
    loco[[length(loco) + 1]] <<- data.frame(candidate = name, cohort = cc, n = sum(te), b = e[["b"]], se = e[["se"]]) }
  if (all(is.finite(S0$NetImmune[te]))) add("C0", S0$NetImmune[te]) else add("C0", rep(0, sum(te)))
  if (all(is.finite(S0$Immune[te]))) add("C1", S0$Immune[te]) else add("C1", rep(0, sum(te)))
  for (nm in names(CAND)) {
    X <- CAND[[nm]]$X; f <- fitpen(X[tr, ], y[tr], coh[tr], CAND[[nm]]$alpha)
    add(nm, as.numeric(X[te, , drop = FALSE] %*% f$b))
  }
  cat(".")
}
cat("\n")
LO <- do.call(rbind, loco); write.csv(LO, file.path(O, "E37_loco_per_cohort.csv"), row.names = FALSE)
PO <- do.call(rbind, lapply(c("C0", "C1", "C2", "C4", "C3"), function(nm) {
  x <- LO[LO$candidate == nm, ]; data.frame(candidate = nm, t(dl_pool(x$b, x$se))) }))
PO$label <- c(C0 = "NetImmune 0.2.1", C1 = "Immune axis 0.2.1", C2 = "module weights (ridge)",
              C4 = "gene weights (ridge)", C3 = "gene weights (elastic net)")[PO$candidate]
print(PO, digits = 3)

## ---- selection: highest pooled held-out log-odds; within 0.02 of a simpler candidate -> the simpler
ord <- c("C0", "C1", "C2", "C4", "C3"); est <- setNames(PO$est, PO$candidate)
best <- names(which.max(est[ord]))
for (nm in ord) { if (nm == best) break; if (est[[nm]] >= est[[best]] - 0.02) { best <- nm; break } }
cat("selected:", best, "\n")

## ---- refit on all 20 and freeze as gene weights
if (best == "C0") {
  w <- c(setNames(rep(0.25 / sapply(MOD[3:6], length), sapply(MOD[3:6], length)), NULL)); stop("C0 selected: freeze NetImmune 0.2.1 weights by hand")
}
if (best == "C1") {
  W <- unlist(lapply(MOD[3:6], function(g) setNames(rep(0.25 / length(g), length(g)), g)))
  lam <- NA; bm <- NULL
} else {
  f <- fitpen(CAND[[best]]$X, y, coh, CAND[[best]]$alpha); lam <- f$lambda
  if (best == "C2") { bm <- f$b; W <- unlist(lapply(names(MOD), function(k) setNames(rep(bm[[k]] / length(MOD[[k]]), length(MOD[[k]])), MOD[[k]]))) }
  else { bm <- NULL; W <- f$b }
}
names(W) <- sub("^[A-Z_]+\\.", "", names(W)); W <- W[W != 0]; W <- W / sum(abs(W))
FW <- data.frame(gene = names(W), weight = unname(W), module = sapply(names(W), function(g) names(MOD)[sapply(MOD, function(v) g %in% v)][1]))
FW <- FW[order(-abs(FW$weight)), ]
write.table(FW, file.path(O, "E37_frozen_weights.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (!is.null(bm)) { cat("module coefficients (all 20 cohorts):\n"); print(round(bm, 4)) }
cat("lambda:", lam, "; genes with non-zero weight:", nrow(FW), "\n"); print(FW, digits = 3, row.names = FALSE)

## ---- development expectations
c2 <- fitpen(M6, y, coh, 0)$b
D1 <- max(est[c("C2", "C3", "C4")]) - est[["C1"]]
EXP <- data.frame(q = c("D1", "D2"),
  expectation = c("no learned candidate exceeds the Immune axis by more than 0.10 (pooled held-out log-odds)",
                  "fitted CAF_ECM weight (C2, all 20 cohorts) at most zero"),
  result = c(sprintf("largest learned minus Immune axis = %+.3f", D1), sprintf("CAF_ECM coefficient %+.4f", c2[["CAF_ECM"]])),
  held = c(D1 <= 0.10, c2[["CAF_ECM"]] <= 0))
write.csv(PO, file.path(O, "E37_loco_pooled.csv"), row.names = FALSE)
write.csv(EXP, file.path(O, "E37_dev_expectations.csv"), row.names = FALSE)
write.csv(data.frame(module = names(c2), C2_coefficient = unname(c2)), file.path(O, "E37_C2_module_coefficients.csv"), row.names = FALSE)
print(EXP)
