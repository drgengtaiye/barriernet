## E15 - LOCKED external validation of the barriernet axes in gastric cohorts downloaded after v16
## (data_download/geo/, series matrices + supplementary outcome files; get_gastric_geo.sh).
##
## PRE-SPECIFIED HYPOTHESES (written before any outcome was examined; nothing below may be tuned):
##   H1  NetImmune: HR per SD < 1 (protective) in each new cohort and pooled.
##   H2  CAF_ECM:   HR per SD > 1 (adverse) in each new cohort and pooled.
##   H3  TAM_SUP:   HR near 1 (95% CI including 1) pooled.
##   H4  Immune modules (APM, IFN, TLS_B, T_CELL) and the Immune axis individually weak
##       (pooled 95% CI including 1, or at least weaker than NetImmune).
##   Reported without expectation: Barrier, CXCL9:SPP1, ESTIMATE, CD8A.
##
## LOCKED METHODS (unchanged from the manuscript):
##   - probes -> gene symbol by highest mean expression (GPL570: GSE62254 annotation, first symbol
##     of '///' lists, as in B1_prep_cohorts.R; GPL8432: GEO platform table Symbol column);
##   - stored scale checked per cohort: if max > 50 the matrix is treated as linear and log2(x+1)
##     is applied once, otherwise it is used as provided;
##   - genes row-z within cohort; module = mean z of its genes present (>= 3 required);
##     Barrier = mean(CAF_ECM, TAM_SUP); Immune = mean(APM, IFN, TLS_B, T_CELL);
##     NetImmune = Immune - Barrier; CXCL9:SPP1 = z(CXCL9) - z(SPP1);
##   - Cox HR per SD (score standardised within cohort); adjusted = + stage (ordinal I-IV) + age
##     where recorded;
##   - pooling by random effects (metafor), reporting the wider of the Hartung-Knapp and
##     DerSimonian-Laird intervals and the larger P;
##   - structure-matched random-gene null for NetImmune: 1,000 axes (4 modules of 8/5/6/7 genes
##     minus 2 of 11) drawn from genes measured in every new cohort, the same gene set in every
##     cohort, credited with its better direction once (on the inverse-variance pooled log HR);
##   - Harrell's C with direction fixed in advance; paired bootstrap (2,000) of Delta C within
##     cohort, pooled by inverse variance.
## Cohorts with no deposited outcome are scored but cannot enter survival analyses; they are
## listed with the reason.
suppressMessages({library(survival); library(metafor); library(matrixStats); library(readxl)})
set.seed(2026)
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
D <- file.path(G, "2/data_download/geo")
O <- file.path(G, "2/analysis/results")
R <- file.path(G, "revision/results")
LOG <- file.path(O, "E15_log.txt"); sink(LOG, split = TRUE)

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

## ---- annotations ----
f570 <- read.csv(file.path(G, "GC_TME_Reprogrammer/data_GEO/GSE62254/GSE62254_fdata.csv"), stringsAsFactors = FALSE)
sym570 <- setNames(sub(" ///.*", "", f570$Gene.Symbol), f570$ID)
a8432 <- read.delim(file.path(D, "platforms/GPL8432_id_symbol.tsv"), stringsAsFactors = FALSE)
sym8432 <- setNames(a8432$Symbol, a8432$ID)

prep <- function(g, symmap) {
  x <- read_matrix(g); m <- x$m
  rng <- range(m, na.rm = TRUE)
  scale_note <- if (rng[2] > 50) "linear (max > 50): log2(x+1) applied once" else "log2 as deposited: used as provided"
  if (rng[2] > 50) m <- log2(pmax(m, 0) + 1)
  m <- m[rowSums(is.na(m)) == 0, , drop = FALSE]
  m <- collapse_max(m, unname(symmap[rownames(m)]))
  cat(sprintf("%-9s stored range %.2f..%.2f -> %s; %d genes x %d samples\n", g, rng[1], rng[2],
              scale_note, nrow(m), ncol(m)))
  list(m = m, meta = x$meta, scale = scale_note, range = rng)
}

## ---- cohorts with outcome ----
cat("== loading cohorts\n")
X <- list()
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
  X[[g]] <- list(m = p$m[, cl$sample], clin = cl, platform = "GPL570 (U133 Plus 2.0)", scale = p$scale)
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
X[["GSE26253"]] <- list(m = p$m[, cl$sample], clin = cl, platform = "GPL8432 (DASL, FFPE)", scale = p$scale)

## ---- cohorts without deposited outcome (scored only; listed) ----
NOOUT <- data.frame(
  cohort = c("GSE35809", "GSE51105", "GSE26901", "GSE13861", "GSE29272", "GSE57303", "GSE22377", "GSE14210"),
  reason = c("no survival deposited (subtype only); same Australian arrays as GSE51105",
             "no survival deposited; CEL files byte-identical in size to GSE35809 (same patients)",
             "clinical files give age, Lauren, AJCC stage and adjuvant chemotherapy only; no survival",
             "clinical file gives TNM, stage, Lauren and grade only; no survival",
             "tumour/normal pairs; no outcome deposited",
             "no outcome deposited",
             "Lauren type only; no outcome deposited",
             "chemotherapy-response labels not deposited (training/validation set labels only)"),
  stringsAsFactors = FALSE)

## ---- overlap checks ----
cat("\n== overlap checks\n")
ids <- lapply(X, function(x) x$clin$pid)
cat("patient-ID overlap GSE15459 vs GSE34942:", length(intersect(ids$GSE15459, ids$GSE34942)), "\n")
O1 <- readRDS(file.path(R, "B1_GSE62254.rds")); O2 <- readRDS(file.path(R, "B1_GSE84437.rds"))
best_match <- function(A, B, label) {
  g <- intersect(rownames(A), rownames(B))
  v <- rowVars(cbind(zrow(A[g, ]), zrow(B[g, ]))); g <- g[order(-v)][1:min(2000, length(g))]
  C <- cor(zrow(A[g, ]), zrow(B[g, ]), method = "spearman")
  top <- apply(C, 1, function(r) sort(r, decreasing = TRUE)[1:2])
  d <- data.frame(pair = label, n_common_genes = length(g), max_r = max(top[1, ]),
                  median_best_r = median(top[1, ]), max_gap = max(top[1, ] - top[2, ]),
                  n_gap_gt_0.10 = sum(top[1, ] - top[2, ] > 0.10))
  print(d); d
}
OV <- rbind(best_match(X$GSE15459$m, X$GSE34942$m, "GSE15459 vs GSE34942"),
            best_match(X$GSE15459$m, O1$expr, "GSE15459 vs GSE62254"),
            best_match(X$GSE34942$m, O1$expr, "GSE34942 vs GSE62254"),
            best_match(X$GSE26253$m, O1$expr, "GSE26253 vs GSE62254"),
            best_match(X$GSE26253$m, O2$expr, "GSE26253 vs GSE84437"))
write.csv(OV, file.path(O, "E15_overlap_check.csv"), row.names = FALSE)

## ---- scores ----
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
S <- Map(score, X, names(X))
write.csv(do.call(rbind, COV), file.path(O, "E15_module_coverage.csv"), row.names = FALSE)
cat("\n== module coverage\n"); print(do.call(rbind, COV)[, c("cohort","module","n_genes","of","missing")], row.names = FALSE)

## ---- Cox per cohort ----
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
PC <- do.call(rbind, lapply(S, function(d) do.call(rbind, c(lapply(SCORES, cox1, d = d, adj = FALSE),
                                                             lapply(SCORES, cox1, d = d, adj = TRUE)))))
write.csv(PC, file.path(O, "E15_per_cohort.csv"), row.names = FALSE)
cat("\n== per-cohort HR per SD (unadjusted)\n")
print(transform(PC[PC$model == "unadjusted" & PC$score %in% c("NetImmune","CAF_ECM","TAM_SUP","Immune"),
                   c("cohort","score","n","events","endpoint","HR","lo","hi","p")], HR = round(HR, 2), lo = round(lo, 2),
                hi = round(hi, 2), p = signif(p, 2)), row.names = FALSE)

## ---- pooling ----
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
orig$score[orig$score == "CXCL9_SPP1"] <- "CXCL9_SPP1"
PL <- list()
for (s in SCORES) {
  for (m in c("unadjusted", "adjusted")) {
    a <- PC[PC$score == s & PC$model == m, ]
    if (nrow(a) >= 2) PL[[length(PL) + 1]] <- cbind(set = "new cohorts", score = s, model = m, pool(a$beta, a$se))
    a2 <- a[a$endpoint == "OS", ]
    if (nrow(a2) >= 2) PL[[length(PL) + 1]] <- cbind(set = "new cohorts, OS only", score = s, model = m, pool(a2$beta, a2$se))
  }
  a <- PC[PC$score == s & PC$model == "unadjusted", ]; o <- orig[orig$score == s, ]
  if (nrow(a) >= 1 && nrow(o) >= 1)
    PL[[length(PL) + 1]] <- cbind(set = "new + original 3", score = s, model = "unadjusted",
                                  pool(c(a$beta, o$beta), c(a$se, o$se)))
  if (nrow(o) >= 2) PL[[length(PL) + 1]] <- cbind(set = "original 3 (reference)", score = s, model = "unadjusted", pool(o$beta, o$se))
}
PL <- do.call(rbind, PL)
write.csv(PL, file.path(O, "E15_pooled.csv"), row.names = FALSE)
cat("\n== pooled (conservative interval)\n")
print(transform(PL[PL$score %in% c("NetImmune","CAF_ECM","TAM_SUP","Barrier","Immune","IFN","CXCL9_SPP1","ESTIMATE","CD8A","APM","TLS_B","T_CELL"),
                   c("set","score","model","k","HR","lo","hi","p","interval","I2")],
                HR = round(HR, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 3), I2 = round(I2)), row.names = FALSE)

## ---- structure-matched random-gene null (NetImmune, new cohorts, unadjusted) ----
cat("\n== random-gene null\n")
Zs <- lapply(X, function(x) zrow(x$m))
gpool <- Reduce(intersect, lapply(Zs, rownames))
gpool <- gpool[sapply(gpool, function(g) all(sapply(Zs, function(z) sd(z[g, ]) > 0)))]
cat("null gene pool (measured in all new cohorts):", length(gpool), "\n")
SIZES <- c(8, 5, 6, 7, 11, 11)
obs <- PC[PC$score == "NetImmune" & PC$model == "unadjusted", ]
fitb <- function(d, x) { co <- summary(coxph(Surv(d$time, d$event) ~ as.numeric(scale(x))))$coefficients; co[1, c(1, 3)] }
iv <- function(b, se) { w <- 1 / se^2; sum(w * b) / sum(w) }
obs_pooled <- iv(obs$beta, obs$se)
NR <- 1000; draws <- matrix(NA, NR, length(X) + 1, dimnames = list(NULL, c(names(X), "pooled")))
for (r in 1:NR) {
  gs <- split(sample(gpool, sum(SIZES)), rep(1:6, SIZES))
  bs <- ses <- numeric(length(X))
  for (i in seq_along(X)) {
    z <- Zs[[i]]; mods <- sapply(gs, function(g) colMeans(z[g, , drop = FALSE]))
    ax <- rowMeans(mods[, 1:4]) - rowMeans(mods[, 5:6])
    ax <- ax[match(S[[i]]$sample, colnames(z))]
    fb <- fitb(S[[i]], ax); bs[i] <- fb[1]; ses[i] <- fb[2]
  }
  pb <- iv(bs, ses)                                        # pooled: direction credited once, on the pooled log HR
  draws[r, ] <- c(abs(bs), abs(pb))                        # per cohort: better direction within the cohort
}
nullp <- data.frame(level = c(names(X), "pooled (inverse variance)"),
                    observed_logHR = c(obs$beta, obs_pooled),
                    null_median_abs = c(apply(abs(draws[, names(X)]), 2, median), median(draws[, "pooled"])),
                    null_p95_abs = c(apply(abs(draws[, names(X)]), 2, quantile, 0.95), quantile(draws[, "pooled"], 0.95)),
                    rank = c(sapply(names(X), function(k) 1 + sum(abs(draws[, k]) >= abs(obs$beta[obs$cohort == k]))),
                             1 + sum(draws[, "pooled"] >= abs(obs_pooled))))
nullp$P <- nullp$rank / (NR + 1)
write.csv(nullp, file.path(O, "E15_null.csv"), row.names = FALSE)
saveRDS(draws, file.path(O, "E15_null_draws.rds"))
print(nullp)
cat("(per-cohort ranks credit each random axis with its better direction within the cohort)\n")

## ---- Harrell's C, paired bootstrap Delta C ----
cat("\n== Harrell's C (direction fixed: higher NetImmune/ESTIMATE/CD8A = lower risk)\n")
cidx <- function(t, e, risk) concordance(Surv(t, e) ~ risk, reverse = TRUE)$concordance
CC <- list(); DC <- list()
for (k in names(S)) {
  d <- S[[k]]; d <- d[complete.cases(d[, c("time", "event", "NetImmune", "ESTIMATE", "CD8A")]), ]
  cs <- sapply(c("NetImmune", "ESTIMATE", "CD8A"), function(s) cidx(d$time, d$event, -d[[s]]))
  CC[[k]] <- data.frame(cohort = k, score = names(cs), C = cs)
  B <- replicate(2000, { i <- sample(nrow(d), replace = TRUE); dd <- d[i, ]
    cn <- cidx(dd$time, dd$event, -dd$NetImmune)
    c(cn - cidx(dd$time, dd$event, -dd$ESTIMATE), cn - cidx(dd$time, dd$event, -dd$CD8A)) })
  for (j in 1:2) DC[[paste(k, j)]] <- data.frame(cohort = k, comparator = c("ESTIMATE", "CD8A")[j],
        dC = cs[1] - cs[j + 1], se = sd(B[j, ]), lo = quantile(B[j, ], 0.025), hi = quantile(B[j, ], 0.975))
}
CC <- do.call(rbind, CC); DC <- do.call(rbind, DC)
for (cmp in c("ESTIMATE", "CD8A")) {
  a <- DC[DC$comparator == cmp, ]; w <- 1 / a$se^2; est <- sum(w * a$dC) / sum(w); se <- sqrt(1 / sum(w))
  DC <- rbind(DC, data.frame(cohort = "pooled (inverse variance)", comparator = cmp, dC = est, se = se,
                             lo = est - 1.96 * se, hi = est + 1.96 * se))
}
DC$p <- 2 * pnorm(-abs(DC$dC / DC$se))
write.csv(CC, file.path(O, "E15_cindex.csv"), row.names = FALSE)
write.csv(DC, file.path(O, "E15_deltaC.csv"), row.names = FALSE)
print(CC, row.names = FALSE); print(transform(DC, dC = round(dC, 4), lo = round(lo, 4), hi = round(hi, 4), p = signif(p, 2)), row.names = FALSE)

## ---- cohort table ----
CT <- do.call(rbind, lapply(names(X), function(k) data.frame(cohort = k, platform = X[[k]]$platform,
        scale = X[[k]]$scale, n = nrow(S[[k]]), events = sum(S[[k]]$event), endpoint = S[[k]]$endpoint[1],
        stage = sum(!is.na(S[[k]]$stage_ord)), age = sum(!is.na(S[[k]]$age)), status = "analysed")))
CT <- rbind(CT, data.frame(cohort = NOOUT$cohort, platform = NA, scale = NA, n = NA, events = NA, endpoint = NA,
                           stage = NA, age = NA, status = paste("excluded:", NOOUT$reason)))
write.csv(CT, file.path(O, "E15_cohort_table.csv"), row.names = FALSE)
cat("\n== cohort table\n"); print(CT, row.names = FALSE)
saveRDS(S, file.path(O, "E15_scores.rds"))
sink()
