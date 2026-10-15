## E11 - prognosis outside gastric cancer.
##
## Question. In gastric cancer CAF_ECM was the strongest prognostic module, but across 11
## checkpoint cohorts from other cancer types (E5) it did not track response. Prognosis
## and response were therefore measured in different cancer types, which confounds the
## comparison. Here the same scores are tested for prognosis (untreated TCGA, OS) in the
## cancer types of the checkpoint cohorts.
##
## Pre-stated expectations (written before any outcome below was computed):
##   CAF_ECM  - if stromal content is prognostic in general, HR per SD > 1 (adverse) pooled
##              over group (i) SKCM, KIRC, LUAD, LUSC, BLCA, and over group (ii) ESCA + STAD,
##              even though it did not predict checkpoint response in those types
##              (E5 pooled log-OR +0.02, 95% CI -0.19 to 0.23).
##   Immune, IFN, NetImmune - HR < 1 (protective) pooled over group (i) and group (ii);
##              KIRC is the expected exception (immune infiltration is adverse in ccRCC).
##   TAM_SUP, C1Q_arm, SPP1_arm - no expected direction (null in gastric cancer).
##   CXCL9_SPP1 - HR < 1 (protective), as it trended in gastric cancer.
##   Dissociation - CAF_ECM carries adverse prognostic information in group (i) with a null
##              response effect; Immune is protective and response-associated.
##
## Data. TCGA pan-cancer RSEM (Xena, log2) for 24 types; TCGA-CDR endpoints [@liu2018cdr].
## Samples: primary tumours (-01), one per patient. SKCM is 78% metastases (-06), so the main
## SKCM analysis uses -01 and -06 (one sample per patient, -01 preferred) with sample type
## as an extra covariate in the adjusted model; SKCM primaries alone are a sensitivity row.
## Scoring as in the manuscript: genes row-z within cancer type (ESCA split: within
## histology), module = mean z of its genes present (>= 3), Barrier = mean(CAF_ECM,
## TAM_SUP), Immune = mean(APM, IFN, TLS_B, T_CELL), NetImmune = Immune - Barrier,
## CXCL9_SPP1 = z(CXCL9) - z(SPP1). MZB1 is stored under its old symbol MGC29506 and is
## renamed; TRAC and TRBC1 are absent from this annotation, so T_CELL uses 5 of 7 genes.
## Effects: Cox HR per SD of the score within cancer type; OS primary, PFI sensitivity;
## unadjusted and adjusted. Adjusted = age + AJCC pathologic stage (I-IV, factor); stage
## from clinical_stage for OV, CESC, UCEC (no AJCC pathologic stage recorded); age + grade
## for LGG; age only for GBM, PRAD, SARC. Stage levels with < 5 patients are merged into
## the adjacent lower level. A type enters a pooled estimate if the endpoint has >= 10 events.
## Pooling: DerSimonian-Laird random effects and Hartung-Knapp (metafor, test = "knha");
## reported interval/P = the wider interval and larger P of the two.
## Groups: (i) SKCM, KIRC, LUAD, LUSC, BLCA; (ii) ESCA adenocarcinoma, ESCA squamous, STAD;
## (iii) all 24 types (ESCA as one type).
## Null: 1,000 random 11-gene sets (CAF_ECM) and 1,000 structure-matched random Immune axes
## (four random modules with the gene counts used here), drawn from genes measured in all
## five group (i) types, the same set in every type, pooled by DL and credited with their
## better direction once (|pooled log HR|).
## METABRIC: the local file carries 489 genes, of which only FN1, STAT1 and STAT3 are module
## genes, and the cBioPortal archive is a failed 243-byte download (no network), so no module
## can be scored. TCGA-BRCA adjusted for age, stage and ESR1-inferred ER status is reported
## in its place.
suppressMessages({library(survival); library(matrixStats); library(metafor)})
set.seed(2611)
D <- "/groups/adv2105_gp/yichen/Yi/com/gastric/pancan"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results"
NRAND <- 1000; MINEV <- 10
MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"),
  C1Q_arm = c("C1QA","C1QB","C1QC"),
  SPP1_arm = c("SPP1","TREM2","MARCO","APOE","LGALS3"))
SC <- c("CAF_ECM","TAM_SUP","APM","IFN","TLS_B","T_CELL","Barrier","Immune","NetImmune",
        "C1Q_arm","SPP1_arm","CXCL9_SPP1")
G1 <- c("SKCM","KIRC","LUAD","LUSC","BLCA"); G2 <- c("ESCA_AC","ESCA_SCC","STAD")
EP <- list(OS = c("OS","OS.time"), PFI = c("PFI","PFI.time"))

zrow <- function(m) { s <- rowSds(m); s[s == 0 | is.na(s)] <- 1; (m - rowMeans(m)) / s }
score_all <- function(Z) {
  used <- c()
  m <- sapply(names(MOD), function(k) { g <- intersect(MOD[[k]], rownames(Z)); used[k] <<- length(g)
    if (length(g) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[g, , drop = FALSE]) })
  S <- as.data.frame(m); rownames(S) <- colnames(Z)
  S$Barrier <- rowMeans(S[, c("CAF_ECM","TAM_SUP")])
  S$Immune <- rowMeans(S[, c("APM","IFN","TLS_B","T_CELL")])
  S$NetImmune <- S$Immune - S$Barrier
  S$CXCL9_SPP1 <- Z["CXCL9", ] - Z["SPP1", ]
  attr(S, "used") <- used; S
}
merge_small <- function(f) {          # merge factor levels with < 5 patients into the level below
  f <- droplevels(f)
  repeat { tb <- table(f); if (length(tb) < 2 || all(tb >= 5)) break
    i <- which(tb < 5)[1]; lv <- levels(f); j <- if (i > 1) i - 1 else 2
    f <- factor(ifelse(as.character(f) == lv[i], lv[j], as.character(f)), levels = lv[-i]) }
  f
}
roman <- function(x) { x <- gsub("^Stage ", "", x); x <- gsub("[ABC][12]?$", "", x)
  factor(ifelse(x %in% c("I","II","III","IV"), x, NA), levels = c("I","II","III","IV")) }

surv <- read.delim(file.path(D, "Survival_SupplementalTable_S1_20171025_xena_sp"), check.names = FALSE,
                   stringsAsFactors = FALSE)
TYPES <- sub("\\.gz$", "", list.files(file.path(D, "expr"), pattern = "\\.gz$"))
stopifnot(length(TYPES) == 24)

## per-type units: expression Z, scores, clinical covariates
make_unit <- function(e, s0, tag) {
  e <- e[rowSums(is.finite(e)) == ncol(e) & rowSds(e) > 0, , drop = FALSE]
  Z <- zrow(e); S <- score_all(Z)
  ca <- s0$`cancer type abbreviation`[1]
  st <- if (ca %in% c("OV","CESC","UCEC")) roman(s0$clinical_stage) else roman(s0$ajcc_pathologic_tumor_stage)
  if (ca == "LGG") st <- factor(ifelse(s0$histological_grade %in% c("G2","G3"), s0$histological_grade, NA))
  if (ca %in% c("GBM","PRAD","SARC")) st <- factor(rep("all", nrow(s0)))
  cl <- data.frame(age = suppressWarnings(as.numeric(s0$age_at_initial_pathologic_diagnosis)), stage = st,
                   met = as.integer(substr(s0$sample, 14, 15) == "06"))
  for (k in names(EP)) { cl[[paste0(k, "_e")]] <- s0[[EP[[k]][1]]]; cl[[paste0(k, "_t")]] <- s0[[EP[[k]][2]]] }
  list(tag = tag, Z = Z, S = S, cl = cl, used = attr(S, "used"))
}
cox_eff <- function(v, cl, ep, adj, tag) {
  d <- cl; d$x <- as.numeric(scale(v)); d$e <- d[[paste0(ep, "_e")]]; d$t <- d[[paste0(ep, "_t")]]
  d <- d[!is.na(d$e) & !is.na(d$t) & d$t > 0 & is.finite(d$x), ]
  f <- "Surv(t, e) ~ x"
  if (adj) {
    d <- d[!is.na(d$age) & !is.na(d$stage), ]; d$stage <- merge_small(d$stage)
    f <- paste(f, "+ age"); if (nlevels(d$stage) > 1) f <- paste(f, "+ stage")
    if (tag == "SKCM" && length(unique(d$met)) > 1) f <- paste(f, "+ met")
  }
  if (sum(d$e) < 3) return(c(n = nrow(d), ev = sum(d$e), b = NA, se = NA, p = NA))
  m <- suppressWarnings(coxph(as.formula(f), data = d)); s <- summary(m)$coefficients
  c(n = nrow(d), ev = sum(d$e), b = s["x", "coef"], se = s["x", "se(coef)"], p = s["x", "Pr(>|z|)"])
}

U <- list(); used_tab <- list()
for (ca in TYPES) {
  e <- as.matrix(read.delim(gzfile(file.path(D, "expr", paste0(ca, ".gz"))), row.names = 1, check.names = FALSE))
  colnames(e) <- gsub("\\.", "-", colnames(e)); rownames(e)[rownames(e) == "MGC29506"] <- "MZB1"
  ty <- substr(colnames(e), 14, 15)
  keep_types <- if (ca == "SKCM") c("01", "06") else "01"
  e <- e[, ty %in% keep_types, drop = FALSE]
  e <- e[, order(substr(colnames(e), 14, 15)), drop = FALSE]           # -01 before -06
  e <- e[, !duplicated(substr(colnames(e), 1, 12)) & colnames(e) %in% surv$sample, drop = FALSE]
  s0 <- surv[match(colnames(e), surv$sample), ]
  if (ca == "ESCA") {
    for (h in c("AC", "SCC")) {
      w <- s0$histological_type %in% c(AC = "Esophagus Adenocarcinoma, NOS", SCC = "Esophagus Squamous Cell Carcinoma")[[h]]
      U[[paste0("ESCA_", h)]] <- make_unit(e[, w, drop = FALSE], s0[w, ], paste0("ESCA_", h))
    }
  }
  if (ca == "SKCM") { w <- substr(colnames(e), 14, 15) == "01"
    U[["SKCM_primary"]] <- make_unit(e[, w, drop = FALSE], s0[w, ], "SKCM_primary") }
  U[[ca]] <- make_unit(e, s0, ca)
  for (u in setdiff(names(U), G1)) U[[u]]$Z <- NULL                     # expression kept for the null only
  cat(sprintf("%-6s n=%4d  genes used: %s\n", ca, ncol(e), paste(names(U[[ca]]$used), U[[ca]]$used, sep = "=", collapse = " ")))
  rm(e); gc(verbose = FALSE)
}
used_tab <- do.call(rbind, lapply(names(U), function(u) data.frame(unit = u, module = names(U[[u]]$used), genes = U[[u]]$used)))

## ---------------- per-type effects ----------------
rows <- list()
for (u in names(U)) for (k in SC) for (ep in names(EP)) for (adj in c(FALSE, TRUE)) {
  r <- cox_eff(U[[u]]$S[[k]], U[[u]]$cl, ep, adj, u)
  rows[[length(rows) + 1]] <- data.frame(unit = u, score = k, endpoint = ep, model = ifelse(adj, "adjusted", "unadjusted"),
    n = r[["n"]], events = r[["ev"]], logHR = r[["b"]], se = r[["se"]], p = r[["p"]])
}
per <- do.call(rbind, rows)
per$HR <- exp(per$logHR); per$lo <- exp(per$logHR - 1.96 * per$se); per$hi <- exp(per$logHR + 1.96 * per$se)
write.csv(per, file.path(O, "E11_per_type.csv"), row.names = FALSE)
write.csv(used_tab, file.path(O, "E11_genes_used.csv"), row.names = FALSE)

## ---------------- pooled by group ----------------
GROUPS <- list(`(i) SKCM KIRC LUAD LUSC BLCA` = G1, `(ii) ESCA-AC ESCA-SCC STAD` = G2,
               `(ii') ESCA STAD` = c("ESCA", "STAD"), `(iii) all 24` = TYPES,
               `(i) with SKCM primaries only` = c("SKCM_primary", G1[-1]))
pool <- function(u) {
  dl <- rma(yi = u$logHR, sei = u$se, method = "DL"); hk <- rma(yi = u$logHR, sei = u$se, method = "DL", test = "knha")
  wide_hk <- (hk$ci.ub - hk$ci.lb) >= (dl$ci.ub - dl$ci.lb)
  data.frame(k = nrow(u), n = sum(u$n), events = sum(u$events), HR = exp(dl$b[1]),
    DL_lo = exp(dl$ci.lb), DL_hi = exp(dl$ci.ub), DL_p = dl$pval, HK_lo = exp(hk$ci.lb), HK_hi = exp(hk$ci.ub), HK_p = hk$pval,
    lo = exp(ifelse(wide_hk, hk$ci.lb, dl$ci.lb)), hi = exp(ifelse(wide_hk, hk$ci.ub, dl$ci.ub)), p = max(dl$pval, hk$pval),
    I2 = dl$I2, n_adverse = sum(u$logHR > 0))
}
pooled <- do.call(rbind, lapply(names(GROUPS), function(g) do.call(rbind, lapply(SC, function(k)
  do.call(rbind, lapply(names(EP), function(ep) do.call(rbind, lapply(c("unadjusted", "adjusted"), function(md) {
    u <- per[per$unit %in% GROUPS[[g]] & per$score == k & per$endpoint == ep & per$model == md &
             per$events >= MINEV & is.finite(per$se), ]
    cbind(data.frame(group = g, score = k, endpoint = ep, model = md), pool(u)) }))))))))
write.csv(pooled, file.path(O, "E11_pooled.csv"), row.names = FALSE)
cat("\n=== pooled HR per SD, OS ===\n")
print(format(pooled[pooled$endpoint == "OS" & pooled$score %in% c("CAF_ECM","TAM_SUP","IFN","Immune","Barrier","NetImmune","CXCL9_SPP1"),
  c("group","score","model","k","HR","lo","hi","p","I2","n_adverse")], digits = 3), row.names = FALSE)

## ---------------- random-gene null, group (i), OS ----------------
Zs <- lapply(G1, function(u) U[[u]]$Z); names(Zs) <- G1
pool_g <- Reduce(intersect, lapply(Zs, rownames))
nIm <- sapply(c("APM","IFN","TLS_B","T_CELL"), function(k) U[["SKCM"]]$used[[k]])
stopifnot(all(sapply(G1, function(u) all(U[[u]]$used[c("APM","IFN","TLS_B","T_CELL")] == nIm))))
dl_b <- function(b, se) { w <- 1 / se^2; fe <- sum(w * b) / sum(w); Q <- sum(w * (b - fe)^2); k <- length(b)
  t2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))); ws <- 1 / (se^2 + t2); sum(ws * b) / sum(ws) }
pooled_b <- function(scores, adj) { r <- sapply(G1, function(u) cox_eff(scores[[u]], U[[u]]$cl, "OS", adj, u)[c("b","se")])
  dl_b(r["b", ], r["se", ]) }
obs <- expand.grid(score = c("CAF_ECM", "Immune"), model = c("unadjusted", "adjusted"), stringsAsFactors = FALSE)
obs$b <- mapply(function(k, md) pooled_b(setNames(lapply(G1, function(u) U[[u]]$S[[k]]), G1), md == "adjusted"), obs$score, obs$model)
draws <- matrix(NA, NRAND, 4, dimnames = list(NULL, paste(obs$score, obs$model, sep = "_")))
for (i in seq_len(NRAND)) {
  g11 <- sample(pool_g, 11); gm <- lapply(nIm, function(n) sample(pool_g, n))
  r11 <- setNames(lapply(G1, function(u) colMeans(Zs[[u]][g11, ])), G1)
  rim <- setNames(lapply(G1, function(u) rowMeans(sapply(gm, function(g) colMeans(Zs[[u]][g, , drop = FALSE])))), G1)
  draws[i, ] <- c(pooled_b(r11, FALSE), pooled_b(rim, FALSE), pooled_b(r11, TRUE), pooled_b(rim, TRUE))
  if (i %% 250 == 0) cat("null draw", i, "\n")
}
obs$null_med <- apply(abs(draws), 2, median); obs$null_p95 <- apply(abs(draws), 2, quantile, 0.95)
obs$rank <- sapply(1:4, function(j) sum(abs(draws[, j]) >= abs(obs$b[j])) + 1)
obs$P_null <- obs$rank / (NRAND + 1)
obs$frac_null_adverse <- apply(draws, 2, function(x) mean(x > 0))
write.csv(obs, file.path(O, "E11_null.csv"), row.names = FALSE)
write.csv(draws, file.path(O, "E11_null_draws.csv"), row.names = FALSE)
cat("\n=== random-gene null, group (i) pooled OS ===\n"); print(format(obs, digits = 3), row.names = FALSE)

## ---------------- dissociation: prognosis (group i, TCGA) vs response (E5, 11 ICB cohorts) ----------------
resp <- read.csv(file.path(O, "E5_icb_pooled.csv")); resp <- resp[resp$set == "all 11", ]
rpc <- read.csv(file.path(O, "E5_icb_per_cohort.csv"))
pg <- pooled[pooled$group == names(GROUPS)[1] & pooled$endpoint == "OS", ]
dis <- do.call(rbind, lapply(SC, function(k) { a <- pg[pg$score == k & pg$model == "adjusted", ]; r <- resp[resp$score == k, ]
  data.frame(score = k, prog_logHR = log(a$HR), prog_lo = log(a$lo), prog_hi = log(a$hi), prog_p = a$p,
             resp_logOR = r$logOR, resp_lo = r$lo, resp_hi = r$hi, resp_p = r$p) }))
TT <- c(melanoma = "SKCM", kidney = "KIRC", lung = "LUAD+LUSC", bladder = "BLCA")
tum <- do.call(rbind, lapply(names(TT), function(t) do.call(rbind, lapply(SC, function(k) {
  rr <- rpc[rpc$tumour == t & rpc$score == k, ]
  rm_ <- if (nrow(rr) > 1) rma(yi = rr$logOR, sei = rr$se, method = "DL") else NULL
  units <- strsplit(TT[[t]], "\\+")[[1]]
  pp <- per[per$unit %in% units & per$score == k & per$endpoint == "OS" & per$model == "adjusted", ]
  pm <- if (nrow(pp) > 1) rma(yi = pp$logHR, sei = pp$se, method = "DL") else NULL
  data.frame(tumour = t, tcga = TT[[t]], score = k, n_icb = nrow(rr),
             resp_logOR = if (is.null(rm_)) rr$logOR else rm_$b[1], resp_se = if (is.null(rm_)) rr$se else rm_$se,
             prog_logHR = if (is.null(pm)) pp$logHR else pm$b[1], prog_se = if (is.null(pm)) pp$se else pm$se) }))))
cor_s <- cor.test(dis$prog_logHR, dis$resp_logOR, method = "spearman", exact = FALSE)
write.csv(dis, file.path(O, "E11_dissociation_scores.csv"), row.names = FALSE)
write.csv(tum, file.path(O, "E11_dissociation_by_tumour.csv"), row.names = FALSE)
cat("\n=== prognosis (group i, adjusted OS) vs response (11 ICB cohorts) ===\n"); print(format(dis, digits = 3), row.names = FALSE)
cat(sprintf("Spearman rho across the 12 scores (log HR vs log OR): %.2f, P = %.3f\n", cor_s$estimate, cor_s$p.value))
cat("\n=== by tumour (TCGA adjusted OS log HR | ICB log OR) ===\n"); print(format(tum[tum$score %in% c("CAF_ECM","Immune","IFN","TAM_SUP"), ], digits = 3), row.names = FALSE)

## ---------------- METABRIC substitute: TCGA-BRCA with ER status inferred from ESR1 ----------------
b <- as.matrix(read.delim(gzfile(file.path(D, "expr", "BRCA.gz")), row.names = 1, check.names = FALSE))
colnames(b) <- gsub("\\.", "-", colnames(b)); esr1 <- b["ESR1", colnames(b) %in% rownames(U$BRCA$S)]; rm(b)
cut <- { dd <- density(esr1); lm_ <- which(diff(sign(diff(dd$y))) == 2) + 1
  lm_ <- lm_[dd$x[lm_] > quantile(esr1, 0.05) & dd$x[lm_] < quantile(esr1, 0.6)]
  if (length(lm_)) dd$x[lm_[which.min(dd$y[lm_])]] else median(esr1) }
cl <- U$BRCA$cl; cl$er <- as.integer(esr1[rownames(U$BRCA$S)] > cut)
brca <- do.call(rbind, lapply(SC, function(k) do.call(rbind, lapply(c("OS", "DSS"), function(ep) {
  d <- cl; d$x <- as.numeric(scale(U$BRCA$S[[k]]))
  d$e <- surv$DSS[match(rownames(U$BRCA$S), surv$sample)]; d$t <- surv$DSS.time[match(rownames(U$BRCA$S), surv$sample)]
  if (ep == "OS") { d$e <- d$OS_e; d$t <- d$OS_t }
  d <- d[complete.cases(d[, c("e", "t", "x", "age", "stage", "er")]) & d$t > 0, ]; d$stage <- merge_small(d$stage)
  m <- coxph(Surv(t, e) ~ x + age + stage + er, data = d); s <- summary(m)$coefficients
  data.frame(score = k, endpoint = ep, n = nrow(d), events = sum(d$e), HR = exp(s["x", 1]),
             lo = exp(s["x", 1] - 1.96 * s["x", 3]), hi = exp(s["x", 1] + 1.96 * s["x", 3]), p = s["x", 5]) }))))
attr(brca, "er_cut") <- cut
write.csv(brca, file.path(O, "E11_brca_substitute.csv"), row.names = FALSE)
cat(sprintf("\n=== TCGA-BRCA (METABRIC substitute), adjusted age + stage + ESR1-inferred ER (cut %.2f; ER+ %.0f%%) ===\n",
            cut, 100 * mean(cl$er, na.rm = TRUE)))
print(format(brca, digits = 3), row.names = FALSE)

## sanity: STAD from this pan-cancer matrix vs the manuscript's TCGA-STAD (CAF_ECM log HR 0.248)
cat(sprintf("\nSTAD CAF_ECM unadjusted OS log HR here %.3f (manuscript 0.248)\n",
            per$logHR[per$unit == "STAD" & per$score == "CAF_ECM" & per$endpoint == "OS" & per$model == "unadjusted"]))
cat("DONE\n")
