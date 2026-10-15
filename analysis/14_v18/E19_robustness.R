## E19 - robustness of the barriernet scores to gene loss and platform (pre-registered:
## E17_E22_prereg.md, section E19; that file is not edited here - deviations are recorded in
## results/E19_README.md).
##
## Data (loaders re-used, not re-implemented):
##   - gastric benchmark (overall survival): TCGA-STAD, GSE84437, GSE62254 from
##     revision/results/B1_<tag>.rds, scored as revision/scripts/B2_scores.R and pooled as
##     E8_prognostic_hardening.R (Cox HR per SD, DerSimonian-Laird and Hartung-Knapp);
##   - PRJEB25780 (pembrolizumab, 44 patients): as E7a_scores_log2.py - the deposited ".TPM"
##     matrix is already log2(TPM + 1) and is used as is; gene z with population SD;
##   - 21 checkpoint cohorts of E12 set "all": the 11 E5 cohorts loaded exactly as
##     E5_icb_split.R (MultiAssayExperiment objects, Ensembl rows, symbol -> first Ensembl id
##     present) and the 10 E12 replication cohorts from results/E12_cohorts.rds scored as
##     E12_icb_expanded.R (scores_sym); inclusion = E12 set "new" in E12_icb_per_cohort.csv.
##     Effect = logistic log-odds of response per SD of the score.
## Scoring (fixed): genes z-scored within cohort on the log2 scale; module = mean z of its genes
## present, >= 3 genes required; Barrier = mean(CAF_ECM, TAM_SUP); Immune = mean(APM, IFN,
## TLS_B, T_CELL); NetImmune = Immune - Barrier; scores standardised to SD 1 within cohort.
## Because a gene's z-score does not depend on which other genes are kept, each cohort's
## 48-gene z matrix is computed once and perturbations only change which rows are averaged.
## Pooling: random effects (metafor); reported interval = the wider of Hartung-Knapp and
## DerSimonian-Laird, reported P = the larger of the two.
##
## Perturbations:
##   - 1,000 random drops each of 10%, 20%, 30% of the 48 genes (5, 10, 14 genes); the same
##     dropped set is applied to every cohort in a draw; a draw is rejected and redrawn if any
##     module that is scorable in a cohort with the full list would fall below 3 genes there;
##   - leave-one-gene-out (48 runs; no redraw possible: a module that falls below 3 genes in a
##     cohort becomes unscorable there and that cohort leaves the affected estimates);
##   - panel subset: the module genes measured in all four targeted-panel replication cohorts of
##     E12 (tech == "panel"), plus post-hoc platform sensitivities (see README).
## Outputs: results/E19_*.csv, results/E19_log.txt.
suppressMessages({library(survival); library(metafor); library(matrixStats); library(parallel)
                  library(MultiAssayExperiment); library(SummarizedExperiment); library(org.Hs.eg.db)})
set.seed(2026)
G <- "/groups/adv2105_gp/yichen/Yi/com/gastric"
O <- file.path(G, "2/analysis/results")
W <- Sys.getenv("E19_OUT", O)          # output directory (results/ unless overridden for a test run)
R <- file.path(G, "revision/results")
ICB <- file.path(G, "icb")
KIM <- file.path(G, "GC_TME_Reprogrammer/validation")
NDRAW <- as.integer(Sys.getenv("E19_NDRAW", "1000")); FRACS <- c(0.1, 0.2, 0.3); NCORES <- 40
sink(file.path(W, "E19_log.txt"), split = TRUE)

MOD <- list(
  CAF_ECM = c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
  TAM_SUP = c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
  APM = c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
  IFN = c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
  TLS_B = c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
  T_CELL = c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"))
GENES <- unlist(MOD, use.names = FALSE); stopifnot(length(GENES) == 48, !anyDuplicated(GENES))
MOD_OF <- setNames(rep(names(MOD), lengths(MOD)), GENES)

## =========================== loaders ===========================
COH <- list()   # each: Z (48-gene z matrix, symbol rows, present genes only), group, outcome
## --- gastric benchmark (B1/B2/E8)
zrow <- function(m) { s <- rowSds(m, na.rm = TRUE); s[s == 0 | is.na(s)] <- 1
                      (m - rowMeans(m, na.rm = TRUE)) / s }
BENCH <- c("TCGA", "GSE84437", "GSE62254")
for (t in BENCH) {
  d <- readRDS(file.path(R, paste0("B1_", t, ".rds")))
  stopifnot(identical(d$clin$sample, colnames(d$expr)))
  COH[[t]] <- list(Z = zrow(d$expr[intersect(GENES, rownames(d$expr)), , drop = FALSE]), group = "benchmark",
                   time = d$clin$time, event = d$clin$event)
}
## --- PRJEB25780 (E7a): matrix already log2(TPM + 1); population SD, zero SD -> 1
b <- read.delim(file.path(KIM, "benchmark_scores.tsv"), row.names = 1, check.names = FALSE)
L <- as.matrix(read.delim(file.path(KIM, "mGC_Pembrolizumab_Kim2018.TPM"), row.names = 1, check.names = FALSE)[, rownames(b)])
L <- L[rowSums(is.finite(L)) == ncol(L), , drop = FALSE]
L <- L[intersect(GENES, rownames(L)), , drop = FALSE]
sd0 <- sqrt(rowMeans((L - rowMeans(L))^2)); sd0[sd0 == 0] <- 1
COH[["PRJEB25780"]] <- list(Z = (L - rowMeans(L)) / sd0, group = "gastric_ici", y = as.integer(b$Responder))
## --- 11 E5 cohorts (E5_icb_split.R load_cohort / scores_of)
sym2ens <- suppressMessages(mapIds(org.Hs.eg.db, keys = GENES, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "list"))
load_cohort <- function(nm) {
  x <- suppressWarnings(readRDS(file.path(ICB, paste0(nm, ".rds"))))
  ex <- experiments(x); an <- names(ex)
  pick <- c("expr_gene_tpm", "expr")[c("expr_gene_tpm", "expr") %in% an][1]
  M <- as.matrix(assay(ex[[pick]])); cd <- as.data.frame(colData(x))
  r <- as.character(cd$response); y <- ifelse(r == "R", 1, ifelse(r == "NR", 0, NA)); names(y) <- rownames(cd)
  keep <- intersect(colnames(M), names(y)[!is.na(y)]); M <- M[, keep, drop = FALSE]; y <- y[keep]
  rownames(M) <- sub("\\..*$", "", rownames(M)); M <- M[!duplicated(rownames(M)), , drop = FALSE]
  M <- M[rowSums(is.finite(M)) == ncol(M) & rowSds(M) > 0, , drop = FALSE]
  ids <- sapply(GENES, function(g) { e <- intersect(sym2ens[[g]], rownames(M)); if (length(e)) e[1] else NA_character_ })
  ids <- ids[!is.na(ids)]
  Mg <- M[ids, , drop = FALSE]; rownames(Mg) <- names(ids)
  list(Z = (Mg - rowMeans(Mg)) / rowSds(Mg), y = unname(y))
}
E5C <- c("ICB_Braun", "ICB_Gide", "ICB_Hugo", "ICB_Jung", "ICB_Limagne1", "ICB_Liu", "ICB_Mariathasan",
         "ICB_Miao1", "ICB_Puch", "ICB_Riaz", "ICB_Van_Allen")
for (nm in E5C) { d <- load_cohort(nm); COH[[nm]] <- list(Z = d$Z, group = "icb", y = d$y) }
## --- 10 E12 replication cohorts (E12_cohorts.rds, scores_sym)
CO <- readRDS(file.path(O, "E12_cohorts.rds"))
aud <- read.csv(file.path(O, "E12_cohort_audit.csv"), stringsAsFactors = FALSE)
per12 <- read.csv(file.path(O, "E12_icb_per_cohort.csv"), stringsAsFactors = FALSE)
NEW <- unique(per12$cohort[per12$set == "new"]); stopifnot(length(NEW) == 10)
for (nm in NEW) {
  M <- CO[[nm]]$M; M <- M[rowSds(M) > 0, , drop = FALSE]; M <- M[intersect(GENES, rownames(M)), , drop = FALSE]
  COH[[nm]] <- list(Z = (M - rowMeans(M)) / rowSds(M), group = "icb", y = unname(CO[[nm]]$y))
}
ICBC <- c(E5C, NEW)
PANEL <- aud$id[aud$id %in% NEW & aud$tech == "panel"]
rm(CO); invisible(gc())

## gene presence (cohort x gene) and full-list module coverage
PRES <- t(sapply(COH, function(c) GENES %in% rownames(c$Z))); colnames(PRES) <- GENES
NFULL <- sapply(names(MOD), function(k) rowSums(PRES[, MOD[[k]], drop = FALSE]))
cat("cohorts:", length(COH), "(benchmark", length(BENCH), "+ PRJEB25780 + ICB", length(ICBC), ")\n")
cat("panel replication cohorts:", PANEL, "\n\nmodule genes present (full list):\n"); print(NFULL)

## =========================== scoring and effects ===========================
score_one <- function(Z, keep) {
  m <- sapply(names(MOD), function(k) {
    g <- intersect(intersect(MOD[[k]], keep), rownames(Z))
    if (length(g) < 3) rep(NA_real_, ncol(Z)) else colMeans(Z[g, , drop = FALSE]) })
  if (is.null(dim(m))) m <- matrix(m, nrow = 1, dimnames = list(NULL, names(MOD)))
  Barrier <- rowMeans(m[, c("CAF_ECM", "TAM_SUP"), drop = FALSE])
  Immune <- rowMeans(m[, c("APM", "IFN", "TLS_B", "T_CELL"), drop = FALSE])
  cbind(m, Barrier = Barrier, Immune = Immune, NetImmune = Immune - Barrier)
}
cox_eff <- function(time, event, x) {
  if (all(is.na(x))) return(c(b = NA, se = NA))
  f <- coxph(Surv(time, event) ~ scale(x)[, 1]); c(b = unname(coef(f)), se = sqrt(vcov(f)[1, 1]))
}
logit_eff <- function(y, s) {
  if (all(is.na(s))) return(c(b = NA, se = NA))
  s <- (s - mean(s)) / sd(s)
  f <- suppressWarnings(glm(y ~ s, family = binomial()))
  c(b = unname(coef(f)[2]), se = unname(sqrt(diag(vcov(f)))[2]))
}
auc <- function(y, v) { if (all(is.na(v))) return(NA_real_)
  r <- rank(v); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
pool <- function(b, se) {
  ok <- is.finite(b) & is.finite(se); b <- b[ok]; se <- se[ok]
  if (length(b) < 2) return(c(est = NA, lo_DL = NA, hi_DL = NA, p_DL = NA, lo_HK = NA, hi_HK = NA, p_HK = NA, lo = NA, hi = NA, p = NA, k = length(b)))
  m <- rma(yi = b, sei = se, method = "DL"); h <- rma(yi = b, sei = se, method = "DL", test = "knha")
  hkw <- (h$ci.ub - h$ci.lb) >= (m$ci.ub - m$ci.lb)
  c(est = m$b[1], lo_DL = m$ci.lb, hi_DL = m$ci.ub, p_DL = m$pval, lo_HK = h$ci.lb, hi_HK = h$ci.ub, p_HK = h$pval,
    lo = if (hkw) h$ci.lb else m$ci.lb, hi = if (hkw) h$ci.ub else m$ci.ub, p = max(m$pval, h$pval), k = length(b))
}

FULL <- lapply(COH, function(c) score_one(c$Z, GENES))
## one run = one kept gene set -> pooled estimates + per-cohort Spearman vs the full score
eval_keep <- function(keep) {
  S <- lapply(COH, function(c) score_one(c$Z, keep))
  hr <- sapply(BENCH, function(t) cox_eff(COH[[t]]$time, COH[[t]]$event, S[[t]][, "NetImmune"]))
  lo <- lapply(c(NetImmune = "NetImmune", Immune = "Immune", CAF_ECM = "CAF_ECM"), function(k)
    sapply(ICBC, function(nm) logit_eff(COH[[nm]]$y, S[[nm]][, k])))
  P <- list(HR_NetImmune = pool(hr["b", ], hr["se", ]),
            OR_Immune = pool(lo$Immune["b", ], lo$Immune["se", ]),
            OR_CAF_ECM = pool(lo$CAF_ECM["b", ], lo$CAF_ECM["se", ]),
            OR_NetImmune = pool(lo$NetImmune["b", ], lo$NetImmune["se", ]))
  out <- unlist(lapply(names(P), function(n) setNames(P[[n]], paste(n, names(P[[n]]), sep = "."))))
  out <- c(out, AUC_NetImmune_PRJEB = auc(COH$PRJEB25780$y, S$PRJEB25780[, "NetImmune"]),
           AUC_Immune_PRJEB = auc(COH$PRJEB25780$y, S$PRJEB25780[, "Immune"]),
           AUC_CAF_ECM_PRJEB = auc(COH$PRJEB25780$y, S$PRJEB25780[, "CAF_ECM"]))
  rho <- function(k) sapply(names(COH), function(nm) { a <- FULL[[nm]][, k]; p <- S[[nm]][, k]
    if (all(is.na(a)) || all(is.na(p))) NA_real_ else cor(a, p, method = "spearman") })
  list(est = out, rho = rho("NetImmune"), rho_imm = rho("Immune"), rho_caf = rho("CAF_ECM"),
       percohort = list(hr = hr, lo = lo))
}

## =========================== 0. reproduction of the published full-score values ===========================
full <- eval_keep(GENES); fe <- full$est
b2 <- readRDS(file.path(R, "B2_scores.rds")); e7 <- read.csv(file.path(O, "E7_scores_log2.csv"), check.names = FALSE)
e8 <- read.csv(file.path(O, "E8_meta_DL_HK.csv")); e12 <- read.csv(file.path(O, "E12_icb_pooled.csv"))
e12a <- e12[e12$set == "all", ]
chk_score <- c(sapply(BENCH, function(t) max(abs(FULL[[t]][, "NetImmune"] - b2[[t]]$NetImmune))),
               PRJEB25780 = max(abs(FULL$PRJEB25780[, "NetImmune"] - e7$NetImmune[match(rownames(b), e7$SEQ_ID)])))
pc <- per12[per12$cohort %in% ICBC & per12$score %in% c("Immune", "CAF_ECM", "NetImmune"), ]
pc_new <- mapply(function(nm, k) full$percohort$lo[[k]]["b", nm], pc$cohort, pc$score)
cat("\nmax |score difference| vs stored scores (B2_scores.rds; E7_scores_log2.csv):\n"); print(signif(chk_score, 3))
cat(sprintf("max |per-cohort log-OR difference| vs E12_icb_per_cohort.csv (%d cohort x score rows): %.2e\n",
            nrow(pc), max(abs(pc_new - pc$logOR))))
repro <- data.frame(
  quantity = c("pooled NetImmune HR per SD, benchmark 3 (OS)", "  wider interval low", "  wider interval high", "  larger P",
               "PRJEB25780 NetImmune AUC",
               "pooled Immune log-OR per SD, 21 ICB cohorts", "  wider interval low", "  wider interval high", "  k",
               "pooled CAF_ECM log-OR per SD, 21 ICB cohorts", "  wider interval low", "  wider interval high", "  k",
               "pooled NetImmune log-OR per SD, 21 ICB cohorts", "  k"),
  published = c(0.74, 0.58, 0.95, NA, 0.797, 0.43, NA, NA, 21, -0.06, NA, NA, 20, NA, 20),
  stored_result = c(e8$HR[e8$score == "NetImmune"], e8$HK_lo[e8$score == "NetImmune"], e8$HK_hi[e8$score == "NetImmune"],
                    max(e8$DL_p[e8$score == "NetImmune"], e8$HK_p[e8$score == "NetImmune"]), NA,
                    e12a$logOR[e12a$score == "Immune"], min(e12a$lo_HK[e12a$score == "Immune"], NA, na.rm = TRUE),
                    e12a$hi_HK[e12a$score == "Immune"], e12a$k[e12a$score == "Immune"],
                    e12a$logOR[e12a$score == "CAF_ECM"], e12a$lo_HK[e12a$score == "CAF_ECM"], e12a$hi_HK[e12a$score == "CAF_ECM"],
                    e12a$k[e12a$score == "CAF_ECM"], e12a$logOR[e12a$score == "NetImmune"], e12a$k[e12a$score == "NetImmune"]),
  reproduced = c(exp(fe[["HR_NetImmune.est"]]), exp(fe[["HR_NetImmune.lo"]]), exp(fe[["HR_NetImmune.hi"]]), fe[["HR_NetImmune.p"]],
                 fe[["AUC_NetImmune_PRJEB"]],
                 fe[["OR_Immune.est"]], fe[["OR_Immune.lo"]], fe[["OR_Immune.hi"]], fe[["OR_Immune.k"]],
                 fe[["OR_CAF_ECM.est"]], fe[["OR_CAF_ECM.lo"]], fe[["OR_CAF_ECM.hi"]], fe[["OR_CAF_ECM.k"]],
                 fe[["OR_NetImmune.est"]], fe[["OR_NetImmune.k"]]))
## PRJEB25780 AUC from the stored E7 scores
repro$stored_result[5] <- auc(e7$Responder, e7$NetImmune)
repro$match <- with(repro, ifelse(is.na(published), abs(reproduced - stored_result) < 1e-6,
                                  abs(reproduced - published) < 0.005 + 1e-9 & (is.na(stored_result) | abs(reproduced - stored_result) < 1e-6)))
write.csv(repro, file.path(W, "E19_reproduction.csv"), row.names = FALSE)
cat("\n=== reproduction of published full-score values ===\n"); print(format(repro, digits = 4), row.names = FALSE)
stopifnot(all(repro$match))

## =========================== 1. random drops ===========================
## a drop set is admissible if every module scorable with the full list (>= 3 genes) in a cohort
## keeps >= 3 genes there
SCORABLE <- NFULL >= 3
admissible <- function(drop) {
  keepP <- PRES[, setdiff(GENES, drop), drop = FALSE]
  n <- sapply(names(MOD), function(k) rowSums(keepP[, intersect(MOD[[k]], colnames(keepP)), drop = FALSE]))
  all(n[SCORABLE] >= 3)
}
## post-hoc sensitivity (literal reading of the prereg wording): a drop set is admissible if each
## module keeps >= 3 of its 48-list genes; a module that then has < 3 genes in a cohort is
## unscorable there and that cohort leaves the affected estimates (as in leave-one-gene-out)
admissible_list <- function(drop) all(sapply(MOD, function(g) sum(!g %in% drop)) >= 3)
draw_set <- function(rule, suffix) {
  out <- list(); log <- list()
  for (f in FRACS) {
    nd <- round(f * length(GENES)); acc <- 0; tried <- 0
    while (acc < NDRAW) {
      d <- sample(GENES, nd); tried <- tried + 1
      if (rule(d)) { acc <- acc + 1; out[[length(out) + 1]] <- list(type = sprintf("drop%02d%s", round(100 * f), suffix), frac = f, drop = sort(d)) }
    }
    log[[length(log) + 1]] <- data.frame(rule = if (suffix == "") "cohort coverage (main)" else "48-list only (sensitivity)",
      fraction = f, genes_dropped = nd, accepted = NDRAW, proposed = tried, acceptance = NDRAW / tried)
  }
  list(runs = out, log = do.call(rbind, log))
}
MAIN <- draw_set(admissible, "")              # drawn first, so the main draws do not depend on the sensitivity
SENS <- draw_set(admissible_list, "_list")
RUNS <- c(MAIN$runs, SENS$runs); drawlog <- rbind(MAIN$log, SENS$log)
cat("\n=== random-drop proposals (rejection sampling on module coverage) ===\n"); print(drawlog, row.names = FALSE)
## genes that can never be dropped (a module in some cohort has exactly 3 of its genes, incl. this one)
locked <- GENES[sapply(GENES, function(g) !admissible(g))]
cat("genes never droppable under the coverage rule:", locked, "\n")
for (g in GENES) RUNS[[length(RUNS) + 1]] <- list(type = "logo", frac = 1 / 48, drop = g)

## =========================== 2. panel subsets ===========================
PLAT <- setNames(aud$platform, aud$id)[PANEL]
panel_strict <- GENES[colSums(!PRES[PANEL, , drop = FALSE]) == 0]
plat_sets <- lapply(split(PANEL, PLAT), function(ids) GENES[colSums(!PRES[ids, , drop = FALSE]) == 0])
## (genes on >= 2 of the 3 panel platforms coincide with the NanoString PanCancer Immune content)
plat_count <- Reduce(`+`, lapply(plat_sets, function(s) GENES %in% s))
cat("genes on >= 2 of 3 panel platforms identical to the PanCancer Immune set:",
    setequal(GENES[plat_count >= 2], plat_sets[["NanoString PanCancer Immune (770 genes)"]]), "\n")
PAN <- c(list(panel_strict = panel_strict),
         setNames(plat_sets, paste0("platform_", c("HTG EdgeSeq OBP (~2.5k genes)" = "HTG_OBP",
                  "NanoString IO360 (750 genes)" = "NanoString_IO360",
                  "NanoString PanCancer Immune (770 genes)" = "NanoString_PanCancerImmune")[names(plat_sets)])))
cat("\n=== panel gene sets ===\n")
for (n in names(PAN)) cat(sprintf("%-34s %2d genes; per module %s\n", n, length(PAN[[n]]),
  paste(names(MOD), sapply(MOD, function(g) sum(g %in% PAN[[n]])), sep = "=", collapse = " ")))
for (n in names(PAN)) RUNS[[length(RUNS) + 1]] <- list(type = n, frac = 1 - length(PAN[[n]]) / 48, drop = setdiff(GENES, PAN[[n]]))

## =========================== evaluate all runs (parallel; deterministic) ===========================
cat(sprintf("\nevaluating %d runs on %d cores ...\n", length(RUNS), NCORES)); t0 <- Sys.time()
EV <- mclapply(RUNS, function(r) { e <- eval_keep(setdiff(GENES, r$drop)); e$percohort <- NULL; e }, mc.cores = NCORES)
stopifnot(!any(sapply(EV, inherits, "try-error")))
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

meta_runs <- data.frame(run = seq_along(RUNS), type = sapply(RUNS, `[[`, "type"), fraction = sapply(RUNS, `[[`, "frac"),
                        n_dropped = sapply(RUNS, function(r) length(r$drop)),
                        dropped = sapply(RUNS, function(r) paste(r$drop, collapse = ";")), stringsAsFactors = FALSE)
EST <- cbind(meta_runs, as.data.frame(do.call(rbind, lapply(EV, `[[`, "est"))))
EST <- rbind(cbind(data.frame(run = 0, type = "full", fraction = 0, n_dropped = 0, dropped = ""), t(full$est)), EST)
EST$HR_NetImmune <- exp(EST$HR_NetImmune.est)
RHO <- cbind(EST[, c("run", "type", "fraction")], rbind(full$rho, do.call(rbind, lapply(EV, `[[`, "rho"))))
RHOI <- cbind(EST[, c("run", "type", "fraction")], rbind(full$rho_imm, do.call(rbind, lapply(EV, `[[`, "rho_imm"))))
RHOC <- cbind(EST[, c("run", "type", "fraction")], rbind(full$rho_caf, do.call(rbind, lapply(EV, `[[`, "rho_caf"))))
write.csv(EST, file.path(W, "E19_runs.csv"), row.names = FALSE)
write.csv(RHO, file.path(W, "E19_spearman_runs.csv"), row.names = FALSE)

## =========================== summaries ===========================
## prespecified directions: NetImmune HR < 1; Immune log-OR > 0; CAF_ECM log-OR <= 0.1
dirs <- function(E) cbind(HR = E$HR_NetImmune.est < 0, Imm = E$OR_Immune.est > 0, CAF = E$OR_CAF_ECM.est <= 0.1)
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
summ <- do.call(rbind, lapply(c("full", "drop10", "drop20", "drop30", "drop10_list", "drop20_list", "drop30_list"), function(ty) {
  E <- EST[EST$type == ty, ]; D <- dirs(E)
  do.call(rbind, lapply(list(c("NetImmune HR per SD (benchmark 3, OS)", "HR_NetImmune.est", "HR_NetImmune.p", "exp", "HR < 1"),
                             c("PRJEB25780 NetImmune AUC", "AUC_NetImmune_PRJEB", NA, "id", "(none prespecified)"),
                             c("Immune log-OR per SD (21 ICB)", "OR_Immune.est", "OR_Immune.p", "id", "log-OR > 0"),
                             c("CAF_ECM log-OR per SD (21 ICB)", "OR_CAF_ECM.est", "OR_CAF_ECM.p", "id", "log-OR <= 0.1"),
                             c("NetImmune log-OR per SD (21 ICB)", "OR_NetImmune.est", "OR_NetImmune.p", "id", "(reported)")), function(v) {
    x <- E[[v[2]]]; tr <- if (v[4] == "exp") exp else identity
    dirok <- switch(v[5], "HR < 1" = mean(D[, "HR"]), "log-OR > 0" = mean(D[, "Imm"]), "log-OR <= 0.1" = mean(D[, "CAF"]),
                    "(reported)" = mean(x > 0), NA)
    kk <- if (is.na(v[3])) NA else E[[sub("\\.p$", ".k", v[3])]]
    data.frame(perturbation = ty, n_runs = nrow(E), estimate = v[1], full_score = tr(EST[[v[2]]][EST$type == "full"]),
               median = tr(median(x)), p2.5 = tr(q(x, 0.025)), p5 = tr(q(x, 0.05)), p95 = tr(q(x, 0.95)), p97.5 = tr(q(x, 0.975)),
               min = tr(min(x)), max = tr(max(x)), direction = v[5], frac_direction_kept = dirok,
               frac_P_lt_0.05 = if (is.na(v[3])) NA else mean(E[[v[3]]] < 0.05),
               k_full = if (is.na(v[3])) NA else EST[[sub("\\.p$", ".k", v[3])]][EST$type == "full"],
               k_min = if (is.na(v[3])) NA else min(kk), frac_k_reduced = if (is.na(v[3])) NA else mean(kk < EST[[sub("\\.p$", ".k", v[3])]][EST$type == "full"]))
  }))
}))
joint <- sapply(c("drop10", "drop20", "drop30", "drop10_list", "drop20_list", "drop30_list"), function(ty) mean(apply(dirs(EST[EST$type == ty, ]), 1, all)))
summ$frac_all_three_directions <- ifelse(summ$perturbation %in% names(joint), joint[summ$perturbation], NA)
write.csv(summ, file.path(W, "E19_summary.csv"), row.names = FALSE)
cat("\n=== pooled estimates under random gene dropout ===\n")
print(format(summ[, c("perturbation", "estimate", "full_score", "median", "p2.5", "p97.5", "direction", "frac_direction_kept", "frac_P_lt_0.05")],
             digits = 3), row.names = FALSE)
cat("\nall three prespecified directions kept jointly:\n"); print(round(joint, 4))

## Spearman with the full NetImmune score, per cohort and dropout fraction
coh_names <- names(COH)
sp <- do.call(rbind, lapply(c("drop10", "drop20", "drop30", "drop10_list", "drop20_list", "drop30_list"), function(ty) {
  Rm <- as.matrix(RHO[RHO$type == ty, coh_names])
  rbind(data.frame(perturbation = ty, cohort = coh_names,
                   group = sapply(COH, `[[`, "group"), tech = ifelse(coh_names %in% PANEL, "panel", "whole transcriptome"),
                   median_rho = apply(Rm, 2, median, na.rm = TRUE), p5_rho = apply(Rm, 2, q, 0.05),
                   min_rho = apply(Rm, 2, function(x) if (all(is.na(x))) NA else min(x, na.rm = TRUE)),
                   frac_rho_ge_0.9 = colMeans(Rm >= 0.9, na.rm = TRUE), frac_undefined = colMeans(is.na(Rm))),
        data.frame(perturbation = ty, cohort = "ALL (cohort x draw)", group = "all", tech = "all",
                   median_rho = median(Rm, na.rm = TRUE), p5_rho = q(Rm, 0.05), min_rho = min(Rm, na.rm = TRUE),
                   frac_rho_ge_0.9 = mean(Rm >= 0.9, na.rm = TRUE), frac_undefined = NA),
        data.frame(perturbation = ty, cohort = "median over cohorts of per-cohort median", group = "all", tech = "all",
                   median_rho = median(apply(Rm, 2, median, na.rm = TRUE), na.rm = TRUE), p5_rho = median(apply(Rm, 2, q, 0.05), na.rm = TRUE),
                   min_rho = NA, frac_rho_ge_0.9 = NA, frac_undefined = NA))
}))
rownames(sp) <- NULL
write.csv(sp, file.path(W, "E19_spearman_summary.csv"), row.names = FALSE)
cat("\n=== Spearman rho, perturbed vs full NetImmune (per cohort) ===\n")
w <- reshape(sp[!grepl("_list", sp$perturbation), c("perturbation", "cohort", "median_rho", "p5_rho")], idvar = "cohort", timevar = "perturbation", direction = "wide")
print(format(w, digits = 3), row.names = FALSE)

## gene-level drop frequency among accepted draws (effect of the coverage rule)
gf <- sapply(c("drop10", "drop20", "drop30", "drop10_list", "drop20_list", "drop30_list"), function(ty) {
  dd <- unlist(strsplit(EST$dropped[EST$type == ty], ";")); as.numeric(table(factor(dd, levels = GENES))) / NDRAW })
gf <- data.frame(gene = GENES, module = MOD_OF, gf, locked = GENES %in% locked)
write.csv(gf, file.path(W, "E19_gene_drop_frequency.csv"), row.names = FALSE)
write.csv(drawlog, file.path(W, "E19_draw_acceptance.csv"), row.names = FALSE)

## leave-one-gene-out
LG <- EST[EST$type == "logo", ]; F0 <- EST[EST$type == "full", ]
logo <- data.frame(gene = LG$dropped, module = MOD_OF[LG$dropped],
  HR_NetImmune = LG$HR_NetImmune, dlogHR = LG$HR_NetImmune.est - F0$HR_NetImmune.est, HR_k = LG$HR_NetImmune.k, HR_p = LG$HR_NetImmune.p,
  AUC_PRJEB = LG$AUC_NetImmune_PRJEB, dAUC = LG$AUC_NetImmune_PRJEB - F0$AUC_NetImmune_PRJEB,
  OR_Immune = LG$OR_Immune.est, dImmune = LG$OR_Immune.est - F0$OR_Immune.est, Immune_k = LG$OR_Immune.k, Immune_p = LG$OR_Immune.p,
  OR_CAF_ECM = LG$OR_CAF_ECM.est, dCAF = LG$OR_CAF_ECM.est - F0$OR_CAF_ECM.est, CAF_k = LG$OR_CAF_ECM.k, CAF_p = LG$OR_CAF_ECM.p,
  OR_NetImmune = LG$OR_NetImmune.est, NetImmune_k = LG$OR_NetImmune.k,
  min_rho_NetImmune = apply(as.matrix(RHO[RHO$type == "logo", coh_names]), 1, min, na.rm = TRUE))
## cohorts in which the dropped gene leaves a module unscorable
logo$cohorts_losing_a_module <- sapply(LG$dropped, function(g) {
  k <- MOD_OF[[g]]; lost <- rownames(PRES)[SCORABLE[, k] & (NFULL[, k] - PRES[, g]) < 3]
  if (length(lost)) paste0(k, ": ", paste(lost, collapse = ";")) else "" })
logo$directions_kept <- with(logo, dlogHR + F0$HR_NetImmune.est < 0 & OR_Immune > 0 & OR_CAF_ECM <= 0.1)
write.csv(logo, file.path(W, "E19_leave_one_gene_out.csv"), row.names = FALSE)
cat("\n=== leave-one-gene-out: largest changes per estimate ===\n")
top <- do.call(rbind, lapply(list(c("dlogHR", "NetImmune log HR (benchmark)"), c("dAUC", "PRJEB25780 NetImmune AUC"),
                                  c("dImmune", "Immune log-OR (21 ICB)"), c("dCAF", "CAF_ECM log-OR (21 ICB)")), function(v) {
  o <- order(-abs(logo[[v[1]]]))[1:5]
  data.frame(estimate = v[2], rank = 1:5, gene = logo$gene[o], module = logo$module[o], change = logo[[v[1]]][o],
             cohorts_losing_a_module = logo$cohorts_losing_a_module[o])
}))
write.csv(top, file.path(W, "E19_leave_one_gene_out_top.csv"), row.names = FALSE)
print(format(top, digits = 3), row.names = FALSE)
cat(sprintf("\nLOGO ranges: HR %.3f-%.3f; AUC %.3f-%.3f; Immune log-OR %+.3f to %+.3f; CAF_ECM log-OR %+.3f to %+.3f; all directions kept %d/48\n",
            min(logo$HR_NetImmune), max(logo$HR_NetImmune), min(logo$AUC_PRJEB), max(logo$AUC_PRJEB), min(logo$OR_Immune), max(logo$OR_Immune),
            min(logo$OR_CAF_ECM), max(logo$OR_CAF_ECM), sum(logo$directions_kept)))

## panel subsets
PS <- EST[EST$type %in% names(PAN), ]
pan <- data.frame(subset = PS$type, n_genes = 48 - PS$n_dropped,
  genes = sapply(PS$type, function(n) paste(PAN[[n]], collapse = ";")),
  modules_scorable_in_subset = sapply(PS$type, function(n) paste(names(MOD)[sapply(MOD, function(g) sum(g %in% PAN[[n]]) >= 3)], collapse = " ")),
  HR_NetImmune = PS$HR_NetImmune, HR_lo = exp(PS$HR_NetImmune.lo), HR_hi = exp(PS$HR_NetImmune.hi), HR_p = PS$HR_NetImmune.p,
  AUC_NetImmune_PRJEB = PS$AUC_NetImmune_PRJEB, AUC_Immune_PRJEB = PS$AUC_Immune_PRJEB,
  OR_Immune = PS$OR_Immune.est, OR_Immune_lo = PS$OR_Immune.lo, OR_Immune_hi = PS$OR_Immune.hi, OR_Immune_p = PS$OR_Immune.p, OR_Immune_k = PS$OR_Immune.k,
  OR_CAF_ECM = PS$OR_CAF_ECM.est, OR_CAF_ECM_lo = PS$OR_CAF_ECM.lo, OR_CAF_ECM_hi = PS$OR_CAF_ECM.hi, OR_CAF_ECM_p = PS$OR_CAF_ECM.p, OR_CAF_ECM_k = PS$OR_CAF_ECM.k,
  OR_NetImmune = PS$OR_NetImmune.est, OR_NetImmune_k = PS$OR_NetImmune.k,
  median_rho_NetImmune = apply(as.matrix(RHO[RHO$type %in% names(PAN), coh_names]), 1, median, na.rm = TRUE),
  min_rho_NetImmune = suppressWarnings(apply(as.matrix(RHO[RHO$type %in% names(PAN), coh_names]), 1, min, na.rm = TRUE)),
  median_rho_Immune = apply(as.matrix(RHOI[RHOI$type %in% names(PAN), coh_names]), 1, median, na.rm = TRUE),
  min_rho_Immune = apply(as.matrix(RHOI[RHOI$type %in% names(PAN), coh_names]), 1, min, na.rm = TRUE),
  median_rho_CAF_ECM = apply(as.matrix(RHOC[RHOC$type %in% names(PAN), coh_names]), 1, median, na.rm = TRUE),
  stringsAsFactors = FALSE)
pan$min_rho_NetImmune[!is.finite(pan$min_rho_NetImmune)] <- NA
pan$directions_kept <- with(pan, ifelse(is.na(HR_NetImmune) | is.na(OR_CAF_ECM), NA, HR_NetImmune < 1 & OR_Immune > 0 & OR_CAF_ECM <= 0.1))
write.csv(pan, file.path(W, "E19_panel_subset.csv"), row.names = FALSE)
pr <- cbind(RHO[RHO$type %in% names(PAN), c("type", coh_names)])
write.csv(pr, file.path(W, "E19_panel_spearman.csv"), row.names = FALSE)
cat("\n=== panel subsets ===\n"); print(format(pan[, -3], digits = 3), row.names = FALSE)

## prespecified expectations (main analysis; the 48-list-rule sensitivity alongside)
expect_of <- function(sfx, label) {
  s30 <- sp[sp$perturbation == paste0("drop30", sfx), ]
  fr <- function(est) sapply(paste0(c("drop10", "drop20", "drop30"), sfx), function(ty) summ$frac_direction_kept[summ$estimate == est & summ$perturbation == ty])
  e <- data.frame(analysis = label,
    expectation = c("median Spearman rho with full NetImmune >= 0.9 at 30% dropout (all cohort x draw)",
                    "  ... and in every cohort (lowest per-cohort median)",
                    sprintf("pooled NetImmune HR < 1 in >= 95%% of draws (%s dropout)", c("10%", "20%", "30%")),
                    sprintf("pooled Immune log-OR > 0 in >= 95%% of draws (%s dropout)", c("10%", "20%", "30%")),
                    sprintf("pooled CAF_ECM log-OR <= 0.1 in >= 95%% of draws (%s dropout)", c("10%", "20%", "30%"))),
    value = c(s30$median_rho[s30$cohort == "ALL (cohort x draw)"], min(s30$median_rho[!grepl("ALL|median over", s30$cohort)], na.rm = TRUE),
              fr("NetImmune HR per SD (benchmark 3, OS)"), fr("Immune log-OR per SD (21 ICB)"), fr("CAF_ECM log-OR per SD (21 ICB)")))
  e$threshold <- c(0.9, 0.9, rep(0.95, 9)); e$held <- e$value >= e$threshold; e
}
expc <- rbind(expect_of("", "main (cohort coverage rule)"), expect_of("_list", "sensitivity (48-list rule)"))
write.csv(expc, file.path(W, "E19_expectations.csv"), row.names = FALSE)
cat("\n=== prespecified expectations ===\n"); print(format(expc, digits = 4), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(W, "E19_sessionInfo.txt"))
sink()
