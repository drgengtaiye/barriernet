## E30 - pooling for Q2 (estimators, randomized part), Q3 (Angio-adjusted estimands a and b), Q4 (endpoint strata)
## Protocol: E30_prereg.md (sha256 in E30_prereg.sha256). Per-comparison inputs: E30_trials.py
## (results/E30_trials_inputs.csv) and the published inputs (E21_meta_inputs.csv, E25_*_inputs.csv, E28 per-comparison
## files). pool() is E21_meta.R's, unchanged (REML point; wider of the REML-HK and DL intervals; larger P); the only
## addition is the REML (Wald) fit for Q2. Step 0 (pooled part): the six-comparison values re-pooled from the E30
## inputs reproduce E21_meta.csv to |difference| < 1e-6 and Table 4 at reported precision; the script stops otherwise.
suppressMessages(library(metafor))
A <- "/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis"; O <- file.path(A, "results")
h <- strsplit(readLines(file.path(A, "E30_prereg.sha256"))[1], " +")[[1]][1]
hh <- system2("sha256sum", file.path(A, "E30_prereg.md"), stdout = TRUE); stopifnot(strsplit(hh, " +")[[1]][1] == h)
cat("E30_prereg.md sha256 matches; metafor", as.character(packageVersion("metafor")), "\n")
options(width = 250)

pool <- function(d) {                                          # E21_meta.R, unchanged
  note <- ""
  re <- withCallingHandlers(rma(yi = d$y, sei = d$se, method = "REML", test = "knha"),
                            warning = function(w) { note <<- gsub("\n", " ", conditionMessage(w)); invokeRestart("muffleWarning") })
  dl <- rma(yi = d$y, sei = d$se, method = "DL")
  dlh <- rma(yi = d$y, sei = d$se, method = "DL", test = "knha")
  fe <- rma(yi = d$y, sei = d$se, method = "FE")
  hk_wider <- (re$ci.ub - re$ci.lb) >= (dl$ci.ub - dl$ci.lb)
  z <- d$y / d$se; zs <- sum(z) / sqrt(nrow(d))
  data.frame(k = nrow(d), k_benefit = sum(d$y > 0), k_benefit_p05 = sum(d$y > 0 & d$p_trial < 0.05),
             k_against_p05 = sum(d$y < 0 & d$p_trial < 0.05),
             est = re$b[1], lo = if (hk_wider) re$ci.lb else dl$ci.lb, hi = if (hk_wider) re$ci.ub else dl$ci.ub,
             p = max(re$pval, dl$pval), interval = if (hk_wider) "HK" else "DL",
             hk_lo = re$ci.lb, hk_hi = re$ci.ub, hk_p = re$pval, tau2_reml = re$tau2, I2 = re$I2,
             dl_est = dl$b[1], dl_lo = dl$ci.lb, dl_hi = dl$ci.ub, dl_p = dl$pval, tau2_dl = dl$tau2,
             dlhk_lo = dlh$ci.lb, dlhk_hi = dlh$ci.ub, dlhk_p = dlh$pval,
             fe_est = fe$b[1], fe_lo = fe$ci.lb, fe_hi = fe$ci.ub, fe_p = fe$pval,
             stouffer_z = zs, stouffer_p = 2 * pnorm(-abs(zs)), reml_warning = note)
}

ORDER <- c("CAF_ECM", "TAM_SUP", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1", "IFN", "APM", "TLS_B", "T_CELL",
           "Barrier", "Immune", "NetImmune", "ESTIMATE", "TMEscore_PC", "TMEscore_z", "CD8A")
Q3S <- c("NetImmune", "Immune", "CAF_ECM", "Barrier", "ESTIMATE", "TMEscore_PC", "TMEscore_z", "CD8A")
LAB <- c(IMM150 = "IMmotion150", ISPY2_pembro = "I-SPY2 pembrolizumab", ISPY2_durva = "I-SPY2 durvalumab-olaparib",
         CM025 = "CheckMate 025", JR101 = "JAVELIN Renal 101", NeoTRIP = "NeoTRIP")
RENAL <- c("IMM150", "CM025", "JR101"); BREAST <- c("ISPY2_pembro", "ISPY2_durva", "NeoTRIP"); ALL6 <- names(LAB)
E21SET <- c(`all six` = "primary: 6 comparisons", renal = "renal only (IMmotion150, CheckMate 025, JAVELIN)",
            breast = "breast only (I-SPY2 x 2, NeoTRIP)", `first three` = "reproduction: 3 published primary",
            `three added` = "3 new primary only")
KS <- list(`all six` = ALL6, renal = RENAL, breast = BREAST, `first three` = ALL6[1:3], `three added` = ALL6[4:6])

IN <- read.csv(file.path(O, "E30_trials_inputs.csv")); IN$p_trial <- IN$p
MI <- read.csv(file.path(O, "E21_meta_inputs.csv")); M21 <- read.csv(file.path(O, "E21_meta.csv"))
U <- IN[IN$analysis == "unadjusted" & IN$key %in% ALL6, ]
pick <- function(d, keys, sc) { d <- d[d$score == sc & d$key %in% keys, ]; d[match(keys, d$key), ] }

## ------------------------------------------------------------------ step 0 (pooled): E21_meta.csv and Table 4
cols <- c("k_benefit", "k_benefit_p05", "est", "lo", "hi", "p", "hk_lo", "hk_hi", "hk_p", "tau2_reml", "dl_est",
          "dl_lo", "dl_hi", "dl_p", "tau2_dl", "dlhk_lo", "dlhk_hi", "dlhk_p", "fe_est", "fe_p", "stouffer_z")
g0 <- do.call(rbind, lapply(names(KS), function(s) do.call(rbind, lapply(ORDER, function(sc) {
  r <- pool(pick(U, KS[[s]], sc)); pb <- M21[M21$set == E21SET[[s]] & M21$score == sc, ]
  data.frame(check = "pooled from E30 inputs vs E21_meta.csv", set = s, score = sc, est = r$est, p = r$p,
             interval_same = r$interval == pb$interval, absdiff = max(abs(unlist(r[cols]) - unlist(pb[cols]))))
}))))
T4 <- data.frame(score = c("CAF_ECM", "TAM_SUP", "Barrier", "IFN", "Immune", "NetImmune", "C1Q_arm", "SPP1_arm", "CXCL9_SPP1",
                           "ESTIMATE", "TMEscore_PC", "TMEscore_z", "CD8A"),
                 est = c(-0.17, 0.05, -0.07, 0.21, 0.20, 0.29, 0.08, -0.03, 0.12, 0.30, 0.30, 0.30, 0.22),
                 lo = c(-0.43, -0.13, -0.28, 0.04, 0.04, 0.03, -0.10, -0.27, -0.04, 0.11, 0.03, 0.05, 0.02),
                 hi = c(0.08, 0.24, 0.15, 0.39, 0.35, 0.55, 0.26, 0.21, 0.28, 0.48, 0.56, 0.56, 0.42),
                 p = c(0.14, 0.48, 0.47, 0.027, 0.021, 0.034, 0.28, 0.73, 0.12, 0.0089, 0.035, 0.029, 0.035),
                 kb = c(0, 4, 2, 6, 5, 6, 4, 3, 5, 6, 5, 6, 4),
                 sp = c(0.0067, 0.82, 0.075, 0.0010, 0.0068, 1.3e-4, 0.47, 0.44, 0.13, 5.5e-5, 9.0e-5, 4.2e-5, 0.018))
t4 <- do.call(rbind, lapply(seq_len(nrow(T4)), function(i) {
  r <- pool(pick(U, ALL6, T4$score[i]))
  data.frame(check = "Table 4 (v31) at reported precision", set = "all six", score = T4$score[i], est = r$est, p = r$p,
             interval_same = NA,
             absdiff = as.numeric(!(abs(r$est - T4$est[i]) <= 0.005 + 1e-9 && abs(r$lo - T4$lo[i]) <= 0.005 + 1e-9 &&
                                    abs(r$hi - T4$hi[i]) <= 0.005 + 1e-9 && signif(r$p, 2) == T4$p[i] &&
                                    r$k_benefit == T4$kb[i] && signif(r$stouffer_p, 2) == T4$sp[i])))
}))
G0 <- rbind(g0, t4); write.csv(G0, file.path(O, "E30_trials_meta_gate.csv"), row.names = FALSE)
cat(sprintf("step 0 (pooled): %d pooled results re-derived from the E30 inputs; largest |difference| from E21_meta.csv %.2e; interval choice identical: %s\n",
            nrow(g0), max(g0$absdiff), all(g0$interval_same)))
cat(sprintf("step 0 (Table 4): %d of %d rows match the published table at reported precision (estimate, CI, P, n favoring, Stouffer P)\n",
            sum(t4$absdiff == 0), nrow(t4)))
stopifnot(max(g0$absdiff) < 1e-6, all(g0$interval_same), all(t4$absdiff == 0))
cat("step 0 pooled gate: PASS\n")

## ------------------------------------------------------------------ Q2: estimators (randomized pools reported, k <= 7)
estimators <- function(d) {
  P <- pool(d); note <- ""
  rz <- withCallingHandlers(rma(yi = d$y, sei = d$se, method = "REML"),
                            warning = function(w) { note <<- gsub("\n", " ", conditionMessage(w)); invokeRestart("muffleWarning") })
  q <- qt(0.975, P$k - 1)
  rbind(data.frame(estimator = "DL", est = P$dl_est, lo = P$dl_lo, hi = P$dl_hi, p = P$dl_p, tau2 = P$tau2_dl, crit = qnorm(0.975)),
        data.frame(estimator = "REML", est = rz$b[1], lo = rz$ci.lb, hi = rz$ci.ub, p = rz$pval, tau2 = rz$tau2, crit = qnorm(0.975)),
        data.frame(estimator = "DL-HK", est = P$dl_est, lo = P$dlhk_lo, hi = P$dlhk_hi, p = P$dlhk_p, tau2 = P$tau2_dl, crit = q),
        data.frame(estimator = "REML-HK", est = P$est, lo = P$hk_lo, hi = P$hk_hi, p = P$hk_p, tau2 = P$tau2_reml, crit = q),
        data.frame(estimator = "manuscript rule", est = P$est, lo = P$lo, hi = P$hi, p = P$p, tau2 = P$tau2_reml,
                   crit = if (P$interval == "HK") q else qnorm(0.975))) |>
    cbind(k = P$k, k_benefit = P$k_benefit, rule_interval = P$interval,
          reml_warning = paste0(P$reml_warning, if (note != "" && note != P$reml_warning) paste0(" [REML Wald] ", note) else ""))
}
MAIN <- list(`reproduction: 3 published primary` = "NetImmune", `3 new primary only` = "NetImmune",
             `primary: 6 comparisons` = T4$score,
             `renal only (IMmotion150, CheckMate 025, JAVELIN)` = c("NetImmune", "CAF_ECM"),
             `breast only (I-SPY2 x 2, NeoTRIP)` = c("NetImmune", "CAF_ECM"),
             `primary, CheckMate 025 OS for PFS` = c("NetImmune", "ESTIMATE", "TMEscore_PC", "TMEscore_z", "CD8A"))
for (s in grep("^primary without", unique(MI$set), value = TRUE)) MAIN[[s]] <- c("NetImmune", "ESTIMATE", "TMEscore_PC", "TMEscore_z", "CD8A")
SETS2 <- setdiff(unique(MI$set), c("+ PRINCE + DUTRENEO (8)", "+ PRINCE + DUTRENEO (>= 1M libraries)"))   # k <= 7
q2 <- do.call(rbind, lapply(SETS2, function(s) do.call(rbind, lapply(ORDER, function(sc) {
  d <- MI[MI$set == s & MI$score == sc, ]; pb <- M21[M21$set == s & M21$score == sc, ]
  e <- estimators(d)
  cbind(family = "score x arm interaction (E21)", set = s, score = sc, e, main_text = sc %in% MAIN[[s]],
        published_rule_est = pb$est, published_rule_lo = pb$lo, published_rule_hi = pb$hi, published_rule_p = pb$p,
        published_source = "E21_meta.csv")
}))))
## E25: within-arm effects and the mutually adjusted interactions; E28: official ESTIMATE, Immune vs CAF_ECM, NetImmune vs ESTIMATE
W <- read.csv(file.path(O, "E25_within_arm_inputs.csv")); W$p_trial <- W$p
J <- read.csv(file.path(O, "E25_joint_inputs.csv")); M25 <- read.csv(file.path(O, "E25_meta.csv"))
CS <- list(`all six` = c("breast", "renal"), breast = "breast", renal = "renal")
fam <- list()
for (sc in c("NetImmune", "CAF_ECM", "Immune")) for (a in c("checkpoint", "control"))
  fam[[length(fam) + 1]] <- list(family = "within-arm effect (E25 Q1)", score = paste(sc, a, "arm"),
                                 d = W[W$score == sc & W$arm == a, c("trial", "cancer", "y", "se", "p_trial")],
                                 pub = M25[M25$question == "Q1 within-arm" & M25$score == sc & M25$arm == a, ],
                                 main = sc %in% c("NetImmune", "CAF_ECM"), src = "E25_meta.csv")
for (v in list(c("Immune, adjusted for CAF_ECM", "immune_adj_y", "immune_adj_se"), c("CAF_ECM, adjusted for Immune", "caf_adj_y", "caf_adj_se")))
  fam[[length(fam) + 1]] <- list(family = "mutually adjusted interaction (E25 Q2)", score = v[1],
                                 d = data.frame(trial = J$trial, cancer = J$cancer, y = J[[v[2]]], se = J[[v[3]]],
                                                p_trial = 2 * pnorm(-abs(J[[v[2]]] / J[[v[3]]]))),
                                 pub = M25[M25$question == "Q2 interaction" & M25$score == v[1], ],
                                 main = TRUE, src = "E25_meta.csv")
cancer_of <- function(tr) ifelse(tr %in% LAB[BREAST], "breast", "renal")
q1 <- read.csv(file.path(O, "E28_Q1_trials_per_comparison.csv")); q1 <- q1[q1$variant == "official, na_min", ]
P28a <- read.csv(file.path(O, "E28_Q1_trials_pooled.csv")); P28b <- read.csv(file.path(O, "E28_Q2_Q6_pooled.csv"))
fam[[length(fam) + 1]] <- list(family = "official ESTIMATE interaction (E28 Q1b)", score = "ESTIMATE official (na_min)",
                               d = data.frame(trial = q1$trial, cancer = cancer_of(q1$trial), y = q1$y, se = q1$se, p_trial = q1$p),
                               pub = P28a[P28a$variant == "na_min", ], main = TRUE,
                               src = "E28_Q1_trials_pooled.csv (numpy rule)")
ja <- read.csv(file.path(O, "E28_Q2a_Q6b_joint.csv")); mb <- read.csv(file.path(O, "E28_Q2b_Q6a_marginal.csv"))
for (v in list(list("Q2a", ja, "Immune minus reversed CAF_ECM, joint model (E28 Q2a)"),
               list("Q6a", mb, "NetImmune minus ESTIMATE, paired marginal (E28 Q6a)"))) {
  dd <- v[[2]][v[[2]]$test == v[[1]], ]
  fam[[length(fam) + 1]] <- list(family = "direct comparison of interactions (E28)", score = v[[3]],
                                 d = data.frame(trial = dd$trial, cancer = cancer_of(dd$trial), y = dd$delta, se = dd$delta_se,
                                                p_trial = dd$delta_p),
                                 pub = P28b[P28b$test == v[[1]], ], main = TRUE,
                                 src = "E28_Q2_Q6_pooled.csv (numpy rule)")
}
q2b <- do.call(rbind, lapply(fam, function(f) do.call(rbind, lapply(names(CS), function(s) {
  d <- f$d[f$d$cancer %in% CS[[s]], ]; pb <- f$pub[f$pub$set == s, ]; stopifnot(nrow(pb) == 1, nrow(d) == if (s == "all six") 6 else 3)
  cbind(family = f$family, set = s, score = f$score, estimators(d), main_text = f$main && s == "all six",
        published_rule_est = pb$est, published_rule_lo = pb$lo, published_rule_hi = pb$hi, published_rule_p = pb$p,
        published_source = f$src)
}))))
Q2 <- rbind(q2, q2b)
rr <- Q2$estimator == "manuscript rule"
Q2$rule_absdiff_vs_published <- ifelse(rr, pmax(abs(Q2$est - Q2$published_rule_est), abs(Q2$lo - Q2$published_rule_lo),
                                                 abs(Q2$hi - Q2$published_rule_hi), abs(Q2$p - Q2$published_rule_p)), NA)
write.csv(Q2, file.path(O, "E30_trials_estimators.csv"), row.names = FALSE)
e21r <- rr & Q2$family == "score x arm interaction (E21)"
cat(sprintf("\nQ2: %d pooled results x 5 estimators; manuscript rule reproduces E21_meta.csv (max |diff| %.1e), E25_meta.csv (%.1e) and the E28 numpy pools (%.1e)\n",
            sum(rr), max(Q2$rule_absdiff_vs_published[e21r]),
            max(Q2$rule_absdiff_vs_published[rr & grepl("E25", Q2$family)]),
            max(Q2$rule_absdiff_vs_published[rr & grepl("E28", Q2$family)])))
stopifnot(max(Q2$rule_absdiff_vs_published[e21r]) < 1e-9, max(Q2$rule_absdiff_vs_published[rr & grepl("E25", Q2$family)]) < 1e-9)
show2 <- Q2[Q2$main_text & Q2$score %in% c("NetImmune", "CAF_ECM", "Immune", "ESTIMATE") &
              Q2$set %in% c("reproduction: 3 published primary", "3 new primary only", "primary: 6 comparisons",
                            "renal only (IMmotion150, CheckMate 025, JAVELIN)", "breast only (I-SPY2 x 2, NeoTRIP)"),
            c("set", "score", "estimator", "k", "est", "lo", "hi", "p", "tau2", "crit")]
show2[, 5:10] <- lapply(show2[, 5:10], signif, 3); print(show2, row.names = FALSE)

## ------------------------------------------------------------------ Q3: Angio-adjusted estimands
ADJ <- IN[IN$analysis == "Angio-adjusted", ]; AXA <- IN[IN$analysis == "Angio x arm (in score model)", ]
ALN <- IN[IN$analysis == "Angio alone", ]
MI6 <- MI[MI$set == "primary: 6 comparisons", ]; MI6$key <- names(LAB)[match(MI6$trial, LAB)]
pooled_rows <- function(d, estimand, set, sc, model) {
  r <- pool(d)
  cbind(level = "pooled", estimand = estimand, set = set, trial = NA, score = sc, model = model, n = NA,
        r[, c("k", "k_benefit", "est", "lo", "hi", "p", "interval", "hk_lo", "hk_hi", "hk_p", "dl_est", "dl_lo", "dl_hi",
              "dl_p", "tau2_reml", "tau2_dl", "I2", "stouffer_z", "stouffer_p", "reml_warning")])
}
q3 <- list()
for (sc in Q3S) {
  q3[[length(q3) + 1]] <- pooled_rows(pick(U, RENAL, sc), "renal, unadjusted (published)", "renal (3)", sc, "unadjusted")
  q3[[length(q3) + 1]] <- pooled_rows(pick(ADJ, RENAL, sc), "(a) renal, Angio-adjusted", "renal (3)", sc, "Angio-adjusted")
  q3[[length(q3) + 1]] <- pooled_rows(pick(MI6, ALL6, sc), "six, unadjusted (published)", "all six", sc, "unadjusted")
  six <- rbind(pick(ADJ, RENAL, sc)[, c("key", "y", "se", "p_trial")], pick(MI6, BREAST, sc)[, c("key", "y", "se", "p_trial")])
  q3[[length(q3) + 1]] <- pooled_rows(six[match(ALL6, six$key), ], "(b) six, renal Angio-adjusted + breast as published",
                                      "all six", sc, "renal Angio-adjusted; breast published")
  q3[[length(q3) + 1]] <- pooled_rows(pick(AXA, RENAL, sc), "Angio x arm, renal (descriptive)", "renal (3)", sc,
                                      "Angio x arm (in score model)")
}
q3[[length(q3) + 1]] <- pooled_rows(pick(ALN, RENAL, "Angio"), "Angio x arm, renal (descriptive)", "renal (3)", "Angio",
                                    "Angio alone (published primary model with Angio as the score)")
Q3P <- do.call(rbind, q3)
tri <- IN[IN$analysis != "unadjusted" | IN$key %in% c(RENAL, "IMM150_mono"), ]
tri <- tri[tri$key %in% c(RENAL, "IMM150_mono") & (tri$score %in% Q3S | tri$analysis == "Angio alone"), ]
Q3T <- data.frame(level = "trial", estimand = ifelse(tri$key == "IMM150_mono", "descriptive: IMmotion150 monotherapy", "per trial"),
                  set = tri$trial, trial = tri$trial, score = tri$score,
                  model = ifelse(tri$analysis == "Angio alone", "Angio alone (published primary model with Angio as the score)", tri$analysis),
                  n = tri$n, k = 1, k_benefit = as.integer(tri$y > 0), est = tri$y, lo = tri$y - qnorm(0.975) * tri$se,
                  hi = tri$y + qnorm(0.975) * tri$se, p = tri$p, interval = "Wald")
Q3T$se <- tri$se
Q3P$se <- NA
unadj <- setNames(tri$y[tri$analysis == "unadjusted"], paste(tri$key, tri$score)[tri$analysis == "unadjusted"])
Q3T$diff_vs_unadjusted <- ifelse(tri$analysis == "Angio-adjusted", tri$y - unadj[paste(tri$key, tri$score)], NA)
Q3T$within_0.05 <- ifelse(tri$analysis == "Angio-adjusted", abs(Q3T$diff_vs_unadjusted) <= 0.05, NA)
Q3P$diff_vs_unadjusted <- NA; Q3P$within_0.05 <- NA
for (sc in Q3S) for (pr in list(c("(a) renal, Angio-adjusted", "renal, unadjusted (published)"),
                                c("(b) six, renal Angio-adjusted + breast as published", "six, unadjusted (published)"))) {
  i <- Q3P$estimand == pr[1] & Q3P$score == sc; j <- Q3P$estimand == pr[2] & Q3P$score == sc
  Q3P$diff_vs_unadjusted[i] <- Q3P$est[i] - Q3P$est[j]
}
allc <- union(names(Q3T), names(Q3P))
fill <- function(d) { for (cc in setdiff(allc, names(d))) d[[cc]] <- NA; d[, allc] }
Q3A <- rbind(fill(Q3T), fill(Q3P))
Q3A <- Q3A[, c("level", "estimand", "set", "trial", "score", "model", "n", "k", "k_benefit", "est", "se", "lo", "hi", "p",
               "interval", "diff_vs_unadjusted", "within_0.05", "hk_lo", "hk_hi", "hk_p", "dl_est", "dl_lo", "dl_hi", "dl_p",
               "tau2_reml", "tau2_dl", "I2", "stouffer_z", "stouffer_p", "reml_warning")]
## pooled unadjusted reference sets must equal the published E21 values
chk <- Q3P[Q3P$estimand %in% c("renal, unadjusted (published)", "six, unadjusted (published)"), ]
ref <- M21[M21$set %in% c(E21SET[["renal"]], E21SET[["all six"]]), ]
ref$estimand <- ifelse(ref$set == E21SET[["renal"]], "renal, unadjusted (published)", "six, unadjusted (published)")
cm <- merge(chk[, c("estimand", "score", "est", "lo", "hi", "p")], ref[, c("estimand", "score", "est", "lo", "hi", "p")], by = c("estimand", "score"))
stopifnot(max(abs(as.matrix(cm[, c("est.x", "lo.x", "hi.x", "p.x")]) - as.matrix(cm[, c("est.y", "lo.y", "hi.y", "p.y")]))) < 1e-9)
write.csv(Q3A, file.path(O, "E30_trials_angio.csv"), row.names = FALSE)
cat("\n=== Q3 pooled (REML; wider of REML-HK / DL interval; larger P) ===\n")
sh <- Q3P[, c("estimand", "score", "k", "k_benefit", "est", "lo", "hi", "p", "interval", "tau2_reml", "diff_vs_unadjusted")]
sh[, c("est", "lo", "hi", "p", "tau2_reml", "diff_vs_unadjusted")] <- lapply(sh[, c("est", "lo", "hi", "p", "tau2_reml", "diff_vs_unadjusted")], signif, 3)
print(sh, row.names = FALSE)
b <- Q3P[Q3P$estimand == "(b) six, renal Angio-adjusted + breast as published" & Q3P$score == "NetImmune", ]
branch <- if (b$est > 0 && b$lo > 0) "1: sign kept, CI excludes zero -> 'the randomized result was not explained by angiogenesis'" else
  if (b$est > 0) "2: sign kept, CI includes zero -> report attenuation; main text: randomized evidence not robust to angiogenesis adjustment" else
    "3: sign reversed -> randomized result reported as explained by angiogenesis"
cat(sprintf("\nQ3 decision (six-comparison NetImmune, estimand b): %+.3f (%.3f to %.3f); P = %.3g -> branch %s\n", b$est, b$lo, b$hi, b$p, branch))

## ------------------------------------------------------------------ Q4: endpoint strata
q4 <- do.call(rbind, lapply(c("renal", "breast", "all six"), function(s) do.call(rbind, lapply(ORDER, function(sc) {
  r <- pool(pick(U, KS[[s]], sc))
  cbind(stratum = s, scale = c(renal = "HR (PFS)", breast = "OR (pCR)", `all six` = "mixed (pre-registered primary)")[[s]],
        score = sc, r[, c("k", "k_benefit", "est", "lo", "hi", "p", "interval", "hk_p", "dl_p", "tau2_reml", "tau2_dl", "I2",
                          "stouffer_z", "stouffer_p", "reml_warning")])
}))))
write.csv(q4, file.path(O, "E30_trials_strata.csv"), row.names = FALSE)
cat("\n=== Q4 strata (randomized rule) ===\n")
w4 <- reshape(q4[, c("stratum", "score", "est", "lo", "hi", "p")], idvar = "score", timevar = "stratum", direction = "wide")
w4[, -1] <- lapply(w4[, -1], signif, 3); print(w4, row.names = FALSE)
cat("\nE30_trials_meta.R done\n")
