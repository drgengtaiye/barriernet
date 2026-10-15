## E9 - re-derive the manuscript statement that the official ssGSEA ESTIMATE immune score
## runs against canonical immune genes on the GSE84437 array (v13: r = -0.40 with PTPRC,
## CD2, CD3E and CD8A). Official scores are EST_immune in B2_scores.rds (B2b_estimate_fix.R).
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
S <- readRDS(file.path(R, "B2_scores.rds"))
for (t in names(S)) {
  e <- readRDS(file.path(R, paste0("B1_", t, ".rds")))$expr
  g <- intersect(c("PTPRC", "CD2", "CD3E", "CD8A"), rownames(e))
  z <- t(scale(t(e[g, , drop = FALSE])))
  m <- colMeans(z, na.rm = TRUE)[S[[t]]$sample]
  cat(sprintf("%-9s genes %d  r(official ESTIMATE immune, mean z of %s) = %+.2f  | r(signature-level ESTIMATE immune) = %+.2f\n",
              t, length(g), paste(g, collapse = "/"), cor(S[[t]]$EST_immune, m, use = "pair"),
              cor(S[[t]]$ESTsig_immune, m, use = "pair")))
}
