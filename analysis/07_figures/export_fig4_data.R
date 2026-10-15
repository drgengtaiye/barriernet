## Export the tables make_fig4.py reads (figures/data/): per-cohort survival with NetImmune and
## Barrier (from revision/results/B2_scores.rds) and the structure-matched null draws (from
## 2/analysis/results/E2_null_draws.rds). These were first written by an inline command;
## this script reproduces them. Usage: Rscript export_fig4_data.R [outdir]
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "/groups/adv2105_gp/yichen/Yi/com/gastric/2/figures/data"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
r <- readRDS("/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results/B2_scores.rds")
for (t in names(r)) write.csv(r[[t]][, c("sample", "time", "event", "NetImmune", "Barrier")],
                              file.path(out, sprintf("km_%s.csv", t)), row.names = FALSE)
x <- readRDS("/groups/adv2105_gp/yichen/Yi/com/gastric/2/analysis/results/E2_null_draws.rds")
w <- data.frame(set = seq_len(dim(x$rnd)[1]), lhr_TCGA = x$rnd[, 1, "b"], lhr_GSE84437 = x$rnd[, 2, "b"],
                lhr_GSE62254 = x$rnd[, 3, "b"], lhr_pooled = x$rnd_pool)
write.csv(w, file.path(out, "null_draws.csv"), row.names = FALSE)
write.csv(data.frame(cohort = colnames(x$obs), lhr = x$obs["b", ]), file.path(out, "null_obs.csv"), row.names = FALSE)
