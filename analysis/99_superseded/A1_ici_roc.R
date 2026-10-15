suppressMessages({library(pROC)})
set.seed(42)
d <- read.delim("/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation/PRJEB25780_benchmark_scores.tsv",
                row.names=1, check.names=FALSE)
cat("n =", nrow(d), " responders =", sum(d$Responder), "\n\n")

# direction: higher score => more likely responder. TIDE is inverted.
markers <- list(
  NetImmune = d$NetImmune_z,
  Barrier   = -d$Barrier_z,
  GEP18     = d$GEP18,
  `PD-L1`   = d$PDL1_expr,
  `TIDE(neg)`= -d$TIDE,
  IPS       = d$IPS,
  `MSI-H`   = as.numeric(d$MSI_H %in% c("True","TRUE",TRUE,"1")),
  `MSI/EBV` = as.numeric(d$MSI_EBV)
)
y <- d$Responder

res <- data.frame()
rocs <- list()
for (nm in names(markers)) {
  r <- roc(y, markers[[nm]], quiet=TRUE, direction="<")
  rocs[[nm]] <- r
  ci <- ci.auc(r, method="bootstrap", boot.n=2000, progress="")
  res <- rbind(res, data.frame(Marker=nm, AUC=as.numeric(auc(r)),
                               CI_low=ci[1], CI_high=ci[3]))
}
res <- res[order(-res$AUC),]
cat("=== AUC with bootstrap 95% CI (2000 resamples) ===\n")
print(format(res, digits=3), row.names=FALSE)

cat("\n=== DeLong test: NetImmune vs each comparator ===\n")
dl <- data.frame()
for (nm in setdiff(names(rocs),"NetImmune")) {
  t <- roc.test(rocs[["NetImmune"]], rocs[[nm]], method="delong")
  dl <- rbind(dl, data.frame(Comparator=nm,
        AUC_NetImmune=as.numeric(auc(rocs[["NetImmune"]])),
        AUC_comparator=as.numeric(auc(rocs[[nm]])),
        Delta=as.numeric(auc(rocs[["NetImmune"]]))-as.numeric(auc(rocs[[nm]])),
        Z=as.numeric(t$statistic), p=t$p.value))
}
dl$p_FDR <- p.adjust(dl$p, method="BH")
print(format(dl, digits=3), row.names=FALSE)
saveRDS(list(res=res, dl=dl), "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results/A1_roc.rds")
write.csv(res, "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results/A1_auc_ci.csv", row.names=FALSE)
write.csv(dl,  "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results/A1_delong.csv", row.names=FALSE)
