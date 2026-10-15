suppressMessages(library(pROC)); set.seed(42)
V <- "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
d <- read.delim(file.path(V,"benchmark_scores.tsv"), row.names=1, check.names=FALSE)
y <- d$Responder
msi <- d$MSI_H %in% c("True","TRUE",TRUE); ebv <- d$EBV_pos %in% c("True","TRUE",TRUE)
cat("n =",nrow(d)," responders =",sum(y),"\n\n")

M <- list(`MSI/EBV`=as.numeric(msi|ebv), IFN=d$IFN, Immune=d$Immune, GEP18=d$GEP18,
          `PD-L1`=d$PDL1_expr, APM=d$APM, NetImmune=d$NetImmune, T_CELL=d$T_CELL,
          TAM_SUP=d$TAM_SUP, `MSI-H`=as.numeric(msi), TLS_B=d$TLS_B,
          `TIDE(neg)`=-d$TIDE, Barrier=-d$Barrier, IPS=d$IPS)

res <- data.frame(); rocs <- list()
for (nm in names(M)) {
  r <- roc(y, M[[nm]], quiet=TRUE, direction="<"); rocs[[nm]] <- r
  ci <- ci.auc(r, method="bootstrap", boot.n=2000)
  res <- rbind(res, data.frame(Marker=nm, AUC=as.numeric(auc(r)), CI_low=ci[1], CI_high=ci[3]))
}
res <- res[order(-res$AUC),]
cat("=== AUC with bootstrap 95% CI (2000 resamples), n=44, 12 responders ===\n")
print(format(res,digits=3), row.names=FALSE)

cat("\n=== DeLong: NetImmune vs each comparator (paired) ===\n")
dl <- data.frame()
for (nm in setdiff(names(rocs),"NetImmune")) {
  t <- roc.test(rocs[["NetImmune"]], rocs[[nm]], method="delong")
  dl <- rbind(dl, data.frame(Comparator=nm, AUC_comp=as.numeric(auc(rocs[[nm]])),
      Delta=as.numeric(auc(rocs[["NetImmune"]]))-as.numeric(auc(rocs[[nm]])),
      Z=as.numeric(t$statistic), p=t$p.value))
}
dl <- dl[order(dl$p),]; dl$p_FDR <- p.adjust(dl$p,"BH")
print(format(dl,digits=3), row.names=FALSE)

write.csv(res,file.path(O,"A2_auc_ci.csv"),row.names=FALSE)
write.csv(dl, file.path(O,"A2_delong.csv"),row.names=FALSE)
cat("\nsaved.\n")
