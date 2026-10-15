suppressMessages(library(estimate))
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
TMP <- "/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-com-gastric/272d872d-f4a5-4666-a44a-998a271df7f6/scratchpad/est2"
dir.create(TMP, showWarnings=FALSE, recursive=TRUE)
res <- readRDS(file.path(O,"B2_scores.rds"))
for (tag in names(res)) {
  m <- readRDS(file.path(O,paste0("B1_",tag,".rds")))$expr
  fin <- file.path(TMP,paste0(tag,".txt"))
  # standard ESTIMATE input: gene symbol in col 1, samples thereafter. NO extra text column.
  out <- cbind(GeneSymbol=rownames(m), as.data.frame(m, check.names=FALSE))
  write.table(out, fin, sep="\t", quote=FALSE, row.names=FALSE)
  fc <- file.path(TMP,paste0(tag,"_c.gct")); fs <- file.path(TMP,paste0(tag,"_e.gct"))
  filterCommonGenes(fin, fc, id="GeneSymbol")
  estimateScore(fc, fs, platform="affymetrix")
  g <- read.delim(fs, skip=2, row.names=1, check.names=FALSE)
  g <- g[, setdiff(colnames(g), c("Description","Description.1")), drop=FALSE]
  cl <- res[[tag]]
  idx <- match(make.names(cl$sample), make.names(colnames(g)))
  cl$EST_immune  <- as.numeric(g["ImmuneScore", idx])
  cl$EST_stromal <- as.numeric(g["StromalScore", idx])
  cl$EST_diff    <- scale(cl$EST_immune)[,1] - scale(cl$EST_stromal)[,1]
  res[[tag]] <- cl
  cat(sprintf("%-10s unmatched=%d  r(ESTimm,Immune)=%.2f  r(ESTstr,Barrier)=%.2f  r(ESTimm,CD8A)=%.2f\n",
      tag, sum(is.na(idx)), cor(cl$EST_immune,cl$Immune,use="pair"),
      cor(cl$EST_stromal,cl$Barrier,use="pair"), cor(cl$EST_immune,cl$CD8A,use="pair")))
}
saveRDS(res, file.path(O,"B2_scores.rds"))
for(n in names(res)) write.csv(res[[n]], file.path(O,paste0("B2_scores_",n,".csv")), row.names=FALSE)
cat("\nsaved.\n")
