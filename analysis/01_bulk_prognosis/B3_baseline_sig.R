suppressMessages({library(estimate); library(matrixStats)})
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
data(SI_geneset, package="estimate", envir=environment())
si <- SI_geneset
cat("SI_geneset class:",class(si)," dim:",paste(dim(si),collapse=" x "),"\n")
imm <- as.character(unlist(si["ImmuneSignature",-1])); imm <- imm[imm!="" & !is.na(imm)]
str_ <- as.character(unlist(si["StromalSignature",-1])); str_ <- str_[str_!="" & !is.na(str_)]
cat("immune genes:",length(imm)," stromal genes:",length(str_),"\n")

zrow <- function(m){ s<-rowSds(m,na.rm=TRUE); s[s==0|is.na(s)]<-1; (m-rowMeans(m,na.rm=TRUE))/s }
res <- readRDS(file.path(O,"B2_scores.rds"))
for (tag in names(res)) {
  m <- readRDS(file.path(O,paste0("B1_",tag,".rds")))$expr
  Z <- zrow(m); cl <- res[[tag]]
  gi <- intersect(imm, rownames(Z)); gs <- intersect(str_, rownames(Z))
  cl$ESTsig_immune  <- colMeans(Z[gi,,drop=FALSE], na.rm=TRUE)
  cl$ESTsig_stromal <- colMeans(Z[gs,,drop=FALSE], na.rm=TRUE)
  cl$ESTsig_diff    <- cl$ESTsig_immune - cl$ESTsig_stromal
  res[[tag]] <- cl
  cat(sprintf("%-10s imm=%d str=%d | r(ESTsig_imm,CD8A)=%.2f  r(ESTsig_str,Barrier)=%.2f  r(ESTsig_diff,NetImmune)=%.2f\n",
      tag, length(gi), length(gs), cor(cl$ESTsig_immune,cl$CD8A,use="pair"),
      cor(cl$ESTsig_stromal,cl$Barrier,use="pair"), cor(cl$ESTsig_diff,cl$NetImmune,use="pair")))
}
saveRDS(res, file.path(O,"B2_scores.rds"))
cat("\nsaved.\n")
