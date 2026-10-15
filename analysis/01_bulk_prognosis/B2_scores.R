suppressMessages({library(GSVA); library(estimate)})
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
TMP <- "/groups/adv2105_gp/yichen/tmp/claude-3874/-groups-adv2105-gp-yichen-Yi-com-gastric/272d872d-f4a5-4666-a44a-998a271df7f6/scratchpad/est"
dir.create(TMP, showWarnings=FALSE, recursive=TRUE)

MOD <- list(
 CAF_ECM=c("COL1A1","COL1A2","LUM","DCN","POSTN","TAGLN","ACTA2","THY1","FAP","FN1","COL3A1"),
 TAM_SUP=c("SPP1","APOE","TREM2","C1QA","C1QB","C1QC","LGALS3","MARCO","IL10","SOCS3","STAT3"),
 APM=c("B2M","HLA-A","HLA-B","HLA-C","TAP1","TAP2","PSMB8","PSMB9"),
 IFN=c("STAT1","IRF1","CXCL9","CXCL10","IFNG"),
 TLS_B=c("CXCL13","LTB","MS4A1","CD79A","CD74","MZB1"),
 T_CELL=c("CD3D","CD3E","TRAC","TRBC1","CD8A","GZMB","NKG7"))

zrow <- function(m){ s<-matrixStats::rowSds(m,na.rm=TRUE); s[s==0|is.na(s)]<-1
                     (m-rowMeans(m,na.rm=TRUE))/s }

score_cohort <- function(tag){
  d <- readRDS(file.path(O,paste0("B1_",tag,".rds"))); m <- d$expr; cl <- d$clin
  Z <- zrow(m)
  mods <- sapply(MOD, function(g){ g<-intersect(g,rownames(Z)); colMeans(Z[g,,drop=FALSE],na.rm=TRUE) })
  mods <- as.data.frame(mods)
  cl$Barrier   <- rowMeans(mods[,c("CAF_ECM","TAM_SUP")])
  cl$Immune    <- rowMeans(mods[,c("APM","IFN","TLS_B","T_CELL")])
  cl$NetImmune <- cl$Immune - cl$Barrier
  for(k in names(MOD)) cl[[k]] <- mods[[k]]

  ## baseline 1: single gene CD8A
  cl$CD8A <- if("CD8A" %in% rownames(Z)) Z["CD8A",] else NA

  ## baseline 2: ssGSEA on the SAME barriernet genesets
  ss <- gsva(ssgseaParam(as.matrix(m), MOD, minSize=3), verbose=FALSE)
  ss <- t(ss)
  cl$Barrier_ss   <- rowMeans(ss[,c("CAF_ECM","TAM_SUP")])
  cl$Immune_ss    <- rowMeans(ss[,c("APM","IFN","TLS_B","T_CELL")])
  cl$NetImmune_ss <- scale(cl$Immune_ss)[,1] - scale(cl$Barrier_ss)[,1]

  ## baseline 3: ESTIMATE immune - stromal
  fin <- file.path(TMP,paste0(tag,".txt")); fout <- file.path(TMP,paste0(tag,"_common.gct"))
  fsc <- file.path(TMP,paste0(tag,"_est.gct"))
  write.table(data.frame(GeneSymbol=rownames(m), Description=rownames(m), m, check.names=FALSE),
              fin, sep="\t", quote=FALSE, row.names=FALSE)
  ok <- TRUE
  tryCatch({ filterCommonGenes(fin, fout, id="GeneSymbol")
             estimateScore(fout, fsc, platform="affymetrix") },
           error=function(e){ ok<<-FALSE; cat("  ESTIMATE failed:",conditionMessage(e),"\n") })
  if(ok){
    g <- read.delim(fsc, skip=2, row.names=1, check.names=FALSE)[,-1,drop=FALSE]
    colnames(g) <- gsub("^X","",colnames(g)); colnames(g) <- gsub("\\.","-",colnames(g))
    idx <- match(gsub("\\.","-",cl$sample), colnames(g))
    cl$EST_immune  <- as.numeric(g["ImmuneScore", idx])
    cl$EST_stromal <- as.numeric(g["StromalScore", idx])
    cl$EST_diff    <- scale(cl$EST_immune)[,1] - scale(cl$EST_stromal)[,1]
  } else { cl$EST_immune<-NA; cl$EST_stromal<-NA; cl$EST_diff<-NA }

  cat(sprintf("%-10s n=%d events=%d  EST_ok=%s\n", tag, nrow(cl), sum(cl$event), ok))
  cl$cohort <- tag; cl
}

suppressMessages(library(matrixStats))
res <- lapply(c("TCGA","GSE84437","GSE62254"), score_cohort)
names(res) <- c("TCGA","GSE84437","GSE62254")
saveRDS(res, file.path(O,"B2_scores.rds"))
for(n in names(res)) write.csv(res[[n]], file.path(O,paste0("B2_scores_",n,".csv")), row.names=FALSE)
cat("\n=== correlation of NetImmune with baselines ===\n")
for(n in names(res)){ d<-res[[n]]
  cat(sprintf("%-10s r(Net,CD8A)=%.2f  r(Net,ssGSEA)=%.2f  r(Net,ESTIMATEdiff)=%.2f\n", n,
    cor(d$NetImmune,d$CD8A,use="pair"), cor(d$NetImmune,d$NetImmune_ss,use="pair"),
    cor(d$NetImmune,d$EST_diff,use="pair"))) }
