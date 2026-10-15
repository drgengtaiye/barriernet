suppressMessages({library(matrixStats)})
P <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
LD <- read.csv(file.path(O,"S3_state_loadings.csv"), row.names=1, check.names=FALSE)
ID <- read.csv(file.path(O,"S3_identity_genes.csv"), row.names=1)
ident <- setNames(strsplit(as.character(ID[,1]), ","), rownames(ID))
zrow <- function(m){ s<-rowSds(m,na.rm=TRUE); s[s==0|is.na(s)]<-1; (m-rowMeans(m,na.rm=TRUE))/s }

build <- function(tag){
  d <- readRDS(file.path(P,paste0("B1_",tag,".rds"))); Z <- zrow(d$expr); cl <- d$clin
  F <- data.frame(row.names=colnames(Z))
  # --- abundance: mean z of identity genes ---
  for (c in names(ident)) {
    g <- intersect(ident[[c]], rownames(Z))
    if (length(g)>=8) F[[paste0("ab_",c)]] <- colMeans(Z[g,,drop=FALSE], na.rm=TRUE)
  }
  # --- state: loading-weighted projection, then residualise on own abundance ---
  for (ax in colnames(LD)) {
    ct <- sub("_S[0-9]+$","",ax); ab <- paste0("ab_",ct)
    L <- LD[[ax]]; names(L) <- rownames(LD); L <- L[!is.na(L)]
    g <- intersect(names(L), rownames(Z))
    if (length(g)<15 || !(ab %in% colnames(F))) next
    raw <- as.numeric(crossprod(Z[g,,drop=FALSE], L[g])) / sum(abs(L[g]))
    F[[paste0("st_",ax)]] <- residuals(lm(raw ~ F[[ab]]))     # orthogonal to abundance
  }
  F <- as.data.frame(scale(F))
  cat(sprintf("%-10s n=%d  features=%d (abundance %d, state %d)\n", tag, nrow(F), ncol(F),
      sum(grepl("^ab_",colnames(F))), sum(grepl("^st_",colnames(F)))))
  list(X=F, clin=cl)
}
res <- lapply(c("TCGA","GSE84437","GSE62254"), build); names(res) <- c("TCGA","GSE84437","GSE62254")
common <- Reduce(intersect, lapply(res, function(r) colnames(r$X)))
cat("\nfeatures common to all 3 cohorts:", length(common), "\n")
for(n in names(res)) res[[n]]$X <- res[[n]]$X[,common,drop=FALSE]
saveRDS(res, file.path(O,"S4_features.rds"))
cat(paste(common, collapse=", "), "\n")
