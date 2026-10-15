suppressMessages({library(SummarizedExperiment)})
R <- "/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer"
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"

collapse_max <- function(mat, sym) {           # collapse probes -> symbol by max mean expr
  keep <- !is.na(sym) & sym!="" & !duplicated(paste0(sym,"_",seq_along(sym)))
  mat <- mat[keep,,drop=FALSE]; sym <- sym[keep]
  o <- order(sym, -rowMeans(mat,na.rm=TRUE))
  mat <- mat[o,,drop=FALSE]; sym <- sym[o]
  mat <- mat[!duplicated(sym),,drop=FALSE]; rownames(mat) <- sym[!duplicated(sym)]
  mat
}

## ---------------- TCGA-STAD ----------------
se <- readRDS(file.path(R,"data/raw/TCGA_STAD/TCGA_STAD_STAR_counts.rds"))
cat("TCGA assays:",paste(assayNames(se),collapse=", "),"\n")
ass <- if("tpm_unstrand" %in% assayNames(se)) "tpm_unstrand" else assayNames(se)[1]
se <- se[, colData(se)$sample_type=="Primary Tumor"]
m  <- assay(se, ass); sym <- as.character(rowData(se)$gene_name)
m  <- log2(m+1); m <- collapse_max(m, sym)
cd <- as.data.frame(colData(se))
tcga <- data.frame(
  sample = colnames(m),
  patient = substr(colnames(m),1,12),
  time  = ifelse(cd$vital_status=="Dead", cd$days_to_death, cd$days_to_last_follow_up),
  event = as.integer(cd$vital_status=="Dead"),
  stage = as.character(cd$ajcc_pathologic_stage),
  age   = as.numeric(cd$age_at_index), stringsAsFactors=FALSE)
tcga$stage_ord <- as.integer(factor(gsub("[ABC]$","",gsub("^Stage ","",tcga$stage)),
                                    levels=c("I","II","III","IV")))
# one sample per patient, drop missing/non-positive survival
tcga <- tcga[!is.na(tcga$time) & tcga$time>0 & !is.na(tcga$event),]
tcga <- tcga[!duplicated(tcga$patient),]
m <- m[, tcga$sample, drop=FALSE]
cat("TCGA:",nrow(m),"genes x",ncol(m),"tumors | events",sum(tcga$event),"\n")
saveRDS(list(expr=m, clin=tcga), file.path(O,"B1_TCGA.rds"))

## ---------------- GSE84437 ----------------
e <- readRDS(file.path(R,"data_GEO/GSE84437/GSE84437_expr.rds"))
f <- read.csv(file.path(R,"data_GEO/GSE84437/GSE84437_fdata.csv"), stringsAsFactors=FALSE)
symcol <- if("Symbol" %in% colnames(f)) "Symbol" else grep("Symbol|ILMN_Gene",colnames(f),value=TRUE)[1]
cat("GSE84437 symbol col:",symcol,"\n")
if (max(e,na.rm=TRUE) > 50) e <- log2(e+1)
sym <- f[[symcol]][match(rownames(e), f$ID)]
e <- collapse_max(e, sym)
p <- read.csv(file.path(R,"data_GEO/GSE84437/GSE84437_pheno.csv"), stringsAsFactors=FALSE)
g84 <- data.frame(sample=p$geo_accession,
  time  = suppressWarnings(as.numeric(gsub(".*: ","",p$characteristics_ch1.6))),
  event = suppressWarnings(as.integer(gsub(".*: ","",p$characteristics_ch1.5))),
  age   = suppressWarnings(as.numeric(gsub(".*: ","",p$characteristics_ch1.1))),
  t_st  = gsub(".*: ","",p$characteristics_ch1.3),
  n_st  = gsub(".*: ","",p$characteristics_ch1.4), stringsAsFactors=FALSE)
g84$stage_ord <- as.integer(gsub("T","",g84$t_st))   # pT as ordinal stage proxy
g84 <- g84[!is.na(g84$time) & g84$time>0 & !is.na(g84$event),]
g84 <- g84[g84$sample %in% colnames(e),]
e <- e[, g84$sample, drop=FALSE]
cat("GSE84437:",nrow(e),"genes x",ncol(e),"| events",sum(g84$event),"\n")
saveRDS(list(expr=e, clin=g84), file.path(O,"B1_GSE84437.rds"))

## ---------------- GSE62254 (ACRG) ----------------
e2 <- readRDS(file.path(R,"data_GEO/GSE62254/GSE62254_expr.rds"))
cat("GSE62254 raw range:",paste(round(range(e2,na.rm=TRUE),2),collapse=" .. "),"\n")
if (max(e2,na.rm=TRUE) < 6) e2 <- e2/log10(2)      # log10 -> log2
f2 <- read.csv(file.path(R,"data_GEO/GSE62254/GSE62254_fdata.csv"), stringsAsFactors=FALSE)
sym2 <- f2$Gene.Symbol[match(rownames(e2), f2$ID)]
sym2 <- sub(" ///.*","",sym2)
e2 <- collapse_max(e2, sym2)
os <- read.csv(file.path(R,"results/bulk_GEO_fromRDS/GSE62254_ACRG_scores_full_OS.csv"), stringsAsFactors=FALSE)
g62 <- data.frame(sample=os$geo_accession, time=as.numeric(os$OS_time),
                  event=as.integer(os$OS_event), age=suppressWarnings(as.numeric(os$age)),
                  pstage=as.character(os$pStage), stringsAsFactors=FALSE)
g62$stage_ord <- as.integer(factor(gsub("[ABC]$","",g62$pstage), levels=c("I","II","III","IV")))
g62 <- g62[!is.na(g62$time) & g62$time>0 & !is.na(g62$event),]
g62 <- g62[g62$sample %in% colnames(e2),]
e2 <- e2[, g62$sample, drop=FALSE]
cat("GSE62254:",nrow(e2),"genes x",ncol(e2),"| events",sum(g62$event),"\n")
saveRDS(list(expr=e2, clin=g62), file.path(O,"B1_GSE62254.rds"))
cat("\nDONE\n")
