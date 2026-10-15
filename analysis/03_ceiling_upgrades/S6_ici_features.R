suppressMessages(library(matrixStats))
O<-"/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
V<-"/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
LD<-read.csv(file.path(O,"S3_state_loadings.csv"),row.names=1,check.names=FALSE)
ID<-read.csv(file.path(O,"S3_identity_genes.csv"),row.names=1)
ident<-setNames(strsplit(as.character(ID[,1]),","),rownames(ID))
zrow<-function(m){s<-rowSds(m,na.rm=TRUE);s[s==0|is.na(s)]<-1;(m-rowMeans(m,na.rm=TRUE))/s}

e<-as.matrix(read.delim(file.path(V,"mGC_Pembrolizumab_Kim2018.TPM"),row.names=1,check.names=FALSE))
e<-log2(e+1); e<-e[rowSums(is.finite(e))==ncol(e),]
cat("ICI expr:",dim(e),"\n")
Z<-zrow(e)
F<-data.frame(row.names=colnames(Z))
for(c in names(ident)){ g<-intersect(ident[[c]],rownames(Z))
  if(length(g)>=8) F[[paste0("ab_",c)]]<-colMeans(Z[g,,drop=FALSE],na.rm=TRUE) }
for(ax in colnames(LD)){ ct<-sub("_S[0-9]+$","",ax); ab<-paste0("ab_",ct)
  L<-LD[[ax]];names(L)<-rownames(LD);L<-L[!is.na(L)]
  g<-intersect(names(L),rownames(Z)); if(length(g)<15||!(ab%in%colnames(F))) next
  raw<-as.numeric(crossprod(Z[g,,drop=FALSE],L[g]))/sum(abs(L[g]))
  F[[paste0("st_",ax)]]<-residuals(lm(raw~F[[ab]])) }
F<-as.data.frame(scale(F))
b<-read.delim(file.path(V,"benchmark_scores.tsv"),row.names=1,check.names=FALSE)
F<-F[rownames(b),,drop=FALSE]
out<-cbind(b[,c("Responder","MSI_H","EBV_pos","NetImmune","Barrier","Immune","GEP18","PDL1_expr","TIDE","IPS")],F)
write.csv(out,file.path(O,"S6_ici_features.csv"))
cat("features:",ncol(F)," samples:",nrow(F),"\n")
