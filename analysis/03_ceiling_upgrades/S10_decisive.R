suppressMessages({library(survival);library(matrixStats);library(pROC)}); set.seed(42)
P<-"/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"; O<-"results"
V<-"/groups/adv2105_gp/yichen/Yi/com/gastric/GC_TME_Reprogrammer/validation"
W<-read.csv(file.path(O,"S9_arrangement_weights.csv"),row.names=1); w<-setNames(W[,1],rownames(W))
zrow<-function(m){s<-rowSds(m,na.rm=TRUE);s[s==0|is.na(s)]<-1;(m-rowMeans(m,na.rm=TRUE))/s}
sas<-function(Z){ g<-intersect(names(w),rownames(Z))
  cat("   genes used:",length(g),"/",length(w),"\n")
  as.numeric(crossprod(Z[g,,drop=FALSE],w[g]))/sum(abs(w[g])) }

SC<-readRDS(file.path(P,"B2_scores.rds")); CT<-c("TCGA","GSE84437","GSE62254")
CLIN<-list(TCGA=c("age","stage_ord"),GSE84437=c("age","stage_ord","n_ord"),GSE62254=c("age","stage_ord"))
cat("=== PROGNOSIS: does the spatial-arrangement score (SAS) add to NetImmune? ===\n")
tab<-data.frame()
for(tg in CT){
  d<-readRDS(file.path(P,paste0("B1_",tg,".rds"))); s<-SC[[tg]]
  if(tg=="GSE84437") s$n_ord<-as.integer(gsub("N","",s$n_st))
  cat(" ",tg,"\n"); s$SAS<-scale(sas(zrow(d$expr)))[,1]
  s$NetImmune_z<-scale(s$NetImmune)[,1]; s$Immune_z<-scale(s$Immune)[,1]
  y<-Surv(s$time,s$event); cv<-CLIN[[tg]]
  u<-summary(coxph(y~SAS,data=s))
  m<-summary(coxph(as.formula(paste("y~SAS+NetImmune_z+",paste(cv,collapse="+"))),data=s))
  c0<-concordance(coxph(as.formula(paste("y~NetImmune_z+",paste(cv,collapse="+"))),data=s))$concordance
  c1<-concordance(coxph(as.formula(paste("y~SAS+NetImmune_z+",paste(cv,collapse="+"))),data=s))$concordance
  tab<-rbind(tab,data.frame(cohort=tg,SAS_uni_HR=u$conf.int[1],SAS_uni_p=u$coefficients[1,5],
    SAS_adj_HR=m$conf.int[1,1],SAS_adj_p=m$coefficients[1,5],
    r_SAS_Net=cor(s$SAS,s$NetImmune_z), C_without=c0,C_with=c1,dC=c1-c0))
  SC[[tg]]<-s
}
print(format(tab,digits=3),row.names=FALSE)
write.csv(tab,file.path(O,"S10_prognosis.csv"),row.names=FALSE)

cat("\n=== ICI: does SAS add to Immune / NetImmune? ===\n")
e<-as.matrix(read.delim(file.path(V,"mGC_Pembrolizumab_Kim2018.TPM"),row.names=1,check.names=FALSE))
e<-log2(e+1); e<-e[rowSums(is.finite(e))==ncol(e),]
b<-read.delim(file.path(V,"benchmark_scores.tsv"),row.names=1,check.names=FALSE)
S<-scale(sas(zrow(e)))[,1]; names(S)<-colnames(e); b$SAS<-S[rownames(b)]
y<-b$Responder
for(v in c("SAS","NetImmune","Immune")){
  r<-roc(y,b[[v]],quiet=TRUE,direction="<"); a<-as.numeric(auc(r)); fl<-a<.5
  if(fl){r<-roc(y,-b[[v]],quiet=TRUE,direction="<");a<-as.numeric(auc(r))}
  ci<-ci.auc(r,method="bootstrap",boot.n=2000)
  cat(sprintf("  %-10s AUC=%.3f (%.3f-%.3f) %s\n",v,a,ci[1],ci[3],ifelse(fl,"[inverted]","")))
}
cat(sprintf("  r(SAS, Immune) = %.3f | r(SAS, NetImmune) = %.3f\n",cor(b$SAS,b$Immune),cor(b$SAS,b$NetImmune)))
m0<-glm(y~scale(Immune),data=b,family=binomial); m1<-glm(y~scale(Immune)+scale(SAS),data=b,family=binomial)
cat(sprintf("  LRT SAS added to Immune: chi2=%.2f p=%.3f\n",2*(logLik(m1)-logLik(m0)),
    pchisq(2*as.numeric(logLik(m1)-logLik(m0)),1,lower.tail=FALSE)))
write.csv(b,file.path(O,"S10_ici_with_SAS.csv"))
