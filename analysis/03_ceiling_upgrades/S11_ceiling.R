suppressMessages({library(survival);library(glmnet);library(matrixStats);library(pROC)}); set.seed(42)
P<-"/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"; O<-"results"
zrow<-function(m){s<-rowSds(m,na.rm=TRUE);s[s==0|is.na(s)]<-1;(m-rowMeans(m,na.rm=TRUE))/s}
CT<-c("TCGA","GSE84437","GSE62254"); SC<-readRDS(file.path(P,"B2_scores.rds"))
E<-lapply(CT,function(t) readRDS(file.path(P,paste0("B1_",t,".rds")))$expr); names(E)<-CT
g<-Reduce(intersect,lapply(E,rownames)); cat("genes common to 3 cohorts:",length(g),"\n")
Z<-lapply(E,function(m) t(zrow(m[g,,drop=FALSE])))
keep<-Reduce(`&`,lapply(Z,function(z) apply(is.finite(z),2,all)))
cat("genes finite in all cohorts:",sum(keep),"\n")
Z<-lapply(Z,function(z) z[,keep,drop=FALSE])
Y<-lapply(CT,function(t) Surv(SC[[t]]$time,SC[[t]]$event)); names(Y)<-CT
cidx<-function(p,y){np<--as.numeric(p); as.numeric(concordance(y~np)$concordance)}

cat("\n=== CEILING TEST: full-transcriptome model vs 6-module score (leave-one-cohort-out) ===\n")
out<-data.frame()
for(te in CT){
  tr<-setdiff(CT,te)
  Xtr<-do.call(rbind,Z[tr]); ytr<-Surv(unlist(lapply(tr,function(t)SC[[t]]$time)),
                                       unlist(lapply(tr,function(t)SC[[t]]$event)))
  for(al in c(1,0.5)){
    cv<-cv.glmnet(Xtr,ytr,family="cox",alpha=al,nfolds=5,standardize=FALSE)
    p<-as.numeric(predict(cv,Z[[te]],s="lambda.min"))
    nz<-sum(coef(cv,s="lambda.min")!=0)
    out<-rbind(out,data.frame(test=te,model=sprintf("full transcriptome EN(a=%.1f)",al),
                              C=cidx(p,Y[[te]]),n_genes=nz))
  }
  out<-rbind(out,data.frame(test=te,model="NetImmune (6 modules, fixed)",
                            C=cidx(-scale(SC[[te]]$NetImmune)[,1],Y[[te]]),n_genes=48))
  out<-rbind(out,data.frame(test=te,model="ESTIMATE (282 genes, fixed)",
                            C=cidx(-scale(SC[[te]]$ESTsig_diff)[,1],Y[[te]]),n_genes=282))
}
W<-reshape(out[,c("test","model","C")],idvar="model",timevar="test",direction="wide")
colnames(W)<-gsub("^C\\.","",colnames(W)); W$mean<-rowMeans(W[,CT]); W<-W[order(-W$mean),]
print(format(W,digits=3),row.names=FALSE)
write.csv(out,file.path(O,"S11_ceiling_prognosis.csv"),row.names=FALSE)
