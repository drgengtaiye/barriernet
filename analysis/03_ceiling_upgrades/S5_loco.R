suppressMessages({library(survival); library(glmnet)}); set.seed(42)
P<-"/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
O<-"/groups/adv2105_gp/yichen/Yi/com/gastric/revision2/results"
FE<-readRDS(file.path(O,"S4_features.rds")); SC<-readRDS(file.path(P,"B2_scores.rds"))
CT<-c("TCGA","GSE84437","GSE62254")
MODS<-c("CAF_ECM","TAM_SUP","APM","IFN","TLS_B","T_CELL")
CLIN<-list(TCGA=c("age","stage_ord"),GSE84437=c("age","stage_ord","n_ord"),GSE62254=c("age","stage_ord"))

D<-list()
for(tg in CT){
  s<-SC[[tg]]; if(tg=="GSE84437") s$n_ord<-as.integer(gsub("N","",s$n_st))
  x<-FE[[tg]]$X; stopifnot(nrow(x)==nrow(s))
  D[[tg]]<-list(new=as.matrix(x), mod=as.matrix(scale(s[,MODS])),
                fixed=data.frame(NetImmune=scale(s$NetImmune)[,1], ESTIMATE=scale(s$ESTsig_diff)[,1],
                                 CD8A=scale(s$CD8A)[,1]),
                y=Surv(s$time,s$event), clin=s[,CLIN[[tg]],drop=FALSE], s=s)
}
fit_en<-function(Xtr,ytr,Xte){
  cv<-cv.glmnet(Xtr,ytr,family="cox",alpha=.5,nfolds=10,standardize=FALSE)
  list(pred=as.numeric(predict(cv,Xte,s="lambda.min")),
       beta=as.matrix(coef(cv,s="lambda.min")))
}
# convention: pred is a RISK score (higher = worse). concordance() expects a protective score.
cidx<-function(pred,y){ np <- -as.numeric(pred); as.numeric(concordance(y ~ np)$concordance) }

SETS<-list(
  `NetImmune (fixed)`  =function(d) -d$fixed$NetImmune,
  `ESTIMATE (fixed)`   =function(d) -d$fixed$ESTIMATE,
  `CD8A (fixed)`       =function(d) -d$fixed$CD8A)
LEARN<-list(
  `EN: 6 barriernet modules`=function(d) d$mod,
  `EN: abundance only (10)` =function(d) d$new[,grepl("^ab_",colnames(d$new)),drop=FALSE],
  `EN: state only (20)`     =function(d) d$new[,grepl("^st_",colnames(d$new)),drop=FALSE],
  `EN: abundance + state (30)`=function(d) d$new)

rows<-list(); betas<-list()
for(te in CT){
  tr<-setdiff(CT,te)
  ytr<-do.call(c,lapply(tr,function(t) D[[t]]$y)); ytr<-Surv(unlist(lapply(tr,function(t) D[[t]]$s$time)),
                                                            unlist(lapply(tr,function(t) D[[t]]$s$event)))
  for(nm in names(SETS)){
    p<-SETS[[nm]](D[[te]])
    rows[[length(rows)+1]]<-data.frame(test=te,model=nm,C=cidx(p,D[[te]]$y),type="fixed")
  }
  for(nm in names(LEARN)){
    Xtr<-do.call(rbind,lapply(tr,function(t) LEARN[[nm]](D[[t]])))
    Xte<-LEARN[[nm]](D[[te]])
    f<-fit_en(Xtr,ytr,Xte)
    rows[[length(rows)+1]]<-data.frame(test=te,model=nm,C=cidx(f$pred,D[[te]]$y),type="learned")
    if(nm=="EN: abundance + state (30)") betas[[te]]<-f$beta
  }
}
R<-do.call(rbind,rows)
W<-reshape(R[,c("test","model","C")],idvar="model",timevar="test",direction="wide")
colnames(W)<-gsub("^C\\.","",colnames(W)); W$mean<-rowMeans(W[,CT])
W<-W[order(-W$mean),]
cat("=== Leave-one-cohort-out C-index (trained on 2 cohorts, tested on the third) ===\n")
print(format(W,digits=3),row.names=FALSE)
write.csv(W,file.path(O,"S5_loco_cindex.csv"),row.names=FALSE)
B<-do.call(cbind,betas); colnames(B)<-names(betas)
B<-B[rowSums(abs(B))>0,,drop=FALSE]
cat("\n=== elastic-net coefficients (abundance+state model), by held-out cohort ===\n")
print(round(B,3))
write.csv(as.matrix(B),file.path(O,"S5_coefs.csv"))
