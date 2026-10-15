suppressMessages({library(survival); library(metafor)})
O <- "/groups/adv2105_gp/yichen/Yi/com/gastric/revision/results"
res <- readRDS(file.path(O,"B2_scores.rds")); set.seed(42)

SCORES <- c(NetImmune="NetImmune", Barrier="Barrier", Immune="Immune",
            CD8A="CD8A", `ESTIMATE(imm-str)`="ESTsig_diff", `NetImmune(ssGSEA)`="NetImmune_ss")
CLIN <- list(TCGA=c("age","stage_ord"), GSE84437=c("age","stage_ord","n_ord"), GSE62254=c("age","stage_ord"))

prep <- function(d,tag){
  if(tag=="GSE84437") d$n_ord <- as.integer(gsub("N","",d$n_st))
  for(s in SCORES) d[[paste0(s,"_z")]] <- as.numeric(scale(d[[s]]))
  d
}
res <- Map(prep, res, names(res))

## ---------- 1. univariate Cox, per SD ----------
uni <- do.call(rbind, lapply(names(res), function(tag){
  d <- res[[tag]]
  do.call(rbind, lapply(names(SCORES), function(lab){
    f <- coxph(Surv(time,event) ~ get(paste0(SCORES[lab],"_z")), data=d)
    s <- summary(f)
    data.frame(cohort=tag, score=lab, HR=s$conf.int[1], lo=s$conf.int[3], hi=s$conf.int[4],
               beta=coef(f)[1], se=sqrt(diag(vcov(f)))[1], p=s$coefficients[1,5])
  }))
}))
cat("=== Univariate Cox (HR per 1 SD) ===\n")
print(format(uni[,c("cohort","score","HR","lo","hi","p")],digits=3), row.names=FALSE)

## ---------- 2. meta-analysis with correct p ----------
cat("\n=== Random-effects meta-analysis (DerSimonian-Laird), correct p-values ===\n")
meta <- do.call(rbind, lapply(names(SCORES), function(lab){
  u <- uni[uni$score==lab,]
  m <- rma(yi=u$beta, sei=u$se, method="DL")
  data.frame(score=lab, HR=exp(m$b[1]), lo=exp(m$ci.lb), hi=exp(m$ci.ub),
             p=m$pval, I2=m$I2, Q=m$QE, Qp=m$QEp)
}))
print(format(meta,digits=3), row.names=FALSE)

## ---------- 3. multivariate Cox ----------
cat("\n=== Multivariate Cox (score + clinical) ===\n")
mv <- do.call(rbind, lapply(names(res), function(tag){
  d <- res[[tag]]; cv <- CLIN[[tag]]
  do.call(rbind, lapply(names(SCORES), function(lab){
    f <- try(coxph(as.formula(paste0("Surv(time,event) ~ ",SCORES[lab],"_z + ",paste(cv,collapse=" + "))), data=d), silent=TRUE)
    if(inherits(f,"try-error")) return(NULL)
    s <- summary(f)
    data.frame(cohort=tag, score=lab, HR=s$conf.int[1,1], lo=s$conf.int[1,3],
               hi=s$conf.int[1,4], p=s$coefficients[1,5])
  }))
}))
print(format(mv,digits=3), row.names=FALSE)

## ---------- 4. C-index / delta C-index ----------
cat("\n=== C-index: clinical baseline vs +score (bootstrap 1000) ===\n")
cidx <- do.call(rbind, lapply(names(res), function(tag){
  d <- res[[tag]]; cv <- CLIN[[tag]]
  keep <- complete.cases(d[,c("time","event",cv,paste0(SCORES,"_z"))]); d <- d[keep,]
  f0 <- coxph(as.formula(paste0("Surv(time,event) ~ ",paste(cv,collapse=" + "))), data=d)
  c0 <- concordance(f0)$concordance
  do.call(rbind, lapply(names(SCORES), function(lab){
    f1 <- coxph(as.formula(paste0("Surv(time,event) ~ ",SCORES[lab],"_z + ",paste(cv,collapse=" + "))), data=d)
    c1 <- concordance(f1)$concordance
    db <- replicate(1000, {
      i <- sample(nrow(d), replace=TRUE); dd <- d[i,]
      a <- try(concordance(coxph(as.formula(paste0("Surv(time,event) ~ ",paste(cv,collapse=" + "))),data=dd))$concordance,silent=TRUE)
      b <- try(concordance(coxph(as.formula(paste0("Surv(time,event) ~ ",SCORES[lab],"_z + ",paste(cv,collapse=" + "))),data=dd))$concordance,silent=TRUE)
      if(inherits(a,"try-error")||inherits(b,"try-error")) NA else b-a })
    db <- db[!is.na(db)]
    data.frame(cohort=tag, score=lab, n=nrow(d), C_clin=c0, C_full=c1, dC=c1-c0,
               dC_lo=quantile(db,.025), dC_hi=quantile(db,.975), p_boot=2*min(mean(db<=0),mean(db>=0)))
  }))
}))
print(format(cidx[,c("cohort","score","n","C_clin","C_full","dC","dC_lo","dC_hi","p_boot")],digits=3), row.names=FALSE)

write.csv(uni,  file.path(O,"B4_univariate.csv"), row.names=FALSE)
write.csv(meta, file.path(O,"B4_meta.csv"), row.names=FALSE)
write.csv(mv,   file.path(O,"B4_multivariate.csv"), row.names=FALSE)
write.csv(cidx, file.path(O,"B4_cindex.csv"), row.names=FALSE)
cat("\nsaved.\n")
