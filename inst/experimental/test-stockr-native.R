# Regression checks for the native-comparison adapter, not stockR modifications.
source("inst/experimental/benchmark-stockr-paired.R")
set.seed(81001)
g <- matrix(rbinom(60*80,2,rep(c(.2,.8),each=30)),60,80)
dimnames(g) <- list(paste0("fish",1:60),paste0("marker",1:80))
for(k in 1:2) {
  set.seed(81002)
  fit <- suppressMessages(stockR::stockSTRUCTURE(t(g),K=k,
    sample.grps=factor(1:60,levels=1:60),start.grps=NULL,control=list(quiet=TRUE)))
  stopifnot(fit$control$method=="Kmeans.EM",fit$control$nKmeanStart==25,
    fit$control$nKPCA==60,fit$control$EM.maxit==100,
    fit$control$EM.minit==5,fit$control$EM.tau.step==.5)
  np <- k-1+nrow(fit$condMeans)*k
  stopifnot(isTRUE(all.equal(fit$AIC,-2*fit$margLogL+2*np)),
    isTRUE(all.equal(fit$BIC,-2*fit$margLogL+log(60)*np)))
  prediction <- .bp_stock_predict(fit,g)
  stopifnot(isTRUE(all.equal(sum(prediction$sample.loglik),fit$margLogL,tolerance=1e-8)))
  # Individual projection cannot silently reorder numerically named sample groups.
  one <- .bp_stock_predict(fit,g[12,,drop=FALSE])
  stopifnot(isTRUE(all.equal(as.numeric(one$membership),as.numeric(prediction$membership[12,]),
    tolerance=1e-8)))
}
cat("Native controls, AIC/BIC and ordered likelihood projection passed.\n")
