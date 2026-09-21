# One-population technical-missingness negative control, fixed pilot seed.
source("R/stock_likelihood_mixture.R")
source("R/stock_mixture_blocks.R")
set.seed(991)
n <- 60; m <- 60
g <- matrix(rbinom(n*m,2,.35),n,m)
x <- array(-Inf,c(n,3,m))
rates <- rep(c(.05,.3,.65),each=20)
for(l in 1:m) for(i in 1:n) {
  x[i,g[i,l]+1,l] <- 0
  if(runif(1)<rates[i]) x[i,,l] <- 0
}
p <- .sm_prepare(x,"log")
reader <- function(j) list(log.likelihood=p$log.likelihood[,,j,drop=FALSE],
  informative=p$informative[,j,drop=FALSE])
info <- rowSums(p$informative)
scores <- matrix(0,4,2)
for(fold in 1:2) {
  weights <- as.numeric(rep(1:2,30)!=fold)
  pooled <- .sm_block_fit(reader,list(1:m),matrix(.5,m,1),1,info,
    300,1e-8,1e-8,weights)
  stopifnot(pooled$converged)
  for(k in 1:4) {
    fit <- .sm_multistart(reader,list(1:m),matrix(.5,m,k),rep(1/k,k),info,
      7,100+fold,300,1e-8,1e-8,weights,start.strategy="mixed",
      prior.center=pooled$frequency[,1],prior.strength=2)
    scores[k,fold] <- sum(fit$sample.loglik[weights==0])
    stopifnot(fit$converged)
  }
}
print(scores)
print(rowSums(scores))
cat("Best predictive K:",which.max(rowSums(scores)),"\n")
