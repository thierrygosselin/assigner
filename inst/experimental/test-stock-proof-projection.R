# Hard-call-limit projection only; not valid for uncertain genotype likelihoods.
source("R/stock_mixture_blocks.R")
x <- matrix(-80*log(10),3,3); diag(x) <- 0
for(p in c(1e-8,.01,.5,.99,1-1e-8)) {
  stopifnot(isTRUE(all.equal(.sm_marker(x,p)$loglik,
    dbinom(0:2,2,p,log=TRUE),tolerance=1e-12)))
}
cat("Near-certain GL/binomial projection equivalence passed.\n")
