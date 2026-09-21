source("R/stock_likelihood_mixture.R")
source("R/stock_mixture_blocks.R")
source("inst/experimental/benchmark-hard-limit.R")
source("inst/experimental/benchmark-streaming-engine.R")
source("inst/experimental/benchmark-stockr-paired.R")
set.seed(991)
g <- matrix(rbinom(45*70,2,.3),45,70)
dimnames(g) <- list(paste0("S",1:45),as.character(1:70))
path <- tempfile(fileext=".gds"); .bp_write_calls(g,path)
stopifnot(isTRUE(all.equal(.hl_read(path),g)))
ctx <- .bs_context(path)
for(k in 1:3) for(strength in c(0,20,50,100,200,500)) {
  p <- matrix(runif(ncol(g)*k,.1,.8),ncol(g),k); pi <- rep(1/k,k)
  fast <- .hl_fit(g,p,pi,strength)
  stopifnot(identical(names(fast$sample.loglik),rownames(g)),identical(ctx$ids,rownames(g)))
  reference <- .sm_block_fit(ctx$reader,ctx$blocks,p,pi,ctx$information,
    300L,1e-8,1e-8,prior.center=colMeans(g)/2,prior.strength=strength)
  for(field in c("frequency","proportion","membership","sample.loglik","objective"))
    stopifnot(isTRUE(all.equal(fast[[field]],reference[[field]],tolerance=1e-7,check.attributes=FALSE)))
}
ctx$close(); unlink(path)
# Guard against silently treating uncertain or missing likelihoods as calls.
uncertain <- list(log.likelihood=array(c(0,-1,-2),c(1,3,1)),informative=matrix(TRUE,1,1))
missing <- list(log.likelihood=array(0,c(1,3,1)),informative=matrix(FALSE,1,1))
stopifnot(inherits(try(.hl_decode(uncertain),silent=TRUE),"try-error"),
  inherits(try(.hl_decode(missing),silent=TRUE),"try-error"))
cat("Hard-limit accelerator matches GDS EM across K and strengths; uncertain input rejected.\n")
