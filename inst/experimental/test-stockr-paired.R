source("R/stock_likelihood_mixture.R")
source("R/stock_mixture_blocks.R")
source("inst/experimental/benchmark-streaming-engine.R")
source("inst/experimental/benchmark-stockr-paired.R")
source("../grur/R/simulation_benchmark.R")
local({
  t<-.grur_structure_truth(c(20,20,20),80,.15,15981)
  g<-t$genotype; g<-g[,apply(g,2,function(x)length(unique(x))>1),drop=FALSE]
  path<-tempfile(fileext=".gds"); on.exit(unlink(path))
  .bp_write_calls(g,path); context<-.bs_context(path,block.snps=17)
  on.exit(context$close(),add=TRUE,after=FALSE)
  for(index in context$blocks) {
    calls<-apply(context$reader(index)$log.likelihood,c(1,3),which.max)-1L
    stopifnot(all(calls==g[,index,drop=FALSE]))
  }
  for(k in 1:3) {
    starts<-.bp_starts(g,k,851)
    stopifnot(identical(starts,.bp_starts(g,k,851)),all(vapply(starts,function(z) length(unique(z))==k,logical(1))))
    a<-.bp_stock(g,k,starts); b<-.bp_gl(context,g,k,starts,0)
    stopifnot(!is.null(a$fit),!is.null(b$fit))
    pred<-.bp_stock_predict(a$fit,g)
    stopifnot(abs(sum(pred$sample.loglik)-a$fit$margLogL)<1e-5,
      max(abs(rowSums(pred$membership)-1))<1e-10)
    one<-.bp_stock_predict(a$fit,g[1,,drop=FALSE])
    stopifnot(isTRUE(all.equal(unname(one$membership[1,]),unname(pred$membership[1,]),tolerance=1e-10)),
      abs(one$sample.loglik-pred$sample.loglik[1])<1e-10)
    if(k==1) stopifnot(abs(b$fit$loglik-a$fit$margLogL)<1e-6)
  }
  capped<-.bp_stock(g,3,list(starts[[1]]),max.iter=1L)
  stopifnot(is.null(capped$fit),!capped$starts$CONVERGED,capped$starts$ITERATIONS==2L)
  stopifnot(context$audit$max.read.cells<=60*17)
})
cat("Matched blind starts, stockR likelihood projection and K=1 equivalence passed\n")
