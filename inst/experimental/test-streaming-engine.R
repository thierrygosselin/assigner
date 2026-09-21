source("R/assign_likelihood_dapc.R")
for(p in c("stock_likelihood_mixture","stock_mixture_blocks","stock_mixture_gds"))
  source(file.path("R",paste0(p,".R")))
source("inst/experimental/benchmark-vectorized-engine.R")
source("inst/experimental/benchmark-streaming-engine.R")
source("../grur/R/simulation_benchmark.R")
set.seed(721)
for(n in c(1L,9L)) for(m in c(1L,13L)) for(e in c("GL","PL","log")) {
  x<-array(if(e=="PL") sample(0:100,n*3*m,TRUE) else -rexp(n*3*m),c(n,3,m))
  x[1,,1]<-NA
  if(n>1) x[2,,1]<-0
  if(n>2&&e!="PL") x[3,1,1]<--Inf
  stopifnot(isTRUE(all.equal(.bs_prepare(x,e),.sm_prepare(x,e),tolerance=1e-12)))
}
for(v in list(c(NA,0,1),c(-Inf,-Inf,-Inf),c(Inf,0,0),c(NaN,0,0)))
  stopifnot(inherits(try(.bs_prepare(array(v,c(1,3,1))),silent=TRUE),"try-error"))
local({
  t<-.grur_structure_truth(c(8,8,8),13,.05,851)
  o<-.grur_observe_structure(t,missing=.1,depth=3,seed=812)
  path<-tempfile(fileext=".gds"); on.exit(unlink(path))
  .grur_write_benchmark_gds(o,t$metadata$INDIVIDUALS,path)
  hash<-tools::md5sum(path); cached<-.bv_context(path)
  for(size in c(1L,4L,13L)) {
    stream<-.bs_context(path,block.snps=size,max.cells=24*size)
    tryCatch({
      for(index in stream$blocks) stopifnot(isTRUE(all.equal(
        stream$reader(index),cached$reader(index),tolerance=1e-12)))
      for(s in c(0,100)) {
        args<-list(frequency=t$frequency,proportion=rep(1/3,3),information=stream$information,
          max.iter=300,tolerance=1e-8,epsilon=1e-8,sample.weights=rep(c(1,1,0),8),
          prior.center=rep(.4,13),prior.strength=s)
        a<-do.call(.bv_fit,c(list(read.block=stream$reader,blocks=stream$blocks),args))
        b<-do.call(.bv_fit,c(list(read.block=cached$reader,blocks=cached$blocks),args))
        for(key in c("frequency","membership","sample.loglik","objective"))
          stopifnot(isTRUE(all.equal(a[[key]],b[[key]],tolerance=1e-7)))
      }
      stopifnot(stream$audit$max.read.cells<=24*size)
      stopifnot(inherits(try(stream$reader(c(2,1)),silent=TRUE),"try-error"))
    },finally=stream$close())
    stream$close() # Closing twice is safe; reading after close must fail.
    stopifnot(inherits(try(stream$reader(1),silent=TRUE),"try-error"))
  }
  sub<-.bs_context(path,sample.ids=t$metadata$INDIVIDUALS[c(2,4,7)],variant.ids=c(2L,5L))
  stopifnot(identical(sub$ids,t$metadata$INDIVIDUALS[c(2,4,7)]),sub$m==2L)
  sub$close(); stopifnot(identical(hash,tools::md5sum(path)))
})
cat("Streaming GL preparation, block-size equivalence, weighted fits, filters, bounds and closure checks passed\n")
