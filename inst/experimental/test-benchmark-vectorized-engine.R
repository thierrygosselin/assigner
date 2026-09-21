source("R/assign_likelihood_dapc.R")
for(p in c("stock_likelihood_mixture","stock_mixture_blocks","stock_mixture_gds","stock_mixture_validation"))
  source(file.path("R",paste0(p,".R")))
source("inst/experimental/benchmark-vectorized-engine.R")
source("../grur/R/simulation_benchmark.R")
set.seed(193)
for(n in c(1,8)) for(m in c(1,7)) {
  x<-array(-rexp(n*3*m),c(n,3,m)); x[1,,1]<-0
  if(n>1) x[2,1,1]<--Inf
  b<-.sm_prepare(x,"log"); p<-seq(.001,.999,length.out=m)
  a<-.bv_moments(b,p)
  for(l in seq_len(m)) {
    z<-.sm_marker(matrix(b$log.likelihood[,,l],n,3),p[l])
    stopifnot(isTRUE(all.equal(a$loglik[,l],z$loglik,tolerance=1e-12)),
      isTRUE(all.equal(a$dosage[,l],z$dosage,tolerance=1e-12)))
  }
}
truth<-.grur_structure_truth(c(8,8,8),12,.1,31)
for(depth in c(Inf,2)) {
  obs<-.grur_observe_structure(truth,missing=.1,depth=depth,seed=21)
  path<-tempfile(fileext=".gds")
  .grur_write_benchmark_gds(obs,truth$metadata$INDIVIDUALS,path)
  context<-.bv_context(path)
  weights<-rep(c(1,1,0),8)
  for(strength in c(0,2,20,100)) {
    old<-.sm_block_fit(context$reader,context$blocks,truth$frequency,rep(1/3,3),
      context$information,300,1e-8,1e-8,weights,rep(.4,12),strength)
    fast<-.bv_fit(context$reader,context$blocks,truth$frequency,rep(1/3,3),
      context$information,300,1e-8,1e-8,weights,rep(.4,12),strength)
    for(key in c("frequency","membership","proportion","sample.loglik","objective","history"))
      stopifnot(isTRUE(all.equal(old[[key]],fast[[key]],tolerance=1e-7)))
    stopifnot(identical(old$converged,fast$converged))
  }
  initial<-matrix(.5,12,3,dimnames=list(context$variants,paste0("STOCK_",1:3)))
  for(strength in c(20,100)) {
  reference<-.sm_fit_gds(path,initial,likelihood="GL",n.starts=7,seed=91,
    sample.weights=setNames(weights,context$ids),start.strategy="mixed",
    regularization="pooled",prior.strength=strength,max.iter=300,block.snps=250)
  original<-.sm_block_fit; .sm_block_fit<-.bv_fit
  accelerated<-.bv_multifit(context,3,strength,weights,seed=91)
  .sm_block_fit<-original
  for(key in c("frequency","membership","sample.loglik","objective"))
    stopifnot(isTRUE(all.equal(unname(reference[[key]]),unname(accelerated[[key]]),tolerance=1e-6)))
  }
  unlink(path)
}
cat("Vectorized moments, weighted EM and seven-start GDS equivalence checks passed\n")
