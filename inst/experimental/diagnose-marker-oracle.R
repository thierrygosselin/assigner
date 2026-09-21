# Privileged known-frequency diagnostic, NOT a deployable estimator or K selector.
# Rscript .../diagnose-marker-oracle.R RUN NEW_OUTPUT
args<-commandArgs(TRUE); stopifnot(length(args)==2,!file.exists(args[2]))
input<-normalizePath(args[1]); output<-args[2]
source("inst/experimental/benchmark-vectorized-engine.R")
source("inst/experimental/benchmark-streaming-engine.R")
source("inst/experimental/diagnose-truth-starts.R")
dir.create(output,recursive=TRUE)
results<-lapply(1:2,function(r) {
  src<-file.path(input,paste0("three_001_",r)); d<-readRDS(file.path(src,"truth.rds"))
  context<-.bs_context(file.path(src,"test.gds")); on.exit(context$close())
  three<-.bv_estep(context$reader,context$blocks,d$frequency,rep(1/3,3),context$n)
  one<-.bv_estep(context$reader,context$blocks,matrix(rowMeans(d$frequency),ncol=1),1,context$n)
  saveRDS(list(three=three,one=one),file.path(output,paste0("oracle_",r,".rds")))
  data.frame(REPLICATE=r,ARI=.tb_partition(three$membership,d$group)$ari,
    LOGLIK_GAIN_VS_TRUE_POOLED=mean(three$sample.loglik-one$sample.loglik))
})
result<-do.call(rbind,results);saveRDS(result,file.path(output,"results.rds"));print(result)
