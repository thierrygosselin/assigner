# Post-primary exploratory initialization sensitivity; never overwrite baseline.
# Rscript .../run-all-marker-sensitivity.R BENCHMARK NEW_OUTPUT
args<-commandArgs(TRUE); stopifnot(length(args)==2,!file.exists(args[2]))
input<-normalizePath(args[1]); output<-args[2]
paths<-c("R/stock_likelihood_mixture.R","R/stock_mixture_blocks.R",
  "inst/experimental/benchmark-streaming-engine.R","inst/experimental/benchmark-stockr-paired.R",
  "inst/experimental/benchmark-all-marker-starts.R","inst/experimental/diagnose-truth-starts.R",
  "../grur/R/simulation_benchmark.R")
for(p in paths) source(p)
dir.create(output,recursive=TRUE)
saveRDS(list(sources=tools::md5sum(c(paths,"inst/experimental/run-all-marker-sensitivity.R")),
  session=sessionInfo(),reason="Exploratory all-marker initialization after observing weak-case failures",
  replicates=1:2,markers=10000L,strength=c(0,100)),file.path(output,"audit.rds"))
manifest<-readRDS(file.path(input,"manifest.rds"))
rows<-manifest[manifest$CONDITION=="three_001"&manifest$MARKERS==10000,]
run<-function(q) {
  row<-rows[q,]; src<-file.path(input,paste(row$CONDITION,row$REPLICATE,sep="_"))
  d<-readRDS(file.path(src,"truth.rds")); g<-d$train
  keep<-which(apply(g,2,function(x)length(unique(x))>1)); g<-g[,keep,drop=FALSE]
  context<-.bs_context(file.path(src,"train.gds"),variant.ids=keep)
  on.exit(context$close())
  starts<-.bp_all_starts(g,3,120000+row$JOB)
  models<-list(stockR=.bp_stock(g,3,starts),GL_0=.bp_gl(context,g,3,starts,0),
    GL_100=.bp_gl(context,g,3,starts,100))
  saveRDS(models,file.path(output,paste0("fits_",q,".rds")))
  test<-.bs_context(file.path(src,"test.gds"),variant.ids=keep)
  on.exit(test$close(),add=TRUE)
  result<-lapply(names(models),function(name) {
    fit<-models[[name]]$fit
    if(is.null(fit)) stop("No converged fit; inspect saved diagnostics.")
    pred<-if(name=="stockR") .bp_stock_predict(fit,d$test[,keep,drop=FALSE]) else
      .sm_block_estep(test$reader,test$blocks,fit$frequency,fit$proportion,test$n)
    saveRDS(pred,file.path(output,paste0("prediction_",q,"_",name,".rds")))
    z<-.tb_partition(pred$membership,d$group)
    data.frame(REPLICATE=row$REPLICATE,METHOD=name,ARI=z$ari,
      TEST_LOGLIK=mean(pred$sample.loglik),FAILED_STARTS=sum(!models[[name]]$starts$CONVERGED))
  })
  do.call(rbind,result)
}
result<-parallel::mclapply(seq_len(nrow(rows)),run,mc.cores=2,mc.set.seed=FALSE)
saveRDS(result,file.path(output,"results.rds"));print(result)
