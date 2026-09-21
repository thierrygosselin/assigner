# Add a true pooled-frequency reference: comparing known K=3 parameters only
# with ESTIMATED K=1 parameters would confound structure and estimation error.
# Rscript .../summarize-truth-starts.R INPUT DIAGNOSTIC NEW_SUMMARY
args<-commandArgs(TRUE)
if(length(args)!=3L) stop("Supply benchmark input, completed diagnostic and NEW summary.")
input<-normalizePath(args[1],mustWork=TRUE)
diagnostic<-normalizePath(args[2],mustWork=TRUE)
output<-args[3]
if(file.exists(output)) stop("Summary output already exists.")
source("R/stock_likelihood_mixture.R")
source("R/stock_mixture_blocks.R")
manifest<-readRDS(file.path(input,"manifest.rds"))
full<-readRDS(file.path(diagnostic,"full-comparisons.rds"))
folds<-readRDS(file.path(diagnostic,"fold-comparisons.rds"))
stopifnot(all(full$TRUTH_START_CONVERGED),all(folds$CONVERGED))
comparisons<-list(); paired<-list()
for(id in full$ID) {
  row<-manifest[manifest$ID==id,]
  truth<-readRDS(file.path(input,"truth",paste0(id,".rds")))
  p<-pmin(pmax(truth$frequency,1e-8),1-1e-8)
  pi<-truth$sample.size/sum(truth$sample.size)
  pooled<-matrix(as.vector(p %*% pi),ncol=1)
  n<-row$SAMPLES; m<-row$MARKERS
  g<-SeqArray::seqOpen(row$GDS,readonly=TRUE)
  raw<-tryCatch(SeqArray::seqGetData(g,"annotation/format/GL"),finally=SeqArray::seqClose(g))
  prepared<-.sm_prepare(array(as.numeric(raw$data),c(n,3,m)),"GL")
  reader<-function(index) list(log.likelihood=prepared$log.likelihood[,,index,drop=FALSE],
    informative=prepared$informative[,index,drop=FALSE])
  k1<-.sm_block_estep(reader,list(seq_len(m)),pooled,1,n)
  k3<-readRDS(file.path(diagnostic,id,"frozen-true-model.rds"))
  f<-folds[folds$ID==id,]
  stopifnot(sum(f$N_TEST)==n)
  comparisons[[id]]<-data.frame(ID=id,
    SEEDED_MINUS_BLIND_K3_PER_SAMPLE=sum(f$TRUTH_START_K3-f$BLIND_K3)/n,
    SEEDED_K3_MINUS_FITTED_K1_PER_SAMPLE=sum(f$TRUTH_START_K3-f$FITTED_K1)/n,
    FROZEN_K3_MINUS_FROZEN_K1_PER_SAMPLE=mean(k3$sample.loglik-k1$sample.loglik),
    FROZEN_K1_MINUS_FITTED_K1_PER_SAMPLE=(sum(k1$sample.loglik)-sum(f$FITTED_K1))/n)
  paired[[id]]<-data.frame(INDIVIDUALS=rownames(truth$genotype),
    FROZEN_K1=k1$sample.loglik,FROZEN_K3=k3$sample.loglik)
}
comparisons<-do.call(rbind,comparisons)
joined<-merge(full,comparisons,by="ID",sort=FALSE)
dir.create(output,recursive=TRUE)
saveRDS(joined,file.path(output,"diagnostic-summary.rds"))
saveRDS(paired,file.path(output,"frozen-paired-scores.rds"))
writeLines(c("Privileged simulation diagnostic; no operational K selection or calibration claim",
  capture.output(print(joined,row.names=FALSE))),file.path(output,"diagnostic-summary.txt"))
print(joined[,c("ID","BLIND_ARI","TRUTH_START_ARI","FROZEN_ARI","OBJECTIVE_GAIN")],row.names=FALSE)
print(comparisons,row.names=FALSE)
cat("Full-fit warnings:",sum(nzchar(full$WARNINGS)),"; fold warnings:",sum(nzchar(folds$WARNINGS)),"\n")
