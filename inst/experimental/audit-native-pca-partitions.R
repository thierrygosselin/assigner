# Confirm partition identity directly, rather than equating identical truth ARI with identical labels.
# Rscript audit-native-pca-partitions.R PCA_BENCHMARK_OUTPUT
args <- commandArgs(TRUE);stopifnot(length(args)==1L)
source("inst/experimental/diagnose-truth-starts.R")
source("inst/experimental/benchmark-stockr-paired.R")
panels <- readRDS(file.path(args[1],"design-locked.rds"))$panels
rows <- lapply(seq_len(nrow(panels)),function(i) {
 r <- panels[i,];p <- file.path(r$SOURCE,r$ID)
 truth <- readRDS(file.path(p,"truth.rds"));keep <- readRDS(file.path(p,"result.rds"))$keep
 native <- readRDS(file.path(p,"stock.rds"))
 do.call(rbind,lapply(seq_along(native$fits),function(s) {
   fit <- native$fits[[s]]
   final <- max.col(.bp_stock_predict(fit,truth$genotype[,keep,drop=FALSE])$membership,ties.method="first")
   data.frame(ID=r$ID,REGIME=r$REGIME,START=s,SELECTED=s==native$selected,
     SAME_PARTITION=abs(.tb_ari(fit$init.grps,final)-1)<1e-10)
 }))
})
result <- do.call(rbind,rows);saveRDS(result,file.path(args[1],"native-partition-audit.rds"))
print(with(subset(result,SELECTED),table(REGIME,SAME_PARTITION)))
print(with(result,table(REGIME,SAME_PARTITION)))
