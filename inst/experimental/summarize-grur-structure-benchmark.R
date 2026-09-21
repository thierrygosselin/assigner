# Read-only analysis of completed pilot and continuation fits.
# Rscript .../summarize-grur-structure-benchmark.R INPUT CONTINUATION NEW_SUMMARY
args<-commandArgs(TRUE)
if(length(args)!=3L) stop("Supply input, continuation, and NEW summary directories.")
input<-normalizePath(args[1],mustWork=TRUE)
extended<-normalizePath(args[2],mustWork=TRUE)
output<-args[3]
if(file.exists(output)) stop("Summary directory already exists.")
manifest<-readRDS(file.path(input,"manifest.rds"))
results<-readRDS(file.path(extended,"all-results.rds"))
stopifnot(identical(manifest$ID,results$ID))
pilot<-readRDS(file.path(input,"pilot-results.rds"))
scores<-diagnostics<-starts<-information<-confidence<-list()
for(i in seq_len(nrow(manifest))) {
  id<-manifest$ID[i]
  path<-if(id %in% names(pilot)) file.path(input,paste0(id,"_CV.rds")) else
    file.path(extended,id,"CV.rds")
  if(!file.exists(path)) next # Errors remain explicit in results, not dropped there.
  cv<-readRDS(path)
  scores[[id]]<-cbind(manifest[rep(i,nrow(cv$summary)),c("ID","SCENARIO","REPLICATE")],cv$summary)
  diagnostics[[id]]<-cbind(ID=id,cv$diagnostics)
  fit.path<-if(id %in% names(pilot)) file.path(input,paste0(id,"_fixed-K-fit.rds")) else
    file.path(extended,id,"fixed-K-fit.rds")
  if(file.exists(fit.path)) {
    fit<-readRDS(fit.path)
    starts[[id]]<-cbind(ID=id,fit$start.diagnostics)
    information[[id]]<-cbind(ID=id,fit$likelihood.information.audit)
    maximum<-apply(fit$membership,1,max)
    confidence[[id]]<-data.frame(ID=id,K=ncol(fit$membership),
      MEAN_MAX_MEMBERSHIP=mean(maximum),FRACTION_AT_LEAST_0.95=mean(maximum>=.95),
      NEAR_BEST_STARTS=fit$near.best.starts)
  }
}
scores<-do.call(rbind,scores); diagnostics<-do.call(rbind,diagnostics)
# Differences within a dataset only: likelihood constants and marker counts
# make raw scores across datasets unsuitable as a performance comparison.
scores$DELTA_PER_SAMPLE<-ave(scores$MEAN_LOGLIK,scores$ID,FUN=function(x) {
  if(any(!is.finite(x))) rep(NA_real_,length(x)) else x-max(x)
})
dir.create(output,recursive=TRUE)
saveRDS(scores,file.path(output,"predictive-scores.rds"))
saveRDS(diagnostics,file.path(output,"fold-diagnostics.rds"))
saveRDS(do.call(rbind,starts),file.path(output,"start-diagnostics.rds"))
saveRDS(do.call(rbind,information),file.path(output,"information-audit.rds"))
saveRDS(do.call(rbind,confidence),file.path(output,"conditional-membership-audit.rds"))
writeLines(c("All cross-validation folds (warnings/errors retained)",
  capture.output(print(diagnostics,row.names=FALSE))),file.path(output,"fold-diagnostics.txt"))
joined<-cbind(manifest[,c("SCENARIO","REPLICATE","SAMPLE_SIZES","MARKERS","UNINFORMATIVE")],results)
writeLines(c("Completed benchmark results; only two replicates per scenario",
  capture.output(print(joined,row.names=FALSE))),file.path(output,"results.txt"))
if(requireNamespace("ggplot2",quietly=TRUE)) {
  figure<-ggplot2::ggplot(scores,ggplot2::aes(K,DELTA_PER_SAMPLE,
    colour=factor(REPLICATE),group=REPLICATE))+
    ggplot2::geom_hline(yintercept=0,colour="grey70")+
    ggplot2::geom_line()+ggplot2::geom_point()+
    ggplot2::facet_wrap(~SCENARIO,scales="free_y",ncol=3)+
    ggplot2::scale_x_continuous(breaks=1:5)+ggplot2::theme_bw(base_size=11)+
    ggplot2::labs(title="Held-out prediction across candidate K",
      subtitle="Zero = best within each dataset; vertical scales differ between panels",
      x="Candidate K",y="Mean held-out log likelihood minus best",colour="Replicate",
      caption="Two simulation replicates; experimental settings, not a power or calibration study")
  ggplot2::ggsave(file.path(output,"predictive-K.png"),figure,width=12,height=11,dpi=160)
}
cat("Dataset errors:",sum(nzchar(results$ERROR)),"\n")
cat("Incomplete CV datasets:",sum(!results$COMPLETE),"\n")
cat("Fold warnings:",sum(!is.na(diagnostics$WARNINGS) & nzchar(diagnostics$WARNINGS)),"\n")
cat("Fold errors:",sum(!is.na(diagnostics$ERROR) & nzchar(diagnostics$ERROR)),"\n")
print(joined[,c("SCENARIO","REPLICATE","ARI","PREDICTIVE_K","CONVERGED","COMPLETE")],row.names=FALSE)
