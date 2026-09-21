# Rscript .../summarize-stockr-marker-benchmark.R RUN NEW_SUMMARY
args<-commandArgs(TRUE)
stopifnot(length(args)==2,!file.exists(args[2]))
input<-normalizePath(args[1]); output<-args[2]
folders<-list.dirs(input,recursive=FALSE)
folders<-folders[startsWith(basename(folders),"fixed_")]
files<-file.path(folders,"result.rds"); files<-files[file.exists(files)]
results<-lapply(files,readRDS)
if(any(!vapply(results,is.data.frame,logical(1)))) stop("Failed jobs: inspect before summarizing.")
results<-do.call(rbind,results)
if(is.null(results)) stop("No completed fixed-K results.")
dir.create(output,recursive=TRUE)
saveRDS(results,file.path(output,"fixed-results.rds"))
writeLines(c("KNOWN K: generating K supplied; not a test of K selection",
  capture.output(print(results,row.names=FALSE))),file.path(output,"fixed-results.txt"))
labels<-c(stockR="stockR EM",GL_0="GL: unpenalized",GL_100="GL: strength 100")
results$METHOD<-factor(results$METHOD,levels=names(labels),labels=labels)
plot.data<-results[results$TRUE_K>1,]
plot.data$CONDITION_LABEL<-paste0("Generating K = ",plot.data$TRUE_K,"; drift = ",plot.data$DRIFT)
plot.data$LINE<-interaction(plot.data$METHOD,plot.data$REPLICATE)
fig<-ggplot2::ggplot(plot.data,ggplot2::aes(MARKERS,ARI,colour=METHOD,group=LINE,
  linetype=factor(REPLICATE)))+ggplot2::geom_line()+ggplot2::geom_point()+
  ggplot2::facet_wrap(~CONDITION_LABEL,ncol=2)+ggplot2::theme_bw(base_size=12)+
  ggplot2::scale_x_log10(breaks=c(1000,3000,5000,10000),labels=c("1,000","3,000","5,000","10,000"))+
  ggplot2::coord_cartesian(ylim=c(-.05,1.02))+
  ggplot2::labs(title="Known-K recovery on independent simulated samples",
    subtitle="300 training + 300 evaluation samples; seven shared blind starts",
    x="Markers (log scale)",y="Adjusted Rand index",colour="Model",linetype="Replicate",
    caption="K supplied, not inferred. Two independent frequency realizations; marker panels nested within each.")
ggplot2::ggsave(file.path(output,"known-K-recovery.png"),fig,width=12,height=8,dpi=150)
cat("Completed panels:",length(unique(results$ID)),"of 40; fitted comparisons:",nrow(results),
  "; converged:",sum(results$CONVERGED),"; failed starts:",sum(results$FAILED_STARTS),"\n")
print(aggregate(ARI~CONDITION+MARKERS+METHOD,results,mean),row.names=FALSE)
