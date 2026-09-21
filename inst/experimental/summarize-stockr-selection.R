# Summarize completed selection jobs; do not omit failed/incomplete jobs silently.
# Rscript .../summarize-stockr-selection.R RUN NEW_OUTPUT
args<-commandArgs(TRUE); stopifnot(length(args)==2,!file.exists(args[2]))
input<-normalizePath(args[1]); output<-args[2]
manifest<-readRDS(file.path(input,"manifest.rds"))
selected<-manifest[manifest$JOB%in%c(7,15,17,23,33,39),]
all.results<-list(); grids<-list(); diagnostics<-list()
for(i in seq_len(nrow(selected))) {
  row<-selected[i,]; d<-file.path(input,paste0("cv_",row$ID))
  x<-readRDS(file.path(d,"result.rds"))
  if(!is.data.frame(x)) stop("Failed selection job: ",row$ID)
  all.results[[i]]<-x
  cv<-do.call(rbind,readRDS(file.path(d,"CV-checkpoint.rds")))
  stopifnot(nrow(cv)==80L,all(is.finite(cv$SCORE)))
  z<-aggregate(SCORE~K+METHOD,cv,mean)
  z$ID<-row$ID; z$DELTA<-z$SCORE-max(z$SCORE)
  z$LABEL<-paste0("K = ",row$TRUE_K,"; drift = ",row$DRIFT,"; markers = ",row$MARKERS)
  grids[[i]]<-z
  counts<-c(starts=0,failed=0,warnings=0)
  for(f in list.files(d,"^fold[0-9]+_K[0-9]+[.]rds$",full.names=TRUE)) {
    for(model in readRDS(f)$models) {
      s<-model$starts
      counts<-counts+c(starts=nrow(s),failed=sum(!s$CONVERGED),
        warnings=if("WARNING"%in%names(s)) sum(nzchar(s$WARNING)) else 0)
    }
  }
  diagnostics[[i]]<-data.frame(ID=row$ID,t(counts))
}
results<-do.call(rbind,all.results); grids<-do.call(rbind,grids)
diagnostics<-do.call(rbind,diagnostics)
dir.create(output,recursive=TRUE)
saveRDS(list(results=results,grid=grids,diagnostics=diagnostics),file.path(output,"selection-summary.rds"))
writeLines(c("Six-panel exploratory predictive K-selection check (not all 40 panels)",
  capture.output(print(results,row.names=FALSE)),capture.output(print(diagnostics,row.names=FALSE))),
  file.path(output,"selection-results.txt"))
fig<-ggplot2::ggplot(grids,ggplot2::aes(K,DELTA,colour=METHOD,group=METHOD))+
  ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::facet_wrap(~LABEL,scales="free_y",ncol=2)+
  ggplot2::scale_x_continuous(breaks=1:4)+ggplot2::theme_bw(base_size=11)+
  ggplot2::labs(title="Five-fold predictive selection: stockR and GL mixture",
    subtitle="Each panel uses 300 training samples; test samples do not select K or strength",
    x="Candidate K",y="Mean held-out log likelihood relative to panel best",colour="Model",
    caption="One realization per displayed condition; GL_2/20/100 denote pooled prior strength. Not native stockR AIC/BIC.")
ggplot2::ggsave(file.path(output,"predictive-K-selection.png"),fig,width=12,height=10,dpi=150)
print(results);print(diagnostics)
