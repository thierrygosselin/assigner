# Rscript .../summarize-joint-validation.R COMPLETED_RUN NEW_SUMMARY
args<-commandArgs(TRUE)
if(length(args)!=2L||file.exists(args[2])) stop("Supply completed run and NEW summary.")
input<-normalizePath(args[1],mustWork=TRUE); output<-args[2]
results<-readRDS(file.path(input,"results.rds"))
manifest<-readRDS(file.path(input,"manifest.rds"))
evaluation<-readRDS(file.path(input,"evaluation.rds"))
grid<-do.call(rbind,lapply(names(results),function(id) data.frame(ID=id,results[[id]]$grid)))
folds<-do.call(rbind,lapply(results,`[[`,"folds"))
grid<-merge(grid,manifest[,c("ID","CONDITION","PANEL","TRAIN_N","MARKERS")],by="ID",sort=FALSE)
grid$DELTA<-ave(grid$SCORE,interaction(grid$ID,grid$V),FUN=function(x) x-max(x))
grid$LABEL<-paste(ifelse(grid$CONDITION=="weak_three","Three populations: drift 0.01","One-population null"),
  paste0("N = ",grid$TRAIN_N,", markers = ",grid$MARKERS),sep="\n")
grid$FOLDS<-paste0(grid$V,"-fold")
dir.create(output,recursive=TRUE)
saveRDS(grid,file.path(output,"cv-grid.rds"))
saveRDS(folds,file.path(output,"cv-folds.rds"))
writeLines(c("Training-selected models evaluated on a separate cohort; no selection using test outcomes",
  capture.output(print(evaluation,row.names=FALSE))),file.path(output,"evaluation.txt"))
cat("CV fits:",nrow(folds),"; converged:",sum(folds$CONVERGED),"\n")
cat("Fold warnings:",sum(nzchar(folds$WARNINGS)),"; errors:",sum(nzchar(folds$ERROR)),"\n")
print(evaluation,row.names=FALSE)
if(requireNamespace("ggplot2",quietly=TRUE)) {
  fig<-ggplot2::ggplot(grid,ggplot2::aes(K,DELTA,colour=factor(STRENGTH),group=STRENGTH))+
    ggplot2::geom_hline(yintercept=0,colour="grey70")+ggplot2::geom_line()+ggplot2::geom_point()+
    ggplot2::facet_grid(LABEL~FOLDS,scales="free_y")+
    ggplot2::scale_x_continuous(breaks=1:4)+ggplot2::theme_bw(base_size=11)+
    ggplot2::labs(title="Training cross-validation: K and regularization",
      subtitle="Scores relative to the best candidate within each panel and fold policy",
      x="Candidate K",y="Mean held-out log likelihood minus best",colour="Strength",
      caption="One population-frequency realization per condition; independent test data not used for these choices")
  ggplot2::ggsave(file.path(output,"training-CV.png"),fig,width=11,height=11,dpi=150)
}
