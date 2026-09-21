# Rscript summarize-pca-kmeans-benchmark.R COMPLETED_PCA_RUN WEAK_STRONG_SUMMARY_DIR
args <- commandArgs(TRUE);stopifnot(length(args)==2L)
out <- normalizePath(args[1],mustWork=TRUE)
results <- readRDS(file.path(out,"results.rds"));audit <- readRDS(file.path(out,"design-locked.rds"))
reference <- readRDS(file.path(args[2],"summary.rds"))$panels
stopifnot(length(results)==80L,!any(vapply(results,inherits,logical(1),"try-error")))
km <- do.call(rbind,lapply(results,`[[`,"kmeans"))
native <- do.call(rbind,lapply(results,`[[`,"native"))
stopifnot(nrow(km)==480L,nrow(native)==240L,setequal(km$ID,reference$ID),
 !anyDuplicated(km[c("ID","COMPONENTS")]),all(km$STATUS=="OK"))
selected <- native[native$SELECTED,];stopifnot(nrow(selected)==80L)
selected$POP_ARI <- reference$POP_ARI[match(selected$ID,reference$ID)]
stopifnot(all(abs(selected$FINAL_ARI-reference$STOCK_ARI[match(selected$ID,reference$ID)])<1e-10))
selected$DELTA <- selected$FINAL_ARI-selected$INIT_ARI
km.summary <- do.call(rbind,lapply(split(km,interaction(km$REGIME,km$N,km$L_NOMINAL,km$COMPONENTS,drop=TRUE)),function(x)
 data.frame(REGIME=x$REGIME[1],N=x$N[1],L=x$L_NOMINAL[1],COMPONENTS=x$COMPONENTS[1],
  REPLICATES=nrow(x),MEAN=mean(x$ARI),SD=sd(x$ARI),MIN=min(x$ARI),MAX=max(x$ARI))))
native.summary <- do.call(rbind,lapply(split(selected,interaction(selected$REGIME,selected$N,selected$L_NOMINAL,drop=TRUE)),function(x)
 data.frame(REGIME=x$REGIME[1],N=x$N[1],L=x$L_NOMINAL[1],REPLICATES=nrow(x),
  INIT_MEAN=mean(x$INIT_ARI),FINAL_MEAN=mean(x$FINAL_ARI),POP_MEAN=mean(x$POP_ARI),
  DELTA=mean(x$DELTA),IMPROVED=sum(x$DELTA>1e-10),WORSENED=sum(x$DELTA< -1e-10),UNCHANGED=sum(abs(x$DELTA)<=1e-10))))
stopifnot(all(km.summary$REPLICATES==10L),all(native.summary$REPLICATES==10L))
saveRDS(list(kmeans=km,kmeans.summary=km.summary,native=native,selected=selected,native.summary=native.summary),
 file.path(out,"summary.rds"))
fmt <- function(x) formatC(x,digits=3,format="f")
lines <- c("# PCA and k-means baseline: the same 80 simulated datasets", "",
 "Principal component analysis (PCA) and k-means were evaluated at six prespecified counts of retained principal components (PCs).",
 "Adjusted Rand index (ARI) measures agreement with true groups, adjusted for chance; one means perfect recovery.",
 "Three groups were supplied throughout. This is training recovery, not selecting the number of groups or independent prediction.",
 "", "## Fixed design", "",
 "Reuse all 40 weak-signal and 40 fixed-generating-F=0.02 panels; ten independent datasets per condition.",
 "Same retained markers and sample order as the native stockR/PopCluster comparison, with no newly simulated data.",
 "Standalone PCA centres genotype dosages and divides by sqrt(2p(1-p)), with p estimated from the fitting genotypes.",
 "Retain 2, 5, 10, 20, 50 or 100 components. Score axes are not whitened. Every setting is reported; truth never selects component count.",
 "Three batches of 25 k-means starts per setting; minimum within-cluster sum of squares selects the batch before scoring truth.",
 "Hartigan-Wong algorithm with iteration cap 100; warnings and failed batches retained. Matching 75 starts is not matching runtime or all optimizer controls.",
 "Native stockR uses uncentred, unscaled rotated genotypes, up to 100 components and 25 k-means starts per native call.",
 "Its saved init.grps allow an exact before/after comparison for all 240 native fits; the displayed native start belongs to the run selected by final likelihood.",
 "That native-start summary is a within-run diagnostic, not a separately selected deployable k-means method.",
 "", "## Standalone PCA + k-means: mean ARI", "",
 "| Regime | Samples | Markers | 2 PCs | 5 PCs | 10 PCs | 20 PCs | 50 PCs | 100 PCs |",
 "|---|---:|---:|---:|---:|---:|---:|---:|---:|")
for(regime in c("WEAK","F002")) for(n in c(120,360)) for(l in c(1000,10000)) {
 x <- subset(km.summary,REGIME==regime & N==n & L==l);x<-x[order(x$COMPONENTS),]
 lines <- c(lines,paste0("| ",paste(c(regime,n,l,fmt(x$MEAN)),collapse=" | ")," |"))
}
lines <- c(lines,"", "PCs are principal components. WEAK uses generating Patterson–Price–Reich score 3, with differentiation varying by condition.",
 "F002 uses fixed generating differentiation 0.02. Markers in these tables are nominal: two F002 panels retain 9,999 of 10,000.",
 "", "## Native stockR: exact saved initialization versus final fit", "",
 "| Regime | Samples | Markers | Native start ARI | Final stockR ARI | PopCluster ARI | Mean change | Improved / worsened / unchanged |",
 "|---|---:|---:|---:|---:|---:|---:|---|")
for(i in order(native.summary$REGIME,native.summary$N,native.summary$L)) {
 r <- native.summary[i,]
 lines <- c(lines,paste0("| ",paste(c(r$REGIME,r$N,r$L,fmt(r$INIT_MEAN),fmt(r$FINAL_MEAN),fmt(r$POP_MEAN),fmt(r$DELTA),
 paste(r$IMPROVED,r$WORSENED,r$UNCHANGED,sep=" / ")),collapse=" | ")," |"))
}
lines <- c(lines,"", "ARI changes measure truth recovery, not optimization of the mixture likelihood.",
 "Differences between standalone PCA-k-means and native stockR cannot be attributed solely to mixture refinement: preprocessing and component counts also differ.",
 "", "![Component sensitivity](component-sensitivity.png)","",
 "Bars show across-dataset standard deviation, not confidence intervals; horizontal lines show mean native-method recovery.",
 "The two initializers are labelled separately to avoid treating them as equivalent.",
 "", "## Audit", "",
 paste("PCA-k-means settings:",nrow(km),"; batches with warnings:",sum(km$WARNINGS),"; failed batches:",sum(km$FAILED_BATCHES)),
 "Direct decomposition and reproducibility tests passed. All 80 sample/marker matches and native final ARIs passed checks.",
 "The exact saved native starting labels were read, not reconstructed from a different random seed.",
 "Inputs, saved fits, component variance, identifiers, hashes, seeds and all-start native results are retained in this artifact directory.")
partition.path <- file.path(out,"native-partition-audit.rds")
if(file.exists(partition.path)) {
 partitions <- readRDS(partition.path)
 stopifnot(nrow(partitions)==240L)
 lines <- c(lines,"", "## Direct native partition audit", "",
  "A separate check compares initial and final partitions directly, allowing arbitrary label numbering (not merely equal ARI against truth).")
 for(regime in c("WEAK","F002")) {
   x <- subset(partitions,REGIME==regime)
   lines <- c(lines,paste0(regime,": identical initial/final partitions in ",
    sum(x$SAME_PARTITION & x$SELECTED),"/",sum(x$SELECTED)," selected runs and ",sum(x$SAME_PARTITION),"/",nrow(x)," total runs."))
 }
}
writeLines(lines,file.path(out,"report.md"))
km.summary$DESIGN <- paste0("N = ",km.summary$N,"\n",format(km.summary$L,trim=TRUE,big.mark=",")," markers")
native.summary$DESIGN <- paste0("N = ",native.summary$N,"\n",format(native.summary$L,trim=TRUE,big.mark=",")," markers")
baselines <- rbind(transform(native.summary,METHOD="Native stockR start",ARI=INIT_MEAN),
 transform(native.summary,METHOD="Native stockR final",ARI=FINAL_MEAN),transform(native.summary,METHOD="PopCluster",ARI=POP_MEAN))
plot <- ggplot2::ggplot(km.summary,ggplot2::aes(x=COMPONENTS,y=MEAN))+
 ggplot2::geom_errorbar(ggplot2::aes(ymin=pmax(-.04,MEAN-SD),ymax=pmin(1.04,MEAN+SD)),width=0,colour="grey60")+
 ggplot2::geom_line(ggplot2::aes(colour="Standalone PCA-k-means",linetype="Standalone PCA-k-means"))+ggplot2::geom_point(ggplot2::aes(colour="Standalone PCA-k-means"),size=1.6)+
 ggplot2::geom_hline(data=baselines,ggplot2::aes(yintercept=ARI,colour=METHOD,linetype=METHOD))+
 ggplot2::facet_grid(DESIGN~REGIME,labeller=ggplot2::labeller(REGIME=c(F002="Generating F = 0.02",WEAK="Weak differentiation")))+ggplot2::scale_x_log10(breaks=audit$components,labels=audit$components)+
 ggplot2::coord_cartesian(ylim=c(-.04,1.04))+
 ggplot2::scale_colour_manual(breaks=c("Standalone PCA-k-means","Native stockR start","Native stockR final","PopCluster"),values=c("Standalone PCA-k-means"="#111827","Native stockR start"="#c87516","Native stockR final"="#16755d","PopCluster"="#3167b0"))+
 ggplot2::scale_linetype_manual(breaks=c("Standalone PCA-k-means","Native stockR start","Native stockR final","PopCluster"),values=c("Standalone PCA-k-means"="solid","Native stockR start"="dotted","Native stockR final"="solid","PopCluster"="dashed"))+
 ggplot2::guides(colour=ggplot2::guide_legend(nrow=2,byrow=TRUE),linetype=ggplot2::guide_legend(nrow=2,byrow=TRUE))+
 ggplot2::theme_bw(base_size=10)+ggplot2::theme(legend.position="bottom",strip.text.y=ggplot2::element_text(size=8))+
 ggplot2::labs(x="Retained principal components (log-spaced axis)",y="Mean adjusted Rand index",colour=NULL,linetype=NULL,
 subtitle="Ten datasets per panel; bars show standard deviation")
ggplot2::ggsave(file.path(out,"component-sensitivity.png"),plot,width=8,height=9,dpi=160)
print(km.summary,row.names=FALSE);print(native.summary,row.names=FALSE)
