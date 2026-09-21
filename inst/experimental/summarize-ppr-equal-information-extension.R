# Rscript inst/experimental/summarize-ppr-equal-information-extension.R PILOT EXTENSION NEW_REPORT_DIR
args <- commandArgs(TRUE);stopifnot(length(args)==3L,!file.exists(args[3]))
inputs <- normalizePath(args[1:2],mustWork=TRUE)
audits <- lapply(inputs,function(p) readRDS(file.path(p,"design-locked.rds")))
stopifnot(identical(unname(audits[[1]]$binary.hash),unname(audits[[2]]$binary.hash)),
 identical(audits[[1]]$stockR,audits[[2]]$stockR),
 identical(audits[[1]]$source.hash,audits[[2]]$source.hash))
load <- function(p,audit) {
 records <- readRDS(file.path(p,"results.rds"))
 stopifnot(!any(vapply(records,inherits,logical(1),"try-error")))
 tab <- do.call(rbind,records)
 design <- audit$design$panels
 stopifnot(nrow(tab)==nrow(design),setequal(tab$ID,design$ID))
 for(i in seq_len(nrow(tab))) {
   row <- tab[i,];d <- design[match(row$ID,design$ID),]
   stopifnot(row$N==d$N,row$K==d$K,row$L_NOMINAL==d$L,
     row$F_GENERATING==d$F,row$REPLICATE==d$REPLICATE)
   dest <- file.path(p,row$ID)
   stopifnot(readRDS(file.path(dest,"execution.rds"))$status==0L)
   stock <- readRDS(file.path(dest,"stock.rds"))
   choice <- readRDS(file.path(dest,"choices-locked.rds"))
   result <- readRDS(file.path(dest,"result.rds"))
   eligible <- which(stock$diagnostics$CONVERGED & is.finite(stock$diagnostics$LOGLIK))
   stopifnot(choice$stock==eligible[which.max(stock$diagnostics$LOGLIK[eligible])],
     choice$pop==which.max(vapply(result$pop,`[[`,numeric(1),"loglik")))
 }
 tab
}
tab <- rbind(load(inputs[1],audits[[1]]),load(inputs[2],audits[[2]]))
design <- rbind(audits[[1]]$design$panels,audits[[2]]$design$panels)
stopifnot(nrow(tab)==40L,!anyDuplicated(tab$ID),!anyDuplicated(design$SEED),
 all(table(tab$N,tab$L_NOMINAL)==10L),
 all(is.finite(tab$POP_ARI)),all(is.finite(tab$STOCK_ARI)))
tab <- tab[order(tab$N,tab$L_NOMINAL,tab$REPLICATE),]
cells <- split(tab,interaction(tab$N,tab$L_NOMINAL,drop=TRUE))
summaries <- lapply(seq_along(cells),function(i) {
 x <- cells[[i]];delta <- x$POP_ARI-x$STOCK_ARI
 set.seed(923000L+i)
 boot <- replicate(10000,mean(sample(delta,length(delta),replace=TRUE)))
 ci <- unname(quantile(boot,c(.025,.975)))
 data.frame(N=x$N[1],N_PER_GROUP=x$N_PER_GROUP[1],L=x$L_NOMINAL[1],REPLICATES=nrow(x),
  F_GENERATING=x$F_GENERATING[1],F_MEASURED_MEAN=mean(x$F_MEASURED),
  S_MEASURED_MIN=min(x$S_MEASURED),S_MEASURED_MAX=max(x$S_MEASURED),
  POP_MEAN=mean(x$POP_ARI),POP_SD=sd(x$POP_ARI),POP_MIN=min(x$POP_ARI),POP_MAX=max(x$POP_ARI),
  STOCK_MEAN=mean(x$STOCK_ARI),STOCK_SD=sd(x$STOCK_ARI),STOCK_MIN=min(x$STOCK_ARI),STOCK_MAX=max(x$STOCK_ARI),
  PAIRED_DIFFERENCE=mean(delta),CI_LOW=ci[1],CI_HIGH=ci[2],POP_HIGHER=sum(delta>0))
})
summary <- do.call(rbind,summaries)
dir.create(args[3],recursive=TRUE);output <- normalizePath(args[3])
saveRDS(list(panels=tab,summary=summary,inputs=inputs,audits=audits,
 bootstrap=list(replicates=10000,seeds=923000L+seq_along(cells),unit="whole paired dataset")),
 file.path(output,"summary.rds"))
f <- function(x) formatC(x,digits=3,format="f")
lines <- c("# Equal-information comparison: ten datasets per condition", "",
 "Forty independent datasets; three true groups (K = 3) supplied to both methods.",
 "Patterson–Price–Reich (PPR) generating score S = N L F² = 3 in every condition.",
 "N counts total samples, with N/3 per group; markers are independent, complete diploid genotypes.",
 "F is the generating parameter, not measured Weir–Cockerham (1984) differentiation.",
 "Adjusted Rand index (ARI) is agreement with true groups, adjusted for chance; one is perfect recovery.",
 "SD denotes standard deviation across independent datasets.", "",
 "## Results", "",
 "| Samples | Per group | Markers | PopCluster mean ARI (SD) | stockR mean ARI (SD) | Paired difference | Exploratory 95% interval | PopCluster higher |",
 "|---:|---:|---:|---:|---:|---:|---:|---:|")
for(i in seq_len(nrow(summary))) {
 r <- summary[i,]
 lines <- c(lines,paste0("| ",paste(c(r$N,r$N_PER_GROUP,r$L,
 paste0(f(r$POP_MEAN)," (",f(r$POP_SD),")"),paste0(f(r$STOCK_MEAN)," (",f(r$STOCK_SD),")"),
 f(r$PAIRED_DIFFERENCE),paste(f(r$CI_LOW),f(r$CI_HIGH),sep=" to "),paste0(r$POP_HIGHER,"/10")),collapse=" | ")," |"))
}
lines <- c(lines,"",
 "Intervals use 10,000 percentile bootstrap resamples of whole paired datasets within each condition.",
 "They are exploratory, not simultaneous confidence intervals across conditions. Replication was expanded after inspecting the pilot.",
 "", "![Paired recovery](paired-recovery.png)","",
 "## Interpretation limits", "",
 "Equal generating S does not imply identical realised differentiation or identical clustering accuracy.",
 "Sample size, marker count and generating differentiation change together; this is not an isolated causal test of any one factor.",
 "Three native runs per method, selected using their own training objective, not truth. Raw likelihoods across methods are not comparable.",
 "The methods see exactly the same genotype dataset within each pair. Ancestral allele frequencies are uniform between 0.05 and 0.95.",
 "This is training recovery with supplied K, not K selection, independent prediction, or a test of the experimental regularized mixture model.",
 "The experiment is not a universal PPR validation or an assessment of low-coverage genotype likelihoods, linkage, imbalance or missingness.",
 "", "## Per-dataset results", "",
 "| ID | Retained markers | Generating F | Measured FST | Measured S | PopCluster ARI | stockR ARI |",
 "|---|---:|---:|---:|---:|---:|---:|")
for(i in seq_len(nrow(tab))) {
 r <- tab[i,]
 lines <- c(lines,paste0("| ",paste(c(r$ID,r$L,formatC(r$F_GENERATING,digits=7,format="f"),
 formatC(r$F_MEASURED,digits=7,format="f"),f(r$S_MEASURED),f(r$POP_ARI),f(r$STOCK_ARI)),collapse=" | ")," |"))
}
lines <- c(lines,"", "## Audit", "",paste("Successful PopCluster executions:",nrow(tab)),
 paste("stockR iteration-limit/ineligible starts recorded:",sum(tab$STOCK_FAILED_STARTS)),
 paste("stockR starts with captured warnings:",sum(tab$STOCK_WARNINGS)),
 "Program binary, stockR version and frozen runner/source hashes match between pilot and extension.",
 "Selected runs were checked against recorded objectives. Genotype/sample order checks were enforced by the runner.",
 "PopCluster's execution also produces admixture output; only its clustering labels are evaluated.",
 "Runtime is not used to compare optimizer efficiency.","",paste("Source artifacts:",paste(inputs,collapse="; ")))
writeLines(lines,file.path(output,"report.md"))
long <- rbind(transform(tab,METHOD="stockR",ARI=STOCK_ARI),transform(tab,METHOD="PopCluster",ARI=POP_ARI))
long$METHOD <- factor(long$METHOD,levels=c("stockR","PopCluster"))
long$SAMPLES <- paste0(long$N," samples (",long$N_PER_GROUP," per group)")
long$MARKERS <- paste0(format(long$L_NOMINAL,trim=TRUE,big.mark=",")," markers")
plot <- ggplot2::ggplot(long,ggplot2::aes(x=METHOD,y=ARI,group=ID))+
 ggplot2::geom_line(alpha=.4,colour="grey50")+ggplot2::geom_point(ggplot2::aes(colour=METHOD),size=2)+
 ggplot2::facet_grid(SAMPLES~MARKERS)+ggplot2::scale_y_continuous(limits=c(min(0,min(long$ARI)),1))+
 ggplot2::theme_bw()+ggplot2::theme(legend.position="none")+
 ggplot2::labs(x=NULL,y="Adjusted Rand index (ARI)",
 title="Recovery with three groups supplied: paired independent datasets",
 subtitle="Generating Patterson–Price–Reich score = 3; ten datasets per condition")
ggplot2::ggsave(file.path(output,"paired-recovery.png"),plot,width=10,height=7,dpi=150)
print(summary,row.names=FALSE)
