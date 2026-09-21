# Rscript inst/experimental/summarize-fixed-f002.R FIXED_OUTPUT WEAK_SUMMARY_DIR NEW_REPORT_DIR
args <- commandArgs(TRUE);stopifnot(length(args)==3L,!file.exists(args[3]))
fixed <- normalizePath(args[1],mustWork=TRUE)
weak <- readRDS(file.path(args[2],"summary.rds"))
audit <- readRDS(file.path(fixed,"design-locked.rds"))
stopifnot(identical(audit$source.hash,weak$audits[[1]]$source.hash),
 identical(unname(audit$binary.hash),unname(weak$audits[[1]]$binary.hash)),
 identical(audit$stockR,weak$audits[[1]]$stockR))
records <- readRDS(file.path(fixed,"results.rds"))
stopifnot(!any(vapply(records,inherits,logical(1),"try-error")))
tab <- do.call(rbind,records);design <- audit$design$panels
stopifnot(nrow(tab)==40L,nrow(design)==40L,setequal(tab$ID,design$ID),
 !anyDuplicated(tab$ID),all(table(tab$N,tab$L_NOMINAL)==10L),
 all(tab$K==3L),all(tab$F_GENERATING==.02),
 all(is.finite(tab$POP_ARI)),all(is.finite(tab$STOCK_ARI)))
for(i in seq_len(nrow(tab))) {
 row <- tab[i,];d <- design[match(row$ID,design$ID),]
 stopifnot(row$N==d$N,row$L_NOMINAL==d$L,row$REPLICATE==d$REPLICATE)
 p <- file.path(fixed,row$ID)
 stopifnot(readRDS(file.path(p,"execution.rds"))$status==0L)
 stock <- readRDS(file.path(p,"stock.rds"));choice <- readRDS(file.path(p,"choices-locked.rds"))
 result <- readRDS(file.path(p,"result.rds"))
 eligible <- which(stock$diagnostics$CONVERGED & is.finite(stock$diagnostics$LOGLIK))
 stopifnot(choice$stock==eligible[which.max(stock$diagnostics$LOGLIK[eligible])],
  choice$pop==which.max(vapply(result$pop,`[[`,numeric(1),"loglik")))
}
all <- rbind(transform(weak$panels,REGIME="Weak: generating PPR score = 3"),
 transform(tab,REGIME="Stronger: generating F = 0.02"))
all$REGIME <- factor(all$REGIME,levels=c("Weak: generating PPR score = 3","Stronger: generating F = 0.02"))
summary <- do.call(rbind,lapply(split(all,interaction(all$REGIME,all$N,all$L_NOMINAL,drop=TRUE)),function(x) {
 data.frame(REGIME=as.character(x$REGIME[1]),N=x$N[1],N_PER_GROUP=x$N_PER_GROUP[1],L=x$L_NOMINAL[1],
  RETAINED_MIN=min(x$L),RETAINED_MAX=max(x$L),REPLICATES=nrow(x),F_GENERATING=x$F_GENERATING[1],
  F_MEASURED=mean(x$F_MEASURED),F_MEASURED_MIN=min(x$F_MEASURED),F_MEASURED_MAX=max(x$F_MEASURED),
  S_GENERATING=x$S_GENERATING[1],STOCK_MEAN=mean(x$STOCK_ARI),STOCK_SD=sd(x$STOCK_ARI),
  STOCK_MIN=min(x$STOCK_ARI),STOCK_MAX=max(x$STOCK_ARI),
  POP_MEAN=mean(x$POP_ARI),POP_SD=sd(x$POP_ARI),POP_MIN=min(x$POP_ARI),POP_MAX=max(x$POP_ARI))
}))
summary <- summary[order(summary$REGIME,summary$N,summary$L),];rownames(summary)<-NULL
dir.create(args[3],recursive=TRUE);out <- normalizePath(args[3])
saveRDS(list(summary=summary,panels=all,audit=audit,weak.sources=weak$inputs),file.path(out,"summary.rds"))
jsonlite::write_json(list(summary=summary,panels=all,
 audit=list(stockR=audit$stockR,popcluster="1.7.0.0 (build 20250715)",
 starts=3,stock.warnings=sum(all$STOCK_WARNINGS),stock.failed.starts=sum(all$STOCK_FAILED_STARTS),
 datasets=nrow(all),source.hashes.match=TRUE)),file.path(out,"summary.json"),dataframe="rows",pretty=TRUE,digits=12,auto_unbox=TRUE)
f <- function(x) formatC(x,digits=3,format="f")
lines <- c("# Native stockR and PopCluster: weak and stronger differentiation", "",
 "Eighty independent simulated datasets, ten per condition. Three true groups (K = 3) supplied to both methods.",
 "Adjusted Rand index (ARI): one means perfect group recovery; near zero means chance-level agreement.",
 "SD denotes standard deviation across datasets. F generating differs from measured Weir–Cockerham (1984) FST.",
 "PPR denotes the Patterson–Price–Reich reference score S = N L F².","",
 "| Regime | Total samples | Per group | Markers | Generating F | Mean measured FST | stockR ARI (SD) | PopCluster ARI (SD) |",
 "|---|---:|---:|---:|---:|---:|---:|---:|")
for(i in seq_len(nrow(summary))) {
 r <- summary[i,]
 lines <- c(lines,paste0("| ",paste(c(r$REGIME,r$N,r$N_PER_GROUP,r$L,
 formatC(r$F_GENERATING,digits=6,format="f"),formatC(r$F_MEASURED,digits=6,format="f"),
 paste0(f(r$STOCK_MEAN)," (",f(r$STOCK_SD),")"),paste0(f(r$POP_MEAN)," (",f(r$POP_SD),")")),collapse=" | ")," |"))
}
lines <- c(lines,"", "![Paired recovery](weak-strong-recovery.png)","",
 "Each method receives the identical complete, independent-marker genotype matrix within each dataset.",
 "Three native runs; selected by each program's training objective before evaluating truth. No truth-guided initialization or run selection.",
 "PopCluster's hard-partition and stockR's marginal mixture objectives differ: raw likelihoods are not compared between methods.",
 "Weak-signal differentiation varies to hold generating S=3; stronger-signal generating F stays fixed at 0.02.",
 "Datasets across regimes are independent, not paired. Lines in the figure connect methods only within one dataset.",
 "This tests group recovery with K supplied, not selecting K, held-out prediction, or genotype-likelihood modelling.",
 "The experimental assigner mixture implementation is not evaluated in these attachments.",
 "The simulator draws ancestral allele frequencies uniformly between 0.05 and 0.95. All three true groups are balanced.",
 "Generating F=0.02 does not force each realised measured FST to equal 0.02.",
 "All executions and selected objectives passed audit; the same binary, stockR version and frozen source hashes were verified across regimes.",
 paste("Captured stockR warnings:",sum(all$STOCK_WARNINGS),"; ineligible/iteration-limit starts:",sum(all$STOCK_FAILED_STARTS)),
 "All nominal marker counts and retained counts are recorded in summary.rds and summary.json.")
writeLines(lines,file.path(out,"report.md"))
long <- rbind(transform(all,METHOD="stockR",ARI=STOCK_ARI),transform(all,METHOD="PopCluster",ARI=POP_ARI))
long$METHOD <- factor(long$METHOD,levels=c("stockR","PopCluster"))
long$DESIGN <- paste0("N = ",long$N,"\n",format(long$L_NOMINAL,trim=TRUE,big.mark=",")," markers")
plot <- ggplot2::ggplot(long,ggplot2::aes(x=METHOD,y=ARI,group=ID))+
 ggplot2::geom_line(colour="grey65",alpha=.65)+ggplot2::geom_point(ggplot2::aes(colour=METHOD),size=1.7,alpha=.8)+
 ggplot2::facet_grid(DESIGN~REGIME)+ggplot2::scale_y_continuous(limits=c(-.035,1.035),breaks=c(0,.5,1))+
 ggplot2::theme_bw(base_size=10)+ggplot2::theme(legend.position="none",strip.text.y=ggplot2::element_text(size=8))+
 ggplot2::labs(x=NULL,y="Adjusted Rand index (ARI)",subtitle="Ten datasets per panel; lines pair methods on identical genotypes")
ggplot2::ggsave(file.path(out,"weak-strong-recovery.png"),plot,width=7.3,height=8.1,dpi=200)
print(summary,row.names=FALSE)
