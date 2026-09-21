# Rscript this_script.R WORKSPACE NEW_RUN
args <- commandArgs(TRUE); stopifnot(length(args)==2L)
workspace <- normalizePath(args[1]); run <- normalizePath(args[2])
old <- readRDS(file.path(workspace,"artifacts/ppr-popcluster-matched-20260921/comparison.rds"))
old <- old[old$F_GENERATING==.001,]
fst <- readRDS(file.path(workspace,"artifacts/stockr-marker-benchmark-20260916/measured-WC84.rds"))
old.fst <- fst$WC84[match(old$ID,fst$ID)]
stopifnot(nrow(old)==2L,all(is.finite(old.fst)))
initial <- data.frame(REPLICATE=1:2,F_MEASURED=old.fst,
  S_MEASURED=300*10000*old.fst^2,POP_ARI=old$POP_TRAIN_ARI,STOCK_ARI=old$STOCK_TRAIN_ARI)
new <- readRDS(file.path(run,"results.rds"))
stopifnot(length(new)==8L,all(vapply(new,is.data.frame,logical(1))))
new <- do.call(rbind,new)
stopifnot(identical(sort(new$REPLICATE),3:10),all(new$N==300L),all(new$L==10000L))
tab <- rbind(initial,new[,names(initial)])
tab$DIFFERENCE <- tab$POP_ARI-tab$STOCK_ARI
set.seed(21092026)
boot <- replicate(10000,mean(sample(tab$DIFFERENCE,replace=TRUE)))
summary <- data.frame(METHOD=c("PopCluster","stockR"),
 MEAN=c(mean(tab$POP_ARI),mean(tab$STOCK_ARI)),
 MEDIAN=c(median(tab$POP_ARI),median(tab$STOCK_ARI)),
 MIN=c(min(tab$POP_ARI),min(tab$STOCK_ARI)),MAX=c(max(tab$POP_ARI),max(tab$STOCK_ARI)))
result <- list(replicates=tab,summary=summary,
 mean.difference=mean(tab$DIFFERENCE),paired.bootstrap.interval=quantile(boot,c(.025,.975)),
 note="Exploratory paired dataset bootstrap, ten replicates, one condition; not broad superiority or K-selection evidence.")
saveRDS(result,file.path(run,"summary.rds"))
print(tab,row.names=FALSE); print(summary,row.names=FALSE)
print(result[c("mean.difference","paired.bootstrap.interval","note")])
