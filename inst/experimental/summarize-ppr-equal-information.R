# Rscript inst/experimental/summarize-ppr-equal-information.R COMPLETED_OUTPUT
args <- commandArgs(TRUE);stopifnot(length(args)==1L)
output <- normalizePath(args[1],mustWork=TRUE)
records <- readRDS(file.path(output,"results.rds"))
stopifnot(!any(vapply(records,inherits,logical(1),"try-error")))
tab <- do.call(rbind,records)
design <- readRDS(file.path(output,"design-locked.rds"))$design$panels
stopifnot(nrow(tab)==nrow(design),setequal(tab$ID,design$ID),
 all(is.finite(tab$POP_ARI)),all(is.finite(tab$STOCK_ARI)))
means <- aggregate(cbind(POP_ARI,STOCK_ARI,F_MEASURED,S_MEASURED)~N+L,tab,mean)
saveRDS(list(per.panel=tab,means=means),file.path(output,"summary.rds"))
fmt <- function(x) formatC(x,digits=5,format="fg",flag="#")
lines <- c("# Equal-information pilot: native stockR and PopCluster",
 "", "K is the supplied number of groups (three). Each condition has only two independent datasets.",
 "ARI is the adjusted Rand index: one is perfect agreement, near zero is chance-level agreement.",
 "PPR is the Patterson–Price–Reich reference score S = N L F², not a recovery guarantee.",
 "N counts total training samples; all three true groups have equal sample sizes.",
 "F generating is not measured Weir–Cockerham (1984) differentiation.",
 "", "| N | Per group | Nominal markers | Retained markers | Replicate | F generating | F measured | S generating | S measured | PopCluster ARI | stockR ARI |",
 "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for(i in seq_len(nrow(tab))) {
 r <- tab[i,]
 lines <- c(lines,paste0("| ",paste(c(r$N,r$N_PER_GROUP,r$L_NOMINAL,r$L,r$REPLICATE,
 fmt(r$F_GENERATING),fmt(r$F_MEASURED),fmt(r$S_GENERATING),fmt(r$S_MEASURED),
 fmt(r$POP_ARI),fmt(r$STOCK_ARI)),collapse=" | ")," |"))
}
lines <- c(lines,"", "Same complete independent-marker genotypes supplied to both programs; three runs each.",
 "Runs selected by their own training objective before evaluation against truth.",
 "This compares training group recovery only, not selection of K or held-out prediction.",
 "No automatic exclusions, uncertain genotypes, linkage or imbalance were tested.",
 "PopCluster execution includes the extra admixture output observed in the earlier pilot; only clustering labels were scored.",
 "Two replicates per condition are a feasibility pilot, not a power curve or general ranking.")
writeLines(lines,file.path(output,"report.md"))
print(tab,row.names=FALSE);print(means,row.names=FALSE)
