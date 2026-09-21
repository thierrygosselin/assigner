args <- commandArgs(TRUE); stopifnot(length(args)==1L)
root <- normalizePath(args[1],mustWork=TRUE)
design <- readRDS(file.path(root,"design-and-audit.rds"))$design
all <- readRDS(file.path(root,"results.rds"))
stopifnot(length(all)==nrow(design),all(vapply(all,is.data.frame,logical(1))))
results <- do.call(rbind,all)
results$CORRECT_K <- results$SELECTED_K==results$TRUE_K
results$UNDER_K <- results$SELECTED_K<results$TRUE_K
results$OVER_K <- results$SELECTED_K>results$TRUE_K
summary <- aggregate(cbind(CORRECT_K,UNDER_K,OVER_K,TEST_ARI)~CONDITION+MARKERS+METHOD,results,mean)
fits <- lapply(file.path(root,design$ID,"fits.rds"),readRDS)
source("inst/experimental/benchmark-stockr-paired.R")
source("inst/experimental/diagnose-truth-starts.R")
# Post-selection diagnostic only: truth selects a supplied-K fit for evaluation,
# never a model used in the unknown-K workflow.
known <- lapply(seq_along(fits),function(i) {
  row <- design[i,]; x <- fits[[i]]; tab <- x$native.diagnostics
  ix <- which(tab$K==row$TRUE_K & tab$VALID)
  best <- ix[which.max(tab$LOGLIK[ix])]
  fit <- x$native[[rownames(tab)[best]]]
  d <- readRDS(file.path(root,paste(row$CONDITION,row$REPLICATE,sep="_"),"truth.rds"))
  g <- d$train[,seq_len(row$MARKERS),drop=FALSE]
  keep <- which(apply(g,2,function(z)length(unique(z))>1))
  pred <- .bp_stock_predict(fit,d$test[,seq_len(row$MARKERS),drop=FALSE][,keep,drop=FALSE])
  data.frame(row,NATIVE_KNOWN_K_ARI=.tb_partition(pred$membership,d$group)$ari)
}) |> do.call(what=rbind)
native <- do.call(rbind,lapply(fits,`[[`,"native.diagnostics"))
gl <- do.call(rbind,lapply(fits,function(x)do.call(rbind,x$cv.diagnostics)))
final <- do.call(rbind,lapply(fits,function(x)x$GL$model$starts))
saveRDS(list(results=results,summary=summary,native=native,GL_cv=gl,GL_final=final,
  native.known.K=known),file.path(root,"summary.rds"))
lines <- c("# Fresh proof-of-concept results", "",
  "Fixed strength 100; predictive K selection versus native stockR AIC/BIC.",
  "New frequencies and samples. Three intermediate replicates; two controls per condition.",
  "300 training and 300 test fish; complete independent loci; near-certain GLs.",
  "Generating drift is not exact measured sample FST. Nested marker counts share frequencies.",
  "",sprintf("%d panels completed. %d GL CV starts (%d nonconverged/failed); %d final GL starts (%d nonconverged/failed); %d native fits (%d invalid; %d with warnings).",
    length(all),nrow(gl),sum(!gl$CONVERGED),nrow(final),sum(!final$CONVERGED),
    nrow(native),sum(!native$VALID),sum(nzchar(native$WARNING))),
  "", "## Summary", "",
  "| Condition | Markers | Workflow | Correct K / replicates | Mean test ARI |",
  "| --- | ---: | --- | ---: | ---: |")
for(i in seq_len(nrow(summary))) {
  r <- summary[i,]; sub <- subset(results,CONDITION==r$CONDITION & MARKERS==r$MARKERS & METHOD==r$METHOD)
  lines <- c(lines,sprintf("| %s | %s | %s | %s / %s | %.4f |",r$CONDITION,r$MARKERS,r$METHOD,
    sum(sub$CORRECT_K),nrow(sub),r$TEST_ARI))
}
lines <- c(lines,"","## Individual results","",
  "| Dataset | Workflow | Selected K | Test ARI |",
  "| --- | --- | ---: | ---: |")
for(i in seq_len(nrow(results))) {
  r <- results[i,]; lines <- c(lines,sprintf("| %s | %s | %s | %.4f |",r$ID,r$METHOD,r$SELECTED_K,r$TEST_ARI))
}
lines <- c(lines,"","## Native stockR with generating K supplied (diagnostic only)","",
  "| Dataset | Test ARI |","| --- | ---: |")
for(i in seq_len(nrow(known))) lines <- c(lines,sprintf("| %s | %.4f |",known$ID[i],known$NATIVE_KNOWN_K_ARI[i]))
lines <- c(lines,"","## Limits","",
  "This is a small clean-data proof of concept, not a precise power estimate or an evaluation of uncertain likelihoods, LD, uneven sampling or demographically independent stocks.",
  "The fixed shrinkage strength was selected in an earlier pilot, not tuned here. Workflows use different K-selection criteria; the comparison does not isolate GL encoding.",
  "No claim of universal superiority or default strength recommendation follows from these simulations.")
writeLines(lines,file.path(root,"proof-results.md")); print(summary,row.names=FALSE)
