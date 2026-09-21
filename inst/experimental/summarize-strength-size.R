args <- commandArgs(TRUE);stopifnot(length(args)==1L)
root <- normalizePath(args[1],mustWork=TRUE)
audit <- readRDS(file.path(root,"audit.rds")); design <- audit$design
all.results <- readRDS(file.path(root,"results.rds"))
stopifnot(length(all.results)==15L,all(vapply(all.results,is.data.frame,logical(1))))
results <- do.call(rbind,all.results)
results$CORRECT_K <- results$K==results$TRUE_K
summary <- aggregate(cbind(CORRECT_K,TEST_ARI)~CONDITION+N+METHOD,results,mean)
fits <- lapply(file.path(root,design$ID,"fits.rds"),readRDS)
cv.diagnostics <- do.call(rbind,lapply(fits,function(x)do.call(rbind,x$diagnostics)))
native.diagnostics <- do.call(rbind,lapply(fits,`[[`,"native.diagnostics"))
final.diagnostics <- do.call(rbind,lapply(fits,function(x)do.call(rbind,lapply(x$final,function(z)z$model$diagnostics))))
grid <- do.call(rbind,lapply(seq_along(fits),function(i) {
  g <- fits[[i]]$grid
  data.frame(design[rep(i,nrow(g)),,drop=FALSE],g,row.names=NULL)
}))
init.warnings <- unlist(lapply(fits,`[[`,"initialization.warnings"))
saveRDS(list(results=results,summary=summary,grid=grid,CV=cv.diagnostics,
  native=native.diagnostics,final=final.diagnostics,initialization.warnings=init.warnings),file.path(root,"summary.rds"))
lines <- c("# Sample-size and regularization results", "",
  "3,000 markers held constant. Fresh realizations; 300 independent test fish.",
  "Structured populations: 100/200/300 fish per group. Nulls use matching total N.",
  "Three structured realizations and two null realizations, nested across sample sizes.",
  "Seven-start initialization unchanged. Fivefold CV; strengths 20/50/100/200/500.",
  "",sprintf("Completed 15 panels: %d CV starts (%d nonconverged/failed), %d final starts (%d nonconverged/failed), %d native fits (%d invalid, %d with warnings); %d recorded CV initialization warnings.",
    nrow(cv.diagnostics),sum(!cv.diagnostics$CONVERGED),nrow(final.diagnostics),sum(!final.diagnostics$CONVERGED),
    nrow(native.diagnostics),sum(!native.diagnostics$VALID),sum(nzchar(native.diagnostics$WARNING)),length(init.warnings)),
  "", "## Recovery", "",
  "| Condition | Total training N | Workflow | Correct K / replicates | Mean test ARI |",
  "| --- | ---: | --- | ---: | ---: |")
for(i in seq_len(nrow(summary))) {
  r <- summary[i,]; x <- subset(results,CONDITION==r$CONDITION & N==r$N & METHOD==r$METHOD)
  lines <- c(lines,sprintf("| %s | %s | %s | %s / %s | %.4f |",r$CONDITION,r$N,r$METHOD,sum(x$CORRECT_K),nrow(x),r$TEST_ARI))
}
lines <- c(lines,"","## Individual choices","",
  "| Dataset | Workflow | K | Strength | Test ARI |",
  "| --- | --- | ---: | ---: | ---: |")
for(i in seq_len(nrow(results))) {
  r <- results[i,];lines<-c(lines,sprintf("| %s | %s | %s | %s | %.4f |",r$ID,r$METHOD,r$K,
    if(is.na(r$STRENGTH)) "—" else r$STRENGTH,r$TEST_ARI))
}
lines<-c(lines,"","## Interpretation limits","",
  "A tuned strength of 500 reaches the grid boundary; it is not an established optimum or a recommended default.",
  "These are small clean-data comparisons, not a power estimate. Nested sample sizes share population frequencies and test cohorts.",
  "The fast calculation is validated only for complete near-certain GLs; it is not a new production likelihood engine.",
  "No model settings were changed during the experiment. Native AIC/BIC and GL predictive selection answer different selection questions.")
writeLines(lines,file.path(root,"strength-size-results.md"))
print(summary,row.names=FALSE);print(subset(results,METHOD=="tuned"),row.names=FALSE)
