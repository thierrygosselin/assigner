args <- commandArgs(TRUE)
stopifnot(length(args)%in%1:2)
root <- normalizePath(args[1],mustWork=TRUE)
all <- readRDS(file.path(root,"results.rds"))
stopifnot(length(all)==40L,!any(vapply(all,inherits,logical(1),"try-error")))
results <- do.call(rbind,all)
files <- file.path(root,unique(results$ID),"result.rds")
fits <- lapply(files,readRDS)
diagnostics <- do.call(rbind,lapply(fits,`[[`,"diagnostics"))
summary <- aggregate(cbind(SELECTED_K,TEST_ARI)~CONDITION+MARKERS+PROTOCOL+CRITERION,
  results,mean)
selection <- subset(results,CRITERION!="known_K")
selection$CORRECT <- selection$SELECTED_K==selection$TRUE_K
selection$UNDER <- selection$SELECTED_K<selection$TRUE_K
selection$OVER <- selection$SELECTED_K>selection$TRUE_K
rates <- aggregate(cbind(CORRECT,UNDER,OVER)~CONDITION+MARKERS+PROTOCOL+CRITERION,
  selection,mean)
stopifnot(nrow(diagnostics)==400L)
saveRDS(list(results=results,diagnostics=diagnostics,summary=summary,rates=rates),
  file.path(root,"summary.rds"))
lines <- c("# Native stockR results", "",
  sprintf("40 panels; %d native fits; %d trace-classified nonconverged/error fits; %d fits with warnings.",
    nrow(diagnostics),sum(!diagnostics$CONVERGED),sum(nzchar(diagnostics$WARNING))),
  "", "Two frequency realizations per condition; panels are nested. No power claims.",
  "Native single = first call. Native best3 = best converged of three native calls per K.",
  "AIC and BIC are returned by stockR itself; each tested at K = 1–4.",
  "", "## All selections and independent-test recovery", "",
  "| Dataset | Protocol | Criterion | Selected K | Test ARI |",
  "| --- | --- | --- | ---: | ---: |")
for(i in seq_len(nrow(results))) {
  r <- results[i,]
  lines <- c(lines,sprintf("| %s | %s | %s | %s | %.4f |",r$ID,r$PROTOCOL,
    r$CRITERION,r$SELECTED_K,r$TEST_ARI))
}
writeLines(lines,file.path(root,"native-results.md"))
if(length(args)==2L) {
  cv <- readRDS(args[2])$results
  gl <- subset(cv,grepl("^GL_",METHOD))
  paired <- merge(subset(results,CRITERION%in%c("AIC","BIC")),
    gl[,c("ID","K","METHOD","ARI")],by="ID")
  stopifnot(nrow(paired)==24L)
  saveRDS(paired,file.path(root,"paired-native-vs-GL-CV.rds"))
  report <- c("# Native stockR versus GL predictive selection", "",
    "Six shared panels only; one realization per condition. Different selection criteria, not an isolated estimator comparison.",
    "GL CV results reuse the earlier frozen run. No GL retuning based on native results.",
    "", "| Dataset | Native protocol | Native criterion | stockR K | stockR test ARI | GL CV K | GL test ARI |",
    "| --- | --- | --- | ---: | ---: | ---: | ---: |")
  for(i in seq_len(nrow(paired))) {
    r <- paired[i,]
    report <- c(report,sprintf("| %s | %s | %s | %s | %.4f | %s | %.4f |",
      r$ID,r$PROTOCOL,r$CRITERION,r$SELECTED_K,r$TEST_ARI,r$K,r$ARI))
  }
  writeLines(report,file.path(root,"paired-native-vs-GL-CV.md"))
}
print(subset(results,MARKERS==10000 & CONDITION=="three_005"),row.names=FALSE)
print(with(diagnostics,table(CONVERGED)))
print(subset(rates,PROTOCOL=="native_best3"),row.names=FALSE)
