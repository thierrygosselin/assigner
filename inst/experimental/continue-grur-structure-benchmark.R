# Run from assigner source. Existing inputs/results are read-only.
# Rscript .../continue-grur-structure-benchmark.R INPUT NEW_OUTPUT [WORKERS]
args <- commandArgs(TRUE)
if(length(args)<2L) stop("Supply the original benchmark and a NEW output directory.")
input <- normalizePath(args[1],mustWork=TRUE)
output <- args[2]
workers <- if(length(args)>2L) as.integer(args[3]) else 2L
if(is.na(workers) || workers<1L || workers>4L) stop("Use 1 to 4 workers.")
if(file.exists(output)) stop("Output already exists; will not overwrite.")
sources <- c("R/assign_likelihood_dapc.R",file.path("R",paste0(c(
  "stock_likelihood_mixture","stock_mixture_blocks","stock_mixture_gds",
  "stock_mixture_validation"),".R")))
for(p in sources) source(p)
manifest <- readRDS(file.path(input,"manifest.rds"))
pilot <- readRDS(file.path(input,"pilot-results.rds"))
stopifnot(is.list(pilot),!is.data.frame(pilot),!anyDuplicated(manifest$ID))
remaining <- which(!manifest$ID %in% names(pilot))
if(!length(remaining)) stop("No remaining datasets.")
stopifnot(all(file.exists(manifest$GDS)))
dir.create(output,recursive=TRUE); output <- normalizePath(output)
saveRDS(list(session=sessionInfo(),source.md5=tools::md5sum(sources),
  input.md5=tools::md5sum(manifest$GDS),remaining=manifest$ID[remaining],
  workers=workers),file.path(output,"run-audit.rds"))

adjusted.rand <- function(a,b) {
  tab<-table(a,b); n<-sum(tab); choose2<-function(x) x*(x-1)/2
  if(n<2L) return(NA_real_)
  expected<-sum(choose2(rowSums(tab)))*sum(choose2(colSums(tab)))/choose2(n)
  maximum<-(sum(choose2(rowSums(tab)))+sum(choose2(colSums(tab))))/2
  if(maximum==expected) return(as.numeric(all(outer(a,a,"==")==outer(b,b,"=="))))
  (sum(choose2(tab))-expected)/(maximum-expected)
}

run.one <- function(j) {
  row<-manifest[j,]; key<-row$ID
  dest<-file.path(output,key); dir.create(dest)
  log<-file(file.path(dest,"progress.txt"),open="wt")
  sink(log); sink(log,type="message")
  on.exit({sink(type="message"); sink(); close(log)})
  started<-Sys.time()
  cat("Starting",key,format(started),"\n")
  result<-tryCatch({
    truth<-readRDS(file.path(input,"truth",paste0(key,".rds")))
    initial<-matrix(.5,row$MARKERS,row$GENERATING_K,
      dimnames=list(as.character(seq_len(row$MARKERS)),paste0("STOCK_",seq_len(row$GENERATING_K))))
    fit<-.sm_fit_gds(row$GDS,initial,likelihood="GL",start.strategy="mixed",
      n.starts=7,seed=36000+j,regularization="pooled",prior.strength=2,max.iter=300)
    saveRDS(fit,file.path(dest,"fixed-K-fit.rds"))
    cat("Known-K fit finished",format(Sys.time()),"\n")
    folds<-setNames(rep(1:2,length.out=row$SAMPLES),rownames(truth$genotype))
    cv<-.sm_cv_gds(row$GDS,K=1:5,folds=folds,n.starts=7,seed=46000+j,
      likelihood="GL",start.strategy="mixed",regularization="pooled",prior.strength=2,max.iter=300)
    saveRDS(cv,file.path(dest,"CV.rds"))
    keep<-!is.na(fit$assignments$ASSIGNED_STOCK)
    data.frame(ID=key,TRUE_K=row$GENERATING_K,TRUTH_WC84=row$TRUTH_WC84,
      CONVERGED=fit$converged,ASSIGNED=sum(keep),ARI=adjusted.rand(truth$group[keep],fit$assignments$ASSIGNED_STOCK[keep]),
      PREDICTIVE_K=cv$best.predictive.K,COMPLETE=all(cv$summary$COMPLETE),ERROR="")
  },error=function(e) data.frame(ID=key,TRUE_K=row$GENERATING_K,
    TRUTH_WC84=row$TRUTH_WC84,CONVERGED=FALSE,ASSIGNED=NA_integer_,ARI=NA_real_,
    PREDICTIVE_K=NA_integer_,COMPLETE=FALSE,ERROR=conditionMessage(e)))
  result$SECONDS<-as.numeric(difftime(Sys.time(),started,units="secs"))
  saveRDS(result,file.path(dest,"result.rds")); print(result)
  result
}

# Independent R processes open their own GDS handles; seeds depend on manifest
# rows, not worker scheduling. No nested parallelism or truth-informed starts.
results<-parallel::mclapply(remaining,run.one,mc.cores=workers,
  mc.preschedule=FALSE,mc.set.seed=FALSE)
if(any(vapply(results,inherits,logical(1),"try-error"))) stop("A worker failed; inspect its progress file.")
results<-do.call(rbind,results)
saveRDS(results,file.path(output,"remaining-results.rds"))
initial<-do.call(rbind,pilot); initial$ERROR<-""; initial$SECONDS<-NA_real_
combined<-rbind(initial,results)
combined<-combined[match(manifest$ID,combined$ID),]
stopifnot(nrow(combined)==nrow(manifest),!anyNA(combined$ID))
saveRDS(combined,file.path(output,"all-results.rds"))
writeLines(c("Same settings as pilot; 24 datasets, two replicates per scenario.",
  capture.output(print(combined,row.names=FALSE))),file.path(output,"all-results.txt"))
audit<-readRDS(file.path(output,"run-audit.rds"))
stopifnot(identical(audit$input.md5,tools::md5sum(manifest$GDS)))
cat("Completed",nrow(results),"additional datasets. Errors:",sum(nzchar(results$ERROR)),"\n")
