# Run from assigner: Rscript inst/experimental/run-stockr-native.R INPUT NEW_OUTPUT
# Native Kmeans.EM, no supplied starts; three independent seeds per K > 1.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
input <- normalizePath(args[1], mustWork = TRUE)
output <- args[2]
if (file.exists(output)) stop("Use a new output directory.")
dir.create(output, recursive = TRUE)
source("inst/experimental/benchmark-stockr-paired.R")
source("inst/experimental/diagnose-truth-starts.R")
manifest <- readRDS(file.path(input, "manifest.rds"))
saveRDS(list(session = sessionInfo(), manifest = manifest,
  sources = lapply(c("inst/experimental/run-stockr-native.R",
    "inst/experimental/benchmark-stockr-paired.R",
    "inst/experimental/diagnose-truth-starts.R"), readLines),
  getIC = stockR:::getIC, defaults = stockR:::set.control,
  design = "K1:4; native defaults except quiet; three seeds; first and best-converged reported separately"),
  file.path(output, "audit.rds"))
run <- function(j) {
  row <- manifest[j, ]; dest <- file.path(output, row$ID); dir.create(dest)
  d <- readRDS(file.path(input, paste(row$CONDITION,row$REPLICATE,sep="_"),"truth.rds"))
  g <- d$train[,seq_len(row$MARKERS),drop=FALSE]
  keep <- which(apply(g,2,function(x) length(unique(x))>1))
  g <- g[,keep,drop=FALSE]
  records <- list(); fits <- list()
  for(k in 1:4) for(s in seq_len(if(k==1) 1 else 3)) {
    seed <- 810000+j*100+k*10+s; set.seed(seed)
    warnings <- character(); before <- proc.time()[3]
    fit <- tryCatch(withCallingHandlers(suppressMessages(stockR::stockSTRUCTURE(
      t(g),K=k,sample.grps=factor(seq_len(nrow(g)),levels=seq_len(nrow(g))),
      start.grps=NULL,control=list(quiet=TRUE))),warning=function(w) {
        warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
      }),error=function(e)e)
    elapsed <- unname(proc.time()[3]-before)
    ok <- !inherits(fit,"error")
    it <- if(!ok) NA_integer_ else if(k==1) 0L else {
      z <- length(fit$margLogLTrace)
      if(z<=fit$control$EM.maxit) z else z-1L
    }
    converged <- ok && (k==1 || it<=fit$control$EM.maxit)
    if(ok) {
      # Check returned native IC values; never substitute another sample size.
      np <- k-1+nrow(fit$condMeans)*k
      stopifnot(isTRUE(all.equal(fit$AIC,-2*fit$margLogL+2*np)),
        isTRUE(all.equal(fit$BIC,-2*fit$margLogL+log(nrow(g))*np)))
    }
    key <- paste(k,s,sep="_"); fits[[key]] <- fit
    records[[key]] <- data.frame(row,K=k,START=s,SEED=seed,CONVERGED=converged,
      ITERATIONS=it,SECONDS=elapsed,RETAINED=length(keep),
      LOGLIK=if(ok) fit$margLogL else NA_real_,AIC=if(ok) fit$AIC else NA_real_,
      BIC=if(ok) fit$BIC else NA_real_,WARNING=paste(warnings,collapse=" | "),
      ERROR=if(ok) "" else conditionMessage(fit))
    saveRDS(records,file.path(dest,"progress.rds"))
  }
  tab <- do.call(rbind,records)
  # Lock choices from training IC before computing any test diagnostics.
  choices <- list()
  for(protocol in c("native_single","native_best3")) {
    ix <- which(tab$CONVERGED & (protocol=="native_best3" | tab$START==1))
    best <- lapply(1:4,function(k) {
      candidates <- ix[tab$K[ix]==k]
      if(!length(candidates)) NA_integer_ else candidates[which.max(tab$LOGLIK[candidates])]
    }) |> unlist()
    for(criterion in c("AIC","BIC","known_K")) {
      chosen <- if(anyNA(best)) NA_integer_ else if(criterion=="known_K")
        best[row$TRUE_K] else best[which.min(tab[[criterion]][best])]
      choices[[paste(protocol,criterion,sep="_")]] <- list(protocol=protocol,
        criterion=criterion,index=chosen)
    }
  }
  saveRDS(choices,file.path(dest,"choices-locked.rds"))
  result <- lapply(choices,function(ch) {
    ix <- ch$index
    pred <- if(is.na(ix)) NULL else .bp_stock_predict(fits[[rownames(tab)[ix]]],
      d$test[,seq_len(row$MARKERS),drop=FALSE][,keep,drop=FALSE])
    data.frame(row,PROTOCOL=ch$protocol,CRITERION=ch$criterion,
      SELECTED_K=if(is.na(ix)) NA else tab$K[ix],
      TEST_ARI=if(is.null(pred)) NA else .tb_partition(pred$membership,d$group)$ari,
      TEST_LOGLIK=if(is.null(pred)) NA else mean(pred$sample.loglik))
  }) |> do.call(what=rbind)
  saveRDS(list(fits=fits,diagnostics=tab,choices=choices,results=result),file.path(dest,"result.rds"))
  cat("Completed",row$ID,"\n")
  result
}
results <- parallel::mclapply(seq_len(nrow(manifest)),run,mc.cores=2L,mc.preschedule=FALSE)
saveRDS(results,file.path(output,"results.rds"))
if(any(vapply(results,inherits,logical(1),"try-error"))) stop("Some panels failed; inspect results.")
print(do.call(rbind,results),row.names=FALSE)
