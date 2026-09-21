# Fresh-replicate proof of concept. Run from assigner with NEW_OUTPUT [workers].
args <- commandArgs(TRUE)
stopifnot(length(args)%in%1:2)
output <- args[1]; workers <- if(length(args)==2) as.integer(args[2]) else 2L
stopifnot(workers%in%1:4)
if(file.exists(output)) stop("Use a new output directory.")
paths <- c("R/stock_likelihood_mixture.R","R/stock_mixture_blocks.R",
  "inst/experimental/benchmark-streaming-engine.R",
  "inst/experimental/benchmark-stockr-paired.R",
  "inst/experimental/diagnose-truth-starts.R","../grur/R/simulation_benchmark.R")
for(p in paths) source(p)
dir.create(output,recursive=TRUE); output <- normalizePath(output)
design <- rbind(data.frame(CONDITION="intermediate",TRUE_K=3L,DRIFT=.005,
  expand.grid(REPLICATE=1:3,MARKERS=c(1000L,3000L,5000L,10000L))),
  data.frame(CONDITION=rep(c("weak","null"),each=2),TRUE_K=rep(c(3L,1L),each=2),
    DRIFT=rep(c(.001,0),each=2),REPLICATE=rep(1:2,2),MARKERS=10000L))
design$ID <- with(design,paste(CONDITION,REPLICATE,MARKERS,sep="_"))
saveRDS(list(design=design,strength=100,K=1:4,folds=5,starts=7,
  train.n=300,test.n=300,session=sessionInfo(),sources=setNames(lapply(c(paths,
    "inst/experimental/run-stock-proof.R"),readLines),c(paths,"inst/experimental/run-stock-proof.R")),
  source.hash=tools::md5sum(c(paths,"inst/experimental/run-stock-proof.R")),
  note="Fresh frequencies; fixed pilot strength100; no tuning; nested panels not independent replicates"),
  file.path(output,"design-and-audit.rds"))
# Generate each realization once. Test individuals share frequencies, not genotypes.
for(condition in unique(design$CONDITION)) for(r in unique(design$REPLICATE[design$CONDITION==condition])) {
  row <- subset(design,CONDITION==condition & REPLICATE==r)[1,]
  seed <- 1920000+match(condition,c("intermediate","weak","null"))*100+r
  truth <- .grur_structure_truth(rep(300/row$TRUE_K,row$TRUE_K),10000,row$DRIFT,seed)
  test <- .grur_benchmark_seed(seed+10000,matrix(rbinom(300*10000,2,
    as.vector(t(truth$frequency[,truth$group,drop=FALSE]))),300,10000))
  dimnames(test) <- list(paste0("T",1:300),colnames(truth$genotype))
  dest <- file.path(output,paste(condition,r,sep="_")); dir.create(dest)
  saveRDS(list(train=truth$genotype,test=test,group=truth$group,seed=seed,
    frequency=truth$frequency),file.path(dest,"truth.rds"))
  .bp_write_calls(truth$genotype,file.path(dest,"train.gds"))
  saveRDS(tools::md5sum(file.path(dest,"train.gds")),file.path(dest,"checksum.rds"))
}
run <- function(j) {
  row <- design[j,]; dest <- file.path(output,row$ID); dir.create(dest)
  log <- file(file.path(dest,"progress.txt"),"wt"); sink(log); sink(log,type="message")
  on.exit({sink(type="message");sink();close(log)})
  result <- tryCatch({
    src <- file.path(output,paste(row$CONDITION,row$REPLICATE,sep="_"))
    d <- readRDS(file.path(src,"truth.rds")); g <- d$train[,seq_len(row$MARKERS),drop=FALSE]
    folds <- rep(1:5,length.out=nrow(g))
    # Labels are not used to make folds, starts or selections.
    saveRDS(folds,file.path(dest,"folds.rds"))
    glfit <- function(train,k,seed) {
      keep <- which(apply(g[train,,drop=FALSE],2,function(x)length(unique(x))>1))
      gt <- g[train,keep,drop=FALSE]
      context <- .bs_context(file.path(src,"train.gds"),sample.ids=rownames(gt),variant.ids=keep)
      on.exit(context$close())
      stopifnot(identical(context$ids,rownames(gt)),identical(context$variants,as.character(keep)))
      starts <- .bp_starts(gt,k,seed)
      list(model=.bp_gl(context,gt,k,starts,100),keep=keep)
    }
    glpredict <- function(model,calls) {
      # Exact near-certain-GL likelihood for these hard-call-limit simulations.
      fit <- model$model$fit; if(is.null(fit)) return(NULL)
      calls <- calls[,model$keep,drop=FALSE]; joint <- matrix(0,nrow(calls),ncol(fit$frequency))
      for(k in seq_len(ncol(joint))) for(index in split(seq_len(ncol(calls)),ceiling(seq_len(ncol(calls))/250))) {
        p <- fit$frequency[index,k]; z <- matrix(dbinom(calls[,index,drop=FALSE],2,
          rep(p,each=nrow(calls)),log=TRUE),nrow(calls))
        joint[,k] <- joint[,k]+rowSums(z)
      }
      joint <- sweep(joint,2,log(fit$proportion),"+"); peak <- apply(joint,1,max)
      score <- peak+log(rowSums(exp(joint-peak)))
      list(membership=exp(joint-score),sample.loglik=score)
    }
    cv <- list(); diagnostics <- list()
    for(f in 1:5) for(k in 1:4) {
      model <- glfit(which(folds!=f),k,1930000+j*100+f*10+k)
      key <- paste(f,k,sep="_"); diagnostics[[key]] <- model$model$starts
      pred <- glpredict(model,g[folds==f,,drop=FALSE])
      cv[[key]] <- data.frame(FOLD=f,K=k,SCORE=if(is.null(pred)) NA_real_ else mean(pred$sample.loglik))
      saveRDS(list(cv=cv,diagnostics=diagnostics),file.path(dest,"CV-checkpoint.rds"))
      cat("Fold",f,"K",k,"done\n");flush(log)
    }
    cv <- do.call(rbind,cv); if(any(!is.finite(cv$SCORE))) stop("Incomplete GL CV grid.")
    scores <- aggregate(SCORE~K,cv,mean)
    chosen <- min(scores$K[scores$SCORE>=max(scores$SCORE)-1e-8])
    keep <- which(apply(g,2,function(x)length(unique(x))>1))
    native <- list(); native.diag <- list()
    for(k in 1:4) for(s in seq_len(if(k==1) 1 else 3)) {
      set.seed(1940000+j*100+k*10+s); warnings <- character(); before <- proc.time()[3]
      fit <- tryCatch(withCallingHandlers(suppressMessages(stockR::stockSTRUCTURE(t(g[,keep,drop=FALSE]),K=k,
        sample.grps=factor(1:nrow(g),levels=1:nrow(g)),control=list(quiet=TRUE))),warning=function(w) {
          warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
      ok <- !inherits(fit,"error"); it <- if(!ok) NA else if(k==1) 0 else {
        n <- length(fit$margLogLTrace); if(n<=fit$control$EM.maxit) n else n-1L
      }
      valid <- ok && (k==1 || it<=fit$control$EM.maxit)
      key <- paste(k,s,sep="_"); native[[key]] <- fit
      native.diag[[key]] <- data.frame(K=k,START=s,VALID=valid,ITERATIONS=it,
        SECONDS=unname(proc.time()[3]-before),WARNING=paste(warnings,collapse=" | "),
        ERROR=if(ok) "" else conditionMessage(fit),LOGLIK=if(ok) fit$margLogL else NA,
        AIC=if(ok) fit$AIC else NA,BIC=if(ok) fit$BIC else NA)
    }
    nt <- do.call(rbind,native.diag)
    best <- vapply(1:4,function(k) {
      ix <- which(nt$K==k & nt$VALID)
      if(!length(ix)) NA_integer_ else ix[which.max(nt$LOGLIK[ix])]
    },integer(1))
    if(anyNA(best)) stop("Incomplete native K grid.")
    native.choices <- c(AIC=best[which.min(nt$AIC[best])],BIC=best[which.min(nt$BIC[best])])
    saveRDS(list(GL_K=chosen,native=native.choices,scores=scores),file.path(dest,"choices-before-test.rds"))
    full <- glfit(1:nrow(g),chosen,1950000+j*10+chosen)
    if(is.null(full$model$fit)) stop("Final GL fit failed.")
    check <- glpredict(full,g)
    stopifnot(isTRUE(all.equal(check$sample.loglik,full$model$fit$sample.loglik,
      tolerance=1e-8,check.attributes=FALSE)))
    test <- d$test[,seq_len(row$MARKERS),drop=FALSE]
    prediction <- glpredict(full,test); if(is.null(prediction)) stop("Final GL fit failed.")
    results <- list(data.frame(row,METHOD="GL_fixed100_CV",SELECTED_K=chosen,
      TEST_ARI=.tb_partition(prediction$membership,d$group)$ari,
      TEST_LOGLIK=mean(prediction$sample.loglik)))
    for(criterion in names(native.choices)) {
      ix <- native.choices[[criterion]]; fit <- native[[rownames(nt)[ix]]]
      pred <- .bp_stock_predict(fit,test[,keep,drop=FALSE])
      results[[criterion]] <- data.frame(row,METHOD=paste0("stockR_",criterion),SELECTED_K=fit$K,
        TEST_ARI=.tb_partition(pred$membership,d$group)$ari,TEST_LOGLIK=mean(pred$sample.loglik))
    }
    stopifnot(identical(readRDS(file.path(src,"checksum.rds")),tools::md5sum(file.path(src,"train.gds"))))
    saveRDS(list(native=native,native.diagnostics=nt,GL=full,GL.prediction=prediction,
      cv=cv,cv.diagnostics=diagnostics,scores=scores),file.path(dest,"fits.rds"))
    do.call(rbind,results)
  },error=function(e) list(ID=row$ID,error=conditionMessage(e)))
  saveRDS(result,file.path(dest,"result.rds")); result
}
results <- parallel::mclapply(seq_len(nrow(design)),run,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
saveRDS(results,file.path(output,"results.rds"))
bad <- !vapply(results,is.data.frame,logical(1))
if(any(bad)) {print(results[bad]);stop("Failures saved.")}
print(do.call(rbind,results),row.names=FALSE)
