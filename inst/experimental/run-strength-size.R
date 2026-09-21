# Frozen factorial experiment: sample size x strength; initialization unchanged.
# Run from assigner: Rscript inst/experimental/run-strength-size.R OUTPUT [WORKERS]
args <- commandArgs(TRUE); stopifnot(length(args)%in%1:2)
output <- args[1]; workers <- if(length(args)==2) as.integer(args[2]) else 2L
stopifnot(workers%in%1:4)
paths <- c("R/stock_likelihood_mixture.R","R/stock_mixture_blocks.R",
  "inst/experimental/benchmark-streaming-engine.R","inst/experimental/benchmark-stockr-paired.R",
  "inst/experimental/benchmark-hard-limit.R","inst/experimental/diagnose-truth-starts.R",
  "../grur/R/simulation_benchmark.R")
for(p in paths) source(p)
hash <- tools::md5sum(c(paths,"inst/experimental/run-strength-size.R"))
design <- expand.grid(REALIZATION=1:5,N=c(300L,600L,900L))
design$CONDITION <- ifelse(design$REALIZATION<=3,"intermediate","null")
design$TRUE_K <- ifelse(design$REALIZATION<=3,3L,1L)
design$REPLICATE <- ifelse(design$REALIZATION<=3,design$REALIZATION,design$REALIZATION-3)
design$ID <- with(design,paste(CONDITION,REPLICATE,N,sep="_"))
strengths <- c(20,50,100,200,500)
dir.create(output,recursive=TRUE,showWarnings=FALSE); output <- normalizePath(output)
audit.path <- file.path(output,"audit.rds")
if(file.exists(audit.path)) {
  audit <- readRDS(audit.path)
  stopifnot(identical(audit$hash,hash),identical(audit$stockR,as.character(packageVersion("stockR"))))
} else saveRDS(list(design=design,strengths=strengths,K=1:4,folds=5,markers=3000,
  train.N=c(300,600,900),test.N=300,starts=7,hash=hash,stockR=as.character(packageVersion("stockR")),
  session=sessionInfo(),sources=setNames(lapply(c(paths,"inst/experimental/run-strength-size.R"),readLines),
    c(paths,"inst/experimental/run-strength-size.R"))),audit.path)
for(r in 1:5) {
  dest <- file.path(output,paste0("data",r)); dir.create(dest,showWarnings=FALSE)
  if(file.exists(file.path(dest,"ready.rds"))) next
  if(any(file.exists(file.path(dest,c("train.gds","test.gds"))))) stop("Partial input generation; inspect files.")
  k <- if(r<=3) 3L else 1L; seed <- 2480000+r
  truth <- .grur_structure_truth(rep(900/k,k),3000,if(k==3) .005 else 0,seed)
  test.groups <- rep(seq_len(k),300/k)
  test <- .grur_benchmark_seed(seed+10000,matrix(rbinom(300*3000,2,
    as.vector(t(truth$frequency[,test.groups,drop=FALSE]))),300,3000))
  dimnames(test) <- list(paste0("T",1:300),colnames(truth$genotype))
  saveRDS(list(group=truth$group,test.group=test.groups,frequency=truth$frequency,seed=seed),file.path(dest,"truth.rds"))
  .bp_write_calls(truth$genotype,file.path(dest,"train.gds")); .bp_write_calls(test,file.path(dest,"test.gds"))
  saveRDS(tools::md5sum(file.path(dest,c("train.gds","test.gds"))),file.path(dest,"ready.rds"))
}
choose <- function(grid) {
  stopifnot(all(is.finite(grid$SCORE)))
  x <- grid[grid$SCORE>=max(grid$SCORE)-1e-8,,drop=FALSE]
  x[order(x$K,x$STRENGTH),,drop=FALSE][1,,drop=FALSE]
}
run <- function(j) {
  row <- design[j,]; dest <- file.path(output,row$ID); dir.create(dest,showWarnings=FALSE)
  result.path <- file.path(dest,"result.rds")
  if(file.exists(result.path)) return(readRDS(result.path))
  log <- file(file.path(dest,"progress.txt"),"wt"); sink(log);sink(log,type="message")
  on.exit({sink(type="message");sink();close(log)})
  result <- tryCatch({
    src <- file.path(output,paste0("data",row$REALIZATION))
    # Nested IDs are determined by simulation design; labels never enter fitting.
    ids <- if(row$TRUE_K==3) unlist(lapply(0:2,function(k) k*300+seq_len(row$N/3))) else seq_len(row$N)
    all.calls <- .hl_read(file.path(src,"train.gds")); g <- all.calls[ids,,drop=FALSE]; rm(all.calls)
    folds <- rep(1:5,length.out=nrow(g)); saveRDS(list(ids=rownames(g),folds=folds),file.path(dest,"samples-folds.rds"))
    cv <- list(); diagnostics <- list(); initialization.warnings <- list()
    for(f in 1:5) {
      train <- which(folds!=f); held <- which(folds==f)
      keep <- which(apply(g[train,,drop=FALSE],2,function(x)length(unique(x))>1))
      gt <- g[train,keep,drop=FALSE]
      for(k in 1:4) {
        warning <- character()
        starts <- withCallingHandlers(.bp_starts(gt,k,2490000+j*100+f*10+k),warning=function(w) {
          warning <<- c(warning,conditionMessage(w));invokeRestart("muffleWarning")})
        initialization.warnings[[paste(f,k,sep="_")]] <- warning
        for(s in strengths) {
          model <- .hl_models(gt,k,starts,s); key <- paste(f,k,s,sep="_")
          diagnostics[[key]] <- model$diagnostics
          pred <- if(is.null(model$fit)) NULL else .hl_predict(g[held,keep,drop=FALSE],model$fit$frequency,model$fit$proportion)
          cv[[key]] <- data.frame(FOLD=f,K=k,STRENGTH=s,SCORE=if(is.null(pred)) NA_real_ else mean(pred$sample.loglik))
        }
        saveRDS(list(cv=cv,diagnostics=diagnostics,warnings=initialization.warnings),file.path(dest,"CV-checkpoint.rds"))
        cat("Fold",f,"K",k,"done\n");flush(log)
      }
    }
    cv <- do.call(rbind,cv); if(any(!is.finite(cv$SCORE))) stop("Incomplete strength/K grid.")
    grid <- aggregate(SCORE~K+STRENGTH,cv,mean)
    choices <- rbind(data.frame(WORKFLOW="fixed100",choose(subset(grid,STRENGTH==100))),
      data.frame(WORKFLOW="tuned",choose(grid)))
    # Native comparator at the same sample sizes, separately from strength effects.
    keep <- which(apply(g,2,function(x)length(unique(x))>1)); gt <- g[,keep,drop=FALSE]
    native <- list(); nd <- list()
    for(k in 1:4) for(s in seq_len(if(k==1) 1 else 3)) {
      set.seed(2500000+j*100+k*10+s); warnings <- character()
      fit <- tryCatch(withCallingHandlers(suppressMessages(stockR::stockSTRUCTURE(t(gt),K=k,
        sample.grps=factor(1:nrow(gt),levels=1:nrow(gt)),control=list(quiet=TRUE))),warning=function(w) {
          warnings <<- c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
      ok <- !inherits(fit,"error"); it <- if(!ok) NA else if(k==1) 0 else {
        n <- length(fit$margLogLTrace);if(n<=fit$control$EM.maxit) n else n-1L}
      key <- paste(k,s,sep="_");native[[key]] <- fit
      nd[[key]] <- data.frame(K=k,START=s,VALID=ok&&(k==1||it<=fit$control$EM.maxit),ITERATIONS=it,
        WARNING=paste(warnings,collapse=" | "),ERROR=if(ok) "" else conditionMessage(fit),
        LOGLIK=if(ok) fit$margLogL else NA,AIC=if(ok) fit$AIC else NA,BIC=if(ok) fit$BIC else NA)
    }
    nd <- do.call(rbind,nd)
    best <- vapply(1:4,function(k) {ix<-which(nd$K==k&nd$VALID);if(!length(ix)) NA_integer_ else ix[which.max(nd$LOGLIK[ix])]},integer(1))
    if(anyNA(best)) stop("Incomplete native K grid.")
    nc <- c(AIC=best[which.min(nd$AIC[best])],BIC=best[which.min(nd$BIC[best])])
    saveRDS(list(GL=choices,native=nc,grid=grid),file.path(dest,"choices-before-test.rds"))
    # Test GDS and truth are accessed only after choices are locked.
    test <- .hl_read(file.path(src,"test.gds")); truth <- readRDS(file.path(src,"truth.rds"))
    rows <- list(); final <- list()
    for(q in 1:nrow(choices)) {
      ch <- choices[q,]; starts <- .bp_starts(gt,ch$K,2510000+j*10+ch$K)
      model <- .hl_models(gt,ch$K,starts,ch$STRENGTH);if(is.null(model$fit)) stop("Final fit failed.")
      # Verify final likelihoods against the actual streaming GDS likelihood engine.
      ctx <- .bs_context(file.path(src,"train.gds"),sample.ids=rownames(gt),variant.ids=keep)
      stopifnot(identical(ctx$ids,rownames(gt)),
        identical(names(model$fit$sample.loglik),rownames(gt)))
      state <- tryCatch(.sm_block_estep(ctx$reader,ctx$blocks,model$fit$frequency,model$fit$proportion,ctx$n),finally=ctx$close())
      stopifnot(isTRUE(all.equal(state$sample.loglik,unname(model$fit$sample.loglik),tolerance=1e-8)))
      pred <- .hl_predict(test[,keep,drop=FALSE],model$fit$frequency,model$fit$proportion)
      final[[ch$WORKFLOW]] <- list(model=model,prediction=pred)
      rows[[ch$WORKFLOW]] <- data.frame(row,METHOD=ch$WORKFLOW,K=ch$K,STRENGTH=ch$STRENGTH,
        TEST_ARI=.tb_partition(pred$membership,truth$test.group)$ari,TEST_LOGLIK=mean(pred$sample.loglik))
    }
    for(criterion in names(nc)) {
      fit <- native[[rownames(nd)[nc[[criterion]]]]];pred <- .bp_stock_predict(fit,test[,keep,drop=FALSE])
      rows[[criterion]] <- data.frame(row,METHOD=paste0("stockR_",criterion),K=fit$K,STRENGTH=NA_real_,
        TEST_ARI=.tb_partition(pred$membership,truth$test.group)$ari,TEST_LOGLIK=mean(pred$sample.loglik))
    }
    stopifnot(identical(readRDS(file.path(src,"ready.rds")),tools::md5sum(file.path(src,c("train.gds","test.gds")))))
    saveRDS(list(cv=cv,grid=grid,diagnostics=diagnostics,initialization.warnings=initialization.warnings,
      native=native,native.diagnostics=nd,final=final),file.path(dest,"fits.rds"))
    do.call(rbind,rows)
  },error=function(e) list(ID=row$ID,error=conditionMessage(e)))
  saveRDS(result,result.path);result
}
results <- parallel::mclapply(seq_len(nrow(design)),run,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
saveRDS(results,file.path(output,"results.rds"))
if(any(!vapply(results,is.data.frame,logical(1)))) {print(results);stop("Inspect saved failures.")}
print(do.call(rbind,results),row.names=FALSE)
