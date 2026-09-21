# Post-hoc development diagnostic on previously missed K, NOT fresh validation.
# Run from assigner with PROOF_ROOT NEW_OUTPUT.
args <- commandArgs(TRUE);stopifnot(length(args)==2L,!file.exists(args[2]))
root <- normalizePath(args[1],mustWork=TRUE); output <- args[2];dir.create(output,recursive=TRUE)
paths <- c("R/stock_mixture_blocks.R","inst/experimental/benchmark-streaming-engine.R",
  "inst/experimental/benchmark-stockr-paired.R","inst/experimental/benchmark-hard-limit.R",
  "inst/experimental/diagnose-truth-starts.R","../grur/R/simulation_benchmark.R")
for(p in paths) source(p)
saveRDS(list(source=tools::md5sum(c(paths,"inst/experimental/revisit-strength-failures.R")),
  session=sessionInfo(),strengths=c(20,50,100,200,500),
  note="Known failures selected after inspection; diagnostic only; original seeds/folds retained"),file.path(output,"audit.rds"))
src <- file.path(root,"intermediate_2"); path <- file.path(src,"train.gds")
calls <- .hl_read(path); all.results <- list()
for(m in c(3000,10000)) {
  j <- if(m==3000) 5L else 11L; id <- paste0("intermediate_2_",m)
  dest <- file.path(output,id);dir.create(dest)
  g <- calls[,seq_len(m),drop=FALSE]; folds <- readRDS(file.path(root,id,"folds.rds"))
  cv <- list();diagnostics <- list()
  for(f in 1:5) for(k in 1:4) {
    train <- which(folds!=f);held <- which(folds==f)
    keep <- which(apply(g[train,,drop=FALSE],2,function(x)length(unique(x))>1))
    gt <- g[train,keep,drop=FALSE]
    starts <- .bp_starts(gt,k,1930000+j*100+f*10+k)
    for(s in c(20,50,100,200,500)) {
      model <- .hl_models(gt,k,starts,s);if(is.null(model$fit))stop("Incomplete diagnostic grid.")
      pred <- .hl_predict(g[held,keep,drop=FALSE],model$fit$frequency,model$fit$proportion)
      key <- paste(f,k,s,sep="_")
      cv[[key]] <- data.frame(FOLD=f,K=k,STRENGTH=s,SCORE=mean(pred$sample.loglik))
      diagnostics[[key]] <- model$diagnostics
    }
    saveRDS(list(cv=cv,diagnostics=diagnostics),file.path(dest,"checkpoint.rds"))
  }
  cv <- do.call(rbind,cv);grid <- aggregate(SCORE~K+STRENGTH,cv,mean)
  old <- readRDS(file.path(root,id,"fits.rds"))$cv
  comparison <- merge(subset(cv,STRENGTH==100),old,by=c("FOLD","K"))
  discrepancy <- max(abs(comparison$SCORE.x-comparison$SCORE.y))
  stopifnot(nrow(comparison)==20L,discrepancy<1e-5)
  eligible <- grid[grid$SCORE>=max(grid$SCORE)-1e-8,]
  chosen <- eligible[order(eligible$K,eligible$STRENGTH),][1,]
  saveRDS(list(choice=chosen,grid=grid,baseline.max.difference=discrepancy),file.path(dest,"choice-before-test.rds"))
  keep <- which(apply(g,2,function(x)length(unique(x))>1));gt <- g[,keep,drop=FALSE]
  model <- .hl_models(gt,chosen$K,.bp_starts(gt,chosen$K,1950000+j*10+chosen$K),chosen$STRENGTH)
  if(is.null(model$fit))stop("Final diagnostic fit failed.")
  d <- readRDS(file.path(src,"truth.rds"))
  pred <- .hl_predict(d$test[,seq_len(m),drop=FALSE][,keep,drop=FALSE],model$fit$frequency,model$fit$proportion)
  result <- data.frame(ID=id,chosen,TEST_ARI=.tb_partition(pred$membership,d$group)$ari,
    TEST_LOGLIK=mean(pred$sample.loglik),BASELINE_SCORE_DIFFERENCE=discrepancy)
  saveRDS(list(result=result,grid=grid,cv=cv,model=model,prediction=pred,diagnostics=diagnostics),file.path(dest,"result.rds"))
  all.results[[id]] <- result;print(result)
}
saveRDS(do.call(rbind,all.results),file.path(output,"results.rds"))
