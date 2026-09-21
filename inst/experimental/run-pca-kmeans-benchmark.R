# Reuse frozen complete-genotype panels; do not generate new data or refit stockR.
# Rscript run-pca-kmeans-benchmark.R PILOT EXTENSION FIXED_F002 NEW_OUTPUT
.pk_scores <- function(g,maximum=100L) {
  stopifnot(is.matrix(g),!anyNA(g),all(g %in% 0:2),nrow(g)>2,ncol(g)>1)
  p <- colMeans(g)/2
  stopifnot(all(p>0 & p<1))
  z <- sweep(sweep(g,2,2*p,"-"),2,sqrt(2*p*(1-p)),"/")
  eig <- eigen(tcrossprod(z),symmetric=TRUE)
  rank <- sum(eig$values>max(eig$values)*1e-10)
  k <- min(maximum,rank,nrow(g)-1L)
  stopifnot(k>0)
  scores <- sweep(eig$vectors[,seq_len(k),drop=FALSE],2,sqrt(eig$values[seq_len(k)]),"*")
  rownames(scores) <- rownames(g)
  list(scores=scores,rank=rank,p=p,variance=eig$values[seq_len(k)]/sum(pmax(eig$values,0)))
}

.pk_kmeans <- function(scores,k,seed,batches=3L,starts=25L) {
  fits <- vector("list",batches); records <- vector("list",batches)
  for(b in seq_len(batches)) {
    set.seed(seed+b);warn <- character();before <- proc.time()[3]
    fit <- tryCatch(withCallingHandlers(stats::kmeans(scores,centers=k,nstart=starts,
      iter.max=100L,algorithm="Hartigan-Wong"),warning=function(w) {
        warn <<- c(warn,conditionMessage(w));invokeRestart("muffleWarning")
      }),error=function(e)e)
    error <- inherits(fit,"error")
    ok <- !error && is.finite(fit$tot.withinss) && (is.null(fit$ifault) || fit$ifault==0L)
    records[[b]] <- data.frame(BATCH=b,SEED=seed+b,OK=ok,
      WITHIN_SS=if(error) NA_real_ else fit$tot.withinss,
      ITERATIONS=if(error) NA_integer_ else fit$iter,
      IFAULT=if(error || is.null(fit$ifault)) NA_integer_ else fit$ifault,
      SECONDS=unname(proc.time()[3]-before),WARNING=paste(warn,collapse=" | "),
      ERROR=if(error) conditionMessage(fit) else "")
    fits[b] <- list(if(error) NULL else fit)
  }
  tab <- do.call(rbind,records);eligible <- which(tab$OK)
  selected <- if(length(eligible)) eligible[which.min(tab$WITHIN_SS[eligible])] else NA_integer_
  list(fits=fits,diagnostics=tab,selected=selected)
}

.pk_main <- function(args) {
  stopifnot(length(args)==4L,!file.exists(args[4]))
  inputs <- normalizePath(args[1:3],mustWork=TRUE)
  source("inst/experimental/diagnose-truth-starts.R")
  source("inst/experimental/benchmark-stockr-paired.R")
  audits <- lapply(inputs,function(p) readRDS(file.path(p,"design-locked.rds")))
  for(a in audits[-1]) stopifnot(identical(a$source.hash,audits[[1]]$source.hash),
    identical(a$stockR,audits[[1]]$stockR),identical(a$binary.hash,audits[[1]]$binary.hash))
  panels <- do.call(rbind,lapply(seq_along(inputs),function(i) {
    x <- audits[[i]]$design$panels;x$SOURCE <- inputs[i]
    x$REGIME <- if(i==3L) "F002" else "WEAK";x
  }))
  stopifnot(nrow(panels)==80L,!anyDuplicated(panels$ID),!anyDuplicated(panels$SEED))
  components <- c(2L,5L,10L,20L,50L,100L)
  dir.create(args[4],recursive=TRUE);out <- normalizePath(args[4])
  own.sources <- c("inst/experimental/run-pca-kmeans-benchmark.R",
    "inst/experimental/diagnose-truth-starts.R","inst/experimental/benchmark-stockr-paired.R")
  saveRDS(list(panels=panels,components=components,K=3L,batches=3L,starts.per.batch=25L,
    kmeans.iter.max=100L,scaling="Centered, divided by sqrt(2*p*(1-p)); PC scores are not whitened",
    initialization.note="Native stockR uses uncentered unscaled rotation, up to 100 components, 25 starts per call",
    selection="Minimum within-cluster sum of squares separately at every PC count; never choose PC count using truth",
    native.selection="Record exact initialization and final labels for every saved stockR run; highlight saved likelihood-selected run",
    session=sessionInfo(),input.audits=audits,source.hash=tools::md5sum(own.sources),
    source.text=lapply(own.sources,readLines)),file.path(out,"design-locked.rds"))
  run <- function(i) {
    row <- panels[i,];p <- file.path(row$SOURCE,row$ID);dest <- file.path(out,row$ID);dir.create(dest)
    cat("Starting",row$ID,"\n")
    truth <- readRDS(file.path(p,"truth.rds"));old <- readRDS(file.path(p,"result.rds"))
    g <- truth$genotype[,old$keep,drop=FALSE]
    stopifnot(nrow(g)==row$N,ncol(g)==old$result$L)
    projected <- .pk_scores(g,max(components))
    fits <- vector("list",length(components));results <- vector("list",length(components))
    for(j in seq_along(components)) {
      q <- components[j];stopifnot(q<=ncol(projected$scores))
      fit <- .pk_kmeans(projected$scores[,seq_len(q),drop=FALSE],3L,9500000L+i*1000L+j*10L)
      # Freeze selection using k-means objective before evaluating truth.
      saveRDS(fit,file.path(dest,paste0("pc",q,"-fit.rds")))
      valid <- !is.na(fit$selected)
      labels <- if(valid) fit$fits[[fit$selected]]$cluster else rep(NA_integer_,nrow(g))
      results[[j]] <- data.frame(ID=row$ID,REGIME=row$REGIME,N=row$N,L=ncol(g),L_NOMINAL=row$L,
        REPLICATE=row$REPLICATE,F_GENERATING=row$F,F_MEASURED=old$result$F_MEASURED,
        COMPONENTS=q,STATUS=if(valid) "OK" else "NO_VALID_RUN",SELECTED_BATCH=fit$selected,
        ARI=if(valid) .tb_ari(labels,truth$group) else NA_real_,
        WITHIN_SS=if(valid) fit$diagnostics$WITHIN_SS[fit$selected] else NA_real_,
        WARNINGS=sum(nzchar(fit$diagnostics$WARNING)),FAILED_BATCHES=sum(!fit$diagnostics$OK))
      fits[[j]] <- labels
    }
    native <- readRDS(file.path(p,"stock.rds"));selected <- readRDS(file.path(p,"choices-locked.rds"))$stock
    native.rows <- lapply(seq_along(native$fits),function(s) {
      fit <- native$fits[[s]]
      stopifnot(fit$control$method=="Kmeans.EM",fit$control$nKPCA==100,
        fit$control$nKmeanStart==25,identical(as.character(fit$data$sample.grps),as.character(seq_len(nrow(g)))),
        identical(names(fit$init.grps),as.character(seq_len(nrow(g)))),
        isTRUE(all.equal(unname(fit$data$SNPdata),unname(t(g)),check.attributes=FALSE)))
      pred <- .bp_stock_predict(fit,g);final <- max.col(pred$membership,ties.method="first")
      stopifnot(all(rowSums(pred$membership==apply(pred$membership,1,max))==1L))
      data.frame(ID=row$ID,REGIME=row$REGIME,N=row$N,L_NOMINAL=row$L,REPLICATE=row$REPLICATE,
        START=s,SELECTED=s==selected,INIT_ARI=.tb_ari(fit$init.grps,truth$group),
        FINAL_ARI=.tb_ari(final,truth$group),LOGLIK=fit$margLogL)
    })
    native.rows <- do.call(rbind,native.rows)
    stopifnot(abs(native.rows$FINAL_ARI[native.rows$SELECTED]-old$result$STOCK_ARI)<1e-10)
    result <- list(kmeans=do.call(rbind,results),native=native.rows,labels=fits,
      variance=projected$variance,frequency=projected$p,sample.ids=rownames(g),marker.ids=colnames(g),
      input.hash=tools::md5sum(file.path(p,c("truth.rds","stock.rds","result.rds","choices-locked.rds"))))
    saveRDS(result,file.path(dest,"result.rds"))
    cat("Completed",row$ID,"\n")
    result[c("kmeans","native")]
  }
  results <- parallel::mclapply(seq_len(nrow(panels)),run,mc.cores=2L,mc.preschedule=FALSE,mc.set.seed=FALSE)
  saveRDS(results,file.path(out,"results.rds"))
  if(any(vapply(results,inherits,logical(1),"try-error"))) stop("Failed panels; inspect output, do not drop them.")
  cat("Completed all",length(results),"panels.\n")
}
if(sys.nframe()==0L) .pk_main(commandArgs(TRUE))
