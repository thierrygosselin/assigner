# Privileged simulation diagnostic; never a deployable initialization method.
# Run from assigner: Rscript .../diagnose-truth-starts.R INPUT EXTENDED NEW_OUTPUT [WORKERS]
.tb_ari <- function(a,b) {
  keep<-!is.na(a)&!is.na(b); a<-a[keep]; b<-b[keep]
  tab<-table(a,b); n<-sum(tab); c2<-function(x) x*(x-1)/2
  if(n<2L) return(NA_real_)
  expected<-sum(c2(rowSums(tab)))*sum(c2(colSums(tab)))/c2(n)
  maximum<-(sum(c2(rowSums(tab)))+sum(c2(colSums(tab))))/2
  if(maximum==expected) return(as.numeric(all(outer(a,a,"==")==outer(b,b,"=="))))
  (sum(c2(tab))-expected)/(maximum-expected)
}
.tb_partition <- function(membership,group) {
  peak<-apply(membership,1,max)
  labels<-max.col(membership,ties.method="first")
  labels[rowSums(membership==peak)!=1L]<-NA_integer_
  list(ari=.tb_ari(group,labels),assigned=sum(!is.na(labels)),
    mean.max.membership=mean(peak))
}

.tb_main <- function(args) {
  if(length(args)<3L) stop("Supply benchmark input, continuation results and NEW output.")
  input<-normalizePath(args[1],mustWork=TRUE)
  extended<-normalizePath(args[2],mustWork=TRUE)
  output<-args[3]; workers<-if(length(args)>3L) as.integer(args[4]) else 2L
  if(is.na(workers)||workers<1L||workers>4L) stop("Use 1 to 4 workers.")
  if(file.exists(output)) stop("Refusing to overwrite diagnostic output.")
  sources<-c("R/assign_likelihood_dapc.R",file.path("R",paste0(c(
    "stock_likelihood_mixture","stock_mixture_blocks","stock_mixture_gds",
    "stock_mixture_validation"),".R")))
  for(p in sources) source(p)
  previous.audit<-readRDS(file.path(extended,"run-audit.rds"))
  stopifnot(identical(previous.audit$source.md5,tools::md5sum(sources)))
  manifest<-readRDS(file.path(input,"manifest.rds"))
  pilot<-readRDS(file.path(input,"pilot-results.rds"))
  selected<-which(manifest$SCENARIO %in% c("balanced_F0.001","balanced_F0.01",
    "balanced_F0.1","more_markers","large_sample"))
  dir.create(output,recursive=TRUE); output<-normalizePath(output)
  input.hash<-tools::md5sum(manifest$GDS[selected])
  saveRDS(list(session=sessionInfo(),source.md5=tools::md5sum(sources),
    input.md5=input.hash,manifest=manifest[selected,],workers=workers,
    interpretation="True frequencies are privileged information, not an ordinary CV procedure"),
    file.path(output,"audit.rds"))
  original.path<-function(id,suffix) {
    if(id %in% names(pilot)) file.path(input,paste0(id,"_",suffix,".rds")) else
      file.path(extended,id,paste0(suffix,".rds"))
  }
  run.one<-function(j) {
    row<-manifest[j,]; id<-row$ID; dest<-file.path(output,id); dir.create(dest)
    log<-file(file.path(dest,"progress.txt"),"wt"); sink(log); sink(log,type="message")
    on.exit({sink(type="message"); sink(); close(log)})
    cat("Starting",id,format(Sys.time()),"\n")
    outcome<-tryCatch({
      truth<-readRDS(file.path(input,"truth",paste0(id,".rds")))
      blind<-readRDS(original.path(id,"fixed-K-fit"))
      blind.cv<-readRDS(original.path(id,"CV"))
      stopifnot(identical(rownames(truth$genotype),blind$sample.ids),
        identical(rownames(truth$frequency),as.character(blind$variant.ids)),
        blind$regularization=="pooled",blind$prior.strength==2)
      n<-row$SAMPLES; m<-row$MARKERS
      p<-pmin(pmax(truth$frequency,blind$epsilon),1-blind$epsilon)
      pi<-truth$sample.size/sum(truth$sample.size)
      # Read GL only; no truth genotypes or population labels enter the likelihood.
      g<-SeqArray::seqOpen(row$GDS,readonly=TRUE)
      raw<-tryCatch(SeqArray::seqGetData(g,"annotation/format/GL"),
        finally=SeqArray::seqClose(g))
      stopifnot(all(raw$length==3L),identical(dim(raw$data),c(as.integer(n),as.integer(3*m))))
      prepared<-.sm_prepare(array(as.numeric(raw$data),c(n,3,m)),"GL")
      reader<-function(index) list(log.likelihood=prepared$log.likelihood[,,index,drop=FALSE],
        informative=prepared$informative[,index,drop=FALSE])
      frozen<-.sm_block_estep(reader,list(seq_len(m)),p,pi,n)
      frozen.summary<-.tb_partition(frozen$membership,truth$group)
      saveRDS(frozen,file.path(dest,"frozen-true-model.rds"))
      # Same fitting objective, bounds and stopping criteria as the saved blind fit.
      # n.starts=1 deliberately tests this supplied start, without another search.
      seeded<-.sm_audit_fit(.sm_fit_gds(row$GDS,p,proportion=pi,likelihood="GL",
        n.starts=1,seed=56000+j,regularization="pooled",prior.strength=2,
        max.iter=300,tolerance=1e-8,epsilon=blind$epsilon))
      if(inherits(seeded$fit,"error")) stop(seeded$fit)
      fit<-seeded$fit; saveRDS(fit,file.path(dest,"truth-start-fit.rds"))
      fit.summary<-.tb_partition(fit$membership,truth$group)
      blind.summary<-.tb_partition(blind$membership,truth$group)
      full<-data.frame(ID=id,TRUTH_WC84=row$TRUTH_WC84,
        BLIND_ARI=blind.summary$ari,TRUTH_START_ARI=fit.summary$ari,FROZEN_ARI=frozen.summary$ari,
        BLIND_OBJECTIVE=blind$objective,TRUTH_START_OBJECTIVE=fit$objective,
        OBJECTIVE_GAIN=fit$objective-blind$objective,LOGLIK_GAIN=fit$loglik-blind$loglik,
        FROZEN_MEAN_MAX=frozen.summary$mean.max.membership,
        FITTED_MEAN_MAX=fit.summary$mean.max.membership,
        TRUTH_START_CONVERGED=fit$converged,ITERATIONS=fit$iterations,WARNINGS=seeded$warnings)
      saveRDS(full,file.path(dest,"full-comparison.rds")); print(full)
      folds<-blind.cv$folds[blind$sample.ids]; fold.rows<-list()
      for(fold in unique(folds)) {
        train<-folds!=fold
        result<-.sm_audit_fit(.sm_fit_gds(row$GDS,p,proportion=pi,likelihood="GL",
          sample.weights=setNames(as.numeric(train),blind$sample.ids),
          n.starts=1,seed=66000+j,regularization="pooled",prior.strength=2,
          max.iter=300,tolerance=1e-8,epsilon=blind$epsilon))
        ok<-!inherits(result$fit,"error")&&isTRUE(result$fit$converged)
        score<-if(ok) sum(result$fit$sample.loglik[!train]) else NA_real_
        saveRDS(result,file.path(dest,paste0("truth-start-fold-",fold,".rds")))
        old<-blind.cv$diagnostics[blind.cv$diagnostics$FOLD==fold,]
        fold.rows[[fold]]<-data.frame(ID=id,FOLD=fold,N_TEST=sum(!train),CONVERGED=ok,
          TRUTH_START_K3=score,BLIND_K3=old$TEST_LOGLIK[old$K==3],
          FITTED_K1=old$TEST_LOGLIK[old$K==1],FROZEN_K3=sum(frozen$sample.loglik[!train]),
          WARNINGS=result$warnings,ERROR=if(inherits(result$fit,"error")) conditionMessage(result$fit) else "")
        cat("Completed fold",fold,format(Sys.time()),"\n")
      }
      fold.rows<-do.call(rbind,fold.rows)
      saveRDS(fold.rows,file.path(dest,"fold-comparison.rds"))
      list(full=full,folds=fold.rows)
    },error=function(e) list(error=conditionMessage(e),ID=id))
    saveRDS(outcome,file.path(dest,"result.rds")); outcome
  }
  results<-parallel::mclapply(selected,run.one,mc.cores=workers,
    mc.preschedule=FALSE,mc.set.seed=FALSE)
  names(results)<-manifest$ID[selected]
  saveRDS(results,file.path(output,"results.rds"))
  stopifnot(identical(input.hash,tools::md5sum(manifest$GDS[selected])))
  errors<-vapply(results,function(x) inherits(x,"try-error")||!is.null(x$error),logical(1))
  if(any(errors)) {print(results[errors]); stop("Diagnostic errors retained; do not drop failed datasets.")}
  full<-do.call(rbind,lapply(results,`[[`,"full"))
  folds<-do.call(rbind,lapply(results,`[[`,"folds"))
  saveRDS(full,file.path(output,"full-comparisons.rds"))
  saveRDS(folds,file.path(output,"fold-comparisons.rds"))
  writeLines(c("Privileged truth diagnostic: not an operational K-selection procedure",
    capture.output(print(full,row.names=FALSE)),capture.output(print(folds,row.names=FALSE))),
    file.path(output,"comparisons.txt"))
  print(full); print(folds)
}
if(sys.nframe()==0L) .tb_main(commandArgs(TRUE))
