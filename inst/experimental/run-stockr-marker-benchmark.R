# Run from assigner: Rscript .../run-stockr-marker-benchmark.R OUTPUT fixed|cv WORKERS [JOB_IDS]
# JOB_IDS is comma-separated row numbers; omit for all. Completed jobs resume.
# Fixed-K is the first stage; CV is a separate, substantially larger stage.
args<-commandArgs(TRUE)
if(length(args)<3||!args[2]%in%c("fixed","cv")) stop("Supply OUTPUT fixed|cv WORKERS [JOB_IDS].")
output<-args[1]; mode<-args[2]; workers<-as.integer(args[3])
stopifnot(workers>=1,workers<=4)
paths<-c("R/stock_likelihood_mixture.R","R/stock_mixture_blocks.R",
  "inst/experimental/benchmark-streaming-engine.R","inst/experimental/benchmark-stockr-paired.R",
  "inst/experimental/diagnose-truth-starts.R","../grur/R/simulation_benchmark.R")
for(p in paths) source(p)
hash<-tools::md5sum(c(paths,"inst/experimental/run-stockr-marker-benchmark.R"))
dir.create(output,recursive=TRUE,showWarnings=FALSE); output<-normalizePath(output)
audit.path<-file.path(output,"audit.rds")
if(file.exists(audit.path)) {
  old<-readRDS(audit.path)
  if(!identical(old$hash,hash)||old$stockR!=as.character(packageVersion("stockR")))
    stop("Sources/package changed; use new output.")
} else saveRDS(list(hash=hash,stockR=as.character(packageVersion("stockR")),session=sessionInfo()),audit.path)
conditions<-data.frame(CONDITION=c("three_001","three_005","three_010","two_010","null"),
  TRUE_K=c(3L,3L,3L,2L,1L),DRIFT=c(.001,.005,.01,.01,0))
design<-list(markers=c(1000L,3000L,5000L,10000L),replicates=1:2,n.train=300L,n.test=300L,
  conditions=conditions,starts=7L,K=1:4,cv.folds=5L,strength=c(2,20,100),
  fixed.strength=c(0,100),max.iter=300L,
  statement="Two independent frequency realizations; nested panels; clean near-certain GLs; exploratory")
saveRDS(design,file.path(output,"design-locked.rds"))
manifest<-merge(expand.grid(REPLICATE=1:2,MARKERS=design$markers),conditions,by=NULL)
manifest$JOB<-seq_len(nrow(manifest)); manifest$ID<-with(manifest,paste(CONDITION,REPLICATE,MARKERS,sep="_"))
saveRDS(manifest,file.path(output,"manifest.rds"))
selected<-if(length(args)>=4) as.integer(strsplit(args[4],",",fixed=TRUE)[[1]]) else manifest$JOB
stopifnot(all(selected%in%manifest$JOB),!anyDuplicated(selected))
# Prepare each generating realization once, reusing its nested marker panels.
for(c in seq_len(nrow(conditions))) for(r in 1:2) {
  id<-paste(conditions$CONDITION[c],r,sep="_"); dest<-file.path(output,id)
  dir.create(dest,showWarnings=FALSE)
  if(!file.exists(file.path(dest,"ready.rds"))) {
    k<-conditions$TRUE_K[c]; seed<-910000+c*100+r
    truth<-.grur_structure_truth(rep(300/k,k),10000,conditions$DRIFT[c],seed)
    g<-truth$genotype
    test<-.grur_benchmark_seed(seed+10000,matrix(stats::rbinom(300*10000,2,
      as.vector(t(truth$frequency[,truth$group,drop=FALSE]))),300,10000))
    dimnames(test)<-list(paste0("T",seq_len(300)),colnames(g))
    saveRDS(list(train=g,test=test,group=truth$group,frequency=truth$frequency,seed=seed),file.path(dest,"truth.rds"))
    .bp_write_calls(g,file.path(dest,"train.gds")); .bp_write_calls(test,file.path(dest,"test.gds"))
    saveRDS(tools::md5sum(file.path(dest,c("train.gds","test.gds"))),file.path(dest,"ready.rds"))
  }
}
run.one<-function(j) {
  row<-manifest[j,]; dest<-file.path(output,paste0(mode,"_",row$ID))
  dir.create(dest,showWarnings=FALSE)
  if(file.exists(file.path(dest,"result.rds"))) return(readRDS(file.path(dest,"result.rds")))
  log<-file(file.path(dest,"progress.txt"),"wt"); sink(log); sink(log,type="message")
  on.exit({sink(type="message");sink();close(log)})
  result<-tryCatch({
    cat("Starting",mode,row$ID,format(Sys.time()),"\n");flush(log)
    src<-file.path(output,paste(row$CONDITION,row$REPLICATE,sep="_"))
    data<-readRDS(file.path(src,"truth.rds")); g<-data$train[,seq_len(row$MARKERS),drop=FALSE]
    # Invariant training markers are omitted for BOTH methods, matching stockR.
    # Selection uses training observations only; no held-out/test genotypes.
    fit.pair<-function(train,k,strengths,seed) {
      keep<-which(apply(g[train,,drop=FALSE],2,function(x)length(unique(x))>1))
      gt<-g[train,keep,drop=FALSE]
      context<-.bs_context(file.path(src,"train.gds"),sample.ids=rownames(gt),variant.ids=keep)
      on.exit(context$close())
      stopifnot(identical(context$ids,rownames(gt)),identical(context$variants,as.character(keep)))
      starts<-.bp_starts(gt,k,seed,design$starts)
      models<-list(stockR=.bp_stock(gt,k,starts))
      for(s in strengths) models[[paste0("GL_",s)]]<-.bp_gl(context,gt,k,starts,s)
      list(models=models,keep=keep,max.block.cells=context$audit$max.read.cells)
    }
    project<-function(pair,path,calls,ids) {
      context<-.bs_context(path,sample.ids=ids,variant.ids=pair$keep)
      on.exit(context$close())
      stopifnot(identical(context$ids,ids))
      lapply(names(pair$models),function(name) {
        fit<-pair$models[[name]]$fit
        if(is.null(fit)) return(NULL)
        if(name=="stockR") .bp_stock_predict(fit,calls[,pair$keep,drop=FALSE]) else
          .sm_block_estep(context$reader,context$blocks,fit$frequency,fit$proportion,context$n)
      }) |> stats::setNames(names(pair$models))
    }
    if(mode=="fixed") {
      paired<-fit.pair(seq_len(nrow(g)),row$TRUE_K,design$fixed.strength,120000+j)
      saveRDS(paired,file.path(dest,"fits.rds"))
      predictions<-project(paired,file.path(src,"test.gds"),data$test,rownames(data$test))
      saveRDS(predictions,file.path(dest,"predictions.rds"))
      rows<-lapply(names(paired$models),function(name) {
        fit<-paired$models[[name]]; pred<-predictions[[name]]
        z<-if(is.null(pred)) list(ari=NA,assigned=0,mean.max.membership=NA) else .tb_partition(pred$membership,data$group)
        data.frame(row,METHOD=name,RETAINED=length(paired$keep),ARI=z$ari,ASSIGNED=z$assigned,
          MEAN_MAX=z$mean.max.membership,TEST_LOGLIK=if(is.null(pred)) NA_real_ else mean(pred$sample.loglik),
          CONVERGED=!is.null(fit$fit),FAILED_STARTS=sum(!fit$starts$CONVERGED),
          FIT_SECONDS=sum(fit$starts$SECONDS),MAX_BLOCK_CELLS=paired$max.block.cells)
      })
      result<-do.call(rbind,rows)
    } else {
      cv<-list(); counter<-0L; folds<-rep(1:5,length.out=nrow(g))
      saveRDS(folds,file.path(dest,"folds.rds"))
      for(fold in 1:5) for(k in design$K) {
        paired<-fit.pair(which(folds!=fold),k,design$strength,130000+j*100+fold*10+k)
        saveRDS(paired,file.path(dest,paste0("fold",fold,"_K",k,".rds")))
        held<-which(folds==fold)
        pred<-project(paired,file.path(src,"train.gds"),g[held,,drop=FALSE],rownames(g)[held])
        for(name in names(pred)) {
          counter<-counter+1L
          cv[[counter]]<-data.frame(FOLD=fold,K=k,METHOD=name,
            SCORE=if(is.null(pred[[name]])) NA_real_ else mean(pred[[name]]$sample.loglik))
        }
        saveRDS(cv,file.path(dest,"CV-checkpoint.rds"))
        cat("Fold",fold,"K",k,"completed",format(Sys.time()),"\n");flush(log)
      }
      cv<-do.call(rbind,cv)
      if(any(!is.finite(cv$SCORE))) stop("Incomplete CV grid; no selection.")
      grid<-aggregate(SCORE~K+METHOD,cv,mean)
      choose<-function(x) {
        x<-x[x$SCORE>=max(x$SCORE)-1e-8,,drop=FALSE]
        strength<-suppressWarnings(as.numeric(sub("GL_","",x$METHOD)))
        x[order(x$K,strength,na.last=TRUE),][1,,drop=FALSE]
      }
      choices<-rbind(choose(grid[grid$METHOD=="stockR",]),choose(grid[grid$METHOD!="stockR",]))
      saveRDS(choices,file.path(dest,"choices-before-test.rds"))
      rows<-list()
      for(q in seq_len(nrow(choices))) {
        ch<-choices[q,]; s<-if(ch$METHOD=="stockR") numeric() else as.numeric(sub("GL_","",ch$METHOD))
        paired<-fit.pair(seq_len(nrow(g)),ch$K,s,140000+j*10+ch$K)
        pred<-project(paired,file.path(src,"test.gds"),data$test,rownames(data$test))[[ch$METHOD]]
        if(is.null(pred)) stop("Selected full fit failed.")
        z<-.tb_partition(pred$membership,data$group)
        rows[[q]]<-data.frame(row,ch,ARI=z$ari,TEST_LOGLIK=mean(pred$sample.loglik))
        saveRDS(list(pair=paired,prediction=pred),file.path(dest,paste0("selected_",q,".rds")))
      }
      result<-do.call(rbind,rows)
    }
    stopifnot(identical(readRDS(file.path(src,"ready.rds")),tools::md5sum(file.path(src,c("train.gds","test.gds")))))
    print(result); result
  },error=function(e) list(ID=row$ID,error=conditionMessage(e)))
  saveRDS(result,file.path(dest,"result.rds")); result
}
results<-parallel::mclapply(selected,run.one,mc.cores=workers,mc.preschedule=FALSE,mc.set.seed=FALSE)
saveRDS(results,file.path(output,paste0(mode,"-last-batch.rds")))
bad<-vapply(results,function(x)!is.data.frame(x),logical(1))
if(any(bad)) {print(results[bad]);stop("Failures retained; inspect before continuing.")}
print(do.call(rbind,results))
