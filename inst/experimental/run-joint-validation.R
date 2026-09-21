# Controlled first-stage experiment: one independent frequency realization per
# condition, nested training panels, independent test cohort, no truth starts.
# Run from assigner: Rscript .../run-joint-validation.R NEW_OUTPUT
args<-commandArgs(TRUE)
if(length(args)!=1L||file.exists(args[1])) stop("Supply a NEW output directory.")
output<-args[1]
source("R/assign_likelihood_dapc.R")
for(p in c("stock_likelihood_mixture","stock_mixture_blocks","stock_mixture_gds","stock_mixture_validation"))
  source(file.path("R",paste0(p,".R")))
source("inst/experimental/benchmark-vectorized-engine.R")
source("inst/experimental/diagnose-truth-starts.R")
source("../grur/R/simulation_benchmark.R")
# Only this R process uses the tested vectorized block update. The production
# sources and installed packages remain unchanged. Same EM and seven starts.
.sm_block_fit<-.bv_fit
dir.create(output,recursive=TRUE); output<-normalizePath(output)
dir.create(file.path(output,"inputs")); dir.create(file.path(output,"truth"))
design<-list(K=1:4,strength=c(2,20,100),fold.count=c(2L,5L),
  primary.policy=5L,n.starts=7L,max.iter=300L,
  sizes=data.frame(PANEL=c("small","joint_large"),N=c(90L,300L),M=c(100L,500L)),
  conditions=c(weak_three=.01,null_one=0),seeds=c(83211L,83212L),
  statement="One frequency realization per condition; exploratory, not a power study")
saveRDS(design,file.path(output,"design-locked.rds"))
paths<-c("R/stock_likelihood_mixture.R","R/stock_mixture_blocks.R","R/stock_mixture_gds.R",
  "R/stock_mixture_validation.R","inst/experimental/benchmark-vectorized-engine.R",
  "inst/experimental/run-joint-validation.R","../grur/R/simulation_benchmark.R")
saveRDS(list(session=sessionInfo(),sources=tools::md5sum(paths)),file.path(output,"run-audit.rds"))
manifest<-list(); z<-0L
for(condition in names(design$conditions)) {
  k<-if(condition=="weak_three") 3L else 1L
  seed<-design$seeds[match(condition,names(design$conditions))]
  master<-.grur_structure_truth(rep(600L/k,k),500,design$conditions[[condition]],seed)
  train.rows<-unlist(lapply(seq_len(k),function(g) which(master$group==g)[seq_len(300L/k)]))
  test.rows<-setdiff(seq_len(600),train.rows)
  saveRDS(master,file.path(output,"truth",paste0(condition,"_master.rds")))
  for(panel in seq_len(nrow(design$sizes))) {
    n<-design$sizes$N[panel]; m<-design$sizes$M[panel]
    selected<-unlist(lapply(seq_len(k),function(g) train.rows[master$group[train.rows]==g][seq_len(n/k)]))
    stopifnot(!length(intersect(selected,test.rows)))
    subset.truth<-function(rows) {
      x<-master; x$genotype<-master$genotype[rows,seq_len(m),drop=FALSE]
      x$frequency<-master$frequency[seq_len(m),,drop=FALSE]
      x$group<-master$group[rows]; x$metadata<-master$metadata[rows,,drop=FALSE]
      x$sample.size<-as.integer(table(x$group)); x
    }
    train<-subset.truth(selected); test<-subset.truth(test.rows)
    id<-paste(condition,design$sizes$PANEL[panel],sep="_"); z<-z+1L
    train.path<-file.path(output,"inputs",paste0(id,"_train.gds"))
    test.path<-file.path(output,"inputs",paste0(id,"_test.gds"))
    .grur_write_benchmark_gds(.grur_observe_structure(train,seed=seed+100),train$metadata$INDIVIDUALS,train.path)
    .grur_write_benchmark_gds(.grur_observe_structure(test,seed=seed+200),test$metadata$INDIVIDUALS,test.path)
    saveRDS(list(train=train,test=test),file.path(output,"truth",paste0(id,".rds")))
    fst<-.grur_benchmark_fst(train$genotype,train$group)$overall
    manifest[[z]]<-data.frame(ID=id,CONDITION=condition,PANEL=design$sizes$PANEL[panel],
      TRAIN_N=n,TEST_N=300L,MARKERS=m,GENERATING_K=k,DRIFT=design$conditions[[condition]],
      TRAIN_WC84=fst,TRAIN_GDS=train.path,TEST_GDS=test.path)
  }
}
manifest<-do.call(rbind,manifest); saveRDS(manifest,file.path(output,"manifest.rds"))
hash<-tools::md5sum(c(manifest$TRAIN_GDS,manifest$TEST_GDS))
saveRDS(hash,file.path(output,"input-hashes.rds"))

# Reproducible tie handling: scores within 1e-8 per individual tie; prefer the
# smaller K, then the weaker regularizer. Selection never sees evaluation data.
choose.model<-function(x) {
  if(!all(x$COMPLETE)||any(!is.finite(x$SCORE))) return(NULL)
  best<-x[x$SCORE>=max(x$SCORE)-1e-8,]
  best[order(best$K,best$STRENGTH),][1,,drop=FALSE]
}
run.one<-function(j) {
  row<-manifest[j,]; id<-row$ID; dest<-file.path(output,id); dir.create(dest)
  log<-file(file.path(dest,"progress.txt"),"wt"); sink(log); sink(log,type="message")
  on.exit({sink(type="message"); sink(); close(log)})
  result<-tryCatch({
    cat("Starting",id,format(Sys.time()),"\n"); flush(log)
    context<-.bv_context(row$TRAIN_GDS)
    cv<-list(); starts<-list(); counter<-0L
    for(v in design$fold.count) {
      folds<-setNames(rep(seq_len(v),length.out=context$n),context$ids)
      saveRDS(folds,file.path(dest,paste0("folds-",v,".rds")))
      for(fold in seq_len(v)) for(strength in design$strength) for(k in design$K) {
        counter<-counter+1L; weights<-as.numeric(folds!=fold)
        seed<-800000+j*100000+v*10000+fold*100+k
        started<-Sys.time()
        attempt<-.sm_audit_fit(.bv_multifit(context,k,strength,weights,seed,
          n.starts=design$n.starts,max.iter=design$max.iter))
        ok<-!inherits(attempt$fit,"error")&&isTRUE(attempt$fit$converged)
        score<-if(ok) sum(attempt$fit$sample.loglik[folds==fold]) else NA_real_
        cv[[counter]]<-data.frame(ID=id,V=v,FOLD=fold,K=k,STRENGTH=strength,
          N_TEST=sum(folds==fold),LOGLIK=score,CONVERGED=ok,SEED=seed,
          WARNINGS=attempt$warnings,ERROR=if(inherits(attempt$fit,"error")) conditionMessage(attempt$fit) else "",
          SECONDS=as.numeric(difftime(Sys.time(),started,units="secs")))
        starts[[counter]]<-if(inherits(attempt$fit,"error")) NULL else attempt$fit$start.diagnostics
        saveRDS(list(cv=cv,starts=starts),file.path(dest,"cv-checkpoint.rds"))
        cat("V",v,"fold",fold,"strength",strength,"K",k,"converged",ok,format(Sys.time()),"\n"); flush(log)
      }
    }
    cv<-do.call(rbind,cv); saveRDS(cv,file.path(dest,"cv-folds.rds"))
    grid<-do.call(rbind,lapply(split(cv,interaction(cv$V,cv$K,cv$STRENGTH,drop=TRUE)),function(x)
      data.frame(V=x$V[1],K=x$K[1],STRENGTH=x$STRENGTH[1],COMPLETE=all(x$CONVERGED),
        SCORE=if(all(x$CONVERGED)) sum(x$LOGLIK)/sum(x$N_TEST) else NA_real_)))
    saveRDS(grid,file.path(dest,"cv-grid.rds"))
    choices<-list()
    for(v in design$fold.count) for(policy in c("tuned","baseline_strength2")) {
      candidates<-grid[grid$V==v & (policy=="tuned"|grid$STRENGTH==2),]
      chosen<-choose.model(candidates)
      if(is.null(chosen)) stop("Incomplete candidate grid; test evaluation withheld.")
      choices[[paste(v,policy)]]<-data.frame(ID=id,POLICY=policy,V=v,
        K=chosen$K,STRENGTH=chosen$STRENGTH,CV_SCORE=chosen$SCORE)
    }
    choices<-do.call(rbind,choices)
    # Lock all choices BEFORE opening any independent evaluation GL or truth.
    saveRDS(choices,file.path(dest,"selection-locked-before-test.rds"))
    test.context<-.bv_context(row$TEST_GDS)
    truth<-readRDS(file.path(output,"truth",paste0(id,".rds")))
    stopifnot(!length(intersect(context$ids,test.context$ids)),
      identical(context$variants,test.context$variants))
    fitted<-list(); evaluations<-list()
    keys<-unique(rbind(choices[,c("K","STRENGTH")],data.frame(K=1,STRENGTH=2)))
    for(q in seq_len(nrow(keys))) {
      k<-keys$K[q]; strength<-keys$STRENGTH[q]; key<-paste(k,strength,sep="_")
      attempt<-.sm_audit_fit(.bv_multifit(context,k,strength,seed=1800000+j*100+k,
        n.starts=design$n.starts,max.iter=design$max.iter))
      if(inherits(attempt$fit,"error")||!attempt$fit$converged) stop("Selected full-training fit failed.")
      fit<-attempt$fit; fitted[[key]]<-fit
      prediction<-.bv_estep(test.context$reader,test.context$blocks,fit$frequency,fit$proportion,test.context$n)
      saveRDS(list(fit=fit,prediction=prediction,warnings=attempt$warnings),file.path(dest,paste0("final-",key,".rds")))
      recovery<-.tb_partition(prediction$membership,truth$test$group)
      evaluations[[key]]<-list(score=mean(prediction$sample.loglik),ari=recovery$ari,
        mean.max=recovery$mean.max.membership,assigned=recovery$assigned)
    }
    rows<-lapply(seq_len(nrow(choices)),function(q) {
      choice<-choices[q,]; e<-evaluations[[paste(choice$K,choice$STRENGTH,sep="_")]]
      data.frame(choice,TEST_MEAN_LOGLIK=e$score,TEST_GAIN_VS_K1=e$score-evaluations[["1_2"]]$score,
        TEST_ARI=e$ari,TEST_MEAN_MAX=e$mean.max,TEST_ASSIGNED=e$assigned,
        AT_K_BOUNDARY=choice$K==max(design$K))
    })
    result<-do.call(rbind,rows); saveRDS(result,file.path(dest,"evaluation.rds")); print(result)
    list(evaluation=result,grid=grid,folds=cv)
  },error=function(e) list(ID=id,error=conditionMessage(e)))
  saveRDS(result,file.path(dest,"result.rds")); result
}
results<-parallel::mclapply(seq_len(nrow(manifest)),run.one,mc.cores=4,
  mc.preschedule=FALSE,mc.set.seed=FALSE)
names(results)<-manifest$ID
saveRDS(results,file.path(output,"results.rds"))
stopifnot(identical(hash,tools::md5sum(c(manifest$TRAIN_GDS,manifest$TEST_GDS))))
bad<-vapply(results,function(x) inherits(x,"try-error")||!is.null(x$error),logical(1))
if(any(bad)) {print(results[bad]); stop("Failures retained; do not silently omit datasets.")}
evaluation<-do.call(rbind,lapply(results,`[[`,"evaluation"))
saveRDS(evaluation,file.path(output,"evaluation.rds"))
writeLines(c("Locked training-CV choices, independent cohort evaluation; exploratory single realization per condition",
  capture.output(print(manifest,row.names=FALSE)),capture.output(print(evaluation,row.names=FALSE))),
  file.path(output,"results.txt"))
print(evaluation)
