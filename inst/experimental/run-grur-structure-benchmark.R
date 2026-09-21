# Run from assigner source: Rscript .../run-grur-structure-benchmark.R OUTPUT
# Default fits the first replicate of the clean three-level and K=1/2 controls.
# All staged datasets are generated. --extended fits every generated dataset.
args <- commandArgs(TRUE)
if(!length(args)) stop("Supply a NEW output directory.")
output <- args[1]
if(file.exists(output)) stop("Output already exists; never overwrite a benchmark.")
source(file.path("..","grur","R","simulation_benchmark.R"))
source("R/assign_likelihood_dapc.R")
for(p in c("stock_likelihood_mixture","stock_mixture_blocks","stock_mixture_gds","stock_mixture_validation"))
  source(file.path("R",paste0(p,".R")))
dir.create(output,recursive=TRUE); output <- normalizePath(output)
dir.create(file.path(output,"truth")); dir.create(file.path(output,"observed"))
saveRDS(sessionInfo(),file.path(output,"session-info.rds"))
scenarios <- list()
add <- function(name,n=c(30,30,30),m=100,fst=.01,missing=0,depth=Inf,pilot=FALSE) {
  scenarios[[name]] <<- list(name=name,n=n,m=m,fst=fst,missing=missing,depth=depth,pilot=pilot)
}
for(f in c(.001,.01,.1)) add(paste0("balanced_F",format(f,scientific=FALSE)),fst=f,pilot=TRUE)
add("small_sample",n=c(10,10,10))
add("large_sample",n=c(100,100,100))
add("more_markers",m=500)
add("unequal",n=c(60,34,6))
add("random_missing",missing=.3)
add("structured_missing",missing=rep(c(.05,.25,.5),each=30))
add("low_depth",depth=3)
add("null_K1",n=90,fst=0,pilot=TRUE)
add("control_K2",n=c(45,45),fst=.1,pilot=TRUE)
saveRDS(scenarios,file.path(output,"scenarios.rds"))
manifest <- list(); z <- 0L
for(s in scenarios) for(replica in 1:2) {
  z <- z+1L; key<-paste0(s$name,"_rep",replica)
  truth <- .grur_structure_truth(s$n,s$m,s$fst,seed=16000+z)
  observation <- .grur_observe_structure(truth,s$missing,s$depth,seed=26000+z)
  fst <- .grur_benchmark_fst(truth$genotype,truth$group)
  masked.fst <- .grur_benchmark_fst(observation$masked.truth,truth$group)
  truth$fst <- fst; truth$masked.truth.fst <- masked.fst
  saveRDS(truth,file.path(output,"truth",paste0(key,".rds")))
  saveRDS(observation[c("deliberate.mask","depth","alt.count","uninformative",
                       "mean.depth","error","seed")],file.path(output,"truth",paste0(key,"_observation-audit.rds")))
  path<-file.path(output,"observed",paste0(key,".gds"))
  .grur_write_benchmark_gds(observation,truth$metadata$INDIVIDUALS,path)
  manifest[[z]] <- data.frame(ID=key,SCENARIO=s$name,REPLICATE=replica,
    SAMPLES=sum(s$n),SAMPLE_SIZES=paste(s$n,collapse="/"),MARKERS=s$m,
    GENERATING_K=length(s$n),DRIFT_FST=s$fst,TRUTH_WC84=fst$overall,
    MASKED_TRUTH_WC84=masked.fst$overall,TRUE_FREQUENCY_GST=truth$frequency.GST,
    UNINFORMATIVE=mean(observation$uninformative),MEAN_DEPTH=s$depth,
    TRUTH_SEED=16000+z,OBSERVATION_SEED=26000+z,
    FIT=(s$pilot && replica==1) || "--extended" %in% args,GDS=path)
  cat("Generated",key,"WC84",fst$overall,"\n")
}
manifest<-do.call(rbind,manifest)
saveRDS(manifest,file.path(output,"manifest.rds"))
writeLines(c("Controlled structure benchmark: generating parameters and measured estimates",
  capture.output(print(manifest[,setdiff(names(manifest),"GDS")],row.names=FALSE))),file.path(output,"manifest.txt"))

# Fit without reading truth labels into the estimator--------------------------
ari <- function(a,b) {
  tab<-table(a,b); choose2<-function(x) x*(x-1)/2; n<-sum(tab)
  expected<-sum(choose2(rowSums(tab)))*sum(choose2(colSums(tab)))/choose2(n)
  maximum<-(sum(choose2(rowSums(tab)))+sum(choose2(colSums(tab))))/2
  if(n<2) return(NA_real_)
  if(maximum==expected) return(as.numeric(all(outer(a,a,"==")==outer(b,b,"=="))))
  (sum(choose2(tab))-expected)/(maximum-expected)
}
results<-list()
for(j in which(manifest$FIT)) {
  row<-manifest[j,]; key<-row$ID
  truth<-readRDS(file.path(output,"truth",paste0(key,".rds")))
  initial<-matrix(.5,row$MARKERS,row$GENERATING_K,
    dimnames=list(as.character(seq_len(row$MARKERS)),paste0("STOCK_",seq_len(row$GENERATING_K))))
  fit<-.sm_fit_gds(row$GDS,initial,likelihood="GL",start.strategy="mixed",
    n.starts=7,seed=36000+j,regularization="pooled",prior.strength=2,max.iter=300)
  saveRDS(fit,file.path(output,paste0(key,"_fixed-K-fit.rds")))
  folds<-setNames(rep(1:2,length.out=row$SAMPLES),rownames(truth$genotype))
  cv<-.sm_cv_gds(row$GDS,K=1:5,folds=folds,n.starts=7,seed=46000+j,
    likelihood="GL",start.strategy="mixed",regularization="pooled",prior.strength=2,max.iter=300)
  saveRDS(cv,file.path(output,paste0(key,"_CV.rds")))
  keep<-!is.na(fit$assignments$ASSIGNED_STOCK)
  results[[key]]<-data.frame(ID=key,TRUE_K=row$GENERATING_K,TRUTH_WC84=row$TRUTH_WC84,
    CONVERGED=fit$converged,ASSIGNED=sum(keep),ARI=ari(truth$group[keep],fit$assignments$ASSIGNED_STOCK[keep]),
    PREDICTIVE_K=cv$best.predictive.K,COMPLETE=all(cv$summary$COMPLETE))
  saveRDS(results,file.path(output,"pilot-results.rds"))
  print(results[[key]])
}
results<-do.call(rbind,results)
writeLines(c("Pilot results: regularized mixed-start fits; not a calibrated stock-count procedure",
  capture.output(print(results,row.names=FALSE))),file.path(output,"pilot-results.txt"))
cat("Outputs:",output,"\n")
