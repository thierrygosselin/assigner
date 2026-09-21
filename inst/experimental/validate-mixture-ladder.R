# Reproducible pilot ladder; not a power or interval-coverage study.
# Run from source root. Optional --stockr runs the original unbalanced example.
source("R/stock_likelihood_mixture.R")
source("R/stock_mixture_blocks.R")
marker.arg <- grep("^--markers=", commandArgs(TRUE), value=TRUE)
marker.count <- if(length(marker.arg)) as.integer(sub("--markers=","",marker.arg[1])) else 60L
stopifnot(length(marker.count)==1, !is.na(marker.count), marker.count>=10)
reader_for <- function(x) {
  p <- .sm_prepare(x, "log")
  list(reader = function(i) list(log.likelihood=p$log.likelihood[,,i,drop=FALSE],
    informative=p$informative[,i,drop=FALSE]), information=rowSums(p$informative),
    blocks=list(seq_len(dim(x)[3])))
}
fit_for <- function(f, m, k, weights, strategy, strength, seed, starts=7L) {
  center <- rep(.5,m)
  if (strength>0) {
    pooled <- .sm_block_fit(f$reader,f$blocks,matrix(.5,m,1),1,f$information,
      300,1e-8,1e-8,weights)
    if (!pooled$converged) stop("Pooled fit incomplete")
    center <- pooled$frequency[,1]
  }
  fit <- .sm_multistart(f$reader,f$blocks,matrix(.5,m,k),rep(1/k,k),f$information,
    starts,seed,300,1e-8,1e-8,weights,start.strategy=strategy,
    prior.center=center,prior.strength=strength)
  failed <- fit$start.diagnostics[!is.na(fit$start.diagnostics$ERROR),]
  if(nrow(failed)) {
    cat("Failed individual starts: K",k,"strength",strength,"seed",seed,"\n")
    print(failed[,c("START","METHOD","ERROR")])
  }
  fit
}
simulate_case <- function(case, seed, m=60L) {
  set.seed(seed)
  sizes <- switch(case, K1=60, K2=c(30,30), unbalanced=c(36,4,20), c(20,20,20))
  labels <- rep(seq_along(sizes),sizes); n<-length(labels); k<-length(sizes)
  pooled <- runif(m,.15,.85)
  scale <- if (case=="weak") .25 else 1.8
  p <- plogis(matrix(qlogis(pooled),m,k)+matrix(rnorm(m*k,sd=scale),m,k))
  g <- matrix(rbinom(n*m,2,as.vector(t(p[,labels,drop=FALSE]))),n,m)
  x <- array(-Inf,c(n,3,m))
  for(l in seq_len(m)) for(i in seq_len(n)) {
    if(case=="low_depth") {
      depth<-rpois(1,2); alt<-rbinom(1,depth,c(.01,.5,.99)[g[i,l]+1])
      x[i,,l]<-dbinom(alt,depth,c(.01,.5,.99),log=TRUE)
    } else x[i,g[i,l]+1,l]<-0
    miss <- if(case=="structured_missing") c(.05,.25,.55)[labels[i]] else .02
    if(runif(1)<miss) x[i,,l]<-0
  }
  # Simulation-only fold stratification ensures each generating group is present.
  folds <- unlist(lapply(sizes,function(s) rep(1:2,length.out=s)))
  list(x=x,labels=labels,folds=folds)
}
if ("--stockr" %in% commandArgs(TRUE)) {
  raw <- get("myData",asNamespace("stockR")); labels<-attr(raw,"grps")
  g<-raw[1:300,]; n<-ncol(g); m<-nrow(g)
  x<-array(-Inf,c(n,3,m))
  for(l in 1:m) for(i in 1:n) {
    if(is.na(g[l,i])) x[i,,l]<-0 else x[i,g[l,i]+1,l]<-0
  }
  f<-reader_for(x); records<-list(); z<-0
  for(fold in 1:3) for(k in 1:3) for(strength in c(0,2)) {
    train<-rep(1:3,length.out=n)!=fold
    fit<-fit_for(f,m,k,as.numeric(train),"mixed",strength,54+fold,13L)
    z<-z+1
    records[[z]]<-data.frame(fold,K=k,strength,converged=fit$converged,
      loglik=fit$loglik,heldout=sum(fit$sample.loglik[!train]),
      selected=fit$start.diagnostics$METHOD[fit$selected.start],
      exact_partition=all(outer(max.col(fit$membership[train,,drop=FALSE]),
        max.col(fit$membership[train,,drop=FALSE]),"==")==outer(labels[train],labels[train],"==")))
    print(records[[z]])
  }
  results<-do.call(rbind,records)
  print(aggregate(heldout~K+strength,results,sum))
} else {
  records<-list(); z<-0
  for(case in c("K1","K2","balanced","unbalanced","weak","low_depth","structured_missing"))
    for(replica in 1:2) {
      sim<-simulate_case(case,1200+replica,m=marker.count); f<-reader_for(sim$x); m<-dim(sim$x)[3]
      for(strength in c(0,2)) for(fold in 1:2) for(k in 1:4) {
        train<-sim$folds!=fold
        fit<-tryCatch(fit_for(f,m,k,as.numeric(train),"mixed",strength,710+fold),error=identity)
        ok<-!inherits(fit,"error") && fit$converged
        z<-z+1
        records[[z]]<-data.frame(case,replica,fold,K=k,strength,converged=ok,
          heldout=if(ok) sum(fit$sample.loglik[!train]) else NA_real_,
          exact_partition=if(ok) all(outer(max.col(fit$membership[train,,drop=FALSE]),
            max.col(fit$membership[train,,drop=FALSE]),"==")==outer(sim$labels[train],sim$labels[train],"==")) else NA,
          error=if(inherits(fit,"error")) conditionMessage(fit) else "")
        if(k==4) cat("Completed fold:",case,"replicate",replica,"strength",strength,
                     "fold",fold,"\n")
      }
      cat("Completed",case,"replicate",replica,"\n")
    }
  results<-do.call(rbind,records)
  print(results,row.names=FALSE)
  totals<-aggregate(heldout~case+replica+K+strength,results,function(x) if(length(x)==2) sum(x) else NA)
  winners<-do.call(rbind,lapply(split(totals,interaction(totals$case,totals$replica,totals$strength)),
    function(d) if(nrow(d)==4 && all(is.finite(d$heldout))) d[which.max(d$heldout),] else NULL))
  cat("Complete-comparison winners (no coverage claim):\n"); print(winners,row.names=FALSE)
}
