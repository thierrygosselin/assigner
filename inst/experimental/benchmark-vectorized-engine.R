# Experimental benchmark accelerator, not an exported or installed backend.
# Same diploid likelihood EM; stores one n-by-marker dosage block per component.
.bv_moments <- function(block, p) {
  n<-dim(block$log.likelihood)[1]; m<-length(p)
  prior<-vapply(p,function(q) stats::dbinom(0:2,2,q,log=TRUE),numeric(3))
  a<-sweep(matrix(block$log.likelihood[,1,],n,m),2,prior[1,],"+")
  b<-sweep(matrix(block$log.likelihood[,2,],n,m),2,prior[2,],"+")
  c<-sweep(matrix(block$log.likelihood[,3,],n,m),2,prior[3,],"+")
  peak<-pmax(a,b,c)
  normalizer<-peak+log(exp(a-peak)+exp(b-peak)+exp(c-peak))
  dosage<-exp(b-normalizer)+2*exp(c-normalizer)
  normalizer[!block$informative]<-0
  list(loglik=normalizer,dosage=dosage)
}
.bv_estep <- function(read.block,blocks,frequency,proportion,n,
                      sample.weights=rep(1,n)) {
  conditional<-matrix(0,n,ncol(frequency))
  for(index in blocks) {
    block<-read.block(index)
    for(j in seq_len(ncol(frequency)))
      conditional[,j]<-conditional[,j]+rowSums(.bv_moments(block,frequency[index,j])$loglik)
  }
  joint<-sweep(conditional,2,log(proportion),"+")
  peak<-apply(joint,1,max)
  normalizer<-peak+log(rowSums(exp(joint-peak)))
  list(membership=exp(joint-normalizer),loglik=sum(normalizer*sample.weights),
    sample.loglik=normalizer)
}
.bv_fit <- function(read.block,blocks,frequency,proportion,information,
                    max.iter,tolerance,epsilon,
                    sample.weights=rep(1,length(information)),
                    prior.center=rep(.5,nrow(frequency)),prior.strength=0) {
  n<-length(information)
  penalty<-function(p) if(prior.strength==0) 0 else
    prior.strength*sum(prior.center*log(p)+(1-prior.center)*log1p(-p))
  state<-.bv_estep(read.block,blocks,frequency,proportion,n,sample.weights)
  history<-state$loglik; objective<-state$loglik+penalty(frequency)
  objective.history<-objective; converged<-FALSE
  for(iteration in seq_len(max.iter)) {
    tau<-state$membership
    fit.weights<-sample.weights*(information>0L)
    proportion<-colSums(tau*fit.weights)/sum(fit.weights)
    if(any(proportion<=0)) stop("Empty component; use a different initialization.")
    for(index in blocks) {
      block<-read.block(index)
      for(j in seq_len(ncol(frequency))) {
        moments<-.bv_moments(block,frequency[index,j])
        weights<-block$informative*(tau[,j]*sample.weights)
        counts<-colSums(weights)
        valid<-counts>0|prior.strength>0
        update<-(colSums(weights*moments$dosage)+prior.strength*prior.center[index])/
          (2*counts+prior.strength)
        frequency[index[valid],j]<-pmin(1-epsilon,pmax(epsilon,update[valid]))
      }
    }
    next.state<-.bv_estep(read.block,blocks,frequency,proportion,n,sample.weights)
    next.objective<-next.state$loglik+penalty(frequency)
    change<-next.objective-objective
    if(change < -1e-7*(1+abs(objective))) stop("EM objective decreased.")
    history<-c(history,next.state$loglik); objective.history<-c(objective.history,next.objective)
    objective<-next.objective; state<-next.state
    if(abs(change)<=tolerance*(1+abs(objective))) {converged<-TRUE; break}
  }
  list(frequency=frequency,proportion=proportion,membership=state$membership,
    loglik=state$loglik,history=history,converged=converged,iterations=iteration,
    epsilon=epsilon,sample.loglik=state$sample.loglik,objective=objective,
    objective.history=objective.history,prior.center=prior.center,prior.strength=prior.strength)
}
.bv_context <- function(gds.path,max.cells=1e6) {
  g<-SeqArray::seqOpen(gds.path,readonly=TRUE); on.exit(SeqArray::seqClose(g))
  ids<-as.character(SeqArray::seqGetData(g,"sample.id"))
  variants<-as.character(SeqArray::seqGetData(g,"variant.id"))
  n<-length(ids); m<-length(variants)
  if(n*m>max.cells) stop("Benchmark cache exceeds max.cells; use the streaming implementation.")
  raw<-SeqArray::seqGetData(g,"annotation/format/GL")
  if(any(raw$length!=3)||!identical(dim(raw$data),c(as.integer(n),as.integer(3*m))))
    stop("Expected diploid GL input.")
  x<-array(as.numeric(raw$data),c(n,3,m)); x[is.nan(x)]<-NA_real_
  prepared<-.sm_prepare(x,"GL")
  reader<-function(index) list(log.likelihood=prepared$log.likelihood[,,index,drop=FALSE],
    informative=prepared$informative[,index,drop=FALSE])
  list(ids=ids,variants=variants,n=n,m=m,reader=reader,
    information=rowSums(prepared$informative),informative=prepared$informative,
    blocks=split(seq_len(m),ceiling(seq_len(m)/250)))
}
.bv_multifit <- function(context,K,strength,weights=rep(1,context$n),seed=1,
                        n.starts=7,max.iter=300) {
  stopifnot(length(weights)==context$n,all(is.finite(weights)),all(weights>=0),sum(weights)>0)
  if(any(colSums(context$informative*weights)==0)) stop("No training information for a marker.")
  pooled<-.bv_fit(context$reader,context$blocks,matrix(.5,context$m,1),1,
    context$information,max.iter,1e-8,1e-8,sample.weights=weights)
  if(!pooled$converged) stop("Pooled fit did not converge.")
  initial<-matrix(.5,context$m,K,dimnames=list(context$variants,paste0("STOCK_",seq_len(K))))
  # .sm_multistart retains the production initialization machinery. This isolated
  # benchmark process aliases its block engine to .bv_fit; installed R is untouched.
  .sm_multistart(context$reader,context$blocks,initial,rep(1/K,K),
    context$information,n.starts,seed,max.iter,1e-8,1e-8,sample.weights=weights,
    start.strategy="mixed",init.max.snps=1000,init.max.cells=1e6,
    prior.center=pooled$frequency[,1],prior.strength=strength)
}
