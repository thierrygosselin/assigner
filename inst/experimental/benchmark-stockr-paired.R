# Shared blind starts for a clean hard-genotype / near-certain-GL comparison.
# Not a replacement for native stockR defaults, and not a GL implementation in stockR.
.bp_starts <- function(g,K,seed,n.starts=7L) {
  if(K==1) return(list(rep(1L,nrow(g))))
  .grur_benchmark_seed(seed,{
    index<-unique(round(seq(1,ncol(g),length.out=min(1000,ncol(g)))))
    z<-scale(g[,index,drop=FALSE],center=TRUE,scale=FALSE)
    sv<-svd(z,nu=min(10,nrow(g)-1,ncol(g)),nv=0)
    rank<-sum(sv$d>max(sv$d)*1e-10)
    if(rank<1) stop("No variable initialization space.")
    scores<-sweep(sv$u[,seq_len(min(10,rank)),drop=FALSE],2,sv$d[seq_len(min(10,rank))],"*")
    lapply(seq_len(n.starts),function(s) {
      if(s<=3) stats::kmeans(scores[,seq_len(min(c(2,5,10)[s],ncol(scores))),drop=FALSE],
        centers=K,nstart=20)$cluster else sample(rep(seq_len(K),length.out=nrow(g)))
    })
  })
}
.bp_stock <- function(g,K,starts,max.iter=300L) {
  diagnostics<-list(); best<-NULL
  for(s in seq_along(starts)) {
    warning<-character(); before<-proc.time()[3]
    fit<-tryCatch(withCallingHandlers(stockR::stockSTRUCTURE(t(g),K=K,
      sample.grps=factor(seq_len(nrow(g)),levels=seq_len(nrow(g))),start.grps=starts[[s]],
      control=list(method="EM",quiet=TRUE,EM.maxit=max.iter,EM.eps=1e-8,EM.tau.step=1)),
      warning=function(w){warning<<-c(warning,conditionMessage(w));invokeRestart("muffleWarning")}),
      error=function(e)e)
    error<-inherits(fit,"error")
    # stockR returns no convergence flag. Beyond maxit the implementation exits
    # even without parameter convergence; retain this limitation explicitly.
    # stockR trims the final trace entry on early convergence, but retains the
    # initial entry at/above the iteration cap. Account for both representations.
    iterations<-if(error) NA_integer_ else if(K==1) 0L else {
      size<-length(fit$margLogLTrace)
      if(size<=max.iter) size else size-1L
    }
    ok<-!error&&is.finite(fit$margLogL)&&iterations<=max.iter
    diagnostics[[s]]<-data.frame(START=s,CONVERGED=ok,ITERATIONS=iterations,
      LOGLIK=if(error) NA_real_ else fit$margLogL,SECONDS=unname(proc.time()[3]-before),
      WARNING=paste(warning,collapse=" | "),ERROR=if(error) conditionMessage(fit) else "")
    if(ok&&(is.null(best)||fit$margLogL>best$margLogL)) best<-fit
  }
  list(fit=best,starts=do.call(rbind,diagnostics))
}
.bp_gl <- function(context,g,K,starts,strength,max.iter=300L) {
  if(!identical(dim(g),c(as.integer(context$n),as.integer(context$m)))) stop("Mismatched genotypes.")
  # Near-certain GL benchmark only: pooled hard-genotype means are the ML limit.
  # Generic uncertain-GL fitting must estimate these from likelihoods instead.
  center<-colMeans(g)/2; diagnostics<-list(); best<-NULL
  for(s in seq_along(starts)) {
    tau<-if(K==1) matrix(1,nrow(g),1) else stockR:::get.init.taus(starts[[s]],2/(K+1),K)
    frequency<-sweep(crossprod(g,tau),2,2*colSums(tau),"/")
    frequency<-pmax(pmin(frequency,1-1e-8),1e-8)
    before<-proc.time()[3]
    fit<-tryCatch(.sm_block_fit(context$reader,context$blocks,frequency,colMeans(tau),
      context$information,max.iter,1e-8,1e-8,prior.center=center,prior.strength=strength),error=function(e)e)
    error<-inherits(fit,"error"); ok<-!error&&isTRUE(fit$converged)
    diagnostics[[s]]<-data.frame(START=s,CONVERGED=ok,
      ITERATIONS=if(error) NA_integer_ else fit$iterations,
      LOGLIK=if(error) NA_real_ else fit$loglik,
      OBJECTIVE=if(error) NA_real_ else fit$objective,
      SECONDS=unname(proc.time()[3]-before),ERROR=if(error) conditionMessage(fit) else "")
    if(ok&&(is.null(best)||fit$objective>best$objective)) best<-fit
  }
  list(fit=best,starts=do.call(rbind,diagnostics))
}
.bp_stock_predict <- function(fit,g) {
  # Reproduce stockR's fitting likelihood bounds, not bounds depending on test N.
  p<-fit$condMeans; floor<-1/(if(fit$K==1) 2 else 10)/fit$data$nFish
  p<-pmax(pmin(p,1-floor),floor); conditional<-matrix(0,nrow(g),fit$K)
  for(k in seq_len(fit$K)) for(index in split(seq_len(ncol(g)),ceiling(seq_len(ncol(g))/250))) {
    z<-matrix(stats::dbinom(g[,index,drop=FALSE],2,rep(p[index,k],each=nrow(g)),log=TRUE),nrow(g))
    conditional[,k]<-conditional[,k]+rowSums(z)
  }
  joint<-sweep(conditional,2,log(fit$pis),"+"); peak<-apply(joint,1,max)
  score<-peak+log(rowSums(exp(joint-peak)))
  list(membership=exp(joint-score),sample.loglik=score)
}
.bp_write_calls <- function(g,path) {
  if(file.exists(path)||anyNA(g)||any(!g %in% 0:2)) stop("Invalid calls or existing output.")
  vcf<-tempfile(fileext=".vcf"); con<-file(vcf,"wt")
  on.exit({try(close(con),silent=TRUE);unlink(vcf)})
  writeLines(c("##fileformat=VCFv4.2",paste0("##contig=<ID=1,length=",ncol(g),">"),
    '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Near-certain simulated likelihood">',
    paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",rownames(g)),collapse="\t")),con)
  for(index in split(seq_len(ncol(g)),ceiling(seq_len(ncol(g))/250))) {
    rows<-vapply(index,function(l) paste(c("1",l,paste0("M",l),"A","C",".","PASS",".","GL",
      c("0,-80,-80","-80,0,-80","-80,-80,0")[g[,l]+1L]),collapse="\t"),character(1))
    writeLines(rows,con)
  }
  close(con); con<-NULL
  SeqArray::seqVCF2GDS(vcf,path,verbose=FALSE)
  invisible(path)
}
