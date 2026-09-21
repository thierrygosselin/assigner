# Bounded-block sample-Gram PCA for an initialization sensitivity experiment.
# Uses all TRAINING markers, with no evaluation samples or population labels.
.bp_all_scores <- function(g,block.snps=250L,max.sample.cells=1e6) {
  n<-nrow(g)
  if(n*n>max.sample.cells) stop("Sample-Gram matrix exceeds memory guard.")
  if(anyNA(g)||any(!is.finite(g))) stop("Complete finite training data required.")
  gram<-matrix(0,n,n)
  for(index in split(seq_len(ncol(g)),ceiling(seq_len(ncol(g))/block.snps))) {
    z<-sweep(g[,index,drop=FALSE],2,colMeans(g[,index,drop=FALSE]),"-")
    gram<-gram+tcrossprod(z)
  }
  eig<-eigen(gram,symmetric=TRUE)
  keep<-which(eig$values>max(eig$values)*1e-10)
  keep<-head(keep,10L)
  if(!length(keep)) stop("No variable initialization space.")
  sweep(eig$vectors[,keep,drop=FALSE],2,sqrt(eig$values[keep]),"*")
}
.bp_all_starts <- function(g,K,seed,n.starts=7L) {
  if(K==1) return(list(rep(1L,nrow(g))))
  baseline<-.bp_starts(g,K,seed,n.starts)
  .grur_benchmark_seed(seed,{
    scores<-.bp_all_scores(g)
    lapply(seq_len(n.starts),function(s) {
      if(s<=3) stats::kmeans(scores[,seq_len(min(c(2,5,10)[s],ncol(scores))),drop=FALSE],
        centers=K,nstart=20)$cluster else baseline[[s]]
    })
  })
}
