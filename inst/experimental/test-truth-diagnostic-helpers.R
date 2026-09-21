source("inst/experimental/diagnose-truth-starts.R")
stopifnot(.tb_ari(c(1,1,2,2),c(2,2,1,1))==1,
  abs(.tb_ari(c(1,1,2,2),c(1,2,1,2))+.5)<1e-12,
  .tb_ari(rep(1,4),rep(2,4))==1,
  is.na(.tb_ari(1,1)))
prob<-rbind(c(.9,.1),c(.8,.2),c(.1,.9),c(.2,.8))
s<-.tb_partition(prob,c(1,1,2,2))
stopifnot(s$ari==1,s$assigned==4,abs(s$mean.max.membership-.85)<1e-12)
s<-.tb_partition(matrix(.5,4,2),c(1,1,2,2))
stopifnot(s$assigned==0,is.na(s$ari))
cat("Truth-diagnostic partition, label-invariance and tie checks passed\n")
source("R/stock_mixture_blocks.R")
set.seed(1)
x<-array(-runif(24),c(4,3,2))
reader<-function(index) list(log.likelihood=x[,,index,drop=FALSE])
one<-.sm_block_estep(reader,list(1:2),matrix(c(.2,.7),2,1),1,4)
three<-.sm_block_estep(reader,list(1:2),matrix(rep(c(.2,.7),3),2,3),rep(1/3,3),4)
stopifnot(isTRUE(all.equal(one$sample.loglik,three$sample.loglik,tolerance=1e-12)))
cat("Identical-frequency frozen K=1/K=3 score equivalence passed\n")
