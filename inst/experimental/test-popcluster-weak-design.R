source("../grur/R/simulation_benchmark.R")
source("inst/experimental/prepare-ppr-popcluster.R")
source("inst/experimental/diagnose-truth-starts.R")
a <- .grur_structure_truth(rep(10,3),20,.001,911003)
b <- .grur_structure_truth(rep(10,3),20,.001,911003)
c <- .grur_structure_truth(rep(10,3),20,.001,911004)
stopifnot(identical(a,b),!identical(a$frequency,c$frequency),
          !identical(a$genotype,c$genotype),all(table(a$group)==10))
stopifnot(identical(unname(ppr_decode(ppr_encode(a$genotype))),unname(a$genotype)))
stopifnot(.tb_ari(c(1,1,2,2),c(2,2,1,1))==1)
if(requireNamespace("mclust",quietly=TRUE)) {
  set.seed(1)
  for(i in 1:100) {
    x <- sample(1:3,300,TRUE); y <- sample(1:3,300,TRUE)
    stopifnot(isTRUE(all.equal(.tb_ari(x,y),mclust::adjustedRandIndex(x,y))))
  }
}
cat("Simulation reproducibility, distinct frequencies, encoding and agreement checks passed.\n")
