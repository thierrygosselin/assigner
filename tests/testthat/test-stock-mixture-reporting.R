sm_reporting_gds <- function(folder) {
  set.seed(150)
  g<-matrix(rbinom(24*30,2,rep(c(.1,.9),each=12)),24,30)
  header<-c("##fileformat=VCFv4.2","##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihood">',
    paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",paste0("s",1:24)),collapse="\t"))
  records<-vapply(1:30,function(l)paste(c("1",l,paste0("v",l),"A","C",".","PASS",".","PL",
    c("0,80,80","80,0,80","80,80,0")[g[,l]+1]),collapse="\t"),character(1))
  vcf<-file.path(folder,"test.vcf");path<-file.path(folder,"test.gds")
  writeLines(c(header,records),vcf);SeqArray::seqVCF2GDS(vcf,path,verbose=FALSE)
  path
}

test_that("partition agreement is label invariant and excludes ties", {
  expect_equal(unname(.sm_partition_agreement(c(1,1,2,2),c(2,2,1,1))),c(1,4))
  expect_equal(unname(.sm_partition_agreement(c(1,NA,2),c(2,1,1))),c(1,2))
  expect_true(is.na(.sm_partition_agreement(c(NA,1),c(1,2))[1]))
  expect_equal(.sm_partition_agreement(rep(1,4),rep(2,4))[1],c(ARI=1))
  expect_lt(.sm_partition_agreement(c(1,1,2,2),c(1,2,1,2))[1],1)
})

test_that("fit reporting exposes ambiguity without changing assignments", {
  skip_if_not_installed("SeqArray")
  folder<-tempfile();dir.create(folder);on.exit(unlink(folder,recursive=TRUE))
  path<-sm_reporting_gds(folder)
  initial<-matrix(rep(c(.1,.9),each=30),30,2,
    dimnames=list(as.character(1:30),c("A","B")))
  fit<-.sm_fit_gds(path,initial,n.starts=2)
  expect_true(all(fit$assignments$NORMALIZED_ENTROPY>=0 & fit$assignments$NORMALIZED_ENTROPY<=1))
  expect_true(all(fit$assignments$PROBABILITY_MARGIN>=0))
  expect_true(all(fit$assignments$ASSIGNMENT_STATUS=="conditional_assignment"))
  expect_equal(fit$start.diagnostics$PARTITION_ARI_VS_SELECTED[fit$start.diagnostics$SELECTED],1)
  expect_equal(fit$start.diagnostics$N_COMPARABLE[fit$start.diagnostics$SELECTED],24)
  one<-.sm_fit_gds(path,initial[,1,drop=FALSE],n.starts=1)
  expect_true(all(is.na(one$assignments$NORMALIZED_ENTROPY)))
  expect_true(all(one$assignments$ASSIGNMENT_STATUS=="single_component_not_evidence_of_panmixia"))
  folds<-setNames(rep(1:3,8),paste0("s",1:24))
  cv<-.sm_cv_gds(path,K=1:2,folds=folds,n.starts=2)
  expect_equal(nrow(cv$fold.support),6)
  expect_equal(min(cv$summary$GAP_FROM_BEST),0)
  expect_true(all(cv$fold.support$FOLD_COMPLETE))
  expect_match(cv$support.interpretation,"no calibrated K confidence")
  limited<-.sm_fit_gds(path,initial,n.starts=4,max.genotypes=72,block.snps=1)
  expect_true(all(is.na(limited$start.diagnostics$PARTITION_ARI_VS_SELECTED)))
  expect_match(limited$start.partition.assessment,"Not assessed")
})
