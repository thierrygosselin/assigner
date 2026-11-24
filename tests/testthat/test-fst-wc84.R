fst_fixture <- function() {
  x <- tidyr::expand_grid(MARKERS = as.character(1:12), INDIVIDUALS = paste0('s',1:8))
  x$STRATA <- rep(rep(c('A','B'),each=4),12)
  x$GT <- rep(c('001001','001001','001002','002002','001002','002002','002002','002002'),12)
  x$CHROM <- rep(rep(c('1','2'),each=6),each=8)
  x$POS <- rep(seq_len(12)*1000,each=8)
  x
}

test_that('current strata API replaces labels and respects the sample whitelist', {
  d <- fst_fixture()
  s <- tibble::tibble(INDIVIDUALS=paste0('s',1:7),STRATA=rep(c('X','Y'),c(4,3)))
  z <- prep_strata(d,s,pop.levels=c('Y','X'),blacklist.id='s1')
  expect_setequal(unique(z$INDIVIDUALS),paste0('s',2:7))
  expect_identical(levels(z$STRATA),c('Y','X'))
  expect_false(any(as.character(z$STRATA) %in% c('A','B')))
})

test_that('non-pairwise WC84 completes, exports and preserves the future plan', {
  skip_if_not_installed('fst')
  work <- tempfile();dir.create(work);withr::local_dir(work)
  future::plan(future::sequential)
  before <- future::plan()
  z <- fst_WC84(fst_fixture(),pairwise=FALSE,ci=FALSE,parallel.core=1,
                linkage.group=TRUE,filename='test',verbose=TRUE)
  expect_true(is.finite(z$fst.overall$FST))
  expect_equal(z$fst.overall$FST,
    sum(z$sigma.loc$lsiga)/sum(z$sigma.loc$lsiga+z$sigma.loc$lsigb+z$sigma.loc$lsigw),tolerance=1e-8)
  expect_setequal(z$fst.linkage.groups$CHROM,c('1','2'))
  expect_true(length(list.files(work,pattern='fst.overall',recursive=TRUE))>0)
  expect_identical(class(future::plan()),class(before))
})

test_that('pairwise output remains available with current metadata API', {
  work <- tempfile();dir.create(work);withr::local_dir(work)
  z <- fst_WC84(fst_fixture(),pairwise=TRUE,ci=FALSE,parallel.core=1,
                linkage.group=TRUE,filename='pair',verbose=FALSE)
  expect_equal(nrow(z$pairwise.fst),1L)
  expect_equal(as.numeric(z$pairwise.fst$FST),z$fst.overall$FST,tolerance=1e-8)
})

test_that('GDS import respects active samples and variants and restores filters', {
  skip_if_not_installed('SeqArray')
  vcf <- tempfile(fileext='.vcf'); gds <- tempfile(fileext='.gds')
  writeLines(c('##fileformat=VCFv4.2',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3',
    '1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1',
    '1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t./.\t1/1\t0/1'),vcf)
  SeqArray::seqVCF2GDS(vcf,gds,verbose=FALSE)
  g <- SeqArray::seqOpen(gds,readonly=TRUE)
  on.exit({SeqArray::seqClose(g);unlink(c(vcf,gds))},add=TRUE)
  testthat::local_mocked_bindings(
    extract_individuals_metadata=function(...) tibble::tibble(INDIVIDUALS=c('s1','s2','s3'),STRATA=c('A','A','B')),
    extract_markers_metadata=function(...) tibble::tibble(VARIANT_ID=1:2),
    .package='genometranslator')
  SeqArray::seqSetFilter(g,sample.id=c('s1','s2'),variant.id=2L,verbose=FALSE)
  before <- SeqArray::seqGetFilter(g)
  d <- import_data(g)
  expect_identical(d$INDIVIDUALS,c('s1','s2'))
  expect_equal(d$GT,c('000000','002002'))
  expect_equal(d$VARIANT_ID,c(2L,2L))
  expect_identical(before,SeqArray::seqGetFilter(g))
})
