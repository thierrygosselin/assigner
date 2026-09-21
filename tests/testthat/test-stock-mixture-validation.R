sm_validation_gds <- function(folder, structured = TRUE) {
  set.seed(150)
  n <- 24; m <- 30
  probability <- if (structured) rep(c(.1, .9), each = n / 2) else rep(.35, n)
  g <- matrix(stats::rbinom(n * m, 2, rep(probability, m)), n, m)
  header <- c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihood">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", paste0("s", 1:n)), collapse = "\t"))
  records <- vapply(1:m, function(l) paste(c("1", l, paste0("v", l), "A", "C", ".",
    "PASS", ".", "PL", c("0,80,80", "80,0,80", "80,80,0")[g[, l] + 1]), collapse = "\t"), character(1))
  vcf <- file.path(folder, "test.vcf"); path <- file.path(folder, "test.gds")
  writeLines(c(header, records), vcf)
  SeqArray::seqVCF2GDS(vcf, path, verbose = FALSE)
  path
}

test_that("exact frequency alignment handles permuted components", {
  x <- matrix(c(.1, .2, .8, .9, .4, .6), 2, 3)
  alignment <- .sm_align(x, x[, c(3, 1, 2)])
  expect_equal(alignment$order, c(2, 3, 1))
  expect_equal(alignment$mean.squared.distance, 0)
})

test_that("regularization and mixed seeds propagate through CV and bootstrap", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive=TRUE))
  path <- sm_validation_gds(folder)
  initial <- matrix(rep(c(.1,.9),each=30),30,2,
    dimnames=list(as.character(1:30),c("A","B")))
  fit <- .sm_fit_gds(path,initial,n.starts=7,start.strategy="mixed",
    regularization="pooled",prior.strength=2)
  boot <- .sm_bootstrap_gds(path,fit,B=2,n.starts=1,seed=91)
  expect_true(all(boot$diagnostics$USABLE))
  expect_false(is.null(boot$intervals))
  folds <- setNames(rep(1:3,8),paste0("s",1:24))
  cv <- .sm_cv_gds(path,K=1:2,folds=folds,n.starts=7,start.strategy="mixed",
    regularization="pooled",prior.strength=2)
  expect_true(all(cv$summary$COMPLETE))
  expect_equal(cv$regularization,"pooled")
  expect_equal(cv$prior.strength,2)
  expect_equal(cv$best.predictive.K,2)
  expect_equal(cv$summary$TOTAL_LOGLIK,
    as.numeric(tapply(cv$heldout$LOGLIK,cv$heldout$K,sum)))
})

test_that("zero-weight holdouts do not estimate mixture parameters", {
  x <- array(-Inf, c(4, 3, 1))
  g <- c(0, 1, 2, 2)
  for (i in 1:4) x[i, g[i] + 1, 1] <- 0
  prepared <- .sm_prepare(x, "log")
  reader <- function(index) prepared
  fit <- .sm_block_fit(reader, list(1), matrix(.5, 1, 1), 1, rep(1, 4),
                       20, 1e-8, 1e-8, c(1, 1, 0, 0))
  expect_equal(as.numeric(fit$frequency), .25)
  expect_equal(fit$sample.loglik, dbinom(g, 2, .25, log = TRUE))
  expect_equal(fit$loglik, sum(dbinom(g[1:2], 2, .25, log = TRUE)))
})

test_that("fixed-K bootstrap returns aligned intervals and preserves state", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  path <- sm_validation_gds(folder)
  initial <- matrix(rep(c(.1, .9), each = 30), 30, 2,
    dimnames = list(as.character(1:30), c("A", "B")))
  fit <- .sm_fit_gds(path, initial, n.starts = 1)
  set.seed(17); before <- .Random.seed
  a <- .sm_bootstrap_gds(path, fit, B = 4, n.starts = 1, seed = 91,
                          frequency.variants = c("1", "2"))
  b <- .sm_bootstrap_gds(path, fit, B = 4, n.starts = 1, seed = 91,
                          frequency.variants = c("1", "2"))
  expect_identical(.Random.seed, before)
  expect_equal(a, b)
  expect_true(all(a$diagnostics$CONVERGED))
  expect_equal(dim(a$intervals$membership), c(3, 24, 2))
  expect_equal(dim(a$intervals$frequency), c(3, 2, 2))
  expect_true(all(a$intervals$membership >= 0 & a$intervals$membership <= 1))
  expect_error(.sm_bootstrap_gds(path, fit, B = 4, max.values = 1), "storage")
})

test_that("held-out scoring distinguishes structured and null smoke fixtures", {
  for (structured in c(TRUE, FALSE)) {
    folder <- tempfile(); dir.create(folder)
    path <- sm_validation_gds(folder, structured)
    folds <- stats::setNames(rep(1:3, 8), paste0("s", 1:24))
    cv <- .sm_cv_gds(path, K = 1:2, folds = folds, n.starts = 3,
                     seed = 22, max.iter = 500)
    expect_true(all(cv$summary$COMPLETE))
    expect_equal(cv$summary$N_TEST, c(24, 24))
    expect_equal(cv$best.predictive.K, if (structured) 2 else 1)
    expect_equal(nrow(cv$heldout), 48)
    unlink(folder, recursive = TRUE)
  }
})

test_that("incomplete resampling and K comparisons do not produce conclusions", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  path <- sm_validation_gds(folder, FALSE)
  initial <- matrix(.5, 30, 1, dimnames = list(as.character(1:30), "A"))
  fit <- .sm_fit_gds(path, initial, n.starts = 1)
  gds <- SeqArray::seqOpen(path)
  on.exit(SeqArray::seqClose(gds), add = TRUE, after = FALSE)
  before <- SeqArray::seqGetFilter(gds)
  expect_warning(boot <- .sm_bootstrap_gds(gds, fit, B = 2, n.starts = 1,
    max.iter = 1, seed = 94), "intervals withheld")
  expect_null(boot$intervals)
  expect_true(any(!boot$diagnostics$USABLE))
  folds <- stats::setNames(rep(1:3, 8), paste0("s", 1:24))
  expect_warning(cv <- .sm_cv_gds(gds, K = 1:2, folds = folds, n.starts = 2,
    max.iter = 1, seed = 82), "no best K")
  expect_true(is.na(cv$best.predictive.K))
  expect_identical(SeqArray::seqGetFilter(gds), before)
  expect_error(.sm_cv_gds(gds, K = 2, folds = folds), "including 1")
})
