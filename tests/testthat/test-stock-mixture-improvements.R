mixture_fixture <- function(x) {
  prepared <- .sm_prepare(x, "log")
  list(reader = function(index) list(
    log.likelihood = prepared$log.likelihood[, , index, drop = FALSE],
    informative = prepared$informative[, index, drop = FALSE]),
    information = rowSums(prepared$informative), blocks = list(seq_len(dim(x)[3])))
}

test_that("Beta MAP updates optimize the stated penalized objective", {
  x <- array(rep(c(0, -Inf, -Inf), each = 2), c(2, 3, 1))
  f <- mixture_fixture(x)
  fit <- .sm_block_fit(f$reader, f$blocks, matrix(1e-8, 1, 1), 1,
    f$information, 100, 1e-10, 1e-8, prior.center = .5, prior.strength = 2)
  expect_equal(as.numeric(fit$frequency), 1/6, tolerance = 1e-9)
  expect_true(all(diff(fit$objective.history) >= -1e-9))
  expect_lt(tail(fit$history,1),fit$history[1])
  expect_equal(fit$objective, fit$loglik + log(1/6) + log(5/6))
  expect_equal(fit$loglik, sum(fit$sample.loglik))
  grid <- seq(.001, .999, length.out = 1000)
  expect_gte(fit$objective, max(5*log1p(-grid)+log(grid)))
})

test_that("uncertain-likelihood MAP matches independent scalar optimization", {
  set.seed(513)
  x <- array(-rexp(8*3*3),c(8,3,3))
  f <- mixture_fixture(x); center <- c(.2,.5,.8)
  fit <- .sm_block_fit(f$reader,f$blocks,matrix(.5,3,1),1,f$information,
    2000,1e-12,1e-8,prior.center=center,prior.strength=3)
  prepared <- .sm_prepare(x,"log")$log.likelihood
  independent <- vapply(1:3,function(l) {
    optimize(function(p) {
      emission <- exp(prepared[,,l]) %*% dbinom(0:2,2,p)
      sum(log(emission))+3*(center[l]*log(p)+(1-center[l])*log1p(-p))
    },c(1e-8,1-1e-8),maximum=TRUE,tol=1e-10)$maximum
  },numeric(1))
  expect_true(fit$converged)
  expect_equal(as.numeric(fit$frequency),independent,tolerance=1e-5)
  expect_true(all(diff(fit$objective.history)>=-1e-9))
})

test_that("mixed starts are reproducible and independent of held-out genotypes", {
  set.seed(874)
  g <- matrix(rbinom(24*24, 2, rep(c(.1, .9), each = 12)), 24, 24)
  x <- array(-20, c(24, 3, 24))
  for (l in 1:24) x[cbind(1:24, g[, l]+1, l)] <- 0
  fit_one <- function(x) {
    f <- mixture_fixture(x)
    suppressWarnings(.sm_multistart(f$reader, f$blocks, matrix(.5, 24, 2), c(.5,.5),
      f$information, 7, 93, 200, 1e-8, 1e-8, c(rep(1, 20), rep(0, 4)),
      start.strategy = "mixed", init.max.snps = 20))
  }
  before <- .Random.seed
  fit <- fit_one(x)
  expect_identical(.Random.seed, before)
  expect_equal(fit$frequency, fit_one(x)$frequency)
  x[21:24, , ] <- x[24:21, c(3,2,1), ]
  altered <- fit_one(x)
  expect_equal(fit$frequency, altered$frequency)
  expect_equal(fit$proportion, altered$proportion)
  expect_equal(fit$start.diagnostics$METHOD,
               c("supplied","pca2","pca5","pca10","random","split","perturb"))
  expect_true(fit$converged)
  expect_true(all(is.na(fit$start.diagnostics$ERROR)))
  expect_identical(fit$selected.start, which.max(fit$start.diagnostics$OBJECTIVE))
})

test_that("initialization respects its cell guard and rejects invariant coordinates", {
  f <- mixture_fixture(array(0, c(5, 3, 4)))
  expect_error(.sm_start_space(f$reader, 4, rep(1,5), rep(1,5), max.cells=4),
               "Insufficient")
  expect_error(.sm_start_space(f$reader, 4, rep(1,5), rep(1,5)), "No training information")
})
