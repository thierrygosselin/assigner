test_that("streamed EM agrees with dense EM on nontrivial likelihoods", {
  set.seed(52)
  x <- array(-rexp(12 * 3 * 7), c(12, 3, 7))
  x[1, , ] <- NA_real_
  x[2, , 4] <- 0
  frequency <- matrix(runif(14, .2, .8), 7, 2)
  prepared <- .sm_prepare(x, "log")
  reader <- function(index) list(
    log.likelihood = prepared$log.likelihood[, , index, drop = FALSE],
    informative = prepared$informative[, index, drop = FALSE])
  information <- rowSums(prepared$informative)
  full <- .sm_fit(x, frequency, encoding = "log", max.iter = 30)
  for (size in c(1, 3, 7)) {
    blocks <- split(1:7, ceiling((1:7) / size))
    streamed <- .sm_block_fit(reader, blocks, frequency, c(.5, .5), information,
                               30, 1e-8, 1e-8)
    expect_equal(streamed$frequency, full$frequency, tolerance = 1e-9)
    expect_equal(streamed$membership, full$membership, tolerance = 1e-9)
    expect_equal(streamed$history, full$history, tolerance = 1e-9)
    expect_true(all(diff(streamed$history) >= -1e-8))
  }
})

test_that("RNG state is restored even when every start fails", {
  set.seed(291)
  before <- .Random.seed
  reader <- function(index) stop("test read failure")
  expect_error(.sm_multistart(reader, list(1L), matrix(.5, 1, 1), 1,
    c(1, 1), 2, 18, 2, 1e-8, 1e-8), "All starts failed")
  expect_identical(.Random.seed, before)
})

test_that("a failed start is retained while subsequent starts can succeed", {
  x <- array(c(0, -5, -10, -10, -5, 0), c(2, 3, 1))
  prepared <- .sm_prepare(x, "log")
  calls <- 0L
  reader <- function(index) {
    calls <<- calls + 1L
    if (calls == 1L) stop("Injected first-start failure")
    prepared
  }
  expect_warning(fit <- .sm_multistart(reader, list(1L), matrix(.5, 1, 1), 1,
    c(1, 1), 3, 27, 500, 1e-8, 1e-8), "starts failed")
  expect_true(fit$converged)
  expect_identical(fit$start.diagnostics$ERROR[1], "Injected first-start failure")
  expect_true(is.na(fit$start.diagnostics$LOGLIK[1]))
  expect_true(fit$selected.start %in% 2:3)
})
