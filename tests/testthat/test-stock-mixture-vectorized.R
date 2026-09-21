test_that("vectorized moments agree with the scalar marker reference", {
  set.seed(903)
  for (n in c(1L, 8L)) for (m in c(1L, 11L)) {
    x <- array(-rexp(n * 3 * m), c(n, 3L, m))
    x[1, , 1] <- 0
    if (n > 1) x[2, 1, 1] <- -Inf
    block <- .sm_prepare(x, "log")
    p <- seq(.001, .999, length.out = m)
    result <- .sm_block_moments(block, p)
    for (l in seq_len(m)) {
      reference <- .sm_marker(matrix(block$log.likelihood[, , l], n, 3), p[l])
      expect_equal(result$loglik[, l], reference$loglik, tolerance = 1e-12)
      expect_equal(result$dosage[, l], reference$dosage, tolerance = 1e-12)
    }
  }
})

test_that("block size preserves weighted penalized EM results", {
  set.seed(907)
  x <- array(-rexp(12 * 3 * 9), c(12, 3, 9))
  x[1, , ] <- NA_real_
  b <- .sm_prepare(x, "GL")
  reader <- function(i) list(log.likelihood = b$log.likelihood[, , i, drop = FALSE],
                            informative = b$informative[, i, drop = FALSE])
  p <- matrix(runif(18, .2, .8), 9, 2)
  for (strength in c(0, 2, 100)) {
    run <- function(blocks) .sm_block_fit(reader, blocks, p, c(.5, .5),
      rowSums(b$informative), 500, 1e-8, 1e-8, rep(c(1, 1, 0), 4),
      prior.center = rep(.4, 9), prior.strength = strength)
    single <- run(as.list(1:9)); full <- run(list(1:9))
    for (key in c("frequency", "membership", "sample.loglik", "objective"))
      expect_equal(single[[key]], full[[key]], tolerance = 1e-7)
    expect_equal(single$converged, full$converged)
  }
})
