test_that("GL and PL encodings give identical mixture fits", {
  pl <- array(c(0, 20, 60, 0, 60, 20, 40, 0, 40, 60, 20, 0), c(4, 3, 1))
  start <- matrix(c(0.2, 0.8), 1, 2)
  a <- .sm_fit(pl, start, encoding = "PL")
  b <- .sm_fit(-pl / 10, start, encoding = "GL")
  expect_equal(a, b)
  expect_true(all(diff(a$history) >= -1e-8))
  expect_equal(rowSums(a$membership), rep(1, 4))
})

test_that("hard calls recover binomial mixture likelihood and frequency", {
  genotype <- c(0, 0, 1, 2)
  x <- array(-Inf, c(4, 3, 1))
  for (i in 1:4) x[i, genotype[i] + 1, 1] <- 0
  p <- matrix(c(0.2, 0.7), 1, 2)
  state <- .sm_estep(x, p, c(0.4, 0.6))
  reference <- cbind(dbinom(genotype, 2, 0.2), dbinom(genotype, 2, 0.7))
  weighted <- sweep(reference, 2, c(0.4, 0.6), "*")
  expect_equal(state$conditional, log(reference))
  expect_equal(state$loglik, sum(log(rowSums(weighted))))
  expect_equal(state$membership, weighted / rowSums(weighted))
  fit <- .sm_fit(x, matrix(0.4, 1, 1), encoding = "log")
  expect_equal(as.numeric(fit$frequency), mean(genotype) / 2)
  expect_true(fit$converged)
})

test_that("flat and missing vectors are neutral to parameter estimation", {
  x <- array(c(0, -2, -8, -8, -2, 0), c(2, 3, 1))
  start <- matrix(c(0.25, 0.75), 1, 2)
  base <- .sm_fit(x, start)
  extra <- array(0, c(3, 3, 1))
  extra[1:2, , ] <- x
  flat <- .sm_fit(extra, start)
  extra[3, , ] <- NA_real_
  missing <- .sm_fit(extra, start)
  expect_equal(flat$frequency, base$frequency)
  expect_equal(flat$proportion, base$proportion)
  expect_equal(flat$membership[3, ], flat$proportion)
  expect_equal(missing, flat)
  expect_error(.sm_fit(array(0, c(2, 3, 1)), start), "Uninformative marker")
})

test_that("normalization and malformed inputs are handled explicitly", {
  x <- array(c(0, -10, -100, -100, -10, 0), c(2, 3, 1))
  start <- matrix(c(0.2, 0.8), 1, 2)
  expect_equal(.sm_fit(x, start), .sm_fit(x - 10000, start))
  x[1, 1, 1] <- NA_real_
  expect_error(.sm_fit(x, start), "partially missing")
  expect_error(.sm_prepare(array(-Inf, c(1, 3, 1))), "Malformed")
  expect_error(.sm_prepare(array(1.5, c(1, 3, 1)), "PL"), "whole numbers")
})
