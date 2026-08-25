test_that("locuswise missingness uses the called samples at each marker", {
  dat <- expand.grid(
    INDIVIDUALS = paste0("i", 1:8),
    MARKERS = c("m1", "m2"),
    stringsAsFactors = FALSE
  )
  dat$POP_ID <- rep(rep(c("p1", "p2"), each = 4), 2)
  dat$GT_BIN <- c(0, 0, 1, 1, 1, 2, 2, 2, 0, NA, 1, 1, 1, 2, 2, 2)

  z <- amova_genomic(dat, hierarchy = "POP_ID", value = "GT_BIN", missing = "locuswise",
                     min.individuals = 2, standardized = FALSE)
  expect_s3_class(z, "assigner_amova")
  expect_equal(z$per_locus$N, c(8, 7))
  expect_true(all(is.finite(z$variance_components)))
})

test_that("locus and block bootstrap return labelled confidence intervals", {
  dat <- expand.grid(INDIVIDUALS = paste0("i", 1:8), MARKERS = paste0("m", 1:6),
                     stringsAsFactors = FALSE)
  dat$STRATA <- rep(rep(c("p1", "p2"), each = 4), 6)
  dat$GT_BIN <- rep(c(0, 0, 1, 1, 1, 2, 2, 2), 6)
  dat$BLOCK <- rep(rep(c("b1", "b2", "b3"), each = 2), each = 8)
  a <- amova_genomic(dat, "STRATA", value = "GT_BIN", standardized = FALSE,
                     resampling = "locus", bootstrap = 19, seed = 3)
  expect_equal(a$uncertainty$marker$replicates, 19)
  expect_equal(a$uncertainty$marker$unit, "locus")
  expect_named(a$uncertainty$marker$intervals,
               c("statistic", "estimate", "lower", "upper", "valid.replicates"))
  b <- amova_genomic(dat, "STRATA", value = "GT_BIN", standardized = FALSE,
                     resampling = "block", block = "BLOCK", bootstrap = 19, seed = 3)
  expect_equal(b$uncertainty$marker$unit, "genomic block")
})

test_that("population jackknife reports leave-one-population-out sensitivity", {
  dat <- expand.grid(INDIVIDUALS = paste0("i", 1:12), MARKERS = paste0("m", 1:3),
                     stringsAsFactors = FALSE)
  dat$STRATA <- rep(rep(c("p1", "p2", "p3"), each = 4), 3)
  dat$GT_BIN <- rep(c(0, 0, 1, 1, 0, 1, 1, 2, 1, 1, 2, 2), 3)
  z <- amova_genomic(dat, "STRATA", value = "GT_BIN", standardized = FALSE,
                     population.jackknife = TRUE)
  expect_setequal(z$uncertainty$population.jackknife$omitted, c("p1", "p2", "p3"))
  expect_true(all(z$uncertainty$population.jackknife$status == "estimable"))
})

test_that("two-level hierarchy is fitted and nesting is checked", {
  dat <- expand.grid(
    INDIVIDUALS = paste0("i", 1:16), MARKERS = c("m1", "m2"),
    stringsAsFactors = FALSE
  )
  dat$REGION <- rep(rep(c("r1", "r2"), each = 8), 2)
  dat$POP_ID <- rep(rep(c("p1", "p2", "p3", "p4"), each = 4), 2)
  dat$GT_BIN <- rep(c(0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 2), 2)
  z <- amova_genomic(dat, c("REGION", "POP_ID"), value = "GT_BIN", standardized = FALSE)
  expect_equal(names(z$variance_components), c("REGION", "POP_ID", "Within"))
  expect_equal(z$global$statistic, c("PHI_ST", "PHI_CT", "PHI_SC"))
})

test_that("identity distance supports standardized statistics", {
  dat <- expand.grid(
    INDIVIDUALS = paste0("i", 1:8), MARKERS = c("m1", "m2"),
    stringsAsFactors = FALSE
  )
  dat$POP_ID <- rep(rep(c("p1", "p2"), each = 4), 2)
  dat$HAP <- rep(c("A", "A", "B", "B", "C", "C", "D", "D"), 2)
  z <- amova_genomic(dat, "POP_ID", value = "HAP", distance = "identity")
  expect_true(any(is.finite(z$global$standardized)))
  expect_true(is.finite(z$global$standardized[z$global$statistic == "PHI_ST"]))
})

test_that("hierarchy-aware permutations return component tests", {
  dat <- expand.grid(
    INDIVIDUALS = paste0("i", 1:16), MARKERS = c("m1", "m2"),
    stringsAsFactors = FALSE
  )
  dat$REGION <- rep(rep(c("r1", "r2"), each = 8), 2)
  dat$POP_ID <- rep(rep(c("p1", "p2", "p3", "p4"), each = 4), 2)
  dat$GT_BIN <- rep(c(0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 2), 2)
  expect_warning(
    amova_genomic(dat, c("REGION", "POP_ID"), value = "GT_BIN", standardized = FALSE,
                  permutations = 9, seed = 11),
    "REGION"
  )
  z <- suppressWarnings(
    amova_genomic(dat, c("REGION", "POP_ID"), value = "GT_BIN", standardized = FALSE,
                  permutations = 9, seed = 11)
  )
  expect_equal(dim(z$permutation$components), c(9, 3))
  expect_true(all(z$permutation$p_value[1:2] >= 0.1))
  expect_true(is.na(z$permutation$p_value[3]))
  expect_equal(z$permutation$design$unique_allocations[1], 3)
  expect_equal(z$permutation$design$minimum_p[1], 1 / 3)
  expect_false(z$permutation$design$alpha_attainable[1])
})

test_that("one master hierarchy permutation is reusable across loci", {
  master <- data.frame(
    .individual = paste0("i", 1:8),
    REGION = rep(c("r1", "r2"), each = 4),
    POP_ID = rep(c("p1", "p2", "p3", "p4"), each = 2)
  )
  set.seed(4)
  permuted <- assigner:::.amova_permute_master(master, c("REGION", "POP_ID"), 1)
  locus_a <- c("i1", "i3", "i5", "i7")
  locus_b <- c("i1", "i2", "i5", "i6")
  expect_identical(
    permuted$REGION[match(intersect(locus_a, locus_b), permuted$.individual)],
    permuted$REGION[match(intersect(locus_b, locus_a), permuted$.individual)]
  )
  expect_equal(length(unique(permuted$REGION)), 2)
})

test_that("loci lacking a required hierarchy level are excluded", {
  dat <- expand.grid(
    INDIVIDUALS = paste0("i", 1:16), MARKERS = c("both", "one_region"),
    stringsAsFactors = FALSE
  )
  dat$REGION <- rep(rep(c("r1", "r2"), each = 8), 2)
  dat$POP_ID <- rep(rep(c("p1", "p2", "p3", "p4"), each = 4), 2)
  dat$GT_BIN <- rep(c(0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 2), 2)
  dat$GT_BIN[dat$MARKERS == "one_region" & dat$REGION == "r2"] <- NA
  z <- amova_genomic(dat, c("REGION", "POP_ID"), value = "GT_BIN", standardized = FALSE)
  expect_equal(z$per_locus$MARKERS, "both")
})

test_that("non-identity distances cannot claim a Meirmans maximum", {
  dat <- data.frame(INDIVIDUALS = paste0("i", 1:4), MARKERS = "m1",
                    POP_ID = rep(c("p1", "p2"), each = 2), GT_BIN = 0:3)
  expect_error(
    amova_genomic(dat, "POP_ID", value = "GT_BIN", standardized = TRUE),
    "requires.*identity"
  )
})

test_that("standard GDS-style numeric GT is selected automatically", {
  dat <- expand.grid(
    INDIVIDUALS = paste0("i", 1:8), MARKERS = c("m1", "m2"),
    stringsAsFactors = FALSE
  )
  dat$POP_ID <- rep(rep(c("p1", "p2"), each = 4), 2)
  dat$GT <- c(0, 0, 1, 1, 1, 2, 2, 2, 0, NA, 1, 1, 1, 2, 2, 2)
  z <- amova_genomic(dat, "POP_ID", standardized = FALSE)
  expect_equal(z$settings$value, ".AMOVA_DOSAGE")
  expect_equal(z$per_locus$N, c(8, 7))
})

test_that("long data-frame hierarchy is preserved without an import round trip", {
  dat <- expand.grid(
    INDIVIDUALS = paste0("i", 1:16), MARKERS = c("m1", "m2"),
    stringsAsFactors = FALSE
  )
  dat$REGION <- rep(rep(c("r1", "r2"), each = 8), 2)
  dat$POP_ID <- rep(rep(c("p1", "p2", "p3", "p4"), each = 4), 2)
  dat$GT <- rep(c(0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 2), 2)
  z <- amova_genomic(dat, c("REGION", "POP_ID"), standardized = FALSE)
  expect_equal(names(z$variance_components), c("REGION", "POP_ID", "Within"))
  expect_false(anyNA(z$per_locus$N))
})

test_that("six-digit calibrated GT is converted to biallelic dosage", {
  dat <- data.frame(
    INDIVIDUALS = paste0("i", 1:8), MARKERS = "m1",
    POP_ID = rep(c("p1", "p2"), each = 4),
    GT = c("001001", "001002", "001002", "002002",
           "001002", "002002", "002002", "000000")
  )
  z <- amova_genomic(dat, "POP_ID", standardized = FALSE)
  expect_equal(z$per_locus$N, 7)
})
