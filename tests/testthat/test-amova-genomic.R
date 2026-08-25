test_that("locuswise missingness uses the called samples at each marker", {
  dat <- expand.grid(
    INDIVIDUALS = paste0("i", 1:8),
    MARKERS = c("m1", "m2"),
    stringsAsFactors = FALSE
  )
  dat$POP_ID <- rep(rep(c("p1", "p2"), each = 4), 2)
  dat$GT_BIN <- c(0, 0, 1, 1, 1, 2, 2, 2, 0, NA, 1, 1, 1, 2, 2, 2)

  z <- amova_genomic(dat, hierarchy = "POP_ID", value = "GT_BIN", missing = "locuswise",
                     min_individuals = 2, standardized = FALSE)
  expect_s3_class(z, "assigner_amova")
  expect_equal(z$per_locus$N, c(8, 7))
  expect_true(all(is.finite(z$variance_components)))
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
  z <- amova_genomic(dat, c("REGION", "POP_ID"), value = "GT_BIN", standardized = FALSE,
                     permutations = 9, seed = 11)
  expect_equal(dim(z$permutation$components), c(9, 3))
  expect_true(all(z$permutation$p_value[1:2] >= 0.1))
  expect_true(is.na(z$permutation$p_value[3]))
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
