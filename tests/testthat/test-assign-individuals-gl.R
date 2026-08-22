make_likelihood_data <- function() {
  ids <- tibble::tibble(
    INDIVIDUALS = c(paste0("A", 1:6), paste0("B", 1:6)),
    STRATA = rep(c("A", "B"), each = 6)
  )
  tidyr::crossing(ids, MARKERS = paste0("m", 1:20)) |>
    dplyr::mutate(
      GT = dplyr::if_else(STRATA == "A", 0L, 2L),
      ID_NUMBER = as.integer(sub("^[AB]", "", INDIVIDUALS)),
      CONFIDENCE = 0.82 + 0.02 * ID_NUMBER,
      GL_HOM_REF = dplyr::if_else(GT == 0L, CONFIDENCE,
        (1 - CONFIDENCE) / 2),
      GL_HET = (1 - CONFIDENCE) / 2,
      GL_HOM_ALT = dplyr::if_else(GT == 2L, CONFIDENCE,
        (1 - CONFIDENCE) / 2),
      PL_HOM_REF = dplyr::if_else(GT == 0L, 0, 100),
      PL_HET = 100,
      PL_HOM_ALT = dplyr::if_else(GT == 2L, 0, 100),
      READ_DEPTH = 4L + ID_NUMBER,
      ALLELE_REF_DEPTH = dplyr::if_else(GT == 0L, READ_DEPTH, 0L),
      ALLELE_ALT_DEPTH = READ_DEPTH - ALLELE_REF_DEPTH
    )
}

test_that("certain GL and PL assign the expected reference population", {
  data <- make_likelihood_data()
  gl <- assign_individuals(data, genotype.method = "GL", verbose = FALSE)
  pl <- assign_individuals(data, genotype.method = "PL", verbose = FALSE)

  expect_true(all(gl$assignment$CORRECT))
  expect_true(all(pl$assignment$CORRECT))
  expect_equal(gl$assignment$INFERRED, pl$assignment$INFERRED)
  expect_true(all(gl$effective.sample.size$EFFECTIVE_N > 0))
})

test_that("unsampled-source score requires allele depth and is standardized", {
  data <- make_likelihood_data()
  no.depth <- dplyr::select(
    data, -ALLELE_REF_DEPTH, -ALLELE_ALT_DEPTH
  )
  expect_error(
    assign_individuals(
      no.depth, genotype.method = "GL", unsampled.source = TRUE,
      verbose = FALSE
    ),
    "requires REF/ALT read depths"
  )

  result <- assign_individuals(
    data, genotype.method = "GL", unsampled.source = TRUE,
    verbose = FALSE
  )
  expect_true(all(is.finite(result$assignment$UNSAMPLED_Z_RAW)))
  expect_true(all(is.finite(result$assignment$UNSAMPLED_Z)))
  expect_equal(
    result$assignment |>
      dplyr::group_by(CURRENT) |>
      dplyr::summarise(MEAN = mean(UNSAMPLED_Z), .groups = "drop") |>
      dplyr::pull(MEAN),
    c(0, 0), tolerance = 1e-8
  )
})

test_that("GL diagnostics retain marker, individual, and block information", {
  data <- make_likelihood_data()
  result <- assign_individuals(
    data, genotype.method = "GL", marker.blocks = 4L, verbose = FALSE
  )
  expect_true(all(c("STRATA", "MARKERS", "EFFECTIVE_N") %in%
    names(result$effective.sample.size.markers)))
  expect_true(all(c("STRATA", "INDIVIDUALS", "EFFECTIVE_N") %in%
    names(result$individual.effective.sample.size)))
  expect_equal(dplyr::n_distinct(result$marker.block.likelihoods$MARKER_BLOCK), 4L)
  expect_true(all(is.finite(result$marker.block.likelihoods$LOG_LIKELIHOOD)))
})

test_that("coverage-matched evaluation GLs do not change reference estimates", {
  data <- make_likelihood_data()
  lower.coverage <- data |>
    dplyr::mutate(
      GL_HOM_REF = 0.55 * (GT == 0L) + 0.225 * (GT != 0L),
      GL_HET = 0.225,
      GL_HOM_ALT = 0.55 * (GT == 2L) + 0.225 * (GT != 2L)
    )
  baseline <- assign_individuals(data, genotype.method = "GL", verbose = FALSE)
  matched <- assign_individuals(
    data, genotype.method = "GL", evaluation.data = lower.coverage,
    verbose = FALSE
  )
  expect_equal(matched$allele.frequencies, baseline$allele.frequencies)
  expect_false(isTRUE(all.equal(
    matched$likelihoods$LOG_LIKELIHOOD,
    baseline$likelihoods$LOG_LIKELIHOOD
  )))
})

test_that("reporting units aggregate collections without replacing them", {
  data <- make_likelihood_data() |>
    dplyr::mutate(REPUNIT = dplyr::if_else(STRATA == "A", "North", "South"))
  result <- assign_individuals(
    data, genotype.method = "GL", reporting.unit = "REPUNIT",
    verbose = FALSE
  )
  expect_true(is.list(result$reporting.unit))
  expect_setequal(result$reporting.unit$assignment$INFERRED, c("North", "South"))
})

test_that("external unsampled scores can be compared without importing WGSassign", {
  data <- make_likelihood_data()
  result <- assign_individuals(
    data, genotype.method = "GL", unsampled.source = TRUE, verbose = FALSE
  )
  external <- dplyr::transmute(
    result$assignment, INDIVIDUALS,
    WGSASSIGN_Z = 0.5 + 1.2 * UNSAMPLED_Z
  )
  comparison <- compare_unsampled_scores(result, external)
  expect_equal(comparison$pearson, 1, tolerance = 1e-10)
  expect_equal(nrow(comparison$scores), nrow(result$assignment))
})

test_that("Monte Carlo assessment keeps holdouts outside the baseline", {
  data <- make_likelihood_data()
  result <- assess_assignment_mc(
    data, repetitions = 2L, mixture.size = 2L, min.remaining = 4L,
    random.seed = 8L, genotype.method = "GT", verbose = FALSE
  )
  expect_equal(nrow(result$holdouts), 4L)
  expect_equal(nrow(result$replicate.summary), 2L)
  expect_true(all(result$assignment$CORRECT))
})
