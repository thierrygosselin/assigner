test_that("Monte Carlo source exclusion distinguishes relative assignment from fit", {
  frequencies <- tidyr::expand_grid(
    STRATA = c("A", "B"), MARKERS = paste0("m", seq_len(40L))
  ) |>
    dplyr::mutate(
      N_GENOTYPED = 20L,
      ALT_FREQ = dplyr::if_else(STRATA == "A", 0.1, 0.9),
      ALT_COUNT = ALT_FREQ * 2 * N_GENOTYPED
    )
  likelihoods <- tidyr::expand_grid(
    INDIVIDUALS = c("a1", "unknown"), CANDIDATE = c("A", "B")
  ) |>
    dplyr::mutate(
      CURRENT = c("A", "A", NA_character_, NA_character_),
      LOG_LIKELIHOOD = c(-15, -120, -200, -210),
      N_MARKERS = 40L
    )
  assignment <- list(
    likelihoods = likelihoods,
    allele.frequencies = frequencies,
    settings = list(engine = "likelihood", genotype.method = "GT")
  )
  result <- test_source_exclusion(
    assignment, simulations = 500L, alpha = 0.01,
    random.seed = 14L, verbose = FALSE
  )
  expect_s3_class(result, "source_exclusion")
  expect_equal(nrow(result$exclusion), 4L)
  expect_true(result$individual.summary$ALL_SOURCES_EXCLUDED[
    result$individual.summary$INDIVIDUALS == "unknown"
  ])
  expect_false(result$settings$exact.missingness.pattern)
})

test_that("source exclusion rejects GL and ML results", {
  result <- list(
    likelihoods = tibble::tibble(), allele.frequencies = tibble::tibble(),
    settings = list(engine = "likelihood", genotype.method = "GL")
  )
  expect_error(
    test_source_exclusion(result, simulations = 100L),
    "native GT likelihood engine"
  )
})
