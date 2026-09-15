test_that("Dlr averages directional log-likelihood ratios", {
  likelihoods <- tibble::tibble(
    INDIVIDUALS = c("A1", "A2", "B1"),
    STRATA = c("A", "A", "B"),
    A = c(5, 7, 1),
    B = c(3, 3, 7)
  )

  # A has ratios 2 and 4; B has ratio 6. Equal population weighting gives 4.5.
  expect_equal(dlr_relative(likelihoods, "A", "B"), 4.5)
})

test_that("Dlr requires finite support from both populations", {
  likelihoods <- tibble::tibble(
    INDIVIDUALS = c("A1", "B1"),
    STRATA = c("A", "B"),
    A = c(2, NA_real_),
    B = c(1, NA_real_)
  )

  expect_error(
    dlr_relative(likelihoods, "A", "B"),
    "at least one finite likelihood ratio per population"
  )
})

test_that("dlr consumes native likelihood tables", {
  likelihoods <- tidyr::crossing(
    INDIVIDUALS = c("A1", "A2", "B1"),
    CANDIDATE = c("A", "B")
  ) |>
    dplyr::mutate(
      CURRENT = dplyr::if_else(grepl("^A", .data$INDIVIDUALS), "A", "B"),
      LOG_LIKELIHOOD = dplyr::case_when(
        .data$INDIVIDUALS == "A1" & .data$CANDIDATE == "A" ~ 5,
        .data$INDIVIDUALS == "A1" ~ 3,
        .data$INDIVIDUALS == "A2" & .data$CANDIDATE == "A" ~ 7,
        .data$INDIVIDUALS == "A2" ~ 3,
        .data$INDIVIDUALS == "B1" & .data$CANDIDATE == "A" ~ 1,
        TRUE ~ 7
      )
    )

  result <- suppressMessages(dlr(likelihoods, verbose = FALSE))
  expect_equal(result$dlr.table$DLR, 4.5)
  expect_setequal(result$assignment$STRATA, c("A", "B"))
})

test_that("Dlr inherits integrated sequencing genotype likelihoods", {
  genomic <- tidyr::crossing(
    INDIVIDUALS = paste0(rep(c("A", "B"), each = 4), rep(1:4, 2)),
    MARKERS = paste0("m", 1:6)
  ) |>
    dplyr::mutate(
      STRATA = substr(.data$INDIVIDUALS, 1, 1),
      GL_HOM_REF = dplyr::if_else(.data$STRATA == "A", 0.90, 0.05),
      GL_HET = 0.05,
      GL_HOM_ALT = dplyr::if_else(.data$STRATA == "B", 0.90, 0.05)
    )
  assignment <- assign_individuals(
    genomic,
    genotype.method = "GL",
    verbose = FALSE
  )
  result <- suppressMessages(dlr(assignment, verbose = FALSE))

  wide <- assignment$likelihoods |>
    dplyr::select("INDIVIDUALS", STRATA = "CURRENT",
                  "CANDIDATE", "LOG_LIKELIHOOD") |>
    tidyr::pivot_wider(names_from = "CANDIDATE",
                       values_from = "LOG_LIKELIHOOD")
  expected <- wide |>
    dplyr::mutate(
      RATIO = dplyr::if_else(.data$STRATA == "A",
                             .data$A - .data$B, .data$B - .data$A)
    ) |>
    dplyr::summarise(MEAN = mean(.data$RATIO), .by = "STRATA") |>
    dplyr::summarise(DLR = mean(.data$MEAN)) |>
    dplyr::pull(.data$DLR)

  expect_equal(result$dlr.table$DLR, round(expected, 2))
})
