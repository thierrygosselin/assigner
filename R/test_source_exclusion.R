#' Monte Carlo exclusion test for candidate source populations
#'
#' Test whether an individual's multilocus genotype is unusually improbable
#' under each candidate source population. This complements the relative ranking
#' returned by [assign_individuals()]: an individual can rank first for a source
#' while still being incompatible with every sampled source.
#'
#' @param result Result returned by [assign_individuals()] using the native
#'   called-genotype likelihood engine.
#' @param data Optional genomic data used for the assignment. Supplying it lets
#'   simulations reproduce the exact non-missing marker set of each individual.
#'   Without it, the function matches only the individual's number of scored
#'   markers by sampling from markers available in the candidate source.
#'   Default: \code{data = NULL}.
#' @param simulations Number of simulated genotypes per candidate source.
#'   Default: \code{simulations = 10000}.
#' @param alpha Lower-tail probability used for exclusion. Default:
#'   \code{alpha = 0.01}.
#' @param frequency.floor Minimum allele frequency used to avoid zero
#'   probabilities. With `NULL`, use `1 / (2 * N + 1)` from the reference count
#'   for each source-marker combination. Default:
#'   \code{frequency.floor = NULL}.
#' @param random.seed Optional simulation seed. Default:
#'   \code{random.seed = NULL}.
#' @param verbose Logical. Display progress messages. Default:
#'   \code{verbose = TRUE}.
#'
#' @return A list containing one row per individual and candidate source with
#' the observed likelihood, Monte Carlo lower-tail probability, exclusion call,
#' and simulation threshold; an individual summary; and simulation settings.
#'
#' @details For each source and marker, genotypes are generated under
#' Hardy-Weinberg proportions using the reference alternate-allele frequency.
#' The same likelihood model used for assignment scores each simulated
#' genotype. The Monte Carlo probability is `(b + 1) / (B + 1)`, where `b` is
#' the number of simulated likelihoods no greater than the observed likelihood
#' and `B` is `simulations`.
#'
#' This is a source-exclusion diagnostic, not proof of origin or migration.
#' Results depend on representative reference sampling, allele-frequency
#' estimates, Hardy-Weinberg assumptions, genotyping error, missingness, linkage
#' disequilibrium, family structure, and whether the true source was sampled.
#' Linked markers make the product likelihood overconfident. Repeat the test on
#' linkage-pruned data or biologically defined marker blocks. Candidate
#' inversion regions should be analysed explicitly because suppressed
#' recombination can contribute many correlated markers.
#'
#' The current simulator is for called diploid genotypes. GL/PL assignment needs
#' coverage-aware read or genotype-likelihood simulation and is not silently
#' approximated here.
#'
#' @references Cornuet JM, Piry S, Luikart G, Estoup A, Solignac M (1999). New
#' methods employing multilocus genotypes to select or exclude populations as
#' origins of individuals. Genetics, 153, 1989-2000.
#'
#' Paetkau D, Slade R, Burden M, Estoup A (2004). Genetic assignment methods for
#' the direct, real-time estimation of migration rate: a simulation-based
#' exploration of accuracy and power. Molecular Ecology, 13, 55-65.
#' \doi{10.1046/j.1365-294X.2004.02008.x}.
#'
#' Manel S, Gaggiotti OE, Waples RS (2005). Assignment methods: matching
#' biological questions with appropriate techniques. Trends in Ecology &
#' Evolution, 20, 136-142. \doi{10.1016/j.tree.2004.12.004}.
#'
#' @export
test_source_exclusion <- function(
    result,
    data = NULL,
    simulations = 10000L,
    alpha = 0.01,
    frequency.floor = NULL,
    random.seed = NULL,
    verbose = TRUE
) {
  if (!is.list(result) ||
      !all(c("likelihoods", "allele.frequencies") %in% names(result))) {
    rlang::abort("`result` must be returned by `assign_individuals()`.")
  }
  simulations <- as.integer(simulations)
  if (length(simulations) != 1L || is.na(simulations) || simulations < 100L) {
    rlang::abort("`simulations` must be a whole number of at least 100.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    rlang::abort("`alpha` must be strictly between 0 and 1.")
  }
  method <- if (is.null(result$settings$genotype.method)) "GT" else
    result$settings$genotype.method
  engine <- if (is.null(result$settings$engine)) "likelihood" else
    result$settings$engine
  if (!identical(method, "GT") || !identical(engine, "likelihood")) {
    rlang::abort("`test_source_exclusion()` currently requires the native GT likelihood engine.")
  }
  if (!is.null(random.seed)) set.seed(random.seed)

  observed.markers <- NULL
  if (!is.null(data)) {
    genome <- genometranslator::read_genome(data = data, import.metadata = TRUE)
    dosage.name <- intersect(c("ALT_DOSAGE", "GT"), names(genome))[1L]
    if (is.na(dosage.name)) {
      rlang::abort("`data` must contain `ALT_DOSAGE` or numeric `GT`.")
    }
    observed.markers <- genome |>
      dplyr::filter(!is.na(.data[[dosage.name]])) |>
      dplyr::distinct(INDIVIDUALS, MARKERS)
  }

  frequencies <- dplyr::mutate(
    result$allele.frequencies,
    STRATA = as.character(STRATA),
    FLOOR = if (is.null(frequency.floor)) {
      1 / (2 * N_GENOTYPED + 1)
    } else frequency.floor,
    P = pmin(1 - .data$FLOOR, pmax(.data$FLOOR, .data$ALT_FREQ))
  ) |>
    dplyr::filter(is.finite(P), N_GENOTYPED > 0L)
  likelihoods <- dplyr::filter(
    result$likelihoods,
    is.finite(LOG_LIKELIHOOD), N_MARKERS > 0L
  )

  if (verbose) {
    message(
      "Running ", simulations, " Monte Carlo genotypes for each individual-source test..."
    )
  }
  tests <- purrr::pmap_dfr(
    likelihoods[c("INDIVIDUALS", "CURRENT", "CANDIDATE", "LOG_LIKELIHOOD", "N_MARKERS")],
    function(INDIVIDUALS, CURRENT, CANDIDATE, LOG_LIKELIHOOD, N_MARKERS) {
      focal.id <- INDIVIDUALS
      p.table <- dplyr::filter(frequencies, STRATA == CANDIDATE)
      if (!is.null(observed.markers)) {
        ids <- dplyr::filter(
          observed.markers, .data$INDIVIDUALS == focal.id
        )$MARKERS
        p.table <- dplyr::filter(p.table, MARKERS %in% ids)
      } else if (nrow(p.table) > N_MARKERS) {
        p.table <- dplyr::slice_sample(p.table, n = N_MARKERS)
      }
      p <- p.table$P
      if (!length(p)) return(tibble::tibble())
      q0 <- (1 - p)^2
      q1 <- 2 * p * (1 - p)
      u <- matrix(stats::runif(simulations * length(p)), nrow = simulations)
      log0 <- 2 * log1p(-p)
      log1 <- log(2) + log(p) + log1p(-p)
      log2 <- 2 * log(p)
      simulated <- rowSums(
        ifelse(u < rep(q0, each = simulations), rep(log0, each = simulations),
          ifelse(
            u < rep(q0 + q1, each = simulations),
            rep(log1, each = simulations), rep(log2, each = simulations)
          )
        )
      )
      probability <- (sum(simulated <= LOG_LIKELIHOOD) + 1) / (simulations + 1)
      tibble::tibble(
        INDIVIDUALS = INDIVIDUALS,
        CURRENT = CURRENT,
        CANDIDATE = CANDIDATE,
        LOG_LIKELIHOOD = LOG_LIKELIHOOD,
        N_MARKERS = length(p),
        EXCLUSION_THRESHOLD = as.numeric(stats::quantile(
          simulated, probs = alpha, names = FALSE, type = 8
        )),
        MONTE_CARLO_P = probability,
        EXCLUDED = probability < alpha
      )
    }
  )
  summary <- tests |>
    dplyr::group_by(INDIVIDUALS, CURRENT) |>
    dplyr::summarise(
      N_SOURCES_TESTED = dplyr::n(),
      N_SOURCES_EXCLUDED = sum(.data$EXCLUDED),
      ALL_SOURCES_EXCLUDED = all(.data$EXCLUDED),
      SOURCES_NOT_EXCLUDED = paste(
        .data$CANDIDATE[!.data$EXCLUDED], collapse = ";"
      ),
      .groups = "drop"
    )
  structure(list(
    exclusion = tests,
    individual.summary = summary,
    settings = list(
      simulations = simulations, alpha = alpha,
      frequency.floor = frequency.floor, exact.missingness.pattern = !is.null(data),
      random.seed = random.seed
    )
  ), class = c("source_exclusion", "list"))
}

#' @export
print.source_exclusion <- function(x, ...) {
  cat("Monte Carlo source-exclusion test\n")
  cat("  Individuals:", nrow(x$individual.summary), "\n")
  cat("  Individuals excluded from every source:",
      sum(x$individual.summary$ALL_SOURCES_EXCLUDED), "\n")
  invisible(x)
}
