# Internal genotype-likelihood assignment helpers -----------------------------

normalize_genotype_likelihoods <- function(data, method) {
  columns <- paste0(method, c("_HOM_REF", "_HET", "_HOM_ALT"))
  missing.columns <- setdiff(columns, names(data))
  if (length(missing.columns)) {
    rlang::abort(paste0(
      "`genotype.method = \"", method, "\"` requires columns: ",
      paste(columns, collapse = ", "), ". Missing: ",
      paste(missing.columns, collapse = ", "), "."
    ))
  }

  x <- as.matrix(data[, columns, drop = FALSE])
  storage.mode(x) <- "double"
  if (method == "PL") {
    x <- 10^(-x / 10)
  } else {
    finite <- x[is.finite(x)]
    probability.scale <- length(finite) && all(finite >= 0 & finite <= 1)
    if (!probability.scale) x <- 10^x
  }

  sums <- rowSums(x, na.rm = TRUE)
  valid <- rowSums(is.finite(x)) == 3L & is.finite(sums) & sums > 0
  x[!valid, ] <- NA_real_
  x[valid, ] <- x[valid, , drop = FALSE] / sums[valid]

  data$GL0 <- x[, 1L]
  data$GL1 <- x[, 2L]
  data$GL2 <- x[, 3L]
  data
}

em_allele_frequency <- function(gl, tolerance = 1e-7, max.iterations = 200L) {
  valid <- stats::complete.cases(gl) & rowSums(gl) > 0
  gl <- gl[valid, , drop = FALSE]
  n <- nrow(gl)
  if (!n) return(list(p = NA_real_, n = 0L, iterations = 0L))

  p <- mean(gl[, 2L] + 2 * gl[, 3L]) / 2
  if (!is.finite(p)) p <- 0.5
  p <- pmin(1 - 1e-8, pmax(1e-8, p))

  for (iteration in seq_len(max.iterations)) {
    prior <- c((1 - p)^2, 2 * p * (1 - p), p^2)
    posterior <- sweep(gl, 2L, prior, `*`)
    denominator <- rowSums(posterior)
    ok <- is.finite(denominator) & denominator > 0
    if (!any(ok)) return(list(p = NA_real_, n = n, iterations = iteration))
    posterior[ok, ] <- posterior[ok, , drop = FALSE] / denominator[ok]
    p.new <- sum(posterior[ok, 2L] + 2 * posterior[ok, 3L]) / (2 * sum(ok))
    p.new <- pmin(1, pmax(0, p.new))
    if (abs(p.new - p) <= tolerance) {
      p <- p.new
      break
    }
    p <- p.new
  }
  list(p = p, n = n, iterations = iteration)
}

observed_information <- function(gl, p) {
  if (!is.finite(p)) return(rep(NA_real_, nrow(gl)))
  g0 <- gl[, 1L]
  g1 <- gl[, 2L]
  g2 <- gl[, 3L]
  u <- g0 * (1 - p)^2 + g1 * 2 * p * (1 - p) + g2 * p^2
  first <- 2 * p * (g0 + g2 - 2 * g1) + 2 * (g1 - g0)
  second <- 2 * (g0 + g2 - 2 * g1)
  information <- (first / u)^2 - second / u
  information[!is.finite(information)] <- NA_real_
  pmax(0, information)
}

expected_read_log_probability <- function(depth, p, sequencing.error) {
  if (!is.finite(depth) || depth < 0 || !is.finite(p)) {
    return(c(mean = NA_real_, variance = NA_real_))
  }
  depth <- as.integer(depth)
  hwe <- c((1 - p)^2, 2 * p * (1 - p), p^2)
  ref.probability <- c(1 - sequencing.error, 0.5, sequencing.error)
  values <- numeric()
  weights <- numeric()
  for (genotype in 0:2) {
    ref.count <- 0:depth
    count.probability <- stats::dbinom(
      ref.count, size = depth, prob = ref.probability[genotype + 1L]
    )
    for (j in seq_along(ref.count)) {
      r <- ref.count[j]
      a <- depth - r
      generated <- c(
        (1 - sequencing.error)^r * sequencing.error^a,
        0.5^depth,
        sequencing.error^r * (1 - sequencing.error)^a
      )
      if (!all(is.finite(generated)) || sum(generated) <= 0) next
      generated <- generated / sum(generated)
      values <- c(values, log(sum(generated * hwe)))
      weights <- c(weights, hwe[genotype + 1L] * count.probability[j])
    }
  }
  if (!length(values) || sum(weights) <= 0) {
    return(c(mean = NA_real_, variance = NA_real_))
  }
  weights <- weights / sum(weights)
  centre <- sum(weights * values)
  c(mean = centre, variance = sum(weights * (values - centre)^2))
}

raw_unsampled_z <- function(data, p, sequencing.error) {
  valid <- stats::complete.cases(
    data[, c("GL0", "GL1", "GL2", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH")]
  ) & is.finite(p)
  if (!any(valid)) return(NA_real_)
  data <- data[valid, , drop = FALSE]
  p <- p[valid]
  observed <- data$GL0 * (1 - p)^2 +
    data$GL1 * 2 * p * (1 - p) + data$GL2 * p^2
  observed <- log(observed)
  expected <- vapply(seq_len(nrow(data)), function(i) {
    expected_read_log_probability(
      depth = data$ALLELE_REF_DEPTH[i] + data$ALLELE_ALT_DEPTH[i],
      p = p[i], sequencing.error = sequencing.error
    )
  }, numeric(2))
  variance <- sum(expected["variance", ], na.rm = TRUE)
  if (!is.finite(variance) || variance <= 0) return(NA_real_)
  (sum(observed, na.rm = TRUE) - sum(expected["mean", ], na.rm = TRUE)) /
    sqrt(variance)
}

calculate_unsampled_scores <- function(
  assignment, genome, allele.frequencies, loo.frequencies,
  leave.one.out, sequencing.error, cutoff
) {
  required <- c("ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH")
  missing.columns <- setdiff(required, names(genome))
  if (length(missing.columns)) {
    rlang::abort(paste0(
      "The unsampled-source score requires REF/ALT read depths in `",
      paste(required, collapse = "` and `"), "`. Missing: ",
      paste(missing.columns, collapse = ", "), "."
    ))
  }

  score_one <- function(id, source, reference.score) {
    individual <- dplyr::filter(genome, INDIVIDUALS == id)
    frequencies <- allele.frequencies |>
      dplyr::filter(as.character(STRATA) == source) |>
      dplyr::select(MARKERS, ALT_FREQ)
    individual <- dplyr::left_join(individual, frequencies, by = "MARKERS")
    p <- individual$ALT_FREQ
    if (reference.score && leave.one.out) {
      loo <- loo.frequencies |>
        dplyr::filter(as.character(STRATA) == source, INDIVIDUALS == id) |>
        dplyr::select(MARKERS, LOO_ALT_FREQ)
      individual <- dplyr::left_join(individual, loo, by = "MARKERS")
      p <- individual$LOO_ALT_FREQ
    }
    raw_unsampled_z(individual, p, sequencing.error)
  }

  assignment$UNSAMPLED_Z_RAW <- vapply(seq_len(nrow(assignment)), function(i) {
    source <- if (is.na(assignment$CURRENT[i])) assignment$INFERRED[i] else
      as.character(assignment$CURRENT[i])
    score_one(
      id = assignment$INDIVIDUALS[i], source = source,
      reference.score = !is.na(assignment$CURRENT[i])
    )
  }, numeric(1))

  calibration <- assignment |>
    dplyr::filter(!is.na(CURRENT), is.finite(UNSAMPLED_Z_RAW)) |>
    dplyr::group_by(CURRENT) |>
    dplyr::summarise(
      Z_MEAN = mean(UNSAMPLED_Z_RAW),
      Z_SD = stats::sd(UNSAMPLED_Z_RAW),
      .groups = "drop"
    ) |>
    dplyr::rename(SOURCE = CURRENT)
  assignment$SOURCE_FOR_Z <- ifelse(
    is.na(assignment$CURRENT), assignment$INFERRED, as.character(assignment$CURRENT)
  )
  assignment <- assignment |>
    dplyr::left_join(calibration, by = c("SOURCE_FOR_Z" = "SOURCE")) |>
    dplyr::mutate(
      UNSAMPLED_Z = dplyr::if_else(
        is.finite(Z_SD) & Z_SD > 0,
        (UNSAMPLED_Z_RAW - Z_MEAN) / Z_SD,
        NA_real_
      ),
      UNSAMPLED_FLAG = dplyr::if_else(
        is.na(UNSAMPLED_Z), NA, UNSAMPLED_Z < cutoff
      )
    ) |>
    dplyr::select(-SOURCE_FOR_Z, -Z_MEAN, -Z_SD)
  list(assignment = assignment, calibration = calibration)
}

estimate_gl_frequency_group <- function(data, tolerance, max.iterations) {
  gl <- as.matrix(data[, c("GL0", "GL1", "GL2")])
  estimate <- em_allele_frequency(gl, tolerance, max.iterations)
  p <- estimate$p
  information <- observed_information(gl, p)
  effective.contribution <- 0.5 * information * p * (1 - p)
  effective.contribution <- pmin(1, pmax(0, effective.contribution))

  tibble::tibble(
    ALT_FREQ = p,
    N_GENOTYPED = estimate$n,
    EFFECTIVE_N = sum(effective.contribution, na.rm = TRUE),
    EM_ITERATIONS = estimate$iterations
  )
}

estimate_gl_loo_group <- function(data, tolerance, max.iterations) {
  gl <- as.matrix(data[, c("GL0", "GL1", "GL2")])
  n <- nrow(gl)
  if (n <= 1L) {
    return(tibble::tibble(
      INDIVIDUALS = data$INDIVIDUALS,
      LOO_ALT_FREQ = NA_real_,
      LOO_N = 0L,
      LOO_EFFECTIVE_N = NA_real_
    ))
  }
  estimates <- lapply(seq_len(n), function(i) {
    loo.gl <- gl[-i, , drop = FALSE]
    estimate <- em_allele_frequency(loo.gl, tolerance, max.iterations)
    information <- observed_information(loo.gl, estimate$p)
    contribution <- 0.5 * information * estimate$p * (1 - estimate$p)
    estimate$effective.n <- sum(pmin(1, pmax(0, contribution)), na.rm = TRUE)
    estimate
  })
  tibble::tibble(
    INDIVIDUALS = data$INDIVIDUALS,
    LOO_ALT_FREQ = vapply(estimates, `[[`, numeric(1), "p"),
    LOO_N = vapply(estimates, `[[`, integer(1), "n"),
    LOO_EFFECTIVE_N = vapply(estimates, `[[`, numeric(1), "effective.n")
  )
}

assign_individuals_gl <- function(
  genome, strata, genotype.method, leave.one.out, frequency.floor,
  effective.n.tolerance, em.tolerance, em.max.iterations,
  unsampled.source, unsampled.cutoff, sequencing.error,
  evaluation.data, marker.blocks, verbose
) {
  genome <- normalize_genotype_likelihoods(genome, genotype.method)

  if (!is.null(strata)) {
    strata.data <- genometranslator::read_strata(strata = strata)$strata
    genome <- dplyr::select(genome, -tidyselect::any_of("STRATA")) |>
      dplyr::left_join(
        dplyr::select(strata.data, INDIVIDUALS, STRATA), by = "INDIVIDUALS"
      )
  }
  if (!"STRATA" %in% names(genome)) {
    rlang::abort("Supply `strata`, or include a `STRATA` column in `data`.")
  }

  individuals <- dplyr::distinct(genome, INDIVIDUALS, STRATA)
  if (anyDuplicated(individuals$INDIVIDUALS)) {
    rlang::abort("Each individual must belong to at most one reference stratum.")
  }
  reference <- dplyr::filter(genome, !is.na(STRATA), nzchar(as.character(STRATA)))
  if (!nrow(reference)) rlang::abort("No reference individuals have a valid `STRATA` value.")

  scoring.genome <- genome
  if (!is.null(evaluation.data)) {
    evaluation <- genometranslator::read_genome(
      data = evaluation.data, import.metadata = TRUE
    ) |>
      normalize_genotype_likelihoods(genotype.method)
    keys <- c("INDIVIDUALS", "MARKERS")
    if (anyDuplicated(dplyr::select(evaluation, dplyr::all_of(keys)))) {
      rlang::abort("`evaluation.data` contains duplicated individual-marker rows.")
    }
    missing.keys <- dplyr::anti_join(
      dplyr::distinct(genome, INDIVIDUALS, MARKERS),
      dplyr::distinct(evaluation, INDIVIDUALS, MARKERS), by = keys
    )
    if (nrow(missing.keys)) {
      rlang::abort("`evaluation.data` must contain every individual-marker row in `data`.")
    }
    scoring.columns <- intersect(
      c("GL0", "GL1", "GL2", "ALLELE_REF_DEPTH", "ALLELE_ALT_DEPTH"),
      names(evaluation)
    )
    scoring.genome <- genome |>
      dplyr::select(-tidyselect::any_of(scoring.columns)) |>
      dplyr::left_join(
        dplyr::select(evaluation, dplyr::all_of(c(keys, scoring.columns))),
        by = keys
      )
  }

  if (isTRUE(verbose)) message("Estimating reference allele frequencies from ", genotype.method)
  allele.frequencies <- reference |>
    dplyr::select(STRATA, MARKERS, INDIVIDUALS, GL0, GL1, GL2) |>
    dplyr::group_by(STRATA, MARKERS) |>
    dplyr::group_modify(~ estimate_gl_frequency_group(
      .x, tolerance = em.tolerance, max.iterations = em.max.iterations
    )) |>
    dplyr::ungroup()

  candidate.strata <- sort(unique(as.character(allele.frequencies$STRATA)))
  common.markers <- allele.frequencies |>
    dplyr::group_by(MARKERS) |>
    dplyr::summarise(N = sum(is.finite(ALT_FREQ)), .groups = "drop") |>
    dplyr::filter(N == length(candidate.strata)) |>
    dplyr::pull(MARKERS)
  allele.frequencies <- dplyr::filter(allele.frequencies, MARKERS %in% common.markers)
  genome <- dplyr::filter(genome, MARKERS %in% common.markers)
  scoring.genome <- dplyr::filter(scoring.genome, MARKERS %in% common.markers)

  marker.blocks <- normalize_marker_blocks(marker.blocks, common.markers)

  effective.sample.size.markers <- dplyr::select(
    allele.frequencies, STRATA, MARKERS, N_GENOTYPED, ALT_FREQ, EFFECTIVE_N
  )

  effective.sample.size <- allele.frequencies |>
    dplyr::filter(ALT_FREQ > 0.05, ALT_FREQ < 0.95) |>
    dplyr::group_by(STRATA) |>
    dplyr::summarise(
      N_REFERENCE = dplyr::n_distinct(reference$INDIVIDUALS[reference$STRATA == dplyr::first(STRATA)]),
      EFFECTIVE_N = mean(EFFECTIVE_N, na.rm = TRUE),
      N_MARKERS = dplyr::n(),
      .groups = "drop"
    )

  individual.effective.sample.size <- reference |>
    dplyr::filter(MARKERS %in% common.markers) |>
    dplyr::left_join(
      dplyr::select(allele.frequencies, STRATA, MARKERS, ALT_FREQ),
      by = c("STRATA", "MARKERS")
    ) |>
    dplyr::group_by(STRATA, INDIVIDUALS) |>
    dplyr::group_modify(~ {
      gl <- as.matrix(.x[, c("GL0", "GL1", "GL2")])
      p <- .x$ALT_FREQ
      contribution <- vapply(seq_len(nrow(.x)), function(i) {
        info <- observed_information(gl[i, , drop = FALSE], p[i])
        pmin(1, pmax(0, 0.5 * info * p[i] * (1 - p[i])))
      }, numeric(1))
      tibble::tibble(
        EFFECTIVE_N = mean(contribution, na.rm = TRUE),
        N_MARKERS = sum(is.finite(contribution))
      )
    }) |>
    dplyr::ungroup()

  if (nrow(effective.sample.size) > 1L) {
    spread <- diff(range(effective.sample.size$EFFECTIVE_N, na.rm = TRUE)) /
      max(effective.sample.size$EFFECTIVE_N, na.rm = TRUE)
    if (is.finite(spread) && spread > effective.n.tolerance) {
      rlang::warn(paste0(
        "Reference-population effective sample sizes differ by ",
        round(100 * spread, 1), "%. Consider exploratory population ",
        "subsampling or coverage equalization before final assignment."
      ))
    }
  }

  loo.frequencies <- NULL
  if (leave.one.out) {
    if (isTRUE(verbose)) message("Estimating leave-one-out reference frequencies")
    loo.frequencies <- reference |>
      dplyr::filter(MARKERS %in% common.markers) |>
      dplyr::select(STRATA, MARKERS, INDIVIDUALS, GL0, GL1, GL2) |>
      dplyr::group_by(STRATA, MARKERS) |>
      dplyr::group_modify(~ estimate_gl_loo_group(
        .x, tolerance = em.tolerance, max.iterations = em.max.iterations
      )) |>
      dplyr::ungroup()
  }

  genotype.data <- scoring.genome |>
    dplyr::select(INDIVIDUALS, CURRENT = STRATA, MARKERS, GL0, GL1, GL2)

  score_candidate <- function(candidate) {
    frequencies <- allele.frequencies |>
      dplyr::filter(as.character(STRATA) == candidate) |>
      dplyr::select(MARKERS, ALT_FREQ, N_GENOTYPED, EFFECTIVE_N)
    scored <- dplyr::left_join(genotype.data, frequencies, by = "MARKERS")

    if (leave.one.out) {
      home.loo <- loo.frequencies |>
        dplyr::filter(as.character(STRATA) == candidate) |>
        dplyr::select(
          MARKERS, INDIVIDUALS, LOO_ALT_FREQ, LOO_N, LOO_EFFECTIVE_N
        )
      scored <- dplyr::left_join(scored, home.loo, by = c("MARKERS", "INDIVIDUALS"))
      use.loo <- !is.na(scored$CURRENT) & as.character(scored$CURRENT) == candidate
      scored$ALT_FREQ[use.loo] <- scored$LOO_ALT_FREQ[use.loo]
      scored$N_GENOTYPED[use.loo] <- scored$LOO_N[use.loo]
      scored$EFFECTIVE_N[use.loo] <- scored$LOO_EFFECTIVE_N[use.loo]
    }

    p <- scored$ALT_FREQ
    floor.value <- if (is.null(frequency.floor)) {
      effective.n <- ifelse(
        is.finite(scored$EFFECTIVE_N) & scored$EFFECTIVE_N > 0,
        scored$EFFECTIVE_N, scored$N_GENOTYPED
      )
      1 / (2 * effective.n + 1)
    } else rep(frequency.floor, length(p))
    p[p == 0] <- floor.value[p == 0]
    p[p == 1] <- 1 - floor.value[p == 1]
    probability <- scored$GL0 * (1 - p)^2 +
      scored$GL1 * 2 * p * (1 - p) + scored$GL2 * p^2
    probability[!is.finite(probability) | probability <= 0] <- NA_real_
    scored$LOG_PROBABILITY <- log(probability)
    summary <- scored |>
      dplyr::group_by(INDIVIDUALS, CURRENT) |>
      dplyr::summarise(
        LOG_LIKELIHOOD = if (all(is.na(LOG_PROBABILITY))) NA_real_ else
          sum(LOG_PROBABILITY, na.rm = TRUE),
        N_MARKERS = sum(!is.na(LOG_PROBABILITY)),
        .groups = "drop"
      ) |>
      dplyr::mutate(CANDIDATE = candidate)
    per.marker <- if (is.null(marker.blocks)) NULL else
      dplyr::transmute(
        scored, INDIVIDUALS, CURRENT, MARKERS,
        CANDIDATE = candidate, LOG_PROBABILITY
      )
    list(summary = summary, per.marker = per.marker)
  }

  scored.candidates <- purrr::map(candidate.strata, score_candidate)
  likelihoods <- purrr::map(scored.candidates, "summary") |>
    purrr::list_rbind() |>
    dplyr::select(INDIVIDUALS, CURRENT, CANDIDATE, LOG_LIKELIHOOD, N_MARKERS) |>
    dplyr::arrange(INDIVIDUALS, dplyr::desc(LOG_LIKELIHOOD), CANDIDATE)
  marker.block.likelihoods <- if (is.null(marker.blocks)) NULL else {
    purrr::map(scored.candidates, "per.marker") |>
      purrr::list_rbind() |>
      summarize_block_likelihoods(marker.blocks)
  }
  ranked <- likelihoods |>
    dplyr::group_by(INDIVIDUALS) |>
    dplyr::arrange(dplyr::desc(LOG_LIKELIHOOD), CANDIDATE, .by_group = TRUE) |>
    dplyr::mutate(RANK = dplyr::row_number()) |>
    dplyr::ungroup()
  best <- ranked |>
    dplyr::filter(RANK == 1L) |>
    dplyr::transmute(INDIVIDUALS, CURRENT, INFERRED = CANDIDATE,
      LOG_LIK_MAX = LOG_LIKELIHOOD, N_MARKERS)
  second <- ranked |>
    dplyr::filter(RANK == 2L) |>
    dplyr::transmute(INDIVIDUALS, SECOND_BEST = CANDIDATE,
      LOG_LIK_SECOND = LOG_LIKELIHOOD)
  home <- likelihoods |>
    dplyr::filter(!is.na(CURRENT), as.character(CURRENT) == CANDIDATE) |>
    dplyr::transmute(INDIVIDUALS, LOG_LIK_HOME = LOG_LIKELIHOOD)
  assignment <- best |>
    dplyr::left_join(second, by = "INDIVIDUALS") |>
    dplyr::left_join(home, by = "INDIVIDUALS") |>
    dplyr::mutate(
      LOG_LIK_RATIO = LOG_LIK_MAX - LOG_LIK_SECOND,
      HOME_LIK_RATIO = LOG_LIK_MAX - LOG_LIK_HOME,
      CORRECT = dplyr::if_else(is.na(CURRENT), NA,
        as.character(CURRENT) == INFERRED)
    ) |>
    dplyr::select(INDIVIDUALS, CURRENT, INFERRED, CORRECT,
      LOG_LIK_HOME, LOG_LIK_MAX, HOME_LIK_RATIO, SECOND_BEST,
      LOG_LIK_SECOND, LOG_LIK_RATIO, N_MARKERS) |>
    dplyr::arrange(CURRENT, INDIVIDUALS)

  unsampled.calibration <- NULL
  if (unsampled.source) {
    if (isTRUE(verbose)) message("Calculating standardized unsampled-source scores")
    unsampled <- calculate_unsampled_scores(
      assignment = assignment, genome = scoring.genome,
      allele.frequencies = allele.frequencies,
      loo.frequencies = loo.frequencies, leave.one.out = leave.one.out,
      sequencing.error = sequencing.error, cutoff = unsampled.cutoff
    )
    assignment <- unsampled$assignment
    unsampled.calibration <- unsampled$calibration
  }

  list(
    assignment = assignment,
    likelihoods = likelihoods,
    allele.frequencies = allele.frequencies,
    effective.sample.size = effective.sample.size,
    effective.sample.size.markers = effective.sample.size.markers,
    individual.effective.sample.size = individual.effective.sample.size,
    marker.block.likelihoods = marker.block.likelihoods,
    unsampled.calibration = unsampled.calibration,
    genotype.method = genotype.method
  )
}
