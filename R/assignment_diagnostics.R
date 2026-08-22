# Assignment validation and diagnostic helpers --------------------------------

normalize_marker_blocks <- function(marker.blocks, markers) {
  if (is.null(marker.blocks)) return(NULL)
  if (length(marker.blocks) == 1L && is.numeric(marker.blocks)) {
    n <- as.integer(marker.blocks)
    if (!is.finite(n) || n < 2L) {
      rlang::abort("A numeric `marker.blocks` value must be at least 2.")
    }
    ordered <- unique(as.character(markers))
    return(tibble::tibble(
      MARKERS = ordered,
      MARKER_BLOCK = paste0("block_", ((seq_along(ordered) - 1L) %% n) + 1L)
    ))
  }
  if (is.atomic(marker.blocks) && !is.null(names(marker.blocks))) {
    marker.blocks <- tibble::tibble(
      MARKERS = names(marker.blocks), MARKER_BLOCK = as.character(marker.blocks)
    )
  } else {
    marker.blocks <- tibble::as_tibble(marker.blocks)
    if (!all(c("MARKERS", "MARKER_BLOCK") %in% names(marker.blocks))) {
      rlang::abort(
        "`marker.blocks` must be a named vector, a MARKERS/MARKER_BLOCK table, or a block count."
      )
    }
    marker.blocks <- dplyr::transmute(
      marker.blocks, MARKERS = as.character(MARKERS),
      MARKER_BLOCK = as.character(MARKER_BLOCK)
    )
  }
  if (anyDuplicated(marker.blocks$MARKERS)) {
    rlang::abort("Each marker can belong to only one marker block.")
  }
  missing <- setdiff(unique(as.character(markers)), marker.blocks$MARKERS)
  if (length(missing)) {
    rlang::abort(paste0(
      "`marker.blocks` does not classify ", length(missing), " marker(s), including: ",
      paste(utils::head(missing, 5L), collapse = ", "), "."
    ))
  }
  marker.blocks
}

summarize_block_likelihoods <- function(per.marker, marker.blocks) {
  if (is.null(marker.blocks)) return(NULL)
  per.marker |>
    dplyr::left_join(marker.blocks, by = "MARKERS") |>
    dplyr::group_by(INDIVIDUALS, CURRENT, CANDIDATE, MARKER_BLOCK) |>
    dplyr::summarise(
      LOG_LIKELIHOOD = if (all(is.na(LOG_PROBABILITY))) NA_real_ else
        sum(LOG_PROBABILITY, na.rm = TRUE),
      N_MARKERS = sum(!is.na(LOG_PROBABILITY)),
      .groups = "drop"
    )
}

#' Compare assigner and WGSassign unsampled-source scores
#'
#' @description Join standardized unsampled-source scores produced by
#'   [assign_individuals()] with independently generated WGSassign scores. This
#'   is a validation aid; it does not call or reproduce WGSassign.
#'
#' @param result An object returned by [assign_individuals()] with
#'   `unsampled.source = TRUE`, or its `assignment` table.
#' @param wgsassign A data frame containing individual identifiers and WGSassign
#'   scores.
#' @param id.column Name of the WGSassign identifier column. Default:
#'   \code{id.column = "INDIVIDUALS"}.
#' @param score.column Name of the WGSassign score column. Default:
#'   \code{score.column = "WGSASSIGN_Z"}.
#'
#' @return A list containing the joined scores, Pearson and Spearman
#'   correlations, and a linear calibration model.
#' @export
compare_unsampled_scores <- function(
  result, wgsassign, id.column = "INDIVIDUALS", score.column = "WGSASSIGN_Z"
) {
  assignment <- if (is.list(result) && !is.null(result$assignment)) {
    result$assignment
  } else result
  if (!all(c("INDIVIDUALS", "UNSAMPLED_Z") %in% names(assignment))) {
    rlang::abort("`result` must contain `INDIVIDUALS` and `UNSAMPLED_Z`.")
  }
  wgsassign <- tibble::as_tibble(wgsassign)
  if (!all(c(id.column, score.column) %in% names(wgsassign))) {
    rlang::abort("The requested WGSassign identifier and score columns were not found.")
  }
  external <- dplyr::transmute(
    wgsassign,
    INDIVIDUALS = as.character(.data[[id.column]]),
    WGSASSIGN_Z = as.numeric(.data[[score.column]])
  )
  joined <- assignment |>
    dplyr::transmute(
      INDIVIDUALS = as.character(INDIVIDUALS),
      ASSIGNER_Z = UNSAMPLED_Z
    ) |>
    dplyr::inner_join(external, by = "INDIVIDUALS") |>
    dplyr::filter(is.finite(ASSIGNER_Z), is.finite(WGSASSIGN_Z))
  if (nrow(joined) < 3L) {
    rlang::abort("At least three individuals with finite scores are required.")
  }
  list(
    scores = joined,
    pearson = stats::cor(joined$ASSIGNER_Z, joined$WGSASSIGN_Z),
    spearman = stats::cor(joined$ASSIGNER_Z, joined$WGSASSIGN_Z, method = "spearman"),
    calibration = stats::lm(WGSASSIGN_Z ~ ASSIGNER_Z, data = joined)
  )
}

#' Monte Carlo assessment of an assignment reference baseline
#'
#' @description Repeatedly remove known reference individuals, treat them as an
#'   artificial mixture of unknown origin, and assign them against the remaining
#'   reference baseline. This evaluates assignment performance without allowing
#'   test individuals to contribute to reference allele frequencies.
#'
#' @param data Genomic data accepted by [assign_individuals()].
#' @param strata Optional strata metadata accepted by
#'   [genometranslator::read_strata()]. Default: \code{strata = NULL}.
#' @param repetitions Number of Monte Carlo replicates. Default:
#'   \code{repetitions = 50}.
#' @param mixture.size Total number of reference individuals held out in every
#'   replicate. Default: \code{mixture.size = 100}.
#' @param min.remaining Minimum individuals retained in every source population.
#'   Default: \code{min.remaining = 5}.
#' @param random.seed Optional simulation seed. Default:
#'   \code{random.seed = NULL}.
#' @param ... Additional arguments passed to [assign_individuals()].
#' @param verbose Logical. Display progress messages. Default:
#'   \code{verbose = TRUE}.
#'
#' @return A list with replicate-level assignments, accuracy summaries, and the
#'   identities held out in every replicate.
#'
#' @details Holdout numbers are allocated approximately in proportion to the
#' reference population sizes, while retaining at least `min.remaining`
#' individuals in every population. Marker selection and other data-dependent
#' preprocessing should be performed independently within each training
#' baseline when they form part of the inferential workflow.
#'
#' @references Moran BM, Anderson EC (2019). Bayesian inference from the
#' conditional genetic stock identification model. Canadian Journal of
#' Fisheries and Aquatic Sciences, 76(4), 551-560.
#' \doi{10.1139/cjfas-2018-0016}.
#' @export
assess_assignment_mc <- function(
  data, strata = NULL, repetitions = 50L, mixture.size = 100L,
  min.remaining = 5L, random.seed = NULL, ..., verbose = TRUE
) {
  repetitions <- as.integer(repetitions)
  mixture.size <- as.integer(mixture.size)
  min.remaining <- as.integer(min.remaining)
  if (repetitions < 1L || mixture.size < 1L || min.remaining < 1L) {
    rlang::abort("`repetitions`, `mixture.size`, and `min.remaining` must be positive integers.")
  }
  genome <- genometranslator::read_genome(data = data, import.metadata = TRUE)
  if (!is.null(strata)) {
    strata.data <- genometranslator::read_strata(strata = strata)$strata
    genome <- dplyr::select(genome, -tidyselect::any_of("STRATA")) |>
      dplyr::left_join(
        dplyr::select(strata.data, INDIVIDUALS, STRATA), by = "INDIVIDUALS"
      )
  }
  if (!"STRATA" %in% names(genome)) {
    rlang::abort("Supply `strata`, or include `STRATA` in `data`.")
  }
  samples <- genome |>
    dplyr::distinct(INDIVIDUALS, STRATA) |>
    dplyr::filter(!is.na(STRATA), nzchar(as.character(STRATA)))
  sizes <- dplyr::count(samples, STRATA, name = "N")
  capacity <- sum(pmax(0L, sizes$N - min.remaining))
  if (mixture.size > capacity) {
    rlang::abort(paste0(
      "`mixture.size` exceeds the available holdout capacity (", capacity,
      ") after retaining `min.remaining` individuals per population."
    ))
  }
  if (!is.null(random.seed)) set.seed(random.seed)

  draw_holdout <- function(iteration) {
    available <- pmax(0L, sizes$N - min.remaining)
    target <- mixture.size * available / sum(available)
    draw <- pmin(available, floor(target))
    remainder <- mixture.size - sum(draw)
    if (remainder > 0L) {
      priority <- order(target - floor(target), decreasing = TRUE)
      for (i in priority) {
        if (!remainder) break
        if (draw[i] < available[i]) {
          draw[i] <- draw[i] + 1L
          remainder <- remainder - 1L
        }
      }
    }
    purrr::map2(sizes$STRATA, draw, function(group, n) {
      ids <- dplyr::filter(samples, STRATA == group)$INDIVIDUALS
      if (!n) character() else sample(ids, n)
    }) |>
      unlist(use.names = FALSE)
  }

  holdouts <- lapply(seq_len(repetitions), draw_holdout)
  results <- lapply(seq_along(holdouts), function(i) {
    if (isTRUE(verbose)) message("Monte Carlo replicate ", i, " of ", repetitions)
    truth <- dplyr::filter(samples, INDIVIDUALS %in% holdouts[[i]]) |>
      dplyr::rename(TRUE_STRATA = STRATA)
    analysis <- genome |>
      dplyr::mutate(STRATA = dplyr::if_else(
        INDIVIDUALS %in% holdouts[[i]], NA_character_, as.character(STRATA)
      ))
    fit <- assign_individuals(analysis, ..., verbose = FALSE)
    fit$assignment |>
      dplyr::filter(INDIVIDUALS %in% holdouts[[i]]) |>
      dplyr::left_join(truth, by = "INDIVIDUALS") |>
      dplyr::mutate(
        REPLICATE = i,
        CORRECT = as.character(INFERRED) == as.character(TRUE_STRATA)
      )
  })
  assignment <- purrr::list_rbind(results)
  summary <- assignment |>
    dplyr::group_by(REPLICATE) |>
    dplyr::summarise(
      N = dplyr::n(), ACCURACY = mean(CORRECT, na.rm = TRUE), .groups = "drop"
    )
  by.stratum <- assignment |>
    dplyr::group_by(TRUE_STRATA) |>
    dplyr::summarise(
      N = dplyr::n(), ACCURACY = mean(CORRECT, na.rm = TRUE), .groups = "drop"
    )
  list(
    assignment = assignment,
    replicate.summary = summary,
    stratum.summary = by.stratum,
    holdouts = tibble::tibble(
      REPLICATE = rep(seq_along(holdouts), lengths(holdouts)),
      INDIVIDUALS = unlist(holdouts, use.names = FALSE)
    )
  )
}
