#' Assign individuals to reference strata
#'
#' @description
#' Estimate allele frequencies in reference strata and calculate the genotype
#' log-likelihood of every individual under every candidate stratum. The native
#' likelihood engine is implemented entirely in R and returns the complete
#' individual-by-stratum score table required for assignment diagnostics and
#' Dlr calculations.
#'
#' @param data A GDS file or object, or tidy genomic data accepted by
#'   [genometranslator::read_genome()]. Genotypes are converted internally to
#'   diploid alternate-allele dosage (`0`, `1`, `2`, or `NA`).
#' @param strata Optional strata file or object accepted by
#'   [genometranslator::read_strata()]. When `NULL`, `STRATA` must be present in
#'   `data`. Default: \code{strata = NULL}.
#' @param engine Assignment engine: native population-genetic likelihood,
#'   elastic-net multinomial regression, XGBoost, or experimental TabPFN.
#'   Default: \code{engine = "likelihood"}.
#' @param leave.one.out Logical. When assigning a reference individual, exclude
#'   that individual from allele-frequency estimation in its home stratum.
#'   Default: \code{leave.one.out = TRUE}.
#' @param frequency.floor Minimum allele frequency used to prevent zero
#'   likelihoods. With `NULL`, the floor is estimated independently for every
#'   stratum and marker as `1 / (2 * N + 1)`, where `N` is the number of
#'   genotyped reference individuals. Default:
#'   \code{frequency.floor = NULL}.
#' @param folds Number of stratified outer validation folds used by machine-
#'   learning engines. Default: \code{folds = 5}.
#' @param random.seed Integer controlling fold creation and stochastic engines.
#'   Default: \code{random.seed = NULL}.
#' @param verbose Logical. Display progress messages. Default:
#'   \code{verbose = TRUE}.
#'
#' @return A list with:
#' \describe{
#'   \item{assignment}{One row per individual with its current stratum,
#'     inferred stratum, home likelihood, maximum likelihood, likelihood-ratio
#'     statistic, and second-best stratum.}
#'   \item{likelihoods}{The complete table containing the log-likelihood of
#'     every individual under every candidate stratum.}
#'   \item{allele.frequencies}{Reference-stratum allele-frequency estimates
#'     before individual leave-one-out adjustment.}
#' }
#'
#' @details
#' For genotype dosage `0`, `1`, and `2`, the Hardy-Weinberg genotype
#' probabilities are `(1-p)^2`, `2p(1-p)`, and `p^2`, respectively, where `p`
#' is the alternate-allele frequency in a candidate stratum. Probabilities are
#' accumulated on the log scale. Missing genotypes do not contribute to the
#' score.
#'
#' The complete `likelihoods` table is intentionally retained. The best and
#' second-best assignments alone are insufficient for Dlr calculations and can
#' conceal assignment ambiguity among additional candidate strata.
#'
#' @section Missing data:
#' Missingness can encode sequencing batch, library, coverage, or processing
#' differences and can therefore create misleading assignment accuracy. The
#' likelihood engine ignores missing focal genotypes, estimates each
#' stratum-marker frequency from called reference genotypes only, and compares
#' candidate strata using the same called-marker set for an individual.
#' Elastic net uses training-fold mean imputation. XGBoost uses its native
#' missing-value routing. TabPFN uses its native preprocessing, including a
#' missingness indicator. All preprocessing and model fitting occur within the
#' training portion of each outer fold. Inspect the returned missingness table
#' and test whether missingness is associated with strata or processing batch.
#'
#' This function performs assignment, not general genomic filtering. Prepare
#' and quality-control the input with
#' \href{https://thierrygosselin.github.io/genometranslator/}{genometranslator}
#' and \href{https://thierrygosselin.github.io/radr/}{radr} first.
#'
#' @examples
#' \dontrun{
#' result <- assigner::assign_individuals(
#'   data = genome,
#'   leave.one.out = TRUE
#' )
#'
#' result$assignment
#' result$likelihoods
#' }
#'
#' @export
assign_individuals <- function(
  data,
  strata = NULL,
  engine = c("likelihood", "elastic_net", "xgboost", "tabpfn"),
  leave.one.out = TRUE,
  frequency.floor = NULL,
  folds = 5L,
  random.seed = NULL,
  verbose = TRUE
) {
  .start <- tgbase::startup(
    f.name = "assign_individuals", package = "assigner", verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  if (missing(data)) rlang::abort("`data` is required.")
  engine <- match.arg(engine)
  if (engine != "likelihood") {
    return(assign_individuals_ml(
      data = data, strata = strata, engine = engine, folds = folds,
      random.seed = random.seed, verbose = verbose
    ))
  }
  if (!is.logical(leave.one.out) || length(leave.one.out) != 1L || is.na(leave.one.out)) {
    rlang::abort("`leave.one.out` must be TRUE or FALSE.")
  }
  if (!is.null(frequency.floor) &&
      (!is.numeric(frequency.floor) || length(frequency.floor) != 1L ||
       is.na(frequency.floor) || frequency.floor <= 0 || frequency.floor >= 0.5)) {
    rlang::abort("`frequency.floor` must be NULL or one number between 0 and 0.5.")
  }

  genome <- genometranslator::read_genome(data = data, import.metadata = TRUE)
  required <- c("INDIVIDUALS", "MARKERS")
  missing.columns <- setdiff(required, names(genome))
  if (length(missing.columns)) {
    rlang::abort(paste0(
      "The genomic data is missing required columns: ",
      paste(missing.columns, collapse = ", "), "."
    ))
  }

  if ("ALT_DOSAGE" %in% names(genome)) {
    genome <- dplyr::mutate(genome, GT = as.integer(ALT_DOSAGE))
  } else if ("GT" %in% names(genome) &&
      all(is.na(genome$GT) | as.character(genome$GT) %in% c("0", "1", "2"))) {
    genome <- dplyr::mutate(genome, GT = as.integer(GT))
  } else if ("GT" %in% names(genome) &&
      all(is.na(genome$GT) | grepl("^[0-9]{6}$", as.character(genome$GT)))) {
    gt.character <- as.character(genome$GT)
    genome$A1 <- substr(gt.character, 1L, 3L)
    genome$A2 <- substr(gt.character, 4L, 6L)
    genome$A1[gt.character == "000000" | is.na(gt.character)] <- NA_character_
    genome$A2[gt.character == "000000" | is.na(gt.character)] <- NA_character_
    marker.alleles <- genome |>
      dplyr::select(MARKERS, A1, A2) |>
      tidyr::pivot_longer(c(A1, A2), values_to = "ALLELE") |>
      dplyr::filter(!is.na(ALLELE)) |>
      dplyr::group_by(MARKERS) |>
      dplyr::summarise(
        REF_ALLELE = sort(unique(ALLELE))[1L],
        N_ALLELES = dplyr::n_distinct(ALLELE),
        .groups = "drop"
      )
    if (any(marker.alleles$N_ALLELES > 2L)) {
      rlang::abort("The native likelihood engine currently requires biallelic markers.")
    }
    genome <- genome |>
      dplyr::left_join(marker.alleles, by = "MARKERS") |>
      dplyr::mutate(
        GT = dplyr::if_else(
          is.na(A1) | is.na(A2), NA_integer_,
          as.integer(A1 != REF_ALLELE) + as.integer(A2 != REF_ALLELE)
        )
      ) |>
      dplyr::select(-A1, -A2, -REF_ALLELE, -N_ALLELES)
  } else {
    rlang::abort(paste0(
      "The native likelihood engine requires `ALT_DOSAGE`, numeric dosage in ",
      "`GT`, or six-digit biallelic genotypes in `GT`."
    ))
  }

  if (!is.null(strata)) {
    strata.data <- genometranslator::read_strata(strata = strata) $strata
    genome <- dplyr::select(genome, -tidyselect::any_of("STRATA")) |>
      dplyr::left_join(
        dplyr::select(strata.data, INDIVIDUALS, STRATA),
        by = "INDIVIDUALS"
      )
  }
  if (!"STRATA" %in% names(genome)) {
    rlang::abort("Supply `strata`, or include a `STRATA` column in `data`.")
  }
  if (any(!is.na(genome$GT) & !genome$GT %in% 0:2)) {
    rlang::abort("`GT` must contain diploid alternate-allele dosages: 0, 1, 2, or NA.")
  }

  individuals <- dplyr::distinct(genome, INDIVIDUALS, STRATA)
  duplicated.ids <- individuals |>
    dplyr::count(INDIVIDUALS, name = "N") |>
    dplyr::filter(N > 1L)
  if (nrow(duplicated.ids)) {
    rlang::abort("Each individual must belong to at most one reference stratum.")
  }

  reference <- dplyr::filter(genome, !is.na(STRATA), nzchar(as.character(STRATA)))
  if (!nrow(reference)) rlang::abort("No reference individuals have a valid `STRATA` value.")

  allele.frequencies <- reference |>
    dplyr::group_by(STRATA, MARKERS) |>
    dplyr::summarise(
      ALT_COUNT = sum(GT, na.rm = TRUE),
      N_GENOTYPED = sum(!is.na(GT)),
      ALT_FREQ = dplyr::if_else(N_GENOTYPED > 0L, ALT_COUNT / (2 * N_GENOTYPED), NA_real_),
      .groups = "drop"
    )

  candidate.strata <- sort(unique(as.character(allele.frequencies$STRATA)))
  common.frequency.markers <- allele.frequencies |>
    dplyr::group_by(MARKERS) |>
    dplyr::summarise(
      N_CANDIDATES = sum(N_GENOTYPED > 0L),
      .groups = "drop"
    ) |>
    dplyr::filter(N_CANDIDATES == length(candidate.strata)) |>
    dplyr::pull(MARKERS)

  genotype.data <- dplyr::select(
    genome, INDIVIDUALS, CURRENT = STRATA, MARKERS, GT
  ) |>
    dplyr::filter(MARKERS %in% common.frequency.markers)

  if (leave.one.out) {
    home.counts <- allele.frequencies |>
      dplyr::select(CURRENT = STRATA, MARKERS, HOME_N_GENOTYPED = N_GENOTYPED)
    genotype.data <- genotype.data |>
      dplyr::left_join(home.counts, by = c("CURRENT", "MARKERS")) |>
      dplyr::filter(
        is.na(CURRENT) | is.na(GT) | is.na(HOME_N_GENOTYPED) |
          HOME_N_GENOTYPED > 1L
      ) |>
      dplyr::select(-HOME_N_GENOTYPED)
  }

  score_candidate <- function(candidate) {
    frequencies <- dplyr::filter(allele.frequencies, as.character(STRATA) == candidate) |>
      dplyr::select(MARKERS, ALT_COUNT, N_GENOTYPED, ALT_FREQ)

    scored <- dplyr::left_join(genotype.data, frequencies, by = "MARKERS")
    use.loo <- leave.one.out & !is.na(scored$CURRENT) &
      as.character(scored$CURRENT) == candidate & !is.na(scored$GT)
    alt.count <- scored$ALT_COUNT - ifelse(use.loo, scored$GT, 0)
    n.genotyped <- scored$N_GENOTYPED - ifelse(use.loo, 1L, 0L)
    p <- ifelse(n.genotyped > 0L, alt.count / (2 * n.genotyped), NA_real_)
    floor.value <- if (is.null(frequency.floor)) {
      ifelse(n.genotyped > 0L, 1 / (2 * n.genotyped + 1), NA_real_)
    } else {
      rep(frequency.floor, length(p))
    }
    p <- pmin(1 - floor.value, pmax(floor.value, p))

    log.probability <- rep(NA_real_, length(p))
    valid <- !is.na(scored$GT) & !is.na(p)
    log.probability[valid & scored$GT == 0] <- 2 * log1p(-p[valid & scored$GT == 0])
    log.probability[valid & scored$GT == 1] <- log(2) + log(p[valid & scored$GT == 1]) +
      log1p(-p[valid & scored$GT == 1])
    log.probability[valid & scored$GT == 2] <- 2 * log(p[valid & scored$GT == 2])

    scored$LOG_PROBABILITY <- log.probability
    scored |>
      dplyr::group_by(INDIVIDUALS, CURRENT) |>
      dplyr::summarise(
        LOG_LIKELIHOOD = if (all(is.na(LOG_PROBABILITY))) NA_real_ else sum(LOG_PROBABILITY, na.rm = TRUE),
        N_MARKERS = sum(!is.na(LOG_PROBABILITY)),
        .groups = "drop"
      ) |>
      dplyr::mutate(CANDIDATE = candidate)
  }

  if (isTRUE(verbose)) {
    message("Calculating likelihoods for ", length(candidate.strata), " candidate strata")
  }
  likelihoods <- purrr::map(candidate.strata, score_candidate) |>
    purrr::list_rbind() |>
    dplyr::select(INDIVIDUALS, CURRENT, CANDIDATE, LOG_LIKELIHOOD, N_MARKERS) |>
    dplyr::arrange(INDIVIDUALS, dplyr::desc(LOG_LIKELIHOOD), CANDIDATE)

  ranked <- likelihoods |>
    dplyr::group_by(INDIVIDUALS) |>
    dplyr::arrange(dplyr::desc(LOG_LIKELIHOOD), CANDIDATE, .by_group = TRUE) |>
    dplyr::mutate(RANK = dplyr::row_number()) |>
    dplyr::ungroup()

  best <- ranked |>
    dplyr::filter(RANK == 1L) |>
    dplyr::transmute(
      INDIVIDUALS, CURRENT, INFERRED = CANDIDATE,
      LOG_LIK_MAX = LOG_LIKELIHOOD, N_MARKERS
    )
  second <- ranked |>
    dplyr::filter(RANK == 2L) |>
    dplyr::transmute(
      INDIVIDUALS, SECOND_BEST = CANDIDATE,
      LOG_LIK_SECOND = LOG_LIKELIHOOD
    )
  home <- likelihoods |>
    dplyr::filter(!is.na(CURRENT), as.character(CURRENT) == CANDIDATE) |>
    dplyr::transmute(INDIVIDUALS, LOG_LIK_HOME = LOG_LIKELIHOOD)

  assignment <- best |>
    dplyr::left_join(second, by = "INDIVIDUALS") |>
    dplyr::left_join(home, by = "INDIVIDUALS") |>
    dplyr::mutate(
      LOG_LIK_RATIO = LOG_LIK_MAX - LOG_LIK_SECOND,
      HOME_LIK_RATIO = LOG_LIK_MAX - LOG_LIK_HOME,
      CORRECT = dplyr::if_else(
        is.na(CURRENT), NA, as.character(CURRENT) == INFERRED
      )
    ) |>
    dplyr::select(
      INDIVIDUALS, CURRENT, INFERRED, CORRECT,
      LOG_LIK_HOME, LOG_LIK_MAX, HOME_LIK_RATIO,
      SECOND_BEST, LOG_LIK_SECOND, LOG_LIK_RATIO, N_MARKERS
    ) |>
    dplyr::arrange(CURRENT, INDIVIDUALS)

  list(
    assignment = assignment,
    likelihoods = likelihoods,
    allele.frequencies = allele.frequencies
  )
}
