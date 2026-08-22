# Pipe operator ----------------------------------------------------------------
#' @title Forward-pipe operator
#' @description magrittr forward-pipe operator
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

# Exposition pipe-operator -----------------------------------------------------
#' @title Exposition pipe-operator
#' @description magrittr Exposition pipe-operator
#' @name %$%
#' @rdname Exposition_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %$%
#' @usage lhs \%$\% rhs
NULL

# compound assignment pipe operator --------------------------------------------
#' @title compound assignment pipe operator
#' @description magrittr compound assignment pipe operator
#' @name %<>%
#' @rdname compound_assignment_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL

# Import used by serialized parallel worker functions.
#' @importFrom carrier crate
NULL

# dplyr n ----------------------------------------------------------------------
# The number of observations in the current group.
#' @title The number of observations in the current group.
#' @description Check dplyr
#' @name n
#' @rdname n
#' @keywords internal
#' @export
#' @importFrom dplyr n
#' @usage n()
NULL

# subsampling_data --------------------------------------------------------------
#' @title subsampling data
#' @description subsampling data
#' @rdname subsampling_data
#' @export
#' @keywords internal


subsampling_data <- function(
  iteration.subsample = 1,
  strata = NULL,
  subsample = NULL,
  random.seed = NULL
) {
  # message(paste0("Creating data subsample: ", iteration.subsample))
  if (is.null(subsample)) {
    subsample.select <- dplyr::mutate(strata, SUBSAMPLE = iteration.subsample)
  } else {

    # Set seed for sampling reproducibility
    if (is.null(random.seed)) {
      random.seed <- sample(x = 1:1000000, size = 1)
    }
    set.seed(random.seed)

    if (subsample > 1) {# integer
      subsample.select <- strata %>%
        dplyr::group_by(STRATA_SEQ) %>%
        dplyr::sample_n(tbl = ., size = subsample, replace = FALSE) %>%
        dplyr::ungroup(.)# sampling individuals for each pop
    }

    subsample.select %<>%
      dplyr::mutate(
        SUBSAMPLE = iteration.subsample,
        RANDOM_SEED_NUMBER = random.seed
      ) %>%
      dplyr::arrange(STRATA_SEQ, ID_SEQ) %>%
      dplyr::ungroup(.)
  }
  return(subsample.select)
} # End subsampling function



# import_subsamples ------------------------------------------------------------
# Write a dataframe containing all the subsample individual assignment

#' @name import_subsamples
#' @title Import individual's assignment results of different subsample folder.
#' @description This function will import all the individual's assignment results
#' of different subsample folder into R.
#' @param dir.path The path to the directory containing the subsample folders.
#' @param imputations (logical) Was the data imputed or not.
#' Default: \code{imputations = FALSE}
#' @return A data frame of all the individual's assignment, with iterations and subsample.

#' @export
#' @rdname import_subsamples

#' @examples
#' \dontrun{
#' subsamples.data <- import_subsamples(
#' dir.path = "assignment_analysis_method_random_imputations_rf_populations",
#' imputations = TRUE
#' )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

import_subsamples <- function(dir.path, imputations = FALSE){

  if (missing(dir.path)) stop("dir.path argument missing")

  sampling.method.files <- list.files(path = dir.path, pattern = "assignment", full.names = FALSE)[1]
  sampling.method <- stringi::stri_detect_fixed(str = sampling.method.files, pattern = "ranked") # looks for ranked

  subsample.folders <- list.files(path = dir.path, pattern = "subsample_", full.names = FALSE)
  data <- list()
  for (i in subsample.folders) {
    sub.name <- stringi::stri_replace_all_fixed(str = i, pattern = "_", replacement = ".", vectorize_all = FALSE)
    if (!sampling.method) {
      if (imputations) {
        filename <- stringi::stri_join(dir.path, "/", i, "/","assignment.random.imputed.results.individuals.iterations.", sub.name, ".tsv")
      } else {
        filename <- stringi::stri_join(dir.path, "/", i, "/","assignment.random.no.imputation.results.individuals.iterations.", sub.name, ".tsv")
      }
    } else {
      if (imputations) {
        filename <- stringi::stri_join(dir.path, "/", i, "/","assignment.ranked.imputed.results.individuals.iterations.", sub.name, ".tsv")
      } else {
        filename <- stringi::stri_join(dir.path, "/", i, "/","assignment.ranked.no.imputation.results.individuals.iterations.", sub.name, ".tsv")
      }
    }
    subsample.data <- readr::read_tsv(file = filename, col_names = TRUE)
    # mutate(SUBSAMPLE = rep(i, n()))
    # filter (MISSING_DATA == 'no.imputation')
    data[[i]] <- subsample.data
  }
  data <- tibble::as_tibble(dplyr::bind_rows(data))
  return(data)
}#End import_subsamples

# import_subsamples_fst---------------------------------------------------------
# Write a dataframe containing all the subsample individual assignment

#' @name import_subsamples_fst
#' @title Import the fst ranking from all the subsample runs inside
#' an assignment folder.
#' @description This function will import all the fst ranking from all the
#' subsample runs inside an assignment folder.
#' @param dir.path The path to the directory containing the subsample folders.
#' @return A data frame of all the Fst and ranking.

#' @export
#' @rdname import_subsamples_fst
#' @examples
#' \dontrun{
#' subsamples.data <- import_subsamples_fst(
#' dir.path = "assignment_analysis_method_ranked_no_imputations_20160425@2321"
#' )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

import_subsamples_fst <- function(dir.path){
  if (missing(dir.path)) stop("dir.path argument missing")

  # sampling.method <- stri_detect_fixed(str = dir.path, pattern = "ranked") # looks for ranked
  # if (sampling.method == FALSE) stop("This function doesn't work for markers sampled randomly")

  subsample.folders <- list.files(path = dir.path, pattern = "subsample_", full.names = FALSE)

  data.subsample <- list()
  for (i in subsample.folders) {
    fst.files.list <- list.files(path = stringi::stri_join(dir.path, "/", i), pattern = "fst.ranked", full.names = FALSE)
    data.fst <- list()
    for (j in fst.files.list) {
      fst.file <- readr::read_tsv(file = stringi::stri_join(dir.path, "/", i, "/", j), col_names = TRUE) %>%
        dplyr::mutate(
          SUBSAMPLE = rep(i, n()),
          ITERATIONS = rep(j, n())
        )
      data.fst[[j]] <- fst.file
    }
    data.fst <- tibble::as_tibble(dplyr::bind_rows(data.fst))
    data.subsample[[i]] <- data.fst
  }
  data <- tibble::as_tibble(dplyr::bind_rows(data.subsample))
  return(data)
}#End import_subsamples_fst


# gsi_sim executable -----------------------------------------------------------

#' Check the gsi_sim installation
#'
#' `gsi_sim`, developed by assigner co-author Eric C. Anderson, is an external
#' command-line program and is not installed by the R package. This function
#' searches the active `PATH` and validates that the executable can be started.
#'
#' @references Anderson EC (2010). Assessing the power of informative subsets
#'   of loci for population assignment: standard methods are upwardly biased.
#'   Molecular Ecology Resources, 10, 701-710.
#'
#' @param error Logical. Stop with installation guidance when `gsi_sim` cannot
#'   be used. Default: \code{error = TRUE}.
#' @param verbose Logical. Display the detected executable. Default:
#'   \code{verbose = TRUE}.
#'
#' @return Invisibly returns the executable path, or an empty string when it is
#'   unavailable and \code{error = FALSE}.
#' @export
check_gsi_sim <- function(error = TRUE, verbose = TRUE) {
  executable <- Sys.which("gsi_sim")
  ok <- nzchar(executable) && file.access(executable, mode = 1L) == 0L

  if (!ok) {
    guidance <- paste0(
      "gsi_sim was not found on PATH. It is an external program and is not ",
      "installed by assigner. Follow the source-installation instructions in ",
      "the assigner get-started vignette, then restart R and run ",
      "assigner::check_gsi_sim().\nSource: https://github.com/eriqande/gsi_sim"
    )
    if (isTRUE(error)) rlang::abort(guidance)
    if (isTRUE(verbose)) message(guidance)
    return(invisible(""))
  }

  status <- suppressWarnings(system2(executable, "--help", stdout = FALSE, stderr = FALSE))
  if (!identical(status, 0L)) {
    guidance <- paste0("gsi_sim was found at ", executable,
      " but could not be started successfully.")
    if (isTRUE(error)) rlang::abort(guidance)
    if (isTRUE(verbose)) message(guidance)
    return(invisible(""))
  }

  if (isTRUE(verbose)) message("gsi_sim: ", executable)
  invisible(unname(executable))
}

#' Locate the gsi_sim executable
#' @return The validated path to `gsi_sim`.
#' @keywords internal
gsi_sim_binary <- function() check_gsi_sim(error = TRUE, verbose = FALSE)




# Internal Matrix Utilities for Pairwise Distance Data -------------------------

#' @title Internal Matrix Utilities for Pairwise Distance Data
#'
#' @description
#' These internal utility functions convert long-format pairwise population
#' data into symmetric matrices, typically used for visualising pairwise
#' FST or genetic distance estimates.
#'
#' @details
#' `make_upper_matrix()` converts long-form pairwise data into a wide upper-triangular matrix.
#' `make_full_symmetric_matrix()` takes that upper matrix and completes the symmetric form.
#'
#' @param data A data frame or tibble with `POP1`, `POP2`, and a value column.
#' @param value.col The unquoted column to use as matrix values.
#' @param fill.value The value to fill missing cells (default: `""`).
#' @param upper.matrix A square matrix (upper triangle filled).
#' @param diagonal.value The value to assign to the diagonal (default: `"0"`).
#'
#' @return A matrix. The first function returns an upper-triangular matrix, and
#' the second a full symmetric matrix with diagonal.
#'
#' @keywords internal
#' @rdname pairwise_matrix_utils
#' @export
#'
# Helper: make upper matrix from long format
make_upper_matrix <- function(data, value.col, fill.value = "") {
  mat <- data %>%
    tidyr::complete(POP1, POP2) %>%
    tidyr::pivot_wider(
      names_from = POP2,
      values_from = {{ value.col }},
      values_fill = fill.value
    ) %>%
    dplyr::rename(POP = POP1)

  rn <- mat$POP
  mat <- as.matrix(mat[, -1])
  rownames(mat) <- rn
  return(mat)
}#End make_upper_matrix

#' @rdname pairwise_matrix_utils
#' @export
# Helper: fill lower triangle and diagonal of matrix
make_full_symmetric_matrix <- function(upper.matrix, diagonal.value = "0") {
  full <- upper.matrix
  full[lower.tri(full)] <- t(full)[lower.tri(full)]
  diag(full) <- diagonal.value
  return(full)
}#End make_full_symmetric_matrix



# match_markers_meta -----------------------------------------------------------
#' @title match_markers_meta
#' @description Integrate markers meta info back into the data
#' @rdname match_markers_meta
#' @export
#' @keywords internal
match_markers_meta <- function(x, markers.meta) {
  x  %<>%
    dplyr::left_join(markers.meta, by = "M_SEQ") %>%
    dplyr::select(-M_SEQ)
  return(x)
}# End match_strata

# fst_write -----------------------------------------------------------
#' @title fst_write
#' @description Write the fst results to working directory or path
#' @rdname fst_write
#' @export
#' @keywords internal
fst_write <- function(x,list.sub, path.folder) {
  readr::write_tsv(
    x = list.sub[[x]],
    file = file.path(path.folder, paste0(x, ".tsv"))
  )
}#End fst_write


# Format mean and range of numeric vector as string ----------------------------
#' @title Format mean and range of numeric vector as string
#' @description Produces a string like "0.123 [0.100 - 0.150]" for a numeric vector.
#' @param x A numeric vector.
#' @param meanx Logical; whether to include the mean value in the output. Defaults to `TRUE`.
#' @param min.max Logical; whether to include the minimum and maximum values in the output. Defaults to `TRUE`.
#' @param digits Number of digits *after the decimal point*.
#' @param scientific Logical; whether to use scientific notation.
#' @return A character string with formatted mean, min, and max values.
#' @keywords internal
format_mean_range <- function(x, meanx = TRUE, min.max = TRUE, digits = 3, scientific = FALSE) {

  # arg for formatC
  fmt.args <- list(
    format = "f",
    digits = digits,
    flag = "",
    big.mark = "",
    decimal.mark = ".",
    drop0trailing = FALSE
  )
  if (scientific) fmt.args$format <- "e"

  # default:
  sep1 <- " ["
  sep2 <- " - "
  sep3 <- "]"

  if (meanx) {
    mean_str <- do.call(formatC, c(list(x = mean(x, na.rm = TRUE)), fmt.args))
  } else {
    mean_str <- ""
    sep1 <- "["
  }


  if (min.max) {
    min_str  <- do.call(formatC, c(list(x = pmax(0, min(x, na.rm = TRUE))),  fmt.args))
    max_str  <- do.call(formatC, c(list(x = max(x, na.rm = TRUE)),  fmt.args))
  } else {
    min_str <- ""
    max_str <- ""
    sep1 <- ""
    sep2 <- ""
    sep3 <- ""
  }
  # Concatenate and return the formatted string
  stringi::stri_c(mean_str, sep1, min_str, sep2, max_str, sep3)
}# END format_mean_range

#' Format numeric using formatC and dynamic arguments
#'
#' @param x A numeric scalar.
#' @param fmt.args A named list of arguments to pass to base::formatC.
#'
#' @return A character string formatted by formatC.
#' @keywords internal

utils_formatC <- function(x, fmt.args) {
  rlang::exec(base::formatC, x = x, !!!fmt.args)
}#END utils_formatC

# Plot the distribution of FST values-------------------------------------------
#' @title Plot the distribution of FST values
#' @description Creates a histogram of FST values from a tibble or data frame.
#'
#' @param data A tibble or data frame containing a numeric `FST` column.
#' @param binwidth A numeric value specifying the histogram bin width. Default is `0.01`.
#'
#' @return A `ggplot2` object showing the histogram of FST values.
#' @keywords internal
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs expand_limits theme element_text
plot_fst_distribution <- function(data, binwidth = 0.01) {
  bold.text <- ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")

  ggplot2::ggplot(data, ggplot2::aes(x = FST, na.rm = TRUE)) +
    ggplot2::geom_histogram(binwidth = binwidth) +
    ggplot2::labs(x = "Fst (overall)") +
    ggplot2::expand_limits(x = 0) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = bold.text,
      axis.title.y = bold.text,
      legend.title = bold.text,
      legend.text = bold.text,
      strip.text.x = bold.text
    )
}# END plot_fst_distribution


# Extract arguments matching --------

#' @title Extract arguments matching a function's formals
#' @rdname extract_matching_args
#' @description
#' Extracts named arguments from an environment (e.g., the calling environment)
#' that match the formal arguments of a target function.
#'
#' Useful when forwarding arguments from a parent function to a subfunction
#' without manually repeating each parameter.
#'
#' @param from.env The environment from which to extract arguments (typically `environment()` or `rlang::caller_env()`).
#' @param to.fn The target function whose formal arguments will be matched.
#' @param .evaluate Logical, whether to force evaluation (default: TRUE).
#' @param .exclude Optional character vector of argument names to exclude.

#' @return A named list of matched and evaluated arguments (excluding `.exclude`).
#'
#' @examples
#' \dontrun{
#' args.for.fst <- extract_matching_args(from.env = environment(), to.fn = assigner::compute_fst)
#' result <- rlang::exec(assigner::compute_fst, x = my_data, !!!args.for.fst)
#' }
#'
#' @keywords internal
#' @export
#'
#' @importFrom rlang env_names env_get set_names
#' @importFrom purrr map
extract_matching_args <- function(from.env, to.fn, .evaluate = TRUE, .exclude = NULL) {
  fn.args <- names(formals(to.fn))
  env.args <- rlang::env_names(from.env)

  matched.args <- env.args[env.args %in% fn.args]

  if (!is.null(.exclude)) {
    matched.args <- setdiff(matched.args, .exclude)
  }

  # Extract the matched arguments and evaluate if requested
  out <- rlang::set_names(matched.args) %>%
    purrr::map(~ rlang::env_get(from.env, .x))

  if (.evaluate) {
    # Force evaluation so they're not lazy expressions/closures
    out <- purrr::map(out, identity)
  }

  out  # Return the evaluated arguments
}
