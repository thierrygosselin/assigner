# assigner_function_header -----------------------------------------------------
#' @title assigner_function_header
#' @description Generate function header
#' @rdname assigner_function_header
#' @keywords internal
#' @export
assigner_function_header <- function(f.name = NULL, start = TRUE, verbose = TRUE) {
  if (is.null(f.name)) invisible(NULL)
  if (start) {
    if (verbose) {
      message(stringi::stri_pad_both(str = "", width = 80L, pad = "#"))
      message(stringi::stri_pad_both(str = paste0(" assigner::", f.name, " "), width = 80L, pad = "#"))
      message(stringi::stri_pad_both(str = "", width = 80L, pad = "#"), "\n")
    }
  } else {
    if (verbose) {
      message(stringi::stri_pad_both(str = paste0(" completed ", f.name, " "), width = 80L, pad = "#"), "\n")
    }
  }
}# End assigner_function_header


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


# GSI_BINARY_SUPPORT -----------------------------------------------------------


#' return the path where gsi_sim should be in the R system paths
#'
#' @keywords internal
#' @name gsi_sim_binary_path
#' @rdname gsi_sim_binary_path
#' @export
gsi_sim_binary_path <- function() {
  file.path(system.file(package = "assigner"), "bin", "gsi_sim")
}

#' return TRUE if gsi_sim exists where it should be
#' @keywords internal
#' @export
#' @name gsi_sim_exists
#' @rdname gsi_sim_exists
gsi_sim_exists <- function() {
  file.exists(gsi_sim_binary_path())
}


#' return TRUE if gsi_sim is executable
#' @keywords internal
#' @export
#' @name gsi_sim_is_executable
#' @rdname gsi_sim_is_executable
gsi_sim_is_executable <- function() {
  NULL #incomplete
}


#' file path to be used in a call to gsi_sim.
#'
#' This version checks to make sure it is there and throws an
#' error with a suggestion of how to get it if it is not there.
#' @export
#' @keywords internal
#' @name gsi_sim_binary
#' @rdname gsi_sim_binary
gsi_sim_binary <- function() {
  if (!gsi_sim_exists()) stop("Can't find the gsi_sim executable where it was expected
                              at ", gsi_sim_binary_path(), ".
                              If you have internet access, you can install it
                              from within R by invoking the function \"install_gsi_sim()\"")

  # then I should check to make sure it is executable

  # if so, return the path
  gsi_sim_binary_path()

}


#' downloads gsi_sim that is appropriate for the operating system
#'
#' If the system is Mac, or Windows, this function will
#' download a precompiled binary from GitHub.  In other cases, or
#' if fromSource == TRUE, this function will attempt to download
#' the source code and compile the program from source and install
#' it.
#'
#' If this function fails, then you can just compile gsi_sim by
#' going to GITHUB_URL and compiling it yourself and naming the
#' executable gsi_sim and putting it at the location specified by the
#' function \code{\link{gsi_sim_binary_path}}.
#' @param commit  The full SHA-1 hash from GitHub from which to get
#' the binary or source
#' @param fromSource If TRUE, download source, even if a binary is available.
#' If FALSE, then it will download a precompiled binary, if available.  If a
#' binary is not available, then it will attempt to download the source.
#' @export
#' @name install_gsi_sim
#' @rdname install_gsi_sim
install_gsi_sim <- function(commit = "080f462c8eff035fa3e9f2fdce26c3ac013e208a", fromSource = FALSE) {

  # make a bin directory
  suppressWarnings(dir.create(file.path(system.file(package = "assigner"), "bin")))

  uname <- Sys.info()["sysname"]
  urlbase <- paste("https://github.com/eriqande/gsi_sim/blob/", commit,
                   "/gsi_sim-", sep = "")

  if (fromSource == FALSE) {
    if (uname == "Darwin") {
      url <- paste(urlbase, "Darwin", sep = "")
    }
    if (uname == "Windows") {
      url <- paste(urlbase, "MINGW32_NT-6.1", sep = "")
    }
    if (uname == "Darwin" || uname == "Windows") {
      message("Downloading file ", url)
      message("And copying to ", gsi_sim_binary_path())
      utils::download.file(url = url, destfile = gsi_sim_binary_path())
      Sys.chmod(gsi_sim_binary_path()) # make it executable
    }
    return(NULL)
  }

  if (uname == "Linux" || fromSource == TRUE) {  # in this case we will just compile from source
    td <- tempdir()

    message("Will be cloning gsi_sim repository to ", td)
    message("Removing any earlier instances of the repository in that temp directory")
    system(paste("cd", td, "; rm -r -f gsi_sim"))
    message("Cloning repository, dealing with submodules, compiling gsi_sim ")
    comm <- paste("cd ", td,
                  "&&  git clone https://github.com/eriqande/gsi_sim.git ",
                  "&&  cd gsi_sim ",
                  "&&  git checkout ", commit,
                  "&& git submodule init && git submodule update ",
                  "&& ./Compile_gsi_sim.sh ", sep = "")
    boing <- system(comm)
    if (boing != 0) {
      stop("Failed trying to clone and compile gsi_sim")
    } else {
      message("Apparently successful compiling gsi_sim.  Now copying to ", gsi_sim_binary_path())
      trycopy <- file.copy(from = paste(td, "/gsi_sim/gsi_sim-", uname, sep = ""),
                           to = gsi_sim_binary_path(),
                           overwrite = TRUE)
      if (trycopy == FALSE) stop("Apparently failed trying to copy ",
                                 paste(td, "/gsi_sim/gsi_sim-", uname, sep = ""),
                                 "to ",
                                 gsi_sim_binary_path())
    }
  }
}
# parallel_core_opt ------------------------------------------------------------
#' @title parallel_core_opt
#' @description Optimization of parallel core argument for radiator
#' @keywords internal
#' @export
parallel_core_opt <- function(parallel.core = NULL, max.core = NULL) {
  # strategy:
  # minimum of 1 core and a maximum of all the core available -2
  # even number of core
  # test
  # parallel.core <- 1
  # parallel.core <- 2
  # parallel.core <- 3
  # parallel.core <- 11
  # parallel.core <- 12
  # parallel.core <- 16
  # max.core <- 5
  # max.core <- 50
  # max.core <- NULL

  # Add-ons options
  # to control the max and min number to use...

  if (is.null(parallel.core)) {
    parallel.core <- parallel::detectCores() - 2
  } else {
    parallel.core <- floor(parallel.core / 2) * 2
    parallel.core <- max(1, min(parallel.core, parallel::detectCores() - 2))
  }

  if (is.null(max.core)) {
    parallel.core.opt <- parallel.core
  } else {
    parallel.core.opt <- min(parallel.core, floor(max.core / 2) * 2)
  }
  return(parallel.core.opt)
}#End parallel_core_opt

# assigner_future: future and future.apply -------------------------------------
#' @name assigner_future
#' @title assigner parallel function
#' @description Updating assigner to use future
# @inheritParams future::plan
# @inheritParams future::availableCores
#' @inheritParams future.apply::future_apply
#' @rdname assigner_future
#' @export
#' @keywords internal
assigner_future <- function(
  .x,
  .f,
  flat.future = c("int", "chr", "dfr", "dfc", "walk", "drop"),
  split.vec = FALSE,
  split.with = NULL,
  split.chunks = 4L,
  parallel.core = parallel::detectCores() - 1,
  forking = FALSE,
  ...
) {
  os <- Sys.info()[['sysname']]
  if (os == "Windows") forking <- FALSE

  opt.change <- getOption("width")
  options(width = 70)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(if (parallel.core > 1L && !forking) future::plan(strategy = "sequential"), add = TRUE)

  # argument for flattening the results
  flat.future <- match.arg(
    arg = flat.future,
    choices = c("int", "chr", "dfr", "dfc", "walk", "drop"),
    several.ok = FALSE
  )

  # splitting into chunks-------------------------------------------------------
  if (split.vec && is.null(split.with)) {
    # d: data, data length, data size
    # sv: split vector
    d <- .x
    df <- FALSE
    if (any(class(d) %in% c("tbl_df","tbl","data.frame"))) {
      d <- nrow(d)
      df <- TRUE
    }
    if (length(d) > 1L) d <- length(d)
    stopifnot(is.integer(d))
    sv <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
    # sv <- as.integer(floor((parallel.core * cpu.rounds * (seq_len(d) - 1) / d) + 1))
    stopifnot(length(sv) == d)

    # split
    if (df) {
      .x$SPLIT_VEC <- sv
      .x %<>% dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% split(x = ., f = sv)
    }
  }

  if (!is.null(split.with)) {
    # check
    if (length(split.with) != 1 || !is.character(split.with)) {
      rlang::abort(message = "Contact author: problem with parallel computation")
    }
    .data <- NULL
    stopifnot(rlang::has_name(.x, split.with))
    if (split.vec) {
      sv <- dplyr::distinct(.x, .data[[split.with]])
      d <- nrow(sv)
      sv$SPLIT_VEC <- as.integer(floor((split.chunks * (seq_len(d) - 1) / d) + 1))
      .x %<>%
        dplyr::left_join(sv, by = split.with) %>%
        dplyr::group_split(.tbl = ., "SPLIT_VEC", .keep = FALSE)
    } else {
      .x %<>% dplyr::group_split(.tbl = ., .data[[split.with]], .keep = TRUE)
    }
  }

  if (parallel.core == 1L) {
    future::plan(strategy = "sequential")
  } else {
    parallel.core <- parallel_core_opt(parallel.core = parallel.core)
    lx <- length(.x)
    if (lx < parallel.core) {
      future::plan(strategy = "multisession", workers = lx)
    } else {
      if (!forking) future::plan(strategy = "multisession", workers = parallel.core)
    }
  }

  # Run the function in parallel and account for dots-dots-dots argument
  if (forking) {
    if (length(list(...)) == 0) {
      rad_map <- switch(
        flat.future,
        int = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            vctrs::vec_c(.ptype = integer())
        },
        chr = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            vctrs::vec_c(.ptype = character())
        },
        dfr = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            dplyr::bind_rows(.)
        },
        dfc = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core) %>%
            dplyr::bind_cols(.)
        },
        walk = {furrr::future_walk},
        drop = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, mc.cores = parallel.core)
        }
      )
    } else {
      rad_map <- switch(
        flat.future,
        int = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            vctrs::vec_c(.ptype = integer())
        },
        chr = {.x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            vctrs::vec_c(.ptype = character())
        },
        dfr = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            dplyr::bind_rows(.)
        },
        dfc = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core) %>%
            dplyr::bind_cols(.)
        },
        walk = {furrr::future_walk},
        drop = {
          .x %<>%
            parallel::mclapply(X = ., FUN = .f, ..., mc.cores = parallel.core)
        }
      )
    }
  } else {
    rad_map <- switch(flat.future,
                      int = {furrr::future_map_int},
                      chr = {furrr::future_map_chr},
                      dfr = {furrr::future_map_dfr},
                      dfc = {furrr::future_map_dfc},
                      walk = {furrr::future_walk},
                      drop = {furrr::future_map}
    )
    # p <- NULL
    p <- progressr::progressor(along = .x)
    opts <- furrr::furrr_options(globals = TRUE, seed = TRUE)
    if (length(list(...)) == 0) {
      .x %<>% rad_map(.x = ., .f = .f, .options = opts)
    } else {
      .x %<>% rad_map(.x = ., .f = .f, ..., .options = opts)
    }
  }
  return(.x)
}#End assigner_future

# PIVOT-GATHER-CAST ------------------------------------------------------------
# rationale for doing this is that i'm tired of using tidyverse or data.table semantics
# tidyr changed from gather/spread to pivot_ functions but their are still very slow compared
# to 1. the original gather/spread and data.table equivalent...

#' @title rad_long
#' @description Gather, melt and pivot_longer
#' @rdname rad_long
#' @keywords internal
#' @export

rad_long <- function(
  x,
  cols = NULL,
  measure_vars = NULL,
  names_to = NULL,
  values_to = NULL,
  variable_factor = TRUE,
  keep_rownames = FALSE,
  tidy = FALSE
){


  # tidyr
  if (tidy) {
    x %>%
      tidyr::pivot_longer(
        data = .,
        cols = -cols,
        names_to = names_to,
        values_to = values_to
      )
  } else {# data.table
    x %>%
      data.table::as.data.table(., keep.rownames = keep_rownames) %>%
      data.table::melt.data.table(
        data = .,
        id.vars = cols,
        measure.vars = measure_vars,
        variable.name = names_to,
        value.name = values_to,
        variable.factor = variable_factor
      ) %>%
      tibble::as_tibble(.)
  }
}#rad_long

#' @title rad_wide
#' @description Spread, dcast and pivot_wider
#' @rdname rad_wide
#' @keywords internal
#' @export
rad_wide <- function(
  x ,
  formula = NULL,
  names_from = NULL,
  values_from = NULL,
  sep = "_",
  fun_aggregate = NULL,
  values_fill = NULL,
  tidy = FALSE
){
  # tidyr
  if (tidy) {
    x %<>%
      tidyr::pivot_wider(
        data = .,
        names_from = names_from,
        values_from = values_from,
        values_fill = values_fill
      )
  } else {# data.table
    if (is.null(fun_aggregate)) {
      x  %>%
        data.table::as.data.table(.) %>%
        data.table::dcast.data.table(
          data = .,
          formula =  formula,
          value.var = values_from,
          sep = sep,
          fill = values_fill
        ) %>%
        tibble::as_tibble(.)
    } else {
      x  %>%
        data.table::as.data.table(.) %>%
        data.table::dcast.data.table(
          data = .,
          formula =  formula,
          value.var = values_from,
          sep = sep,
          fun.aggregate = fun_aggregate,
          fill = values_fill
        ) %>%
        tibble::as_tibble(.)
    }
  }
}#rad_wide


# assigner_clock ---------------------------------------------------------------
#' @title assigner_tic
#' @description assigner tictoc function
#' @rdname assigner_tic
#' @keywords internal
#' @export
assigner_tic <- function(timing = proc.time()) {
  invisible(timing)
}# End assigner_tic

#' @title assigner_toc
#' @description assigner tictoc function
#' @rdname assigner_toc
#' @keywords internal
#' @export
assigner_toc <- function(
  timing,
  end.message = "Computation time, overall:",
  verbose = TRUE
) {
  if (verbose) message("\n", end.message, " ", round((proc.time() - timing)[[3]]), " sec")
}# End assigner_toc


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



#' @title change_matrix_strata
#' @description Integrate strata info back into the matrix
#' @rdname change_matrix_strata
#' @export
#' @keywords internal
change_matrix_strata <- function(x, s) {
  # x
  # class(colnames(x))
  # class(rownames(x))

  if (length(dim(x)) > 1) {
    if (identical(colnames(x), as.character(s$STRATA_SEQ))) {
      colnames(x) <- as.character(s$STRATA)
    } else {
      rlang::abort(message = "Contact author, problem with strata levels in fst")
    }
    if (identical(rownames(x), as.character(s$STRATA_SEQ))) {
      rownames(x) <- as.character(s$STRATA)
    } else {
      rlang::abort(message = "Contact author, problem with strata levels in fst")
    }
  }
  return(x)
}#change_matrix_strata

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

#' Format a numeric vector using formatC and flexible arguments
#'
#' @param x A numeric vector.
#' @param fmt.args A named list of arguments to pass to base::formatC.
#'
#' @return A character vector with formatted values.
#' @keywords internal
utils_formatC_vec <- function(x, fmt.args) {
  # rlang::exec() is not vectorised, so we use vapply + utils_formatC
  vapply(
    x,
    FUN = function(val) rlang::exec(base::formatC, x = val, !!!fmt.args),
    FUN.VALUE = character(1)
  )
}# END utils_formatC_vec

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




