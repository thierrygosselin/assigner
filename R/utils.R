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

# Exposition pipe-operator
#' @title Exposition pipe-operator
#' @description magrittr Exposition pipe-operator
#' @name %$%
#' @rdname Exposition_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %$%
#' @usage lhs \%$\% rhs
NULL

# compound assignment pipe operator
#' @title compound assignment pipe operator
#' @description magrittr compound assignment pipe operator
#' @name %<>%
#' @rdname compound_assignment_pipe_operator
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL


# subsampling_data --------------------------------------------------------------
#' @title subsampling data
#' @description subsampling data
#' @rdname subsampling_data
#' @export
#' @keywords internal
#' @importFrom dplyr mutate group_by ungroup arrange sample_n sample_frac


subsampling_data <- function(
  iteration.subsample = 1,
  ind.pop.df = NULL,
  subsample = NULL,
  random.seed = NULL
) {
  # message(paste0("Creating data subsample: ", iteration.subsample))
  if (is.null(subsample)) {
    subsample.select <- ind.pop.df %>% 
      dplyr::mutate(SUBSAMPLE = rep(iteration.subsample, n()))
  } else {
    
    # Set seed for sampling reproducibility
    if (is.null(random.seed)) {
      random.seed <- sample(x = 1:1000000, size = 1)
      set.seed(random.seed)
    } else {
      set.seed(random.seed)
    }
    
    if (subsample > 1) {# integer
      subsample.select <- ind.pop.df %>%
        dplyr::group_by(POP_ID) %>%
        dplyr::sample_n(tbl = ., size = subsample, replace = FALSE)# sampling individuals for each pop
    }

    subsample.select <- subsample.select %>% 
      dplyr::mutate(
        SUBSAMPLE = rep(iteration.subsample, n()),
        RANDOM_SEED_NUMBER = rep(random.seed, n())
      ) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS) %>% 
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
#' @importFrom dplyr bind_rows
#' @importFrom tibble as_data_frame
#' @importFrom stringi stri_join stri_detect_fixed stri_replace_all_fixed
#' @importFrom readr read_tsv

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
  data <- tibble::as_data_frame(dplyr::bind_rows(data))
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
#' @importFrom stringi stri_join
#' @importFrom dplyr bind_rows
#' @importFrom tibble as_data_frame
#' @importFrom readr read_tsv

#' @examples
#' \dontrun{
#' subsamples.data <- import_subsamples_fst(
#' dir.path = "assignment_analysis_method_ranked_no_imputations_20160425@2321"
#' )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

import_subsamples_fst <- function(dir.path){
  if (missing (dir.path)) stop("dir.path argument missing")
  
  # sampling.method <- stri_detect_fixed(str = dir.path, pattern = "ranked") # looks for ranked
  # if (sampling.method == FALSE) stop("This function doesn't work for markers sampled randomly")
  
  subsample.folders <- list.files(path = dir.path, pattern = "subsample_", full.names = FALSE)
  
  data.subsample <- list()
  for (i in subsample.folders) {
    fst.files.list <- list.files(path = stringi::stri_join(dir.path, "/", i), pattern = "fst.ranked", full.names = FALSE)
    data.fst <- list()
    for (j in fst.files.list) {
      fst.file <- readr::read_tsv(file = stringi::stri_join(dir.path, "/", i, "/", j), col_names = TRUE) %>% 
        mutate(
          SUBSAMPLE = rep(i, n()),
          ITERATIONS = rep(j, n())
        )
      data.fst[[j]] <- fst.file
    }
    data.fst <- tibble::as_data_frame(dplyr::bind_rows(data.fst))
    data.subsample[[i]] <- data.fst
  }
  data <- tibble::as_data_frame(dplyr::bind_rows(data.subsample))
  return(data)
}#End import_subsamples_fst


# GSI_BINARY_SUPPORT -----------------------------------------------------------


#' return the path where gsi_sim should be in the R system paths
#' 
#' @keywords internal
gsi_sim_binary_path <- function() {
  file.path(system.file(package = "assigner"), "bin", "gsi_sim")
}

#' return TRUE if gsi_sim exists where it should be
#' @keywords internal
gsi_sim_exists <- function() {
  file.exists(gsi_sim_binary_path())
}


#' return TRUE if gsi_sim is executable
#' @keywords internal
gsi_sim_is_executable <- function() {
  NULL #incomplete
}


#' file path to be used in a call to gsi_sim.
#' 
#' This version checks to make sure it is there and throws an
#' error with a suggestion of how to get it if it is not there.
#' @export
#' @keywords internal
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
# @keywords internal
#' @importFrom utils download.file
#' 
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
    message("")
    message("Removing any earlier instances of the repository in that temp directory")
    message("")
    system(paste("cd", td, "; rm -r -f gsi_sim"))
    message("Cloning repository, dealing with submodules, compiling gsi_sim ")
    message("")
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



