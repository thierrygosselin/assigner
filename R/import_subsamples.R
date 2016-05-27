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
#' @import dplyr
#' @import stringi


#' @examples
#' \dontrun{
#' subsamples.data <- import_subsamples(
#' dir.path = "assignment_analysis_method_random_imputations_rf_populations",
#' imputations = TRUE
#' )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

import_subsamples <- function(dir.path, imputations){
  if (missing (dir.path)) stop("dir.path argument missing")
  if (missing (imputations)) imputations <- FALSE
  
  sampling.method.files <- list.files(path = dir.path, pattern = "assignment", full.names = FALSE)[1]
  sampling.method <- stri_detect_fixed(str = sampling.method.files, pattern = "ranked") # looks for ranked
  
  subsample.folders <- list.files(path = dir.path, pattern = "subsample_", full.names = FALSE)
  data <- list()
  for (i in subsample.folders) {
    sub.name <- stri_replace_all_fixed(str = i, pattern = "_", replacement = ".", vectorize_all = FALSE)
    if (sampling.method == FALSE){
      if (imputations == TRUE){
        filename <- stri_paste(dir.path, "/", i, "/","assignment.random.imputed.results.individuals.iterations.", sub.name, ".tsv")
      } else {
        filename <- stri_paste(dir.path, "/", i, "/","assignment.random.no.imputation.results.individuals.iterations.", sub.name, ".tsv")
      }
    } else {
      if (imputations == TRUE){
        filename <- stri_paste(dir.path, "/", i, "/","assignment.ranked.imputed.results.individuals.iterations.", sub.name, ".tsv")
      } else {
        filename <- stri_paste(dir.path, "/", i, "/","assignment.ranked.no.imputation.results.individuals.iterations.", sub.name, ".tsv")
      }
    }
    subsample.data <- read_tsv(file = filename, col_names = TRUE) 
    # mutate(SUBSAMPLE = rep(i, n()))
      # filter (MISSING_DATA == 'no.imputation')
    data[[i]] <- subsample.data
  }
  data <- as_data_frame(bind_rows(data))
  return(data)
}
