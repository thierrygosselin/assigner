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
#' @import dplyr
#' @import stringi


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
    fst.files.list <- list.files(path = stri_paste(dir.path, "/", i), pattern = "fst.ranked", full.names = FALSE)
    data.fst <- list()
    for (j in fst.files.list) {
      fst.file <- read_tsv(file = stri_paste(dir.path, "/", i, "/", j), col_names = TRUE) %>% 
        mutate(
          SUBSAMPLE = rep(i, n()),
          ITERATIONS = rep(j, n())
        )
      data.fst[[j]] <- fst.file
    }
    data.fst <- as_data_frame(bind_rows(data.fst))
    data.subsample[[i]] <- data.fst
  }
  data <- as_data_frame(bind_rows(data.subsample))
  return(data)
}
