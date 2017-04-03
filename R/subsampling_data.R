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
    if (subsample < 1) { # proportion
      subsample.select <- ind.pop.df %>%
        dplyr::group_by(POP_ID) %>%
        dplyr::sample_frac(tbl = ., size = subsample, replace = FALSE)# sampling individuals for each pop
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
