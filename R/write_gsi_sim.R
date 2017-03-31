# write a gsi_sim file

#' @name write_gsi_sim

#' @title Write a gsi_sim file from a data frame (wide or long/tidy).

#' @description Write a gsi_sim file from a data frame (wide or long/tidy). 
#' Used internally in \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy genomic data set in the working directory tidy formats.
#' \emph{How to get a tidy data frame ?}
#' Look for \pkg{stackr} \code{\link{tidy_genomic_data}}.

#' @param pop.levels (option, string) This refers to the levels in a factor. In this 
#' case, the id of the pop.
#' Use this argument to have the pop ordered your way instead of the default 
#' alphabetical or numerical order. e.g. \code{pop.levels = c("QUE", "ONT", "ALB")} 
#' instead of the default \code{pop.levels = c("ALB", "ONT", "QUE")}. 
#' Default: \code{pop.levels = NULL}. If you find this too complicated, there is also the
#' \code{strata} argument that can do the same thing, see below.

#' @param pop.labels (optional, string) Use this argument to rename/relabel
#' your pop or combine your pop. e.g. To combine \code{"QUE"} and \code{"ONT"} 
#' into a new pop called \code{"NEW"}:
#' (1) First, define the levels for your pop with \code{pop.levels} argument: 
#' \code{pop.levels = c("QUE", "ONT", "ALB")}. 
#' (2) then, use \code{pop.labels} argument: 
#' \code{pop.levels = c("NEW", "NEW", "ALB")}.#' 
#' To rename \code{"QUE"} to \code{"TAS"}:
#' \code{pop.labels = c("TAS", "ONT", "ALB")}.
#' Default: \code{pop.labels = NULL}. If you find this too complicated, there is also the
#' \code{strata} argument that can do the same thing, see below.
#' 
#' @param strata (optional) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. 
#' Default: \code{strata = NULL}. Use this argument to rename or change 
#' the populations id with the new \code{STRATA} column.
#' The \code{STRATA} column can be any hierarchical grouping.

#' @param filename The name of the file written to the working directory.

# @param ... other parameters passed to the function.

#' @return A gsi_sim input file is saved to the working directory. 
#' @export
#' @rdname write_gsi_sim

#' @importFrom stackr change_pop_names

#' @importFrom data.table fread dcast.data.table as.data.table
#' @importFrom tibble as_data_frame
#' @importFrom tidyr separate gather unite 
#' @importFrom dplyr n_distinct rename mutate select left_join arrange
#' @importFrom stringi stri_replace_all_regex stri_paste stri_replace_all_fixed
#' @importFrom utils write.table count.fields


#' @references Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008)
#' An improved method for predicting the accuracy of genetic stock identification.
#' Canadian Journal of Fisheries and Aquatic Sciences 65, 7:1475-1486.
#' @references Anderson, E. C. (2010) Assessing the power of informative subsets of
#' loci for population assignment: standard methods are upwardly biased.
#' Molecular ecology resources 10, 4:701-710.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_gsi_sim <- function(
  data, 
  pop.levels = NULL, 
  pop.labels = NULL, 
  strata = NULL,
  filename = "gsi_sim.unname.txt"
  ) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file necessary to write the gsi_sim file is missing")
  
  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...
  # pop.levels, pop.labels, pop.select, strata, etc
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  if (!is.null(pop.labels)) {
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # Import data
  if (is.vector(data)) {
    input <- stackr::tidy_wide(data = data, import.metadata = FALSE)
  } else {
    input <- data
  }  
  colnames(input) <- stringi::stri_replace_all_fixed(
    str = colnames(input), 
    pattern = "GENOTYPE", 
    replacement = "GT", 
    vectorize_all = FALSE)
  
  # remove space in POP_ID
  input$POP_ID <- stringi::stri_replace_all_fixed(input$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)
  
  # Info for gsi_sim input -----------------------------------------------------
  n.individuals <- dplyr::n_distinct(input$INDIVIDUALS)  # number of individuals
  
  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  # n.pop <- dplyr::n_distinct(input$POP_ID)
  n.markers <- dplyr::n_distinct(input$MARKERS)          # number of markers
  list.markers <- unique(input$MARKERS)                  # list of markers
  
  
  # Spread/dcast in wide format ------------------------------------------------------
  input <- dplyr::select(input, MARKERS, POP_ID, INDIVIDUALS, GT) %>%
    tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
    tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% 
    dplyr::arrange(MARKERS) %>%
    tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELES, sep = "_") %>%
    dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>%
    data.table::as.data.table(x = .) %>% 
    data.table::dcast.data.table(data = .,
                                 formula = POP_ID + INDIVIDUALS ~ MARKERS_ALLELES, 
                                 value.var = "GT") %>% 
    tibble::as_data_frame(.)
  
  # change sep in individual name
  input$INDIVIDUALS <- stringi::stri_replace_all_fixed(
    str = input$INDIVIDUALS, 
    pattern = c("_", ":"), 
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )
  
  # population levels and strata------------------------------------------------
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      # message("strata file: yes")
      number.columns.strata <- max(utils::count.fields(strata, sep = "\t"))
      col.types <- stringi::stri_join(rep("c", number.columns.strata), collapse = "")
      suppressMessages(strata.df <- readr::read_tsv(
        file = strata,
        col_names = TRUE,
        col_types = col.types
      ) %>% 
        dplyr::rename(POP_ID = STRATA))
    } else {
      # message("strata object: yes")
      colnames(strata) <- stringi::stri_replace_all_fixed(
        str = colnames(strata), 
        pattern = "STRATA", 
        replacement = "POP_ID", 
        vectorize_all = FALSE
      )
      strata.df <- strata
    }
    
    # Remove potential whitespace in pop_id
    strata.df$POP_ID <- stringi::stri_replace_all_fixed(
      strata.df$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE
    )
    
    strata.df$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = strata.df$INDIVIDUALS, 
      pattern = c("_", ":"), 
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )
    
    input <- dplyr::select(.data = input, -POP_ID) %>% 
      dplyr::left_join(strata.df, by = "INDIVIDUALS")
  }
  
  # using pop.levels and pop.labels info if present
  input <- stackr::change_pop_names(data = input,
                                    pop.levels = pop.levels,
                                    pop.labels = pop.labels)
  
  input <- dplyr::arrange(.data = input, POP_ID, INDIVIDUALS)
  
  # write gsi_sim file ---------------------------------------------------------
  
  # open the connection to the file
  filename.connection <- file(filename, "w") 
  
  # Line 1: number of individuals and the number of markers
  writeLines(text = stringi::stri_join(n.individuals, n.markers, sep = " "), con = filename.connection, sep = "\n")
  
  # Line 2 and + : List of markers
  writeLines(text = stringi::stri_join(list.markers, sep = "\n"), con = filename.connection, sep = "\n")
  
  # close the connection to the file
  close(filename.connection) # close the connection
  
  # remaining lines, individuals and genotypes
  pop <- input$POP_ID # Create a vector with the population
  input <- suppressWarnings(dplyr::select(.data = input, -POP_ID))  # remove pop id
  gsi_sim.split <- split(input, pop)  # split gsi_sim by populations
  pop.string <- as.character(unique(pop))
  for (k in pop.string) {
    utils::write.table(x = as.data.frame(stringi::stri_join("pop", k, sep = " ")), file = filename, append = TRUE, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
    # readr::write_delim(x = as.data.frame(stringi::stri_join("pop", k, sep = " ")), path = filename, delim = "\n", append = TRUE, col_names = FALSE)
    # readr::write_delim(x = gsi_sim.split[[k]], path = filename, delim = " ", append = TRUE, col_names = FALSE)
    utils::write.table(x = gsi_sim.split[[k]], file = filename, append = TRUE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  }
  
  gsi_sim.split <- input <-pop <- pop.string <- NULL
  return(filename)
} # End write_gsi function

