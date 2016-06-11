# write a gsi_sim file

#' @name write_gsi_sim

#' @title Write a gsi_sim file from a data frame (wide or long/tidy).

#' @description Write a gsi_sim file from a data frame (wide or long/tidy). 
#' Used internally in \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. To import, the function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}}. See details for more info.
#' 
#' \emph{How to get a tidy data frame ?}
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{tidy_genomic_data}} can transform 6 genomic data formats 
#' in a tidy data frame.

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

#' @param ... other parameters passed to the function.

#' @return A gsi_sim input file is saved to the working directory. 
#' @export
#' @rdname write_gsi_sim
#' @import stackr
#' @import dplyr
#' @import stringi
#' @importFrom data.table fread
#' @importFrom data.table dcast.data.table
#' @importFrom data.table as.data.table

#' @details \strong{Input data:}
#'  
#' To discriminate the long from the wide format, 
#' the function \pkg{stackr} \code{\link[stackr]{read_long_tidy_wide}} searches 
#' for \code{MARKERS or LOCUS} in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns: 
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals), 
#' the remaining columns are the markers in separate columns storing genotypes.
#' 
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info. 
#' (e.g. from a VCF see \pkg{stackr} \code{\link{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID}, 
#' \code{MARKERS or LOCUS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#' 
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 of 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 of 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#' 
#' \emph{How to get a tidy data frame ?}
#' \pkg{stackr} \code{\link{tidy_genomic_data}} can transform 6 genomic data formats 
#' in a tidy data frame.


#' @references Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008)
#' An improved method for predicting the accuracy of genetic stock identification.
#' Canadian Journal of Fisheries and Aquatic Sciences 65, 7:1475-1486.
#' @references Anderson, E. C. (2010) Assessing the power of informative subsets of
#' loci for population assignment: standard methods are upwardly biased.
#' Molecular ecology resources 10, 4:701-710.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_gsi_sim <- function (
  data, 
  pop.levels = NULL, 
  pop.labels = NULL, 
  strata = NULL,
  filename = "gsi_sim.unname.txt", 
  ...) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file necessary to write the gsi_sim file is missing")
  if (missing(pop.levels)) pop.levels <- NULL
  if (missing(pop.labels)) pop.labels <- NULL
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (missing(strata)) strata <- NULL
  if (missing(filename)) filename <- "gsi_sim.unname.txt"
  
  # data <- "skipjack.filtered_tidy.tsv" #long
  # data <- "skipjack.wide.test.tsv" #wide
  # data <- data.select
  
  # Import data
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = FALSE)
  
  # Info for gsi_sim input -----------------------------------------------------
  n.individuals <- n_distinct(input$INDIVIDUALS)  # number of individuals
  if ("LOCUS" %in% colnames(input)) input <- rename(.data = input, MARKERS = LOCUS)
  n.markers <- n_distinct(input$MARKERS)          # number of markers
  list.markers <- unique(input$MARKERS)           # list of markers
  
  # Spread/dcast in wide format ------------------------------------------------------
  input <- input %>%
    tidyr::separate(data = ., col = GT, into = .(A1, A2), sep = 3, remove = TRUE) %>% 
    tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% 
    arrange(MARKERS) %>%
    tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELES, sep = "_") %>%
    arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES)
  
  input <- data.table::dcast.data.table(
    as.data.table(input), 
    formula = POP_ID + INDIVIDUALS ~ MARKERS_ALLELES, 
    value.var = "GT"
  ) %>% 
    as_data_frame()
  
  # change sep in individual name
  input$INDIVIDUALS <- stri_replace_all_fixed(
    str = input$INDIVIDUALS, 
    pattern = c("_", ":"), 
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )
  
  # population levels ----------------------------------------------------------
  if (is.null(strata)){ # no strata
    if(is.null(pop.levels)) { # no pop.levels
      if (is.factor(input$POP_ID)) {
        input$POP_ID <- droplevels(x = input$POP_ID)
      } else {
        input$POP_ID <- factor(input$POP_ID)
      }
    } else { # with pop.levels
      input <- input %>%
        mutate( # Make population ready
          POP_ID = factor(
            stri_replace_all_regex(
              POP_ID, 
              stri_paste("^", pop.levels, "$", sep = ""), 
              pop.labels,
              vectorize_all = FALSE), 
            levels = unique(pop.labels), 
            ordered = TRUE
          )
        )
    }
  } else { # Make population ready with the strata provided
    if (is.vector(strata)) {
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA)
    } else {
      strata.df <- strata
    }
    # change sep in individual name
    strata.df$INDIVIDUALS <- stri_replace_all_fixed(
      str = strata.df$INDIVIDUALS, 
      pattern = c("_", ":"), 
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )
    
    if(is.null(pop.levels)) { # no pop.levels
      input <- input %>%
        select(-POP_ID) %>% 
        mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        left_join(strata.df, by = "INDIVIDUALS") %>% 
        mutate(POP_ID = factor(POP_ID))
    } else { # with pop.levels
      input <- input %>%
        select(-POP_ID) %>% 
        mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        left_join(strata.df, by = "INDIVIDUALS") %>%
        mutate(
          POP_ID = factor(
            stri_replace_all_regex(
              POP_ID, 
              stri_paste("^", pop.levels, "$", sep = ""), 
              pop.labels, 
              vectorize_all = FALSE
            ),
            levels = unique(pop.labels), ordered = TRUE
          )
        )
    }
  }
  
  # write gsi_sim file ---------------------------------------------------------
  
  # open the connection to the file
  filename.connection <- file(filename, "w") 
  
  # Line 1: number of individuals and the number of markers
  writeLines(text = stri_join(n.individuals, n.markers, sep = " "), con = filename.connection, sep = "\n")
  
  # Line 2 and + : List of markers
  writeLines(text = stri_paste(list.markers, sep = "\n"), con = filename.connection, sep = "\n")
  
  # close the connection to the file
  close(filename.connection) # close the connection
  
  # remaining lines, individuals and genotypes
  pop <- input$POP_ID  # Create a vector with the population
  input <- suppressWarnings(input %>% select(-POP_ID))  # remove pop id
  gsi_sim.split <- split(input, pop)  # split gsi_sim by populations
  
  for (k in levels(pop)) {
    write_delim(x = as.data.frame(stri_join("pop", k, sep = " ")), path = filename, delim = "\n", append = TRUE, col_names = FALSE)
    write_delim(x = gsi_sim.split[[k]], path = filename, delim = " ", append = TRUE, col_names = FALSE)
  }
  
  gsi_sim.split<- NULL
  input <- NULL
  return(filename)
} # End write_gsi function

