# write a gsi_sim file

#' @name write_gsi_sim
#' @title Used internally in assigner to write a gsi_sim file
#' @description Write a gsi_sim file
#' @param data A file or object in the global environment containing at least 
#' these 2 id columns: 
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals.), 
#' the remaining columns are the markers.
#' @param markers.names (character) Name of the markers.
#' @param i Iterations.
#' @param m Markers number used.
#' @param directory.subsample The name of the subsample directory.
#' @param filename The name of the file written to the working directory.
#' @param ... other parameters passed to the function.
#' @return A gsi_sim file is saved to the working directory. 
#' @export
#' @rdname write_gsi_sim
#' @import dplyr
#' @import stringi

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


write_gsi_sim <- function (data, 
                           markers.names = NULL, 
                           directory.subsample = NULL, 
                           filename = "gsi_sim.unname.txt", 
                           i = NULL, 
                           m = NULL, 
                           ...) {
  # data <- data.select
  # markers.names = markers.names
  # filename = filename
  # i = i
  # m = m
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file necessary to write the gsi_sim file is missing")
  if (missing(markers.names)) markers.names <- NULL
  if (missing(i)) i <- NULL
  if (missing(m)) m <- NULL
  if (missing(directory.subsample)) directory.subsample <- NULL
  if (missing(filename)) filename <- "gsi_sim.unname.txt"
  
   
  data$POP_ID <- droplevels(x = data$POP_ID)
  n.individuals <- n_distinct(data$INDIVIDUALS)  # number of individuals
  pop <- data$POP_ID  # Create a vector with the population ordered by levels
  data <- suppressWarnings(data %>% select(-POP_ID))  # remove pop id
  gsi_sim.split <- split(data, pop)  # split gsi_sim by populations
  filename <- filename  # gsi_sim filename
  
  # Line 1: number of individuals and the number of markers
  # line1_gsi_sim <- as.data.frame(stri_join(n.individuals, m, sep = " "))
  # write.table(line1_gsi_sim, file = paste0(directory.subsample, filename), col.names = FALSE, row.names = FALSE, quote = FALSE)
  filename.connection <- file(paste0(directory.subsample, filename), "w") # open the connection to the file
  writeLines(text = stri_join(n.individuals, m, sep = " "), con = filename.connection, sep = "\n")

  # Markers names
  # loci.table <- as.data.frame(markers.names)
  # write_delim(x = loci.table, path = paste0(directory.subsample, filename), delim = "\n", append = TRUE, col_names = FALSE)
  writeLines(text = stri_paste(markers.names, sep = "\n"), con = filename.connection, sep = "\n")
  close(filename.connection) # close the connection
  
  # remaining lines, individuals and genotypes
  for (k in levels(pop)) {
    # pop.line <- as.data.frame(stri_join("pop", k, sep = " "))
    write_delim(x = as.data.frame(stri_join("pop", k, sep = " ")), path = paste0(directory.subsample, filename), delim = "\n", append = TRUE, col_names = FALSE)
    write_delim(x = gsi_sim.split[[k]], path = paste0(directory.subsample, filename), delim = " ", append = TRUE, col_names = FALSE)
  }
  return(filename)
} # End write_gsi function
