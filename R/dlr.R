#' @name dlr
#' @title Genotype likelihood ratio distance (Dlr)
#' @description The function computes Paetkau's et al. (1997) genotype likelihood
#' ratio distance (Dlr).
#' @param data The output assignment file (home likelihood or
#' likelihood ratio statistics) from GENODIVE.
#' @param strata A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} column is used here as the populations id of your sample. 

#' @param filename (optional) Name of the file prefix for
#' the matrix and the table written in the working directory. 
#' @return A list with 3 objects of class: table ($dlr.table), dist (a lower
#' diagonal matrix, $dlr.dist), data.frame (a mirrored matrix, $dlr.matrix).

#' @importFrom dplyr mutate select filter group_by ungroup filter_ mutate_ summarise ungroup left_join rename
#' @importFrom readr read_delim write_tsv read_table read_tsv
#' @importFrom lazyeval interp
#' @importFrom stringi stri_dup stri_join stri_replace_all_fixed stri_sub
#' @importFrom stats as.dist dist
#' @importFrom utils combn
#' @importFrom tibble rownames_to_column data_frame
#' @importFrom purrr discard flatten_dbl

#' @export 
#' @rdname dlr

#' @references Paetkau D, Slade R, Burden M, Estoup A (2004) 
#' Genetic assignment methods for the direct, real-time estimation of 
#' migration rate: a simulation-based exploration of accuracy and power. 
#' Molecular Ecology, 13, 55-65.
#' @references Paetkau D, Waits LP, Clarkson PL, Craighead L, Strobeck C (1997)
#'  An empirical evaluation of genetic distance statistics using microsatellite
#'   data from bear (Ursidae) populations. Genetics, 147, 1943-1957.
#' @references Meirmans PG, Van Tienderen PH (2004) genotype and genodive: 
#' two programs for the analysis of genetic diversity of asexual organisms. 
#' Molecular Ecology Notes, 4, 792-794.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

dlr <- function(data, strata, filename = NULL) {
  cat("#######################################################################\n")
  cat("########################### assigner::Dlr #############################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  
  if (missing(data)) stop("GenoDive file missing")
  if (missing(strata)) stop("Strata file missing")
  
  # import and modify the assignment file form GenoDive-------------------------
  message("Importing GenoDive assignment results")
  temp.file <- suppressWarnings(suppressMessages(readr::read_table(file = data)))
  max.lines = nrow(temp.file) - 1
  skip.number <- which(stringi::stri_detect_fixed(str = temp.file$X1,
                                                  pattern = "membership")) + 2
  assignment <- suppressMessages(
    suppressWarnings(
      readr::read_delim(
        data,
        delim = "\t",
        skip = skip.number,
        n_max = max.lines,
        col_names = TRUE,
        progress = interactive())) %>% 
      dplyr::filter(!is.na(Current)) %>%
      dplyr::rename(
        INDIVIDUALS = Individual, POP_ID = Current, INFERRED = Inferred,
        LIK_MAX = Lik_max, LIK_HOME = Lik_home, LIK_RATIO = Lik_ratio))
  
  temp.file <- skip.number <- NULL
  
  message("Importing strata file")
  strata.df <- readr::read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
    dplyr::rename(POP_ID = STRATA)
  
  # check that same number of individuals and pop...
  if (!identical(sort(assignment$INDIVIDUALS), sort(strata.df$INDIVIDUALS))) {
    stop("Assignment file and strata don't have the same individuals")
  }
  
  assignment.pop <- ncol(assignment) - 6
  strata.pop <- dplyr::n_distinct(strata.df$POP_ID)
  if (assignment.pop != strata.pop) stop("Assignment file and strata don't have the same number of populations")
  
  fixed.header <- c("INDIVIDUALS", "POP_ID", "INFERRED", "LIK_MAX", "LIK_HOME", "LIK_RATIO")
  header.pop <- purrr::discard(.x = colnames(assignment), .p = colnames(assignment) %in% fixed.header)
  
  get.pop <- strata.df %>% 
    dplyr::filter(INDIVIDUALS %in% header.pop)
  
  strata.df <- NULL
  
  # Change header for the real pop names
  colnames(assignment) <- stringi::stri_replace_all_fixed(
    str = colnames(assignment), pattern = get.pop$INDIVIDUALS, replacement = get.pop$POP_ID, vectorize_all = FALSE)
  
  # Change POP_ID and inferred for the real pop name
  assignment <- assignment %>% 
    dplyr::mutate(
      POP_ID = stringi::stri_replace_all_fixed(
        str = POP_ID, pattern = get.pop$INDIVIDUALS, replacement = get.pop$POP_ID, vectorize_all = FALSE),
      INFERRED = stringi::stri_replace_all_fixed(
        str = INFERRED, pattern = get.pop$INDIVIDUALS, replacement = get.pop$POP_ID, vectorize_all = FALSE)
    )
  
  
  # Dlr relative for one combination of pop-------------------------------------
  dlr_relative <- function(pop1, pop2){
    
    dlr <- suppressWarnings(
      assignment %>%
        dplyr::filter_(lazyeval::interp(~ POP_ID == as.name(pop1) | POP_ID == as.name(pop2))) %>%
        dplyr::group_by(INDIVIDUALS) %>%
        dplyr::mutate_(
          RATIO1 = lazyeval::interp(~pop1 - pop2, pop1 = as.name(pop1), pop2 = as.name(pop2)),
          RATIO2 = lazyeval::interp(~pop2 - pop1, pop1 = as.name(pop1), pop2 = as.name(pop2)),
          RATIO = lazyeval::interp(~ifelse(POP_ID == pop1, c.RATIO1, c.RATIO2),
                                   POP_ID = quote(POP_ID), pop1 = as.name("pop1"),
                                   c.RATIO1 = quote(RATIO1), c.RATIO2 = quote(RATIO2))) %>% 
        dplyr::group_by(POP_ID) %>%
        dplyr::summarise(DLR_RELATIVE = (sum(RATIO)/length(RATIO)^2)) %>%
        dplyr::ungroup(.) %>%
        dplyr::summarise(DLR_RELATIVE = sum(DLR_RELATIVE)/2)
    )
    return(dlr)
  }#End dlr_relative
  
  # All combination of populations----------------------------------------------
  message("Calculating Dlr...")
  pop.pairwise <- utils::combn(unique(get.pop$POP_ID), 2)
  pop.pairwise <- matrix(pop.pairwise, nrow = 2)
  
  # Dlr for all pairwise populations--------------------------------------------
  dlr.all.pop <- as.numeric()
  for (i in 1:ncol(pop.pairwise)) {
    dlr.all.pop[i] <- dlr_relative(pop1 = pop.pairwise[1,i], 
                                   pop2 = pop.pairwise[2,i])
  }
  dlr.all.pop <- as.numeric(dlr.all.pop)
  pop.pairwise <- NULL
  
  # Table with Dlr--------------------------------------------------------------
  names.pairwise <- utils::combn(unique(get.pop$POP_ID), 2, paste, collapse = '-')
  
  dlr.table <- tibble::data_frame(PAIRWISE_POP = names.pairwise, DLR = dlr.all.pop) %>%
    dplyr::mutate(DLR = round(as.numeric(DLR), 2))
  
  
  # Dist and Matrix-------------------------------------------------------------
  dlr.dist <- stats::dist(1:length(unique(get.pop$POP_ID)))
  dlr.dist.matrix <- dlr.all.pop
  attributes(dlr.dist.matrix) <- attributes(dlr.dist)
  dlr.dist.matrix <- as.matrix(dlr.dist.matrix)
  colnames(dlr.dist.matrix) <- rownames(dlr.dist.matrix) <- unique(get.pop$POP_ID)
  dlr.dist.matrix <- stats::as.dist(dlr.dist.matrix)
  
  dlr.matrix <- as.data.frame(as.matrix(dlr.dist.matrix)) %>%
    tibble::rownames_to_column(df = ., var = "POP")
  
  get.pop <- NULL
  cat("############################### RESULTS ###############################\n")
  # Results---------------------------------------------------------------------
  res <- list(
    assignment = assignment,
    dlr.table = dlr.table,
    dlr.dist = dlr.dist.matrix,
    dlr.matrix = dlr.matrix
    )
  
  # Write file to working directory --------------------------------------------
  if (is.null(filename)) {
    message("Writing files to directory: no")
  } else {
    # saving table
    filename.table <- stringi::stri_join(filename, "table.tsv", sep = ".") 
    readr::write_tsv(dlr.table, filename.table)
    
    # saving matrix
    filename.matrix <- stringi::stri_join(filename, "matrix.tsv", sep = ".") 
    readr::write_tsv(dlr.matrix, filename.matrix)
    message("Writing files to directory: yes")
    message("Filenames : ", "\n", filename.table, "\n", filename.matrix)
  }
  timing <- proc.time() - timing
  message("Computation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res)
}

# Dlr absolute
# dlr.absolute <- assignment %>%
#   group_by(Individuals) %>%
#   mutate_(
#     RATIO1 = lazyeval::interp(~pop1 - pop2, pop1 = as.name(pop1), pop2 = as.name(pop2)),
#     RATIO2 = lazyeval::interp(~pop2 - pop1, pop1 = as.name(pop1), pop2 = as.name(pop2)),
#     RATIO = lazyeval::interp(~ifelse(Populations == pop1, c.RATIO1, c.RATIO2), Populations = quote(Populations), pop1 = as.name("pop1"), c.RATIO1 = quote(RATIO1), c.RATIO2 = quote(RATIO2))) %>% 
#     group_by(Populations) %>%
#     summarise(DLR_ABSOLUTE = sum(RATIO)/length(RATIO)) %>%
#     ungroup %>%
#     summarise(DLR_ABSOLUTE = sum(DLR_ABSOLUTE)/2)
