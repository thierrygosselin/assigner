#' @name dlr
#' @title Genotype likelihood ratio distance (Dlr)
#' @description The function computes Paetkau's et al. (1997) genotype likelihood
#' ratio distance (Dlr).
#' @param data The output assignment file (home likelihood or
#' likelihood ratio statistics) from GENODIVE.
#' @param strata A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} column is used here as the populations id of your sample.

#' @param plots (optional) Generate Dlr plots for all the pairwise populations
#' in the dataset. The plots are \code{ggplot2} objects that can be modified with
#' the proper \code{ggplot2} syntax. Default: \code{plots = FALSE}.
#' @param filename (optional) Name of the file prefix for
#' the matrix and the table written in the working directory.
#' @param parallel.core (optional) The number of core for parallel computation.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @return A list with 5 objects:
#' \enumerate{
#' \item the assignment results ($assignment),
#' \item the dlr pairwise table ($dlr.table),
#' \item the lower diagonal dlr distance matrix ($dlr.dist),
#' \item a data.frame with the dlr distance mirrored ($dlr.matrix),
#' \item the list of dlr plots ($dlr.plots)
#' }

#' @importFrom dplyr mutate select filter group_by ungroup filter_ mutate_ summarise ungroup left_join rename
#' @importFrom readr read_delim write_tsv read_table read_tsv
#' @importFrom stringi stri_dup stri_join stri_replace_all_fixed stri_sub
#' @importFrom stats as.dist dist
#' @importFrom utils combn
#' @importFrom tibble rownames_to_column tibble
#' @importFrom purrr discard flatten_dbl
#' @importFrom parallel detectCores
#' @importFrom rlang .data UQ

#' @examples
#' \dontrun{
#' dlr <- assigner::dlr(
#' data = "assignment.gdv", strata = "my.strata.tsv", plots = TRUE)
#'
#' # to get the plots list:
#' plot.list <- dlr$dlr.plots
#' # access and isolate in different object a plot with $
#' }


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

dlr <- function(
  data, strata, plots = FALSE, filename = NULL,
  parallel.core = parallel::detectCores() - 1) {
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

  # All combination of populations----------------------------------------------
  message("Calculating Dlr...")
  pop.pairwise <- utils::combn(unique(get.pop$POP_ID), 2)
  # pop.pairwise <- matrix(pop.pairwise, nrow = 2)

  # Dlr for all pairwise populations--------------------------------------------
  dlr.all.pop <- as.numeric()
  for (i in 1:ncol(pop.pairwise)) {
    dlr.all.pop[i] <- dlr_relative(data = assignment,
                                   pop1 = pop.pairwise[1,i],
                                   pop2 = pop.pairwise[2,i])
  }
  dlr.all.pop <- as.numeric(dlr.all.pop)
  # pop.pairwise <- NULL

  # Table with Dlr--------------------------------------------------------------
  names.pairwise <- utils::combn(unique(get.pop$POP_ID), 2, paste, collapse = '-')

  dlr.table <- tibble::tibble(PAIRWISE_POP = names.pairwise, DLR = dlr.all.pop) %>%
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

  # get.pop <- NULL

  # create list of results -----------------------------------------------------
  res <- list(
    assignment = assignment,
    dlr.table = dlr.table,
    dlr.dist = dlr.dist.matrix,
    dlr.matrix = dlr.matrix
  )
  # Dlr plots ------------------------------------------------------------------
  if (plots) {
    # all combination of individual pair
    pop.pairwise <- utils::combn(unique(get.pop$POP_ID), 2, simplify = FALSE)

    # get the number of pairwise comp.
    number.pairwise <- length(pop.pairwise)

    # Optimizing cpu usage
    # if (number.pairwise <= 50) {
    #   round.cpu <- floor(number.pairwise / parallel.core)
    # } else {
    #   round.cpu <- floor(number.pairwise / (50 * parallel.core))
    # }
    # # as.integer is usually twice as light as numeric vector...
    # split.vec <- as.integer(floor((parallel.core * round.cpu * (1:number.pairwise - 1) / number.pairwise) + 1))
    # pop.pairwise <- split(x = pop.pairwise, f = split.vec)

    message("Generating ", number.pairwise, " Dlr plots...")
    # round.cpu <- split.vec <- number.pairwise <- NULL

    dlr.plots <- list()
    dlr.plots <- .assigner_parallel(
      X = pop.pairwise,
      FUN = plot_dlr,
      mc.preschedule = FALSE,
      mc.silent = FALSE,
      mc.cleanup = TRUE,
      mc.cores = parallel.core,
      data = assignment
    ) %>% purrr::flatten(.)
    res$dlr.plots <- dlr.plots
  }#End Dlr plots

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
#   dplyr::group_by(INDIVIDUALS) %>%
#   dplyr::mutate(RATIO = dplyr::if_else(
#     POP_ID == rlang::UQ(pop1),
#     rlang::.data[[pop1]] - rlang::.data[[pop2]],
#     rlang::.data[[pop2]] - rlang::.data[[pop1]])
#   ) %>%
#   dplyr::group_by(POP_ID) %>%
#   dplyr::summarise(DLR_ABSOLUTE = sum(RATIO)/length(RATIO)) %>%
#   dplyr::ungroup(.) %>%
#   dplyr::summarise(DLR_ABSOLUTE = sum(DLR_ABSOLUTE)/2)

# Internal nested functions: ---------------------------------------------------

# dlr relative -----------------------------------------------------------------
#' @title dlr_relative
#' @description Calculate relative Dlr between 2 pops
#' @rdname dlr_relative
#' @export
#' @keywords internal
dlr_relative <- function(data, pop1, pop2){
  dlr <- data %>%
    dplyr::filter(POP_ID == rlang::UQ(pop1) | POP_ID == rlang::UQ(pop2)) %>%
    dplyr::group_by(INDIVIDUALS) %>%
    dplyr::mutate(RATIO = dplyr::if_else(
      POP_ID == rlang::UQ(pop1),
      rlang::.data[[pop1]] - rlang::.data[[pop2]],
      rlang::.data[[pop2]] - rlang::.data[[pop1]])
    ) %>%
    dplyr::group_by(POP_ID) %>%
    dplyr::summarise(DLR_RELATIVE = (sum(RATIO)/length(RATIO)^2)) %>%
    dplyr::ungroup(.) %>%
    dplyr::summarise(DLR_RELATIVE = sum(DLR_RELATIVE)/2)
  return(dlr)
}#End dlr_relative

# plot_dlr ---------------------------------------------------------------------
#' @name plot_dlr
#' @title Assignment plot of genotype likelihood (Dlr).
#' @description The function generate a figure similar to Paetkau's et al. (2004)
#' Fig 6.
#' @param data The assigment object from
#' \code{\link[assigner]{dlr}} output.
#' @param pop.pairwise List of pairwise populations to generate the Dlr plot.
#' @param dlr (optional) Character string with Dlr value.
#' @param x.dlr (optional) Position to the x-axis of the Dlr value.
#' @param y.dlr (optional) Position to the y-axis of the Dlr value.
#' @param fst (optional) Character string with Fst value.
#' @param x.fst (optional) Position to the x-axis of the Fst value.
#' @param y.fst (optional) Position to the y-axis of the Fst value.
#' @param filename (optional) Name of the figure written
#' in the working directory.
#' @param plot.width (optional) Width in cm of the figure.
#' @param plot.height (optional) height in cm of the figure.
#' @param plot.dpi (optional) Number of dpi for the figure (e.g 600).
#' @return A list with the assignment table and the assignment plot.
#' @export
#' @rdname plot_dlr
#' @keywords internal

#' @importFrom stringi stri_join stri_sub stri_replace_all_fixed stri_detect_fixed stri_replace_na
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs sample_n sample_frac one_of
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid annotate geom_abline scale_shape_manual scale_fill_manual geom_abline aes_string geom_point scale_x_continuous scale_y_continuous
#' @importFrom readr read_table read_delim read_tsv
#' @importFrom purrr discard flatten_dbl

#' @references Paetkau D, Slade R, Burden M, Estoup A (2004)
#' Genetic assignment methods for the direct, real-time estimation of migration
#' rate: a simulation-based exploration of accuracy and power.
#' Molecular Ecology, 13, 55-65.
#' @references Meirmans PG, Van Tienderen PH (2004) genotype and genodive:
#' two programs for the analysis of genetic diversity of asexual organisms.
#' Molecular Ecology Notes, 4, 792-794.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

plot_dlr <- function(
  data,
  pop.pairwise,
  dlr = NULL,
  x.dlr = NULL,
  y.dlr = NULL,
  fst = NULL,
  x.fst = NULL,
  y.fst = NULL,
  filename = NULL,
  plot.width = NULL,
  plot.height = NULL,
  plot.dpi = NULL) {

  if (missing(data)) stop("GenoDive file missing")
  if (missing(pop.pairwise)) stop("Missing pop.pairwise argument")

  # import data ----------------------------------------------------------------
  assignment.select <- dplyr::rename(data, Populations = POP_ID)

  POPA <- pop.pairwise[1]
  POPB <- pop.pairwise[2]

  # Check that POPA and POPB are found in the data
  if (!POPA %in% colnames(assignment.select)) stop("POPA value is not in the assignment file")
  if (!POPB %in% colnames(assignment.select)) stop("POPB value is not in the assignment file")


  pop.select <- c(POPA, POPB)
  fixed.header <- c("INDIVIDUALS", "Populations", "INFERRED", "LIK_MAX", "LIK_HOME", "LIK_RATIO")

  dlr.plot.name <- stringi::stri_join("dlr_plot_pop_", stringi::stri_join(pop.select, collapse = "_"))

  assignment.select <- suppressWarnings(
    assignment.select %>%
      dplyr::filter(Populations %in% pop.select) %>%
      dplyr::select(dplyr::one_of(c(fixed.header, pop.select))))
  
  scale.temp <- suppressWarnings(
    unique(
      dplyr::select(assignment.select,
                    dplyr::one_of(pop.select)) %>%
        purrr::flatten_dbl(.)
    ))
  min.scale <- min(scale.temp)
  max.scale <- max(scale.temp)

  assignment.plot  <- ggplot2::ggplot(assignment.select, ggplot2::aes_string(x = POPA, y = POPB)) +
    ggplot2::geom_point(ggplot2::aes(fill = Populations, shape = Populations), na.rm = TRUE, alpha = 0.8, size = 4) +
    ggplot2::geom_abline(slope = 1) +
    ggplot2::scale_x_continuous(name = stringi::stri_join("Log (genotype likelihood) population: ", POPA), limits = c(min.scale, max.scale)) +
    ggplot2::scale_y_continuous(name = stringi::stri_join("Log (genotype likelihood) population: ", POPB), limits = c(min.scale, max.scale)) +
    ggplot2::scale_shape_manual(values = c(21, 24)) +
    ggplot2::scale_fill_manual(values = c("black", NA)) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.y = ggplot2::element_text(angle = 0, size = 10, family = "Helvetica", face = "bold")
    )

  if (is.null(dlr)) {
    assignment.plot <- assignment.plot
  } else {
    assignment.plot <- assignment.plot + ggplot2::annotate("text", x = x.dlr, y = y.dlr,
                                                           label = dlr, colour = "black")
  }

  if (is.null(fst)) {
    assignment.plot <- assignment.plot
  } else {
    assignment.plot <- assignment.plot + ggplot2::annotate("text", x = x.fst, y = y.fst,
                                                           label = fst, colour = "black")
  }


  if (is.null(filename)) {
    # message("Saving was not selected...")
  } else {
    message("Saving the figure was selected, the filename: ", filename)

    if (is.null(plot.width)) stop("plot.width argument is required")
    if (is.null(plot.height)) stop("plot.height argument is required")
    if (is.null(plot.dpi)) stop("plot.dpi argument is required")
    ggplot2::ggsave(filename, width = plot.width, height = plot.height,
                    dpi = plot.dpi, units = "cm", useDingbats = FALSE)
  }
  res <- list(assignment.plot)
  names(res) <- dlr.plot.name
  return(res)
}#plot_dlr
