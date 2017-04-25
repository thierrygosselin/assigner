#' @name plot_dlr
#' @title Assignment plot of genotype likelihood (Dlr).
#' @description The function generate a figure similar to Paetkau's et al. (2004)
#' Fig 6.
#' @param data The assigment object from 
#' \code{\link[assigner]{dlr}} output.
#' @param POPA First population to compare.
#' @param POPB Second population to compare (with A).
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

#' @importFrom stringi stri_join stri_sub stri_replace_all_fixed stri_detect_fixed stri_replace_na
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs sample_n sample_frac one_of
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid annotate geom_abline scale_shape_manual scale_fill_manual geom_abline aes_string geom_point scale_x_continuous scale_y_continuous
#' @importFrom readr read_table read_delim read_tsv
#' @importFrom purrr discard flatten_dbl

#' @examples
#' \dontrun{
#' # if you didn't run assigner::dlr before:
#' dlr.fig <- assigner::dlr(
#' data = "assignment.gdv", strata = "my.strata.tsv")$assignment %>%
#' assigner::plot_dlr(data = ., POPA = "mypop1", POPB = "mypop2")
#' 
#' # to get the plot:
#' dlr$assignment.plot
#' 
#' # If you already have the results of assigner::dlr:
#' assigner::plot_dlr(data = results, POPA = "mypop1", POPB = "mypop2")$assignment.plot
#' }

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
  POPA,
  POPB,
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
  if (missing(POPA)) stop("Missing POPA argument")
  if (missing(POPB)) stop("Missing POPB argument")
  
  # import data ----------------------------------------------------------------
  assignment.select <- dplyr::rename(data, Populations = POP_ID)
  
  # Check that POPA and POPB are found in the data
  if (!POPA %in% colnames(assignment.select)) stop("POPA value is not in the assignment file")
  if (!POPB %in% colnames(assignment.select)) stop("POPB value is not in the assignment file")
  
  pop.select <- c(POPA, POPB)
  fixed.header <- c("INDIVIDUALS", "Populations", "INFERRED", "LIK_MAX", "LIK_HOME", "LIK_RATIO")
  
  assignment.select <- suppressWarnings(assignment.select %>% 
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
    message("Saving was not selected...")
  } else {
    message("Saving the figure was selected, the filename: ", filename)
    
    if (is.null(plot.width)) stop("plot.width argument is required")
    if (is.null(plot.height)) stop("plot.height argument is required")
    if (is.null(plot.dpi)) stop("plot.dpi argument is required")
    ggplot2::ggsave(filename, width = plot.width, height = plot.height,
                    dpi = plot.dpi, units = "cm", useDingbats = FALSE)
  }
  
  res <- list(assignment.data.plot = assignment.select, assignment.plot = assignment.plot)
  return(res)
}#plot_dlr
