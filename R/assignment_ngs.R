# Assignment analysis of massive parallel sequencing data
#' @name assignment_ngs

#' @title Assignment analysis tailored for RADseq data

#' @description
#' The arguments in the \code{assignment_ngs} function were tailored for the
#' reality of GBS/RADseq data for assignment analysis while
#' maintaining a reproducible workflow. Assignment are conducted using
#' \href{https://github.com/eriqande/gsi_sim}{gsi_sim} or 
#' \code{\link[adegenet]{adegenet}}. 
#' 
#' \itemize{
#'   \item \strong{Input file:} various file format are supported (see \code{data} argument below)
#'   \item \strong{Filters:} genotypes, markers, individuals and populations can be 
#'   filtered and/or selected in several ways using blacklist,
#'   whitelist and other arguments (see details).
#'   \item \strong{Cross-Validations:} Markers can be randomly selected for a classic LOO (Leave-One-Out)
#'   assignment or chosen based on ranked Fst for a thl
#'   (Training, Holdout, Leave-one-out) assignment analysis
#'   \item \strong{Imputations:} Map-independent imputation of missing genotype/alleles
#'   using Random Forest or the most frequent category.
#'   \item \strong{Assignment analysis:} conducted in 
#'   \href{https://github.com/eriqande/gsi_sim}{gsi_sim}, a tool 
#'   for doing and simulating genetic stock identification and 
#'   developed by Eric C. Anderson, or 
#' \href{https://github.com/thibautjombart/adegenet}{adegenet}, 
#' an R package developed by Thibaul Jombart
#'   \item \strong{Parallel:} The assignment can be conduncted on multiple CPUs
#'   \item \strong{Results:} Assignment results in raw or processed tables and figures
#' }

#' @inheritParams radiator::tidy_genomic_data 

#' @param strata (optional/required) Optional for file format with population 
#' grouping integrated (e.g. vcf is not population-wise and requires a strata file).
#' 
#' The strata file is a tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The strata column is cleaned of a white spaces that interfere with some
#' packages or codes: space is changed to an underscore \code{_}.
#' Default: \code{strata = NULL}.

#' @param pop.levels (optional) Documented in 
#' \pkg{radiator} \code{\link[radiator]{read_strata}}.
#' Default: \code{pop.levels = NULL}.

#' @param assignment.analysis (character) Assignment analysis conducted with 
#' \code{assignment.analysis = "gsi_sim"} or 
#' \code{assignment.analysis = "adegenet"}.

#' @param sampling.method (character) Should the markers be randomly selected
#' \code{sampling.method == "random"} for a classic Leave-One-Out (LOO) assignment or
#' chosen based on ranked Fst \code{sampling.method == "ranked"}, used in a
#' Training-Holdout-Leave One Out (thl) assignment. 
#' \emph{See details}.

#' @param adegenet.dapc.opt (optional, character) \strong{Argument available only when 
#' using:
#' \code{assignment.analysis = "adegenet"} with
#' \code{sampling.method == "random"}}.
#' 
#' The assignment using dapc can use the optimized alpha score 
#' \code{adegenet.dapc.opt == "optim.a.score"} or 
#' cross-validation \code{adegenet.dapc.opt == "xval"}
#' for stability of group membership probabilities. 
#' For fine tuning the trade-off between power of discrimination and over-fitting.
#' See \pkg{adegenet} documentation for more details.
#' \code{adegenet.dapc.opt == "xval"} doesn't work with missing data, so it's 
#' only available with \strong{imputed data} (i.e. imputation.method == "rf" or "max").
#' With non imputed data or the default: \code{adegenet.dapc.opt == "optim.a.score"}.


#' @param adegenet.n.rep (optional, integer) 
#' When \code{adegenet.dapc.opt == "xval"}, the number of replicates to be 
#' carried out at each level of PC retention. 
#' Default: \code{adegenet.n.rep = 30}.
#' See \pkg{adegenet} documentation for more details.
#' @param adegenet.training (optional, numeric) 
#' When \code{adegenet.dapc.opt == "xval"}, the proportion of data (individuals) 
#' to be used for the training set.
#' Default: \code{adegenet.training = 0.9}, if all groups have >= 10 members. 
#' Otherwise, training.set scales automatically to the largest proportion 
#' that still ensures all groups will be present in both training 
#' and validation sets.
#' See \pkg{adegenet} documentation for more details.

#' @param thl (character, integer, proportion) For \code{sampling.method = "ranked"} only.
#' Default \code{thl = 1}, 1 individual sample is used as holdout. This individual is not
#' participating in the markers ranking. For each marker number,
#' the analysis will be repeated with all the indiviuals in the data set
#' (e.g. 500 individuals, 500 times 500, 1000, 2000 markers).
#' If a proportion is used e.g. \code{thl = 0.15}, 15 percent of individuals in each
#' populations are chosen randomly as holdout individuals.
#' With \code{thl = "all"} all individuals are used for ranking (not good) and
#' the argument \code{iteration.method = 1} is set by default.
#' For the other thl values, you can create different holdout individuals lists
#' with the \code{iteration.method} argument below (bootstrap).

#' @param iteration.method (integer) With random marker selection the iterations argument =
#' the number of iterations to repeat marker resampling. 
#' Default: \code{iteration.method = 10}.
#' With \code{marker.number = c(500, 1000)} and default iterations setting,
#' 500 markers will be randomly chosen 10 times and 1000 markers will be randomly
#' chosen 10 times.
#' 
#' \strong{Notes:} If all the markers are used, with \code{marker.number = "all"}
#' or in a series of marker number groupings \code{marker.number = c(200, 500, "all")}, 
#' the number of iteration is automatically set to 1. The remaining groupings
#' are unaffected.
#' 
#' For the ranked method, using \code{thl = 1}, the analysis
#' will be repeated for each individuals in the data set for every
#' \code{marker.number} selected. With a proportion argument \code{thl = 0.15},
#' 15 percent of individuals in each populations are chosen randomly as holdout
#' individuals and this process is reapeated the number of times chosen by the
#' \code{iteration.method} value.

#' @param subsample (Integer or Character, optional) 
#' With \code{subsample = 36}, 36 individuals in each populations are chosen
#' randomly to represent the dataset. This integer as to be smaller than the
#' population with min sample size, if higher, the minimum sample size found 
#' in the data will be used as default. In doubt, use \code{subsample = "min"},
#' the function will use the smallest population sample size found in the data.
#' Default: \code{subsample = NULL} (no subsampling).

#' @param iteration.subsample (Integer) The number of iterations to repeat 
#' subsampling.
#' With \code{subsample = 20} and \code{iteration.subsample = 10},
#' 20 individuals/populations will be randomly chosen 10 times.
#' Default: \code{iteration.subsample = 1}.

#' @param marker.number (Integer or string of number or "all") Calculations with
#' fixed or subsample of your markers.
#' e.g. To test 500, 1000, 2000 and all  the markers:
#' \code{marker.number = c(500, 1000, 2000, "all")}.
#' To use only 500 makers \code{marker.number = 500}.
#' Default = \code{marker.number = "all"}.

#' @param folder (optional) The name of the folder created in the working 
#' directory to save the files/results. Default: \code{folder = NULL} will create
#' the folder for you with data and time.

#' @param filename (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. 
#' Several info will be appended to the name of the file.
#' Default \code{assignment_data.txt}.

#' @param keep.gsi.files (logical, optional) With the default, 
#' the intermediate input and output gsi_sim files will be deleted from the 
#' directory when finished processing. I you decide to keep the files, with 
#' \code{keep.gsi.files = TRUE}, remember to allocate a large chunk of the disk 
#' space for the analysis.
#' Default: \code{keep.gsi.files = FALSE} 

#' @param random.seed (integer, optional) For reproducibility, set an integer
#' that will be used inside function that requires randomness. With default,
#' a random number is generated and printed in the appropriate output.
#' Default: \code{random.seed = NULL}.

#' @param ... (optional) To pass further argument for fine-tuning the 
#' function (filters). See details.

#' @details 
#' \strong{Input files:} see \pkg{radiator} \code{\link[radiator]{tidy_genomic_data}}
#' for detailed information about supported file format.
#' 
#' \strong{Available filters:}
#' 
#' Further arguments can be passed via the \emph{dots-dots-dots}: 
#' \itemize{
#' \item whitelist.markers
#' }
#' For argument documentation see \pkg{radiator} \code{\link[radiator]{tidy_genomic_data}}.
#' 
#' 
#' \strong{Imputations:}
#' 
#' The imputations using Random Forest requires more time to compute
#' and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals
#' will require 15 min. This is using multiple CPUs. To have your computer ready
#' for parallel computing during imputations follow the steps in the
# \href{https://github.com/thierrygosselin/radiator/blob/master/vignettes/vignette_imputations_parallel.Rmd}{vignette}
#' (~10 min)
#' 
#' \strong{THL, Ranking and Fst:}
#' 
#' With \code{sampling.method = "ranked"}, the markers are first 
#' arranged by \emph{decreasing} values of Fst.
#' The Fst is computed with \code{\link{fst_WC84}} function, that uses a fast 
#' implementation of Weir and Cockerham 1984 Fst/Theta equations. 

#' @return Depending on arguments selected, several files are written to the your
#' working directory or \code{folder}
#' The output in your global environment is a list. To view the assignment results
#' \code{$assignment} to view the ggplot2 figure \code{$plot.assignment}. 
#' See example below.

#' @note \code{assignment_ngs} assumes that the command line version of 
#' \href{https://github.com/eriqande/gsi_sim}{gsi_sim} 
#' is properly installed into \code{file.path(system.file(package = "assigner"), "bin", "gsi_sim")}.
#' Things are set up so that it will try running gsi_sim, and if it does not find it, the 
#' program will throw an error and ask the user to run \code{\link{install_gsi_sim}}
#' which will do its best to put a usable copy of gsi_sim where it is needed. 
#' To do so, you must be connected to the internet. If that doesn't work, you will
#' need to compile the program yourself, or get it yourself, and the manually copy
#' it to \code{file.path(system.file(package = "assigner"), "bin", "gsi_sim")}.
#' To compile \href{https://github.com/eriqande/gsi_sim}{gsi_sim}, follow the 
#' instruction here: \url{https://github.com/eriqande/gsi_sim}.

#' @export
#' @rdname assignment_ngs
#' @importFrom parallel detectCores
#' @importFrom stringi stri_join stri_sub stri_replace_all_fixed stri_detect_fixed stri_replace_na
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs sample_n sample_frac
#' @importFrom radiator tidy_genomic_data change_pop_names write_genind filter_maf detect_genomic_format filter_common_markers filter_monomorphic
#' @importFrom stats var median quantile
#' @importFrom purrr map flatten keep discard
#' @importFrom adegenet genind
#' @importFrom readr read_delim read_tsv
#' @importFrom adegenet optim.a.score dapc xvalDapc indNames pop predict.dapc
#' @importFrom tibble as_data_frame data_frame
#' @importFrom tidyr separate
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_summary labs theme element_blank element_text geom_jitter scale_colour_manual scale_y_reverse theme_light geom_bar facet_grid


#' @examples
#' \dontrun{
#' assignment.treefrog <- assignment_ngs(
#'     data = "batch_1.vcf", strata = "strata.treefrog.tsv",
#'     assignment.analysis = "gsi_sim",
#'      marker.number = c(500, 5000, "all"),
#'      sampling.method = "ranked", thl = 0.3,
#'      subsample = 25, iteration.subsample = 10)
#' 
#' Since the 'folder' argument is missing, it will be created automatically
#' inside your working directory.
#' 
#' To create a dataframe with the assignment results: 
#' assignment <- assignment.treefrog$assignment.
#' 
#' To plot the assignment using ggplot2 and facet 
#' (with subsample by current pop):
#' assignment.treefrog$plot.assignment + ggplot2::facet_grid(SUBSAMPLE~CURRENT).
#' 
#' To view the full range of y values = Assignment success(%): 
#' assignment.treefrog$plot.assignment + 
#' ggplot2::facet_grid(SUBSAMPLE~CURRENT) + 
#' ggplot2::scale_y_continuous(limits = c(0,100)) 
#' To save the plot:
#' ggplot2::ggsave("assignment.pdf", height = 35, width = 60,dpi = 600,
#' units = "cm", useDingbats = FALSE)
#' 
#' # If you want to remove underscore in population names that contained white space:
#' facet_names <- c(
#' `some_pop` = "Some POP", 
#' `some_other_pop` = "This is what I want", 
#' `OVERALL` = "Overall")
#' 
#' # use the labeller in the facet_grid or facet_wrap call:
#' assignment.treefrog$plot.assignment + 
#' ggplot2::facet_grid(~CURRENT, ggplot2::labeller = ggplot2::as_labeller(facet_names)) + 
#' ggplot2::scale_y_continuous(limits = c(0,100)) 
#' figure # this one should be ok with custom naming in the facets.
#' }


#' @references Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008)
#' An improved method for predicting the accuracy of genetic stock identification.
#' Canadian Journal of Fisheries and Aquatic Sciences 65, 7:1475-1486.
#' @references Anderson, E. C. (2010) Assessing the power of informative subsets of
#' loci for population assignment: standard methods are upwardly biased.
#' Molecular ecology resources 10, 4:701-710.
#' @references Weir BS, Cockerham CC (1984) Estimating F-Statistics for the
#' Analysis of Population Structure. Evolution, 38, 1358–1370.
#' @references Jombart T, Devillard S, Balloux F. 
#' Discriminant analysis of principal components: a new method for the analysis 
#' of genetically structured populations. 
#' BMC Genet. 2010:11: 94. doi:10.1186/1471-2156-11-94
#' @references Jombart T, Ahmed I. adegenet 1.3-1: new tools for the analysis 
#' of genome-wide SNP data. 
#' Bioinformatics. 2011:27: 3070–3071. doi:10.1093/bioinformatics/btr521

#' @seealso \href{https://github.com/eriqande/gsi_sim}{gsi_sim} and 
#' \href{https://github.com/eriqande/rubias}{rubias}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and Eric C. Anderson

assignment_ngs <- function(
  data,
  assignment.analysis = c("gsim_sim", "adegenet"),
  sampling.method = c("ranked", "random"),
  adegenet.dapc.opt = "optim.a.score",
  adegenet.n.rep = 30,
  adegenet.training = 0.9,
  thl = 1,
  iteration.method = 10,
  subsample = NULL,
  iteration.subsample = 1,
  marker.number = "all",
  strata = NULL,
  pop.levels = NULL,
  verbose = FALSE,
  folder = NULL,
  filename = "assignment_data.txt",
  keep.gsi.files = FALSE,
  random.seed = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  
  ## testing 
  # whitelist.markers = NULL
  
  
  cat("#######################################################################\n")
  cat("###################### assigner::assignment_ngs #######################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  res.list <- list() # results to keep stored in this list
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  if (missing(assignment.analysis)) stop("assignment.analysis argument missing")
  if (missing(sampling.method)) stop("sampling.method argument missing")
  if (assignment.analysis == "gsi_sim" & !gsi_sim_exists()) {
    stop("Can't find the gsi_sim executable where it was expected at ", gsi_sim_binary_path(), ".  
         If you have internet access, you can install it
         from within R by invoking the function \"install_gsi_sim(fromSource = TRUE)\"")
  }
  
  assignment.analysis <- match.arg(
    arg = assignment.analysis, 
    choices = c("gsi_sim", "adegenet"), 
    several.ok = FALSE)
  
  sampling.method <- match.arg(
    arg = sampling.method, 
    choices = c("ranked", "random"), 
    several.ok = FALSE) 
  
  if (assignment.analysis == "gsi_sim") message("Assignment analysis with gsi_sim")
  if (assignment.analysis == "adegenet") message("Assignment analysis with adegenet")
  
  # Correct iteration.method default if using only "all" in marker.number
  if (length(marker.number) == 1 && marker.number == "all" && sampling.method == "random") {
    iteration.method <- 1
  }
  
  if ("all" %in% marker.number && sampling.method == "random") {
    manage.all <- TRUE
  } else {
    manage.all <- FALSE
  }
  
  if (thl == "all") {
    message("Note: with thl == \"all\", automatically setting iteration.method = 1")
    iteration.method <- 1
  }
  
  # dotslist -------------------------------------------------------------------
  dotslist <- list(...)
  
  want <- c("whitelist.markers")
  unknowned_param <- setdiff(names(dotslist), want)

  if (length(unknowned_param) > 0) {
    stop("Unknowned \"...\" parameters ",
         stringi::stri_join(unknowned_param, collapse = " "))
  }
  assigner.dots <- dotslist[names(dotslist) %in% want]

  whitelist.markers <- assigner.dots[["whitelist.markers"]]
  
  # POP levels -----------------------------------------------------------------
  
  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...
  # pop.levels, pop.labels, pop.select, strata, etc
  # if (!is.null(pop.levels) & is.null(pop.labels)) {
  #   pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  #   pop.labels <- pop.levels
  # }
  # if (!is.null(pop.labels)) {
  #   pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  # }
  # if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  # if (!is.null(pop.select)) {
  #   pop.select <- stringi::stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  # }
  
  # store function call
  res.list$call <- match.call()
  
  # Create a folder based on filename to save the output files -----------------
  if (is.null(folder)) {
    # Get date and time to have unique filenaming
    file.date <- format(Sys.time(), "%Y%m%d@%H%M")
    imputation.method <- NULL
    # if (is.null(imputation.method)) {
      message("Map-independent imputations: no")
      directory <- stringi::stri_join(getwd(),"/", "assignment_analysis_", "method_", sampling.method, "_no_imputations_", file.date, "/", sep = "")
      dir.create(file.path(directory))
    # } else {
    #   message("Map-independent imputations: yes")
    #   directory <- stringi::stri_join(getwd(),"/","assignment_analysis_", "method_", sampling.method, "_imputations_", imputation.method,"_", hierarchical.levels, "_", file.date, "/", sep = "")
    #   dir.create(file.path(directory))
    # }
    message("Folder: ", directory)
    file.date <- NULL #unused object
  } else {
    directory <- stringi::stri_join(getwd(), "/", folder, "/", sep = "")
    dir.create(file.path(directory))
    message("Folder: ", directory)
  }
  
  # File type detection --------------------------------------------------------
  data.type <- radiator::detect_genomic_format(data)
  
  # Strata argument required for VCF and haplotypes files ----------------------
  if (data.type == "haplo.file" | data.type == "vcf.file") {
    if (is.null(strata)) stop("strata argument is required")
  }
  
  # Import input ---------------------------------------------------------------
  input <- radiator::tidy_genomic_data(
    data = data, 
    whitelist.markers = whitelist.markers, 
    strata = strata,
    filename = NULL,
    verbose = FALSE
  )
  
  # need to remove space in POP_ID name to work in gsi_sim
  input$POP_ID <- stringi::stri_replace_all_fixed(input$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)
  
  # input <- radiator::change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)
  pop.levels <- pop.labels <- NULL
  input <- radiator::change_pop_names(data = input, pop.levels = unique(pop.labels), pop.labels = unique(pop.labels))
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
  # subsampling data------------------------------------------------------------
  # create the subsampling list
  ind.pop.df <-  dplyr::distinct(.data = input, POP_ID, INDIVIDUALS)
  
  if (!is.null(subsample)) {
    min.pop.n <- min(dplyr::count(x = ind.pop.df, POP_ID, sort = TRUE)$n)
    # Control the subsample argument, replace if > than min sample size
    if (rlang::is_bare_numeric(subsample)) {
      if (subsample <= min.pop.n) {
        subsample <- subsample
      } else {
        message("Warning: value of subsample argument is higher than the min sample size found in the data")
        message("Replacing subsample value by: ", min.pop.n)
        subsample <- min.pop.n
      }
    } else {
      # replace min by the min sample size found in the data
      if (subsample == "min") {
        subsample <- min.pop.n
        message("Using subsample size of: ", subsample)
      } else {
        stop("Wrong subsample value")
      }
    }
  }
  
  subsample.list <- purrr::map(
    .x = 1:iteration.subsample,
    .f = subsampling_data,
    ind.pop.df = ind.pop.df,
    subsample = subsample,
    random.seed = random.seed
  )
  
  # keep track of subsampling individuals and write to directory
  if (is.null(subsample)) {
    message("Subsampling: not selected")
  } else {
    message("Subsampling: selected")
    subsampling.individuals <- dplyr::bind_rows(subsample.list)
    readr::write_tsv(
      x = subsampling.individuals, 
      path = paste0(directory, "subsampling_individuals.tsv"), 
      col_names = TRUE, 
      append = FALSE
    )
    res.list$subsampling.individuals <- subsampling.individuals
  } # End subsampling
  
  # unused objects
  subsampling.individuals <- ind.pop.df <- NULL
  
  # assignment analysis --------------------------------------------------------
  res <- purrr::map(
    .x = subsample.list,
    .f = assignment_function,
    input = input,
    subsample = subsample,
    assignment.analysis = assignment.analysis,
    marker.number = marker.number,
    sampling.method = sampling.method,
    iteration.method = iteration.method,
    filename = filename,
    directory = directory,
    keep.gsi.files = keep.gsi.files,
    verbose = verbose,
    parallel.core = parallel.core,
    manage.all = manage.all,
    random.seed = random.seed,
    res.list = res.list,
    thl = thl,
    adegenet.dapc.opt = adegenet.dapc.opt,
    adegenet.n.rep = adegenet.n.rep,
    adegenet.training = adegenet.training
  )
  res <- dplyr::bind_rows(res)
  
  # if (is.null(imputation.method)) {
    filename.res <- stringi::stri_join("assignment", sampling.method, "no.imputation.results.summary.stats.subsample", "tsv", sep = ".")
  # } else {# with imputations
  #   filename.res <- stringi::stri_join("assignment", sampling.method, "imputed.results,summary.stats.subsample", "tsv", sep = ".")
  # }
  readr::write_tsv(
    x = res, path = paste0(directory,filename.res),
    col_names = TRUE, append = FALSE
  )
  
  # Summary of the subsampling iterations
  if (iteration.subsample > 1) {
    res.pop <- dplyr::filter(.data = res, CURRENT != "OVERALL") %>%
      dplyr::group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
      dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
      dplyr::summarise(
        MEAN = round(mean(ASSIGNMENT_PERC), 2),
        SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
        MIN = round(min(ASSIGNMENT_PERC), 2),
        MAX = round(max(ASSIGNMENT_PERC), 2),
        MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
        QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
        QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
      ) %>% 
      dplyr::mutate(SUBSAMPLE = rep("OVERALL", n())) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(CURRENT, MARKER_NUMBER)
    
    res.overall <- res.pop %>% 
      dplyr::group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
      dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
      dplyr::summarise(
        MEAN = round(mean(ASSIGNMENT_PERC), 2),
        SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
        MIN = round(min(ASSIGNMENT_PERC), 2),
        MAX = round(max(ASSIGNMENT_PERC), 2),
        MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
        QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
        QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
      ) %>%
      dplyr::mutate(
        CURRENT = rep("OVERALL", n()),
        SUBSAMPLE = rep("OVERALL", n())
      ) %>% 
      dplyr::ungroup(.) %>% 
      dplyr::arrange(CURRENT, MARKER_NUMBER)
    
    res.pop.overall <- suppressWarnings(
      dplyr::bind_rows(res.pop, res.overall) %>%
        dplyr::mutate(CURRENT = factor(CURRENT, levels = levels(res.pop$CURRENT), ordered = TRUE)) %>%
        dplyr::arrange(CURRENT, MARKER_NUMBER) %>%
        dplyr::mutate(
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE
        ) %>%
        dplyr::select(CURRENT, MARKER_NUMBER, MEAN, MEDIAN, SE, MIN, MAX, QUANTILE25, QUANTILE75, SE_MIN, SE_MAX, METHOD, MISSING_DATA, SUBSAMPLE)
    )
    
    res <- dplyr::bind_rows(
      res %>% 
        dplyr::mutate(SUBSAMPLE = as.character(SUBSAMPLE)), 
      res.pop.overall) %>% 
      dplyr::mutate(SUBSAMPLE = factor(SUBSAMPLE, levels = c(1:iteration.subsample, "OVERALL"), ordered = TRUE)) %>% 
      dplyr::arrange(CURRENT, MARKER_NUMBER, SUBSAMPLE)
    
    # if (is.null(imputation.method)) {
      filename.assignment.sum.subsample <- stringi::stri_join("assignment", sampling.method, "no.imputation.results.summary.stats.subsample.overall", "tsv", sep = ".")
    # } else {# with imputations
    #   filename.assignment.sum.subsample <- stringi::stri_join("assignment", sampling.method, "imputed.results.summary.stats.subsample.overall", "tsv", sep = ".")
    # }
    readr::write_tsv(
      x = res, path = paste0(directory, filename.assignment.sum.subsample),
      col_names = TRUE, append = FALSE
    )
    
    # unused objects
    res.pop.overall <- res.overall <- res.pop <- NULL
  } # End summary of the subsampling iterations
  
  # Assignment plot ------------------------------------------------------------
  # if (is.null(imputation.method)) { # no imputation
    plot.assignment <- ggplot2::ggplot(res, ggplot2::aes(x = factor(MARKER_NUMBER), y = MEAN)) +
      ggplot2::geom_point(size = 2, alpha = 0.5, na.rm = TRUE) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = SE_MIN, ymax = SE_MAX), width = 0.3, na.rm = TRUE) +
      ggplot2::scale_y_continuous(breaks = c(0, 10, 20 ,30, 40, 50, 60, 70, 80, 90, 100)) +
      ggplot2::labs(x = "Marker number",
                    y = "Assignment success (%)") +
      ggplot2::theme(
        legend.position = "bottom",      
        panel.grid.minor.x = ggplot2::element_blank(), 
        panel.grid.major.y = ggplot2::element_line(colour = "grey60", linetype = "dashed"), 
        axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", face = "bold", angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
      ) +
      ggplot2::theme_bw()
  # } else {#with imputations
  #   
  #   plot.assignment <- ggplot2::ggplot(res, ggplot2::aes(x = factor(MARKER_NUMBER), y = MEAN)) +
  #     ggplot2::geom_point(ggplot2::aes(colour = MISSING_DATA), size = 2, alpha = 0.8) +
  #     ggplot2::geom_errorbar(ggplot2::aes(ymin = SE_MIN, ymax = SE_MAX), width = 0.3) +
  #     ggplot2::scale_colour_manual(name = "Missing data", values = c("gray33", "dodgerblue")) +
  #     ggplot2::scale_y_continuous(breaks = c(0, 10, 20 ,30, 40, 50, 60, 70, 80, 90, 100)) +
  #     ggplot2::labs(x = "Marker number",
  #                   y = "Assignment success (%)") +
  #     ggplot2::theme(
  #       legend.position = "bottom",      
  #       panel.grid.minor.x = ggplot2::element_blank(), 
  #       panel.grid.major.y = ggplot2::element_line(colour = "grey60", linetype = "dashed"), 
  #       axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
  #       axis.text.x = ggplot2::element_text(size = 8, family = "Helvetica", face = "bold", angle = 90, hjust = 1, vjust = 0.5), 
  #       axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"), 
  #       axis.text.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
  #     ) +
  #     ggplot2::theme_bw()
  # } # end plot
  
  # results --------------------------------------------------------------------
  res.list$assignment <- res
  res.list$plot.assignment <- plot.assignment
  timing <- proc.time() - timing
  message("Computation time: ", round(timing[[3]]), " sec")
  cat("############################## completed ##############################\n")
  return(res.list)
} # End assignment_ngs

# Internal Nested Functions -----------------------------------------------------------

# assignment_gsi_sim------------------------------------------------------------
#' @title assignment with gsi_sim
#' @description assignment with gsi_sim
#' @rdname assignment_gsi_sim
#' @export
#' @keywords internal

assignment_gsi_sim <- function(
  input = NULL,
  strata.df = NULL,
  select.markers = NULL,
  markers.names = NULL,
  missing.data = NULL,
  i = NULL,
  m = NULL,
  holdout = NULL,
  filename = NULL,
  directory.subsample = NULL,
  keep.gsi.files = FALSE,
  pop.labels = NULL,
  sampling.method = "random",
  thl = 1
) {
  # missing.data <- "no.imputation" #test
  data.select <- suppressWarnings(
    dplyr::semi_join(input, select.markers, by = "MARKERS") %>%
      dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  )
  
  # Write gsi_sim input file to directory
  input.gsi <- radiator::write_gsi_sim(
    data = data.select, 
    pop.levels = unique(pop.labels), 
    pop.labels = unique(pop.labels), 
    strata = NULL, 
    filename = filename
  )
  
  # Run gsi_sim ------------------------------------------------------------
  output.gsi <- stringi::stri_replace_all_fixed(input.gsi, pattern = "txt", replacement = "output.txt")
  setwd(directory.subsample)
  system(paste(gsi_sim_binary(), "-b", input.gsi, "--self-assign > ", output.gsi))
  
  # Option remove the input file from directory to save space
  if (!keep.gsi.files) file.remove(input.gsi)
  
  # Get Assignment results -------------------------------------------------
  
  assignment <- suppressWarnings(
    suppressMessages(readr::read_delim(output.gsi, col_names = "ID", delim = "\t")) %>%
      tidyr::separate(ID, c("KEEPER", "ASSIGN"), sep = ":/", extra = "warn") %>%
      dplyr::filter(KEEPER == "SELF_ASSIGN_A_LA_GC_CSV") %>%
      tidyr::separate(ASSIGN, c("INDIVIDUALS", "ASSIGN"), sep = ";", extra = "merge") %>%
      tidyr::separate(ASSIGN, c("INFERRED", "OTHERS"), sep = ";", convert = TRUE, numerals = "no.loss", extra = "merge") %>%
      tidyr::separate(OTHERS, c("SCORE", "OTHERS"), sep = ";;", convert = TRUE, numerals = "no.loss", extra = "merge") %>%
      tidyr::separate(OTHERS, c("SECOND_BEST_POP", "OTHERS"), sep = ";", convert = TRUE, numerals = "no.loss", extra = "merge") %>%
      tidyr::separate(OTHERS, c("SECOND_BEST_SCORE", "OTHERS"), sep = ";;", convert = TRUE, numerals = "no.loss")
  )
  assignment <- suppressWarnings(
    dplyr::mutate(.data = assignment, INDIVIDUALS = as.character(INDIVIDUALS)) %>% 
      dplyr::left_join(strata.df, by = "INDIVIDUALS") %>%
      dplyr::rename(CURRENT = POP_ID) %>% 
      dplyr::mutate(
        INFERRED = factor(INFERRED, levels = unique(pop.labels), ordered = TRUE),
        INFERRED = droplevels(INFERRED),
        SECOND_BEST_POP = factor(SECOND_BEST_POP, levels = unique(pop.labels), ordered = TRUE),
        SECOND_BEST_POP = droplevels(SECOND_BEST_POP),
        SCORE = round(SCORE, 2),
        SECOND_BEST_SCORE = round(SECOND_BEST_SCORE, 2),
        MARKER_NUMBER = as.numeric(rep(m, n())), # m: Number of markers
        METHOD = rep(sampling.method, n()),
        MISSING_DATA = rep(missing.data, n())
      ) %>%
      dplyr::select(INDIVIDUALS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
      dplyr::arrange(CURRENT)
  )
  
  if (sampling.method == "random") {
    assignment <- dplyr::mutate(.data = assignment, ITERATIONS = rep(i, n())) %>% 
      dplyr::select(INDIVIDUALS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA, ITERATIONS) %>%
      dplyr::arrange(CURRENT)
  }
  
  if (sampling.method == "ranked" & thl != "all") {
    assignment <- dplyr::filter(.data = assignment, INDIVIDUALS %in% holdout$INDIVIDUALS)
  }
  
  # Option remove the output file from directory to save space
  if (!keep.gsi.files) file.remove(output.gsi)
  
  # saving preliminary results
  if (sampling.method == "ranked") {
    
    assignment <- assignment %>%
      dplyr::mutate(METHOD = rep(stringi::stri_join("ranked_thl_", thl) , n()))
    
    if (thl != 1 & thl != "all") {
      assignment <- dplyr::mutate(
        .data = assignment,
        CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
        CURRENT = droplevels(CURRENT),
        ITERATIONS = rep(i, n())
      )
    }
  }
  return(assignment)
} # End assignment_gsi_sim function

# assignment_adegenet-----------------------------------------------------------
#' @title assignment with adegenet
#' @description assignment with adegenet
#' @rdname assignment_adegenet
#' @export
#' @keywords internal

assignment_adegenet <- function(
  data = NULL,
  select.markers = NULL,
  markers.names = NULL,
  missing.data = NULL,
  adegenet.dapc.opt = "optim.a.score",
  adegenet.n.rep = 30,
  adegenet.training = 0.9, 
  i = NULL,
  m = NULL,
  holdout = NULL,
  sampling.method = "random",
  thl = 1,
  subsample.id = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  # data <- genind.object #test
  # missing.data <- "no.imputation" #test
  data.select <- data[loc = select.markers$MARKERS]
  
  # Run adegenet -----------------------------------------------------------
  pop.data <- data.select@pop
  pop.data <- droplevels(pop.data)
  
  
  if (sampling.method == "random") {
    # DAPC optimized a-score 
    if (adegenet.dapc.opt == "optim.a.score") {
      dapc.best.optim.a.score <- adegenet::optim.a.score(
        adegenet::dapc(data.select, 
                       n.da = length(levels(pop.data)), 
                       n.pca = round((length(adegenet::indNames(data.select))/3) - 1, 0)), 
        pop = pop.data, 
        plot = FALSE
      )$best
      message("a-score optimisation for iteration: ", i) # message not working in parallel...
      
      # DAPC
      dapc.assignment <- adegenet::dapc(data.select, n.da = length(levels(pop.data)), n.pca = dapc.best.optim.a.score, pop = pop.data)
      message("DAPC iteration: ", i)
      message("DAPC marker group: ", m)
    }
    
    # DAPC with Cross-Validation
    if (adegenet.dapc.opt == "xval") {
      dapc.assignment <- adegenet::xvalDapc(
        x = data.select@tab, 
        grp = adegenet::pop(data.select),
        n.da = length(levels(pop.data)),
        n.pca.max = round((length(adegenet::indNames(data.select))/3) - 1, 0), 
        n.rep = adegenet.n.rep , 
        training.set = adegenet.training, 
        result = "groupMean", 
        center = TRUE, 
        scale = FALSE, 
        xval.plot = FALSE, 
        parallel = "multicore", 
        ncpus = parallel.core
      )$DAPC
      
      message("DAPC iteration: ", i)
      message("DAPC marker group: ", m)
    }
  }
  
  if (sampling.method == "ranked") {
    
    # Alpha-Score DAPC training data
    training.data <- data.select[!adegenet::indNames(data.select) %in% holdout$INDIVIDUALS] # training dataset
    pop.training <- training.data@pop
    pop.training <- droplevels(pop.training)
    
    dapc.best.optim.a.score <- adegenet::optim.a.score(
      adegenet::dapc(training.data, 
                     n.da = length(levels(pop.training)),
                     n.pca = round(((length(adegenet::indNames(training.data))/3) - 1), 0)),
      pop = pop.training,
      plot = FALSE
    )$best
    message("a-score optimisation for iteration: ", i)
    
    dapc.training <- adegenet::dapc(training.data,
                                    n.da = length(levels(pop.training)),
                                    n.pca = dapc.best.optim.a.score,
                                    pop = pop.training)
    message("DAPC of training data set for iteration: ", i)
    
    # DAPC holdout individuals
    holdout.data <- data.select[adegenet::indNames(data.select) %in% holdout$INDIVIDUALS] # holdout dataset
    pop.holdout <- holdout.data@pop
    pop.holdout <- droplevels(pop.holdout)
    assignment.levels <- levels(pop.holdout) # for figure
    rev.assignment.levels <- rev(assignment.levels)  # for figure 
    
    dapc.predict.holdout <- adegenet::predict.dapc(dapc.training, newdata = holdout.data)
    message("Assigning holdout data for iteration: ", i)
  }
  
  
  # Get Assignment results -------------------------------------------------
  if (sampling.method == "random") {
    assignment <- tibble::data_frame(ASSIGNMENT_PERC = summary(dapc.assignment)$assign.per.pop*100) %>% 
      dplyr::bind_cols(tibble::data_frame(POP_ID = levels(pop.data))) %>%
      dplyr::mutate(ASSIGNMENT_PERC = round(ASSIGNMENT_PERC, 2)) %>% 
      dplyr::select(POP_ID, ASSIGNMENT_PERC)
  }        
  if (sampling.method == "ranked") {
    assignment <- data.frame(
      INDIVIDUALS = adegenet::indNames(holdout.data), 
      CURRENT = pop.holdout,
      INFERRED = dapc.predict.holdout$assign, dapc.predict.holdout$posterior
    ) %>%
      dplyr::mutate(
        CURRENT = factor(CURRENT, levels = rev.assignment.levels, ordered = TRUE),
        INFERRED = factor(INFERRED, levels = assignment.levels, ordered = TRUE)
      )
  }
  
  assignment <- dplyr::mutate(.data = assignment,
                              METHOD = rep(sampling.method, n()),
                              ITERATIONS = rep(i, n()),
                              MARKER_NUMBER = as.numeric(rep(m, n())),
                              MISSING_DATA = rep(missing.data, n()))
  
  return(assignment)
} # End assignment_adegenet function


# marker_selection-----------------------------------------------------------
#' @title marker_selection 
#' @description Random selection of marker function + iteration.method
#' @rdname marker_selection
#' @export
#' @keywords internal
marker_selection <- function(
  iteration.method,
  m = NULL,
  random.seed = NULL,
  unique.markers = NULL
) {
  m <- as.numeric(m)
  
  # Set seed for random sampling
  if (is.null(random.seed)) {
    random.seed <- sample(x = 1:1000000, size = 1)
    set.seed(random.seed)
  } else {
    set.seed(random.seed)
  }
  
  select.markers <- dplyr::sample_n(tbl = unique.markers, size = m, replace = FALSE) %>%
    dplyr::arrange(MARKERS) %>%
    dplyr::mutate(
      ITERATIONS = rep(iteration.method, n()),
      MARKER_NUMBER = rep(m, n()),
      RANDOM_SEED_NUMBER = rep(random.seed, n())
    )
  return(select.markers)
}#End marker_selection



# assignment_random--------------------------------------------------------------
#' @title assignment_random
#' @description assignment_random
#' @rdname assignment_random
#' @export
#' @keywords internal

assignment_random <- function(
  x, # the list of dataframe with random markers previously marker.random.list
  assignment.analysis = "gsi_sim",
  input = NULL,
  genind.object = NULL,
  strata.df = NULL,
  directory.subsample = NULL,
  keep.gsi.files = FALSE,
  sampling.method = "random",
  subsample.id = NULL,
  filename = NULL,
  adegenet.n.rep = 30,
  adegenet.training = 0.9,
  holdout = NULL
) {
  x <- tibble::as_data_frame(x)
  # x <- marker.random.list[[1]]
  i <- as.integer(unique(x$ITERATIONS))      # iteration
  m <- as.numeric(unique(x$MARKER_NUMBER))   # number of marker selected
  
  select.markers <- dplyr::ungroup(x) %>% 
    dplyr::select(MARKERS) %>% 
    dplyr::arrange(MARKERS)
  
  # get the list of loci after filter
  markers.names <- unique(select.markers$MARKERS)
  
  # if (is.null(imputation.method)) {
    filename.imp <- "no_imputation.txt"
    missing.data <- "no.imputation"
  # } else {
  #   filename.imp <- "imputed.txt"
  #   
  #   if (imputation.method == "rf") {
  #     if (hierarchical.levels == "populations") {
  #       missing.data <- "imputed RF populations"
  #     } else {
  #       missing.data <- "imputed RF global"
  #     }
  #   } else {
  #     if (hierarchical.levels == "populations") {
  #       missing.data <- "imputed max populations"
  #     } else {
  #       missing.data <- "imputed max global"
  #     }
  #   }
  # }
  
  # Modify filename
  filename <- stringi::stri_join(directory.subsample, filename, sep = "")
  filename <- stringi::stri_replace_all_fixed(
    filename,
    pattern = ".txt",
    replacement = stringi::stri_join(
      "", "iteration", i, "markers", m, 
      filename.imp, sep = "_"
    )
  )
  
  
  if (assignment.analysis == "gsi_sim") {
    assignment <- assignment_gsi_sim(
      input = input,
      strata.df = strata.df,
      select.markers = select.markers,
      markers.names = markers.names,
      missing.data = missing.data, 
      i = i, 
      m = m,
      holdout = NULL,
      filename = filename,
      directory.subsample = directory.subsample,
      keep.gsi.files = keep.gsi.files,
      pop.labels = pop.labels,
      sampling.method = sampling.method
    )
  }
  if (assignment.analysis == "adegenet") {
    assignment <- assignment_adegenet(
      data = genind.object,
      select.markers = select.markers,
      adegenet.dapc.opt = "optim.a.score",
      adegenet.n.rep = adegenet.n.rep, 
      adegenet.training = adegenet.training,
      markers.names = markers.names,
      missing.data = "no.imputation", 
      i = i, 
      m = m,
      holdout = NULL
    )
  }
  
  # unused objects
  x <- i <- m <- select.markers <- markers.names <- filename <- NULL
  filename.imp <- missing.data <- NULL
  return(assignment)
}#End assignment_random

# assignment_marker_loop--------------------------------------------------------
#' @title assignment_marker_loop
#' @description assignment_marker_loop
#' @rdname assignment_marker_loop
#' @export
#' @keywords internal

assignment_marker_loop <- function(
  x, # marker.number
  assignment.analysis = "gsi_sim",
  fst.ranked = NULL,
  i = NULL,
  input = NULL,
  genind.object = NULL,
  strata.df = NULL,
  sampling.method = NULL,
  subsample.id = NULL,
  holdout = NULL,
  thl = 1,
  adegenet.dapc.opt = NULL,
  adegenet.n.rep = NULL,
  adegenet.training = NULL,
  directory.subsample = NULL,
  filename = NULL,
  keep.gsi.files = FALSE
) {
  # x <- marker.number[1]
  x <- as.numeric(x)
  message("Marker number: ", x)
  
  select.markers <- dplyr::filter(.data = fst.ranked, RANKING <= x) %>%
    dplyr::select(MARKERS)
  
  # get the list of markers after filter
  markers.names <- unique(select.markers$MARKERS)
  
  # if (is.null(imputation.method)) {
    filename.imp <- "no_imputation.txt"
    missing.data <- "no.imputation"
  # } else {
  #   filename.imp <- "imputed.txt"
  #   
  #   if (imputation.method == "rf") {
  #     if (hierarchical.levels == "populations") {
  #       missing.data <- "imputed RF populations"
  #     } else {
  #       missing.data <- "imputed RF global"
  #     }
  #   } else {
  #     if (hierarchical.levels == "populations") {
  #       missing.data <- "imputed max populations"
  #     } else {
  #       missing.data <- "imputed max global"
  #     }
  #   }
  # }
  
  
  # Modify filename
  filename <- stringi::stri_join(directory.subsample, filename, sep = "")
  filename <- stringi::stri_replace_all_fixed(
    filename,
    pattern = ".txt",
    replacement = stringi::stri_join(
      "", "iteration", i, "markers", x, filename.imp, sep = "_")
  )
  
  if (assignment.analysis == "gsi_sim") {
    assignment <- assignment_gsi_sim(
      input = input,
      strata.df = strata.df,
      select.markers = select.markers,
      markers.names = markers.names,
      missing.data = missing.data,
      i = i,
      m = x,
      holdout = holdout,
      filename = filename,
      directory.subsample = directory.subsample,
      keep.gsi.files = keep.gsi.files,
      sampling.method = sampling.method,
      thl = thl
    )
  }
  
  
  if (assignment.analysis == "adegenet") {
    assignment <- assignment_adegenet(
      data = genind.object,
      select.markers = select.markers,
      markers.names = markers.names,
      missing.data = missing.data,
      adegenet.dapc.opt = adegenet.dapc.opt,
      adegenet.n.rep = adegenet.n.rep,
      adegenet.training = adegenet.training,
      i = i, 
      m = x,
      sampling.method = sampling.method,
      subsample.id = subsample.id,
      holdout = holdout
    )
  }
  # unused objects
  # x <- select.markers <- markers.names <- NULL
  # filename.imp <- missing.data <- filename <- NULL
  
  # x <- as.character(x)
  # assignment.marker <- list()
  # assignment.marker[[x]] <- assignment
  # return(assignment.marker)
  return(assignment)
}#End assignment_marker_loop

# assignment_ranking-----------------------------------------------------------
#' @title assignment_ranking
#' @description assignment_ranking
#' @rdname assignment_ranking
#' @export
#' @keywords internal
assignment_ranking <- function(
  iterations.list = NULL,
  thl = NULL,
  input = NULL,
  holdout.individuals = NULL,
  directory.subsample = NULL,
  marker.number = NULL,
  assignment.analysis = NULL,
  genind.object = NULL,
  strata.df = NULL,
  sampling.method = NULL,
  subsample.id = NULL,
  adegenet.dapc.opt = NULL,
  adegenet.n.rep = NULL,
  adegenet.training = NULL,
  filename = NULL,
  keep.gsi.files = NULL
) {
  # i <- iterations.list[[1]]
  i <- iterations.list
  
  # Ranking Fst with training dataset (keep holdout individuals out)
  message("Ranking markers based on Fst")
  # THL = "all"
  if (thl == "all") {
    holdout <- NULL
    fst.ranked <- assigner::fst_WC84(
      data = input, 
      pop.levels = unique(pop.labels), pop.labels = unique(pop.labels), strata = NULL, 
      holdout.samples = NULL
    )$fst.ranked
    # if (!is.null(imputation.method)) {
    #   fst.ranked.imp <- assigner::fst_WC84(
    #     data = input.imp, 
    #     pop.levels = unique(pop.labels), pop.labels = unique(pop.labels), strata = NULL, 
    #     holdout.samples = NULL
    #   )$fst.ranked
    # }
  }
  
  if (thl == 1) {
    holdout <- tibble::data_frame(INDIVIDUALS = i)
    fst.ranked <- assigner::fst_WC84(
      data = input, 
      pop.levels = unique(pop.labels), pop.labels = unique(pop.labels), strata = NULL, 
      holdout.samples = holdout$INDIVIDUALS
    )$fst.ranked
    # if (!is.null(imputation.method)) {
    #   fst.ranked.imp <- assigner::fst_WC84(
    #     data = input.imp, 
    #     pop.levels = unique(pop.labels), pop.labels = unique(pop.labels), strata = NULL, 
    #     holdout.samples = holdout$INDIVIDUALS
    #   )$fst.ranked
    # }
  }
  
  # thl proportion or > 1
  if (thl != 1 && thl != "all") {
    holdout <- dplyr::filter(.data = holdout.individuals, ITERATIONS %in% i)
    fst.ranked <- assigner::fst_WC84(
      data = input, 
      pop.levels = unique(pop.labels), pop.labels = unique(pop.labels), strata = NULL, 
      holdout.samples = holdout$INDIVIDUALS
    )$fst.ranked
    # if (!is.null(imputation.method)) {
    #   fst.ranked.imp <- assigner::fst_WC84(
    #     data = input.imp, 
    #     pop.levels = unique(pop.labels), pop.labels = unique(pop.labels), strata = NULL, 
    #     holdout.samples = holdout$INDIVIDUALS
    #   )$fst.ranked
    # }
  }
  
  readr::write_tsv(
    x = fst.ranked, 
    path = stringi::stri_join(directory.subsample,
                              "fst.ranked_", i, ".tsv", sep = ""),
    col_names = TRUE, 
    append = FALSE
  )
  
  # if (!is.null(imputation.method)) {  # With imputations
  #   readr::write_tsv(
  #     x = fst.ranked.imp, 
  #     path = stringi::stri_join(directory.subsample,
  #                               "fst.ranked_", i, "_imputed",".tsv", sep = ""),
  #     col_names = TRUE, 
  #     append = FALSE
  #   )
  # }
  
  # looping through the markers
  message("Going throught the marker.number")
  assignment.marker <- list()
  assignment.marker <- suppressWarnings(
    purrr::map(
      .x = marker.number, 
      .f = assignment_marker_loop,
      assignment.analysis = assignment.analysis,
      fst.ranked = fst.ranked,
      i = i,
      input = input,
      genind.object = genind.object,
      strata.df = strata.df,
      pop.labels = unique(pop.labels),
      sampling.method = sampling.method,
      subsample.id = subsample.id,
      holdout = holdout,
      thl = thl,
      adegenet.dapc.opt = adegenet.dapc.opt,
      adegenet.n.rep = adegenet.n.rep,
      adegenet.training = adegenet.training,
      directory.subsample = directory.subsample,
      filename = filename,
      keep.gsi.files = keep.gsi.files
    ) %>%
      dplyr::bind_rows(.)
  )
  
  # if (!is.null(imputation.method)) {
  #   assignment.marker.imp <- list()
  #   assignment.marker.imp <- suppressWarnings(
  #     purrr::map(
  #       .x = marker.number, 
  #       .f = assignment_marker_loop,
  #       assignment.analysis = assignment.analysis,
  #       fst.ranked = fst.ranked.imp,
  #       i = i,
  #       input = input.imp,
  #       genind.object = genind.object.imp,
  #       strata.df = strata.df,
  #       sampling.method = sampling.method,
  #       subsample.id = subsample.id,
  #       holdout = holdout,
  #       thl = thl,
  #       adegenet.dapc.opt = adegenet.dapc.opt,
  #       adegenet.n.rep = adegenet.n.rep,
  #       adegenet.training = adegenet.training,
  #       directory.subsample = directory.subsample,
  #       filename = filename,
  #       keep.gsi.files = keep.gsi.files
  #     ) %>%
  #       dplyr::bind_rows(.)
  #   )
  #   
  #   assignment.res.summary <- suppressWarnings(
  #     dplyr::bind_rows(assignment.marker, assignment.marker.imp)
  #   )
  # } else {
    assignment.res.summary <- assignment.marker
  # }
  assignment.marker.imp <- assignment.marker <- NULL
  
  message("Summarizing the assignment analysis results by iterations and marker group")
  
  readr::write_tsv(
    x = assignment.res.summary,
    path = stringi::stri_join(directory.subsample,
                              "assignment_", i, ".tsv", sep = ""), 
    col_names = TRUE, 
    append = FALSE
  )
  holdout <- fst.ranked <- fst.ranked.imp <- i <- NULL
  return(assignment.res.summary)
}  # End assignment_ranking


# assignment_function-----------------------------------------------------------
#' @title assignment_function
#' @description assignment_function
#' @rdname assignment_function
#' @export
#' @keywords internal

assignment_function <- function(
  x,
  input = NULL,
  subsample = NULL,
  assignment.analysis = "gsi_sim",
  marker.number = NULL,
  sampling.method = "random",
  iteration.method = 10,
  filename = "assignment_data.txt",
  directory = NULL,
  keep.gsi.files = FALSE,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1,
  manage.all = NULL,
  random.seed = NULL,
  base.filename = NULL,
  res.list = NULL,
  thl = 1,
  adegenet.dapc.opt = NULL,
  adegenet.n.rep = NULL,
  adegenet.training = NULL
) {
  
  pop.labels <- NULL
  
  # x <- subsample.list[[1]] # test
  subsample.id <- unique(x$SUBSAMPLE)
  
  if (!is.null(subsample)) {
    message("Analyzing subsample: ", subsample.id)
  }
  
  # Updating directories for subsampling
  if (is.null(subsample)) {
    directory.subsample <- directory
  } else {
    directory.subsample <- paste0(directory, "subsample_", subsample.id, "/")
    dir.create(file.path(directory.subsample))
  }
  
  # Keep only the subsample
  input <- dplyr::semi_join(input, x, by = c("POP_ID", "INDIVIDUALS"))
  
  # unused object
  x <- NULL
  
  # Markers in common between all populations (optional) ---------------------
  # if (common.markers) { # keep only markers present in all pop
  #   input <- radiator::keep_common_markers(data = input)$input
  # } # End common markers
  # 
  # Minor Allele Frequency filter --------------------------------------------
  # if (!is.null(maf.thresholds)) {
  #   input <- radiator::filter_maf(
  #   data = input,
  #   interactive.filter = FALSE,
  #   maf.thresholds = maf.thresholds,
  #   parallel.core = parallel.core,
  #   verbose = FALSE)$tidy.filtered.maf
  # } # End of MAF filters
  
  # Keep a new strata df -------------------------------------------------------
  strata.df <- dplyr::distinct(.data = input, INDIVIDUALS, POP_ID)
  
  # Adegenet no imputations --------------------------------------------------
  if (assignment.analysis == "adegenet") {
    message("Creating genind object")
    genind.object <- radiator::write_genind(data = input)
  } else {
    genind.object <- NULL
  }
  
  # Imputations --------------------------------------------------------------
  # if (!is.null(imputation.method)) {
  #   message("Preparing the data for imputations")
  #   input.imp <- radiator::radiator_imputations_module(
  #     data = input, 
  #     imputation.method = imputation.method, 
  #     hierarchical.levels = hierarchical.levels, 
  #     verbose = verbose, 
  #     parallel.core = parallel.core, 
  #     filename = stringi::stri_join(directory, "dataset.imputed.tsv")
  #   )
  #   res.list$tidy.data.imputed <- input.imp
  #   # prepare the imputed dataset for gsi_sim or adegenet
  #   message("Preparing imputed data set for assignement analysis")
  #   
  #   # adegenet
  #   if (assignment.analysis == "adegenet") {
  #     genind.object.imp <- radiator::write_genind(data = input.imp)
  #   } else {
  #     genind.object.imp <- NULL
  #   } # end adegenet
  #   
  # } else {
  #   input.imp <- NULL
  #   genind.object.imp <- NULL
  # } # End imputations
  
  # Sampling of markers ------------------------------------------------------
  # unique list of markers after all the filtering
  unique.markers <- dplyr::distinct(.data = input, MARKERS)
  
  # if "all" is present in the list, change to the maximum number of markers
  marker.number <- as.numeric(
    stringi::stri_replace_all_fixed(str = marker.number, 
                                    pattern = "all", 
                                    replacement = nrow(unique.markers), 
                                    vectorize_all = TRUE)
  )
  
  # In marker.number, remove marker group higher than the max number of markers
  removing.marker <- purrr::keep(.x = marker.number, .p = marker.number > nrow(unique.markers))
  
  if (length(removing.marker) > 0) {
    message(
      "Removing marker.number higher than the max number of markers: ", 
      stringi::stri_join(removing.marker, collapse = ", ")
    )
    marker.number <- purrr::discard(.x = marker.number, .p = marker.number > nrow(unique.markers))
  }
  
  # Random method ------------------------------------------------------------
  if (sampling.method == "random") {
    message("Conducting Assignment analysis with markers selected randomly")
    # Number of times to repeat the sampling of markers
    iterations.list <- seq(iteration.method)
    # iterations.list <- 1:10 # test
    
    # Go through the function with the marker number selected
    message("Making a list containing all the markers combinations")
    
    if (manage.all) {# only 1 iterations when random method + max markers
      marker.number.mod <- purrr::discard(
        .x = marker.number, .p = marker.number == nrow(unique.markers)
      )
      marker.random.list <- list() # marker random list
      for (m in marker.number.mod) {
        res <- purrr::map(
          .x = 1:iteration.method,
          .f = marker_selection,
          m = m,
          random.seed = random.seed,
          unique.markers = unique.markers
        )
        marker.random.list[[m]] <- res
      }
      m <- nrow(unique.markers)
      marker.random.list[[m]] <- purrr::map(
        .x = 1:1,
        .f = marker_selection,
        m = m,
        random.seed = random.seed,
        unique.markers = unique.markers
      )
      
      marker.random.list <- purrr::flatten(marker.random.list)
      
    } else {# all the iterations requested by user
      marker.random.list <- list() # marker random list
      for (m in marker.number) {
        res <- purrr::map(
          .x = 1:iteration.method,
          .f = marker_selection,
          m = m,
          random.seed = random.seed,
          unique.markers = unique.markers
        )
        
        marker.random.list[[m]] <- res
      }
      marker.random.list <- purrr::flatten(marker.random.list)
    }
    # Write the results to the working directory
    marker.random.df <- tibble::as_data_frame(dplyr::bind_rows(marker.random.list))
    readr::write_tsv(
      x = marker.random.df,
      path = paste0(directory.subsample, "markers.random.tsv"),
      col_names = TRUE, append = FALSE
    )
    
    res.list$marker.random.df <- marker.random.df
    
    message("Starting parallel computations for the assignment analysis: random markers")
    message("For progress: monitor activity in the folder...")
    
    assignment.res <- NULL
    assignment.res <- suppressWarnings(
      .assigner_parallel(
        X = marker.random.list,
        FUN = assignment_random,
        assignment.analysis = assignment.analysis,
        input = input,
        genind.object = genind.object,
        strata.df = strata.df,
        directory.subsample = directory.subsample,
        keep.gsi.files = keep.gsi.files,
        sampling.method = sampling.method,
        subsample.id = subsample.id,
        filename = filename,
        adegenet.n.rep = adegenet.n.rep,
        adegenet.training = adegenet.training,
        holdout = NULL,
        mc.preschedule = FALSE,
        mc.silent = FALSE,
        mc.cleanup = TRUE,
        mc.cores = parallel.core
      ) %>% 
        dplyr::bind_rows(assignment.res)
    )
    
    # with imputations
    # if (!is.null(imputation.method)) {
    #   assignment.res.imp <- NULL
    #   assignment.res.imp <- suppressWarnings(
    #     .assigner_parallel(
    #       X = marker.random.list, 
    #       FUN = assignment_random,
    #       assignment.analysis = assignment.analysis,
    #       input = input.imp,
    #       genind.object = genind.object.imp,
    #       strata.df = strata.df,
    #       directory.subsample = directory.subsample,
    #       keep.gsi.files = keep.gsi.files,
    #       sampling.method = sampling.method,
    #       subsample.id = subsample.id,
    #       filename = filename,
    #       adegenet.n.rep = adegenet.n.rep,
    #       adegenet.training = adegenet.training,
    #       holdout = NULL,
    #       mc.preschedule = FALSE, 
    #       mc.silent = FALSE,
    #       mc.cleanup = TRUE,
    #       mc.cores = parallel.core
    #     )
    #   )
    #   assignment.res <- suppressWarnings(dplyr::bind_rows(assignment.res.imp) %>%
    #                                        dplyr::bind_rows(assignment.res))
    # }
    # Compiling the results
    message("Compiling results")
    if (assignment.analysis == "adegenet") {
      assignment.res <- suppressWarnings(
        dplyr::rename(.data = assignment.res,CURRENT = POP_ID) %>% 
          dplyr::mutate(
            SUBSAMPLE = rep(subsample.id, n()),
            CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
            CURRENT = droplevels(CURRENT)
          ) %>% 
          dplyr::arrange(CURRENT, MARKER_NUMBER, MISSING_DATA, ITERATIONS)
      )
    } else {
      assignment.res <- suppressWarnings(
        dplyr::mutate(.data = assignment.res, SUBSAMPLE = rep(subsample.id, n())) %>% 
          dplyr::arrange(CURRENT, INDIVIDUALS, MARKER_NUMBER, MISSING_DATA, ITERATIONS)
      )
    }
    
    # Write the tables to directory
    # assignment results
    if (assignment.analysis == "gsi_sim") {
      if (is.null(subsample)) {
        # if (is.null(imputation.method)) {
          filename.assignment.res <- stringi::stri_join(
            "assignment", sampling.method, "no.imputation", "results", 
            "individuals", "iterations", "tsv", sep = "."
          )
        # } else {# with imputations
        #   filename.assignment.res <- stringi::stri_join(
        #     "assignment", sampling.method, "imputed", "results", "individuals", 
        #     "iterations", "tsv", sep = "."
        #   )
        # }
      } else {# with subsampling
        # if (is.null(imputation.method)) {
          filename.assignment.res <- stringi::stri_join(
            "assignment", sampling.method, "no.imputation", "results", "individuals",
            "iterations", "subsample", subsample.id, "tsv", sep = "."
          )
        # } else {# with imputations
        #   filename.assignment.res <- stringi::stri_join(
        #     "assignment", sampling.method, "imputed", "results", "individuals", 
        #     "iterations", "subsample", subsample.id, "tsv", sep = "."
        #   )
        # }
      }
      readr::write_tsv(
        x = assignment.res, 
        path = paste0(directory.subsample,filename.assignment.res), 
        col_names = TRUE, 
        append = FALSE
      )
    } else {# with adegenet
      if (is.null(subsample)) {
        # if (is.null(imputation.method)) {
          filename.assignment.res <- stringi::stri_join(
            "assignment", sampling.method, "no.imputation", "results", 
            "iterations", "tsv", sep = "."
          )
        # } else {# with imputations
          # filename.assignment.res <- stringi::stri_join(
          #   "assignment", sampling.method, "imputed", "results", "iterations", 
          #   "tsv", sep = "."
          # )
        # }
      } else {# with subsampling
        # if (is.null(imputation.method)) {
          filename.assignment.res <- stringi::stri_join(
            "assignment", sampling.method, "no.imputation", "results", 
            "iterations", "subsample", subsample.id, "tsv", sep = "."
          )
        # } else {# with imputations
        #   filename.assignment.res <- stringi::stri_join(
        #     "assignment", sampling.method, "imputed", "results", "iterations",
        #     "subsample", subsample.id, "tsv", sep = "."
        #   )
        # }
      }
      readr::write_tsv(
        x = assignment.res, 
        path = paste0(directory.subsample,filename.assignment.res), 
        col_names = TRUE, append = FALSE
      )
    }
    
    if (assignment.analysis == "gsi_sim") {
      assignment.stats.pop <- suppressWarnings(
        assignment.res %>%
          dplyr::group_by(CURRENT, INFERRED, MARKER_NUMBER, MISSING_DATA, ITERATIONS, METHOD) %>%
          dplyr::tally(.) %>%
          dplyr::group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, ITERATIONS, METHOD) %>%
          dplyr::mutate(TOTAL = sum(n)) %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate(MEAN_i = round(n/TOTAL*100, 0)) %>%
          dplyr::filter(as.character(CURRENT) == as.character(INFERRED)) %>%
          dplyr::select(CURRENT, MEAN_i, MARKER_NUMBER, MISSING_DATA, ITERATIONS, METHOD) %>%
          dplyr::mutate(
            CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = T),
            CURRENT = droplevels(CURRENT)
          ) %>%
          dplyr::group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
          dplyr::summarise(
            MEAN = round(mean(MEAN_i), 2),
            SE = round(sqrt(stats::var(MEAN_i)/length(MEAN_i)), 2),
            MIN = round(min(MEAN_i), 2),
            MAX = round(max(MEAN_i), 2),
            MEDIAN = round(stats::median(MEAN_i), 2),
            QUANTILE25 = round(stats::quantile(MEAN_i, 0.25), 2),
            QUANTILE75 = round(stats::quantile(MEAN_i, 0.75), 2)
          ) %>%
          dplyr::ungroup(.) %>% 
          dplyr::arrange(CURRENT, MARKER_NUMBER)
      )
    }
    
    if (assignment.analysis == "adegenet") {
      assignment.stats.pop <- suppressWarnings(
        assignment.res %>%
          dplyr::group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
          dplyr::summarise(
            MEAN = round(mean(ASSIGNMENT_PERC), 2),
            SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
            MIN = round(min(ASSIGNMENT_PERC), 2),
            MAX = round(max(ASSIGNMENT_PERC), 2),
            MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
            QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
            QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
          ) %>%
          dplyr::ungroup(.) %>% 
          dplyr::mutate(
            CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
            CURRENT = droplevels(CURRENT)
          ) %>%
          dplyr::arrange(CURRENT, MARKER_NUMBER)
      )
    }
    
    # Next step is common for gsi_sim and adegenet
    pop.levels.assignment.stats.overall <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
    
    assignment.stats.overall <- assignment.stats.pop %>%
      dplyr::group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
      dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
      dplyr::summarise(
        MEAN = round(mean(ASSIGNMENT_PERC), 2),
        SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
        MIN = round(min(ASSIGNMENT_PERC), 2),
        MAX = round(max(ASSIGNMENT_PERC), 2),
        MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
        QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
        QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
      ) %>%
      dplyr::mutate(CURRENT = rep("OVERALL", n())) %>%
      dplyr::ungroup(.) %>% 
      dplyr::arrange(CURRENT, MARKER_NUMBER)
    
    assignment.summary.stats <- suppressWarnings(
      dplyr::bind_rows(assignment.stats.pop, assignment.stats.overall) %>%
        dplyr::mutate(CURRENT = factor(CURRENT, levels = pop.levels.assignment.stats.overall, ordered = TRUE)) %>%
        dplyr::arrange(CURRENT, MARKER_NUMBER) %>%
        dplyr::mutate(
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE,
          ITERATIONS = rep(iteration.method, n())
        ) %>%
        dplyr::select(CURRENT, MARKER_NUMBER, MEAN, MEDIAN, SE, MIN, MAX, QUANTILE25, QUANTILE75, SE_MIN, SE_MAX, METHOD, MISSING_DATA, ITERATIONS)
    )
    
    
    # Write the tables to directory
    # assignment summary stats
    if (is.null(subsample)) {
      # if (is.null(imputation.method)) {
        filename.assignment.sum <- stringi::stri_join(
          "assignment", sampling.method, "no.imputation", "results", 
          "summary.stats", "tsv", sep = "."
        )
      # } else {# with imputations
      #   filename.assignment.sum <- stringi::stri_join(
      #     "assignment", sampling.method, "imputed", "results", "summary.stats",
      #     "tsv", sep = "."
      #   )
      # }
    } else {# with subsampling
      # if (is.null(imputation.method)) {
        filename.assignment.sum <- stringi::stri_join(
          "assignment", sampling.method, "no.imputation", "results", 
          "summary.stats", "subsample", subsample.id, "tsv", sep = "."
        )
      # } else {# with imputations
      #   filename.assignment.sum <- stringi::stri_join(
      #     "assignment", sampling.method, "imputed", "results", "summary.stats",
      #     "subsample", subsample.id, "tsv", sep = "."
      #   )
      # }
    }
    readr::write_tsv(
      x = assignment.summary.stats, 
      path = paste0(directory.subsample,filename.assignment.sum), 
      col_names = TRUE, append = FALSE
    )
    
  } # End method random
  
  # Ranked method ------------------------------------------------------------
  if (sampling.method == "ranked") {
    message("Conducting Assignment analysis with ranked markers")
    
    # individuals and pop df
    ind.pop.df <- dplyr::ungroup(input) %>% dplyr::distinct(POP_ID, INDIVIDUALS)
    
    # thl selection
    message("Using thl method, ranking Fst with training samples...")
    
    # Will go through the individuals in the list one by one.
    if (thl == 1) {
      iterations.list <- ind.pop.df$INDIVIDUALS
      holdout.individuals <- ind.pop.df %>%
        dplyr::mutate(ITERATIONS = stringi::stri_join("HOLDOUT", seq(1:nrow(ind.pop.df)), sep = "_"))
    }
    
    # no holdout for that one
    if (thl == "all") {
      iterations.list <- iteration.method
      holdout.individuals <- NULL
      message("Warning: using all the individuals for ranking markers based on Fst\nNo holdout samples")
      message("Recommended reading: \nAnderson, E. C. (2010) Assessing the power of informative subsets of
loci for population assignment: standard methods are upwardly biased.\nMolecular ecology resources 10, 4:701-710.")
    }
    
    if (thl != 1 && thl != "all") {
      # Set seed for random sampling
      if (is.null(random.seed)) {
        random.seed <- sample(x = 1:1000000, size = 1)
        set.seed(random.seed)
      } else {
        set.seed(random.seed)
      }
      
      
      # Create x (iterations) list of y (thl) proportion of individuals per pop.
      if (stringi::stri_detect_fixed(thl, ".") && thl < 1) {
        holdout.individuals.list <- list()
        iterations.list <- 1:iteration.method
        for (x in iterations.list) {
          holdout.individuals <- ind.pop.df %>%
            dplyr::group_by(POP_ID) %>%
            dplyr::sample_frac(thl, replace = FALSE) %>%  # sampling fraction for each pop
            dplyr::arrange(POP_ID, INDIVIDUALS) %>%
            dplyr::ungroup(.) %>%
            dplyr::select(INDIVIDUALS) %>%
            dplyr::mutate(
              ITERATIONS = rep(x, n()),
              RANDOM_SEED_NUMBER = rep(random.seed, n())
            )
          
          holdout.individuals.list[[x]] <- holdout.individuals
        }
        holdout.individuals <- dplyr::bind_rows(holdout.individuals.list)
        holdout.individuals.list <- NULL
      }
      
      # Create x (iterations) list of y (thl) individuals per pop.
      if (thl > 1 && thl != "all") {
        holdout.individuals.list <- list()
        iterations.list <- 1:iteration.method
        for (x in iterations.list) {
          holdout.individuals <- ind.pop.df %>%
            dplyr::group_by(POP_ID) %>%
            dplyr::sample_n(tbl = ., size = thl, replace = FALSE) %>% # sampling individuals for each pop
            dplyr::arrange(POP_ID, INDIVIDUALS) %>%
            dplyr::ungroup(.) %>%
            dplyr::select(INDIVIDUALS) %>%
            dplyr::mutate(
              ITERATIONS = rep(x, n()),
              RANDOM_SEED_NUMBER = rep(random.seed, n())
            )
          holdout.individuals.list[[x]] <- holdout.individuals
        }
        holdout.individuals <- dplyr::bind_rows(holdout.individuals.list)
        holdout.individuals.list <- NULL
      }
      x <- NULL
    } # End tracking holdout individuals
    
    ind.pop.df <- NULL
    
    readr::write_tsv(
      x = holdout.individuals, 
      path = paste0(directory.subsample,"holdout.individuals.tsv"), 
      col_names = TRUE, 
      append = FALSE
    )
    res.list$holdout.individuals <- holdout.individuals
    message("Holdout samples saved in your folder")
    
    # Going through the loop of holdout individuals
    message("Starting parallel computations for the assignment analysis: ranked markers")
    message("For progress: monitor activity in the folder...")
    assignment.res <- list()
    assignment.res <- .assigner_parallel(
      # assignment.res <- parallel::mclapply(
      X = iterations.list, 
      FUN = assignment_ranking, 
      thl = thl,
      input = input,
      holdout.individuals = holdout.individuals,
      directory.subsample = directory.subsample,
      marker.number = marker.number,
      assignment.analysis = assignment.analysis,
      genind.object = genind.object,
      strata.df = strata.df,
      sampling.method = sampling.method,
      subsample.id = subsample.id,
      adegenet.dapc.opt = adegenet.dapc.opt,
      adegenet.n.rep = adegenet.n.rep,
      adegenet.training = adegenet.training,
      filename = filename,
      keep.gsi.files = keep.gsi.files,
      mc.preschedule = FALSE, 
      mc.silent = FALSE,
      mc.cleanup = FALSE,
      mc.cores = parallel.core
    )
    
    # Compiling the results
    message("Compiling results")
    assignment.res.summary <- suppressWarnings(
      dplyr::bind_rows(assignment.res) %>% 
        dplyr::mutate(SUBSAMPLE = rep(subsample.id, n())) %>% 
        # arrange(CURRENT, INDIVIDUALS, MARKER_NUMBER, MISSING_DATA, ITERATIONS)
        dplyr::arrange(CURRENT, INDIVIDUALS, MARKER_NUMBER, MISSING_DATA)
    )
    
    # assignment results
    if (is.null(subsample)) {
      # if (is.null(imputation.method)) {
        filename.assignment.res <- stringi::stri_join("assignment", sampling.method, "no.imputation", "results", "individuals", "iterations", "tsv", sep = ".")
      # } else {# with imputations
      #   filename.assignment.res <- stringi::stri_join("assignment", sampling.method, "imputed", "results", "individuals", "iterations", "tsv", sep = ".")
      # }
    } else {# with subsampling
      # if (is.null(imputation.method)) {
        filename.assignment.res <- stringi::stri_join("assignment", sampling.method, "no.imputation", "results", "individuals","iterations", "subsample", subsample.id, "tsv", sep = ".")
      # } else {# with imputations
      #   filename.assignment.res <- stringi::stri_join("assignment", sampling.method, "imputed", "results", "individuals", "iterations", "subsample", subsample.id, "tsv", sep = ".")
      # }
    }
    readr::write_tsv(
      x = assignment.res.summary,
      path = paste0(directory.subsample,filename.assignment.res),
      col_names = TRUE, append = FALSE
    )
    
    
    if (thl == 1 | thl == "all") {
      assignment.stats.pop <- assignment.res.summary %>%
        dplyr::mutate(
          CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
          CURRENT = droplevels(CURRENT)
        ) %>% 
        dplyr::group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
        dplyr::summarise(
          n = length(CURRENT[as.character(CURRENT) == as.character(INFERRED)]),
          TOTAL = length(CURRENT)
        ) %>%
        dplyr::ungroup(.) %>% 
        dplyr::mutate(MEAN = round(n/TOTAL*100, 0)) %>% 
        dplyr::select(-n, -TOTAL)
      
      pop.levels.assignment.stats.overall <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
      
      assignment.stats.overall <- assignment.stats.pop %>% 
        dplyr::group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
        dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
        dplyr::summarise(
          MEAN = round(mean(ASSIGNMENT_PERC), 2),
          SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
          MIN = round(min(ASSIGNMENT_PERC), 2),
          MAX = round(max(ASSIGNMENT_PERC), 2),
          MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
          QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
          QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
        ) %>% 
        dplyr::mutate(CURRENT = rep("OVERALL", n())) %>% 
        dplyr::arrange(CURRENT, MARKER_NUMBER)
      
      assignment.summary.stats <- suppressWarnings(
        dplyr::bind_rows(assignment.stats.pop, assignment.stats.overall) %>%
          dplyr::mutate(CURRENT = factor(CURRENT, levels = pop.levels.assignment.stats.overall, ordered = TRUE)) %>%
          dplyr::arrange(CURRENT, MARKER_NUMBER) %>%
          dplyr::mutate(
            SE_MIN = MEAN - SE,
            SE_MAX = MEAN + SE
          )
      )
      
    } else {
      # if (assignment.analysis == "gsi_sim") {
      assignment.res.summary.prep <- assignment.res.summary %>% 
        dplyr::group_by(CURRENT, MARKER_NUMBER, METHOD, MISSING_DATA, ITERATIONS) %>%
        dplyr::summarise(
          n = length(CURRENT[as.character(CURRENT) == as.character(INFERRED)]),
          TOTAL = length(CURRENT)
        ) %>%
        dplyr::ungroup(.) %>% 
        dplyr::mutate(ASSIGNMENT_PERC = round(n/TOTAL*100, 0)) %>% 
        dplyr::select(-n, -TOTAL)
      # }
      
      if (is.null(subsample)) {
        # if (is.null(imputation.method)) {
          filename.assignment.res.sum <- stringi::stri_join("assignment", sampling.method, "no.imputation", "results", "summary", "tsv", sep = ".")
        # } else {# with imputations
        #   filename.assignment.res.sum <- stringi::stri_join("assignment", sampling.method, "imputed", "results", "summary", "tsv", sep = ".")
        # }
      } else {# with subsampling
        # if (is.null(imputation.method)) {
          filename.assignment.res.sum <- stringi::stri_join("assignment", sampling.method, "no.imputation", "results", "summary", "subsample", subsample.id, "tsv", sep = ".")
        # } else {# with imputations
        #   filename.assignment.res.sum <- stringi::stri_join("assignment", sampling.method, "imputed", "results", "summary", "subsample", subsample.id, "tsv", sep = ".")
        # }
      }
      readr::write_tsv(
        x = assignment.res.summary.prep,
        path = paste0(directory.subsample,filename.assignment.res.sum),
        col_names = TRUE, append = FALSE
      )
      
      assignment.stats.pop <- assignment.res.summary.prep %>%
        # assignment.stats.pop <- assignment.res.summary %>%
        dplyr::mutate(
          CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
          CURRENT = droplevels(CURRENT)
        ) %>%
        dplyr::group_by(CURRENT, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
        dplyr::summarise(
          MEAN = round(mean(ASSIGNMENT_PERC), 2),
          SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
          MIN = round(min(ASSIGNMENT_PERC), 2),
          MAX = round(max(ASSIGNMENT_PERC), 2),
          MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
          QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
          QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2),
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE
        ) %>%
        dplyr::ungroup(.) %>% 
        dplyr::arrange(CURRENT, MARKER_NUMBER)
      
      pop.levels.assignment.stats.overall <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
      
      assignment.stats.overall <- assignment.stats.pop %>%
        dplyr::group_by(MARKER_NUMBER, METHOD, MISSING_DATA) %>%
        dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
        dplyr::summarise(
          MEAN = round(mean(ASSIGNMENT_PERC), 2),
          SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
          MIN = round(min(ASSIGNMENT_PERC), 2),
          MAX = round(max(ASSIGNMENT_PERC), 2),
          MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
          QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
          QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2),
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE
        ) %>%
        dplyr::mutate(CURRENT = rep("OVERALL", n())) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(CURRENT, MARKER_NUMBER)
      
      assignment.summary.stats <- suppressWarnings(
        dplyr::bind_rows(assignment.stats.pop, assignment.stats.overall) %>%
          dplyr::mutate(CURRENT = factor(CURRENT, levels = pop.levels.assignment.stats.overall, ordered = TRUE)) %>%
          dplyr::arrange(CURRENT, MARKER_NUMBER)
      )
    } # End thl != 1
    
    # Write the tables to directory
    # assignment summary stats
    if (is.null(subsample)) {
      # if (is.null(imputation.method)) {
        filename.assignment.sum <- stringi::stri_join("assignment", sampling.method, "no.imputation", "results", "summary.stats", "tsv", sep = ".")
      # } else {# with imputations
      #   filename.assignment.sum <- stringi::stri_join("assignment", sampling.method, "imputed", "results", "summary.stats", "tsv", sep = ".")
      # }
    } else {# with subsampling
      # if (is.null(imputation.method)) {
        filename.assignment.sum <- stringi::stri_join("assignment", sampling.method, "no.imputation", "results", "summary.stats", "subsample", subsample.id, "tsv", sep = ".")
      # } else {# with imputations
      #   filename.assignment.sum <- stringi::stri_join("assignment", sampling.method, "imputed", "results", "summary.stats", "subsample", subsample.id, "tsv", sep = ".")
      # }
    }
    readr::write_tsv(
      x = assignment.summary.stats,
      path = paste0(directory.subsample,filename.assignment.sum),
      col_names = TRUE, append = FALSE
    )
  } # End of ranked thl method
  
  # update the assignment with subsampling iterations id
  assignment.summary.stats <- assignment.summary.stats %>% 
    dplyr::mutate(SUBSAMPLE = rep(subsample.id, n()))
  return(assignment.summary.stats)
}# End assignment_function
