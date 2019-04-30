# Assignment analysis of massive parallel sequencing data
#' @name assignment_ngs

#' @title Assignment analysis tailored for RADseq data

#' @description
#' The arguments in the \code{assignment_ngs} function were tailored for the
#' reality of RADseq data for assignment analysis while
#' maintaining a reproducible workflow.
#' Assignment assumptions are listed in the section below.
#' 
#' \itemize{
#'   \item \strong{Input file:} various file format are supported (see \code{data} argument below).
#'   \item \strong{Cross-Validations:} Markers can be randomly selected for a classic LOO (Leave-One-Out)
#'   assignment or chosen based on ranked Fst for a thl
#'   (Training, Holdout, Leave-one-out) assignment analysis.
#'   \item \strong{Assignment analysis:} conducted in 
#'   \href{https://github.com/eriqande/gsi_sim}{gsi_sim}, a tool 
#'   for doing and simulating genetic stock identification and 
#'   developed by Eric C. Anderson, or 
#'   \href{https://github.com/thibautjombart/adegenet}{adegenet}, 
#'  an R package developed by Thibaul Jombart.
#'   \item \strong{Parallel:} The assignment can be conduncted on multiple CPUs.
#'   The R GUI is unstable with this functions, I recommend using 
#'   \href{https://www.rstudio.com/products/rstudio/download/}{RStudio}. 
#' }
#' @param data Several input format are accepted. assigner uses \pkg{radiator}
#' \code{\link[radiator]{tidy_genomic_data}} module to import the data.
#' See function documentation for more details.

#' @inheritParams radiator::tidy_genomic_data 
#' @inheritParams radiator::read_strata 

#' @param assignment.analysis (character) Assignment analysis conducted with 
#' \code{assignment.analysis = "gsi_sim"} or 
#' \code{assignment.analysis = "adegenet"}.
#' See \strong{Details} section below for installing
#' \href{https://github.com/eriqande/gsi_sim}{gsi_sim}.
#' 

#' @param marker.number (Integer or string of number or "all") The assignment
#' analysis can use all your markers (default) or a subset of your markers.
#' e.g. To test 500, 1000, 2000 and all the markers:
#' \code{marker.number = c(500, 1000, 2000, "all")}.
#' To use only 500 makers \code{marker.number = 500}. How those markers are sampled
#' is determined with the argument \code{markers.sampling}, next.
#' Default = \code{marker.number = "all"}.

#' @param markers.sampling (character) 2 options for markers selection:
#' \enumerate{
#' \item \code{markers.sampling == "random"} the subset of markers are selected 
#' randomly, this is the classic \emph{Leave-One-Out} (LOO) assignment.
#' \item \code{markers.sampling == "ranked"} the subset of markers are first 
#' ranked based on an overall \emph{decreasing} Fst values.
#' The Fst is computed with \code{\link{fst_WC84}} function, that uses a fast 
#' implementation of Weir and Cockerham 1984 Fst/Theta equations. This selection
#' method is used during \emph{Training-Holdout-Leave One Out} (thl)
#' assignment. How many markers are selected is controlled with argument \code{thl}.
#' }

#' @param thl (character, integer, proportion) For \code{markers.sampling = "ranked"} only.
#' Several options are available:
#' \enumerate{
#' \item \code{thl = 1}: Only 1 individual sample is used as holdout. This individual is not
#' participating in the markers ranking. For each marker number,
#' the analysis will be repeated with all the indiviuals in the data set
#' (e.g. 500 individuals, 500 times 500, 1000, 2000 markers). This is the default. 
#' \item \code{proportion}: e.g. \code{thl = 0.15}, 15 percent of the individuals,
#' in each strata, are chosen randomly as holdout individuals.
#' \item \code{thl = "all"}: all individuals are used for ranking (not good) and
#' the argument \code{iteration.method = 1} is set by default. This is for testing
#' only.
#' }
#' Different lists of holdout individuals are generated when the argument
#' \code{iteration.method} (bootstrap) is used.

#' @param iteration.method (integer) 
#' With \strong{random markers selection} method, the iterations argument =
#' the number of iterations to repeat marker resampling. 
#' Default: \code{iteration.method = 10}.
#'
#' With \code{marker.number = c(500, 1000)} and default iterations setting,
#' 500 markers will be randomly chosen 10 times and 1000 markers will be randomly
#' chosen 10 times.
#' 
#' \strong{Notes:} If all the markers are used, with \code{marker.number = "all"}
#' or in a series of marker number groupings \code{marker.number = c(200, 500, "all")}, 
#' the number of iteration is automatically set to 1. The remaining groupings
#' are unaffected.
#' 
#' With \strong{ranked makers selection} method, using \code{thl = 1}, the analysis
#' will be repeated for each individuals in the data set for every
#' \code{marker.number} selected. With a proportion argument \code{thl = 0.15},
#' 15 percent of individuals in each populations are chosen randomly as holdout
#' individuals and this process is reapeated the number of times chosen by the
#' \code{iteration.method} value.

#' @param subsample (Integer or Character, optional) 
#' This argument subsample individuals.
#' With \code{subsample = 36}, 36 individuals in each populations are chosen
#' randomly to represent the dataset. This integer as to be smaller than the
#' population with min sample size, if higher, the minimum sample size found 
#' in the data will be used as default. In doubt, use \code{subsample = "min"},
#' the function will use the smallest population sample size found in the data.
#' The number of times this process is repeated is controlled by the argument
#' \code{iteration.subsample}.
#' Default: \code{subsample = NULL} (no subsampling).

#' @param iteration.subsample (Integer) The number of iterations to repeat 
#' subsampling of individuals.
#' With \code{subsample = 20} and \code{iteration.subsample = 10},
#' 20 individuals/populations will be randomly chosen 10 times.
#' Default: \code{iteration.subsample = 1}.

#' @param ... (optional) To pass further argument for fine-tuning the 
#' function (see advanced section below).

# ... -------------------------- -----------------------------------------------
#' @section Advance mode:
#'
#' Ideally, forget about this section.
#' For advance users, through \emph{dots-dots-dots ...} you can pass several
#' arguments for fine-tuning the function:
#' \enumerate{
#' \item \code{adegenet.dapc.opt} (optional, character) Argument available only when 
#' using:
#' \code{assignment.analysis = "adegenet"} with
#' \code{markers.sampling == "random"}.
#' 
#' The assignment using dapc can use the optimized alpha score 
#' \code{adegenet.dapc.opt == "optim.a.score"} or 
#' cross-validation \code{adegenet.dapc.opt == "xval"}
#' for stability of group membership probabilities. 
#' For fine tuning the trade-off between power of discrimination and over-fitting.
#' See \pkg{adegenet} documentation for more details.
#' \code{adegenet.dapc.opt == "xval"} doesn't work with missing data, so it's 
#' only available with full dataset or \strong{imputed dataset}.
#' With non imputed data or the default: \code{adegenet.dapc.opt == "optim.a.score"}.
#' 
#' \item \code{adegenet.n.rep}: (optional, integer) 
#' When \code{adegenet.dapc.opt == "xval"}, the number of replicates to be 
#' carried out at each level of PC retention. 
#' Default: \code{adegenet.n.rep = 30}.
#' See \pkg{adegenet} documentation for more details.
#' 
#' 
#' \item \code{adegenet.training}: (optional, numeric) 
#' When \code{adegenet.dapc.opt == "xval"}, the proportion of data (individuals) 
#' to be used for the training set.
#' Default: \code{adegenet.training = 0.9}, if all groups have >= 10 members. 
#' Otherwise, training.set scales automatically to the largest proportion 
#' that still ensures all groups will be present in both training 
#' and validation sets.
#' See \pkg{adegenet} documentation for more details.
#' 
#' \item \code{folder}: (optional) The name of the folder created in the working 
#' directory to save the files/results. Default: \code{folder = NULL} will create
#' the folder for you with data and time.
#' 
#' \item \code{filename}: (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. 
#' Several info will be appended to the name of the file.
#' Default \code{assignment_data.txt}.
#' 
#' \item \code{keep.gsi.files}: (logical, optional) With the default, 
#' the intermediate input and output gsi_sim files will be deleted from the 
#' directory when finished processing. I you decide to keep the files, with 
#' \code{keep.gsi.files = TRUE}, remember to allocate a large chunk of the disk 
#' space for the analysis.
#' Default: \code{keep.gsi.files = FALSE} 
#' 
#' \item \code{random.seed}: (integer, optional) For reproducibility, set an integer
#' that will be used inside function that requires randomness. With default,
#' a random number is generated and printed in the appropriate output.
#' Default: \code{random.seed = NULL}.
#' }

#' @section Assumptions:
#' \enumerate{
#' \item \strong{Individuals QC}: Bad individual QC \strong{will bias}
#' the assignment results.
#' \itemize{
#' \item \strong{remove duplicates samples: } when found within the same strata,
#' duplicates generate a false population signal, when they are found 
#' between strata (yes, I've seen it),
#' it's generating noise and the core population signal is diluted.
#' \item \strong{remove individual with outlier heterozygosity: } 
#' unchecked, outlier individuals based on heterozygosity will generate false population
#' signal when the sample as lower heterozygosity and match against several
#' strata (week assignment) when the sample as higher heterozygosity.
#' \item \strong{remove individuals with too many missing: } these individuals
#' will exacerbate or dilute the signal, depending on correlation with heterozygosity
#' and presence of pattern of missingness.
#' }
#' \item \strong{Markers QC}: Bad markers QC \strong{will bias}
#' the assignment results.
#' \itemize{
#' \item \strong{low MAC}: improper Minor Allele Count filtering generate noise. 
#' The LOO and THL methods, both removes samples during model construction,
#' if MAC is too low, the population core signature is greatly impacted at each iteration.
#' \item \strong{Linkage disequilibrium}: remove highly linked markers.
#' \item \strong{HWE}: remove markers in very strong Hardy-Weinberg disequilibrium 
#' likely artefactual and/or genotyping errors.
#' }
#' \item \strong{Strata}: bad K selection will result in poor assingment results.
#' \item \strong{filtered data:} Don't expect to filter your data here.
#' \pkg{radiator} was designed for this, and \code{filter_rad} will handle 
#' the issues mentioned above. By default, the function will only remove
#' monomorphic markers and keep markers in common between strata.
#' }

# Life cycle -------- ----------------------------------------------------------
#' @section Life cycle:
#'
#' Map-independent imputation of missing genotype is avaible in my other R
#' package called \href{https://github.com/thierrygosselin/grur}{grur}.
#'
#' Use \href{https://github.com/thierrygosselin/grur}{grur} to :
#' \enumerate{
#' \item \strong{Visualize your missing data: } before imputing your genotypes,
#' visualize your missing data.
#' Several visual tools are available inside \href{https://github.com/thierrygosselin/grur}{grur} to
#' help you decide the best strategy after.
#' \item \strong{Optimize: }
#' use \href{https://github.com/thierrygosselin/grur}{grur} imputation module
#' and other functions to optimize the imputations of your dataset.
#' You need to test arguments. Failing to conduct tests and adjust imputations arguments
#' will \strong{generate artifacts} and/or \strong{exacerbate bias}.
#' Using defaults is not optional here...
#' \item \strong{genomic_converter: }
#' use the output argument inside \href{https://github.com/thierrygosselin/grur}{grur}
#' imputation module to generate the required formats for assigner (e.g. a tidy dataset)
#' }
#' 
#' 
#' \strong{Deprecated arguments: }
#' 
#' \itemize{
#' \item \code{sampling.method}: renamed \code{markers.sampling}.
#' }
#' 
#' @details Using \href{https://github.com/eriqande/gsi_sim}{gsi_sim}:
#'
#' \code{assignment_ngs} assumes that the command line version of 
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

# Return ------------------------ ----------------------------------------------
#' @return Depending on arguments selected, several folders and files are written:
#' \enumerate{
#' \item \code{assigner_assignment_ngs_args_date@time.tsv}: For reproducibility,
#' the function call, arguments and values used along the default arguments.
#' \item \code{assignment.plot.pdf}: The assignment figure.
#' \item \code{assignment.results.summary.stats.tsv}: tibble
#' of summarized assignment statistics.
#' \item \code{assignment.ranked.results.summary.stats.all.subsamples.tsv}:
#' When subsampling is used, this file contains the raw results of all subsample before 
#' summarizing.

#' }
#' 
#' \strong{THL: Training, Holdout, Leave-one-out}:
#'
#' Intermediate files are written, you can inspect specific iterations
#' and/or subsample:
#' \enumerate{
#' \item \code{assignment.ranked.results.iterations.raw.tsv}: tibble
#' with raw intermediate results of all iterations.
#' \item \code{assignment.ranked.results.iterations.summary.tsv}: tibble 
#' with intermediate summary over iterations.
#' \item \code{holdout.individuals.tsv}: tibble with the holdout individuals and
#' associated iteration and random seed number.
#' }
#' 
#' 
#' \strong{LOO: Leave-one-out}:
#'
#' Intermediate files are written, you can inspect specific iterations
#' and/or subsample:
#' \enumerate{
#' \item \code{assignment.random.results.iterations.raw.tsv}: tibble
#' with raw intermediate results of all iterations.
#' \item \code{markers.random.tsv}: tibble with the random markers selected for
#' each iteration with associated random seed number.
#' }
#' 
#' 
#' \strong{Folders}
#' 
#' The names have the different iterations \emph{i}
#' starting with \code{assignment_}\emph{i} contains:
#' \itemize{
#' \item \code{assignment_}\emph{i}\code{.tsv}: the assignment result, for the 
#' iteration.
#' \item \code{fst.ranked_}\emph{i}\code{.tsv}: for THL method, the ranked Fst 
#' per markers, for the iteration.
#' \item \code{gsi_sim_seeds}: the \code{gsi_sim} random seeds when this program 
#' is used, for the iteration.
#' }
#' The output in your global environment is a list. To view the assignment results
#' \code{$assignment} to view the ggplot2 figure \code{$plot.assignment}. 
#' See example below.


#' @export
#' @rdname assignment_ngs


#' @examples
#' \dontrun{
#' assignment.treefrog <- assignment_ngs(
#'     data = "batch_1.vcf", strata = "strata.treefrog.tsv",
#'     assignment.analysis = "gsi_sim",
#'      marker.number = c(500, 5000, "all"),
#'      markers.sampling = "ranked", thl = 0.3
#'      )
#' 
#' # To create a dataframe with the assignment results: 
#' assignment <- assignment.treefrog$assignment
#' 
#' # To plot the assignment using ggplot2 and facet 
#' fig <- assignment.treefrog$plot.assignment
#' 
#' # To view the full range of y values = Assignment success(%): 
#' fig + ggplot2::scale_y_continuous(limits = c(0,100)) 
#' 
#' # If you want to remove underscore in population names that contained white space:
#' facet_names <- c(
#'     `some_pop` = "Some POP", 
#'     `some_other_pop` = "This is what I want", 
#'     `OVERALL` = "Overall")
#' 
#' # use the labeller in the facet_grid or facet_wrap call:
#' fig + 
#'     ggplot2::facet_grid(
#'         SUBSAMPLE ~ CURRENT, 
#'         ggplot2::labeller = ggplot2::as_labeller(facet_names)
#'         ) + 
#'     ggplot2::scale_y_continuous(limits = c(0,100)) 
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

#' @importFrom radiator tidy_genomic_data change_pop_names write_genind detect_genomic_format
#' @importFrom adegenet optim.a.score dapc xvalDapc indNames pop predict.dapc

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com} and Eric C. Anderson

assignment_ngs <- function(
  data,
  strata = NULL,
  pop.levels = NULL,
  assignment.analysis = c("gsim_sim", "adegenet"),
  markers.sampling = c("ranked", "random"),
  marker.number = "all",
  thl = 1,
  iteration.method = 10,
  subsample = NULL,
  iteration.subsample = 1,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  ## testing
  # verbose = TRUE
  # parallel.core = parallel::detectCores() - 1
  # strata = NULL
  # pop.levels = NULL
  # adegenet.dapc.opt = "optim.a.score"
  # adegenet.n.rep = 30
  # adegenet.training = 0.9
  # folder = NULL
  # filename = "assignment_data.txt"
  # keep.gsi.files = FALSE
  # random.seed = NULL
  # whitelist.markers = NULL
  cat("################################################################################\n")
  cat("########################## assigner::assignment_ngs ############################\n")
  cat("################################################################################\n")
  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date/time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  res.list <- list() # results to keep stored in this list
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(cat("########################## assignment_ngs completed ############################\n"), add = TRUE)
  
  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator::radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE), 
    keepers = c("adegenet.dapc.opt", "adegenet.n.rep", "adegenet.training",
                "folder", "filename", "keep.gsi.files", "random.seed",
                "whitelist.markers"),
    deprecated = "sampling.method",
    verbose = FALSE
  )
  if (is.null(adegenet.dapc.opt)) adegenet.dapc.opt <- "optim.a.score"
  if (is.null(adegenet.n.rep)) adegenet.n.rep <- 30L
  if (is.null(adegenet.training)) adegenet.training <- 0.9
  if (is.null(filename)) filename <- "assignment_data.txt"
  if (is.null(keep.gsi.files)) keep.gsi.files <- FALSE
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  if (missing(assignment.analysis)) stop("assignment.analysis argument missing")
  if (missing(markers.sampling)) stop("markers.sampling argument missing")
  if (assignment.analysis == "gsi_sim" & !gsi_sim_exists()) {
    rlang::abort(
      "Can't find the gsi_sim executable where it was expected.\nFollow instruction for gsi_sim on assigner webpage")  
  }
  
  assignment.analysis <- match.arg(
    arg = assignment.analysis, 
    choices = c("gsi_sim", "adegenet"), 
    several.ok = FALSE)
  
  markers.sampling <- match.arg(
    arg = markers.sampling, 
    choices = c("ranked", "random"), 
    several.ok = FALSE) 
  
  message("Assignment analysis with ", assignment.analysis)
  
  # Correct iteration.method default if using only "all" in marker.number
  if (length(marker.number) == 1 && marker.number == "all" && markers.sampling == "random") {
    iteration.method <- 1
  }
  
  if (thl == "all") {
    message("Note: with thl == \"all\", automatically setting iteration.method = 1")
    iteration.method <- 1
  }
  
  # Folder----------------------------------------------------------------------
  directory <- radiator::generate_folder(
    f = folder,
    rad.folder = stringi::stri_join("assignment_analysis_method_", markers.sampling, sep = ""),
    prefix_int = FALSE,
    internal = FALSE,
    file.date = file.date,
    verbose = verbose)
  
  # write the dots file
  radiator::write_rad(
    data = rad.dots,
    path = directory,
    filename = stringi::stri_join(
      "assigner_assignment_ngs_args_", file.date, ".tsv"),
    tsv = TRUE,
    internal = FALSE,
    write.message = "Function call and arguments stored in: ",
    verbose = FALSE
  )
  
  # Import input ---------------------------------------------------------------
  input <- radiator::tidy_genomic_data(
    data = data, 
    whitelist.markers = whitelist.markers, 
    strata = strata,
    filename = NULL,
    verbose = FALSE,
    path.folder = directory
  )
  
  # Strata and pop levels ------------------------------------------------------
  input <- radiator::change_pop_names(data = input, pop.levels = pop.levels)
  strata <-  radiator::generate_strata(data = input, pop.id = TRUE)
  
  # Subsampling ----------------------------------------------------------------
  subsample.list <- generate_subsamples(
    subsample = subsample,
    iteration.subsample = iteration.subsample, 
    strata = strata,
    random.seed = random.seed,
    path.folder = directory
  )
  
  # Analysis -------------------------------------------------------------------
  res <- suppressWarnings(
    purrr::map_df(
      .x = subsample.list,
      .f = assignment_subsamples,
      input = input,
      subsample = subsample,
      assignment.analysis = assignment.analysis,
      marker.number = marker.number,
      markers.sampling = markers.sampling,
      iteration.method = iteration.method,
      filename = filename,
      directory = directory,
      keep.gsi.files = keep.gsi.files,
      verbose = verbose,
      parallel.core = parallel.core,
      random.seed = random.seed,
      thl = thl,
      adegenet.dapc.opt = adegenet.dapc.opt,
      adegenet.n.rep = adegenet.n.rep,
      adegenet.training = adegenet.training
    )
  )
  
  if (iteration.subsample > 1) {
    filename.res <- stringi::stri_join(
      "assignment.", markers.sampling, ".results.summary.stats.all.subsamples.tsv"
    )
  } else {
    filename.res <- "assignment.results.summary.stats.tsv"
  }
  
  readr::write_tsv(
    x = res, path = file.path(directory, filename.res),
    col_names = TRUE, append = FALSE
  )
  
  # Summary of the subsampling iterations
  if (iteration.subsample > 1) {
    res.pop <- dplyr::filter(.data = res, CURRENT != "OVERALL") %>%
      dplyr::group_by(CURRENT, MARKER_NUMBER, METHOD) %>%
      dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
      assigner_stats(x = .) %>% 
      dplyr::mutate(SUBSAMPLE = rep("OVERALL", n())) %>%
      dplyr::ungroup(.)
    
    res.overall <- suppressWarnings(
      res.pop %>% 
        dplyr::group_by(MARKER_NUMBER, METHOD) %>%
        dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
        assigner_stats(x = .) %>%
        dplyr::ungroup(.) %>% 
        dplyr::mutate(
          CURRENT = "OVERALL",
          SUBSAMPLE = "OVERALL"
        ) %>% 
        dplyr::bind_rows(res.pop) %>%
        dplyr::mutate(
          CURRENT = factor(CURRENT, levels = levels(res.pop$CURRENT), ordered = TRUE),
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE
        ) %>%
        dplyr::arrange(CURRENT, MARKER_NUMBER) %>%
        dplyr::select(CURRENT, MARKER_NUMBER, MEAN, MEDIAN, SE, MIN, MAX, 
                      QUANTILE25, QUANTILE75, SE_MIN, SE_MAX, METHOD, SUBSAMPLE)
    )
    res.pop <- NULL
    
    suppressWarnings(
      res %<>% 
        dplyr::mutate(SUBSAMPLE = as.character(SUBSAMPLE)) %>% 
        dplyr::bind_rows(res.overall) %>% 
        dplyr::mutate(
          SUBSAMPLE = factor(
            x = SUBSAMPLE, 
            levels = c(1:iteration.subsample, "OVERALL"), 
            ordered = TRUE
          )
        ) %>% 
        dplyr::arrange(CURRENT, MARKER_NUMBER, SUBSAMPLE)
    )
    res.overall <- NULL # unused objects
    filename <- "assignment.results.summary.stats.tsv"
    readr::write_tsv(
      x = res, path = file.path(directory, filename),
      col_names = TRUE, append = FALSE
    )
  } # End summary of the subsampling iterations
  
  # results in the list
  res.list$assignment <- res
  
  # Assignment plot ------------------------------------------------------------
  res.list$plot.assignment <- plot_assignment(x = res, path.folder = directory)
  
  # Return results --------------------------------------------------------------------
  return(res.list)
} # End assignment_ngs

# Internal Nested Functions ----------------    --------------------------------
# generate_subsamples ----------------------------------------------------------
#' @title generate_subsamples
#' @description generate_subsamples
#' @rdname generate_subsamples
#' @export
#' @keywords internal
generate_subsamples <- function(
  subsample, 
  iteration.subsample, 
  strata, 
  random.seed = NULL, 
  path.folder = NULL
) {
  if (!is.null(subsample)) {
    min.pop.n <- min(dplyr::count(x = strata, POP_ID, sort = TRUE)$n)
    # Control the subsample argument, replace if > than min sample size
    if (rlang::is_bare_numeric(subsample)) {
      if (subsample > min.pop.n) {
        message("Warning: subsample argument value > min sample size of some strata")
        message("Replacing subsample value by: ", min.pop.n)
        subsample <- min.pop.n
      }
    } else {
      # replace min by the min sample size found in the data
      if (subsample == "min") {
        subsample <- min.pop.n
        message("Using subsample size of: ", subsample)
      } else {
        rlang::abort("Wrong subsample value")
      }
    }
  }
  
  subsample.list <- purrr::map(
    .x = 1:iteration.subsample,
    .f = subsampling_data, # in utils.R
    strata = strata,
    subsample = subsample,
    random.seed = random.seed
  )
  
  # keep track of subsampling individuals and write to directory
  if (is.null(subsample)) {
    message("Subsampling: not selected")
  } else {
    message("Subsampling: selected")
    readr::write_tsv(
      x = dplyr::bind_rows(subsample.list), 
      path = file.path(path.folder, "subsampling_individuals.tsv"), 
      col_names = TRUE, 
      append = FALSE
    )
  } # End subsampling
  return(subsample.list)
}#End generate_subsamples

# assignment_subsamples-----------------------------------------------------------
#' @title assignment_subsamples
#' @description assignment_subsamples
#' @rdname assignment_subsamples
#' @export
#' @keywords internal
assignment_subsamples <- function(
  x,
  input = NULL,
  subsample = NULL,
  assignment.analysis = "gsi_sim",
  marker.number = NULL,
  markers.sampling = "random",
  iteration.method = 10,
  filename = "assignment_data.txt",
  directory = NULL,
  keep.gsi.files = FALSE,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1,
  random.seed = NULL,
  base.filename = NULL,
  thl = 1,
  adegenet.dapc.opt = NULL,
  adegenet.n.rep = NULL,
  adegenet.training = NULL
) {
  
  # x <- subsample.list[[1]] # test
  subsample.id <- unique(x$SUBSAMPLE)
  directory.subsample <- directory
  
  # Updating directories for subsampling
  if (!is.null(subsample)) {
    message("Analyzing subsample: ", subsample.id)
    directory.subsample <- file.path(directory, stringi::stri_join("subsample_", subsample.id))
    dir.create(directory.subsample)
  }
  
  # Keep only the subsample
  input %<>% dplyr::filter(INDIVIDUALS %in% x$INDIVIDUALS)
  
  # unused object
  x <- NULL
  
  # Keep a new strata
  # strata.df <- dplyr::distinct(.data = input, INDIVIDUALS, POP_ID)
  strata <- radiator::generate_strata(data = input, pop.id = TRUE)
  
  # If Adegenet, generate genind object
  genind.object <- NULL
  if (assignment.analysis == "adegenet") {
    message("Creating genind object")
    genind.object <- radiator::write_genind(data = input)
  }
  
  # Sampling of markers
  # unique list of markers after all the filtering
  unique.markers <- dplyr::distinct(.data = input, MARKERS)
  
  # clean and prepare marker.number argument
  marker.number %<>% clean_marker_number(x = ., unique.markers = unique.markers)
  
  # Random method ------------------------------------------------------------
  if (markers.sampling == "random") {
    message("Conducting Assignment analysis with markers selected randomly")
    # Number of times to repeat the sampling of markers
    iterations.list <- seq(iteration.method)
    # iterations.list <- 1:10 # test
    
    # Go through the function with the marker number selected
    message("Making a list containing all the markers combinations")
    if ("all" %in% marker.number) {# only 1 iterations when random method + max markers
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
    readr::write_tsv(
      x = tibble::as_tibble(dplyr::bind_rows(marker.random.list)),
      path = file.path(directory.subsample, "markers.random.tsv"),
      col_names = TRUE, append = FALSE
    )
    message("Starting parallel computations, for progress monitor activity in folder...")
    assignment.res <- NULL
    assignment.res <- suppressWarnings(
      .assigner_parallel_mc(#.assigner_parallel(
        X = marker.random.list,
        FUN = assignment_random,
        assignment.analysis = assignment.analysis,
        input = input,
        genind.object = genind.object,
        strata = strata,
        directory.subsample = directory.subsample,
        keep.gsi.files = keep.gsi.files,
        markers.sampling = markers.sampling,
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
    
    # Compiling the results
    message("Compiling results")
    if (assignment.analysis == "adegenet") {
      assignment.res <- suppressWarnings(
        dplyr::rename(.data = assignment.res, CURRENT = POP_ID) %>% 
          dplyr::mutate(SUBSAMPLE = subsample.id) %>% 
          dplyr::arrange(CURRENT, MARKER_NUMBER, ITERATIONS)
      )
    } else {
      assignment.res <- suppressWarnings(
        dplyr::mutate(.data = assignment.res, SUBSAMPLE = subsample.id) %>% 
          dplyr::arrange(CURRENT, INDIVIDUALS, MARKER_NUMBER, ITERATIONS)
      )
    }
    
    # Write the tables to directory
    # assignment results
    if (assignment.analysis == "gsi_sim") {
      if (is.null(subsample)) {
        filename.assignment.res <- stringi::stri_join(
          "assignment.", markers.sampling, ".results.iterations.raw.tsv"
        )
      } else {# with subsampling
        filename.assignment.res <- stringi::stri_join(
          "assignment.", markers.sampling, ".results.iterations.raw.subsample.", subsample.id, ".tsv"
        )
      }
      readr::write_tsv(
        x = assignment.res, 
        path = file.path(directory.subsample,filename.assignment.res), 
        col_names = TRUE, 
        append = FALSE
      )
    } else {# with adegenet
      if (is.null(subsample)) {
        filename.assignment.res <- stringi::stri_join(
          "assignment.", markers.sampling, ".results.iterations.tsv"
        )
      } else {# with subsampling
        filename.assignment.res <- stringi::stri_join(
          "assignment.", markers.sampling, ".results.iterations.subsample.", subsample.id, ".tsv"
        )
      }
      readr::write_tsv(
        x = assignment.res, 
        path = file.path(directory.subsample,filename.assignment.res), 
        col_names = TRUE, append = FALSE
      )
    }
    
    if (assignment.analysis == "gsi_sim") {
      assignment.stats.pop <- suppressWarnings(
        assignment.res %>%
          dplyr::group_by(CURRENT, INFERRED, MARKER_NUMBER, ITERATIONS, METHOD) %>%
          dplyr::tally(.) %>%
          dplyr::group_by(CURRENT, MARKER_NUMBER, ITERATIONS, METHOD) %>%
          dplyr::mutate(TOTAL = sum(n)) %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate(ASSIGNMENT_PERC = round(n/TOTAL*100, 0)) %>%
          # here ASSIGNMENT_PERC replaces MEAN_i used previously
          dplyr::filter(as.character(CURRENT) == as.character(INFERRED)) %>%
          dplyr::select(CURRENT, ASSIGNMENT_PERC, MARKER_NUMBER, ITERATIONS, METHOD) %>%
          dplyr::group_by(CURRENT, MARKER_NUMBER, METHOD) %>%
          assigner_stats(x = .) %>% 
          dplyr::ungroup(.) %>% 
          dplyr::arrange(CURRENT, MARKER_NUMBER)
      )
    }
    
    if (assignment.analysis == "adegenet") {
      assignment.stats.pop <- suppressWarnings(
        assignment.res %>%
          dplyr::group_by(CURRENT, MARKER_NUMBER, METHOD) %>%
          assigner_stats(x = .) %>%
          dplyr::ungroup(.) %>%
          dplyr::arrange(CURRENT, MARKER_NUMBER)
      )
    }
    
    # Next step is common for gsi_sim and adegenet
    assigner.levels <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
    
    assignment.stats.overall <- assignment.stats.pop %>%
      dplyr::group_by(MARKER_NUMBER, METHOD) %>%
      dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
      assigner_stats(x = .) %>%
      dplyr::mutate(CURRENT = rep("OVERALL", n())) %>%
      dplyr::ungroup(.) %>% 
      dplyr::arrange(CURRENT, MARKER_NUMBER)
    
    assignment.summary.stats <- suppressWarnings(
      dplyr::bind_rows(assignment.stats.pop, assignment.stats.overall) %>%
        dplyr::mutate(CURRENT = factor(CURRENT, levels = assigner.levels, ordered = TRUE)) %>%
        dplyr::arrange(CURRENT, MARKER_NUMBER) %>%
        dplyr::mutate(
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE,
          ITERATIONS = rep(iteration.method, n())
        ) %>%
        dplyr::select(CURRENT, MARKER_NUMBER, MEAN, MEDIAN, SE, MIN, MAX, 
                      QUANTILE25, QUANTILE75, SE_MIN, SE_MAX, METHOD, ITERATIONS)
    )
    assignment.stats.pop <- assignment.stats.overall <- NULL
    # update the assignment with subsampling iterations id
    assignment.summary.stats %<>% dplyr::mutate(SUBSAMPLE = subsample.id)
    # assignment summary stats
    if (!is.null(subsample)) {
      filename.assignment.sum <- stringi::stri_join(
        "assignment.", markers.sampling, ".results.iterations.summary.subsample.", subsample.id, ".tsv"
      )
      readr::write_tsv(
        x = assignment.summary.stats,
        path = file.path(directory.subsample, filename.assignment.sum),
        col_names = TRUE, append = FALSE
      )
    }
  } # End method random
  
  # Ranked method ------------------------------------------------------------
  if (markers.sampling == "ranked") {
    message("Conducting Assignment analysis using Training, Holdout, Leave-one-out")
    message("Using training samples to rank markers based on Fst")
    
    # individuals and pop df
    strata <- radiator::generate_strata(data = input, pop.id = TRUE)
    
    # Will go through the individuals in the list one by one.
    hs <- generate_holdhout(
      thl = thl, 
      strata = strata, 
      iteration.method = iteration.method, 
      random.seed = random.seed, 
      path.folder = directory.subsample
    )
    holdout.individuals <- hs$holdout.individuals
    iterations.list <- hs$iterations.list
    hs <- NULL
    message("Holdout samples saved in your folder")
    
    # Going through the loop of holdout individuals
    message("Starting parallel computations, for progress monitor activity in folder...")
    assignment.res.summary <- list()
    # assignment.res <- .assigner_parallel(
    assignment.res.summary <- suppressWarnings(
      .assigner_parallel_mc(
        X = iterations.list, 
        FUN = assignment_ranking, 
        thl = thl,
        input = input,
        holdout.individuals = holdout.individuals,
        directory.subsample = directory.subsample,
        marker.number = marker.number,
        assignment.analysis = assignment.analysis,
        genind.object = genind.object,
        strata = strata,
        markers.sampling = markers.sampling,
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
      ) %>% 
        dplyr::bind_rows() %>% 
        dplyr::mutate(SUBSAMPLE = subsample.id) %>% 
        dplyr::arrange(CURRENT, INDIVIDUALS, MARKER_NUMBER)
    )
    
    # assignment results
    if (is.null(subsample)) {
      filename.assignment.res <- stringi::stri_join(
        "assignment.", 
        markers.sampling,
        ".results.iterations.raw.tsv")
    } else {# with subsampling
      filename.assignment.res <- stringi::stri_join(
        "assignment.", 
        markers.sampling,
        ".results.iterations.raw.subsample.",
        subsample.id, ".tsv")
    }
    readr::write_tsv(
      x = assignment.res.summary,
      path = file.path(directory.subsample, filename.assignment.res),
      col_names = TRUE, append = FALSE
    )
    
    
    if (thl == 1 | thl == "all") {
      assignment.stats.pop <- assignment.res.summary %>%
        dplyr::group_by(CURRENT, MARKER_NUMBER, METHOD) %>%
        dplyr::summarise(
          n = length(CURRENT[as.character(CURRENT) == as.character(INFERRED)]),
          TOTAL = length(CURRENT)
        ) %>%
        dplyr::ungroup(.) %>% 
        dplyr::mutate(MEAN = round(n / TOTAL * 100, 0)) %>% 
        dplyr::select(-n, -TOTAL)
      
      assigner.levels <- c(unique(assignment.stats.pop$CURRENT), "OVERALL")
      
      assignment.stats.overall <- assignment.stats.pop %>% 
        dplyr::group_by(MARKER_NUMBER, METHOD) %>%
        dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
        assigner_stats(x = .) %>% 
        dplyr::mutate(CURRENT = rep("OVERALL", n())) %>% 
        dplyr::arrange(CURRENT, MARKER_NUMBER)
      
      assignment.summary.stats <- suppressWarnings(
        dplyr::bind_rows(assignment.stats.pop, assignment.stats.overall) %>%
          dplyr::mutate(CURRENT = factor(CURRENT, levels = assigner.levels, ordered = TRUE)) %>%
          dplyr::arrange(CURRENT, MARKER_NUMBER) %>%
          dplyr::mutate(
            SE_MIN = MEAN - SE,
            SE_MAX = MEAN + SE
          )
      )
      
    } else {
      assignment.res.summary.prep <- assignment.res.summary %>% 
        dplyr::group_by(CURRENT, MARKER_NUMBER, METHOD, ITERATIONS) %>%
        dplyr::summarise(
          n = length(CURRENT[as.character(CURRENT) == as.character(INFERRED)]),
          TOTAL = length(CURRENT)
        ) %>%
        dplyr::ungroup(.) %>% 
        dplyr::mutate(ASSIGNMENT_PERC = round(n/TOTAL*100, 0)) %>% 
        dplyr::select(-n, -TOTAL)
      
      if (is.null(subsample)) {
        filename.assignment.res.sum <- stringi::stri_join(
          "assignment.", markers.sampling, ".results.iterations.summary.tsv")
      } else {# with subsampling
        filename.assignment.res.sum <- stringi::stri_join(
          "assignment.", markers.sampling,
          ".results.iterations.summary.subsample.", subsample.id, ".tsv")
      }
      readr::write_tsv(
        x = assignment.res.summary.prep,
        path = file.path(directory.subsample,filename.assignment.res.sum),
        col_names = TRUE, append = FALSE
      )
      
      assignment.stats.pop <- assignment.res.summary.prep %>%
        dplyr::group_by(CURRENT, MARKER_NUMBER, METHOD) %>%
        assigner_stats(x = .) %>%
        dplyr::ungroup(.) %>% 
        dplyr::arrange(CURRENT, MARKER_NUMBER)
      assignment.res.summary.prep <- NULL
      assigner.levels <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
      
      assignment.stats.overall <- assignment.stats.pop %>%
        dplyr::group_by(MARKER_NUMBER, METHOD) %>%
        dplyr::rename(ASSIGNMENT_PERC = MEAN) %>%
        assigner_stats(x = .) %>%
        dplyr::mutate(CURRENT = rep("OVERALL", n())) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(CURRENT, MARKER_NUMBER)
      
      assignment.summary.stats <- suppressWarnings(
        dplyr::bind_rows(assignment.stats.pop, assignment.stats.overall) %>%
          dplyr::mutate(CURRENT = factor(CURRENT, levels = assigner.levels, ordered = TRUE)) %>%
          dplyr::arrange(CURRENT, MARKER_NUMBER)
      )
      assignment.stats.overall <- assignment.stats.pop <- NULL
    } # End thl != 1
    
    # update the assignment with subsampling iterations id
    assignment.summary.stats %<>% dplyr::mutate(SUBSAMPLE = subsample.id)
    # assignment summary stats
    if (!is.null(subsample)) {
      # filename.assignment.sum <- stringi::stri_join("assignment", markers.sampling, "results", "summary.stats", "tsv", sep = ".")
      # } else {# with subsampling
      filename.assignment.sum <- stringi::stri_join(
        "assignment.", markers.sampling, ".results.summary.stats.subsample.",
        subsample.id, ".tsv")
      readr::write_tsv(
        x = assignment.summary.stats,
        path = file.path(directory.subsample,filename.assignment.sum),
        col_names = TRUE, append = FALSE
      )
    }
  } # End of ranked thl method
  
  return(assignment.summary.stats)
}# End assignment_subsamples

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
  strata = NULL,
  markers.sampling = NULL,
  subsample.id = NULL,
  adegenet.dapc.opt = NULL,
  adegenet.n.rep = NULL,
  adegenet.training = NULL,
  filename = NULL,
  keep.gsi.files = NULL
) {
  i <- iterations.list
  # i <- iterations.list[[1]]
  
  # folder of iterations
  i.folder <- radiator::generate_folder(
    f = directory.subsample,
    rad.folder = stringi::stri_join("assignment_", i),
    prefix_int = FALSE,
    internal = FALSE, 
    append.date = FALSE,
    verbose = FALSE)
  
  # Ranking Fst with training dataset (keep holdout individuals out)
  message("Ranking markers based on Fst")
  if (thl == "all") holdout <- NULL
  if (thl == 1) holdout <- tibble::tibble(INDIVIDUALS = i) %$% INDIVIDUALS
  
  # thl proportion or > 1
  if (thl != 1 && thl != "all") {
    holdout <- dplyr::filter(.data = holdout.individuals, ITERATIONS %in% i) %$% 
      INDIVIDUALS
  }
  
  fst.ranked <- assigner::fst_WC84(
    data = input, 
    strata = NULL, 
    holdout.samples = holdout
  ) %$%
    fst.ranked
  
  readr::write_tsv(
    x = fst.ranked, 
    path = file.path(i.folder, stringi::stri_join("fst.ranked_", i, ".tsv", sep = "")),
    col_names = TRUE, 
    append = FALSE
  )
  
  # looping through the markers
  message("Going through the marker.number")
  assignment.res <- list()
  assignment.res <- suppressWarnings(
    purrr::map_df(
      .x = marker.number, 
      .f = assignment_marker_loop,
      assignment.analysis = assignment.analysis,
      fst.ranked = fst.ranked,
      i = i,
      input = input,
      genind.object = genind.object,
      strata = strata,
      markers.sampling = markers.sampling,
      subsample.id = subsample.id,
      holdout = holdout,
      thl = thl,
      adegenet.dapc.opt = adegenet.dapc.opt,
      adegenet.n.rep = adegenet.n.rep,
      adegenet.training = adegenet.training,
      directory.subsample = i.folder,
      filename = filename,
      keep.gsi.files = keep.gsi.files
    )
  )
  message("Summarizing the assignment analysis results by iterations and marker group")
  readr::write_tsv(
    x = assignment.res,
    path = file.path(i.folder, stringi::stri_join("assignment_", i, ".tsv")), 
    col_names = TRUE, 
    append = FALSE
  )
  holdout <- fst.ranked <- fst.ranked.imp <- i <- NULL
  return(assignment.res)
}  # End assignment_ranking

# assignment_marker_loop--------------------------------------------------------
#' @title assignment_marker_loop
#' @description assignment_marker_loop
#' @rdname assignment_marker_loop
#' @export
#' @keywords internal
assignment_marker_loop <- function(
  m, # marker.number
  assignment.analysis = "gsi_sim",
  fst.ranked = NULL,
  i = NULL,
  input = NULL,
  genind.object = NULL,
  strata = NULL,
  markers.sampling = NULL,
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
  # m <- marker.number[1]
  m <- as.numeric(m)
  message("Marker number: ", m)
  
  select.markers <- dplyr::filter(.data = fst.ranked, RANKING <= m) %$% MARKERS
  
  # Modify filename
  filename <- file.path(directory.subsample, filename)
  filename <- stringi::stri_replace_all_fixed(
    filename,
    pattern = ".txt",
    replacement = stringi::stri_join(
      "", "_iteration_", i, "_markers_", m, ".txt")
  )
  
  if (assignment.analysis == "gsi_sim") {
    assignment <- assignment_gsi_sim(
      input = input,
      strata = strata,
      select.markers = select.markers,
      i = i,
      m = m,
      holdout = holdout,
      filename = filename,
      directory.subsample = directory.subsample,
      keep.gsi.files = keep.gsi.files,
      markers.sampling = markers.sampling,
      thl = thl
    )
  }
  
  
  if (assignment.analysis == "adegenet") {
    assignment <- assignment_adegenet(
      data = genind.object,
      select.markers = select.markers,
      adegenet.dapc.opt = adegenet.dapc.opt,
      adegenet.n.rep = adegenet.n.rep,
      adegenet.training = adegenet.training,
      i = i, 
      m = m,
      markers.sampling = markers.sampling,
      subsample.id = subsample.id,
      holdout = holdout
    )
  }
  return(assignment)
}#End assignment_marker_loop

# assignment_gsi_sim------------------------------------------------------------
#' @title assignment with gsi_sim
#' @description assignment with gsi_sim
#' @rdname assignment_gsi_sim
#' @export
#' @keywords internal
assignment_gsi_sim <- function(
  input = NULL,
  strata = NULL,
  select.markers = NULL,
  i = NULL,
  m = NULL,
  holdout = NULL,
  filename = NULL,
  directory.subsample = NULL,
  keep.gsi.files = FALSE,
  markers.sampling = "random",
  thl = 1
) {
  # Write gsi_sim input file to directory
  input.gsi <- radiator::write_gsi_sim(
    data = input %>% 
      dplyr::filter(MARKERS %in% select.markers) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS), 
    filename = filename
  )
  
  # Run gsi_sim
  output.gsi <- stringi::stri_replace_all_fixed(
    str = input.gsi, pattern = "txt", replacement = "output.txt")
  setwd(directory.subsample)
  system(
    paste(assigner::gsi_sim_binary(), "-b", input.gsi, "--self-assign > ", 
          output.gsi)
  )
  
  # Option remove the input file from directory to save space
  if (!keep.gsi.files) file.remove(input.gsi)
  
  # Import GSI_SIM results
  assignment <- suppressWarnings(
    suppressMessages(readr::read_delim(output.gsi, col_names = "ID", delim = "\t")) %>%
      tidyr::separate(ID, c("KEEPER", "ASSIGN"), sep = ":/", extra = "warn") %>%
      dplyr::filter(KEEPER == "SELF_ASSIGN_A_LA_GC_CSV") %>%
      tidyr::separate(ASSIGN, c("INDIVIDUALS", "ASSIGN"), sep = ";", extra = "merge") %>%
      tidyr::separate(ASSIGN, c("INFERRED", "OTHERS"), sep = ";", convert = TRUE, 
                      numerals = "no.loss", extra = "merge") %>%
      tidyr::separate(OTHERS, c("SCORE", "OTHERS"), sep = ";;", convert = TRUE, 
                      numerals = "no.loss", extra = "merge") %>%
      tidyr::separate(OTHERS, c("SECOND_BEST_POP", "OTHERS"), sep = ";", 
                      convert = TRUE, numerals = "no.loss", extra = "merge") %>%
      tidyr::separate(OTHERS, c("SECOND_BEST_SCORE", "OTHERS"), sep = ";;", 
                      convert = TRUE, numerals = "no.loss") %>% 
      dplyr::mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>% 
      dplyr::left_join(strata, by = "INDIVIDUALS") %>%
      dplyr::rename(CURRENT = POP_ID) %>% 
      dplyr::mutate(
        SCORE = round(SCORE, 2),
        SECOND_BEST_SCORE = round(SECOND_BEST_SCORE, 2),
        MARKER_NUMBER = as.numeric(rep(m, n())), # m: Number of markers
        METHOD = rep(markers.sampling, n())
      ) %>%
      dplyr::select(INDIVIDUALS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, 
                    SECOND_BEST_SCORE, MARKER_NUMBER, METHOD) %>%
      dplyr::arrange(CURRENT)
  )
  
  if (markers.sampling == "random") {
    assignment %<>% 
      dplyr::mutate(ITERATIONS = i) %>% 
      dplyr::select(INDIVIDUALS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, 
                    SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, ITERATIONS) %>%
      dplyr::arrange(CURRENT)
  }
  
  if (markers.sampling == "ranked" & thl != "all") {
    assignment %<>% dplyr::filter(INDIVIDUALS %in% holdout)
  }
  
  # Option remove the output file from directory to save space
  if (!keep.gsi.files) file.remove(output.gsi)
  
  # saving preliminary results
  if (markers.sampling == "ranked") {
    
    assignment %<>% 
      dplyr::mutate(METHOD = rep(stringi::stri_join("ranked_thl_", thl) , n()))
    
    if (thl != 1 & thl != "all") {
      assignment %<>% dplyr::mutate(ITERATIONS = rep(i, n()))
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
  adegenet.dapc.opt = "optim.a.score",
  adegenet.n.rep = 30,
  adegenet.training = 0.9, 
  i = NULL,
  m = NULL,
  holdout = NULL,
  markers.sampling = "random",
  thl = 1,
  subsample.id = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  # data <- genind.object #test
  data.select <- data[loc = select.markers]
  
  # Run adegenet
  pop.data <- data.select@pop
  pop.data <- droplevels(pop.data)
  
  # check for missing data 
  if (anyNA(data.select@tab)) {
    message("Missing data detected: optim.a.score automatically selected")
    adegenet.dapc.opt <- "optim.a.score"
  }
  
  if (markers.sampling == "random") {
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
      dapc.assignment <- adegenet::dapc(
        data.select, n.da = length(levels(pop.data)),
        n.pca = dapc.best.optim.a.score, 
        pop = pop.data)
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
  
  if (markers.sampling == "ranked") {
    
    # Alpha-Score DAPC training data
    training.data <- data.select[!adegenet::indNames(data.select) %in% holdout] # training dataset
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
    holdout.data <- data.select[adegenet::indNames(data.select) %in% holdout] # holdout dataset
    pop.holdout <- holdout.data@pop
    pop.holdout <- droplevels(pop.holdout)
    assignment.levels <- levels(pop.holdout) # for figure
    rev.assignment.levels <- rev(assignment.levels)  # for figure 
    
    dapc.predict.holdout <- adegenet::predict.dapc(dapc.training, newdata = holdout.data)
    message("Assigning holdout data for iteration: ", i)
  }
  
  
  # Get Assignment results
  if (markers.sampling == "random") {
    assignment <- tibble::tibble(ASSIGNMENT_PERC = summary(dapc.assignment)$assign.per.pop*100) %>% 
      dplyr::bind_cols(tibble::tibble(POP_ID = levels(pop.data))) %>%
      dplyr::mutate(ASSIGNMENT_PERC = round(ASSIGNMENT_PERC, 2)) %>% 
      dplyr::select(POP_ID, ASSIGNMENT_PERC)
  }        
  if (markers.sampling == "ranked") {
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
  
  assignment %<>%
    dplyr::mutate(
      METHOD = rep(markers.sampling, n()),
      ITERATIONS = rep(i, n()),
      MARKER_NUMBER = as.numeric(rep(m, n()))
    )
  
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
  if (is.null(random.seed)) random.seed <- sample(x = 1:1000000, size = 1)
  set.seed(random.seed)
  
  select.markers <- dplyr::sample_n(tbl = unique.markers, size = m, replace = FALSE) %>%
    dplyr::arrange(MARKERS) %>%
    dplyr::mutate(
      ITERATIONS = iteration.method,
      MARKER_NUMBER = m,
      RANDOM_SEED_NUMBER = random.seed
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
  strata = NULL,
  directory.subsample = NULL,
  keep.gsi.files = FALSE,
  markers.sampling = "random",
  subsample.id = NULL,
  filename = NULL,
  adegenet.n.rep = 30,
  adegenet.training = 0.9,
  holdout = NULL
) {
  x <- tibble::as_tibble(x)
  # x <- marker.random.list[[1]]
  i <- as.integer(unique(x$ITERATIONS))      # iteration
  m <- as.numeric(unique(x$MARKER_NUMBER))   # number of marker selected
  
  select.markers <- dplyr::ungroup(x) %$% MARKERS
  
  # folder of iterations
  i.folder <- radiator::generate_folder(
    f = directory.subsample,
    rad.folder = stringi::stri_join("assignment_", i),
    prefix_int = FALSE,
    internal = FALSE, 
    append.date = FALSE,
    verbose = FALSE)
  
  # Modify filename
  filename <- file.path(i.folder, filename)
  filename <- stringi::stri_replace_all_fixed(
    filename,
    pattern = ".txt",
    replacement = stringi::stri_join(
      "", "_iteration_", i, "_markers_", m, ".txt")
  )
  
  if (assignment.analysis == "gsi_sim") {
    assignment <- assignment_gsi_sim(
      input = input,
      strata = strata,
      select.markers = select.markers,
      i = i, 
      m = m,
      holdout = NULL,
      filename = filename,
      directory.subsample = i.folder,
      keep.gsi.files = keep.gsi.files,
      markers.sampling = markers.sampling
    )
  }
  if (assignment.analysis == "adegenet") {
    assignment <- assignment_adegenet(
      data = genind.object,
      select.markers = select.markers,
      adegenet.dapc.opt = "optim.a.score",
      adegenet.n.rep = adegenet.n.rep, 
      adegenet.training = adegenet.training,
      i = i, 
      m = m,
      holdout = NULL
    )
  }
  message("Summarizing the assignment analysis results by iterations and marker group")
  readr::write_tsv(
    x = assignment,
    path = file.path(i.folder, stringi::stri_join("assignment_", i, ".tsv")), 
    col_names = TRUE, 
    append = FALSE
  )
  
  # unused objects
  x <- i <- m <- select.markers <- filename <- NULL
  return(assignment)
}#End assignment_random

# assigner_stats--------------------------------------------------------------
#' @title assigner_stats
#' @description assigner_stats
#' @rdname assigner_stats
#' @export
#' @keywords internal
assigner_stats <- function(x) {
  x %<>%
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
    )
}

# plot_assignment --------------------------------------------------------------
#' @title plot_assignment
#' @description plot_assignment
#' @rdname plot_assignment
#' @export
#' @keywords internal
plot_assignment <- function(x, path.folder = NULL) {
  if (is.null(path.folder)) path.folder <- getwd()
  
  plot.assignment <- ggplot2::ggplot(
    x, ggplot2::aes(x = factor(MARKER_NUMBER), y = MEAN)) +
    ggplot2::geom_point(size = 2, alpha = 0.5, na.rm = TRUE) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = SE_MIN, ymax = SE_MAX), width = 0.3, na.rm = TRUE) +
    ggplot2::scale_y_continuous(
      breaks = c(0, 10, 20 ,30, 40, 50, 60, 70, 80, 90, 100)) +
    ggplot2::labs(x = "Marker number",
                  y = "Assignment success (%)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",      
      panel.grid.minor.x = ggplot2::element_blank(), 
      panel.grid.major.y = ggplot2::element_line(colour = "grey60", 
                                                 linetype = "dashed"), 
      axis.title.x = ggplot2::element_text(size = 10, face = "bold"), 
      axis.text.x = ggplot2::element_text(size = 8, face = "bold", angle = 90, 
                                          hjust = 1, vjust = 0.5), 
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"), 
      axis.text.y = ggplot2::element_text(size = 10, face = "bold")
    ) +
    ggplot2::facet_grid(SUBSAMPLE~CURRENT)
  
  n.pop <- length(unique(x$CURRENT))
  n.sub <- length(unique(x$SUBSAMPLE))
  
  # To save the plot:
  ggplot2::ggsave(
    plot = plot.assignment, 
    filename = file.path(path.folder, "assignment.plot.pdf"), 
    height = 10 * n.sub, width = 5 * n.pop, dpi = 300, 
    units = "cm", useDingbats = FALSE, limitsize = FALSE
  )
  return(plot.assignment)
}

# clean_marker_number ----------------------------------------------------------
#' @title clean_marker_number
#' @description clean_marker_number
#' @rdname clean_marker_number
#' @export
#' @keywords internal
clean_marker_number <- function(x, unique.markers) {
  # x is the marker.number
  # if "all" is present in the list, change to the maximum number of markers
  x <- as.numeric(
    stringi::stri_replace_all_fixed(str = x, 
                                    pattern = "all", 
                                    replacement = nrow(unique.markers), 
                                    vectorize_all = TRUE)
  )
  
  # In marker.number, remove marker group higher than the max number of markers
  removing.marker <- purrr::keep(.x = x, .p = x > nrow(unique.markers))
  
  if (length(removing.marker) > 0) {
    message(
      "Removing marker.number higher than the max number of markers: ", 
      stringi::stri_join(removing.marker, collapse = ", ")
    )
    x <- purrr::discard(.x = x, .p = x > nrow(unique.markers))
  }
  return(x)
}#End clean_marker_number

# generate_holdhout ----------------------------------------------------------
#' @title generate_holdhout
#' @description Generate holdhout samples based on thl and iterations.list
#' @rdname generate_holdhout
#' @export
#' @keywords internal
generate_holdhout <- function(
  thl, 
  strata, 
  iteration.method, 
  random.seed = NULL, 
  path.folder = NULL
) {
  res <- list()
  # Set seed for random sampling
  if (is.null(random.seed)) random.seed <- sample(x = 1:1000000, size = 1)
  set.seed(random.seed)
  
  # required func
  pick_samples <- function(x, strata, random.seed = NULL, thl, frac = FALSE) {
    holdout.individuals <- dplyr::group_by(strata, POP_ID)
    
    if (frac) {
      holdout.individuals %<>%
        dplyr::sample_frac(tbl = ., size = thl, replace = FALSE)
    } else {
      holdout.individuals %<>%
        dplyr::sample_n(tbl = ., size = thl, replace = FALSE)
    }
    
    holdout.individuals %<>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS) %>%
      dplyr::select(INDIVIDUALS) %>%
      dplyr::mutate(
        ITERATIONS = x,
        RANDOM_SEED_NUMBER = random.seed
      )
    return(holdout.individuals)
  }#End pick_samples
  
  # Generate holdout and iterations.list
  if (rlang::is_bare_integerish(thl)) {
    if (thl == 1) {
      res$iterations.list <- strata$INDIVIDUALS
      res$holdout.individuals <- strata %>%
        dplyr::mutate(ITERATIONS = stringi::stri_join("HOLDOUT", seq(1:nrow(strata)), sep = "_"))
    } else {# Generate holdhout samples using integer
      res$iterations.list <- 1:iteration.method
      res$holdout.individuals <- purrr::map_dfr(
        .x = res$iterations.list, 
        .f = pick_samples,
        strata = strata, 
        random.seed = random.seed, 
        thl = thl, 
        frac = FALSE
      )
    }
  } else {
    # no holdout for that one
    if (thl == "all") {
      res$iterations.list <- iteration.method
      res$holdout.individuals <- NULL
      message("Warning: using all the individuals for ranking markers based on Fst\nNo holdout samples")
      message("Recommended reading: \nAnderson, E. C. (2010) Assessing the power of informative subsets of
              loci for population assignment: standard methods are upwardly biased.\nMolecular ecology resources 10, 4:701-710.")
    } else {# Generate holdhout samples using proportion
      res$iterations.list <- 1:iteration.method
      res$holdout.individuals <- purrr::map_dfr(
        .x = res$iterations.list, 
        .f = pick_samples,
        strata = strata, 
        random.seed = random.seed, 
        thl = thl, 
        frac = TRUE
      )
    }
  }# End tracking holdout individuals
  
  readr::write_tsv(
    x = res$holdout.individuals, 
    path = file.path(path.folder, "holdout.individuals.tsv"), 
    col_names = TRUE, 
    append = FALSE
  )
  return(res)
}#End generate_holdhout
