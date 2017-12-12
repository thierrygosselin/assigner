#' @name assignment_mixture

#' @title Mixture/Baseline assignment analysis tailored for RADseq data

#' @description
#' The arguments in the \code{assignment_mixture} function were tailored for the
#' reality of GBS/RADseq data for assignment analysis of mixture or unknown samples, 
#' while maintaining a reproducible workflow. Assignment are conducted
#' using \code{gsi_sim} or \code{\link[adegenet]{adegenet}}.
#' 
#' \itemize{
#'   \item \strong{Input file:} various file format are supported (see \code{data} argument below)
#'   \item \strong{Filters:} genotypes, markers, individuals and populations can be 
#'   filtered and/or selected in several ways using blacklist,
#'   whitelist and other arguments.
#'   \item \strong{Cross-Validations:} Markers can be randomly selected for a classic LOO (Leave-One-Out)
#'   assignment or chosen based on ranked Fst for a thl
#'   (Training, Holdout, Leave-one-out) assignment analysis.
#'   Baseline = training set, Mixture/unknown = holdout/test set.
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
#' @inheritParams assignment_ngs 

#' @param strata (optional/required) Required for VCF and haplotypes files, 
#' optional for the other file formats supported.
#' 
#' \strong{mixture samples}: use \code{"mixture"} or \code{"unknown"} inside the
#' \code{STRATA} column to identify samples not in the baseline,
#' but see details on the other ways to identify your unknown or mixture samples.
#' 
#' The strata file is a tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping. 
#' To create a strata file see \code{\link[radiator]{individuals2strata}}.
#' If you have already run 
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data, 
#' the strata file is similar to a stacks `population map file`, make sure you 
#' have the required column names  (\code{INDIVIDUALS} and \code{STRATA}).
#' Default: \code{strata = NULL}.

#' @param mixture (optional) A file in the working directory (e.g. "mixture.txt")
#' with mixture individual ID. The column header is \code{INDIVIDUALS}.
#' Default: \code{mixture = NULL}, but see details on the different ways to 
#' identify your unknown or mixture samples.


#' @param sampling.method (character) Should the markers be randomly selected
#' \code{sampling.method = "random"}, for a classic Leave-One-Out (LOO)
#' assignment or chosen based on ranked Fst \code{sampling.method = "ranked"},
#' using the baseline samples for the training and the mixture samples as holdout. 
#' Classic Leave-one-out is then used to assign individual mixture samples.

#' @param iteration.method With random marker selection the iterations argument =
#' the number of iterations to repeat marker resampling.
#' With \code{marker.number = c(500, 1000)} and default iterations setting,
#' 500 markers will be randomly chosen 10 times and 1000 markers will be randomly
#' chosen 10 times.
#' 
#' \strong{Notes:} If all the markers are used, with \code{marker.number = "all"}
#' or in a series of marker number groupings \code{marker.number = c(200, 500, "all")},
#' the number of iteration is automatically set to 1. The remaining groupings
#' are unaffected.
#' Default: \code{iteration.method = 10}.

#' @param folder (optional) The name of the folder created in the working 
#' directory to save the files/results. Default: \code{folder = NULL}.

#' @param filename (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. 
#' Default \code{filename = assignment_data.txt}.
#' The number of markers used will be appended to the name of the file.

#' @param subsample (Integer or Character, optional)
#' With \code{subsample = 36}, 36 individuals in the \strong{baseline} of
#' each populations are chosen
#' randomly to represent the dataset. This integer as to be smaller than the
#' population with min sample size, if higher, the minimum sample size found
#' in the data will be used as default. In doubt, use \code{subsample = "min"},
#' the function will use the smallest population sample size found in the data.
#' Default: \code{subsample = NULL} (no subsampling).


#' @param iteration.subsample (Integer) The number of iterations to repeat 
#' subsampling, default: \code{iteration.subsample = 1}.
#' With \code{subsample = 20} and \code{iteration.subsample = 10},
#' 20 individuals/populations will be randomly chosen 10 times.


#' @inheritParams radiator::radiator_imputations_module 

#' @param impute.mixture (Logical) Imputations of mixture samples.
#' Default: \code{impute.mixture = FALSE}. For no imputation. 
#' For \code{impute.mixture = TRUE} the hierarchical.levels (see below)
#' for the mixture samples is automatically set to 
#' \code{hierarchical.levels = "global"}. Warning: bias could be introduced by imputing
#' missing genotype in the mixture samples.

#' @details
#' 
#' \strong{2 options to identify your unknown or mixture samples, using:}
#' \itemize{
#'   \item \strong{strata} argument with "mixture" or "unknown" 
#'   instead of the grouping id (e.g. populations) for the other samples.
#'    \item \strong{mixture} argument using a file in the working directory 
#'    (e.g. "mixture.txt") with mixture individual ID. 
#'    The column header is \code{INDIVIDUALS}.
#'   }
#' 
#' \strong{Input files:}
#' \strong{Input files:} see \pkg{radiator} \code{\link[radiator]{tidy_genomic_data}}
#' for detailed information about supported file format.
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
#' \href{https://github.com/thierrygosselin/radiator/blob/master/vignettes/vignette_imputations_parallel.Rmd}{vignette}
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
#' \code{$assignment} or \code{$assignment.mixture.summary.results}

#' @note \code{assignment_mixture} assumes that the command line version of gsi_sim 
#' is properly installed into \code{file.path(system.file(package = "assigner"), "bin", "gsi_sim")}.
#' Things are set up so that it will try running gsi_sim, and if it does not find it, the 
#' program will throw an error and ask the user to run \code{\link{install_gsi_sim}}
#' which will do its best to put a usable copy of gsi_sim where it is needed. To do 
#' so, you must be connected to the internet. If that doesn't work, you will
#' need to compile the program yourself, or get it yourself, and the manually copy
#' it to \code{file.path(system.file(package = "assigner"), "bin", "gsi_sim")}.
#' To compile gsi_sim, follow the 
#' instruction here: \url{https://github.com/eriqande/gsi_sim}.

#' @export
#' @rdname assignment_mixture
#' @importFrom parallel detectCores
#' @importFrom stringi stri_join stri_sub stri_replace_all_fixed stri_detect_fixed stri_replace_na
#' @importFrom adegenet optim.a.score predict.dapc dapc indNames
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join funs sample_n sample_frac intersect
#' @importFrom stats var median quantile
#' @importFrom purrr map map2 flatten flatten_df keep discard
#' @importFrom radiator tidy_genomic_data radiator_imputations_module write_genind snp_ld keep_common_markers radiator_maf_module detect_genomic_format
#' @importFrom tibble as_data_frame data_frame
#' @importFrom tidyr spread gather separate

#' @examples
#' \dontrun{
#' # with adegenet DAPC for the assignment and sampling.method = "random":
#' assignment.treefrog <- assignment_mixture(
#' data = "batch_1.vcf",
#' strata = "strata.tsv",
#' mixture = "mixture.treefrog.tsv",
#' assignment.analysis = "adegenet",
#' whitelist.markers = "whitelist.vcf.txt",
#' snp.ld = NULL,
#' common.markers = TRUE,
#' marker.number = c(500, 5000, "all"),
#' sampling.method = "random",
#' blacklist.id = "blacklist.id.tsv",
#' subsample = 25,
#' iteration.subsample = 5
#' filename = "treefrog.txt",
#' keep.gsi.files = FALSE,
#' pop.levels = c("PAN", "COS", "mixture")
#' imputation.method = NULL,
#' parallel.core = 12
#' )
#' # with gsi_sim for the mixture assignment and sampling.method = "ranked"
#' # Here I also want to impute the genotypes of the data (baseline and mixture) 
#' # using random forest:
#' assignment.tuna <- assignment_mixture(
#' data = "data.frame.genotypes.tuna.tsv",
#' mixture = "cohort.tuna.tsv",
#' assignment.analysis = "gsi_sim",
#' common.markers = TRUE,
#' marker.number = c(100, 200, 300),
#' sampling.method = "ranked",
#' subsample = 25,
#' iteration.subsample = 5
#' filename = "tuna.txt",
#' keep.gsi.files = FALSE,
#' pop.levels = c("BAJ", "IND", "mixture"),
#' imputation.method = "rf", 
#' impute.mixture = TRUE, 
#' hierarchical.levels = "populations", 
#' verbose = FALSE,
#' parallel.core = 12
#' )
#' 
#' Since the 'folder' argument is missing, it will be created automatically
#' inside your working directory.
#' 
#' use $ to access the data frames in the list created.
#' }


#' @references Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008)
#' An improved method for predicting the accuracy of genetic stock identification.
#' Canadian Journal of Fisheries and Aquatic Sciences 65, 7:1475-1486.
#' @references Anderson, E. C. (2010) Assessing the power of informative subsets of
#' loci for population assignment: standard methods are upwardly biased.
#' Molecular ecology resources 10, 4:701-710.
#' @references Catchen JM, Amores A, Hohenlohe PA et al. (2011)
#' Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences.
#' G3, 1, 171-182.
#' @references Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013)
#' Stacks: an analysis tool set for population genomics.
#' Molecular Ecology, 22, 3124-3140.
#' @references Weir BS, Cockerham CC (1984) Estimating F-Statistics for the
#' Analysis of Population Structure. Evolution, 38, 1358–1370.
#' @references Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival,
#' Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841--860.
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
#' Bender D, et al. 
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. 
#' American Journal of Human Genetics. 2007; 81: 559–575. doi:10.1086/519795

#' @seealso \code{gsi_sim} development page is available here: \url{https://github.com/eriqande/gsi_sim}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}


assignment_mixture <- function(
  data,
  strata = NULL,
  mixture = NULL,
  assignment.analysis,
  sampling.method,
  iteration.method = 10,
  subsample = NULL,
  iteration.subsample = 1,
  marker.number = "all",
  blacklist.id = NULL,
  blacklist.genotype = NULL,
  whitelist.markers = NULL,
  monomorphic.out = TRUE,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  max.marker = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  imputation.method = NULL,
  impute.mixture = FALSE,
  hierarchical.levels = "populations",
  verbose = FALSE,
  folder = NULL,
  filename = "assignment_data.txt",
  keep.gsi.files = FALSE,
  random.seed = NULL,
  parallel.core = parallel::detectCores() - 1
) {
  cat("#######################################################################\n")
  cat("#################### assigner::assignment_mixture #####################\n")
  cat("#######################################################################\n")
  timing <- proc.time()
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file missing")
  if (missing(assignment.analysis)) stop("assignment.analysis argument missing")
  if (assignment.analysis == "gsi_sim" & !gsi_sim_exists()) {
    stop("Can't find the gsi_sim executable where it was expected at ", gsi_sim_binary_path(), ".  
          If you have internet access, you can install it
          from within R by invoking the function \"install_gsi_sim(fromSource = TRUE)\"")
  }
  
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
 
  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...
  # pop.levels, pop.labels, pop.select, strata, etc
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stringi::stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  if (!is.null(pop.labels)) {
    if (length(pop.labels) != length(pop.levels)) stop("pop.labels and pop.levels must have the same length")
    pop.labels <- stringi::stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  if (!is.null(pop.select)) {
    pop.select <- stringi::stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  # store function call
  function.call <- match.call()
  
  # File type detection ********************************************************
  data.type <- radiator::detect_genomic_format(data)
  
  
  if (data.type == "haplo.file") {
    message("With stacks haplotype file the maf.approach is automatically set to: haplotype")
    maf.approach <- "SNP"
    # confusing, but because the haplotpe file doesn't have snp info, only locus info
    # it's treated as markers/snp info and filtered the same way as the approach by SNP.
    # but it's really by haplotype
  }
  
  if (maf.approach == "haplotype") {
    if (data.type != "vcf.file" | data.type != "haplo.file") {
      stop("The haplotype approach during MAF filtering is for VCF and
           stacks haplotypes file, only. Use the snp approach for the other file types")
    }
  }
  
  # Strata argument required for VCF and haplotypes files-----------------------
  if (data.type == "haplo.file" | data.type == "vcf.file") {
    if (is.null(strata)) stop("strata argument is required")
  }
  
  # Base filename
  base.filename <- filename # filename will change from time to time in the function
  
  # Create a folder based on filename to save the output files -----------------
  if (is.null(folder)) {
    # Get date and time to have unique filenaming
    file.date <- stringi::stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stringi::stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stringi::stri_sub(file.date, from = 1, to = 13)
    
    if (is.null(imputation.method)) {
      message("Map-imputation: no")
      directory <- stringi::stri_join(getwd(),"/", "assignment_mixture_analysis", "_no_imputations_", file.date, "/", sep = "")
      dir.create(file.path(directory))
    } else {
      message("Map-imputation: yes")
      directory <- stringi::stri_join(getwd(),"/","assignment_mixture_analysis", "_imputations_", imputation.method,"_", hierarchical.levels, "_", file.date, "/", sep = "")
      dir.create(file.path(directory))
    }
    message(stringi::stri_join("Folder: ", directory))
    file.date <- NULL #unused object
  } else {
    directory <- stringi::stri_join(getwd(), "/", folder, "/", sep = "")
    dir.create(file.path(directory))
    message(stringi::stri_join("Folder: ", directory))
  }
  
  # Import input ---------------------------------------------------------------
  input <- radiator::tidy_genomic_data(
    data = data, 
    vcf.metadata = FALSE,
    blacklist.id = blacklist.id, 
    blacklist.genotype = blacklist.genotype, 
    whitelist.markers = whitelist.markers, 
    monomorphic.out = monomorphic.out, 
    max.marker = max.marker,
    snp.ld = snp.ld,
    common.markers = FALSE, 
    maf.thresholds = maf.thresholds,
    maf.pop.num.threshold = maf.pop.num.threshold,
    maf.approach = maf.approach,
    maf.operator = maf.operator,
    strata = strata, 
    pop.levels = pop.levels, 
    pop.labels = pop.labels, 
    pop.select = pop.select,
    filename = NULL
  )
  
  # check
  # levels(input$POP_ID)
  # dplyr::n_distinct(input$MARKERS)
  
  # change "unknown" to "mixture" for simplicity of pop_id below ---------------
  input$POP_ID <- stringi::stri_replace_all_fixed(
    str = input$POP_ID, pattern = "unknnown", replacement = "mixture",
    vectorize_all = FALSE
  )
  # need to remove space in POP_ID name to work in gsi_sim
  input$POP_ID <- stringi::stri_replace_all_fixed(
    input$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE
  )
  
  # mixture data ---------------------------------------------------------------
  if ("mixture" %in% unique(input$POP_ID)) {
    mixture.df <- input %>% 
      dplyr::filter(POP_ID == "mixture") %>% 
      dplyr::distinct(INDIVIDUALS)
    
    message(stringi::stri_join("Number of samples identified as unknown/mixture: ", length(mixture.df$INDIVIDUALS)))
  }
  
  if (!is.null(mixture)) {
    mixture.df <- readr::read_tsv(file = mixture, col_names = TRUE, col_types = "c")
    input <- dplyr::mutate(.data = input, POP_ID = dplyr::if_else(INDIVIDUALS %in% mixture.df$INDIVIDUALS, "mixture", as.character(POP_ID)))
    message(stringi::stri_join("Number of samples identified as unknown/mixture: ", length(mixture.df$INDIVIDUALS)))
  }
  
  # keep a new pop.levels and pop.labels  --------------------------------------
  if (is.null(pop.levels)) {
    pop.levels <- pop.labels <- levels(factor(input$POP_ID))
  }
  
  # subsampling data -----------------------------------------------------------
  # create the subsampling list
  ind.pop.df <- dplyr::distinct(.data = input, POP_ID, INDIVIDUALS) %>% dplyr::arrange(POP_ID, INDIVIDUALS)
  
  if (!is.null(subsample)) {
    min.pop.n <- min(dplyr::count(x = dplyr::filter(ind.pop.df, POP_ID != "mixture"), POP_ID, sort = TRUE)$n)
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
      } else {
        stop("Wrong subsample value")
      }
    }
  }
  
  subsample.list <- purrr::map(
    .x = 1:iteration.subsample,
    .f = subsampling_baseline,#internal function below
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
  } # End subsampling
  
  # unused objects
  ind.pop.df <- subsampling.individuals <- NULL
  
  # assignment analysis --------------------------------------------------------
  res <- purrr::map(
    .x = subsample.list,
    .f = assignment_function,
    input = input,
    subsample = subsample,
    assignment.analysis = assignment.analysis,
    mixture.df = mixture.df,
    snp.ld = snp.ld,
    common.markers = common.markers,
    maf.thresholds = maf.thresholds,
    maf.pop.num.threshold = maf.pop.num.threshold,
    maf.approach = maf.approach,
    maf.operator = maf.operator,
    marker.number = marker.number,
    pop.levels = pop.levels,
    pop.labels = pop.labels,
    sampling.method = sampling.method,
    iteration.method = iteration.method,
    filename = filename,
    directory = directory,
    keep.gsi.files = keep.gsi.files,
    imputation.method = imputation.method,
    impute.mixture = impute.mixture,
    hierarchical.levels = hierarchical.levels,
    verbose = verbose,
    parallel.core = parallel.core,
    manage.all = manage.all,
    random.seed = random.seed,
    base.filename = base.filename
  )
  res <- dplyr::bind_rows(res)
  
  if (!is.null(subsample)) {
    readr::write_tsv(x = res, path = paste0(directory, "assignment.mixture.results.tsv"), col_names = TRUE, append = FALSE)
  }
  
  # Summary of the subsampling iterations
  if (sampling.method == "random") {
    if (assignment.analysis == "gsi_sim") {
      assignment.mixture.summary.subsample <- res %>% 
        dplyr::group_by(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, INFERRED, SUBSAMPLE) %>%
        dplyr::summarise(
          MEAN_MARKERS_COMMON = mean(MARKERS_COMMON, na.rm = TRUE),
          NUMBER_ITERATIONS = length(ITERATIONS),
          MEAN_ITERATIONS = round((NUMBER_ITERATIONS/iteration.method)*100, 2),
          MEAN = round(mean(SCORE), 2),
          SE = round(sqrt(stats::var(SCORE)/length(SCORE)), 2),
          MIN = round(min(SCORE), 2),
          MAX = round(max(SCORE), 2),
          MEDIAN = round(stats::median(SCORE), 2),
          QUANTILE25 = round(stats::quantile(SCORE, 0.25), 2),
          QUANTILE75 = round(stats::quantile(SCORE, 0.75), 2)
        ) %>% 
        dplyr::mutate(
          TOTAL_ITERATIONS = rep(iteration.method, n())
        ) %>% 
        dplyr::select(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MEAN_MARKERS_COMMON, MISSING_DATA, SUBSAMPLE, INFERRED, NUMBER_ITERATIONS, TOTAL_ITERATIONS, MEAN_ITERATIONS, MEAN, SE, MIN, MAX, MEDIAN, QUANTILE25, QUANTILE75) %>% 
        dplyr::arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE) %>% 
        dplyr::ungroup(.)
    }
    if (assignment.analysis == "adegenet") {
      assignment.mixture.summary.subsample <- res %>% 
        # dplyr::select(-X1, -X2) %>% 
        dplyr::ungroup(.) %>%
        dplyr::mutate(CURRENT = factor(CURRENT)) %>% 
        dplyr::group_by(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE) %>%
        dplyr::summarise(
          NUMBER_ITERATIONS = length(ITERATIONS),
          MEAN_ITERATIONS = round((NUMBER_ITERATIONS/iteration.method)*100, 2)
        ) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE, CURRENT, INFERRED) %>% 
        dplyr::group_by(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, MISSING_DATA) %>%
        dplyr::summarise(
          MEAN_SUBSAMPLE = round(mean(MEAN_ITERATIONS), 2),
          SE = round(sqrt(stats::var(MEAN_ITERATIONS)/length(MEAN_ITERATIONS)), 2),
          MIN = round(min(MEAN_ITERATIONS), 2),
          MAX = round(max(MEAN_ITERATIONS), 2),
          MEDIAN = round(stats::median(MEAN_ITERATIONS), 2),
          QUANTILE25 = round(stats::quantile(MEAN_ITERATIONS, 0.25), 2),
          QUANTILE75 = round(stats::quantile(MEAN_ITERATIONS, 0.75), 2)
        ) %>%
        dplyr::ungroup(.) %>%
        dplyr::arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, CURRENT, INFERRED)
    }
  } # end random
  
  if (sampling.method == "ranked") {
    if (assignment.analysis == "gsi_sim") {
      assignment.mixture.summary.subsample <- res %>% 
        dplyr::group_by(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA, INFERRED) %>%
        dplyr::summarise(
          MEAN_MARKERS_COMMON = mean(MARKERS_COMMON, na.rm = TRUE),
          NUMBER_SUBSAMPLE = length(SUBSAMPLE),
          MEAN_SUBSAMPLE = round((NUMBER_SUBSAMPLE/iteration.subsample)*100, 2),
          MEAN = round(mean(SCORE), 2),
          SE = round(sqrt(stats::var(SCORE)/length(SCORE)), 2),
          MIN = round(min(SCORE), 2),
          MAX = round(max(SCORE), 2),
          MEDIAN = round(stats::median(SCORE), 2),
          QUANTILE25 = round(stats::quantile(SCORE, 0.25), 2),
          QUANTILE75 = round(stats::quantile(SCORE, 0.75), 2)
        ) %>% 
        dplyr::mutate(TOTAL_SUBSAMPLE = rep(iteration.subsample, n())) %>% 
        dplyr::select(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MEAN_MARKERS_COMMON, METHOD, MISSING_DATA, INFERRED, NUMBER_SUBSAMPLE, TOTAL_SUBSAMPLE, MEAN_SUBSAMPLE, MEAN, SE, MIN, MAX, MEDIAN, QUANTILE25, QUANTILE75) %>% 
        dplyr::arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA) %>% 
        dplyr::ungroup(.)
    }
    if (assignment.analysis == "adegenet") {
      assignment.mixture.summary.subsample <- res %>%
        dplyr::select(-ITERATIONS) %>%
        dplyr::ungroup(.) %>%
        dplyr::mutate(CURRENT = factor(CURRENT)) %>%
        dplyr::group_by(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
        dplyr::summarise(
          NUMBER_SUBSAMPLE = length(SUBSAMPLE),
          MEAN_SUBSAMPLE = round((NUMBER_SUBSAMPLE/iteration.subsample)*100, 2)) %>% 
        dplyr::ungroup(.) %>% 
        dplyr::mutate(TOTAL_SUBSAMPLE = rep(iteration.subsample, n())) %>% 
        dplyr::select(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA, NUMBER_SUBSAMPLE, TOTAL_SUBSAMPLE, MEAN_SUBSAMPLE) %>% 
        dplyr::arrange(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA)
    }
  } # end ranked
  
  
  # assignment summary results -------------------------------------------------
  if (is.null(imputation.method)) {
    filename.assignment.sum <- stringi::stri_join("assignment.mixture.summary.results", "no.imputation", sampling.method, "tsv", sep = ".")
  } else {# with imputations
    filename.assignment.sum <- stringi::stri_join("assignment.mixture.summary.results", "imputed", sampling.method, "tsv", sep = ".")
  }
  readr::write_tsv(x = assignment.mixture.summary.subsample, path = paste0(directory,filename.assignment.sum), col_names = TRUE, append = FALSE)
  
  # results
  res.list <- list(call = function.call, assignment = res, assignment.mixture.summary.results = assignment.mixture.summary.subsample)
  
  timing <- proc.time() - timing
  message(stringi::stri_join("Computation time: ", round(timing[[3]]), " sec"))
  cat("############################## completed ##############################\n")
  return(res.list)
} # End assignment_mixture

# Internal Nested Functions -----------------------------------------------------------

# subsampling_baseline --------------------------------------------------------------
#' @title subsampling baseline data
#' @description subsampling baseline data
#' @rdname subsampling_baseline
#' @export
#' @keywords internal
subsampling_baseline <- function(
  iteration.subsample = 1,
  ind.pop.df = NULL,
  subsample = NULL,
  random.seed = NULL,
  ...
) {
  if (is.null(subsample)) {
    subsample.select <- ind.pop.df %>% 
      dplyr::mutate(SUBSAMPLE = rep(iteration.subsample, n()))
  } else {
    # separate all the mixture samples
    mixture.select <- ind.pop.df %>% dplyr::filter(POP_ID == "mixture")
    
    # Set seed for sampling reproducibility
    if (is.null(random.seed)) {
      random.seed <- sample(x = 1:1000000, size = 1)
      set.seed(random.seed)
    } else {
      set.seed(random.seed)
    }
    
    # subsample the baseline
    if (subsample > 1) { # integer
      subsample.select <- ind.pop.df %>%
        dplyr::filter(POP_ID != "mixture") %>% 
        dplyr::group_by(POP_ID) %>%
        dplyr::sample_n(tbl = ., size = subsample, replace = FALSE)# sampling individuals for each pop
    }
    
    # Join baseline and mixture back in 1 dataset
    subsample.select <- dplyr::bind_rows(subsample.select, mixture.select) %>% 
      dplyr::mutate(
        SUBSAMPLE = rep(iteration.subsample, n()),
        RANDOM_SEED_NUMBER = rep(random.seed, n())
      ) %>%
      dplyr::arrange(POP_ID, INDIVIDUALS) %>% 
      dplyr::ungroup(.)
  }
  return(subsample.select)
} # End subsampling function

# mixture_baseline_intersect----------------------------------------------------
#' @title mixture_baseline_intersect
#' @description mixture_baseline_intersect
#' @rdname mixture_baseline_intersect
#' @export
#' @keywords internal
mixture_baseline_intersect <- function(
  mixture.id,
  inferred.baseline,
  input.mixture = NULL,
  input.baseline = NULL
) {
  mixture.id.data <- dplyr::ungroup(input.mixture) %>%
    dplyr::filter(GT != "000000") %>% 
    dplyr::filter(INDIVIDUALS == mixture.id) %>% 
    dplyr::distinct(MARKERS)
  
  inferred.baseline.data <- dplyr::ungroup(input.baseline) %>% 
    dplyr::filter(GT != "000000") %>% 
    dplyr::filter(POP_ID == inferred.baseline) %>% 
    dplyr::distinct(MARKERS)
  
  marker.in.common <- nrow(dplyr::intersect(x = mixture.id.data, y = inferred.baseline.data))
  # diff <- dplyr::setdiff(x = mixture.id.data, y = inferred.baseline.data)
  # marker.in.common.df <- tibble::data_frame(INDIVIDUALS = mixture.id, INFERRED = inferred.baseline, MARKERS_COMMON = marker.in.common)
  return(marker.in.common)
}# End mixture_baseline_intersect

# assignment_gsi_sim------------------------------------------------------------
#' @title assignment with gsi_sim
#' @description assignment with gsi_sim
#' @rdname assignment_gsi_sim
#' @export
#' @keywords internal

assignment_gsi_sim <- function(
  data = NULL,
  strata.df = NULL,
  select.markers = NULL,
  markers.names = NULL,
  missing.data = NULL,
  i = NULL,
  m = NULL,
  filename = NULL,
  directory.subsample = NULL,
  keep.gsi.files = FALSE,
  pop.labels = NULL,
  sampling.method = "random",
  ...
) {
  # data <- input #test
  # missing.data <- "no.imputation" #test
  data.select <- suppressWarnings(
    dplyr::semi_join(data, select.markers, by = "MARKERS") %>%
      dplyr::arrange(POP_ID, INDIVIDUALS, MARKERS)
  )
  
  input.baseline <- dplyr::filter(.data = data.select, POP_ID != "mixture")
  input.mixture <- dplyr::filter(.data = data.select, POP_ID == "mixture")
  
  # Baseline filename
  baseline.filename <- stringi::stri_replace_all_fixed(
    filename, 
    pattern = ".txt",
    replacement = "_baseline.txt", 
    vectorize_all = FALSE
  )
  baseline.input.gsi <- stringi::stri_join(directory.subsample, baseline.filename)
  
  # save input file to directory
  assigner::write_gsi_sim(
    data = input.baseline, 
    # markers.names = markers.names, 
    filename = baseline.input.gsi
  )
  
  # Mixture filename
  mixture.filename <- stringi::stri_replace_all_fixed(
    filename, 
    pattern = ".txt",
    replacement = "_mixture.txt",
    vectorize_all = FALSE
  )
  mixture.input.gsi <- stringi::stri_join(directory.subsample, mixture.filename)
  
  # save input file to directory
  assigner::write_gsi_sim(
    data = input.mixture, 
    markers.names = markers.names, 
    # directory.subsample = directory.subsample, 
    # i = i, 
    # m = m,
    filename = mixture.input.gsi
  )
  
  # Run gsi_sim ------------------------------------------------------------
  output.gsi <- stringi::stri_replace_all_fixed(mixture.input.gsi, pattern = "txt", replacement = "output.txt")
  setwd(directory.subsample)
  system(paste("gsi_sim -b ", baseline.input.gsi, "-t ", mixture.input.gsi, " > ", output.gsi))
  
  # Option remove the input file from directory to save space
  if (!keep.gsi.files) {
    file.remove(baseline.input.gsi)
    file.remove(mixture.input.gsi)
  }
  
  # Get Assignment results -------------------------------------------------
  assignment <- suppressWarnings(suppressMessages(
    readr::read_delim(output.gsi, col_names = "ID", delim = "\t") %>%
      tidyr::separate(ID, c("KEEPER", "ASSIGN"), sep = ":/", extra = "warn") %>%
      dplyr::filter(KEEPER == "GMA_FULL_EM_INDIVS_CSV") %>%
      tidyr::separate(ASSIGN, c("INDIVIDUALS", "ASSIGN"), sep = ";", extra = "merge") %>%
      tidyr::separate(ASSIGN, c("INFERRED", "OTHERS"), sep = ";", convert = TRUE, numerals = "no.loss", extra = "merge") %>%
      tidyr::separate(OTHERS, c("SCORE", "OTHERS"), sep = ";;", convert = TRUE, numerals = "no.loss", extra = "merge") %>%
      tidyr::separate(OTHERS, c("SECOND_BEST_POP", "OTHERS"), sep = ";", convert = TRUE, numerals = "no.loss", extra = "merge") %>%
      tidyr::separate(OTHERS, c("SECOND_BEST_SCORE", "OTHERS"), sep = ";;", convert = TRUE, numerals = "no.loss")
  ))
  
  assignment <- suppressWarnings(
    assignment %>%
      dplyr::mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>% 
      dplyr::left_join(strata.df, by = "INDIVIDUALS") %>%
      dplyr::rename(CURRENT = POP_ID) %>% 
      dplyr::mutate(
        ANALYSIS = rep("mixture", n()),
        CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
        # CURRENT = factor(CURRENT, levels = pop.levels, labels = pop.labels, ordered = TRUE),
        CURRENT = droplevels(CURRENT),
        INFERRED = factor(INFERRED, levels = unique(pop.labels), ordered = TRUE),
        INFERRED = droplevels(INFERRED),
        SECOND_BEST_POP = factor(SECOND_BEST_POP, levels = unique(pop.labels), ordered = TRUE),
        SECOND_BEST_POP = droplevels(SECOND_BEST_POP),
        SCORE = round(SCORE, 2),
        SECOND_BEST_SCORE = round(SECOND_BEST_SCORE, 2),
        MARKER_NUMBER = as.numeric(rep(m, n())),
        METHOD = rep(sampling.method, n()),
        MISSING_DATA = rep(missing.data, n())
      ) %>%
      dplyr::select(INDIVIDUALS, ANALYSIS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
      dplyr::arrange(CURRENT)
  )
  
  if (sampling.method == "random") {
    assignment <- assignment %>% 
      dplyr::mutate(ITERATIONS = rep(i, n())) %>% 
      dplyr::select(INDIVIDUALS, ANALYSIS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA, ITERATIONS) %>%
      dplyr::arrange(CURRENT)
  }
  
  # Markers between mixture and inferred pop -----------------------------------
  
  mixture.baseline.itersect <- dplyr::ungroup(assignment) %>% 
    dplyr::select(INDIVIDUALS, INFERRED) %>% 
    dplyr::mutate(
      MARKERS_COMMON = unlist(
        purrr::map2(
          .x = INDIVIDUALS, .y = INFERRED, .f = mixture_baseline_intersect,
          input.mixture = input.mixture, input.baseline = input.baseline
        )
      )
    )
  
  assignment <- assignment %>% 
    dplyr::full_join(mixture.baseline.itersect, by = c("INDIVIDUALS", "INFERRED")) %>% 
    dplyr::ungroup(.)
  
  # Option remove the output file from directory to save space
  if (!keep.gsi.files) file.remove(output.gsi)
  
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
  i = NULL,
  m = NULL,
  holdout = NULL,
  sampling.method = "random",
  subsample.id = NULL,
  ...
) {
  # data <- genind.object #test
  # missing.data <- "no.imputation" #test
  data.select <- data[loc = select.markers$MARKERS]
  
  # Run adegenet *********************************************************
  pop.data <- data.select@pop
  pop.data <- droplevels(pop.data)
  
  # # Alpha-Score DAPC
  # # When all the individuals are accounted for in the model construction
  # dapc.best.optim.a.score <- optim.a.score(dapc(data.select, n.da = length(levels(pop.data)), n.pca = round((length(indNames(data.select))/3)-1, 0)), pop = pop.data, plot = FALSE)$best
  # message(stringi::stri_join("a-score optimisation for iteration:", i, sep = " ")) # message not working in parallel...
  # 
  # # DAPC with all the data
  # dapc.all <- dapc(data.select, n.da = length(levels(pop.data)), n.pca = dapc.best.optim.a.score, pop = pop.data)
  # message(stringi::stri_join("DAPC iteration:", i, sep = " "))
  # message(stringi::stri_join("DAPC marker group:", m, sep = " "))
  
  # Alpha-Score DAPC training data
  training.data <- data.select[data.select@strata$POP_ID != "mixture"]
  # training.data <- data.select[!indNames(data.select) %in% holdout$INDIVIDUALS] # training dataset
  # indNames(training.data)
  # training.data@strata
  pop.training <- training.data@pop
  pop.training <- droplevels(pop.training)
  
  dapc.best.optim.a.score <- adegenet::optim.a.score(
    adegenet::dapc(
      training.data, 
      n.da = length(levels(pop.training)),
      n.pca = round(((length(adegenet::indNames(training.data))/3) - 1), 0)
    ),
    pop = pop.training, plot = FALSE)$best
  message(stringi::stri_join("a-score optimisation for iteration:", i, sep = " "))
  
  dapc.training <- adegenet::dapc(
    training.data,
    n.da = length(levels(pop.training)), 
    n.pca = dapc.best.optim.a.score, 
    pop = pop.training)
  message(stringi::stri_join("DAPC of training data set for iteration:", i, sep = " "))
  
  # DAPC holdout individuals
  holdout.data <- data.select[data.select@strata$POP_ID == "mixture"]
  # indNames(holdout.data)
  # holdout.data@strata
  # holdout.data <- data.select[indNames(data.select) %in% holdout$INDIVIDUALS] # holdout dataset
  pop.holdout <- holdout.data@pop
  pop.holdout <- droplevels(pop.holdout)
  
  assignment.levels <- levels(pop.data)
  rev.assignment.levels <- rev(assignment.levels)
  
  dapc.predict.holdout <- adegenet::predict.dapc(dapc.training, newdata = holdout.data)
  message(stringi::stri_join("Assigning holdout data for iteration:", i, sep = " "))
  
  
  # Get Assignment results -----------------------------------------------
  if (sampling.method == "ranked") {
    i <- "not available with sampling.method = ranked"
  }
  
  assignment <- data.frame(INDIVIDUALS = adegenet::indNames(holdout.data), POP_ID = pop.holdout, ASSIGN = dapc.predict.holdout$assign, dapc.predict.holdout$posterior) %>% 
    dplyr::rename(CURRENT = POP_ID, INFERRED = ASSIGN) %>%
    dplyr::mutate(
      ANALYSIS = rep("mixture", n()),
      MARKER_NUMBER = as.numeric(rep(m, n())),
      METHOD = rep(sampling.method, n()),
      MISSING_DATA = rep(missing.data, n()),
      SUBSAMPLE = rep(subsample.id, n()),
      CURRENT = factor(CURRENT, levels = rev.assignment.levels, ordered = TRUE),
      INFERRED = factor(INFERRED, levels = assignment.levels, ordered = TRUE),
      ITERATIONS = rep(i, n())
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
  base.filename = NULL,
  input = NULL,
  genind.object = NULL,
  strata.df = NULL,
  imputation.method = NULL,
  hierarchical.levels = "populations",
  directory.subsample = NULL,
  keep.gsi.files = FALSE,
  sampling.method = "random",
  pop.labels = NULL,
  subsample.id = NULL,
  ...
) {
  x <- tibble::as_data_frame(x)
  # x <- x[1] # test
  i <- as.numeric(unique(x$ITERATIONS))      # iteration
  m <- as.numeric(unique(x$MARKER_NUMBER))   # number of marker selected
  
  select.markers <- dplyr::ungroup(x) %>% dplyr::select(MARKERS) %>% 
    dplyr::arrange(MARKERS)
  
  # get the list of loci after filter
  markers.names <- unique(select.markers$MARKERS)
  
  # Assignment analysis
  
  if (is.null(imputation.method)) {
    filename.imp <- "no_imputation.txt"
    missing.data <- "no.imputation"
  } else {
    filename.imp <- "imputed.txt"
    
    if (imputation.method == "rf") {
      if (hierarchical.levels == "populations") {
        missing.data <- "imputed RF populations"
      } else {
        missing.data <- "imputed RF global"
      }
    } else {
      if (hierarchical.levels == "populations") {
        missing.data <- "imputed max populations"
      } else {
        missing.data <- "imputed max global"
      }
    }
  }
  
  # Modify filename
  filename <- stringi::stri_replace_all_fixed(
    base.filename,
    pattern = ".txt",
    replacement = stringi::stri_join(
      "", "iteration", i, "markers", m, 
      filename.imp, sep = "_"
    )
  )
  
  if (assignment.analysis == "gsi_sim") {
    assignment <- assignment_gsi_sim(
      data = input,
      strata.df = strata.df,
      select.markers = select.markers,
      markers.names = markers.names,
      missing.data = missing.data, 
      i = i, 
      m = m,
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
      markers.names = markers.names,
      missing.data = missing.data, 
      i = i, 
      m = m,
      sampling.method = sampling.method,
      subsample.id = subsample.id
    )
  }
  
  # unused objects
  x <- i <- m <- select.markers <- markers.names <- filename <- NULL
  filename.imp <- missing.data <- NULL
  
  return(assignment)
} #End assignment_random


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
  pop.levels = NULL,
  pop.labels = NULL,
  sampling.method = NULL,
  iteration.method = NULL,
  base.filename = NULL,
  keep.gsi.files = FALSE,
  imputation.method = NULL,
  hierarchical.levels = NULL,
  parallel.core = parallel::detectCores() - 1,
  directory.subsample = NULL,
  subsample.id = NULL,
  ...
) {
  message("Marker number: ", x)
  # x <- 200 # test
  # x <- marker.number
  x <- as.numeric(x)
  
  select.markers <- dplyr::filter(.data = fst.ranked, RANKING <= x) %>%
    dplyr::select(MARKERS)
  
  # get the list of markers after filter
  markers.names <- unique(select.markers$MARKERS)
  
  if (is.null(imputation.method)) {
    filename.imp <- "no_imputation.txt"
    missing.data <- "no.imputation"
  } else {
    filename.imp <- "imputed.txt"
    
    if (imputation.method == "rf") {
      if (hierarchical.levels == "populations") {
        missing.data <- "imputed RF populations"
      } else {
        missing.data <- "imputed RF global"
      }
    } else {
      if (hierarchical.levels == "populations") {
        missing.data <- "imputed max populations"
      } else {
        missing.data <- "imputed max global"
      }
    }
  }
  
  # Modify filename
  filename <- stringi::stri_replace_all_fixed(
    base.filename,
    pattern = ".txt",
    replacement = stringi::stri_join(
      "", "markers", x, 
      filename.imp, sep = "_"
    )
  )
  
  if (assignment.analysis == "gsi_sim") {
    assignment <- assignment_gsi_sim(
      data = input,
      strata.df = strata.df,
      select.markers = select.markers,
      markers.names = markers.names,
      missing.data = missing.data,
      i = NULL,
      m = x,
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
      markers.names = markers.names,
      missing.data = missing.data, 
      i = i, 
      m = x,
      sampling.method = sampling.method,
      subsample.id = subsample.id
    )
  }
  # unused objects
  x <- select.markers <- markers.names <- NULL
  filename.imp <- missing.data <- filename <- NULL
  
  return(assignment)
}#End assignment_marker_loop

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
  mixture.df = NULL,
  snp.ld = NULL,
  common.markers = TRUE,
  maf.thresholds = NULL,
  maf.pop.num.threshold = 1,
  maf.approach = "SNP",
  maf.operator = "OR",
  marker.number = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  sampling.method = "random",
  iteration.method = 10,
  filename = "assignment_data.txt",
  directory = NULL,
  keep.gsi.files = FALSE,
  imputation.method = NULL,
  impute.mixture = FALSE,
  hierarchical.levels = "populations",
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1,
  manage.all = NULL, random.seed = NULL, base.filename = NULL,
  ...
) {
  # x <- subsample.list[[1]] # test
  subsampling.individuals <- x
  subsample.id <- unique(subsampling.individuals$SUBSAMPLE)
  
  if (!is.null(subsample)) {
    message(paste("Analyzing subsample: ", subsample.id))
  }
  
  # Updating directories for subsampling
  if (is.null(subsample)) {
    directory.subsample <- directory
  } else {
    directory.subsample <- paste0(directory, "subsample_", subsample.id, "/")
    dir.create(file.path(directory.subsample))
  }
  
  # Keep only the subsample
  input <- dplyr::semi_join(input, subsampling.individuals, by = c("POP_ID", "INDIVIDUALS"))
  
  # unused object
  subsampling.individuals <- x <- NULL
  
  # LD control... keep only 1 SNP per haplotypes/reads (optional) ------------
  if (!is.null(snp.ld)) {
    input <- radiator::snp_ld(data = input, snp.ld = snp.ld)
  } # End of snp.ld control
  
  
  # Markers in common between all populations (optional) ---------------------
  if (common.markers) { # keep only markers present in all pop
    input <- radiator::keep_common_markers(data = input)$input
  } # End common markers
  
  # Minor Allele Frequency filter --------------------------------------------
  if (!is.null(maf.thresholds)) {
    message("Note that MAF filtering is conducted on baseline data only")
    # MAF
    mixture.data <- dplyr::filter(.data = input, POP_ID == "mixture")
    baseline.data <- dplyr::filter(.data = input, POP_ID != "mixture")
    
    # maf.thresholds <- c(0.05, 0.05) # test
    baseline.data <- radiator::radiator_maf_module(
      data = baseline.data,
      maf.thresholds = maf.thresholds,
      maf.pop.num.threshold = maf.pop.num.threshold,
      maf.approach = maf.approach,
      maf.operator = maf.operator
    )
    
    input <- dplyr::bind_rows(mixture.data, baseline.data$input) %>% 
      dplyr::arrange(POP_ID, INDIVIDUALS)
    
    mixture.data <- baseline.data <- NULL # unused objects
  } # End of MAF filters
  
  # Keep a strata df --------------------------------------------------------
  strata.df <- dplyr::distinct(.data = input, INDIVIDUALS, POP_ID)
  
  # Adegenet no imputations --------------------------------------------------
  if (assignment.analysis == "adegenet") {
    message("Creating genind object")
    genind.object <- radiator::write_genind(data = input)
  } else {
    genind.object <- NULL
  }
  
  input.baseline <- dplyr::filter(.data = input, POP_ID != "mixture")
  input.mixture <- dplyr::filter(.data = input, POP_ID == "mixture")
  
  # Imputations --------------------------------------------------------------
  if (!is.null(imputation.method)) {
    message("Preparing the data for imputations")
    
    # imputation for the mixture samples, if selected, is always conducted globally. 
    message("Imputations of baseline samples ...")
    input.baseline.imp <- radiator::radiator_imputations_module(
      data = input.baseline, 
      imputation.method = imputation.method, 
      hierarchical.levels = hierarchical.levels, 
      verbose = verbose, 
      parallel.core = parallel.core, 
      filename = NULL
    )
    
    # combine the mixture (no imputation) + the imputed baseline
    input.imp <- suppressWarnings(
      dplyr::bind_rows(input.baseline.imp, input.mixture)
    )
    
    if (impute.mixture) {
      # impute globally the mixture samples
      message("Imputations computed globally for mixture samples:")
      input.imp <- radiator::radiator_imputations_module(
        data = input.imp, 
        imputation.method = imputation.method, 
        hierarchical.levels = "global", 
        verbose = verbose, 
        parallel.core = parallel.core, 
        filename = NULL
      )
    } # End impute.mixture
    
    # input.baseline.imp <- dplyr::filter(.data = input.imp, POP_ID != "mixture")
    # input.mixture.imp <- dplyr::filter(.data = input.imp, POP_ID == "mixture")
    
    # test <- input.imp %>% dplyr::filter(GT == "000000")
    # prep. adegenet
    if (assignment.analysis == "adegenet") {
      genind.object.imp <- radiator::write_genind(data = input.imp)
    } else {
      genind.object.imp <- NULL
    } # end adegenet
  } # end imputations
  
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
  
  # Random method --------------------------------------------------------------
  if (sampling.method == "random") {
    message("Conducting Assignment analysis with markers selected randomly")
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
    readr::write_tsv(
      x = tibble::as_data_frame(dplyr::bind_rows(marker.random.list)),
      path = paste0(directory.subsample, "markers.random.tsv"),
      col_names = TRUE, append = FALSE
    )
    
    message("Starting parallel computations for the assignment analysis: random markers")
    message("For progress: monitor activity in the folder...")
    
    assignment.res <- NULL
    assignment.res <- suppressWarnings(
      .assigner_parallel(
        X = marker.random.list, 
        FUN = assignment_random, 
        mc.preschedule = FALSE, 
        mc.silent = FALSE,
        mc.cleanup = TRUE,
        mc.cores = parallel.core,
        assignment.analysis = assignment.analysis,
        base.filename = base.filename,
        input = input,
        genind.object = genind.object,
        strata.df = strata.df,
        imputation.method = NULL,
        hierarchical.levels = NULL,
        directory.subsample = directory.subsample,
        keep.gsi.files = keep.gsi.files,
        sampling.method = sampling.method,
        pop.labels = pop.labels,
        subsample.id = subsample.id
      ) %>% 
        dplyr::bind_rows(assignment.res)
    )
    
    # with imputations
    if (!is.null(imputation.method)) {
      assignment.res.imp <- NULL
      assignment.res.imp <- suppressWarnings(
        .assigner_parallel(
          X = marker.random.list, 
          FUN = assignment_random, 
          mc.preschedule = FALSE, 
          mc.silent = FALSE,
          mc.cleanup = TRUE,
          mc.cores = parallel.core,
          assignment.analysis = assignment.analysis,
          base.filename = base.filename,
          input = input.imp,
          genind.object = genind.object.imp,
          strata.df = strata.df,
          imputation.method = imputation.method,
          hierarchical.levels = hierarchical.levels,
          directory.subsample = directory.subsample,
          keep.gsi.files = keep.gsi.files,
          sampling.method = sampling.method,
          pop.labels = pop.labels,
          subsample.id = subsample.id
        )
      )
      assignment.res <- dplyr::bind_rows(assignment.res.imp) %>%
        dplyr::bind_rows(assignment.res)
    }
    
    # Compiling the results
    message("Compiling results")
    
    assignment.res <- suppressWarnings(
      dplyr::mutate(.data = assignment.res, SUBSAMPLE = rep(subsample.id, n())) %>% 
        dplyr::arrange(INDIVIDUALS, MARKER_NUMBER, MISSING_DATA, ITERATIONS)
    )
    
    # Write to the directory assignment results
    if (is.null(imputation.method)) {
      filename.assignment.res <- stringi::stri_join("assignment.mixture", "no.imputation", sampling.method, "tsv", sep = ".")
    } else {# with imputations
      filename.assignment.res <- stringi::stri_join("assignment.mixture", "imputed", sampling.method, "tsv", sep = ".")
    }
    readr::write_tsv(x = assignment.res, path = paste0(directory.subsample, filename.assignment.res), col_names = TRUE, append = FALSE)
    
    if (assignment.analysis == "gsi_sim") {
      assignment.mixture.summary.stats <- assignment.res %>% 
        dplyr::group_by(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, INFERRED, SUBSAMPLE) %>%
        dplyr::summarise(
          MEAN_MARKERS_COMMON = mean(MARKERS_COMMON, na.rm = TRUE),
          NUMBER_ITERATIONS = length(ITERATIONS),
          MEAN_ITERATIONS = round((NUMBER_ITERATIONS/iteration.method)*100, 2),
          MEAN = round(mean(SCORE), 2),
          SE = round(sqrt(stats::var(SCORE)/length(SCORE)), 2),
          MIN = round(min(SCORE), 2),
          MAX = round(max(SCORE), 2),
          MEDIAN = round(stats::median(SCORE), 2),
          QUANTILE25 = round(stats::quantile(SCORE, 0.25), 2),
          QUANTILE75 = round(stats::quantile(SCORE, 0.75), 2)
        ) %>% 
        dplyr::mutate(TOTAL_ITERATIONS = rep(iteration.method, n())) %>% 
        dplyr::select(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MEAN_MARKERS_COMMON, MISSING_DATA, SUBSAMPLE, INFERRED, NUMBER_ITERATIONS, TOTAL_ITERATIONS, MEAN_ITERATIONS, MEAN, SE, MIN, MAX, MEDIAN, QUANTILE25, QUANTILE75) %>% 
        dplyr::arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE)
    }
    
    if (assignment.analysis == "adegenet") {
      assignment.mixture.summary.stats <- suppressWarnings(
        assignment.res %>%
          dplyr::ungroup(.) %>%
          dplyr::mutate(CURRENT = factor(CURRENT)) %>% 
          dplyr::group_by(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE) %>%
          dplyr::summarise(
            NUMBER_ITERATIONS = length(ITERATIONS),
            MEAN_ITERATIONS = round((NUMBER_ITERATIONS/iteration.method)*100, 2)
          ) %>%
          dplyr::ungroup(.) %>%
          dplyr::arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE, CURRENT, INFERRED)
      )
    }
    
    # Next step is common for gsi_sim and adegenet
    # Write the tables to directory
    # assignment summary stats
    if (is.null(imputation.method)) {
      filename.assignment.sum <- stringi::stri_join("assignment.mixture.summary.results", "no.imputation", sampling.method, "tsv", sep = ".")
    } else {# with imputations
      filename.assignment.sum <- stringi::stri_join("assignment.mixture.summary.results", "imputed", sampling.method, "tsv", sep = ".")
    }
    readr::write_tsv(x = assignment.mixture.summary.stats, path = paste0(directory.subsample,filename.assignment.sum), col_names = TRUE, append = FALSE)
  } # End method random
  
  # Ranked method ------------------------------------------------------------
  if (sampling.method == "ranked") {
    message("Conducting Assignment analysis with ranked markers")
    
    # List of all individuals
    # ind.pop.df <- dplyr::ungroup(input) %>% dplyr::distinct(POP_ID, INDIVIDUALS)
    
    message("Ranking Fst with training samples...")
    holdout.individuals <- mixture.df
    
    readr::write_tsv(
      x = holdout.individuals, 
      path = paste0(directory.subsample,"holdout.individuals.tsv"), 
      col_names = TRUE, append = FALSE
    )
    message("Holdout samples = mixture samples: saved in your folder")
    
    # Going through the loop of holdout individuals
    message("Starting parallel computations for the assignment analysis: ranked markers")
    message("For progress: monitor activity in the folder...")
    
    # Ranking Fst with training dataset (keep holdout individuals out)
    message("Ranking markers based on Fst with training samples")
    fst.ranked <- assigner::fst_WC84(
      data = input,
      holdout.samples = holdout.individuals$INDIVIDUALS
    )$fst.ranked
    
    readr::write_tsv(
      x = fst.ranked, 
      path = paste0(directory.subsample, "fst_ranked.tsv"), 
      col_names = TRUE, 
      append = FALSE
    )
    if (!is.null(imputation.method)) {
      fst.ranked.imp <- assigner::fst_WC84(
        data = input.imp, 
        holdout.samples = holdout.individuals$INDIVIDUALS
      )$fst.ranked
      
      readr::write_tsv(
        x = fst.ranked.imp, 
        path = paste0(directory.subsample, "fst_ranked_imputed.tsv"), 
        col_names = TRUE, 
        append = FALSE
      )
    }
    
    # Markers numbers loop function
    message("Going throught the marker.number")
    
    if (length(marker.number) > 1) {
      # no imputation
      assignment.res <- list()
      assignment.res <- .assigner_parallel(
        X = marker.number, 
        FUN = assignment_marker_loop,
        mc.preschedule = FALSE, 
        mc.silent = FALSE, 
        mc.cleanup = TRUE,
        mc.cores = parallel.core,
        assignment.analysis = assignment.analysis,
        fst.ranked = fst.ranked,
        i = NULL,
        input = input,
        genind.object = genind.object,
        strata.df = strata.df,
        pop.levels = pop.levels,
        pop.labels = pop.labels,
        sampling.method = sampling.method,
        iteration.method = iteration.method,
        base.filename = base.filename,
        keep.gsi.files = keep.gsi.files,
        imputation.method = NULL,
        hierarchical.levels = NULL,
        parallel.core = parallel.core,
        directory.subsample = directory.subsample,
        subsample.id = subsample.id
      )
      assignment.res <- suppressWarnings(dplyr::bind_rows(assignment.res))
      
      
      # with imputations
      # todo: replace with do.call or purrr::invoke
      if (!is.null(imputation.method)) {
        assignment.res.imp <- list()
        assignment.res.imp <- .assigner_parallel(
          X = marker.number, 
          FUN = assignment_marker_loop,
          mc.preschedule = FALSE, 
          mc.silent = FALSE, 
          mc.cleanup = TRUE,
          mc.cores = parallel.core,
          assignment.analysis = assignment.analysis,
          fst.ranked = fst.ranked.imp,
          i = NULL,
          input = input.imp,
          genind.object = genind.object.imp,
          strata.df = strata.df,
          pop.levels = pop.levels,
          pop.labels = pop.labels,
          sampling.method = sampling.method,
          iteration.method = iteration.method,
          base.filename = base.filename,
          keep.gsi.files = keep.gsi.files,
          imputation.method = imputation.method,
          hierarchical.levels = hierarchical.levels,
          parallel.core = parallel.core,
          directory.subsample = directory.subsample,
          subsample.id = subsample.id
        )
        assignment.res.summary <- suppressWarnings(
          dplyr::bind_rows(assignment.res.imp) %>% 
            dplyr::bind_rows(assignment.res)
        )
      }
      assignment.res.summary <- assignment.res
    } else {
      assignment.res <- assignment_marker_loop(
        x = marker.number,
        assignment.analysis = assignment.analysis,
        fst.ranked = fst.ranked,
        i = NULL,
        input = input,
        genind.object = genind.object,
        strata.df = strata.df,
        pop.levels = pop.levels,
        pop.labels = pop.labels,
        sampling.method = sampling.method,
        iteration.method = iteration.method,
        base.filename = base.filename,
        keep.gsi.files = keep.gsi.files,
        imputation.method = NULL,
        hierarchical.levels = NULL,
        parallel.core = parallel.core,
        directory.subsample = directory.subsample,
        subsample.id = subsample.id
      )
      
      if (!is.null(imputation.method)) {
        assignment.res.imp <- assignment_marker_loop(
          x = marker.number,
          assignment.analysis = assignment.analysis,
          fst.ranked = fst.ranked.imp,
          i = NULL,
          input = input.imp,
          genind.object = genind.object.imp,
          strata.df = strata.df,
          pop.levels = pop.levels,
          pop.labels = pop.labels,
          sampling.method = sampling.method,
          iteration.method = iteration.method,
          base.filename = base.filename,
          keep.gsi.files = keep.gsi.files,
          imputation.method = imputation.method,
          hierarchical.levels = hierarchical.levels,
          parallel.core = parallel.core,
          directory.subsample = directory.subsample,
          subsample.id = subsample.id
        )
        assignment.res.summary <- suppressWarnings(
          dplyr::bind_rows(assignment.res, assignment.res.imp)
          )
      }
      assignment.res.summary <- assignment.res
    }
    assignment.res.imp <- assignment.res <- NULL
    
    # Compiling the results
    message("Compiling results")
    assignment.res.summary <- dplyr::mutate(
      .data = assignment.res.summary,
      SUBSAMPLE = rep(subsample.id, n())
    ) %>% 
      dplyr::arrange(INDIVIDUALS, MARKER_NUMBER, MISSING_DATA)
    
    # Write to the directory assignment results
    if (is.null(imputation.method)) {
      filename.assignment.res <- stringi::stri_join(
        "assignment.mixture", "no.imputation", sampling.method, "tsv", sep = ".")
    } else {# with imputations
      filename.assignment.res <- stringi::stri_join(
        "assignment.mixture", "imputed", sampling.method, "tsv", sep = ".")
    }
    readr::write_tsv(x = assignment.res.summary,
                     path = paste0(directory.subsample, filename.assignment.res),
                     col_names = TRUE, append = FALSE)
  } # End of ranked thl method
  
  return(assignment.res.summary)
} # End assignment_function
