# Compute Weir and Cockerham (1984) Fst

#' @name fst_WC84

#' @title A fast implementation of Weir and Cockerham (1984) Fst/Theta
#' (overall and paiwise estimates)

#' @description The function calculates Weir and Cockerham (1984)
#' Fst for diploid genomes. Both overall and pairwise Fst can be estimated with
#' confidence intervals based on bootstrap of markers (resampling with replacement).
#' The function gives identical results \emph{at the 9th decimal} when tested
#' against \code{genet.dist} in \code{hierfstat}. Using the
#' argument \code{snprelate = TRUE} will compute the Fst with
#' \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}. This implementation
#' gives slightly upward bias values but provided the fastest computations I know,
#' but it doesn't compute confidence intervals, for now.
#' For an R implementation, \code{\link{fst_WC84}} is very fast.
#' The computations takes advantage of \pkg{dplyr}, \pkg{tidyr}, \pkg{purrr},
#' \pkg{parallel} and \pkg{SNPRelate}.
#' The impact of unbalanced design on estimates can be tested by using the
#' subsample argument (see advance mode section).
#'
#' \emph{Special concerns for genome-wide estimate and filtering bias}
#'
#' During computation, the function first starts by keeping only the polymorphic
#' markers in common between the populations. Keep this in mind when filtering
#' your markers to use this function characteristic strategically to get
#' better genome-wide estimate. This is even more important when your project
#' involves more than 2 populations that evolved more by neutral processes
#' (e.g. genetic drift) than by natural selection (see the vignette for more details).

#' @note \strong{Negative Fst} are technical artifact of the computation
#' (see Roesti el al. 2012) and are automatically replaced with zero inside
#' this function.
#'
#' \strong{Why no p-values ?}
#'
#' There is no null hypothesis testing with \emph{P}-values and its rarely if ever
#' the appropriate model in population genomics, despite its popularity
#' with molecular ecologists interested in population differentiation.
#'
#' "The important scientific question is the real magnitude of the differentiation,
#' not the smallness of the \emph{P} value" -Lou Jost
#'
#' Confidence intervals provided with the \emph{F}-statistics
#' enables more reliable conclusions about the biological trends in the data.
#' A confidence interval describes the statistical uncertainty of the \emph{F}-statistics
#' estimate. If the confidence interval include zero, then the null hypothesis
#' cannot be rejected. If the confidence interval does not include zero, the null
#' hypothesis can be rejected and you can also have an appreciation of the real
#' magnitude of the statistical differentiation whether its large or small.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link[radiator]{tidy_genomic_data}}.
#' You can also use this function to filter your dataset using
#' whitelist of markers, blacklist of individuals and genotypes.

#' @param snprelate (optional, logical) Use \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}
#' to compute the Fst.
#' It's the fastest computation I've seen so far!
#'
#' However, testing with different RADseq datasets as shown several upward bias
#' with \code{SNPRelate::snpgdsFst} (last version tested was v.1.16.0).
#' I compared the results with assigner, hierfstat and strataG
#' (results available upon request).
#' The SNPRelate author as not given me good reason to belive the issue is fully
#' resolved, consequently, the option is no longer available, until further notice.
#' Default: \code{snprelate = FALSE}

#' @param pop.levels (optional, string) This refers to the levels in a factor. In this
#' case, the id of the pop.
#' Use this argument to have the pop ordered your way instead of the default
#' alphabetical or numerical order. e.g. \code{pop.levels = c("QUE", "ONT", "ALB")}
#' instead of the default \code{pop.levels = c("ALB", "ONT", "QUE")}.
#' Default: \code{pop.levels = NULL}.

#' @param strata (optional, data frame) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' If a \code{strata} file is specified, the strata file will have
#' precedence over any grouping found data file (\code{data}).
#' The \code{STRATA} column can be any hierarchical grouping.
#' Default: \code{strata = NULL}.


#' @param pairwise (optional, logical) With \code{pairwise = TRUE}, the
#' pairwise WC84 Fst is calculated between populations.
#' Default: \code{pairwise = FALSE}.

#' @param ci (optional, logical) Compute bootstrapped confidence intervals.
#' Default: \code{ci = FALSE}.

#' @param iteration.ci (optional, integer) The number of iterations for
#' the boostraps (resampling with replacement of markers).
#' Default: \code{iteration.ci = 100}.

#' @param quantiles.ci (optional, double)
#' The quantiles for the bootstrapped confidence intervals.
#' Default: \code{quantiles.ci = c(0.025,0.975)}.

#' @param heatmap.fst (logical) Generate a heatmap with the Fst values in
#' lower matrix and CI in the upper matrix.
#' The heatmap can also be generated separately after the Fst
#' analysis using the separate function: \code{\link{heatmap_fst}}.
#' Default: \code{heatmap.fst = FALSE}.

#' @param digits (optional, integer) The number of decimal places to be used in
#' results.
#' Default: \code{digits = 9}.

#' @param parallel.core (optional, integer) The number of core for parallel computation
#' of pairwise Fst. See also the advance mode section below.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.


#' @param filename (optional, character) Give filename prefix for the output
#' directory, this will trigger saving results.
#' Default: \code{filename = "fst_WC84"}.

#' @param ... other parameters passed to the function.

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#'
#' \item \code{filter.monomorphic} (logical, optional) By default monomorphic
#' markers present in the dataset are removed (and it should stay that way...).
#' Default: \code{filter.monomorphic = TRUE}.
#'
#' \item \code{holdout.samples} (optional, data frame) Samples that don't participate in the Fst
#' computation (supplementary). Data frame with one column \code{INDIVIDUALS}.
#' This argument is used inside assignment analysis.
#' Default: \code{holdout.samples = NULL}.
#'
#' \item \code{subsample} (Integer or character)
#' With \code{subsample = 36}, 36 individuals in each populations are chosen
#' randomly to represent the dataset. With \code{subsample = "min"}, the
#' minimum number of individual/population found in the data is used automatically.
#' Default is no subsampling, \code{subsample = NULL}.
#' \item \code{iteration.subsample} (Integer) The number of iterations to repeat
#' subsampling.
#' With \code{subsample = 20} and \code{iteration.subsample = 10},
#' 20 individuals/populations will be randomly chosen 10 times.
#' Default: \code{iteration.subsample = 1}.
#'
#'
#' \item \code{calibrate.alleles} (logical)
#' Un-calibrated alleles can bias estimate and by default the function expect that
#' the REF/ALT alleles are calibrated. Using \code{calibrate.alleles = TRUE},
#' can take a bit more time.
#' Default: \code{calibrate.alleles = FALSE}.
#' }



#' @return The function returns a list with several objects.
#' When sumsample is selected the objects end with \code{.subsample}.
#' \itemize{
#'  \item \code{$subsampling.individuals}: the combinations of individuals and subsamples,
#'  \item \code{$sigma.loc}: the variance components per locus, with
#'       (\code{lsiga}: among populations,
#'       \code{lsigb}: among individuals within populations,
#'       \code{lsigw}: within individuals)
#'  \item \code{$fst.markers}: the fst by markers,
#'  \item \code{$fst.ranked}: the fst ranked,
#'  \item \code{$fst.overall}: the mean fst overall markers and the number of markers
#'  \item \code{$fis.markers}: the fis by markers,
#'  \item \code{$fis.overall}: the mean fis overall markers and the number of markers,
#'  \item \code{$fst.plot}: the histogram of the overall Fst per markers,
#'  \item \code{$pairwise.fst}: the pairwise fst in long/tidy data frame and the number of markers ,
#'  \item \code{$pairwise.fst.upper.matrix}: the pairwise fst in a upper triangle matrix,
#'  \item \code{$pairwise.fst.full.matrix}: the pairwise fst matrix (duplicated upper and lower triangle),
#'  \item \code{$pairwise.fst.ci.matrix}: matrix with pairwise fst in the upper triangle
#'  and the confidence intervals in the lower triangle.
#'  \item when subsample is selected \code{$pairwise.fst.subsample.mean} is a summary
#'  of all pairwise comparisons subsample. The mean is calculated accross summary
#'  statistics.
#' }

#' @export
#' @rdname fst_WC84
#' @examples
#' \dontrun{
#' wombat.fst.pairwise <- fst_WC84(
#'     data = "wombat.filtered.tidy.tsv",
#'     pop.levels = c("ATL", "MLE", "BIS", "PMO", "SOL", "TAS", "ECU"),
#'     pairwise = TRUE,
#'     ci = TRUE,
#'     iteration.ci = 10000,
#'     quantiles.ci = c(0.025,0.975),
#'     parallel.core = 8,
#'     verbose = TRUE,
#'     filename = "wombat",
#'     heatmap.fst = TRUE
#' )
#'
#' # To get the overall Fst estimate:
#' wombat.fst.pairwise$fst.overall
#'
#' # To get the Fst plot:
#' wombat.fst.pairwise$fst.plot
#'
#' #To get the pairwise Fst values with confidence intervals in a data frame:
#' df <- wombat.fst.pairwise$pairwise.fst
#' }

#' @references Excoffier L, Smouse PE, Quattro JM.
#' Analysis of molecular variance inferred from metric distances among
#' DNA haplotypes: application to human mitochondrial DNA restriction data.
#' Genetics. 1992;131: 479-491.
#' @references Jost L, D vs. GST: Response to Heller and Siegismund (2009) and
#' Ryman and Leimar (2009).
#' Molecular Ecology. 2009; 18:10 2088-2091.
#' https://doi.org/10.1111/j.1365-294x.2009.04186.x
#' @references Meirmans PG, Van Tienderen PH (2004) genotype and genodive:
#' two programs for the analysis of genetic diversity of asexual organisms.
#' Molecular Ecology Notes, 4, 792-794.
#' @references Michalakis Y, Excoffier L.
#' A generic estimation of population
#' subdivision using distances between alleles with special reference for
#' microsatellite loci.
#' Genetics. 1996;142: 1061-1064.
#' @references Weir BS, Cockerham CC (1984) Estimating F-Statistics for the
#' Analysis of Population Structure.
#' Evolution, 38, 1358-1370.
#' @references Roesti M, Salzburger W, Berner D. (2012)
#' Uninformative polymorphisms bias genome scans for signatures of selection.
#' BMC Evol Biol., 12:94. doi:10.1111/j.1365-294X.2012.05509.x
#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606

#' @seealso
#' From \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive} manual:
#' \emph{'In general, rather than to test differentiation between all pairs of
#' populations,
#' it is advisable to perform an overall test of population differentiation,
#' possibly using a hierarchical population structure, (see AMOVA)'}.
#' To compute an AMOVA, use \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive}
#' or \code{Phi_st_Meirmans} in \code{mmod}.
#'
#' \href{https://github.com/jgx65/hierfstat/}{hierfstat}
#'
#' For Fisher's exact test and p-values per markers
#' see \code{mmod} \code{diff_test}.
#'
#' \strong{Vignette for this function:} \href{https://www.dropbox.com/s/tiq4yenzmgzc2f5/fst_confidence_intervals.html?dl=0}{how to do the pairwise and overall Fst with confidence intervals and build the phylogenetic tree}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# Fst function: Weir & Cockerham 1984
fst_WC84 <- function(
    data,
    snprelate = FALSE,
    strata = NULL,
    pop.levels = NULL,
    pairwise = FALSE,
    ci = FALSE,
    iteration.ci = 100,
    quantiles.ci = c(0.025,0.975),
    heatmap.fst = FALSE,
    digits = 9,
    filename = "fst_WC84",
    parallel.core = parallel::detectCores() - 2,
    verbose = FALSE,
    ...
) {
  ## test
  # data
  # snprelate = FALSE
  # pop.levels = NULL
  # strata = NULL
  # holdout.samples = NULL
  # pairwise = FALSE
  # ci = FALSE
  # iteration.ci = 100
  # quantiles.ci = c(0.025, 0.975)
  # subsample = NULL
  # iteration.subsample = 1
  # digits = 9
  # parallel.core = parallel::detectCores() - 1
  # verbose = TRUE
  # filename = "coral_fst"
  # heatmap.fst = FALSE
  # calibrate.alleles = FALSE

  # Cleanup---------------------------------------------------------------------
  assigner_function_header(f.name = "fst_WC84", verbose = verbose)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date/time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  opt.digits <- getOption("digits")
  options(digits = digits, scipen = 999)
  timing <- assigner_tic()
  res <- list() # where the results will be stored
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(options(digits = opt.digits), add = TRUE)
  on.exit(assigner_toc(timing), add = TRUE)
  on.exit(assigner_function_header(f.name = "fst_WC84", start = FALSE, verbose = verbose), add = TRUE)

  # Function call and dotslist -------------------------------------------------
  rad.dots <- assigner::assigner_dots(
    func.name = as.list(sys.call())[[1]], # Get name of the calling function
    fd = rlang::fn_fmls_names(),# Grab formal argument names of the function
    args.list = as.list(environment()),# Capture current environment's arguments
    dotslist = rlang::dots_list(
      ..., # Capture the `...` content
      .homonyms = "error", # Error on duplicate names
      .check_assign = TRUE), # Ensure arguments can be assigned
    keepers = c("holdout.samples", "subsample",
                "iteration.subsample", "blacklist.id",
                "calibrate.alleles"),# Expected dot arguments
    verbose = FALSE
  )

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data is missing")
  if (!ci && heatmap.fst) {
    heatmap.fst <- FALSE
    if (verbose) message("\nconfidence intervals not selected, heatmap.fst: FALSE\n")
  }

  # filename & folder ----------------------------------------------------------
  if (is.null(filename)) filename <- "fst_WC84"

  path.folder <- radiator::generate_folder(
    rad.folder = filename,
    path.folder = NULL,
    file.date = file.date,
    verbose = verbose)

  # write the dots file
  radiator::write_radiator_tsv(
    data = rad.dots,
    path.folder = path.folder,
    filename = "assigner_fst_WC84_args",
    date = TRUE,
    internal = FALSE,
    write.message = "Function call and arguments stored in: ",
    verbose = verbose
  )

  snprelate <- FALSE
  if (snprelate) {
    # Check that snprelate is installed
    if (!"SNPRelate" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install SNPRelate for this option:\n
                 install.packages("BiocManager")
                 BiocManager::install("SNPRelate")')
    }
    rlang::abort("Until the bias observed with SNPRelate is resolved, the option is unavailable.")
  }


  # plan -----------------------------------------------------------------------

  # 1. Detecting and importing the data.
  # 2. Prepare strata and population levels
  # 3. Stripping and selecting data.
  # 4. Create subsample list and write to directory
  # 5. Write temp data set
  # 6. Compute global Fst
  # 7. Compute pairwise Fst
  # 8. Compiling results
  # 9. Write and plot


  # 1. Detecting and importing the data ----------------------------------------
  data <- import_data(
    data = data,
    calibrate.alleles = calibrate.alleles,
    verbose = verbose
  )

  # 2. Prepare strata and population levels ------------------------------------
  data <- prep_strata(
    data = data,
    strata = strata,
    verbose = verbose
  )

  # 3. Stripping and selecting data---------------------------------------------
  env.arg <- rlang::current_env()
  data %<>%
    radiator::strip_rad(
      x = .,
      env.arg = env.arg,
      keep.strata = TRUE,
      verbose = FALSE
    ) %>%
    dplyr::select(tidyselect::one_of(c("ID_SEQ", "STRATA_SEQ", "M_SEQ", "GT")))

  # The function will make the bk object...
  # strata.bk
  # markers.meta.bk
  # genotypes.meta.bk
  # pop.levels.bk

  pop.levels <- unique(data$STRATA_SEQ) # integers
  # pop.levels <- as.character(unique(data$STRATA_SEQ))

  # 4. Generating subsample list and write to directory ----------------------------
  subsample.list <- create_subsample_list(
    subsample = subsample,
    strata.bk = strata.bk,
    iteration.subsample = iteration.subsample,
    path.folder = path.folder,
    verbose = verbose
  )

  # decide if plot.fst will be generated in compute_fst (to manage the iterations)
  plot.fst <- TRUE
  if (iteration.subsample > 1) plot.fst <- FALSE


  # 5. Write temp data set -------------------------------------------------------
  if (verbose) cli::cli_progress_step(msg = "Writing temporary data")

  args.for.prep <- extract_matching_args(
    from.env = environment(),
    to.fn = write_temp_data_fst,
    .evaluate = TRUE,
    .exclude = c("data", "strata", "subsample.list")
  )
  # names(args.for.prep)

  options(cli.progress_handlers = "progressr", cli.progress_bar_width = 80)
  progressr::handlers(global = TRUE)
  # progressr::handlers(progressr::handler_txtprogressbar(char = "\U1F41F"))
  progressr::handlers(progressr::handler_txtprogressbar(char = " ><((('> "))


  fst.temp.files <- purrr::map(
    .x = subsample.list,
    .f = function(subsample) {
      rlang::exec(
        assigner::write_temp_data_fst,
        subsample.list = subsample,  # Explicitly pass the current element of subsample.list
        data = data,                 # Explicitly pass the data
        strata = strata.bk,          # Explicitly pass the strata
        !!!args.for.prep
      )
    },
    .progress = list(
      type = "iterator",
      format = "{cli::pb_name} Writing temp data: {cli::pb_bar} {cli::pb_percent}",
      clear = TRUE
    )
  )
  args.for.prep <- NULL

  options(cli.progress_handlers = NULL)

  if (pairwise) {
    if (!is.factor(data$STRATA_SEQ)) data$STRATA_SEQ <- factor(data$STRATA_SEQ)
    # pop list
    pop.list <- levels(data$STRATA_SEQ)
  }

  data <- NULL # no longer needed ... we have the files...

  # 6. Compute global Fst ------------------------------------------------------
  if (verbose) cli::cli_progress_step(msg = "Computing global Fst")

  # Using purrr::map vs furrr::future_map

  # Extract only the matching arguments from the current environment
  # args.for.fst <- extract_matching_args(
  #   from.env = environment(),
  #   to.fn = assigner::compute_fst
  # )
  # names(args.for.fst)

  # global.fst <- purrr::map(
  #   .x = fst.temp.files,
  #   .f = function(file) {
  #     rlang::exec(
  #       assigner::compute_fst,
  #       x = file,
  #       !!!args.for.fst
  #     )
  #   },
  #   .progress = list(
  #     type = "iterator",
  #     format = "{cli::pb_name} Writing temp data: {cli::pb_bar} {cli::pb_percent}",
  #     clear = TRUE
  #   )
  # )

  # furrr::future_map
  # args.for.fst doesnt function with furrr...

  # Set up parallel processing plan
  future::plan(future::multisession, workers = parallel::detectCores() - 1)

  # Ensure the plan is reset after this block (optional but recommended)
  on.exit(future::plan(future::sequential), add = TRUE)

  res <- list()
  global <- furrr::future_map(
    .x = fst.temp.files,
    .f = function(file) {
      assigner::compute_fst(
        x = file,
        ci = ci,
        iteration.ci = iteration.ci,
        quantiles.ci = quantiles.ci,
        digits = digits,
        path.folder = path.folder,
        plot.fst = plot.fst
      )
    },
    .options = furrr::furrr_options(seed = TRUE),
    .progress = TRUE
  )
  res <- append(res, global)
  global <- NULL

  # 7. Compute pairwise Fst -------------------------------------------------------

  # by default
  res$pairwise.fst <- list(pairwise.fst = "pairwise fst not selected")
  if (pairwise) {
    if (verbose) cli::cli_progress_step(msg = "Computing pairwise Fst")

    # all combination of populations
    pop.pairwise <- utils::combn(pop.list, 2, simplify = FALSE)
    pop.list <- NULL
    # if (verbose) message("Number of pairwise computations: ", length(pop.pairwise))
    # if (verbose) cli::cli_progress_message("Number of pairwise computations: {length(pop.pairwise)}")

    future::plan(future::multisession, workers = parallel::detectCores() - 1)
    on.exit(future::plan(future::sequential), add = TRUE)

    pairwise.fst <- furrr::future_map(
      .x = fst.temp.files,
      .f = function(file) {
        pairwise_fst(
          x = file,
          pop.pairwise = pop.pairwise,
          ci = ci,
          iteration.ci = iteration.ci,
          quantiles.ci = quantiles.ci,
          digits = digits,
          path.folder = path.folder
        )
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    ) %>%
      dplyr::bind_rows()
    res$pairwise.fst <- list(pairwise.fst = pairwise.fst)
    pop.pairwise <- NULL
  }# END pairwise

  options(cli.progress_handlers = NULL)

  # delete the temp files
  temp.files.paths <- unlist(fst.temp.files)
  purrr::walk(temp.files.paths[file.exists(temp.files.paths)], unlink)
  if (verbose) cli::cli_progress_message("Removing temporary files")

  # # str(res)
  # # check that the calculations did the exact number ...
  # if (length(res) != iteration.subsample) {
  #   rlang::abort(message("The function didn't complete the number of iterations, contact author"))
  # }
  # subsample.list <- NULL

  # 8. Compiling results--------------------------------------------------------
  if (verbose) message("Preparing results")
  nms <- res %>% purrr::map(names) %>% purrr::reduce(union)

  # These are the objects in nms:
  # sigma.loc
  # fst.markers
  # fst.ranked
  # fst.overall
  # fis.markers
  # fis.overall
  # fst.plot

  # when pairwise is selected:
  # pairwise.fst
  # pairwise.fst
  # pairwise.fst.upper.matrix
  # pairwise.fst.full.matrix
  # pairwise.fst.ci.matrix

  # defaults:
  pairwise.fst <- "pairwise fst not selected"
  upper.mat.fst <- "pairwise fst not selected"
  full.mat.fst <- "pairwise fst not selected"
  pairwise.fst.ci.matrix <- "pairwise fst not selected"

  if (pairwise) {
    # Tibble of pairwise Fst
    pairwise.fst <- fst_stats(
      x = "pairwise.fst",
      l = res,
      digits = digits,
      m = markers.meta.bk,
      s = dplyr::distinct(strata.bk, STRATA, STRATA_SEQ),
      iteration.subsample = iteration.subsample
    ) %>%
      dplyr::bind_rows()

    # FST upper matrix
    pairwise.fst.upper.matrix <- make_upper_matrix(
      data = dplyr::select(pairwise.fst, POP1, POP2, FST),
      value.col = FST,
      fill.value = " "
    )

    # Full symmetric matrix
    # merge upper and lower matrix

    pairwise.fst.full.matrix <- make_full_symmetric_matrix(
      upper.matrix = pairwise.fst.upper.matrix,
      diagonal.value = "0"
    )



    # default
    pairwise.fst.ci.matrix <- "confidence intervals not selected"

    # CI matrix logic
    if (ci) {
      upper.mat.ci <- pairwise.fst %>%
        dplyr::select(POP1, POP2, CI_LOW, CI_HIGH) %>%
        tidyr::unite(CI, CI_LOW, CI_HIGH, sep = " - ") %>%
        make_upper_matrix(
          data = .,
          value.col = CI,
          fill.value = ""
        )

      pairwise.fst.ci.matrix <- pairwise.fst.full.matrix
      pairwise.fst.ci.matrix[upper.tri(pairwise.fst.ci.matrix)] <-
        upper.mat.ci[upper.tri(upper.mat.ci)]
      upper.mat.ci <- NULL
    }
  }

  # update res list:
  res$pairwise.fst <- pairwise.fst
  res$pairwise.fst.upper.matrix <- pairwise.fst.upper.matrix
  res$pairwise.fst.full.matrix <- pairwise.fst.full.matrix
  res$pairwise.fst.ci.matrix <- pairwise.fst.ci.matrix


  # remove pairwise.fst and manage the rest
  nms %<>% purrr::discard(.p = nms %in% "pairwise.fst")

  # TODO
  # "sigma.loc"
  # "fst.markers"
  # "fst.ranked"
  # "fst.overall"
  # "fis.markers"
  # "fis.overall"
  # "fst.plot"

  res.temp <- purrr::map(
    .x = nms,
    .f = fst_stats,
    l = res,
    digits = digits,
    m = markers.meta.bk,
    s = dplyr::distinct(strata.bk, STRATA, STRATA_SEQ),
    iteration.subsample = iteration.subsample
  ) %>%
    purrr::flatten(.)

  res <- append(res, res.temp)
  res.temp <- NULL # unsused object

  # 9. Write and plot-----------------------------------------------------------

  # plot fst for overall iterations.
  # the plot is not generated with iterations = 1
  if (iteration.subsample > 1) {
    if (!is.null(res$fst.markers) && "MEAN" %in% colnames(res$fst.markers)) {
      res$fst.plot <- assigner::plot_fst_distribution(
        data = dplyr::rename(.data = res$fst.markers, FST = MEAN)
      )
    } else {
      message("res$fst.markers is NULL or does not contain the 'MEAN' column.")
    }

  }

  # write results
  if (!is.null(filename)) {
    purrr::walk(
      .x = list("sigma.loc", "fst.markers", "fst.ranked", "fst.overall",
                "fis.markers", "fis.overall", "pairwise.fst"),
      .f = fst_write, list.sub = res, path.folder = path.folder
    )
    # fst.plot
    ggplot2::ggsave(
      filename = file.path(path.folder, "fst.plot.pdf"),
      plot = res$fst.plot,
      width = 15, height = 10,
      dpi = 300, units = "cm", device = "pdf", limitsize = FALSE,
      useDingbats = FALSE
    )
    saveRDS(
      object = res$pairwise.fst.upper.matrix,
      file = file.path(path.folder, "pairwise.fst.upper.matrix.RData"))
    saveRDS(
      object = res$pairwise.fst.full.matrix,
      file = file.path(path.folder, "pairwise.fst.full.matrix.RData"))
    saveRDS(
      object = res$pairwise.fst.ci.matrix,
      file = file.path(path.folder, "pairwise.fst.ci.matrix.RData"))
  }

  # heatmap.fst
  if (heatmap.fst) {
    res$heatmap.fst <- heatmap_fst(
      pairwise.fst.full.matrix = res$pairwise.fst.full.matrix,
      pairwise.fst.ci.matrix = res$pairwise.fst.ci.matrix,
      digits = digits,
      path.folder = path.folder,
      filename = filename)
  } else {
    res$heatmap.fst <- "confidence intervals not selected, heatmap.fst => FALSE"
  }

  # End -------------------------------------------------------------------
  if (verbose) {
    message(stringi::stri_pad_both(str = "FST RESULTS", width = 80L, pad = "#"))

    # default:
    sep1 <- ""
    sep2 <- ""
    sep3 <- ""

    mean.fst <- res$fst.overall[[if (is.null(subsample)) "FST" else "MEAN"]]

    if (ci) {
      sep1 <- " ["
      sep2 <- " - "
      sep3 <- "]"
      min.fst <- res$fst.overall$MIN
      max.fst <- res$fst.overall$MAX
    }
    fst.message <- stringi::stri_c(mean.fst, sep1, min.fst, sep2, max.fst, sep3)
    message("Fst (overall): ", fst.message)
  }
  # if (verbose) cli::cli_progress_done()
  return(res)
}#END fst_WC84

# Internal Nested Functions to compute WC84 Fst --------------------------------

# Import and Standardize Genomic Data ------------------------------------------
#' @title Import and Standardize Genomic Data
#' @description Import and format genomic data into a tidy format compatible with downstream analyses in \pkg{radiator}. This function detects the input data type (tidy data frame, GDS file, or SeqVarGDSClass) and performs the appropriate transformation.
#' @param data Genomic dataset. Can be a tibble/data frame in wide format, a GDS file path, or a \code{SeqVarGDSClass} object.
#' @param calibrate.alleles Logical. If \code{TRUE}, force recalibration of REF/ALT allele genotypes using \code{\link[radiator]{calibrate_alleles}}. Default: \code{FALSE}.
#' @param verbose Logical. Show detailed messages during import. Default: \code{FALSE}.
#' @return A tidy data frame (tibble) with standardized genotype format and required metadata (e.g., GT column).
#' @seealso \code{\link[radiator]{detect_genomic_format}}, \code{\link[radiator]{tidy_wide}}, \code{\link[radiator]{calibrate_alleles}}, \code{\link[radiator]{read_rad}}, \code{\link[radiator]{gds2tidy}}
#' @rdname import_data
#' @export
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # From tidy data
#' tidy_data <- import_data(data = my_tibble)
#'
#' # From GDS file
#' tidy_data <- import_data(data = "myfile.gds")
#'
#' # With calibration
#' tidy_data <- import_data(data = my_tibble, calibrate.alleles = TRUE)
#' }

import_data <- function(data, calibrate.alleles = FALSE, verbose = FALSE) {
  data.type <- radiator::detect_genomic_format(data)

  if (verbose) cli::cli_progress_step(msg = "Importing data...")

  if (!data.type %in% c("tbl_df", "SeqVarGDSClass", "gds.file")) {
    rlang::abort("Input not supported for this function: read function documentation")
  }
  # GDS or Tidy
  if (data.type == "tbl_df") {
    data %<>% radiator::tidy_wide(data = ., import.metadata = TRUE)

    if (!rlang::has_name(data, "GT") || calibrate.alleles) {
      if (verbose) cli::cli_progress_step(msg = "generating the necessary genotype format")
      data %<>% radiator::calibrate_alleles(data = ., gt = TRUE, verbose = TRUE) %$% input
    }

  } else if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    radiator::radiator_packages_dep(package = "SeqArray", cran = FALSE, bioc = TRUE)

    if (data.type == "gds.file") {
      data %<>% radiator::read_rad(data = ., verbose = TRUE)
    }

    data %<>% radiator::gds2tidy(gds = ., pop.id = FALSE, gt = TRUE)
  } else {
    rlang::abort("Input not supported for this function: read function documentation")
  }
  # cli::cli_progress_done()
  return(data)
}#END import_data

# Prepare stratification data ------------------------------------------

#' @title Prepare stratification data
#' @description Read and join strata information to the main dataset, standardising population levels.
#' @rdname prep_strata
#'
#' @param data A tibble containing genetic data.
#' @param strata Path to the strata file or a tibble containing strata information.
#' @param pop.levels (optional) A character vector specifying the population levels to use for factor levels of `STRATA`. Default: \code{NULL}.
#' @param blacklist.id (optional) A vector of individual IDs to blacklist. Default: \code{NULL}.
#' @param verbose (logical) Show progress messages. Default: \code{FALSE}.
#'
#' @return A tibble with updated `STRATA` column and joined strata info.
#'
#' @keywords internal
#' @export

prep_strata <- function(data, strata, pop.levels = NULL, blacklist.id = NULL, verbose = FALSE) {
  if (verbose) cli::cli_progress_step(msg = "Preparing the stratification")

  # Strata
  strata <- radiator::read_strata(
    strata = strata,
    pop.id = TRUE,
    blacklist.id = blacklist.id,
    pop.levels = NULL,
    verbose = verbose) %$%
    strata

  # population levels and strata
  if (!is.null(strata)) {
    data <- radiator::join_strata(
      data = data,
      strata = strata,
      pop.id = FALSE,
      verbose = FALSE
    )
  }

  # no longer use POP_ID
  if (rlang::has_name(data, "POP_ID") && !rlang::has_name(data, "STRATA")) {
    data %<>% dplyr::rename(STRATA = POP_ID)
  }

  check.levels <- as.character(unique(data$STRATA))
  same.levels <- FALSE
  if (!is.null(pop.levels)) {
    same.levels <- identical(pop.levels, check.levels)
  } else {
    pop.levels <- check.levels
  }

  if (!same.levels) {
    data %<>%
      dplyr::mutate(STRATA = factor(x = STRATA, levels = pop.levels)) %>%
      dplyr::arrange(STRATA)
  }
  same.levels <- check.levels <- NULL
  # if (verbose) cli::cli_progress_done()
  return(data)
} #END prep_strata


# Create a list of subsampled individuals --------------------------------------
#' @title Create a list of subsampled individuals
#' @description Generate a list of subsampled individuals from strata for replicated subsampling. Optionally writes the output to file.
#' @rdname create_subsample_list
#'
#' @param subsample Integer or `"min"`: number of individuals to subsample per stratum, or `"min"` to use the minimum stratum size. Default: \code{NULL}.
#' @param strata.bk A tibble of the original strata information, with \code{STRATA_SEQ} as a grouping column.
#' @param iteration.subsample Number of subsampling replicates. Default: \code{1L}.
#' @param path.folder Path to the directory where subsample individuals file will be written. Default: \code{NULL}.
#' @param verbose (logical) Show messages during execution. Default: \code{FALSE}.
#'
#' @return A list of tibbles, each containing subsampled individuals for one iteration.
#'
#' @keywords internal
#' @export


create_subsample_list <- function(subsample, strata.bk, iteration.subsample, path.folder = NULL, verbose = FALSE) {

  # create the subsampling list
  if (!is.null(subsample) && !is.numeric(subsample)) {
    if (verbose) cli::cli_progress_step(msg = "Generating subsample list")

    if (subsample == "min") {
      subsample <- NULL
      subsample <- strata.bk %>%
        dplyr::group_by(STRATA_SEQ) %>%
        dplyr::tally(.) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::ungroup(.) %>%
        dplyr::pull(n) %>%
        unique
    }
  }

  # Generate subsample list
  subsample.list <- purrr::map(
    .x = 1:iteration.subsample,
    .f = subsampling_data,
    strata = strata.bk,
    subsample = subsample,
    random.seed = NULL
  )

  # keep track of subsampling individuals and write to directory
  if (!is.null(subsample)) {
    # if (verbose) message("Subsampling size: ", subsample)
    # subsampling.size.mes <- paste0("Subsampling size: ", subsample)
    # if (verbose) cli::cli_progress_message(subsampling.size.mes)
    if (verbose) cli::cli_progress_message("Subsampling size: {subsample}")
    subsample.list %>%
      dplyr::bind_rows() %>%
      readr::write_tsv(
        x = .,
        file = file.path(path.folder, "subsampling.individuals.tsv")
      )
  } # End subsampling
  # if (verbose) cli::cli_progress_done()
  # message("it works")
  return(subsample.list)
}#End create_subsample_list

# Write Temporary Data to Parquet ----------------------------------------------
#' @rdname write_temp_data_fst
#' @title Write Temporary Data to Parquet
#'
#' @description
#' This function filters the input data based on subsample information, removes missing genotypes, and optionally removes holdout individuals. It then writes the resulting data to a Parquet file.
#'
#' @param subsample.list A data frame or tibble containing information about the subsamples. It should include a column `SUBSAMPLE` used to identify each subsample.
#' @param data A data frame or tibble containing the dataset that will be filtered based on the subsample and holdout samples.
#' @param holdout.samples A vector of IDs or individual identifiers to be excluded from the dataset (optional).
#' @param strata A data frame or tibble containing information about the strata, used to identify and filter holdout samples.
#' @param file.date A character string representing the date (or timestamp) to be included in the filename for the Parquet file.
#' @param path.folder The folder path where the resulting Parquet file will be saved.
#' @param verbose A logical indicating whether to display progress messages. Defaults to `TRUE`.
#' @keywords internal
#' @export

#' @return A character string representing the path to the saved Parquet file.
#'
#' @details
#' The function filters the `data` by the `ID_SEQ` values found in the `subsample.list`. Missing genotypes (denoted by `GT == "000000"`) are removed. If `holdout.samples` are provided, those individuals are removed from the `data` based on matching IDs in the `strata` data frame.
#' The filtered data is then written to a Parquet file with the specified filename and path.
#'
#'
#' @examples
#' \dontrun{
#' # Example usage
#' write_temp_data_fst(
#'   subsample.list = subsample_info,
#'   data = my_data
#' )
#' }

write_temp_data_fst <- function(
    subsample.list = NULL,
    data,
    holdout.samples = NULL,
    strata = NULL,
    file.date = NULL,
    path.folder = NULL,
    verbose = TRUE
) {
  i <- unique(subsample.list$SUBSAMPLE)
  # message("Subsample: ", i)
  # cli::cli_progress_update()
  # cli::cli_progress_message("subsample: {i}")

  # Keep only subsample and remove missing genotypes
  data %<>%
    dplyr::filter(ID_SEQ %in% subsample.list$ID_SEQ) %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::mutate(ITERATIONS = i)
  subsample.list <- NULL #unused object

  # Remove holdout individuals if provided
  if (!is.null(holdout.samples)) {
    # if (verbose) cli::cli_progress_message("Removing holdout individuals")
    # message("Removing holdout individuals\nFst computation...")

    match.col <- if (is.numeric(holdout.samples) || any(holdout.samples %in% unique(data$ID_SEQ))) {
      "ID_SEQ"
    } else {
      "INDIVIDUALS"
    }

    holdout.samples <- strata %>%
      dplyr::filter(.data[[match.col]] %in% holdout.samples) %>%
      dplyr::pull(ID_SEQ)

    data %<>% dplyr::filter(!ID_SEQ %in% holdout.samples)
  }
  holdout.samples <- match.col <- strata <- NULL
  filename <- paste0("assigner_temp_file", file.date, "_", i, ".arrow.parquet") %>%
    file.path(path.folder, .)
  arrow::write_parquet(x = data, sink = filename)
  return(filename)
}#End write_temp_data_fst










# compute_fst------------------------------------------------------------------

#' @title compute_fst
#' @description main function
#' @rdname compute_fst
#' @export
#' @keywords internal

compute_fst <- carrier::crate(function(
    x,
    ci = FALSE,
    iteration.ci = 100,
    quantiles.ci = c(0.025,0.975),
    digits = 9,
    path.folder = NULL,
    plot.fst = TRUE,
    ...
) {
  # TEST
  # ci = FALSE
  # iteration.ci = 100
  # quantiles.ci = c(0.025,0.975)
  # digits = 9
  # path.folder = NULL
  # plot = TRUE

  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
  `%$%` <- magrittr::`%$%`
  # cli::cli_progress_step(msg = "Global fst")

  # Check if input_data is a file path or an object
  if (inherits(x, "character") && file.exists(x)) {
    # If it's a file path, read the file (assuming it's a Parquet file)
    # message("Reading data from file: ", input_data)
    x <- arrow::read_parquet(x)
  }

  i <- 1L #default
  if (rlang::has_name(x, "ITERATIONS")) {
    i <- unique(x$ITERATIONS)
    x %<>% dplyr::select(-ITERATIONS)
  }


  # Removing monomorphic markers
  bl <- x %>%
    dplyr::filter(GT != "000000") %>%
    dplyr::transmute(
      M_SEQ,
      A1 = stringi::stri_sub(GT, 1, 3),
      A2 = stringi::stri_sub(GT, 4, 6)
    ) %>%
    tidyr::pivot_longer(
      cols = c(A1, A2),
      names_to = NULL,
      values_to = "ALLELES"
    ) %>%
    dplyr::distinct(M_SEQ, ALLELES) %>%
    dplyr::count(M_SEQ) %>%
    dplyr::filter(n == 1L) %>%
    dplyr::pull(M_SEQ)

  # Remove the markers from the dataset
  n.markers.removed <- length(bl)
  if (n.markers.removed > 0) {
    # message("Number of monomorphic markers removed: ", n.markers.removed)
    x %<>%  dplyr::filter(!M_SEQ %in% bl)
  }
  bl <- NULL

  # number of marker used for computation
  n.markers <- length(unique(x$M_SEQ))

  # Fst prep: 3 parts
  pop.select <- unique(x$STRATA_SEQ)

  fst.prep <- x %>%
    dplyr::mutate(
      `1` = stringi::stri_sub(GT, 1,3),
      `2` = stringi::stri_sub(GT, 4,6),
      GT = NULL,
      HET = dplyr::if_else(`1` != `2`, 1L, 0L)
    ) %>%
    tidyr::pivot_longer(
      cols = -c("M_SEQ", "STRATA_SEQ", "ID_SEQ", "HET"),
      names_to = "ALLELE_GROUP",
      values_to = "ALLELES"
    ) %>%
    dplyr::mutate(ALLELE_GROUP = as.integer(ALLELE_GROUP))

  # alleles
  fa <- fst.prep %>%
    dplyr::group_by(M_SEQ, STRATA_SEQ, ALLELES) %>%
    dplyr::summarise(
      n = dplyr::n(),
      MHO = length(HET[HET == 1]),
      .groups = "drop"
    ) %>%
    tidyr::complete(
      STRATA_SEQ,
      tidyr::nesting(M_SEQ, ALLELES),
      fill = list(MHO = 0, n = 0)
    ) %>%
    dplyr::group_by(M_SEQ, STRATA_SEQ) %>%
    dplyr::mutate(
      NAPL = sum(n), # Number of alleles per locus
      FREQ_APL = n / NAPL # Frequency of alleles per pop and locus
    ) %>%
    dplyr::ungroup(.)


  # n per locus
  npl <- x %>%
    dplyr::distinct(M_SEQ, STRATA_SEQ) %>%
    dplyr::mutate(NPL = 1L) %>%
    tidyr::pivot_wider(
      names_from = STRATA_SEQ,
      values_from = NPL,
      values_fill = list(NPL = 0L)
    )


  nil <- x %>%
    dplyr::count(M_SEQ, STRATA_SEQ, name = "NIL") %>%
    tidyr::pivot_wider(
      names_from = STRATA_SEQ,
      values_from = NIL,
      values_fill = list(NIL = 0L)
    )

  # using Matrix could speed up... need more test...

  # nil <- x %>%
  #   dplyr::count(MARKERS, STRATA, name = "NIL") %>%
  #   dplyr::mutate(
  #     MARKERS = as.factor(MARKERS),
  #     STRATA = as.factor(STRATA)
  #   )

  # Create the sparse matrix
  # Convert to dense matrix (optional depending on size)
  # Convert to tibble

  # nil  <- Matrix::sparseMatrix(
  #     i = as.integer(nil$MARKERS),
  #     j = as.integer(nil$STRATA),
  #     x = nil$NIL,
  #     dims = c(length(levels(nil$MARKERS)), length(levels(nil$STRATA))),
  #     dimnames = list(levels(nil$MARKERS), levels(nil$STRATA))
  #   ) %>%
  #   as.matrix() %>%
  #   as.data.frame() %>%
  #   tibble::rownames_to_column(var = "MARKERS") %>%
  #   tibble::as_tibble()


  count.locus <- npl %>%
    dplyr::select(M_SEQ) %>%
    dplyr::mutate(
      NPL = rowSums(x = npl[-1], na.rm = TRUE),# much longer using rowwise
      NIL = rowSums(x = nil[-1], na.rm = TRUE),
      NIPL = rowSums(x = nil[-1]^2, na.rm = TRUE), # read nipl square sum
      NC = (NIL - NIPL / NIL) / (NPL - 1)#correction
    )
  npl <- nil <- NULL
  x <- NULL # no longer required

  # integrating markers and alleles (like ref and alt for bi-allelic data)
  fst.prep %<>%
    dplyr::distinct(M_SEQ, ALLELES) %>%
    dplyr::left_join(count.locus, by = "M_SEQ")
  count.locus <- NULL

  # Not so bad part ----
  fst.prep %<>%
    dplyr::full_join(
      fa %<>%
        dplyr::group_by(M_SEQ, ALLELES) %>%
        dplyr::mutate(FREQ_AL = sum(n) / sum(NAPL)) %>%
        dplyr::ungroup(.)
      , by = c("M_SEQ", "ALLELES")
    )
  fa <- NULL

  fst.prep %<>%
    dplyr::mutate(
      NIPL = NAPL / 2,
      MHOM = round(((NAPL * FREQ_APL - MHO) / 2), 0),
      dum = NIPL * (FREQ_APL - 2 * FREQ_APL ^ 2) + MHOM
    ) %>%
    dplyr::group_by(M_SEQ, ALLELES) %>%
    dplyr::mutate(
      SSi = sum(dum, na.rm = TRUE),
      dum1 = NIPL * (FREQ_APL - FREQ_AL) ^ 2,
      SSP = 2 * sum(dum1, na.rm = TRUE)
    ) %>%
    dplyr::group_by(M_SEQ, STRATA_SEQ) %>%
    dplyr::mutate(SSG = NIPL * FREQ_APL - MHOM) %>%
    dplyr::group_by(M_SEQ, ALLELES) %>%
    dplyr::mutate(
      sigw = round(sum(SSG, na.rm = TRUE), 2) / NIL,# ntal -> NIL
      MSP = SSP / (NPL - 1),
      MSI = SSi / (NIL - NPL),
      sigb = 0.5 * (MSI - sigw),
      siga = 0.5 / NC * (MSP - MSI)
    ) %>%
    dplyr::ungroup(.)


  # variance components of allele frequencies for each allele
  # siga: among populations
  # sigb: among individuals within/between populations
  # sigw: within individuals

  # Fast part ------
  fst.prep %<>%
    dplyr::group_by(M_SEQ, ALLELES) %>%
    dplyr::summarise(
      siga = mean(siga, na.rm = TRUE),
      sigb = mean(sigb, na.rm = TRUE),
      sigw = mean(sigw, na.rm = TRUE),
      .groups = "drop"
    )

  # variance components per locus
  # lsiga: among populations
  # lsigb: among individuals within/between populations
  # lsigw: within individuals

  sigma.loc <- fst.prep %>%
    dplyr::group_by(M_SEQ) %>%
    dplyr::summarise(
      lsiga = round(sum(siga, na.rm = TRUE), digits),
      lsigb = round(sum(sigb, na.rm = TRUE), digits),
      lsigw = round(sum(sigw, na.rm = TRUE), digits),
      .groups = "drop"
    )

  fst.fis.markers <- sigma.loc %>%
    dplyr::group_by(M_SEQ) %>%
    dplyr::summarise(
      FST = round(lsiga/(lsiga + lsigb + lsigw), digits),
      FIS = round(lsigb/(lsigb + lsigw), digits),
      .groups = "keep"
    ) %>%
    dplyr::mutate(FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0)) %>%
    dplyr::ungroup(.)

  fst.fis.overall <- fst.prep %>%
    dplyr::summarise(# not big change if using colSums...
      tsiga = sum(siga, na.rm = TRUE),
      tsigb = sum(sigb, na.rm = TRUE),
      tsigw = sum(sigw, na.rm = TRUE)
    ) %>%
    dplyr::summarise(
      FST = round(tsiga / (tsiga + tsigb + tsigw), digits),
      FIS = round(tsigb / (tsigb + tsigw), digits)
    ) %>%
    # automatically remove negative Fst
    dplyr::mutate(FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0))
  # add new column with number of markers
  fst.fis.overall$N_MARKERS <- n.markers

  # Confidence Intervals -----------------------------------------------------
  # over loci for the overall Fst estimate
  # could do something similar for Fis...

  if (ci) {
    # the function:
    boot.fst.list <- purrr::map(
      .x = 1:iteration.ci,
      .f = assigner::boot_ci,
      fst.prep = fst.prep,
      digits = digits
    )
    boot.fst <- dplyr::bind_rows(boot.fst.list)
    boot.fst.summary <- boot.fst %>%
      dplyr::summarise(
        CI_LOW = round(
          stats::quantile(FST, probs = quantiles.ci[1], na.rm = TRUE), digits),
        CI_HIGH = round(
          stats::quantile(FST, probs = quantiles.ci[2], na.rm = TRUE),digits)
      ) %>%
      dplyr::mutate(
        CI_LOW = dplyr::if_else(CI_LOW < 0, true = 0, false = CI_LOW, missing = 0),
        CI_HIGH = dplyr::if_else(CI_HIGH < 0, true = 0, false = CI_HIGH, missing = 0)
      )
  }

  # Fst markers  -------------------------------------------------------------
  fst.markers <- fst.fis.markers %>%
    dplyr::select(M_SEQ, FST) %>%
    dplyr::arrange(M_SEQ)

  # Ranked fst   -------------------------------------------------------------
  n.row.fst <- nrow(fst.markers)
  fst.ranked <- fst.markers %>%
    dplyr::arrange(dplyr::desc(FST)) %>%
    dplyr::select(M_SEQ, FST) %>%
    dplyr::mutate(
      RANKING = seq(from = 1, to = n.row.fst),
      QUARTILE = dplyr::ntile(FST,10)
    )

  # Fst overall  -------------------------------------------------------------
  if (ci) {
    fst.overall <- fst.fis.overall %>%
      dplyr::select(FST, N_MARKERS) %>%
      dplyr::bind_cols(boot.fst.summary)
  } else {
    fst.overall <- fst.fis.overall %>%
      dplyr::select(FST, N_MARKERS)
  }
  fst.overall$N_MONOMORPHIC_BL <- n.markers.removed

  # Fis markers  -------------------------------------------------------------
  fis.markers <- dplyr::select(.data = fst.fis.markers, M_SEQ, FIS) %>%
    dplyr::arrange(M_SEQ)

  # Fis overall   ------------------------------------------------------------
  fis.overall <- dplyr::select(.data = fst.fis.overall, FIS, N_MARKERS)

  # Plot -----------------------------------------------------------------------
  fst.plot <- "Fst plot distribution not selected"
  if (plot.fst) fst.plot <- assigner::plot_fst_distribution(data = fst.markers)

  # ITERATIONS -----------------------------------------------------------------
  # Manage the iterations
  sigma.loc$ITERATIONS <- i
  fst.markers$ITERATIONS <- i
  fst.ranked$ITERATIONS <- i
  fis.markers$ITERATIONS <- i
  fis.overall$ITERATIONS <- i
  fst.overall$ITERATIONS <- i

  # Results --------------------------------------------------------------------
  res <- list(
    sigma.loc = sigma.loc,
    fst.markers = fst.markers,
    fst.ranked = fst.ranked,
    fst.overall = fst.overall,
    fis.markers = fis.markers,
    fis.overall = fis.overall,
    fst.plot = fst.plot
  )
  return(res)
})# End compute_fst function

# pairwise_fst------------------------------------------------------------------
#' @title pairwise_fst
#' @description main function
#' @rdname pairwise_fst
#' @export
#' @keywords internal

pairwise_fst <- carrier::crate(function(
    x,
    pop.pairwise = NULL,
    ci = FALSE,
    iteration.ci = 100,
    quantiles.ci = c(0.025,0.975),
    digits = 9,
    path.folder = NULL,
    ...
) {
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
  `%$%` <- magrittr::`%$%`


  # x <- fst.temp.files[[1]]

  # read the data
  if (inherits(x, "character") && file.exists(x)) {
    # If it's a file path, read the file (assuming it's a Parquet file)
    # message("Reading data from file: ", input_data)
    x <- arrow::read_parquet(x)
  }

  i <- 1L #default
  if (rlang::has_name(x, "ITERATIONS")) {
    i <- unique(x$ITERATIONS)
    x %<>% dplyr::select(-ITERATIONS)
  }

  fst_pair <- function(
    pop.pairwise = NULL,
    x = NULL,
    i = NULL,
    ci = NULL,
    iteration.ci = NULL,
    quantiles.ci = NULL,
    digits = NULL,
    path.folder = NULL
  ) {

    x %<>%
      dplyr::filter(STRATA_SEQ %in% pop.pairwise)# %>%
    # dplyr::mutate(STRATA_SEQ = droplevels(x = STRATA_SEQ))

    # common markers
    # message("Tailoring markers between strata pair")
    common.set <- intersect(
      unique(x$M_SEQ[x$STRATA_SEQ == pop.pairwise[1]]),
      unique(x$M_SEQ[x$STRATA_SEQ == pop.pairwise[2]])
    )

    x %<>%
      dplyr::filter(M_SEQ %in% common.set) %>%
      assigner::compute_fst(
        x = .,
        ci = ci,
        iteration.ci = iteration.ci,
        quantiles.ci = quantiles.ci,
        digits = digits,
        path.folder = path.folder,
        plot = FALSE
      ) %>%
      purrr::pluck("fst.overall") %>%
      dplyr::mutate(
        POP1 = pop.pairwise[1],
        POP2 = pop.pairwise[2],
        .before = 1  # insert at the beginning
      ) %>%
      dplyr::mutate(ITERATIONS = i)
    pop.pairwise <- common.set <- NULL #unused
    return(x)
  }#END inside_pairwise

  res <- purrr::map(
    .x = pop.pairwise,
    .f = function(y) {
      fst_pair(
        pop.pairwise = y,
        x = x,
        i = i,
        ci = ci,
        iteration.ci = iteration.ci,
        quantiles.ci = quantiles.ci,
        digits = digits,
        path.folder = path.folder
      )
    }
  ) %>%
    dplyr::bind_rows()
  return(res)
}) # End pairwise_fst

# pairwise_fst_snprelate--------------------------------------------------------

# @title pairwise_fst_snprelate
# @description Pairwise Fst function with SNPRelate
# @rdname pairwise_fst_snprelate
# @export
# @keywords internal

# pairwise_fst_snprelate <- function(pop.pairwise, data, strata, unique.markers.pop) {
#
#   strata <- dplyr::filter(.data = strata, POP_ID %in% pop.pairwise) %>% # filter the pop
#     dplyr::mutate(POP_ID = droplevels(POP_ID)) # remove unnecessary factors
#
#   # markers in common between pair of pop
#   set1 <- unique.markers.pop %>%
#     dplyr::filter(POP_ID == pop.pairwise[1]) %>%
#     dplyr::select(MARKERS)
#   set2 <- unique.markers.pop %>%
#     dplyr::filter(POP_ID == pop.pairwise[2]) %>%
#     dplyr::select(MARKERS)
#   common.set <- dplyr::intersect(set1, set2) %>%
#     dplyr::arrange(MARKERS)
#
#   # fst.snprelate <- NULL
#   fst.snprelate <- SNPRelate::snpgdsFst(
#     gdsobj = data,
#     population = strata$POP_ID, # factors required
#     sample.id = strata$INDIVIDUALS,
#     snp.id = common.set$MARKERS,
#     method = "W&C84",
#     remove.monosnp = TRUE,
#     maf = NaN,
#     missing.rate = NaN,
#     autosome.only = FALSE,
#     with.id = FALSE,
#     verbose = FALSE
#   )
#   return(fst.snprelate)
# }

# boot_ci-----------------------------------------------------------------------

#' @title boot_ci
#' @description Confidence interval function
#' @rdname boot_ci
#' @export
#' @keywords internal

boot_ci <- carrier::crate(function(x, fst.prep, digits = 9){
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
  `%$%` <- magrittr::`%$%`

  markers.list <- fst.prep %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(M_SEQ) %>%
    dplyr::arrange(M_SEQ)

  subsample.markers <- markers.list %>%
    dplyr::sample_n(tbl = ., size = nrow(markers.list), replace = TRUE) %>%
    dplyr::arrange(M_SEQ)

  markers.list <- NULL

  fst.fis.overall.iterations <- fst.prep %>%
    dplyr::right_join(subsample.markers, by = "M_SEQ") %>%
    dplyr::ungroup(.) %>%
    dplyr::summarise(
      tsiga = sum(siga, na.rm = TRUE),
      tsigb = sum(sigb, na.rm = TRUE),
      tsigw = sum(sigw, na.rm = TRUE)
    ) %>%
    dplyr::summarise(
      FST = round(tsiga/(tsiga + tsigb + tsigw), digits),
      FIS = round(tsigb/(tsigb + tsigw), digits)
    ) %>%
    dplyr::mutate(
      ITERATIONS = rep(x, dplyr::n()),
      FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0)
    )
  return(fst.fis.overall.iterations)
}) # End boot_ci function



# heatmap_fst-------------------------------------------------------------------
#' @title heatmap_fst
#' @description Function that generate an Heatmap of Fst and CI values
#' @param pairwise.fst (object or path). Tibble with all the info generated by
#' \code{\link{fst_WC84}}. If this is not used, \code{pairwise.fst.full.matrix}
#' and \code{pairwise.fst.ci.matrix} are both required.
#' @param pairwise.fst.full.matrix (object or path).
#' @param pairwise.fst.ci.matrix (object or path).
#' @param pop.levels (optional, character, string) If not supplied,
#' the order is set from the colnames of the full fst matrix.
#' Default: \code{pop.levels = NULL}.
#' @param n.s (optional, logical) To have an * when the Fst value is not
#' significant (0 is the lower bound of the CI).
#' Default: \code{n.s = TRUE}.
#' @param digits (optional, integer) The number of digits showed in the heatmap.
#' Default: \code{digits = 5}.
#' @param color.low (optional, character) Color of lower bound.
#' Default: \code{color.low = "blue"}.
#' @param color.mid (optional, character) Mid color value.
#' Default: \code{color.mid = "yellow"}.
#' @param color.high (optional, character) Color of higher bound.
#' Default: \code{color.high = "red"}.
#' @param text.size (optional, integer) Size of the values.
#' Default: \code{text.size = 2}.
#' @param plot.size (optional, integer) By default the size is
#' \code{the number of strata * 2} in cm.
#' Default: \code{plot.size = NULL}.
#' @param path.folder (optional, character)
#' Default: \code{path.folder = NULL}. Default will use the working directory.
#' @param filename (optional, character) Name of the plot to write.
#' Default: \code{filename = NULL}. With default, the plot is not written to disk.
#' @rdname heatmap_fst
#' @export
# @keywords internal
heatmap_fst <- function(
    pairwise.fst = NULL,
    pairwise.fst.full.matrix = NULL,
    pairwise.fst.ci.matrix = NULL,
    pop.levels = NULL,
    n.s = TRUE,
    digits = 5,
    color.low = "blue",
    color.mid = "yellow",
    color.high = "red",
    text.size = 4,
    plot.size = NULL,
    path.folder = NULL,
    filename = NULL
) {

  # ## test
  # pairwise.fst
  # ## pairwise.fst.full.matrix
  # ## pairwise.fst.ci.matrix
  # n.s = TRUE
  # digits = 5
  # color.low = "blue"
  # color.mid = "yellow"
  # color.high = "red"
  # text.size = 4
  # plot.size = 40
  # filename = NULL
  # path.folder = NULL
  # pop.levels = NULL


  # if (
  #   missing(pairwise.fst) ||
  #   missing(pairwise.fst.full.matrix) ||
  #   missing(pairwise.fst.ci.matrix)
  # ) {
  #   rlang::abort("pairwise.fst.full.matrix and/or pairwise.fst.ci.matrix are missing")
  # }
  cli::cli_progress_step(msg = "Generating heatmap")

  pairwise.tib <- FALSE
  if (is.null(pairwise.fst)) {
    if (
      missing(pairwise.fst.full.matrix) ||
      missing(pairwise.fst.ci.matrix)
    ) {
      rlang::abort("pairwise.fst.full.matrix and/or pairwise.fst.ci.matrix are missing")
    }

    if (is.vector(pairwise.fst.full.matrix)) {
      data.fst <- readRDS(pairwise.fst.full.matrix)
    } else {
      data.fst <- pairwise.fst.full.matrix
    }

    if (is.vector(pairwise.fst.ci.matrix)) {
      data.ci <- readRDS(pairwise.fst.ci.matrix)
    } else {
      data.ci <- pairwise.fst.ci.matrix
    }
    if (is.null(pop.levels)) {
      pop.levels <- colnames(data.fst)
    }

  } else {
    if (missing(pairwise.fst.full.matrix)) {
      rlang::abort("pairwise.fst.full.matrix and/or pairwise.fst.ci.matrix are missing")
    }
    pairwise.tib <- TRUE
    if (is.vector(pairwise.fst)) {
      data.fst <- readr::read_tsv(file = pairwise.fst)
    } else {
      data.fst <- pairwise.fst
    }
    if (is.null(pop.levels)) {
      if (is.factor(data.fst$POP1)) {
        pop.levels <- levels(data.fst$POP1)
      } else {
        pop.levels <- as.character(unique(c(data.fst$POP1, data.fst$POP2)))
      }
    }
    data.fst %<>%
      dplyr::select(POP1, POP2, FST, CI_LOW, CI_HIGH, FST_HEATMAP) %>%
      tidyr::unite(data = ., col = CI, c("CI_LOW", "CI_HIGH"), sep = "\n")

    data.fst <- dplyr::bind_rows(
      data.fst %>% dplyr::select(POP1, POP2, FST, CI),
      data.fst %>% dplyr::select(POP1 = POP2, POP2 = POP1, FST, CI = FST_HEATMAP)
    ) %>%
      dplyr::mutate(FST = as.numeric(FST))
  }

  inv.levels <- rev(pop.levels)

  if (!pairwise.tib) {
    data.fst %<>%
      radiator::distance2tibble(
        x = .,
        remove.diag = FALSE,
        na.diag = TRUE,
        remove.lower = FALSE,
        relative = FALSE,
        pop.levels = pop.levels
      ) %>%
      magrittr::set_colnames(x = ., value = c("POP1", "POP2", "FST"))

    rounder <- function(x, digits) round(as.numeric(x), digits)
    if (max(stringi::stri_length(data.fst$FST), na.rm = TRUE) - 1 != digits) {
      round.num <- TRUE
    } else {
      round.num <- FALSE
    }

    if (n.s) round.num <- TRUE

    data.ci %<>%
      radiator::distance2tibble(
        x = .,
        remove.diag = FALSE,
        na.diag = TRUE,
        remove.lower = FALSE,
        relative = FALSE,
        distance.class.double = FALSE,
        pop.levels = pop.levels
      ) %>%
      magrittr::set_colnames(x = ., value = c("POP1", "POP2", "CI"))

    # test1 <- data.fst %>% dplyr::left_join(data.ci, by = c("POP1", "POP2"))

    data.fst %<>% dplyr::left_join(data.ci, by = c("POP1", "POP2"))

    # Verify no negative ...
    # data without CI for the lower diag
    if (n.s || round.num) {
      data.ci <- dplyr::filter(data.fst, stringi::stri_detect_fixed(str = CI, pattern = " - ")) %>%
        tidyr::separate(data = ., col = CI, into = c("LOW", "HIGH"), sep = " - ") %>%
        dplyr::mutate(dplyr::across(.cols = c("FST", "LOW", "HIGH"), .fns = rounder, digits = digits)) %>%
        dplyr::mutate(
          dplyr::across(
            .cols = c("FST", "LOW", "HIGH"),
            .fns = ~ pmax(.x, 0)# remove negative values...
          )
        ) %>%
        dplyr::mutate(NS = dplyr::if_else(LOW == 0, TRUE, FALSE)) %>%
        dplyr::mutate(dplyr::across(.cols = c("LOW", "HIGH"), .fns = format, scientific = FALSE)) %>%
        tidyr::unite(data = ., col = CI, c("LOW", "HIGH"), sep = "\n")
      ns <- dplyr::distinct(data.ci, POP1, POP2, NS) %>%
        dplyr::rename(POP3 = POP2, POP2 = POP1) %>%
        dplyr::rename(POP1 = POP3) %>%
        dplyr::mutate(POP1 = as.character(POP1), POP2 = as.character(POP2))
      data.ci %<>% dplyr::select(-NS)
      # data.fst
      suppressWarnings(
        data.fst %<>%
          dplyr::filter(!stringi::stri_detect_fixed(str = CI, pattern = " - ") | is.na(CI)) %>%
          dplyr::mutate(dplyr::across(.cols = c("FST", "CI"), .fns = rounder, digits = digits)) %>%
          dplyr::mutate(dplyr::across(.cols = "CI", .fns = format, scientific = FALSE)) %>%
          dplyr::left_join(ns, by = c("POP1", "POP2"))
      )

      if (n.s) {
        suppressWarnings(
          data.fst %<>%
            dplyr::mutate(CI = dplyr::if_else(NS, paste0(CI, "*"), CI), NS = NULL)
        )
      }
      suppressWarnings(data.fst %<>% dplyr::bind_rows(data.ci))
      data.ci <- ns <- NULL
    }
  }

  # Factors
  data.fst$POP2 <- factor(x = as.character(data.fst$POP2), levels = inv.levels, ordered = TRUE)
  data.fst$POP1 <- factor(x = as.character(data.fst$POP1), levels = pop.levels, ordered = TRUE)

  heatmap.fst <- ggplot2::ggplot(
    data = data.fst, ggplot2::aes(x = POP1, y = POP2, fill = FST)
  ) +
    ggplot2::geom_tile(color = "white", alpha = 0.7) +
    ggplot2::geom_text(ggplot2::aes(x = POP1, y = POP2, label = CI),
                       color = "black", size = text.size, na.rm = TRUE) +
    ggplot2::scale_fill_viridis_c(name = "Fst", option = "H") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    ) +
    ggplot2::coord_equal()

  print(heatmap.fst)

  if (!is.null(filename)) {
    if (is.null(path.folder)) path.folder <- getwd()

    if (digits > 5) {
      mult.fact <- 3
    } else {
      mult.fact <- 2
    }
    if (is.null(plot.size)) plot.size <- max(20, length(levels(data.fst$POP1)) * mult.fact)

    heatmap.pdf <- stringi::stri_join(filename, "_heatmap.fst.pdf")
    heatmap.png <- stringi::stri_join(filename, "_heatmap.fst.png")

    ggplot2::ggsave(
      filename = file.path(path.folder, heatmap.png),
      plot = heatmap.fst,
      width = plot.size,
      height = plot.size,
      dpi = 200,
      units = "cm",
      device = "png",
      limitsize = FALSE)

    ggplot2::ggsave(
      filename = file.path(path.folder, heatmap.pdf),
      plot = heatmap.fst,
      width = plot.size,
      height = plot.size,
      dpi = 300,
      units = "cm",
      device = "pdf",
      limitsize = FALSE,
      useDingbats = FALSE)

  }
  return(heatmap.fst)
}#End heatmap_fst

# fst_stats-----------------------------------------------------------
#' @title fst_stats
#' @description Generate useful stats
#' @rdname fst_stats
#' @export
#' @keywords internal

fst_stats <- function(
    x,
    l,
    digits = 9L,
    m = NULL,
    s = NULL,
    iteration.subsample = NULL
    ) {
  # x = "pairwise.fst" # test
  # x = "sigma.loc"
  # x = "fst.markers"
  # x = "fst.ranked"
  # x = "fst.overall"
  # x = "fis.markers"
  # x = "fis.overall"

  subsampling <- iteration.subsample > 1

  res.stats <- list()
  want <- c("sigma.loc", "fst.markers", "fst.ranked", "fst.overall",
            "fis.markers", "fis.overall", "pairwise.fst")

  #1. flatten the tibble
  # res1 <- purrr::pluck(l, 1, x)

  res1 <- purrr::map(l, x)

  if (x %in% want) {
    res1 %<>% dplyr::bind_rows(.)
    # res1 <- purrr::pluck(l, 1, x)
  }

  #2. generate the stats
  want <- c("fst.markers","fis.markers", "fst.overall", "fis.overall")
  if (x %in% want && subsampling) {
    # message("Stats: ", x)

    # defaults:
    group.by <- NULL
    overall <- TRUE
    v <- "FIS"

    # mod
    if (x %in% c("fst.markers","fis.markers")) {
      group.by <- "M_SEQ"
      overall <- FALSE
    }

    if (x %in% c("fst.markers","fst.overall")) v <- "FST"

    if (!is.null(group.by)) {
      res1 %<>% dplyr::group_by(.data[[group.by]])
    }

    if (overall) {
      n.markers.mean <- res1 %>%
        dplyr::summarise(N_MARKERS_MEAN = mean(N_MARKERS))
    }

    res1 %<>%
      dplyr::summarise(
        MEAN = mean(.data[[v]], na.rm = TRUE),
        SE = sqrt(stats::var(.data[[v]], na.rm = TRUE) / length(.data[[v]])),
        SD = stats::sd(.data[[v]], na.rm = TRUE),
        MEDIAN = stats::median(.data[[v]], na.rm = TRUE),
        Q25 = stats::quantile(.data[[v]], 0.25, na.rm = TRUE),
        Q75 = stats::quantile(.data[[v]], 0.75, na.rm = TRUE),
        IQR = stats::IQR(.data[[v]], na.rm = TRUE),
        MIN = min(.data[[v]], na.rm = TRUE),
        MAX = max(.data[[v]], na.rm = TRUE),
        ITERATIONS = dplyr::n(),
        .groups = "drop"
      )

    if (overall) {
      res1$N_MARKERS_MEAN <- n.markers.mean$N_MARKERS_MEAN
    }

    res1 %<>%
      dplyr::mutate(
        dplyr::across(
          .cols = where(is.numeric),
          .fns = round,
          digits = digits
        )
      )

  }

  #3. adjust some objects
  # pairwise.fst
  if (x == "pairwise.fst") {
    res1 %<>%
      dplyr::mutate(dplyr::across(c(POP1, POP2), as.integer)) %>%
      dplyr::mutate(
        dplyr::across(
          .cols = c("FST", "CI_LOW", "CI_HIGH"),
          .fns = ~ pmax(.x, 0)
        )
      )

    if (subsampling) {
      fmt.args <- list(format = "f", digits = digits, decimal.mark = ".")

      res1 %<>%
        dplyr::group_by(POP1, POP2) %>%
        dplyr::summarise(
          FST_RANGE = format_mean_range(FST, meanx = FALSE, digits = digits),
          FST = format_mean_range(FST, min.max = FALSE, digits = digits),
          N_MARKERS = format_mean_range(N_MARKERS, digits = 0),
          N_MONOMORPHIC_BL = format_mean_range(N_MONOMORPHIC_BL, digits = 0),
          CI_LOW = utils_formatC(min(CI_LOW, na.rm = TRUE), fmt.args = fmt.args),
          CI_HIGH = utils_formatC(max(CI_HIGH, na.rm = TRUE), fmt.args = fmt.args),
          .groups = "drop"
        )
    }

    res1 %<>%
      dplyr::mutate(
        NS = dplyr::if_else(as.numeric(CI_LOW) == 0, TRUE, FALSE),
        FST_HEATMAP = dplyr::if_else(NS, paste0(FST, "*"), FST),
        ITERATIONS = rep(iteration.subsample, dplyr::n())
      ) %>%
      dplyr::left_join(s, by = c("POP1" = "STRATA_SEQ")) %>%
      dplyr::select(-POP1, POP1 = STRATA) %>%
      dplyr::left_join(s, by = c("POP2" = "STRATA_SEQ")) %>%
      dplyr::select(-POP2, POP2 = STRATA) %>%
      dplyr::select(
        tidyselect::any_of(
          c("POP1", "POP2", "FST", "FST_RANGE",
            "CI_LOW", "CI_HIGH",  "N_MARKERS", "N_MONOMORPHIC_BL",
            "NS", "FST_HEATMAP", "ITERATIONS")
          )
        )
  }


  #4. Change strata and markers
  want <- c("sigma.loc", "fst.markers", "fst.ranked", "fis.markers")
  if (x %in% want) {
    # iterations.string <- res1 %>% dplyr::count(M_SEQ)
    res1 %<>% match_markers_meta(x = ., markers.meta = m)
  }

  # 5. fst.plot
  if (x %in% "fst.plot" && !subsampling) res1 <- res1[[1]]

  res.stats[[x]] <- res1
  return(res.stats)
}#End fst_stats

# Minimize arguments list-------------------------------------------------------
#' @title Minimize arguments list for parallel execution
#' @description Removes unused or large objects from a list of arguments.
#' Designed to reduce memory overhead when passing arguments to parallel workers.
#' @param args A named list of arguments (e.g., args.for.fst)
#' @param keep A character vector of argument names to keep
#' @param verbose Show what is being dropped (default: TRUE)
#' @return A slimmed-down list with only relevant arguments
#' @keywords internal
#' @export
minimize_args_list <- function(args, keep, verbose = TRUE) {
  if (!is.list(args) || is.null(names(args))) {
    rlang::abort("Input `args` must be a named list.")
  }

  drop <- setdiff(names(args), keep)

  if (length(drop) > 0 && verbose) {
    message("Dropping unused arguments from args list:\n    ",
            stringi::stri_paste(drop, collapse = "\n    "))
  }

  args[keep]
}





