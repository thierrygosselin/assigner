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
#' There is no null hypothesis testing with \emph{P}-values.
#' Confidence intervals provided with the \emph{F}-statistics
#' enables more reliable conclusions about the biological trends in the data.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{radiator} \code{\link{tidy_genomic_data}}.
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

#' @param digits (optional, integer) The number of decimal places to be used in
#' results.
#' Default: \code{digits = 9}.

#' @param parallel.core (optional, integer) The number of core for parallel computation
#' of pairwise Fst.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.


#' @param filename (optional, character) Give filename prefix, this will trigger
#' saving results in a directory.
#' Default: \code{filename = NULL}.

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
#' \item \code{heatmap.fst} to generate an heatmap with the Fst values in
#' lower matrix and CI in the upper matrix.
#' Default: \code{heatmap.fst = FALSE}.
#' The heatmap can also be generated and more finetuned separately after the Fst
#' analysis using \code{\link{heatmap_fst}}.
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
  digits = 9,
  filename = NULL,
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
  # filter.monomorphic=TRUE
  # calibrate.alleles = FALSE

  # Cleanup---------------------------------------------------------------------
  assigner_function_header(f.name = "fst_WC84", verbose = verbose)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date/time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- assigner_tic()
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(assigner_toc(timing), add = TRUE)
  on.exit(assigner_function_header(f.name = "fst_WC84", start = FALSE, verbose = verbose), add = TRUE)

  # Function call and dotslist -------------------------------------------------
  rad.dots <- assigner::assigner_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    keepers = c("filter.monomorphic", "holdout.samples", "subsample",
                "iteration.subsample", "heatmap.fst", "blacklist.id",
                "calibrate.alleles"),
    verbose = FALSE
  )
  dots.filename <- stringi::stri_join("assigner_fst_WC84_args_", file.date, ".tsv")
  # currently not saved

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("data is missing")
  if (!ci && heatmap.fst) {
    heatmap.fst <- FALSE
    if (verbose) message("\nconfidence intervals not selected, heatmap.fst: FALSE\n")
  }
  if (!filter.monomorphic) {
    message("filter.monomorphic = FALSE... not a good idea, but lets do it...")
  }
  # filename & folder ----------------------------------------------------------
  path.folder <- NULL
  if (!is.null(filename)) {
    filename <- stringi::stri_join(filename, "_fst_WC84")
    path.folder <- radiator::generate_folder(
      f = filename,
      file.date = file.date,
      verbose = verbose)
  }

  if (snprelate) {
    # Check that snprelate is installed
    if (!"SNPRelate" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install SNPRelate for this option:\n
                 install.packages("BiocManager")
                 BiocManager::install("SNPRelate")')
    }
    rlang::abort("Until the bias observed with SNPRelate is resolved, the option is unavailable.")
  }

  # Import data ---------------------------------------------------------------
  if (verbose) message("Importing data")
  data %<>% radiator::tidy_wide(data = ., import.metadata = TRUE)

  if (!rlang::has_name(data, "GT") || calibrate.alleles) {
    data %<>%
      radiator::calibrate_alleles(data = ., parallel.core = parallel.core) %$%
      input
  }

  # Strata----------------------------------------------------------------------
  strata <- radiator::read_strata(
    strata = strata,
    pop.id = TRUE,
    blacklist.id = blacklist.id,
    pop.levels = NULL,
    verbose = verbose) %$%
    strata

  # population levels and strata------------------------------------------------
  if (!is.null(strata)) {
    data <- radiator::join_strata(
      data = data, strata = strata, pop.id = TRUE, verbose = FALSE)
  }

  if (!rlang::has_name(data, "POP_ID") && rlang::has_name(data, "STRATA")) {
    data %<>% dplyr::rename(POP_ID = STRATA)
  }

  pop.levels.bk <- pop.levels
  if (is.null(pop.levels)) pop.levels.bk <- unique(data$POP_ID)

  data %<>%
    dplyr::mutate(POP_ID = factor(x = POP_ID, levels = pop.levels.bk)) %>%
    dplyr::arrange(POP_ID)

  # strip the data -------------------------------------------------------------
  strata.bk <- markers.meta.bk <- genotypes.meta.bk <- NULL
  env.arg <- rlang::current_env()
  data %<>%
    radiator::strip_rad(
      x = .,
      m = c("VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"),
      env.arg = env.arg,
      keep.strata = TRUE,
      verbose = FALSE
    ) %>%
    dplyr::select(tidyselect::any_of(c("M_SEQ", "STRATA_SEQ", "ID_SEQ", "GT"))) %>%
    dplyr::rename(MARKERS = M_SEQ, STRATA = STRATA_SEQ, INDIVIDUALS = ID_SEQ)

  pop.levels <- unique(data$STRATA)

  # subsampling data------------------------------------------------------------
  # create the subsampling list
  if (!is.null(subsample) && !is.numeric(subsample)) {
    heatmap.fst <- FALSE
    if (subsample == "min") {
      subsample <- strata.bk %>%
        dplyr::group_by(STRATA_SEQ) %>%
        dplyr::tally(.) %>%
        dplyr::filter(n == min(n)) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(n) %>%
        purrr::flatten_int(.)
    }
  }


  subsample.list <- purrr::map(
    .x = 1:iteration.subsample,
    .f = subsampling_data,
    strata = strata.bk,
    subsample = subsample
  )

  # keep track of subsampling individuals and write to directory
  if (!is.null(subsample)) {
    if (verbose) message("Subsampling: selected")
    res$subsample$subsampling.individuals <- subsample.list %>%
      dplyr::bind_rows() %>%
      readr::write_tsv(
        x = .,
        file = file.path(path.folder, "subsampling.individuals.tsv")
      )
  } # End subsampling

  # Calculations ----------------------------------------------------------------
  subsample.fst <- purrr::map(
    .x = subsample.list,
    .f = fst_subsample,
    data = data,
    snprelate = snprelate,
    strata = strata.bk,
    holdout.samples = holdout.samples,
    pairwise = pairwise,
    ci = ci,
    iteration.ci = iteration.ci,
    quantiles.ci = quantiles.ci,
    digits = digits,
    subsample = subsample,
    path.folder = path.folder,
    parallel.core = parallel.core,
    verbose = verbose
  )
  subsample.list <- NULL

  # Compiling results-----------------------------------------------------------
  if (verbose) message("Generating statistics...")
  # no subsampling --------------------------
  if (is.null(subsample)) {

    # These are the objects:
    # sigma.loc
    # fst.markers
    # fst.ranked
    # fst.overall
    # fis.markers
    # fis.overall
    # fst.plot
    # pairwise.fst & pairwise.fst.mean
    # pairwise.fst.upper.matrix & pairwise.fst.upper.matrix.mean
    # pairwise.fst.full.matrix & pairwise.fst.full.matrix.mean
    # merge upper and lower matrix
    # pairwise.fst.ci.matrix

    # subsample.fst <- purrr::flatten(subsample.fst)
    # res <- purrr::prepend(x = res, values = purrr::flatten(subsample.fst))
    # change strata --------
    nms <- subsample.fst %>% purrr::map(names) %>% purrr::reduce(union)
    res <- purrr::map(
      .x = nms,
      .f = fst_stats,
      l = subsample.fst,
      digits = digits,
      m = markers.meta.bk,
      s = dplyr::distinct(strata.bk, POP_ID, STRATA_SEQ),
      subsample = FALSE
    ) %>%
      purrr::flatten(.)

    # test1 <- res$pairwise.fst
    # test2 <- res$pairwise.fst.upper.matrix


    # pairwise.fst.upper.matrix
    # pairwise.fst.full.matrix
    # pairwise.fst.ci.matrix
    # write results --------
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
  }# end of compiling results NO SUBSAMPLE

  # Compile subsampling results --------------
  if (!is.null(subsample)) {
    nms <- subsample.fst %>% purrr::map(names) %>% purrr::reduce(union)
    res$subsample <- purrr::map(
      .x = nms,
      fst_stats,
      l = subsample.fst,
      digits = digits,
      m = markers.meta.bk,
      s = dplyr::distinct(strata.bk, POP_ID, STRATA_SEQ),
      subsample = TRUE
      ) %>%
      purrr::flatten(.)

    # These are the objects:
    # sigma.loc
    # fst.markers
    # fst.ranked
    # fst.overall
    # fis.markers
    # fis.overall
    # fst.plot
    # pairwise.fst & pairwise.fst.mean
    # pairwise.fst.upper.matrix & pairwise.fst.upper.matrix.mean
    # pairwise.fst.full.matrix & pairwise.fst.full.matrix.mean
    # merge upper and lower matrix
    # pairwise.fst.ci.matrix

    # test1 <- res$subsample$sigma.loc
    # test2 <- res$subsample$pairwise.fst
    # test3 <- res$subsample$pairwise.fst.full.matrix
    # test3[3]

    # Work on the matrix of FST-------
    res$subsample$pairwise.fst.upper.matrix.mean <- res$subsample$pairwise.fst %>%
      dplyr::select(POP1, POP2, FST) %>%
      tidyr::complete(data = ., POP1, POP2) %>%
      assigner::rad_wide(x = ., formula = "POP1 ~ POP2", values_from = "FST", values_fill = "") %>%
      dplyr::rename(POP = POP1)
    rn <- res$subsample$pairwise.fst.upper.matrix.mean$POP # rownames
    res$subsample$pairwise.fst.upper.matrix.mean <- as.matrix(res$subsample$pairwise.fst.upper.matrix.mean[,-1])# make matrix without first column
    rownames(res$subsample$pairwise.fst.upper.matrix.mean) <- rn

    # pairwise.fst.full.matrix & pairwise.fst.full.matrix.mean
    res$subsample$pairwise.fst.full.matrix.mean <- res$subsample$pairwise.fst.upper.matrix.mean # bk of upper.mat.fst
    lower.mat.fst <- t(res$subsample$pairwise.fst.full.matrix.mean) # transpose

    # merge upper and lower matrix
    res$subsample$pairwise.fst.full.matrix.mean[lower.tri(res$subsample$pairwise.fst.full.matrix.mean)] <- lower.mat.fst[lower.tri(lower.mat.fst)]
    diag(res$subsample$pairwise.fst.full.matrix.mean) <- "0"


    # write results --------

    if (!is.null(filename)) {
      purrr::walk(
        .x = list("sigma.loc", "fst.markers", "fst.ranked", "fst.overall",
                  "fis.markers", "fis.overall", "pairwise.fst"),
        .f = fst_write, list.sub = res$subsample, path.folder = path.folder
      )
      saveRDS(
        object = res$subsample$pairwise.fst.upper.matrix.mean,
        file = file.path(path.folder, "pairwise.fst.upper.matrix.RData")
      )
      saveRDS(
        object = res$subsample$pairwise.fst.full.matrix.mean,
        file = file.path(path.folder, "pairwise.fst.full.matrix.RData"))
    }

    # CI ----------
    # defaults
    res$subsample$pairwise.fst.ci.matrix <-
      res$subsample$pairwise.fst.ci.matrix.mean <-
      "confidence intervals not selected"

    if (ci) {
      # pairwise.fst.ci.matrix
      # pairwise.fst.ci.matrix.mean
      lower.mat.ci.sub <- res$subsample$pairwise.fst %>%
        dplyr::select(POP1, POP2, CI_LOW, CI_HIGH) %>%
        tidyr::unite(data = ., CI, CI_LOW, CI_HIGH, sep = " - ") %>%
        tidyr::complete(data = ., POP1, POP2) %>%
        assigner::rad_wide(x = ., formula = "POP1 ~ POP2", values_from = "CI", values_fill = "") %>%
        dplyr::rename(POP = POP1)

      cn <- colnames(lower.mat.ci.sub) # bk of colnames
      lower.mat.ci.sub <- t(lower.mat.ci.sub[,-1]) # transpose
      colnames(lower.mat.ci.sub) <- cn[-1] # colnames - POP
      lower.mat.ci.sub = as.matrix(lower.mat.ci.sub) # matrix

      # merge upper and lower matrix
      pairwise.fst.ci.matrix.sub <- res$subsample$pairwise.fst.upper.matrix.mean # bk upper.mat.fst
      pairwise.fst.ci.matrix.sub[lower.tri(pairwise.fst.ci.matrix.sub)] <- lower.mat.ci.sub[lower.tri(lower.mat.ci.sub)]
      res$subsample$pairwise.fst.ci.matrix.mean <- pairwise.fst.ci.matrix.sub
      pairwise.fst.ci.matrix.sub <- NULL
      if (!is.null(filename)) {
        saveRDS(
          object = res$subsample$pairwise.fst.ci.matrix.mean,
          file = file.path(path.folder, "pairwise.fst.ci.matrix.RData")
          )
      }
    }
  } # end of compiling subsample results

  # heatmap.fst ----------------------------------------------------------------
  if (heatmap.fst) {
    if (verbose) message("Generating heatmap...")
    res$heatmap.fst <- heatmap_fst(
      pairwise.fst.full.matrix = res$pairwise.fst.full.matrix,
      pairwise.fst.ci.matrix = res$pairwise.fst.ci.matrix,
      digits = digits,
      path.folder = path.folder,
      filename = filename)
  }

  # End -------------------------------------------------------------------
  if (verbose) {
    cat("################################### RESULTS ####################################\n")
    if (is.null(subsample)) {
      if (ci) {
        message("Fst (overall): ", res$fst.overall$FST, " [", res$fst.overall$CI_LOW, " - ", res$fst.overall$CI_HIGH, "]")
      } else{
        message("Fst (overall): ", res$fst.overall$FST)
      }
    } else {
      message("Fst (overall): ", res$subsample$fst.overall$MEAN)
    }
  }
  return(res)
}

# Internal Nested Functions to compute WC84 Fst --------------------------------

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
  parallel.core = parallel::detectCores(.) - 2
) {
  # TEST
  # ci = FALSE
  # iteration.ci = 100
  # quantiles.ci = c(0.025,0.975)
  # digits = 9
  # path.folder = NULL
  ## x = data
  # x <- dplyr::filter(data, GT != "000000")
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
  `%$%` <- magrittr::`%$%`

  # Removing monomorphic markers
  x %<>%
    radiator::filter_monomorphic(data = ., internal = TRUE, path.folder = path.folder)

  # number of marker used for computation
  n.markers <- length(unique(x$MARKERS))
  count.locus <- dplyr::bind_cols(
    dplyr::distinct(x, MARKERS, STRATA) %>% dplyr::count(MARKERS, name = "NPL"),# number of populations per locus
    # dplyr::distinct(x, MARKERS, INDIVIDUALS) %>% dplyr::count(MARKERS, name = "NIL") %>% dplyr::select(-MARKERS)# number of individuals per locus
    dplyr::distinct(x, MARKERS1 = MARKERS, INDIVIDUALS) %>% dplyr::count(MARKERS1, name = "NIL")# number of individuals per locus
  )

  if (!identical(count.locus$MARKERS1, count.locus$MARKERS)) {
    rlang::abort("Problem...contact author")
  } else {
    count.locus %<>% dplyr::select(-MARKERS1)
  }

  # faster:
  count.locus %<>%
    dplyr::bind_cols(
      dplyr::count(x, STRATA, MARKERS, name = "NIPL") %>%
        dplyr::mutate(NIPL_SQ = NIPL ^ 2) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::summarise(NIPL_SQ_SUM = sum(NIPL_SQ, na.rm = TRUE), .groups = "drop") %>%
        dplyr::rename(MARKERS1 = MARKERS)
    )

  if (!identical(count.locus$MARKERS1, count.locus$MARKERS)) {
    rlang::abort("contact author")
  } else {
    count.locus %<>%
      dplyr::select(-MARKERS1) %>%
      dplyr::mutate(NC = (NIL - NIPL_SQ_SUM / NIL) / (NPL - 1))#correction
  }

  # prep data: allele per locus and frequency of alleles
  fst.prep <- x %>%
    dplyr::mutate(
      `1` = stringi::stri_sub(GT, 1,3),
      `2` = stringi::stri_sub(GT, 4,6),
      GT = NULL,
      HET = dplyr::if_else(`1` != `2`, 1L, 0L)
    ) %>%
    assigner::rad_long(
      x = .,
      cols = c("MARKERS", "STRATA", "INDIVIDUALS", "HET"),
      names_to = "ALLELE_GROUP",
      values_to = "ALLELES",
      variable_factor = TRUE
    ) %>%
    dplyr::mutate(ALLELE_GROUP = as.integer(ALLELE_GROUP))

  count.locus %<>%
    dplyr::full_join(dplyr::distinct(fst.prep, MARKERS, ALLELES), by = "MARKERS")

  fa <- dplyr::count(x = fst.prep, MARKERS, ALLELES, STRATA) %>%
    tidyr::complete(data = ., STRATA, tidyr::nesting(MARKERS, ALLELES), fill = list(n = 0)) %>%
    dplyr::group_by(MARKERS, STRATA) %>%
    dplyr::mutate(
      NAPL = sum(n),
      FREQ_APL = n / NAPL # Frequency of alleles per pop and locus
    ) %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::mutate(FREQ_AL = sum(n) / sum(NAPL)) %>% #Frequency of alleles per locus
    dplyr::full_join(count.locus, by = c("MARKERS", "ALLELES")) %>%
    dplyr::arrange(MARKERS, STRATA)

  count.locus <- NULL

  fst.prep %<>%
    dplyr::select(-ALLELE_GROUP) %>%
    dplyr::group_by(MARKERS, STRATA, ALLELES) %>%
    dplyr::summarise(MHO = length(HET[HET == 1]), .groups = "drop") %>%
    tidyr::complete(data = ., STRATA, tidyr::nesting(MARKERS, ALLELES), fill = list(MHO = 0)) %>%
    dplyr::arrange(MARKERS, ALLELES, STRATA) %>%
    dplyr::full_join(fa, by = c("STRATA", "MARKERS", "ALLELES")) %>%
    dplyr::mutate(
      NIPL = NAPL / 2,
      MHOM = round(((NAPL * FREQ_APL - MHO) / 2), 0),
      dum = NIPL * (FREQ_APL - 2 * FREQ_APL ^ 2) + MHOM
    ) %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::mutate(
      SSi = sum(dum, na.rm = TRUE),
      dum1 = NIPL * (FREQ_APL - FREQ_AL) ^ 2,
      SSP = 2 * sum(dum1, na.rm = TRUE)
    ) %>%
    dplyr::group_by(MARKERS, STRATA) %>%
    dplyr::mutate(SSG = NIPL * FREQ_APL - MHOM) %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::mutate(
      sigw = round(sum(SSG, na.rm = TRUE), 2) / NIL,# ntal -> NIL
      MSP = SSP / (NPL - 1),
      MSI = SSi / (NIL - NPL),
      sigb = 0.5 * (MSI - sigw),
      siga = 0.5 / NC * (MSP - MSI)
    ) %>%
    dplyr::ungroup(.)
  fa <- NULL

  # variance components of allele frequencies for each allele
  # siga: among populations
  # sigb: among individuals within/between populations
  # sigw: within individuals
  sigma.loc.alleles <- fst.prep %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::summarise(
      siga = mean(siga, na.rm = TRUE),
      sigb = mean(sigb, na.rm = TRUE),
      sigw = mean(sigw, na.rm = TRUE),
      .groups = "drop"
    )
  fst.prep <- NULL

  # variance components per locus
  # lsiga: among populations
  # lsigb: among individuals within/between populations
  # lsigw: within individuals

  sigma.loc <- sigma.loc.alleles %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::summarise(
      lsiga = round(sum(siga, na.rm = TRUE), digits),
      lsigb = round(sum(sigb, na.rm = TRUE), digits),
      lsigw = round(sum(sigw, na.rm = TRUE), digits),
      .groups = "drop"
    )

  fst.fis.markers <- sigma.loc %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::summarise(
      FST = round(lsiga/(lsiga + lsigb + lsigw), digits),
      FIS = round(lsigb/(lsigb + lsigw), digits),
      .groups = "keep"
    ) %>%
    dplyr::mutate(FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0)) %>%
    dplyr::ungroup(.)

  fst.fis.overall <- sigma.loc.alleles %>%
    dplyr::summarise(
      tsiga = sum(siga, na.rm = TRUE),
      tsigb = sum(sigb, na.rm = TRUE),
      tsigw = sum(sigw, na.rm = TRUE)
    ) %>%
    dplyr::summarise(
      FST = round(tsiga / (tsiga + tsigb + tsigw), digits),
      FIS = round(tsigb / (tsigb + tsigw), digits)
    ) %>%
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
      sigma.loc.alleles = sigma.loc.alleles,
      digits = digits
    )
    boot.fst <- dplyr::bind_rows(boot.fst.list)
    boot.fst.summary <- boot.fst %>%
      dplyr::summarise(
        CI_LOW = round(
          stats::quantile(FST, probs = quantiles.ci[1], na.rm = TRUE), digits),
        CI_HIGH = round(
          stats::quantile(FST, probs = quantiles.ci[2], na.rm = TRUE),digits)
      )
  }

  # Fst markers  -------------------------------------------------------------
  fst.markers <- fst.fis.markers %>%
    dplyr::select(MARKERS, FST) %>%
    dplyr::arrange(MARKERS)

  # Ranked fst   -------------------------------------------------------------
  fst.ranked <- fst.markers %>%
    dplyr::arrange(dplyr::desc(FST)) %>%
    dplyr::select(MARKERS, FST) %>%
    dplyr::mutate(
      RANKING = seq(from = 1, to = dplyr::n()),
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

  # Fis markers  -------------------------------------------------------------
  fis.markers <- dplyr::select(.data = fst.fis.markers, MARKERS, FIS) %>%
    dplyr::arrange(MARKERS)

  # Fis overall   ------------------------------------------------------------
  fis.overall <- dplyr::select(.data = fst.fis.overall, FIS, N_MARKERS)

  # Plot -----------------------------------------------------------------------
  fst.plot <- ggplot2::ggplot(fst.markers, ggplot2::aes(x = FST, na.rm = TRUE)) +
    ggplot2::geom_histogram(binwidth = 0.01) +
    ggplot2::labs(x = "Fst (overall)") +
    ggplot2::expand_limits(x = 0) +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.title = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      legend.text = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      strip.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"))

  # Results ------------------------------------------------------------------
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
}) # End compute_fst function

# pairwise_fst------------------------------------------------------------------
#' @title pairwise_fst
#' @description Pairwise Fst function
#' @rdname pairwise_fst
#' @export
#' @keywords internal

pairwise_fst <- carrier::crate(function(
  list.pair,
  pop.pairwise = NULL,
  data = NULL,
  ci = FALSE,
  iteration.ci = 100,
  quantiles.ci = c(0.025,0.975),
  digits = 9,
  path.folder = path.folder
) {
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
  `%$%` <- magrittr::`%$%`

  pop.select <- stringi::stri_join(purrr::flatten(pop.pairwise[list.pair]))
  data.select <- data %>%
    dplyr::filter(STRATA %in% pop.select) %>%
    dplyr::mutate(STRATA = droplevels(x = STRATA))

  # common markers
  common.set <- intersect(
    unique(data.select$MARKERS[data.select$STRATA == pop.select[1]]),
    unique(data.select$MARKERS[data.select$STRATA == pop.select[2]])
  )

  data.select %<>% dplyr::filter(MARKERS %in% common.set)
  common.set <- NULL

  fst.select <- tibble::tibble(POP1 = pop.select[1], POP2 = pop.select[2]) %>%
    dplyr::bind_cols(
      assigner::compute_fst(
        x = data.select,
        ci = ci,
        iteration.ci = iteration.ci,
        quantiles.ci = quantiles.ci,
        digits = digits,
        path.folder = path.folder) %$%
        fst.overall
    )
  data.select <- pop.select <- NULL
  return(fst.select)
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

boot_ci <- carrier::crate(function(x, sigma.loc.alleles, digits = 9){
  `%>%` <- magrittr::`%>%`
  `%<>%` <- magrittr::`%<>%`
  `%$%` <- magrittr::`%$%`
  markers.list <- sigma.loc.alleles %>%
    dplyr::ungroup(.) %>%
    dplyr::distinct(MARKERS) %>%
    dplyr::arrange(MARKERS)

  subsample.markers <- markers.list %>%
    dplyr::sample_n(tbl = ., size = nrow(markers.list), replace = TRUE) %>%
    dplyr::arrange(MARKERS)

  fst.fis.overall.iterations <- sigma.loc.alleles %>%
    dplyr::right_join(subsample.markers, by = "MARKERS") %>%
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

# fst_subsample-----------------------------------------------------------------
#' @title fst_subsample
#' @description Function that link all with subsampling
#' @rdname fst_subsample
#' @export
#' @keywords internal

fst_subsample <- function(
  x,
  data,
  snprelate = FALSE,
  strata = NULL,
  holdout.samples = NULL,
  pairwise = FALSE,
  ci = FALSE,
  iteration.ci = 100,
  quantiles.ci = c(0.025,0.975),
  digits = 9,
  subsample = NULL,
  path.folder = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE
) {
  # x <- subsample.list[[1]] # test

  res <- list()# create list to store results
  # Managing subsampling -------------------------------------------------------
  subsample.id <- unique(x$SUBSAMPLE)
  if (!is.null(subsample) && (verbose)) message("Analyzing subsample: ", subsample.id)

  # genotyped data and holdout sample ------------------------------------------
  data %<>%
    dplyr::filter(INDIVIDUALS %in% x$ID_SEQ) %>%  # Keep only the subsample
    dplyr::filter(GT != "000000")
  x <- NULL #unused object


  # if holdout set, removes individuals
  if (!is.null(holdout.samples)) {
    message("Removing holdout individuals\nFst computation...")
    holdout.samples <- strata %>%
      dplyr::filter(INDIVIDUALS %in% holdout.samples) %$%
      ID_SEQ

    data %<>%
      dplyr::filter(!INDIVIDUALS %in% holdout.samples)
  }

  # Compute global Fst ---------------------------------------------------------
  if (verbose) message("Global fst...")
  global.res <- compute_fst(
    x = data,
    ci = ci,
    iteration.ci = iteration.ci,
    quantiles.ci = quantiles.ci,
    digits = digits,
    path.folder = path.folder
  )

  res <- append(res, global.res)
  global.res <- NULL # unsused object

  # Compute pairwise Fst -------------------------------------------------------
  if (pairwise) {
    if (verbose) message("Pairwise fst...")
    if (!is.factor(data$STRATA)) data$STRATA <- factor(data$STRATA)
    pop.list <- levels(data$STRATA) # pop list
    # all combination of populations
    pop.pairwise <- utils::combn(pop.list, 2, simplify = FALSE)
    npp <- length(pop.pairwise)
    if (verbose) message("Number of pairwise computations: ", npp)

    # Fst for all pairwise populations
    list.pair <- seq_len(npp)
    fst.all.pop <- assigner::assigner_future(
      .x = list.pair,
      .f = pairwise_fst,
      flat.future = "dfr",
      split.vec = FALSE,
      split.with = NULL,
      parallel.core = min(10L, parallel.core),
      pop.pairwise = pop.pairwise,
      data = data,
      ci = ci,
      iteration.ci = iteration.ci,
      quantiles.ci = quantiles.ci,
      path.folder = path.folder
    )

    # Table with Fst
    pairwise.fst <- fst.all.pop %>%
      dplyr::mutate(
        POP1 = factor(POP1, levels = pop.list),
        POP2 = factor(POP2, levels = pop.list)
        # N_MARKERS = as.integer(N_MARKERS)
      ) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), .fns = round, digits = digits))
    # }#End pairwise Fst

    # Matrix--------------------------------------------------------------------
    upper.mat.fst <- pairwise.fst %>%
      dplyr::select(POP1, POP2, FST) %>%
      tidyr::complete(data = ., POP1, POP2) %>%
      assigner::rad_wide(x = ., formula = "POP1 ~ POP2", values_from = "FST", values_fill = "") %>%
      dplyr::rename(POP = POP1)
    rn <- upper.mat.fst$POP # rownames
    upper.mat.fst <- as.matrix(upper.mat.fst[,-1])# make matrix without first column
    rownames(upper.mat.fst) <- rn

    # get the full matrix with identical lower and upper diagonal
    # the diagonal is filled with 0

    full.mat.fst <- upper.mat.fst # bk of upper.mat.fst
    lower.mat.fst <- t(full.mat.fst) # transpose
    # merge upper and lower matrix
    full.mat.fst[lower.tri(full.mat.fst)] <- lower.mat.fst[lower.tri(lower.mat.fst)]
    diag(full.mat.fst) <- "0"

    if (ci) {
      # bind upper and lower diagonal of matrix
      lower.mat.ci <- pairwise.fst %>%
        dplyr::select(POP1, POP2, CI_LOW, CI_HIGH) %>%
        tidyr::unite(data = ., CI, CI_LOW, CI_HIGH, sep = " - ") %>%
        tidyr::complete(data = ., POP1, POP2) %>%
        assigner::rad_wide(x = ., formula = "POP1 ~ POP2", values_from = "CI", values_fill = "") %>%
        dplyr::rename(POP = POP1)

      cn <- colnames(lower.mat.ci) # bk of colnames
      lower.mat.ci <- t(lower.mat.ci[,-1]) # transpose
      colnames(lower.mat.ci) <- cn[-1] # colnames - POP
      lower.mat.ci = as.matrix(lower.mat.ci) # matrix

      # merge upper and lower matrix
      pairwise.fst.ci.matrix <- upper.mat.fst # bk upper.mat.fst
      pairwise.fst.ci.matrix[lower.tri(pairwise.fst.ci.matrix)] <- lower.mat.ci[lower.tri(lower.mat.ci)]
    } else {
      pairwise.fst.ci.matrix <- "confidence intervals not selected"
    }

  } else {
    pairwise.fst <- "pairwise fst not selected"
    upper.mat.fst <- "pairwise fst not selected"
    full.mat.fst <- "pairwise fst not selected"
    pairwise.fst.ci.matrix <- "pairwise fst not selected"
  }

  res$pairwise.fst <- pairwise.fst
  res$pairwise.fst.upper.matrix <- upper.mat.fst
  res$pairwise.fst.full.matrix <- full.mat.fst
  res$pairwise.fst.ci.matrix <- pairwise.fst.ci.matrix
  return(res)
}#End fst_subsample

# heatmap_fst-------------------------------------------------------------------
#' @title heatmap_fst
#' @description Function that generate an Heatmap of Fst and CI values
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
  pairwise.fst.full.matrix,
  pairwise.fst.ci.matrix,
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

  if (missing(pairwise.fst.full.matrix) || missing(pairwise.fst.ci.matrix)) {
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
  # else {
    # if (length(pop.levels) != length(union(rownames(data.fst), colnames(data.fst)))) {
      # rlang::abort(message = "Contact author, problem with strata levels in heatmap")
    # }
  # }

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

  inv.levels <- rev(levels(data.fst$POP2))
  # pop.levels <- levels(data.fst$POP2)
  rounder <- function(x, digits) round(as.numeric(x), digits)
  if (max(stringi::stri_length(data.fst$FST), na.rm = TRUE) != digits) {
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


  data.fst %<>% dplyr::left_join(data.ci, by = c("POP1", "POP2"))
  # median.fst <- median(x = data.fst$FST, na.rm = TRUE)
  mean.fst <- mean(x = data.fst$FST, na.rm = TRUE)
  min.fst <- min(x = data.fst$FST, na.rm = TRUE)
  max.fst <- max(x = data.fst$FST, na.rm = TRUE)
  data.fst$POP2 <- factor(x = as.character(data.fst$POP2),
                          levels = inv.levels, ordered = TRUE)

  # data without CI for the lower diag
  if (n.s || round.num) {
    data.ci <- dplyr::filter(data.fst, stringi::stri_detect_fixed(str = CI, pattern = " - ")) %>%
      tidyr::separate(data = ., col = CI, into = c("LOW", "HIGH"), sep = " - ") %>%
      dplyr::mutate(dplyr::across(.cols = c("FST", "LOW", "HIGH"), .fns = rounder, digits = digits)) %>%
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
        # dplyr::mutate_at(.tbl = ., .vars = "FST", .funs = format, scientific = FALSE) %>%
        # dplyr::mutate(FST = dplyr::if_else(NS, paste0(FST, "*"), FST), NS = NULL)
      )
    }
    suppressWarnings(data.fst %<>% dplyr::bind_rows(data.ci))
    data.ci <- ns <- NULL
    data.fst$POP2 <- factor(x = as.character(data.fst$POP2), levels = inv.levels, ordered = TRUE)
    data.fst$POP1 <- factor(x = as.character(data.fst$POP1), levels = pop.levels, ordered = TRUE)
  }

  heatmap.fst <- ggplot2::ggplot(
    data = data.fst, ggplot2::aes(x = POP1, y = POP2, fill = FST)
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(x = POP1, y = POP2, label = CI),
                       color = "black", size = text.size, na.rm = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = color.low,
      mid = color.mid,
      high = color.high,
      # midpoint = median.fst,
      midpoint = mean.fst,
      na.value = "white",
      limit = NULL,
      # limit = c(min.fst, max.fst),
      space = "Lab"
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10, family = "Helvetica", face = "bold")
    )
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

# assigner_fst_stats -----------------------------------------------------------
#' @title assigner_fst_stats
#' @description Generate useful stats
#' @rdname assigner_fst_stats
#' @export
#' @keywords internal

assigner_fst_stats <- function(
  x,
  v,
  group.by = NULL,
  outliers = FALSE,
  overall = FALSE,
  keep.groups = "drop",
  digits = NULL
) {


  if (outliers) {
    message("not implemented yet")
    # OUTLIERS_LOW = Q25 - (1.5 * IQR),
    # OUTLIERS_HIGH = Q75 + (1.5 * IQR),
    # OUTLIERS_LOW = ifelse(OUTLIERS_LOW < MIN, MIN, OUTLIERS_LOW), # don'T use dplyr::if_else here... you don't want to preserve types
    # OUTLIERS_LOW_N = length(.data[[v]][.data[[v]] < OUTLIERS_LOW]),
    # OUTLIERS_HIGH = ifelse(OUTLIERS_HIGH > MAX, MAX, OUTLIERS_HIGH),
    # OUTLIERS_HIGH_N = length(.data[[v]][.data[[v]] > OUTLIERS_HIGH]),
  }

  # keep.groups <- match.arg(arg = keep.groups, choices = c("drop", "keep"))

  if (!is.null(group.by)) {
    x %<>% dplyr::group_by(.data[[group.by]])
  }

  if (overall) {
    x %<>%
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
        N_MARKERS_MEAN = mean(N_MARKERS),
        .groups = keep.groups
      )
  } else {
    x %<>%
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
        .groups = keep.groups
      )
  }

  if (!is.null(digits)) {
    x %<>% dplyr::mutate(dplyr::across(.cols = where(is.numeric), .fns = round, digits = digits))
  }

  return(x)
}#End assigner_fst_stats

# fst_stats-----------------------------------------------------------
#' @title fst_stats
#' @description Generate useful stats
#' @rdname fst_stats
#' @export
#' @keywords internal

fst_stats <- function(x, l, digits = 9L, m = NULL, s = NULL, subsample = FALSE) {
  res <- list()
  message(x)
  want <- c("sigma.loc", "fst.markers", "fst.ranked", "fst.overall",
            "fis.markers", "fis.overall", "pairwise.fst")

  #1. flatten the tibble
  if (x %in% want) {
    res1 <- purrr::map(l, x) %>% dplyr::bind_rows(.)
  } else {
    res1 <- purrr::map(l, x)
  }

  #2. generate the stats
  want <- c("fst.markers","fis.markers", "fst.overall", "fis.overall")
  if (x %in% want && subsample) {
    # message("Stats: ", x)

    # defaults:
    group.by <- NULL
    overall <- TRUE
    v <- "FIS"

    # mod
    if (x %in% c("fst.markers","fis.markers")) {
      group.by <- "MARKERS"
      overall <- FALSE
    }

    if (x %in% c("fst.markers","fst.overall")) v <- "FST"

    res1 %<>%
      assigner_fst_stats(
        x = .,
        v = v,
        group.by = group.by,
        overall = overall,
        digits = digits
      )

    # res$subsample  <- purrr::modify_in(
    #   .x = res$subsample,
    #   .where = list(i),
    #   .f = assigner_fst_stats,
    #   v = v,
    #   group.by = group.by,
    #   overall = overall,
    #   digits = digits
    # )
  }

  #3. adjust some objects
  # pairwise.fst
  if (x == "pairwise.fst" && subsample) {
    res1 %<>%
      dplyr::group_by(POP1, POP2) %>%
      dplyr::summarise(
        dplyr::across(
          tidyselect::everything(), .fns = mean, na.rm = TRUE
        ), .groups = "drop"
      ) %>%
      dplyr::mutate(
        ITERATIONS = rep(iteration.subsample, dplyr::n()),
        dplyr::across(.cols = c("POP1", "POP2"), .fns = as.integer)
      ) %>%
      dplyr::left_join(s, by = c("POP1" = "STRATA_SEQ")) %>%
      dplyr::select(-POP1, POP1 = POP_ID) %>%
      dplyr::left_join(s, by = c("POP2" = "STRATA_SEQ")) %>%
      dplyr::select(-POP2, POP2 = POP_ID) %>%
      dplyr::select(POP1, POP2, FST, N_MARKERS, tidyselect::everything())
  }

  #4. Change strata and markers
  if (x %in% c("sigma.loc", "fst.markers", "fst.ranked", "fis.markers")
  ) {
    res1 %<>% match_markers_meta(x = ., markers.meta = m)
  }

  want <- c("fst.plot", "pairwise.fst.upper.matrix",
            "pairwise.fst.full.matrix", "pairwise.fst.ci.matrix")
  if (x %in% want && !subsample) res1 <- res1[[1]]

  want <- c("pairwise.fst.upper.matrix", "pairwise.fst.full.matrix", "pairwise.fst.ci.matrix")
  if (x %in% want) {
    if (subsample) {
      res1 <- purrr::set_names(x = res1, nm = seq_len(length(res1)))
      res1 <- purrr::map(
        .x = res1,
        .f = change_matrix_strata,
        s = s
      )
    } else {
      res1 %<>% change_matrix_strata(x = ., s = s)
    }
  }

  if (x == "pairwise.fst" && !subsample) {
    res1 %<>%
      dplyr::mutate(dplyr::across(.cols = c("POP1", "POP2"), .fns = as.integer)) %>%
      dplyr::left_join(s, by = c("POP1" = "STRATA_SEQ")) %>%
      dplyr::select(-POP1, POP1 = POP_ID) %>%
      dplyr::left_join(s, by = c("POP2" = "STRATA_SEQ")) %>%
      dplyr::select(-POP2, POP2 = POP_ID) %>%
      dplyr::select(POP1, POP2, FST, N_MARKERS, tidyselect::everything())
  }

  res[[x]] <- res1
  return(res)
}#End fst_stats

#' @title change_matrix_strata
#' @description Integrate strata info back into the matrix
#' @rdname change_matrix_strata
#' @export
#' @keywords internal
change_matrix_strata <- function(x, s) {
  # x
  # class(colnames(x))
  # class(rownames(x))
  colnames(x) %<>%
    as.character(.) %>%
    stringi::stri_replace_all_regex(
      str = .,
      pattern = paste0("^", as.character(s$STRATA_SEQ)),
      replacement = as.character(s$POP_ID),
      vectorize_all = FALSE
    )
  # colnames(x)
  rownames(x) %<>%
    as.character(.) %>%
    stringi::stri_replace_all_regex(
      str = .,
      pattern = paste0("^", as.character(s$STRATA_SEQ)),
      replacement = as.character(s$POP_ID),
      vectorize_all = FALSE
    )
  #
  return(x)
}#change_matrix_strata

# match_markers_meta -----------------------------------------------------------
#' @title match_markers_meta
#' @description Integrate markers meta info back into the data
#' @rdname match_markers_meta
#' @export
#' @keywords internal
match_markers_meta <- function(x, markers.meta) {
    x  %<>%
      dplyr::rename(M_SEQ = MARKERS) %>%
      dplyr::left_join(markers.meta, by = "M_SEQ") %>%
      dplyr::select(-M_SEQ)
  return(x)
}# End match_strata

# fst_write -----------------------------------------------------------
#' @title fst_write
#' @description Write the fst results to working directory or path
#' @rdname fst_write
#' @export
#' @keywords internal
fst_write <- function(x,list.sub, path.folder) {
  readr::write_tsv(
    x = list.sub[[x]],
    file = file.path(path.folder, paste0(x, ".tsv"))
  )
}#End fst_write
