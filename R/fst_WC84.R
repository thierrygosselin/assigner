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
#' precedence over any grouping found input file (\code{data}). 
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
#' @importFrom radiator tidy_wide change_pop_names detect_biallelic_markers generate_folder radiator_dots
#' @importFrom tidyr separate gather spread unite
#' @importFrom purrr map flatten transpose flatten_int
#' @importFrom dplyr mutate mutate_if mutate_all summarise group_by ungroup select rename full_join left_join anti_join right_join semi_join filter n_distinct distinct arrange sample_n bind_rows bind_cols ntile desc n
#' @importFrom stats quantile
#' @importFrom utils count.fields combn
# @importFrom SNPRelate snpgdsOpen snpgdsClose snpgdsFst snpgdsCreateGeno
#' @importFrom stringi stri_replace_all_regex stri_join stri_replace_na stri_sub
#' @importFrom readr read_tsv
#' @importFrom parallel detectCores
#' @importFrom ggplot2 ggplot aes expand_limits geom_histogram labs theme element_blank element_text scale_colour_manual  facet_grid


#' @examples
#' \dontrun{
#' wombat.fst.pairwise <- fst_WC84(
#' data = "wombat.filtered.tidy.tsv", 
#' pop.levels = c("ATL", "MLE", "BIS", "PMO", "SOL", "TAS", "ECU"),
#' holdout.samples = NULL,
#' pairwise = TRUE,
#' ci = TRUE, 
#' iteration.ci = 10000, 
#' quantiles.ci = c(0.025,0.975),
#' parallel.core = 8,
#' verbose = TRUE,
#' filename = "wombat",
#' heatmap.fst = TRUE
#' )
#' To get the overall Fst estimate:
#' wombat.fst.pairwise$fst.overall
#' To get the Fst plot:
#' wombat.fst.pairwise$fst.plot
#' To get the pairwise Fst values with confidence intervals in a data frame:
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
#' it is adviseable to perform an overall test of population differentiation, 
#' possibly using a hierarchical population structure, (see AMOVA)'}
#' 
#' To compute an AMOVA, use \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive}
#' or \code{Phi_st_Meirmans} in \code{mmod}.
#' 
#' \code{hierfstat} is available on 
#' CRAN \url{http://cran.r-project.org/web/packages/hierfstat/} and 
#' github \url{https://github.com/jgx65/hierfstat/}
#' 
#' Link for \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive}
#' 
#' For Fisher's exact test and p-values per markers 
#' see \code{mmod} \code{diff_test}.
#' 
#' \code{\link[radiator]{tidy_genomic_data}} to transform numerous genomic data 
#' format in tidy data frames.

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
  
  # fst.snprelate <- NULL # remove after bias test
  # gds.file.connection <- NULL
  if (verbose) {
    cat("################################################################################\n")
    cat("############################# assigner::fst_WC84 ###############################\n")
    cat("################################################################################\n")
  }
  
  # Cleanup---------------------------------------------------------------------
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date/time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(width = 70)
  timing <- proc.time()# for timing
  res <- list()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(timing <- proc.time() - timing, add = TRUE)
  on.exit(if (verbose) message("\nComputation time, overall: ", round(timing[[3]]), " sec"), add = TRUE)
  on.exit(if (verbose) cat("############################# fst_WC84 completed ###############################\n"), add = TRUE)
  
  # Function call and dotslist -------------------------------------------------
  rad.dots <- radiator::radiator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE), 
    keepers = c("filter.monomorphic", "holdout.samples", "subsample",
                "iteration.subsample", "heatmap.fst", "blacklist.id"),
    verbose = verbose
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
  if (!is.null(filename)) filename <- stringi::stri_join(filename, "_fst_WC84")
  path.folder <- radiator::generate_folder(
    f = filename,
    file.date = file.date,
    verbose = verbose)
  
  
  if (snprelate) {
    # Check that snprelate is installed
    if (!"SNPRelate" %in% utils::installed.packages()[,"Package"]) {
      rlang::abort('Please install SNPRelate for this option:\n
                 install.packages("BiocManager")
                 BiocManager::install("SNPRelate")')
    }
    rlang::abort("Until the bias observed with SNPRelate is resolved, the option is unavailable.")
  }
  #   # snprelate <- FALSE
  #   stop("Until the bias observed with SNPRelate is resolved, the option will remain unavailable.")
  #   message("Fst computations with SNPRelate")
  #   if (ci) {
  #     message("Confidence Intervals are not implemented with SNPRelate, for now...")
  #     ci <- FALSE
  #   }
  # } else {
  #   message("Fst computations with assigner built-in function")
  # }
  
  # Import data ---------------------------------------------------------------
  if (verbose) message("Importing data")
  want <- c("MARKERS", "POP_ID", "STRATA", "INDIVIDUALS", "GT")
  input <- suppressWarnings(
    radiator::tidy_wide(data = data, import.metadata = FALSE) %>% 
  dplyr::select(dplyr::one_of(want))# remove unnecessary columns
  )
  
  # Strata----------------------------------------------------------------------
  strata.df <- radiator::read_strata(
    strata = strata,
    pop.id = TRUE,
    blacklist.id = blacklist.id,
    pop.levels = pop.levels,
    verbose = verbose) %$% strata
  
  # population levels and strata------------------------------------------------
  if (!is.null(strata)) {
    input <- radiator::join_strata(
      data = input, strata = strata.df, pop.id = TRUE, verbose = FALSE)
  } else {
    strata <- strata.df <- radiator::generate_strata(data = input)
  }
  
  # subsampling data------------------------------------------------------------
  # create the subsampling list
  if (is.null(subsample)) {
    iteration.subsample <- 1
  } else {
    if (subsample == "min") {
      subsample <- strata.df %>% 
        dplyr::group_by(POP_ID) %>% 
        dplyr::tally(.) %>% 
        dplyr::filter(n == min(n)) %>% 
        dplyr::ungroup(.) %>% 
        dplyr::select(n) %>% 
        purrr::flatten_int(.)
    }
  }
  
  subsample.list <- purrr::map(# map_df ?
    .x = 1:iteration.subsample,
    .f = subsampling_data,
    ind.pop.df = strata.df,
    subsample = subsample
  )
  
  # keep track of subsampling individuals and write to directory
  if (!is.null(subsample)) {
    if (verbose) message("Subsampling: selected")
    subsampling.individuals <- dplyr::bind_rows(subsample.list)
    readr::write_tsv(
      x = subsampling.individuals, 
      path = "assigner_fst_WC84_subsampling_individuals.tsv", 
      col_names = TRUE, 
      append = FALSE
    )
    res$subsampling.individuals <- subsampling.individuals
    
    # Note to myself: is this necessary ? Simplify...
    if (!is.null(filename)) {
      readr::write_tsv(
        x = res$subsampling.individuals, 
        path = file.path(path.folder, "subsampling.individuals.tsv"))
    }
  } # End subsampling
  
  # unused objects
  subsampling.individuals <- NULL
  
  # Calculations ----------------------------------------------------------------
  subsample.fst <- purrr::map(
    .x = subsample.list,
    .f = fst_subsample,
    input = input,
    snprelate = snprelate,
    pop.levels = pop.levels, 
    strata = strata,
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
  
  # Compile subsampling results ------------------------------------------------
  if (is.null(subsample)) {
    subsample.fst <- purrr::flatten(subsample.fst)
    
    # sigma.loc
    res$sigma.loc <- subsample.fst$sigma.loc
    if (!is.null(filename) && is.data.frame(res$sigma.loc)) { 
      readr::write_tsv(
        x = res$sigma.loc, 
        path = file.path(path.folder, "sigma.loc.tsv"))
    }
    
    # fst.markers
    res$fst.markers <- subsample.fst$fst.markers
    if (!is.null(filename) && is.data.frame(res$fst.markers)) { 
      readr::write_tsv(
        x = res$fst.markers, 
        path = file.path(path.folder, "fst.markers.tsv"))
    }
    # fst.ranked
    res$fst.ranked <- subsample.fst$fst.ranked
    if (!is.null(filename) && is.data.frame(res$fst.ranked)) { 
      readr::write_tsv(
        x = res$fst.ranked, 
        path = file.path(path.folder, "fst.ranked.tsv"))
    }
    
    # fst.overall
    res$fst.overall <- subsample.fst$fst.overall
    if (!is.null(filename) && is.data.frame(res$fst.overall)) { 
      readr::write_tsv(
        x = res$fst.overall, 
        path = file.path(path.folder, "fst.overall.tsv"))
    }
    # fis.markers
    res$fis.markers <- subsample.fst$fis.markers
    if (!is.null(filename) && is.data.frame(res$fis.markers)) {
      readr::write_tsv(
        x = res$fis.markers, 
        path = file.path(path.folder, "fis.markers.tsv"))
    }
    # fis.overall
    res$fis.overall <- subsample.fst$fis.overall
    if (!is.null(filename) && is.data.frame(res$fis.overall)) {  
      readr::write_tsv(
        x = res$fis.overall, 
        path = file.path(path.folder, "fis.overall.tsv"))
    }
    # fst.plot
    res$fst.plot <- subsample.fst$fst.plot
    if (!is.null(filename)) {  
      ggplot2::ggsave(
        filename = file.path(path.folder, "fst.plot.pdf"), 
        plot = res$fst.plot,
        width = 15, height = 10,
        dpi = 300, units = "cm", device = "pdf", limitsize = FALSE,
        useDingbats = FALSE)
    }
    
    # pairwise.fst
    res$pairwise.fst <- subsample.fst$pairwise.fst
    if (!is.null(filename) && is.data.frame(res$pairwise.fst)) {  
      readr::write_tsv(
        x = res$pairwise.fst, 
        path = file.path(path.folder, "pairwise.fst.tsv"))
    }
    # pairwise.fst.upper.matrix
    res$pairwise.fst.upper.matrix <- subsample.fst$pairwise.fst.upper.mat
    if (!is.null(filename)) {
      pairwise.fst.upper.matrix <- res$pairwise.fst.upper.matrix
      saveRDS(pairwise.fst.upper.matrix, file = file.path(path.folder, "pairwise.fst.upper.matrix.RData"))
      pairwise.fst.upper.matrix <- NULL
    }
    # pairwise.fst.full.matrix
    res$pairwise.fst.full.matrix <- subsample.fst$pairwise.fst.full.mat
    if (!is.null(filename)) {
      pairwise.fst.full.matrix <- res$pairwise.fst.full.matrix
      saveRDS(pairwise.fst.full.matrix, file = file.path(path.folder, "pairwise.fst.full.matrix.RData"))
      pairwise.fst.full.matrix <- NULL
    }
    # pairwise.fst.ci.matrix
    res$pairwise.fst.ci.matrix <- subsample.fst$pairwise.fst.ci.matrix
    if (!is.null(filename)) { 
      pairwise.fst.ci.matrix <- res$pairwise.fst.ci.matrix
      saveRDS(pairwise.fst.ci.matrix, file = file.path(path.folder, "pairwise.fst.ci.matrix.RData"))
      pairwise.fst.ci.matrix <- NULL
    }
  } else {
    # With SUBSAMPLING --------
    # test <- subsample.fst[1]
    subsample.fst.transposed <- purrr::transpose(subsample.fst)
    # names(subsample.fst.transposed)
    
    # sigma.loc
    res$sigma.loc.subsample <- subsample.fst.transposed[["sigma.loc"]]
    if (!is.null(filename) && is.data.frame(res$sigma.loc.subsample)) {
      readr::write_tsv(
        x = dplyr::bind_rows(res$sigma.loc.subsample), 
        path = file.path(path.folder, "sigma.loc.tsv"))
    }
    
    # fst.markers
    res$fst.markers.subsample <- dplyr::bind_rows(subsample.fst.transposed[["fst.markers"]]) %>% 
      dplyr::group_by(MARKERS) %>% 
      dplyr::summarise(
        MEAN = mean(FST),
        SE = sqrt(stats::var(FST)/length(FST)),
        MIN = min(FST),
        MAX = max(FST),
        MEDIAN = stats::median(FST),
        QUANTILE25 = stats::quantile(FST, 0.25),
        QUANTILE75 = stats::quantile(FST, 0.75),
        ITERATIONS = length(FST)
      ) %>% 
      dplyr::mutate_if(.tbl = ., .predicate =  is.numeric, .funs = round, digits = digits)
    
    if (!is.null(filename) && is.data.frame(res$fst.markers.subsample)) {
      readr::write_tsv(
        x = res$fst.markers.subsample, 
        path = file.path(path.folder, "fst.markers.tsv"))
    }
    # fst.ranked
    res$fst.ranked.subsample <- dplyr::bind_rows(subsample.fst.transposed[["fst.ranked"]])
    if (!is.null(filename) && is.data.frame(res$fst.ranked.subsample)) {
      readr::write_tsv(
        x = res$fst.ranked.subsample, 
        path = file.path(path.folder, "fst.ranked.tsv"))
    }
    # fst.overall
    res$fst.overall.subsample <- dplyr::bind_rows(subsample.fst.transposed[["fst.overall"]]) %>%
      dplyr::summarise(
        MEAN = mean(FST),
        SE = sqrt(stats::var(FST)/length(FST)),
        MIN = min(FST),
        MAX = max(FST),
        MEDIAN = stats::median(FST),
        QUANTILE25 = stats::quantile(FST, 0.25),
        QUANTILE75 = stats::quantile(FST, 0.75),
        ITERATIONS = length(FST),
        N_MARKERS_MEAN = mean(N_MARKERS)
      ) %>% 
      dplyr::mutate_all(.tbl = ., .funs = round, digits = digits)
    
    if (!is.null(filename) && is.data.frame(res$fst.overall.subsample)) {
      readr::write_tsv(
        x = res$fst.overall.subsample, 
        path = file.path(path.folder, "fst.overall.tsv"))
    }
    # fis.markers
    res$fis.markers.subsample <- dplyr::bind_rows(subsample.fst.transposed[["fis.markers"]]) %>% 
      dplyr::group_by(MARKERS) %>% 
      dplyr::summarise(
        MEAN = mean(FIS),
        SE = sqrt(stats::var(FIS)/length(FIS)),
        MIN = min(FIS),
        MAX = max(FIS),
        MEDIAN = stats::median(FIS),
        QUANTILE25 = stats::quantile(FIS, 0.25),
        QUANTILE75 = stats::quantile(FIS, 0.75),
        ITERATIONS = length(FIS)
      ) %>% 
      dplyr::mutate_if(.tbl = ., .predicate =  is.numeric, .funs = round, digits = digits)
    if (!is.null(filename) && is.data.frame(res$fis.markers.subsample)) {
      readr::write_tsv(
        x = res$fis.markers.subsample, 
        path = file.path(path.folder, "fis.markers.tsv"))
    }
    # fis.overall
    res$fis.overall.subsample <- dplyr::bind_rows(subsample.fst.transposed[["fis.overall"]]) %>%
      dplyr::summarise(
        MEAN = mean(FIS),
        SE = sqrt(stats::var(FIS)/length(FIS)),
        MIN = min(FIS),
        MAX = max(FIS),
        MEDIAN = stats::median(FIS),
        QUANTILE25 = stats::quantile(FIS, 0.25),
        QUANTILE75 = stats::quantile(FIS, 0.75),
        ITERATIONS = length(FIS),
        N_MARKERS_MEAN = mean(N_MARKERS)
      ) %>% 
      dplyr::mutate_all(.tbl = ., .funs = round, digits = digits)
    if (!is.null(filename) && is.data.frame(res$fis.overall.subsample)) {
      readr::write_tsv(
        x = res$fis.overall.subsample, 
        path = file.path(path.folder, "fis.overall.tsv"))
    }
    
    # fst.plot
    res$fst.plot.subsample <- subsample.fst.transposed[["fst.plot"]]
    
    # pairwise.fst
    res$pairwise.fst.subsample <- subsample.fst.transposed[["pairwise.fst"]]
    
    # pairwise.fst.mean
    res$pairwise.fst.subsample.mean <- dplyr::bind_rows(res$pairwise.fst.subsample) %>% 
      dplyr::group_by(POP1, POP2) %>% 
      dplyr::summarise_all(.tbl = ., .funs = mean, na.rm = TRUE) %>% 
      dplyr::mutate(ITERATIONS = rep(iteration.subsample, n()))
    
    if (!is.null(filename) && is.data.frame(res$pairwise.fst.subsample.mean)) {
      readr::write_tsv(
        x = res$pairwise.fst.subsample.mean, 
        path = file.path(path.folder, "pairwise.fst.tsv"))
    }
    # pairwise.fst.upper.matrix
    res$pairwise.fst.upper.matrix.subsample <- subsample.fst.transposed[["pairwise.fst.upper.matrix"]]
    
    # pairwise.fst.upper.matrix.mean
    res$pairwise.fst.upper.matrix.subsample.mean <- res$pairwise.fst.subsample.mean %>% 
      dplyr::select(POP1, POP2, FST) %>% 
      tidyr::spread(data = ., POP2, FST, fill = "", drop = FALSE) %>% 
      dplyr::rename(POP = POP1)
    rn <- res$pairwise.fst.upper.matrix.subsample.mean$POP # rownames
    res$pairwise.fst.upper.matrix.subsample.mean <- as.matrix(res$pairwise.fst.upper.matrix.subsample.mean[,-1])# make matrix without first column
    rownames(res$pairwise.fst.upper.matrix.subsample.mean) <- rn
    if (!is.null(filename)) {
      pairwise.fst.upper.matrix.subsample.mean <- res$pairwise.fst.upper.matrix.subsample.mean
      saveRDS(pairwise.fst.upper.matrix.subsample.mean, file = file.path(path.folder, "pairwise.fst.upper.matrix.RData"))
      pairwise.fst.upper.matrix.subsample.mean <- NULL
    }
    # pairwise.fst.full.matrix
    res$pairwise.fst.full.matrix.subsample <- subsample.fst.transposed[["pairwise.fst.full.matrix"]]
    
    # pairwise.fst.full.matrix.mean
    res$pairwise.fst.full.matrix.subsample.mean <- res$pairwise.fst.upper.matrix.subsample.mean # bk of upper.mat.fst
    lower.mat.fst <- t(res$pairwise.fst.full.matrix.subsample.mean) # transpose
    
    
    # merge upper and lower matrix
    res$pairwise.fst.full.matrix.subsample.mean[lower.tri(res$pairwise.fst.full.matrix.subsample.mean)] <- lower.mat.fst[lower.tri(lower.mat.fst)] 
    diag(res$pairwise.fst.full.matrix.subsample.mean) <- "0"
    
    if (!is.null(filename)) {
      pairwise.fst.full.matrix.subsample.mean <- res$pairwise.fst.full.matrix.subsample.mean
      saveRDS(pairwise.fst.full.matrix.subsample.mean, file = file.path(path.folder, "pairwise.fst.full.matrix.RData"))
      pairwise.fst.full.matrix.subsample.mean <- NULL
    }
    if (ci) {
      # pairwise.fst.ci.matrix
      res$pairwise.fst.ci.matrix.subsample <- subsample.fst.transposed[["pairwise.fst.ci.matrix"]]
      
      # pairwise.fst.ci.matrix.mean
      lower.mat.ci.sub <- res$pairwise.fst.subsample.mean %>% 
        dplyr::select(POP1, POP2, CI_LOW, CI_HIGH) %>% 
        tidyr::unite(data = ., CI, CI_LOW, CI_HIGH, sep = " - ") %>% 
        tidyr::spread(data = ., POP2, CI, fill = "", drop = FALSE) %>% 
        dplyr::rename(POP = POP1)
      
      cn <- colnames(lower.mat.ci.sub) # bk of colnames
      lower.mat.ci.sub <- t(lower.mat.ci.sub[,-1]) # transpose
      colnames(lower.mat.ci.sub) <- cn[-1] # colnames - POP
      lower.mat.ci.sub = as.matrix(lower.mat.ci.sub) # matrix
      
      # merge upper and lower matrix
      pairwise.fst.ci.matrix.sub <- res$pairwise.fst.upper.matrix.subsample.mean # bk upper.mat.fst
      pairwise.fst.ci.matrix.sub[lower.tri(pairwise.fst.ci.matrix.sub)] <- lower.mat.ci.sub[lower.tri(lower.mat.ci.sub)]
      res$pairwise.fst.ci.matrix.subsample.mean <- pairwise.fst.ci.matrix.sub
      if (!is.null(filename)) {
        pairwise.fst.ci.matrix.subsample.mean <- res$pairwise.fst.ci.matrix.subsample.mean
        saveRDS(pairwise.fst.ci.matrix.subsample.mean, file = file.path(path.folder, "pairwise.fst.ci.matrix.RData"))
        pairwise.fst.ci.matrix.subsample.mean <- NULL
      }
    } else {
      res$pairwise.fst.ci.matrix.subsample <- "confidence intervals not selected"
      res$pairwise.fst.ci.matrix.subsample.mean <- "confidence intervals not selected"
    }
  }
  
  # heatmap.fst ----------------------------------------------------------------
  if (heatmap.fst) {
    if (verbose) message("Generating heatmap fst...")
    res$heatmap.fst <- heatmap_fst(
      pairwise.fst.full.matrix = res$pairwise.fst.full.matrix, 
      pairwise.fst.ci.matrix = res$pairwise.fst.ci.matrix, 
      digits = digits, 
      pop.levels = pop.levels,
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
      message("Fst (overall): ", res$fst.overall.subsample$MEAN)
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

compute_fst <- function(
  x,
  ci = FALSE,
  iteration.ci = 100,
  quantiles.ci = c(0.025,0.975),
  digits = 9,
  path.folder = NULL,
  parallel.core = parallel::detectCores() - 2
) {
  # TEST
  # ci = FALSE
  # iteration.ci = 100
  # quantiles.ci = c(0.025,0.975)
  # digits = 9
  # path.folder = NULL
  ## x = data.genotyped
  # x <- dplyr::filter(data, GT != "000000")
  
  # Removing monomorphic markers------------------------------------------------
  x <- radiator::filter_monomorphic(data = x, internal = TRUE, verbose = FALSE, path.folder = path.folder)
  
  # number of marker used for computation
  n.markers <- length(unique(x$MARKERS))
  
  # count.locus <- dplyr::group_by(.data = x, MARKERS) %>%
  #   dplyr::summarise(
  #     NPL = length(unique(as.character(POP_ID))),# number of populations per locus
  #     NIL = n() # number of individuals per locus
  #   )
  count.locus <- dplyr::bind_cols(
    dplyr::distinct(.data = x, MARKERS, POP_ID) %>% dplyr::count(MARKERS, name = "NPL"),# number of populations per locus
    dplyr::distinct(.data = x, MARKERS, INDIVIDUALS) %>% dplyr::count(MARKERS, name = "NIL")# number of individuals per locus
  )
  if (!identical(count.locus$MARKERS, count.locus$MARKERS1)) {
    rlang::abort("contact author")
  } else {
    count.locus %<>% dplyr::select(-MARKERS1)
  }
  
  # count.locus.pop <- dplyr::count(x, POP_ID, MARKERS, name = "NIPL") %>%
  #   dplyr::mutate(NIPL_SQ = NIPL ^ 2) %>%
  #   dplyr::group_by(MARKERS) %>%
  #   dplyr::summarise(NIPL_SQ_SUM = sum(NIPL_SQ, na.rm = TRUE)) %>%
  #   dplyr::full_join(count.locus, by = "MARKERS") %>%
  #   dplyr::mutate(NC = (NIL - NIPL_SQ_SUM / NIL) / (NPL - 1))#correction
  
  
  # faster:
  count.locus.pop <- dplyr::bind_cols(
    count.locus,
    dplyr::count(x, POP_ID, MARKERS, name = "NIPL") %>%
      dplyr::mutate(NIPL_SQ = NIPL ^ 2) %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(NIPL_SQ_SUM = sum(NIPL_SQ, na.rm = TRUE))
  )
  if (!identical(count.locus.pop$MARKERS, count.locus.pop$MARKERS1)) {
    rlang::abort("contact author")
  } else {
    count.locus.pop %<>%
      dplyr::select(-MARKERS1) %>%
      dplyr::mutate(NC = (NIL - NIPL_SQ_SUM / NIL) / (NPL - 1))#correction
  }
  count.locus <- NULL
  
  # numbers corrected
  if (n.markers > 30000) {
    allele.locus <- radiator::separate_gt(
      x = dplyr::select(x, MARKERS, POP_ID, INDIVIDUALS, GT),
      sep = 3,
      gt = "GT",
      gather = TRUE,
      haplotypes = FALSE,
      exclude = c("MARKERS", "POP_ID", "INDIVIDUALS"),
      parallel.core = parallel.core
    )
  } else {
    allele.locus <- radiator::separate_gt(
      x = dplyr::select(x, MARKERS, POP_ID, INDIVIDUALS, GT),
      sep = 3,
      gt = "GT",
      gather = TRUE,
      haplotypes = FALSE,
      exclude = c("MARKERS", "POP_ID", "INDIVIDUALS"),
      parallel.core = 1
    )
  }
  
  correction <- dplyr::distinct(.data = allele.locus, MARKERS, ALLELES) %>%
    dplyr::full_join(count.locus.pop, by = "MARKERS") %>%
    dplyr::arrange(MARKERS, ALLELES)
  count.locus.pop <- NULL
  
  freq.alleles <- dplyr::count(x = allele.locus, MARKERS, ALLELES, POP_ID) %>%
    tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, ALLELES), fill = list(n = 0)) %>%
    dplyr::group_by(MARKERS, POP_ID) %>%
    dplyr::mutate(
      NAPL = sum(n),
      FREQ_APL = n / NAPL # Frequency of alleles per pop and locus
    ) %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::mutate(FREQ_AL = sum(n) / sum(NAPL)) %>% #Frequency of alleles per locus
    dplyr::full_join(correction, by = c("MARKERS", "ALLELES")) %>%
    dplyr::arrange(MARKERS, POP_ID)
  correction <- allele.locus <- NULL
  
  fst.stats.prep <- dplyr::mutate(
    .data = x,
    het = dplyr::if_else(stringi::stri_sub(GT, 1, 3) != stringi::stri_sub(GT, 4, 6), 1, 0)
  )
  
  if (n.markers > 30000) {
    fst.stats.prep %<>% radiator::separate_gt(
      x =.,
      sep = 3,
      gt = "GT",
      gather = TRUE,
      haplotypes = FALSE,
      exclude = c("MARKERS", "POP_ID", "INDIVIDUALS", "het"),
      parallel.core = parallel.core
    )
  } else {
    fst.stats.prep %<>% radiator::separate_gt(
      x =.,
      sep = 3,
      gt = "GT",
      gather = TRUE,
      haplotypes = FALSE,
      exclude = c("MARKERS", "POP_ID", "INDIVIDUALS", "het"),
      parallel.core = 1
    )
  }
  
  fst.stats.prep %<>%
    dplyr::select(-ALLELE_GROUP) %>%
    dplyr::group_by(MARKERS, POP_ID, ALLELES) %>%
    dplyr::summarise(MHO = length(het[het == 1])) %>%
    dplyr::ungroup(.) %>%
    tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, ALLELES), fill = list(MHO = 0)) %>%
    dplyr::arrange(MARKERS, ALLELES, POP_ID) %>%
    dplyr::full_join(freq.alleles, by = c("POP_ID", "MARKERS", "ALLELES")) %>%
    dplyr::mutate(
      NIPL = NAPL/2,
      MHOM = round(((NAPL * FREQ_APL - MHO)/2), 0),
      dum = NIPL * (FREQ_APL - 2 * FREQ_APL^2) + MHOM
    ) %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::mutate(
      SSi = sum(dum, na.rm = TRUE),
      dum1 = NIPL * (FREQ_APL - FREQ_AL)^2,
      SSP = 2 * sum(dum1, na.rm = TRUE)
    ) %>%
    dplyr::group_by(MARKERS, POP_ID) %>%
    dplyr::mutate(SSG = NIPL * FREQ_APL - MHOM) %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::mutate(
      sigw = round(sum(SSG, na.rm = TRUE), 2)/NIL,# ntal -> NIL
      MSP = SSP/(NPL - 1),
      MSI = SSi/(NIL - NPL),
      sigb = 0.5 * (MSI - sigw),
      siga = 1/2/NC * (MSP - MSI)
    ) %>%
    dplyr::ungroup(.)
  freq.alleles <- NULL
  
  # variance components of allele frequencies for each allele
  # siga: among populations
  # sigb: among individuals within/between populations
  # sigw: within individuals
  sigma.loc.alleles <- fst.stats.prep %>%
    dplyr::group_by(MARKERS, ALLELES) %>%
    dplyr::summarise(
      siga = mean(siga, na.rm = TRUE),
      sigb = mean(sigb, na.rm = TRUE),
      sigw = mean(sigw, na.rm = TRUE)
    ) %>%
    dplyr::ungroup(.)
  fst.stats.prep <- NULL
  # variance components per locus
  # lsiga: among populations
  # lsigb: among individuals within/between populations
  # lsigw: within individuals
  
  sigma.loc <- sigma.loc.alleles %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::summarise(
      lsiga = round(sum(siga, na.rm = TRUE), digits),
      lsigb = round(sum(sigb, na.rm = TRUE), digits),
      lsigw = round(sum(sigw, na.rm = TRUE), digits)
    ) %>%
    dplyr::ungroup(.)
  
  fst.fis.markers <- sigma.loc %>%
    dplyr::group_by(MARKERS) %>%
    dplyr::summarise(
      FST = round(lsiga/(lsiga + lsigb + lsigw), digits),
      FIS = round(lsigb/(lsigb + lsigw), digits)
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
  if (ci) {
    # the function:
    boot.fst.list <- purrr::map(
      .x = 1:iteration.ci,
      .f = boot_ci,
      sigma.loc.alleles = sigma.loc.alleles,
      digits = digits
    )
    boot.fst <- dplyr::bind_rows(boot.fst.list)
    boot.fst.summary <- boot.fst %>%
      dplyr::summarise(
        CI_LOW = round(stats::quantile(FST,
                                       probs = quantiles.ci[1],
                                       na.rm = TRUE),
                       digits),
        CI_HIGH = round(stats::quantile(FST,
                                        probs = quantiles.ci[2],
                                        na.rm = TRUE),
                        digits)
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
      RANKING = seq(from = 1, to = n()),
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
} # End compute_fst function

# pairwise_fst------------------------------------------------------------------
#' @title pairwise_fst
#' @description Pairwise Fst function
#' @rdname pairwise_fst
#' @export
#' @keywords internal

pairwise_fst <- function(
  list.pair, 
  pop.pairwise = NULL,
  # unique.markers.pop = NULL,
  data.genotyped = NULL,
  ci = FALSE, 
  iteration.ci = 100, 
  quantiles.ci = c(0.025,0.975), 
  digits = 9,
  path.folder = path.folder
) {
  
  pop.select <- stringi::stri_join(purrr::flatten(pop.pairwise[list.pair]))
  data.select <- data.genotyped %>%
    dplyr::filter(POP_ID %in% pop.select) %>%
    dplyr::mutate(POP_ID = droplevels(x = POP_ID))
  
  # common markers
  common.set <- intersect(
    unique(data.select$MARKERS[data.select$POP_ID == pop.select[1]]),
    unique(data.select$MARKERS[data.select$POP_ID == pop.select[2]])
  )
  
  data.select %<>% dplyr::filter(MARKERS %in% common.set)
  common.set <- NULL
  
  fst.select <- tibble::tibble(POP1 = pop.select[1], POP2 = pop.select[2]) %>% 
    dplyr::bind_cols(
      compute_fst(
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
} # End pairwise_fst

# pairwise_fst_snprelate--------------------------------------------------------

# @title pairwise_fst_snprelate
# @description Pairwise Fst function with SNPRelate
# @rdname pairwise_fst_snprelate
# @export
# @keywords internal

# pairwise_fst_snprelate <- function(pop.pairwise, data, strata, unique.markers.pop) {
#   
#   strata.df <- dplyr::filter(.data = strata, POP_ID %in% pop.pairwise) %>% # filter the pop
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
#     population = strata.df$POP_ID, # factors required
#     sample.id = strata.df$INDIVIDUALS,
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

boot_ci <- function(x, sigma.loc.alleles, digits = 9){
  
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
      ITERATIONS = rep(x, n()),
      FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0)
    )
  return(fst.fis.overall.iterations)
} # End boot_ci function

# fst_subsample-----------------------------------------------------------------
#' @title fst_subsample
#' @description Function that link all with subsampling
#' @rdname fst_subsample
#' @export
#' @keywords internal

fst_subsample <- function(
  x,
  input,
  snprelate = FALSE,
  pop.levels = NULL, 
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
  res <- list()# create list to store results
  
  # Managing subsampling -------------------------------------------------------
  # x <- subsample.list[[1]] # test
  subsample.id <- unique(x$SUBSAMPLE)
  
  if (!is.null(subsample)) {
    if (verbose) message("Analyzing subsample: ", subsample.id)
  }
  
  # Keep only the subsample
  input <- dplyr::semi_join(input, x, by = c("POP_ID", "INDIVIDUALS"))
  x <- NULL #unused object
  
  
  # SNPRelate and other prep ---------------------------------------------------
  # if (snprelate) {
  #   message("Detect if data is biallelic")
  #   biallelic <- radiator::detect_biallelic_markers(data = input)
  #   
  #   if (!biallelic) {
  #     warning("Data is not biallelic = cannot run assigner with SNPRelate Fst function")
  #     message("Continuing the Fst computations with built in function...")
  #     snprelate <- FALSE
  #   }
  #   
  #   # if holdout set, removes individuals
  #   if (!is.null(holdout.samples)) {
  #     message("Removing holdout individuals\nFst computation...")
  #     input <- dplyr::filter(.data = input, !INDIVIDUALS %in% holdout.samples)
  #   }
  #   
  #   unique.markers.pop <- input %>% 
  #     dplyr::filter(!is.na(GT_BIN)) %>%
  #     dplyr::distinct(MARKERS, POP_ID)
  #   
  # } else {
  # genotyped data and holdout sample
  data.genotyped <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
    dplyr::filter(GT != "000000")
  
  # if holdout set, removes individuals
  if (!is.null(holdout.samples)) {
    message("Removing holdout individuals\nFst computation...")
    data.genotyped <- dplyr::filter(.data = data.genotyped, !INDIVIDUALS %in% holdout.samples)
  }
  #No longer necessary: duplicate the dataset during parallel process...
  # unique.markers.pop <- dplyr::distinct(.data = data.genotyped, MARKERS, POP_ID)
  
  # Compute global Fst ---------------------------------------------------------
  # if (snprelate) {
  #   if (verbose) message("Generating GDS format...")
  #   
  #   strata.df <- dplyr::distinct(gds.genotypes, POP_ID, INDIVIDUALS) %>%
  #     dplyr::mutate(POP_ID = factor(POP_ID))
  #   
  #   snp.id <- dplyr::distinct(.data = gds.genotypes, MARKERS) %>%
  #     dplyr::arrange(MARKERS) %>%
  #     purrr::flatten_chr(.)
  #   
  #   gds.genotypes <- suppressWarnings(
  #     dplyr::select(.data = gds.genotypes, MARKERS, INDIVIDUALS, GT_BIN) %>%
  #       dplyr::group_by(MARKERS) %>% 
  #       tidyr::spread(data = ., key = INDIVIDUALS, value = GT_BIN) %>% 
  #       dplyr::arrange(MARKERS) %>%
  #       tibble::column_to_rownames(df = ., var = "MARKERS") %>% 
  #       data.matrix(.)
  #   )  
  #   
  #   # gds.data <- NULL
  #   SNPRelate::snpgdsCreateGeno(
  #     gds.fn = "assigner.gds",
  #     genmat = gds.genotypes,
  #     sample.id = strata.df$INDIVIDUALS,
  #     snp.id = snp.id,
  #     snp.rs.id = NULL,
  #     snp.chromosome = NULL,
  #     snp.position = NULL,
  #     snp.allele = NULL,
  #     snpfirstdim = TRUE,
  #     compress.annotation = "ZIP_RA.max",
  #     compress.geno = "",
  #     other.vars = NULL
  #   )
  #   
  #   # Compute the global Fst
  #   if (verbose) message("Global fst calculation")
  #   gds.file.connection <- SNPRelate::snpgdsOpen("assigner.gds")
  #   fst.snprelate <- SNPRelate::snpgdsFst(
  #     gdsobj = gds.file.connection,
  #     population = strata.df$POP_ID,
  #     sample.id = strata.df$INDIVIDUALS,
  #     snp.id = NULL,
  #     method = "W&C84",
  #     remove.monosnp = TRUE,
  #     maf = NaN,
  #     missing.rate = NaN,
  #     autosome.only = FALSE,
  #     with.id = FALSE,
  #     verbose = FALSE
  #   )$Fst
  #   
  #   fst.overall <- tibble::tibble(FST = round(fst.snprelate, digits)) %>% 
  #     dplyr::mutate(FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0))
  #   
  #   res$fst.overall <- fst.overall
  #   
  #   # close SNPRelate connection if no more computation with SNPRelate
  #   if (!pairwise) {
  #     SNPRelate::snpgdsClose(gds.file.connection)
  #   }
  
  
  # } else {
  if (verbose) message("Global fst calculation")
  global.res <- compute_fst(
    x = data.genotyped, 
    ci = ci, 
    iteration.ci = iteration.ci, 
    quantiles.ci = quantiles.ci, 
    digits = digits,
    path.folder = path.folder
  )
  res$sigma.loc <- global.res$sigma.loc
  res$fst.markers <- global.res$fst.markers
  res$fst.ranked <- global.res$fst.ranked
  res$fst.overall <- global.res$fst.overall
  res$fis.markers <- global.res$fis.markers
  res$fis.overall <- global.res$fis.overall
  res$fst.plot <- global.res$fst.plot
  global.res <- NULL # unsused object
  # }
  
  # Compute pairwise Fst -------------------------------------------------------
  if (pairwise) {
    if (verbose) message("Pairwise fst calculation")
    
    if (!is.factor(input$POP_ID)) input$POP_ID <- factor(input$POP_ID)
    pop.list <- levels(input$POP_ID) # pop list
    # all combination of populations
    pop.pairwise <- utils::combn(pop.list, 2, simplify = FALSE) 
    
    # if (snprelate) {
    #   fst.all.pop <- purrr::map(
    #     .x = pop.pairwise,
    #     .f = pairwise_fst_snprelate,
    #     data = gds.file.connection,
    #     strata = strata.df,
    #     unique.markers.pop = unique.markers.pop
    #   )
    #   # Table with Fst
    #   pairwise.fst <- t(data.frame(pop.pairwise)) %>% 
    #     tibble::as_tibble(.) %>%
    #     dplyr::bind_cols(dplyr::bind_rows(fst.all.pop)) %>% 
    #     dplyr::rename(POP1 = V1, POP2 = V2, FST = Fst) %>% 
    #     dplyr::mutate(
    #       POP1 = factor(POP1, levels = pop.list, ordered = TRUE),
    #       POP2 = factor(POP2, levels = pop.list, ordered = TRUE),
    #       FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0),
    #       FST = round(FST, digits)
    #     )
    #   
    #   # close SNPRelate connection
    #   SNPRelate::snpgdsClose(gds.file.connection)
    
    # } else {
    # Fst for all pairwise populations
    list.pair <- 1:length(pop.pairwise)
    # list.pair <- 5 #  test
    fst.all.pop <- .assigner_parallel(
      # fst.all.pop <- mclapply(
      X = list.pair, 
      FUN = pairwise_fst, 
      mc.preschedule = FALSE, 
      mc.silent = FALSE, 
      mc.cores = parallel.core,
      pop.pairwise = pop.pairwise,
      # unique.markers.pop = unique.markers.pop,
      data.genotyped = data.genotyped,
      ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci,
      path.folder = path.folder
    )
    # Table with Fst
    pairwise.fst <- dplyr::bind_rows(fst.all.pop) %>% 
      dplyr::mutate(
        POP1 = factor(POP1, levels = pop.list, ordered = TRUE),
        POP2 = factor(POP2, levels = pop.list, ordered = TRUE)
      )
    # }#End pairwise Fst
    
    # Matrix--------------------------------------------------------------------
    upper.mat.fst <- pairwise.fst %>% 
      dplyr::select(POP1, POP2, FST) %>% 
      tidyr::spread(data = ., POP2, FST, fill = "", drop = FALSE) %>% 
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
        tidyr::spread(data = ., POP2, CI, fill = "", drop = FALSE) %>% 
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
#' @param n.s (optional, logical) To have an * when the Fst value is not 
#' significative (0 is the lower bound of the CI).
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
#' @param plot.size (optional, integer) By default the size of the plot is set
#' to 40 cm x 40 cm.
#' Default: \code{plot.size = 40}.
#' @param pop.levels (optional, character) the pop.levels to have the pop ordered as desired.
#' Default: \code{pop.levels = NULL}.
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
  n.s = TRUE, 
  digits = 5, 
  color.low = "blue", 
  color.mid = "yellow", 
  color.high = "red",
  text.size = 4,
  plot.size = 40,
  pop.levels = NULL,
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
  # pop.levels = NULL
  # filename = NULL
  # path.folder = NULL
  
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
  
  
  data.fst %<>%
    radiator::distance2tibble(x = ., remove.diag = FALSE, na.diag = TRUE, remove.lower = FALSE, relative = FALSE, pop.levels = pop.levels) %>% 
    magrittr::set_colnames(x = ., value = c("POP1", "POP2", "FST"))
  inv.levels <- rev(levels(data.fst$POP2))
  pop.levels <- levels(data.fst$POP2)
  if (max(stringi::stri_length(data.fst$FST), na.rm = TRUE) != digits) {
    round.num <- TRUE
    rounder <- function(x, digits) round(as.numeric(x), digits)
  } else {
    round.num <- FALSE
  }
  
  if (n.s) round.num <- TRUE
  
  data.ci %<>%
    radiator::distance2tibble(
      x = ., remove.diag = FALSE, na.diag = TRUE,
      remove.lower = FALSE, relative = FALSE,
      distance.class.double = FALSE, pop.levels = pop.levels) %>% 
    magrittr::set_colnames(x = ., value = c("POP1", "POP2", "CI"))
  
  data.fst %<>% dplyr::left_join(data.ci, by = c("POP1", "POP2"))
  # median.fst <- median(x = data.fst$FST, na.rm = TRUE)
  mean.fst <- mean(x = data.fst$FST, na.rm = TRUE)
  min.fst <- min(x = data.fst$FST, na.rm = TRUE)
  max.fst <- max(x = data.fst$FST, na.rm = TRUE)
  data.fst$POP2 <- factor(x = as.character(data.fst$POP2), 
                          levels = inv.levels, ordered = TRUE)
  
  # data without CI
  if (n.s || round.num) {
    data.ci <- dplyr::filter(data.fst, stringi::stri_detect_fixed(str = CI, pattern = " - ")) %>% 
      tidyr::separate(data = ., col = CI, into = c("LOW", "HIGH"), sep = " - ") %>% 
      dplyr::mutate_at(.tbl = ., .vars = c("LOW", "HIGH"), .funs = rounder, digits = 5) %>% 
      dplyr::mutate(NS = dplyr::if_else(LOW == 0, TRUE, FALSE)) %>% 
      dplyr::mutate_at(.tbl = ., .vars = c("LOW", "HIGH"), .funs = format, scientific = FALSE) %>% 
      tidyr::unite(data = ., col = CI, c("LOW", "HIGH"), sep = "\n")
    ns <- dplyr::distinct(data.ci, POP1, POP2, NS) %>% 
      dplyr::rename(POP3 = POP2, POP2 = POP1) %>% 
      dplyr::rename(POP1 = POP3)
    data.ci %<>% dplyr::select(-NS)
    # data.fst 
    suppressWarnings(
      data.fst %<>% 
        dplyr::filter(!stringi::stri_detect_fixed(str = CI, pattern = " - ") | is.na(CI)) %>%
        dplyr::mutate_at(.tbl = ., .vars = c("CI"), .funs = rounder, digits = 5) %>%
        dplyr::mutate_at(.tbl = ., .vars = c("CI"), .funs = format, scientific = FALSE) %>% 
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
    data.fst$POP2 <- factor(x = as.character(data.fst$POP2), levels = inv.levels, ordered = TRUE)
    data.fst$POP1 <- factor(x = as.character(data.fst$POP1), levels = pop.levels, ordered = TRUE)
  }
  
  heatmap.fst <- ggplot2::ggplot(
    data = data.fst, 
    ggplot2::aes(x = POP1, y = POP2, fill = FST)) + 
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
      limit = c(min.fst, max.fst), space = "Lab") +
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
    heatmap.n <- stringi::stri_join(filename, "_heatmap.fst.pdf")
    ggplot2::ggsave(
      filename = file.path(path.folder, heatmap.n),
      plot = heatmap.fst,
      width = plot.size, height = plot.size,
      dpi = 300, units = "cm", device = "pdf", limitsize = FALSE,
      useDingbats = FALSE)
  }
  return(heatmap.fst)
}#End heatmap_fst

