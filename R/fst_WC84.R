# Compute Weir and Cockerham (1984) Fst

#' @name fst_WC84

#' @title A fast implementation of Weir and Cockerham (1984) Fst/Theta 
#' (overall and paiwise estimates)

#' @description The function computes Weir and Cockerham (1984) 
#' Fst for diploid genomes. Both overall and pairwise Fst can be estimated with 
#' confidence intervals based on bootstrap of markers (resampling with replacement). 
#' The function gives identical results \emph{at the 9th decimal} when tested 
#' against \code{\link[hierfstat]{genet.dist}} in \pkg{hierfstat}. Using the 
#' argument \code{snprelate = TRUE} will compute the Fst with 
#' \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}. This implementation
#' gives slightly upward bias values but provided the fastest computations I know,
#' but it doesn't compute confidence intervals, for now. 
#' For an R implementation, \code{\link{fst_WC84}} is very fast. 
#' The computations takes advantage of \pkg{dplyr}, \pkg{tidyr}, \pkg{purrr}, 
#' \pkg{data.table}, \pkg{parallel} and \pkg{SNPRelate}.
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


#' @param data A file in the working directory or object in the global environment 
#' in wide or long (tidy) formats. To import, the function uses internally
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}}. See details for more info.
#' 
#' \emph{How to get a tidy data frame ?}
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{tidy_genomic_data}} can transform 6 genomic data formats 
#' in a tidy data frame. You can also use this function to filter your dataset using
#' whitelist of markers, blacklist of individuals and genotypes.

#' @param snprelate (logical) Use \href{https://github.com/zhengxwen/SNPRelate}{SNPRelate}
#' to compute the Fst. Testing as shown an upward bias with \code{SNPRelate::snpgdsFst} function.
#' The values are 0.99 correlated with my built in codes and with \code{Hierfstat}.
#' But it's the fastest computation I've seen so far!
#' Default: \code{snprelate = FALSE}

#' @param pop.levels (optional, string) This refers to the levels in a factor. In this 
#' case, the id of the pop.
#' Use this argument to have the pop ordered your way instead of the default 
#' alphabetical or numerical order. e.g. \code{pop.levels = c("QUE", "ONT", "ALB")} 
#' instead of the default \code{pop.levels = c("ALB", "ONT", "QUE")}. 
#' Default: \code{pop.levels = NULL}.

#' @param pop.labels (optional, string) Use this argument to rename/relabel
#' your pop or combine your pop. e.g. To combine \code{"QUE"} and \code{"ONT"} 
#' into a new pop called \code{"NEW"}:
#' (1) First, define the levels for your pop with \code{pop.levels} argument: 
#' \code{pop.levels = c("QUE", "ONT", "ALB")}. 
#' (2) then, use \code{pop.labels} argument: 
#' \code{pop.levels = c("NEW", "NEW", "ALB")}.#' 
#' To rename \code{"QUE"} to \code{"TAS"}:
#' \code{pop.labels = c("TAS", "ONT", "ALB")}.
#' Default: \code{pop.labels = NULL}. If you find this too complicated, there is also the
#' \code{strata} argument that can do the same thing, see below.

#' @param strata (optional, data frame) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' If a \code{strata} file is specified, the strata file will have
#' precedence over any grouping found input file (\code{data}). 
#' The \code{STRATA} column can be any hierarchical grouping.
#' Default: \code{strata = NULL}.

#' @param holdout.samples (optional, data frame) Samples that don't participate in the Fst 
#' computation (supplementary). Data frame with one column \code{INDIVIDUALS}.
#' Default: \code{holdout.samples = NULL}.

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

#' @param ... other parameters passed to the function.

#' @return With pairwise comparison computed, the function returns a list with 
#' 11 objects:
#' \itemize{
#'   \item \code{$sigma.loc}: the variance components per locus 
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
#' }

#' @details \strong{Input data:}
#'  
#' To discriminate the long from the wide format, 
#' the function \pkg{stackr} \code{\link[stackr]{read_long_tidy_wide}} searches 
#' for \code{MARKERS or LOCUS} in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns: 
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals), 
#' the remaining columns are the markers in separate columns storing genotypes.
#' 
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info. 
#' (e.g. from a VCF see \pkg{stackr} \code{\link{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID}, 
#' \code{MARKERS or LOCUS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#' 
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 of 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 of 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#' 
#' \emph{How to get a tidy data frame ?}
#' \pkg{stackr} \code{\link{tidy_genomic_data}} can transform 6 genomic data formats 
#' in a tidy data frame.

#' @export
#' @rdname fst_WC84
#' @import parallel
#' @import ggplot2
#' @importFrom stackr read_long_tidy_wide discard_monomorphic_markers keep_common_markers change_pop_names detect_biallelic_markers
#' @importFrom tidyr separate gather spread unite
#' @importFrom purrr map flatten
#' @importFrom dplyr mutate summarise group_by ungroup select rename full_join left_join anti_join right_join semi_join filter n_distinct distinct arrange sample_n bind_rows bind_cols
#' @importFrom stats quantile
#' @importFrom utils count.fields combn
#' @importFrom SNPRelate snpgdsOpen snpgdsClose snpgdsFst snpgdsCreateGeno
#' @importFrom tibble data_frame column_to_rownames has_name as_data_frame
#' @importFrom stringi stri_replace_all_regex stri_join stri_replace_na stri_sub
#' @importFrom readr read_tsv

#' @examples
#' \dontrun{
#' wombat.fst.pairwise <- fst_WC84(
#' data = "wombat.filtered.tidy.tsv", 
#' sep = "/", 
#' pop.levels = c("ATL", "MLE", "BIS", "PMO", "SOL", "TAS", "ECU"),
#' holdout.samples = NULL,
#' pairwise = TRUE,
#' ci = TRUE, 
#' iteration.ci = 10000, 
#' quantiles.ci = c(0.025,0.975),
#' parallel.core = 8,
#' verbose = TRUE
#' )
#' To get the overall Fst estimate:
#' wombat.fst.pairwise$fst.overall
#' To get the Fst plot:
#' wombat.fst.pairwise$fst.plots
#' To get the pairwise Fst values with confidence intervals in a data frame:
#' wombat.fst.pairwise$pairwise.fst
#' To get the pairwise fst and ci matrix in a data frame:
#' # rename, data frame, put rownames in column
#' pairwise.fst.ci.df <- data.frame(pairwise.fst.ci.matrix) %>% add_rownames("POP")
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
#' or \code{\link[mmod]{Phi_st_Meirmans}} 
#' in \code{\link[mmod]{mmod}}.
#' 
#' \code{hierfstat} is available on 
#' CRAN \url{http://cran.r-project.org/web/packages/hierfstat/} and 
#' github \url{https://github.com/jgx65/hierfstat/}
#' 
#' Link for \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive}
#' 
#' For Fisher's exact test and p-values per markers 
#' see \code{mmod} \code{\link[mmod]{diff_test}}.
#' 
#' \code{\link[stackr]{tidy_genomic_data}} to transform numerous genomic data 
#' format in tidy data frames.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# Fst function: Weir & Cockerham 1984
fst_WC84 <- function(
  data,
  snprelate = FALSE,
  pop.levels = NULL, 
  pop.labels = NULL, 
  strata = NULL,
  holdout.samples = NULL,
  pairwise = FALSE,
  ci = FALSE,
  iteration.ci = 100,
  quantiles.ci = c(0.025,0.975),
  digits = 9,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...) {
  
  if (verbose) {
    cat("#######################################################################\n")
    cat("######################### assigner::fst_WC84 ##########################\n")
    cat("#######################################################################\n")
    timing <- proc.time()
  }
  
  if (snprelate) {
    message("Fst computations with SNPRelate")
    if (ci) {
      message("Confidence Intervals are not implemented with SNPRelate, for now...")
      ci <- FALSE
    }
  } else {
    message("Fst computations with assigner built-in function")
  }
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file is missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # Import data ---------------------------------------------------------------
  if (verbose) message("Importing data")
  input <- stackr::read_long_tidy_wide(data = data, import.metadata = TRUE)
  
  # For long tidy format, switch LOCUS to MARKERS column name, if found MARKERS not found
  if (tibble::has_name(input, "LOCUS") && !tibble::has_name(input, "MARKERS")) {
    input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  }
  
  # population levels and strata------------------------------------------------
  if (!is.null(strata)) {
    if (is.vector(strata)) {
      # message("strata file: yes")
      number.columns.strata <- max(utils::count.fields(strata, sep = "\t"))
      col.types <- stringi::stri_join(rep("c", number.columns.strata), collapse = "")
      suppressMessages(strata.df <- readr::read_tsv(file = strata, col_names = TRUE, col_types = col.types) %>% 
                         dplyr::rename(POP_ID = STRATA))
    } else {
      # message("strata object: yes")
      colnames(strata) <- stringi::stri_replace_all_fixed(
        str = colnames(strata), 
        pattern = "STRATA", 
        replacement = "POP_ID", 
        vectorize_all = FALSE
      )
      strata.df <- strata
    }
    
    # Remove potential whitespace in pop_id
    strata.df$POP_ID <- stringi::stri_replace_all_fixed(strata.df$POP_ID, pattern = " ", replacement = "_", vectorize_all = FALSE)
    
    strata.df$INDIVIDUALS <- stringi::stri_replace_all_fixed(
      str = strata.df$INDIVIDUALS, 
      pattern = c("_", ":"), 
      replacement = c("-", "-"),
      vectorize_all = FALSE
    )
    
    input <- input %>%
      dplyr::select(-POP_ID) %>% 
      dplyr::left_join(strata.df, by = "INDIVIDUALS")
  }
  
  # using pop.levels and pop.labels info if present
  input <- stackr::change_pop_names(data = input, pop.levels = pop.levels, pop.labels = pop.labels)
  
  # SNPRelate prep -------------------------------------------------------------
  if (snprelate) {
    message("Detect if data is biallelic")
    biallelic <- stackr::detect_biallelic_markers(data = input)
    
    if (!biallelic) {
      warning("Data is not biallelic = cannot run assigner with SNPRelate Fst function")
      message("Continuing the Fst computations with built in function...")
      snprelate <- FALSE
    }
  }
  
  if (snprelate) {
    # if holdout set, removes individuals
    if (!is.null(holdout.samples)) {
      message("removing holdout individuals")
      input <- dplyr::filter(.data = input, !INDIVIDUALS %in% holdout.samples)
    }
    
    unique.markers.pop <- input %>% 
      dplyr::filter(!is.na(GT_BIN)) %>%
      dplyr::distinct(MARKERS, POP_ID)
    
  } else {
    # genotyped data and holdout sample
    data.genotyped <- dplyr::select(.data = input, MARKERS, POP_ID, INDIVIDUALS, GT) %>% 
      dplyr::filter(GT != "000000")
    
    # if holdout set, removes individuals
    if (!is.null(holdout.samples)) {
      message("removing holdout individuals")
      data.genotyped <- dplyr::filter(.data = data.genotyped, !INDIVIDUALS %in% holdout.samples)
    }
    
    unique.markers.pop <- dplyr::distinct(.data = data.genotyped, MARKERS, POP_ID)
    
  }
  
  # results stored in this list:
  res <- list()
  
  # Function to compute WC84 Fst ----------------------------------------------
  
  # fst function
  compute_fst <- function(x, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci) {
    # x = data.genotyped # test
    
    # Removing monomorphic markers------------------------------------------------
    # mono.markers <- x %>%
    #   dplyr::select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
    #   tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
    #   tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
    #   dplyr::filter(GT != "000") %>%
    #   dplyr::group_by(MARKERS, GT) %>% 
    #   dplyr::tally(.) %>%
    #   dplyr::ungroup(.) %>% 
    #   dplyr::select(MARKERS) %>% 
    #   dplyr::group_by(MARKERS) %>% 
    #   dplyr::tally(.) %>% 
    #   dplyr::filter(n == 1) %>% 
    #   dplyr::select(MARKERS)
    # 
    # # Remove the markers from the dataset
    # if (length(mono.markers$MARKERS) > 0) {
    #   x <- dplyr::anti_join(x, mono.markers, by = "MARKERS")
    # }
    
    x <- stackr::discard_monomorphic_markers(data = x, verbose = FALSE)$input
    
    # number of marker used for computation 
    n.markers <- dplyr::n_distinct(x$MARKERS)
    
    count.locus <- dplyr::group_by(.data = x, MARKERS) %>%
      dplyr::summarise(
        NPL = n_distinct(POP_ID),# number of populations per locus
        NIL = n() # number of individuals per locus
      )
    
    count.locus.pop <- dplyr::group_by(.data = x, POP_ID, MARKERS) %>%
      dplyr::tally(.) %>%
      dplyr::rename(NIPL = n) %>%
      dplyr::mutate(NIPL_SQ = NIPL^2) %>% 
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(NIPL_SQ_SUM = sum(NIPL_SQ, na.rm = TRUE)) %>%
      dplyr::full_join(count.locus, by = "MARKERS") %>%
      dplyr::mutate(NC = (NIL - NIPL_SQ_SUM/NIL)/(NPL - 1))#correction
    
    count.locus <- NULL
    
    # numbers corrected
    allele.locus <- x %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(GT, 1, 3),
        A2 = stringi::stri_sub(GT, 4,6)
      ) %>% 
      dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>% 
      tidyr::gather(key = ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, POP_ID, MARKERS))
    
    correction <- dplyr::distinct(.data = allele.locus, MARKERS, ALLELES) %>%
      dplyr::full_join(count.locus.pop, by = "MARKERS") %>% 
      dplyr::arrange(MARKERS, ALLELES)
    
    count.locus.pop <- NULL
    
    freq.alleles <- allele.locus %>%
      dplyr::group_by(MARKERS, ALLELES, POP_ID) %>% 
      dplyr::tally(.) %>%
      dplyr::ungroup(.) %>%
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
    
    fst.stats.prep <- x %>%
      dplyr::mutate(
        het = ifelse(stringi::stri_sub(GT, 1, 3) != stringi::stri_sub(GT, 4, 6), 1, 0),
        AL1 = stringi::stri_sub(GT, 1, 3),
        AL2 = stringi::stri_sub(GT, 4, 6)
      ) %>% 
      dplyr::select(-GT) %>%
      tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -c(INDIVIDUALS, MARKERS, POP_ID, het)) %>%
      dplyr::select(-ALLELES_GROUP) %>% 
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
      )
    
    
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
      ) 
    
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
      )
    
    fst.fis.markers <- sigma.loc %>% 
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        FST = round(lsiga/(lsiga + lsigb + lsigw), digits),
        FIS = round(lsigb/(lsigb + lsigw), digits)
      ) %>% 
      dplyr::mutate(FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0))
    
    fst.fis.overall <- dplyr::ungroup(sigma.loc.alleles) %>%
      dplyr::summarise(
        tsiga = sum(siga, na.rm = TRUE),
        tsigb = sum(sigb, na.rm = TRUE),
        tsigw = sum(sigw, na.rm = TRUE)
      ) %>% 
      dplyr::summarise(
        FST = round(tsiga/(tsiga + tsigb + tsigw), digits),
        FIS = round(tsigb/(tsigb + tsigw), digits)
      ) %>% 
      dplyr::mutate(FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0))
    # add new column with number of markers
    fst.fis.overall$N_MARKERS <- n.markers
    
    # Confidence Intervals -----------------------------------------------------
    # over loci for the overall Fst estimate
    if (ci) {
      # the function:
      boot.fst.list <- purrr::map(.x = 1:iteration.ci, .f = boot_ci, sigma.loc.alleles = sigma.loc.alleles)
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
      dplyr::arrange(desc(FST)) %>%
      dplyr::select(MARKERS, FST) %>%
      dplyr::mutate(
        RANKING = seq(from = 1, to = n()),
        QUARTILE = ntile(FST,10)
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
    fst.plot <- ggplot(fst.markers, aes(x = FST, na.rm = T)) +
      geom_histogram(binwidth = 0.01) +
      labs(x = "Fst (overall)") +
      expand_limits(x = 0) +
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold"))
    
    # Results ------------------------------------------------------------------
    res$sigma.loc <- sigma.loc
    res$fst.markers <- fst.markers
    res$fst.ranked <- fst.ranked
    res$fst.overall <- fst.overall
    res$fis.markers <- fis.markers
    res$fis.overall <- fis.overall
    res$fst.plot <- fst.plot
    
    return(res)
  } # End compute_fst function
  
  # Pairwise Fst function
  pairwise_fst <- function(list.pair, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci) {
    
    pop.select <- stringi::stri_join(purrr::flatten(pop.pairwise[list.pair]))
    
    # data.select <- stackr::keep_common_markers(data = data.select) # longer than below
    # common markers
    set1 <- unique.markers.pop %>%
      dplyr::filter(POP_ID == pop.select[1]) %>%
      dplyr::select(MARKERS)
    
    set2 <- unique.markers.pop %>%
      dplyr::filter(POP_ID == pop.select[2]) %>%
      dplyr::select(MARKERS)
    
    common.set <- dplyr::intersect(set1, set2) %>%
      dplyr::arrange(MARKERS)
    
    data.genotyped <- suppressWarnings(dplyr::semi_join(data.genotyped, common.set, by = "MARKERS"))
    
    data.select <- data.genotyped %>% 
      dplyr::filter(POP_ID %in% pop.select) %>% 
      dplyr::mutate(POP_ID = droplevels(x = POP_ID))
    
    fst.select <- compute_fst(x = data.select, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci)
    df.select <- tibble::data_frame(POP1 = pop.select[1], POP2 = pop.select[2])
    df.select <- dplyr::bind_cols(df.select, fst.select$fst.overall) 
    fst.select <- NULL
    return(df.select)
  } # End pairwise_fst
  
  # Pairwise Fst function with SNPRelate
  pairwise_fst_snprelate <- function(pop.pairwise, data, strata, unique.markers.pop) {
    
    strata.df <- dplyr::filter(.data = strata, POP_ID %in% pop.pairwise) %>% # filter the pop
      dplyr::mutate(POP_ID = droplevels(POP_ID)) # remove unnecessary factors
    
    # markers in common between pair of pop
    set1 <- unique.markers.pop %>% 
      dplyr::filter(POP_ID == pop.pairwise[1]) %>% 
      dplyr::select(MARKERS)
    set2 <- unique.markers.pop %>% 
      dplyr::filter(POP_ID == pop.pairwise[2]) %>% 
      dplyr::select(MARKERS)
    common.set <- dplyr::intersect(set1, set2) %>% 
      dplyr::arrange(MARKERS)
    
    fst.snprelate <- SNPRelate::snpgdsFst(
      gdsobj = data,
      population = strata.df$POP_ID, # factors required
      sample.id = strata.df$INDIVIDUALS,
      snp.id = common.set$MARKERS,
      method = "W&C84",
      remove.monosnp = TRUE,
      maf = NaN,
      missing.rate = NaN,
      autosome.only = FALSE,
      with.id = FALSE,
      verbose = FALSE
    )
    return(fst.snprelate)
  }
  
  # Confidence interval function
  boot_ci <- function(x, sigma.loc.alleles){
    
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
  
  # Compute global Fst ---------------------------------------------------------
  if (snprelate) {
    if (verbose) message("Generating GDS format...")
    # keep markers in common
    gds.genotypes <- suppressMessages(stackr::keep_common_markers(data = input))
    
    strata.df <- dplyr::distinct(gds.genotypes, POP_ID, INDIVIDUALS) %>%
      dplyr::mutate(POP_ID = factor(POP_ID))
    
    snp.id <- dplyr::distinct(.data = gds.genotypes, MARKERS) %>%
      dplyr::arrange(MARKERS) %>%
      purrr::flatten_chr(.)
    
    gds.genotypes <- suppressWarnings(
      dplyr::select(.data = gds.genotypes, MARKERS, INDIVIDUALS, GT_BIN) %>%
        dplyr::group_by(MARKERS) %>% 
        tidyr::spread(data = ., key = INDIVIDUALS, value = GT_BIN) %>% 
        dplyr::arrange(MARKERS) %>%
        tibble::column_to_rownames(df = ., var = "MARKERS") %>% 
        data.matrix(.)
    )  
    
    gds.data <- NULL
    gds.data <- SNPRelate::snpgdsCreateGeno(
      gds.fn = "assigner.gds",
      genmat = gds.genotypes,
      sample.id = strata.df$INDIVIDUALS,
      snp.id = snp.id,
      snp.rs.id = NULL,
      snp.chromosome = NULL,
      snp.position = NULL,
      snp.allele = NULL,
      snpfirstdim = TRUE,
      compress.annotation = "ZIP_RA.max",
      compress.geno = "",
      other.vars = NULL
    )
    
    # Compute the global Fst
    if (verbose) message("Computing global fst")
    gds.file.connection <- SNPRelate::snpgdsOpen("assigner.gds")
    fst.snprelate <- SNPRelate::snpgdsFst(
      gdsobj = gds.file.connection,
      population = strata.df$POP_ID,
      sample.id = strata.df$INDIVIDUALS,
      snp.id = NULL,
      method = "W&C84",
      remove.monosnp = TRUE,
      maf = NaN,
      missing.rate = NaN,
      autosome.only = FALSE,
      with.id = FALSE,
      verbose = FALSE
    )$Fst
    
    fst.overall <- tibble::data_frame(FST = round(fst.snprelate, digits)) %>% 
      dplyr::mutate(FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0))
    
    res$fst.overall <- fst.overall
    
    # close SNPRelate connection if no more computation with SNPRelate
    if (!pairwise) {
      SNPRelate::snpgdsClose(gds.file.connection)
    }
    
    
  } else {
    if (verbose) message("Computing global fst")
    res <- compute_fst(x = data.genotyped, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci)
  }
  
  # Compute pairwise Fst -------------------------------------------------------
  if (pairwise) {
    if (verbose) message("Computing paiwise fst")
    
    pop.list <- levels(input$POP_ID) # pop list
    # all combination of populations
    pop.pairwise <- utils::combn(pop.list, 2, simplify = FALSE) 
    
    if (snprelate) {
      fst.all.pop <- purrr::map(
        .x = pop.pairwise,
        .f = pairwise_fst_snprelate,
        data = gds.file.connection,
        strata = strata.df,
        unique.markers.pop = unique.markers.pop
      )
      # Table with Fst
      pairwise.fst <- t(data.frame(pop.pairwise)) %>% 
        tibble::as_data_frame(.) %>%
        dplyr::bind_cols(dplyr::bind_rows(fst.all.pop)) %>% 
        dplyr::rename(POP1 = V1, POP2 = V2, FST = Fst) %>% 
        dplyr::mutate(
          POP1 = factor(POP1, levels = pop.list, ordered = TRUE),
          POP2 = factor(POP2, levels = pop.list, ordered = TRUE),
          FST = dplyr::if_else(FST < 0, true = 0, false = FST, missing = 0),
          FST = round(FST, digits)
        )
      
      # close SNPRelate connection
      SNPRelate::snpgdsClose(gds.file.connection)
      
    } else {
      # Fst for all pairwise populations
      list.pair <- 1:length(pop.pairwise)
      # list.pair <- 5 #  test
      fst.all.pop <- mclapply(
        X = list.pair, 
        FUN = pairwise_fst, 
        mc.preschedule = FALSE, 
        mc.silent = FALSE, 
        mc.cores = parallel.core, 
        ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci
      )
      # Table with Fst
      pairwise.fst <- dplyr::bind_rows(fst.all.pop) %>% 
        dplyr::mutate(
          POP1 = factor(POP1, levels = pop.list, ordered = TRUE),
          POP2 = factor(POP2, levels = pop.list, ordered = TRUE)
        )
    }#End pairwise Fst
    
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
      pairwise.fst.ci.matrix <- "pairwise fst not selected"
    }
    
    
  } else {
    pairwise.fst <- "pairwise fst not selected"
    upper.mat.fst <- "pairwise fst not selected"
    full.mat.fst <- "pairwise fst not selected"
    pairwise.fst.ci.matrix <- "pairwise fst not selected"
  }
  
  # messages -------------------------------------------------------------------
  if (verbose) {
    cat("############################### RESULTS ###############################\n")
    if (ci) {
      message(stringi::stri_join("Fst (overall): ", res$fst.overall$FST, " [", res$fst.overall$CI_LOW, " - ", res$fst.overall$CI_HIGH, "]"))
    } else{
      message(stringi::stri_join("Fst (overall): ", res$fst.overall$FST))
    }
    timing <- proc.time() - timing
    message(stringi::stri_join("Computation time: ", round(timing[[3]]), " sec"))
    cat("#######################################################################\n")
  }
  # Results pairwise -----------------------------------------------------------
  res$pairwise.fst <- pairwise.fst
  res$pairwise.fst.upper.matrix <- upper.mat.fst
  res$pairwise.fst.full.matrix <- full.mat.fst
  res$pairwise.fst.ci.matrix <- pairwise.fst.ci.matrix
  return(res)
}
