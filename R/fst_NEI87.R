# Compute Nei's 1987 Fst

#' @name fst_NEI87

#' @title A fast implementation of Nei's 1987 Fst (overall and paiwise estimates)

#' @description 
#' The function computes Nei's 1987 Fst for diploid genomes. 
#' 2 version are offered: the classic Nei's Gst and Nei's G'st (prime), 
#' that comes with a correction for the bias that stems from sampling a limited 
#' number of populations.
#' Both overall and pairwise Fst can be estimated with 
#' confidence intervals based on bootstrap of markers (resampling with replacement). 
#' The function should give identical results \emph{at the 4th decimal} when tested 
#' against \code{\link[hierfstat]{genet.dist}} in \pkg{hierfstat} and with 
#' the Fst computed in \code{Calculate Distances} or 
#' \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive}.
#' The fastest computation is still 
#' \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive}, 
#' but here, the R solution computes confidence intervals and is very fast. 
#' The computations takes advantage of \pkg{tidyverse} packages, 
#' \pkg{data.table} and \pkg{parallel}.
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
#' \code{\link[stackr]{tidy_wide}}. See details for more info.
#' 
#' \emph{How to get a tidy data frame ?}
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{tidy_genomic_data}} can transform 6 genomic data formats 
#' in a tidy data frame. You can also use this function to filter your dataset using
#' whitelist of markers, blacklist of individuals and genotypes.

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

#' @param strata (optional) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' If a \code{strata} file is specified, the strata file will have
#' precedence over any grouping found input file (\code{data}). 
#' The \code{STRATA} column can be any hierarchical grouping.
#' Default: \code{strata = NULL}.

#' @param holdout.samples (optional) Samples that don't participate in the Fst 
#' computation (supplementary). Data frame with one column \code{INDIVIDUALS}.
#' Default: \code{holdout.samples = NULL}.

#' @param pairwise (logical, optional) With \code{pairwise = TRUE}, the 
#' pairwise NEI87 Fst is calculated between populations. 
#' Default: \code{pairwise = FALSE}.

#' @param ci (logical, optional) Compute bootstrapped confidence intervals. 
#' Default: \code{ci = FALSE}.

#' @param iteration.ci (integer, optional) The number of iterations for 
#' the boostraps (resampling with replacement of markers). 
#' Default: \code{iteration.ci = 100}.

#' @param quantiles.ci (double, optional) 
#' The quantiles for the bootstrapped confidence intervals. 
#' Default: \code{quantiles.ci = c(0.025,0.975)}.

#' @param digits (optional, integer) The number of decimal places to be used in 
#' results.
#' Default: \code{digits = 9}.

#' @param parallel.core (optional) The number of core for parallel computation 
#' of pairwise Fst. 
#' If not selected \code{detectCores() - 1} is used as default.

#' @param verbose (logical, optional) \code{verbose = TRUE} to be chatty 
#' during execution. 
#' Default: \code{verbose = FALSE}.

#' @param ... other parameters passed to the function.

#' @return With pairwise comparison computed, the function returns a list with 
#' 10 objects:
#' \itemize{
#'  \item \code{$fst.markers}: Nei's Gst, Nei's G'st and Jost's D by markers,
#'  \item \code{$fst.ranked}: Nei's Gst, Nei's G'st and Jost's D by markers ranked by Nei's Gst,
#'  \item \code{$fst.overall}: the overall Nei's Gst, Nei's G'st and Jost's D by markers with confidence intervals.
#'  \item \code{$fis.markers}: the fis by markers,
#'  \item \code{$fis.overall}: the mean fis overall markers with confidence intervals and the number of markers,
#'  \item \code{$fst.plot}: the histogram of the overall G'st per markers,
#'  \item \code{$pairwise.fst}: pairwise Nei's Gst, Nei's G'st and Jost's D in long/tidy data frame and the number of markers ,
#'  \item \code{$pairwise.fst.upper.matrix}: the pairwise fst prime in a upper triangle matrix,
#'  \item \code{$pairwise.fst.full.matrix}: the pairwise fst prime matrix (duplicated upper and lower triangle),
#'  \item \code{$pairwise.fst.ci.matrix}: matrix with pairwise fst prime in the upper triangle
#'  and the confidence intervals in the lower triangle.
#' }

#' @details \strong{Input data:}
#'  
#' To discriminate the long from the wide format, 
#' the function \pkg{stackr} \code{\link[stackr]{tidy_wide}} searches 
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
#' @rdname fst_NEI87
#' @importFrom dplyr select distinct n_distinct group_by ungroup rename arrange tally filter if_else mutate summarise left_join inner_join right_join anti_join semi_join full_join summarise_each_ funs summarise_if mutate_if count bind_rows bind_cols ntile desc n
#' @importFrom stackr tidy_wide change_pop_names
#' @importFrom tidyr spread gather unite separate complete nesting
#' @importFrom stringi stri_replace_all_regex stri_sub stri_join
#' @importFrom purrr map
#' @importFrom data.table fread melt.data.table as.data.table
#' @importFrom tibble as_data_frame data_frame
#' @importFrom readr read_tsv
#' @importFrom utils count.fields combn
#' @importFrom stats quantile
#' @importFrom parallel detectCores

#' @examples
#' \dontrun{
#' wombat.fst.pairwise <- fst_NEI87(
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

#' @references Nei M. (1987) Molecular Evolutionary Genetics.
#' Columbia University Press
#' @references Meirmans PG, Van Tienderen PH (2004) genotype and genodive: 
#' two programs for the analysis of genetic diversity of asexual organisms. 
#' Molecular Ecology Notes, 4, 792-794.
#' @references Roesti M, Salzburger W, Berner D. (2012)
#' Uninformative polymorphisms bias genome scans for signatures of selection. 
#' BMC Evol Biol., 12:94. doi:10.1111/j.1365-294X.2012.05509.x

#' @seealso
#' From \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive} manual:
#' \emph{'In general, rather than to test differentiation between all pairs of 
#' populations,
#' it is adviseable to perform an overall test of population differentiation, 
#' possibly using a hierarchical population structure, (see AMOVA)'}
#' 
#' To compute an AMOVA, use \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive}
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

# Fst function: Nei's 1987
fst_NEI87 <- function(
  data,
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
    cat("######################## assigner::fst_NEI87 ##########################\n")
    cat("#######################################################################\n")
  }
  
  message("WARNING: This function is still under testing, use with caution, compare the results with GENODIVE and report bugs")
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file is missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # Import data ---------------------------------------------------------------
  if (verbose) message("Importing data")
  input <- stackr::tidy_wide(data = data)
  
  # Change individuals names containing special character
  input$INDIVIDUALS <- stringi::stri_replace_all_fixed(
    str = input$INDIVIDUALS, 
    pattern = c("_", ":"), 
    replacement = c("-", "-"),
    vectorize_all = FALSE
  )
  # switch LOCUS to MARKERS if found
  if ("LOCUS" %in% colnames(input)) input <- dplyr::rename(.data = input, MARKERS = LOCUS)
  
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
  
  
  # genotyped data and holdout sample  -----------------------------------------
  if (is.null(holdout.samples)) { # use all the individuals
    data.genotyped <- input %>%
      dplyr::filter(GT != "000000") # Check for df and plink...
  } else {# if holdout set, removes individuals
    message("removing holdout individuals")
    data.genotyped <- input %>%
      dplyr::filter(GT != "000000") %>% # remove missing genotypes
      # remove supplementary individual before ranking markers with Fst
      dplyr::filter(!INDIVIDUALS %in% holdout.samples) 
  }
  
  # Function to compute Nei's 1987 Fst ----------------------------------------------
  
  # Confidence interval function
  boot_ci <- function(x, fst.data){
    # x <- 1
    markers.list <- fst.data %>% 
      dplyr::ungroup(.) %>% 
      dplyr::distinct(MARKERS) %>% 
      dplyr::arrange(MARKERS)
    
    subsample.markers <- markers.list %>% 
      sample_n(tbl = ., size = nrow(markers.list), replace = TRUE) %>% 
      dplyr::arrange(MARKERS)
    
    fst.data.overall.iterations <- fst.data %>%
      dplyr::right_join(subsample.markers, by = "MARKERS") %>% 
      dplyr::ungroup(.) %>%
      dplyr::mutate(
        HS = MN / (MN - 1) * (1 - MSP2 - HO / 2 / MN),
        HT = 1 - MP2 + HS / MN / NP - HO / 2 / MN / NP,
        FIS = 1 - HO / HS,
        DST = HT - HS,
        DST_P = NP / (NP - 1) * DST,
        HT_P = HS + DST_P,
        NEI_FST = DST / HT,
        NEI_FST = dplyr::if_else(NEI_FST < 0, true = 0, false = NEI_FST, missing = 0),
        NEI_FST_P = DST_P / HT_P,
        NEI_FST_P = dplyr::if_else(NEI_FST_P < 0, true = 0, false = NEI_FST_P, missing = 0),
        JOST_D = DST_P / (1 - HS)
      ) %>% 
      dplyr::summarise_if(is.numeric, funs(mean(., na.rm = TRUE))) %>% 
      dplyr::mutate(
        NEI_FST = DST / HT,
        NEI_FST = dplyr::if_else(NEI_FST < 0, true = 0, false = NEI_FST, missing = 0),
        NEI_FST_P = DST_P / HT_P,
        NEI_FST_P = dplyr::if_else(NEI_FST_P < 0, true = 0, false = NEI_FST_P, missing = 0),
        FIS = 1 - (HO / HS),
        JOST_D = DST_P / (1 - HS)
      ) %>% 
      dplyr::ungroup(.) %>%
      dplyr::mutate_if(is.numeric, funs( round(x = ., digits = digits))) %>%
      dplyr::select(HO, HS, HT, DST, HT_P, DST_P, NEI_FST, NEI_FST_P, FIS, JOST_D) %>% 
      dplyr::mutate(ITERATIONS = rep(x, n()))
    return(fst.data.overall.iterations)
  } # End boot_ci function
  
  # fst function
  compute_fst <- function(x, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci) {
    # x = data.select
    # x = data.genotyped # test
    # Markers in common between all populations ********************************
    pop.number <- dplyr::n_distinct(x$POP_ID)
    
    pop.filter <- x %>% 
      dplyr::group_by(MARKERS) %>%
      dplyr::filter(dplyr::n_distinct(POP_ID) == pop.number) %>%
      dplyr::arrange(MARKERS) %>%
      dplyr::distinct(MARKERS)
    
    # number of marker used for computation 
    n.markers <- dplyr::n_distinct(pop.filter$MARKERS)
    
    #Filter
    x <- suppressWarnings(dplyr::semi_join(x, pop.filter, by = "MARKERS"))
    
    # ununsed objects
    pop.filter <- NULL
    pop.number <- NULL
    
    # Removing monomorphic markers------------------------------------------------
    mono.markers <- x %>%
      dplyr::select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
      tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
      tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
      dplyr::filter(GT != "000") %>%
      dplyr::group_by(MARKERS, GT) %>% 
      dplyr::tally(.) %>%
      dplyr::ungroup(.) %>% 
      dplyr::select(MARKERS) %>% 
      dplyr::group_by(MARKERS) %>% 
      dplyr::tally(.) %>% 
      dplyr::filter(n == 1) %>% 
      dplyr::select(MARKERS)
    
    
    # Remove the markers from the dataset
    if (length(mono.markers$MARKERS) > 0) {
      x <- dplyr::anti_join(x, mono.markers, by = "MARKERS")
    }
    
    mono.markers <- NULL
    
    # split the data per alleles and melt
    x <- x %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(GT, 1, 3),
        A2 = stringi::stri_sub(GT, 4,6)
      ) %>% 
      dplyr::select(-GT) %>% 
      data.table::as.data.table() %>% 
      data.table::melt.data.table(
        data = ., 
        id.vars = c("MARKERS", "INDIVIDUALS", "POP_ID"), 
        variable.name = "ALLELES",
        variable.factor = FALSE,
        value.name = "GT"
      ) %>% 
      tibble::as_data_frame()
      
    
    # frequency per markers, alleles, pop
    p <- x %>%
      dplyr::group_by(MARKERS, POP_ID) %>%
      dplyr::count(GT) %>% 
      dplyr::mutate(P = n / sum(n)) %>% 
      dplyr::select(-n) %>% 
      dplyr::arrange(MARKERS, POP_ID, GT) %>%  #%>% complete(data = ., POP_ID, nesting(MARKERS, GT), fill = list(P = 0)) %>%
      dplyr::ungroup(.)
    
    # mp: mean frequency per markers
    # mp2: sum of square mean frequency per markers
    mean.p2 <- p %>% 
      tidyr::complete(data = ., POP_ID, tidyr::nesting(MARKERS, GT), fill = list(P = 0)) %>%
      dplyr::group_by(MARKERS, GT) %>% 
      dplyr::summarise(MP = mean(P, na.rm = TRUE)) %>% 
      dplyr::group_by(MARKERS) %>% 
      dplyr::summarise(MP2 = sum(MP^2)) %>% 
      dplyr::ungroup(.)
    
    # msp2 mean frequency per markers per pop
    mean.frequency.markers <- p %>%
      dplyr::group_by(MARKERS, POP_ID) %>% 
      dplyr::summarise(SP2 = sum(P^2)) %>% 
      dplyr::group_by(MARKERS) %>% 
      dplyr::summarise(MSP2 = mean(SP2, na.rm = TRUE)) %>% 
      dplyr::ungroup(.)
    
    # For diploid-------------------------------------------------------------------
    # Mean heterozygosity observed per pop and markers
    # mean heterozygosity across all markers
    mean.het.obs.markers <- x %>%
      dplyr::group_by(POP_ID, MARKERS, INDIVIDUALS) %>% 
      dplyr::mutate(HO = dplyr::if_else(GT[ALLELES == "A1"] != GT[ALLELES == "A2"], 1, 0)) %>% 
      dplyr::group_by(POP_ID, MARKERS) %>% 
      dplyr::summarise(HO = mean(HO)) %>% 
      dplyr::group_by(MARKERS) %>% 
      dplyr::summarise(HO = mean(HO)) %>% 
      dplyr::ungroup(.)
    
    # mn: corrected mean number of individuals per markers
    #n: number of individuals, per pop and markers
    fst.data <- x %>%
      dplyr::group_by(POP_ID, MARKERS) %>%
      dplyr::distinct(INDIVIDUALS) %>% 
      dplyr::tally(.) %>% 
      dplyr::mutate(N_INV = 1 / n) %>% 
      dplyr::ungroup(.) %>% 
      dplyr::group_by(MARKERS) %>% 
      dplyr::summarise(
        NP = sum(!is.na(n)), # number of pop per markers
        MN = NP / sum(N_INV, na.rm = TRUE)
      ) %>% 
      dplyr::ungroup(.) %>% 
      dplyr::distinct(MARKERS, NP, MN) %>% 
      dplyr::full_join(mean.het.obs.markers, by = "MARKERS") %>% 
      dplyr::full_join(mean.frequency.markers, by = "MARKERS") %>%
      dplyr::full_join(mean.p2, by = "MARKERS") %>% 
      dplyr::mutate(
        HS = MN / (MN - 1) * (1 - MSP2 - HO / 2 / MN),#Expected Heterozygosity within populations
        HT = 1 - MP2 + HS / MN / NP - HO / 2 / MN / NP,# Total Gene diversity
        FIS = 1 - HO / HS,
        DST = HT - HS,
        DST_P = NP / (NP - 1) * DST,
        HT_P = HS + DST_P,
        NEI_FST = DST / HT,
        NEI_FST = dplyr::if_else(NEI_FST < 0, true = 0, false = NEI_FST, missing = 0),
        NEI_FST_P = DST_P / HT_P,
        NEI_FST_P = dplyr::if_else(NEI_FST_P < 0, true = 0, false = NEI_FST_P, missing = 0),
        JOST_D = DST_P / (1 - HS)
      )
    
    fst.data.select <- fst.data %>% 
      dplyr::select(MARKERS, HO, HS, HT, DST, HT_P, DST_P, NEI_FST, NEI_FST_P, FIS, JOST_D)
    
    mean.p2 <- mean.het.obs.markers <- mean.frequency.markers <- NULL
    
    overall <- fst.data.select %>% 
      dplyr::summarise_if(is.numeric, funs(mean(., na.rm = TRUE))) %>% 
      dplyr::mutate(
        MARKERS = "OVERALL",
        NEI_FST = DST / HT,
        NEI_FST = dplyr::if_else(NEI_FST < 0, true = 0, false = NEI_FST, missing = 0),
        NEI_FST_P = DST_P / HT_P,
        NEI_FST_P = dplyr::if_else(NEI_FST_P < 0, true = 0, false = NEI_FST_P, missing = 0),
        FIS = 1 - (HO / HS),
        JOST_D = DST_P / (1 - HS)
      ) %>% 
      dplyr::ungroup(.) %>%
      dplyr::mutate_if(is.numeric, dplyr::funs( round(x = ., digits = digits))) %>%
      dplyr::select(MARKERS, HO, HS, HT, DST, HT_P, DST_P, NEI_FST, NEI_FST_P, FIS, JOST_D)
    # add new column with number of markers
    
    overall$N_MARKERS <- n.markers
    
    # Confidence Intervals -----------------------------------------------------
    # over loci for the overall Fst estimate
    if (ci) {
      # the function:
      boot.fst.list <- purrr::map(.x = 1:iteration.ci, .f = boot_ci, fst.data = fst.data)
      boot.fst.summary <- dplyr::bind_rows(boot.fst.list) %>% 
        dplyr::summarise(
          FIS_CI_LOW = stats::quantile(FIS, probs = quantiles.ci[1], na.rm = TRUE),
          FIS_CI_HIGH = stats::quantile(FIS, probs = quantiles.ci[2], na.rm = TRUE),
          DST_CI_LOW = stats::quantile(DST, probs = quantiles.ci[1], na.rm = TRUE),
          DST_CI_HIGH = stats::quantile(DST, probs = quantiles.ci[2], na.rm = TRUE),
          DST_P_CI_LOW = stats::quantile(DST_P, probs = quantiles.ci[1], na.rm = TRUE),
          DST_P_CI_HIGH = stats::quantile(DST_P, probs = quantiles.ci[2], na.rm = TRUE),
          NEI_FST_CI_LOW = stats::quantile(NEI_FST, probs = quantiles.ci[1], na.rm = TRUE),
          NEI_FST_CI_HIGH = stats::quantile(NEI_FST, probs = quantiles.ci[2], na.rm = TRUE),
          NEI_FST_P_CI_LOW = stats::quantile(NEI_FST_P, probs = quantiles.ci[1], na.rm = TRUE),
          NEI_FST_P_CI_HIGH = stats::quantile(NEI_FST_P, probs = quantiles.ci[2], na.rm = TRUE),
          JOST_D_CI_LOW = stats::quantile(JOST_D, probs = quantiles.ci[1], na.rm = TRUE),
          JOST_D_CI_HIGH = stats::quantile(JOST_D, probs = quantiles.ci[2], na.rm = TRUE)
        ) %>% 
        dplyr::mutate_if(is.numeric, funs( round(x = ., digits = digits)))
      fst.data <- NULL
    } else {
      fst.data <- NULL
    }
    
    # Fst markers  -------------------------------------------------------------
    fst.markers <- fst.data.select %>% 
      dplyr::select(MARKERS, NEI_FST, NEI_FST_P, JOST_D) %>% 
      dplyr::arrange(MARKERS)
    
    # Ranked fst   -------------------------------------------------------------
    fst.ranked <- fst.markers %>%
      dplyr::arrange(dplyr::desc(NEI_FST)) %>%
      dplyr::mutate(
        RANKING = seq(from = 1, to = n()),
        QUARTILE = dplyr::ntile(NEI_FST,10)
      )
    
    # Fst overall  -------------------------------------------------------------
    if (ci) {
      fst.overall <- overall %>% 
        dplyr::bind_cols(boot.fst.summary) %>%
        dplyr::select(NEI_FST, NEI_FST_CI_LOW, NEI_FST_CI_HIGH, NEI_FST_P, NEI_FST_P_CI_LOW, NEI_FST_P_CI_HIGH, JOST_D, JOST_D_CI_LOW, JOST_D_CI_HIGH, N_MARKERS)
    } else {
      fst.overall <- overall %>% 
        dplyr::select(NEI_FST, NEI_FST_P, JOST_D, N_MARKERS)
    }
    
    # Fis markers  -------------------------------------------------------------
    fis.markers <- fst.data.select %>% 
      dplyr::select(MARKERS, FIS) %>% 
      dplyr::arrange(MARKERS)
    
    # Fis overall   ------------------------------------------------------------
    fis.overall <- overall %>% dplyr::select(FIS, N_MARKERS)
    if (ci) {
      fis.overall <- overall %>% 
        dplyr::bind_cols(boot.fst.summary) %>%
        dplyr::select(FIS, FIS_CI_LOW, FIS_CI_HIGH, N_MARKERS)
    } else {
      fis.overall <- overall %>% dplyr::select(FIS, N_MARKERS)
    }
    
    # Plot -----------------------------------------------------------------------
    fst.plot <- ggplot(fst.markers, aes(x = NEI_FST_P, na.rm = T)) +
      geom_histogram(binwidth = 0.01) +
      labs(x = "Nei's G'st (overall)") +
      expand_limits(x = 0) +
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold"))
    
    # Results ------------------------------------------------------------------
    res <- list()
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
    # list.pair <- 2
    pop.select <- stringi::stri_join(flatten(pop.pairwise[list.pair]))
    data.select <- data.genotyped %>%
      dplyr::filter(POP_ID %in% pop.select) %>% 
      dplyr::mutate(POP_ID = droplevels(x = POP_ID))
    fst.select <- compute_fst(x = data.select, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci)
    # if (ci){
    df.select <- tibble::data_frame(POP1 = pop.select[1], POP2 = pop.select[2])
    df.select <- dplyr::bind_cols(df.select, fst.select$fst.overall) 
    # %>% dplyr::rename(FST = fst.overall)
    # } else {
    # df.select <- tibble::data_frame(POP1 = pop.select[1], 
    # POP2 = pop.select[2], 
    # FST = fst.select$fst.overall
    # )
    # }
    fst.select <- NULL
    return(df.select)
  } # End pairwise_fst
  
  
  # Compute global Fst ---------------------------------------------------------
  if (verbose) message("Computing global fst")
  res <- compute_fst(x = data.genotyped, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci)
  
  # Compute pairwise Fst -------------------------------------------------------
  if (pairwise) {
    if (verbose) message("Computing paiwise fst")
    pop.list <- levels(input$POP_ID) # pop list
    # all combination of populations
    pop.pairwise <- utils::combn(unique(pop.list), 2, simplify = FALSE) 
    # Fst for all pairwise populations
    list.pair <- 1:length(pop.pairwise)
    # list.pair <- 5 #  test
    
    fst.all.pop <- .assigner_parallel(
      X = list.pair, 
      FUN = pairwise_fst, 
      mc.preschedule = FALSE, 
      mc.silent = FALSE, 
      mc.cores = parallel.core, 
      ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci
    )
    
    # Table with Fst------------------------------------------------------------
    pairwise.fst <- dplyr::bind_rows(fst.all.pop) %>% 
      dplyr::mutate(
        POP1 = factor(POP1, levels = pop.list, ordered = TRUE),
        POP2 = factor(POP2, levels = pop.list, ordered = TRUE)
      )
    # Matrix--------------------------------------------------------------------
    upper.mat.fst <- pairwise.fst %>% 
      dplyr::select(POP1, POP2, NEI_FST_P) %>% 
      tidyr::spread(data = ., POP2, NEI_FST_P, fill = "", drop = FALSE) %>% 
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
        dplyr::select(POP1, POP2, NEI_FST_P_CI_LOW, NEI_FST_P_CI_HIGH) %>% 
        tidyr::unite(data = ., CI, NEI_FST_P_CI_LOW, NEI_FST_P_CI_HIGH, sep = " - ") %>% 
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
      message(stringi::stri_join("Nei's G'st (overall): ", res$fst.overall$NEI_FST_P, " [", res$fst.overall$NEI_FST_P_CI_LOW, " - ", res$fst.overall$NEI_FST_P_CI_HIGH, "]"))
    } else{
      message(stringi::stri_join("Nei's G'st (overall): ", res$fst.overall$NEI_FST_P))
    }
    cat("#######################################################################\n")
  }
  # Results pairwise -----------------------------------------------------------
  res$pairwise.fst <- pairwise.fst
  res$pairwise.fst.upper.matrix <- upper.mat.fst
  res$pairwise.fst.full.matrix <- full.mat.fst
  res$pairwise.fst.ci.matrix <- pairwise.fst.ci.matrix
  return(res)
}
