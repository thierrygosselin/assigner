# Compute Weir and Cockerham (1984) Fst

#' @name fst_WC84

#' @title A fast implementation of Weir and Cockerham (1984) Fst/Theta 
#' (overall and paiwise estimates)

#' @description The function computes Weir and Cockerham (1984) 
#' Fst for diploid genomes. Both overall and pairwise Fst can be estimated with 
#' confidence intervals based on bootstrap of markers (resampling with replacement). 
#' The function gives identical results \emph{at the 9th decimal} when tested 
#' against \code{\link[hierfstat]{genet.dist}} in \pkg{hierfstat} and with 
#' the Fst computed in \code{Calculate Distances} or \code{Paiwise Differentiation} 
#' options in \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive}, 
#' that uses the Analysis of Molecular Variance 
#' (AMOVA, Excoffier et al., 1992; Michalakis and Excoffier, 1996). 
#' The fastest computation is still 
#' \href{http://www.bentleydrummer.nl/software/software/GenoDive.html}{GenoDive}, 
#' but it does compute confidence intervals. For an R implementation, \code{\link{fst_WC84}} is very fast. 
#' The computations takes advantage of \pkg{dplyr}, \pkg{tidyr}, \pkg{purrr}, 
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
#' 
#' @note \strong{This function is dedicated to Louis Bernatchez!
#' In the hope that your students find this function fast and usefull.}

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
#' pairwise WC84 Fst is calculated between populations. 
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
#' If not selected \code{detectCores()-1} is used as default.

#' @param verbose (logical, optional) \code{verbose = TRUE} to be chatty 
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
#' @import stringi
#' @import dplyr
#' @import utils
#' @import stackr
#' @importFrom purrr map
#' @importFrom data.table fread

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

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1"){
  utils::globalVariables(
    c("GENOTYPE", "FIS", "POP1", "POP2", "tsiga", "tsigb", "tsigw", "CI", 
      "CI_LOW", "CI_HIGH", "N_MARKERS","sigma.loc.alleles")
  )
}


# Fst function: Weir & Cockerham 1984
fst_WC84 <- function(data,
                     pop.levels = NULL, 
                     pop.labels = NULL, 
                     strata = NULL,
                     holdout.samples = NULL,
                     pairwise = FALSE,
                     ci = FALSE,
                     iteration.ci = 100,
                     quantiles.ci = c(0.025,0.975),
                     digits = 9,
                     parallel.core = detectCores()-1,
                     verbose = FALSE,
                     ...) {
  
  if(verbose){
    cat("#######################################################################\n")
    cat("######################### assigner::fst_WC84 ##########################\n")
    cat("#######################################################################\n")
  }
  
  
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) stop("Input file necessary to write the genepop file is missing")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # Import data ---------------------------------------------------------------
  if(verbose) message("Importing data")
  input <- stackr::read_long_tidy_wide(data = data)
  
  # switch LOCUS to MARKERS if found
  if ("LOCUS" %in% colnames(input)) input <- rename(.data = input, MARKERS = LOCUS)
  
  # population levels and strata  ----------------------------------------------
  if (is.null(strata)){ # no strata
    if(is.null(pop.levels)) { # no pop.levels
      if (is.factor(input$POP_ID)) {
        input$POP_ID <- droplevels(x = input$POP_ID)
      } else {
        input$POP_ID <- factor(input$POP_ID)
      }
    } else { # with pop.levels
      input <- input %>%
        mutate( # Make population ready
          POP_ID = factor(
            stri_replace_all_regex(
              POP_ID, 
              stri_paste("^", pop.levels, "$", sep = ""), 
              pop.labels,
              vectorize_all = FALSE), 
            levels = unique(pop.labels), 
            ordered = TRUE
          )
        )
    }
  } else { # Make population ready with the strata provided
    if (is.vector(strata)) {
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA)
    } else {
      strata.df <- strata
    }
    if(is.null(pop.levels)) { # no pop.levels
      input <- input %>%
        select(-POP_ID) %>% 
        mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        left_join(strata.df, by = "INDIVIDUALS") %>% 
        mutate(POP_ID = factor(POP_ID))
    } else { # with pop.levels
      input <- input %>%
        select(-POP_ID) %>% 
        mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        left_join(strata.df, by = "INDIVIDUALS") %>%
        mutate(
          POP_ID = factor(
            stri_replace_all_regex(
              POP_ID, 
              stri_paste("^", pop.levels, "$", sep = ""), 
              pop.labels, 
              vectorize_all = FALSE
            ),
            levels = unique(pop.labels), ordered = TRUE
          )
        )
    }
  }
  
  # Get the number of pop  -----------------------------------------------------
  # pop.number <- n_distinct(input$POP_ID)
  
  # genotyped data and holdout sample  -----------------------------------------
  if (is.null(holdout.samples)) { # use all the individuals
    data.genotyped <- input %>%
      filter(GT != "000000") # Check for df and plink...
  } else { # if holdout set, removes individuals
    message("removing holdout individuals")
    data.genotyped <- input %>%
      filter(GT != "000000") %>% # remove missing genotypes
      # remove supplementary individual before ranking markers with Fst
      filter(!INDIVIDUALS %in% holdout.samples) 
  }
  
  # Function to compute WC84 Fst ----------------------------------------------
  
  # Confidence interval function
  boot_ci <- function(x, sigma.loc.alleles){
    
    markers.list <- sigma.loc.alleles %>% 
      ungroup() %>% 
      distinct(MARKERS) %>% 
      arrange(MARKERS)
    
    subsample.markers <- markers.list %>% 
      sample_n(tbl = ., size = nrow(markers.list), replace = TRUE) %>% 
      arrange(MARKERS)
    
    fst.fis.overall.iterations <- sigma.loc.alleles %>% 
      right_join(subsample.markers, by = "MARKERS") %>% 
      ungroup %>%
      summarise(
        tsiga = sum(siga, na.rm = TRUE),
        tsigb = sum(sigb, na.rm = TRUE),
        tsigw = sum(sigw, na.rm = TRUE)
      ) %>% 
      summarise(
        FST = round(tsiga/(tsiga + tsigb + tsigw), digits),
        FIS = round(tsigb/(tsigb + tsigw), digits)
      ) %>% 
      mutate(ITERATIONS = rep(x, n()))
    return(fst.fis.overall.iterations)
  } # End boot_ci function
  
  # fst function
  compute_fst <- function(x, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci) {
    # x = data.genotyped # test
    
    # Markers in common between all populations ********************************
    pop.number <- n_distinct(x$POP_ID)
    
    pop.filter <- x %>% 
      group_by(MARKERS) %>%
      filter(n_distinct(POP_ID) == pop.number) %>%
      arrange(MARKERS) %>%
      distinct(MARKERS)
    
    # number of marker used for computation 
    n.markers <- n_distinct(pop.filter$MARKERS)
    
    #Filter
    x <- suppressWarnings(x %>% semi_join(pop.filter, by = "MARKERS"))
    
    # ununsed objects
    pop.filter <- NULL
    pop.number <- NULL
    
    # Removing monomorphic markers------------------------------------------------
    mono.markers <- x %>%
      select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
      tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
      tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
      filter(GT != "000") %>%
      group_by(MARKERS, GT) %>% 
      tally %>%
      ungroup() %>% 
      select(MARKERS) %>% 
      group_by(MARKERS) %>% 
      tally %>% 
      filter(n == 1) %>% 
      select(MARKERS)
    
    
    # Remove the markers from the dataset
    if (length(mono.markers$MARKERS) > 0) {
      x <- anti_join(x, mono.markers, by = "MARKERS")
    }
    
    # The similar hierfstat steps ----------------------------------------------
    n.pop.locus <- x %>%
      select(MARKERS, POP_ID) %>%
      group_by(MARKERS) %>%
      distinct(POP_ID, .keep_all = TRUE) %>%
      tally %>%
      rename(npl = n)
    
    ind.count.locus <- x %>%
      group_by(MARKERS) %>%
      tally
    
    ind.count.locus.pop <- x %>%
      group_by(POP_ID, MARKERS) %>%
      tally %>%
      rename(nal = n) %>%
      mutate(
        nal_sq = nal^2,
        N_IND_GENE = nal*2
      )
    
    freq.al.locus <- x %>%
      tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3) %>%
      # separate the genotypes into alleles
      tidyr::gather(key = ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, POP_ID, MARKERS))
    
    # freq.al.locus2 <- tidyr::separate(data = x, col = GT, into = c("A1", "A2"), sep = 3)
    # system.time(
    #   freq.al.locus2 <- data.table::melt.data.table(
    #     data = as.data.table(freq.al.locus2), 
    #     id.vars = c("INDIVIDUALS", "POP_ID", "MARKERS"), 
    #     variable.name = "ALLELES_GROUP",
    #     variable.factor = FALSE,
    #     value.name = "ALLELES"
    #   ) %>% 
    #     as_data_frame()
    # )
    
    # identical(freq.al.locus, freq.al.locus2)
    
    #pop
    freq.al.locus.pop <- suppressWarnings(
      freq.al.locus %>%
        group_by(POP_ID, MARKERS, ALLELES) %>%
        tally %>%
        full_join(ind.count.locus.pop, by = c("POP_ID", "MARKERS")) %>%
        mutate(P = n/N_IND_GENE) %>% # Freq. Allele per pop
        select(POP_ID, MARKERS, ALLELES, P) %>%
        group_by(MARKERS, ALLELES) %>%
        tidyr::spread(data = ., key = POP_ID, value = P) %>%
        tidyr::gather(key = POP_ID, value = P, -c(MARKERS, ALLELES)) %>%
        mutate(P = as.numeric(stri_replace_na(str = P, replacement = 0))) %>%
        full_join(ind.count.locus.pop, by = c("POP_ID", "MARKERS"))
    )    
    
    freq.al.locus.global <- suppressWarnings(
      freq.al.locus %>%
        group_by(MARKERS, ALLELES) %>%
        tally %>%
        full_join(
          ind.count.locus %>%
            rename(N = n), 
          by = "MARKERS"
        ) %>%
        mutate(pb = n/(2*N)) %>% # Global Freq. Allele
        select(MARKERS, ALLELES, pb)
    )    
    
    mean.n.pop.corrected.per.locus <- suppressWarnings(
      ind.count.locus.pop %>%
        group_by(MARKERS) %>%
        summarise(nal_sq_sum = sum(nal_sq, na.rm = TRUE)) %>%
        full_join(ind.count.locus, by = "MARKERS") %>%
        mutate(nal_sq_sum_nt = (n - nal_sq_sum/n)) %>%
        full_join(n.pop.locus, by = "MARKERS") %>%
        mutate(ncal = nal_sq_sum_nt/(npl-1)) %>%
        select(MARKERS, ncal)
    )
    
    ncal <- suppressWarnings(
      freq.al.locus %>%
        # select(MARKERS, ALLELES) %>%
        # group_by(MARKERS, ALLELES) %>%
        distinct(MARKERS, ALLELES) %>%
        full_join(ind.count.locus, by = "MARKERS") %>%
        rename(ntal = n) %>%
        full_join(mean.n.pop.corrected.per.locus, by = "MARKERS")
    )
    
    data.genotyped.het <- x %>%
      mutate(het = ifelse(stri_sub(GT, 1, 3) != stri_sub(GT, 4, 6), 1, 0))
    
    
    fst.stats.prep <- suppressWarnings(
      data.genotyped.het %>%
        group_by(MARKERS, POP_ID) %>%
        # = the number of heterozygote individuals per pop and markers
        summarise(mho = sum(het, na.rm = TRUE)) %>%  
        group_by(MARKERS) %>%
        tidyr::spread(data = ., key = POP_ID, value = mho) %>%
        tidyr::gather(key = POP_ID, value = mho, -MARKERS) %>%
        mutate(mho = as.numeric(stri_replace_na(str = mho, replacement = 0))) %>%
        full_join(freq.al.locus.pop, by = c("POP_ID", "MARKERS")) %>%
        mutate(
          mhom = round(((2 * nal * P - mho)/2), 0),
          dum = nal * (P - 2 * P^2) + mhom
        ) %>%
        group_by(MARKERS, ALLELES) %>%
        full_join(freq.al.locus.global, by = c("MARKERS", "ALLELES")) %>%
        mutate(
          SSi = sum(dum, na.rm = TRUE),
          dum1 = nal * (P - pb)^2
        ) %>%
        group_by(MARKERS, ALLELES) %>%
        mutate(SSP = 2 * sum(dum1, na.rm = TRUE)) %>%
        group_by(MARKERS, POP_ID) %>%
        mutate(SSG = nal * P - mhom) %>%
        group_by(MARKERS, ALLELES) %>%
        full_join(ncal, by = c("MARKERS", "ALLELES")) %>%
        full_join(n.pop.locus, by = "MARKERS") %>%
        rename(ntalb = npl) %>%
        mutate(
          sigw = round(sum(SSG, na.rm = TRUE), 2)/ntal,
          MSP = SSP/(ntalb - 1),
          MSI = SSi/(ntal - ntalb),
          sigb = 0.5 * (MSI - sigw),
          siga = 1/2/ncal * (MSP - MSI)
        )
    )
    
    # variance components of allele frequencies for each allele
    # siga: among populations
    # sigb: among individuals within populations
    # sigw: within individuals
    sigma.loc.alleles <- fst.stats.prep %>%
      group_by(MARKERS, ALLELES) %>% 
      summarise(
        siga = mean(siga, na.rm = TRUE),
        sigb = mean(sigb, na.rm = TRUE),
        sigw = mean(sigw, na.rm = TRUE)
      ) 
    
    # variance components per locus
    # lsiga: among populations
    # lsigb: among individuals within populations
    # lsigw: within individuals
    
    sigma.loc <- sigma.loc.alleles %>% 
      group_by(MARKERS) %>%
      summarise(
        lsiga = round(sum(siga, na.rm = TRUE), digits),
        lsigb = round(sum(sigb, na.rm = TRUE), digits),
        lsigw = round(sum(sigw, na.rm = TRUE), digits)
      )
    
    fst.fis.markers <- sigma.loc %>% 
      group_by(MARKERS) %>%
      summarise(
        FST = round(lsiga/(lsiga + lsigb + lsigw), digits),
        FIS = round(lsigb/(lsigb + lsigw), digits)
      )
    
    fst.fis.overall <-sigma.loc.alleles %>% 
      ungroup %>%
      summarise(
        tsiga = sum(siga, na.rm = TRUE),
        tsigb = sum(sigb, na.rm = TRUE),
        tsigw = sum(sigw, na.rm = TRUE)
      ) %>% 
      summarise(
        FST = round(tsiga/(tsiga + tsigb + tsigw), digits),
        FIS = round(tsigb/(tsigb + tsigw), digits)
      )
    # add new column with number of markers
    fst.fis.overall$N_MARKERS <- n.markers
    
    # Confidence Intervals -----------------------------------------------------
    # over loci for the overall Fst estimate
    if (ci){
      # the function:
      boot.fst.list <- purrr::map(.x = 1:iteration.ci, .f = boot_ci, sigma.loc.alleles = sigma.loc.alleles)
      boot.fst <- bind_rows(boot.fst.list)
      boot.fst.summary <- boot.fst %>% 
        summarise(
          CI_LOW = round(quantile(FST, 
                                  probs = quantiles.ci[1], 
                                  na.rm = TRUE), 
                         digits),
          CI_HIGH = round(quantile(FST, 
                                   probs = quantiles.ci[2], 
                                   na.rm = TRUE), 
                          digits)
        )
    }
    
    # Fst markers  -------------------------------------------------------------
    fst.markers <- fst.fis.markers %>% 
      select(MARKERS, FST) %>% 
      arrange(MARKERS)
    
    # Ranked fst   -------------------------------------------------------------
    fst.ranked <- fst.markers %>%
      arrange(desc(FST)) %>%
      select(MARKERS, FST) %>%
      mutate(
        RANKING = seq(from = 1, to = n()),
        QUARTILE = ntile(FST,10)
      )
    
    # Fst overall  -------------------------------------------------------------
    if (ci){
      fst.overall <- fst.fis.overall %>% 
        select(FST, N_MARKERS) %>% 
        bind_cols(boot.fst.summary)
    } else {
      fst.overall <- fst.fis.overall %>% 
        select(FST, N_MARKERS)
    }
    
    # Fis markers  -------------------------------------------------------------
    fis.markers <- fst.fis.markers %>% 
      select(MARKERS, FIS) %>% 
      arrange(MARKERS)
    
    # Fis overall   ------------------------------------------------------------
    fis.overall <- fst.fis.overall %>% select(FIS, N_MARKERS)
    
    # Plot -----------------------------------------------------------------------
    fst.plot <- ggplot(fst.markers, aes(x = FST, na.rm = T))+
      geom_histogram(binwidth = 0.01)+
      labs(x = "Fst (overall)")+
      expand_limits(x = 0)+
      theme(
        legend.position = "none",
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"),
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.title = element_text(size = 10, family = "Helvetica", face = "bold"),
        legend.text = element_text(size = 10, family = "Helvetica", face = "bold"),
        strip.text.x = element_text(size = 10, family = "Helvetica", face = "bold"))
    
    # Results ------------------------------------------------------------------
    res <- list()
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
    pop.select <- stri_paste(flatten(pop.pairwise[list.pair]))
    data.select <- data.genotyped %>% 
      filter(POP_ID %in% pop.select) %>% 
      mutate(POP_ID = droplevels(x = POP_ID))
    fst.select <- compute_fst(x = data.select, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci)
    # if (ci){
    df.select <- data_frame(POP1 = pop.select[1], POP2 = pop.select[2])
    df.select <- bind_cols(df.select, fst.select$fst.overall) 
    # %>% rename(FST = fst.overall)
    # } else {
    # df.select <- data_frame(POP1 = pop.select[1], 
    # POP2 = pop.select[2], 
    # FST = fst.select$fst.overall
    # )
    # }
    fst.select <- NULL
    return(df.select)
  } # End pairwise_fst
  
  
  # Compute global Fst ---------------------------------------------------------
  if(verbose) message("Computing global fst")
  res <- compute_fst(x = data.genotyped, ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci)
  
  # Compute pairwise Fst -------------------------------------------------------
  if (pairwise) {
    if (verbose) message("Computing paiwise fst")
    pop.list <- levels(input$POP_ID) # pop list
    # all combination of populations
    pop.pairwise <- combn(unique(pop.list), 2, simplify = FALSE) 
    # Fst for all pairwise populations
    list.pair <- 1:length(pop.pairwise)
    # list.pair <- 5 #  test
    fst.all.pop <- parallel::mclapply(
      X = list.pair, 
      FUN = pairwise_fst, 
      mc.preschedule = FALSE, 
      mc.silent = FALSE, 
      mc.cores = parallel.core, 
      ci = ci, iteration.ci = iteration.ci, quantiles.ci = quantiles.ci
    )
    
    # Table with Fst------------------------------------------------------------
    pairwise.fst <- bind_rows(fst.all.pop) %>% 
      mutate(
        POP1 = factor(POP1, levels = pop.list, ordered = TRUE),
        POP2 = factor(POP2, levels = pop.list, ordered = TRUE)
      )
    # Matrix--------------------------------------------------------------------
    upper.mat.fst <- pairwise.fst %>% 
      select (POP1, POP2, FST) %>% 
      tidyr::spread(data = ., POP2, FST, fill = "", drop = FALSE) %>% 
      rename(POP = POP1)
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
    
    if (ci){
      # bind upper and lower diagonal of matrix
      lower.mat.ci <- pairwise.fst %>% 
        select (POP1, POP2, CI_LOW, CI_HIGH) %>% 
        tidyr::unite(data = ., CI, CI_LOW, CI_HIGH, sep = " - ") %>% 
        tidyr::spread(data = ., POP2, CI, fill = "", drop = FALSE) %>% 
        rename(POP = POP1)
      
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
  if(verbose){
    cat("############################### RESULTS ###############################\n")
    if (ci) {
      message(stri_paste("Fst (overall): ", res$fst.overall$FST, " [", res$fst.overall$CI_LOW, " - ", res$fst.overall$CI_HIGH, "]"))
    } else{
      message(stri_paste("Fst (overall): ", res$fst.overall$FST))
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
