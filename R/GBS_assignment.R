# Write a gsi_sim file from STACKS VCF file

#' @name GBS_assignment
#' @title Assignment analysis in gsi_sim with GBS data produced by STACKS workflow
#' @description \code{gsi_sim} is a tool for doing and simulating genetic stock
#' identification and developed by Eric C. Anderson.
#' The arguments in the \code{GBS_assignment} function were tailored for the
#' reality of GBS data for assignment analysis while
#' maintaining a reproducible workflow.
#' The input data is a VCF file produced by STACKS. Individuals, populations and
#' markers can be filtered and/or selected in several ways using blacklist,
#' whitelist and other arguments. Map-independent imputation of missing genotype
#' using Random Forest or the most frequent category is also available.
#' Markers can be randomly selected for a classic LOO (Leave-One-Out)
#' assignment or chosen based on ranked Fst for a THL
#' (Training, Holdout, Leave-one-out) assignment analysis.

#' @param vcf.file The VCF file created by STACKS.
#' @param whitelist.markers (optional) A whitelist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by CHROM and/or locus and/or by snp.
#' The whitelist is in the working directory (e.g. "whitelist.txt").
#' de novo CHROM column with 'un' need to be changed to 1. 
#' Default \code{NULL} for no whitelist of markers.

#' @param blacklist.genotype (optional) Useful to erase genotype with below 
#' average quality, e.g. genotype with more than 2 alleles in diploid likely 
#' sequencing errors or genotypes with poor genotype likelihood or coverage. 
#' The blacklist as a minimum of 2 column headers (markers and individuals). 
#' Markers can be 1 column (CHROM or LOCUS or POS), 
#' a combination of 2 (e.g. CHROM and POS or CHROM and LOCUS or LOCUS and POS) or 
#' all 3 (CHROM, LOCUS, POS) The markers columns must be designated: CHROM (character
#' or integer) and/or LOCUS (integer) and/or POS (integer). The id column designated
#' INDIVIDUALS (character) columns header. The blacklist must be in the working 
#' directory (e.g. "blacklist.genotype.txt"). For de novo VCF, CHROM column 
#' with 'un' need to be changed to 1. Default \code{NULL} for no blacklist of 
#' genotypes to erase.

#' @param snp.LD (optional) Minimize linkage disequilibrium (LD) by choosing
#' among these 3 options: \code{"random"} selection, \code{"first"} or
#' \code{"last"} SNP on the same read/haplotype. Default = \code{NULL}.
#' @param common.markers (optional) Logical. Default = \code{FALSE}.
#' With \code{TRUE}, will keep markers genotyped in all the populations.
#' @param maf.local.threshold (double) (optional) Filter local/populations maf with a
#' threshold before conduncting the assignment analysis. Default = \code{NULL}.
#' @param maf.global.threshold (double) (optional) Filter global/overall maf with a
#' threshold before conduncting the assignment analysis. Default = \code{NULL}.
#' @param maf.pop.num.threshold (integer) When maf threshold is used,
#' this argument is for the number of pop required to pass the maf thresholds
#' to keep the locus. Default is \code{maf.pop.num.threshold = 1}
#' @param maf.approach Character. By \code{maf.approach = "SNP"} or by \code{maf.approach = "haplotype"}.
#' The function will consider the SNP or ID/LOCUS/haplotype/read MAF statistics to filter the marker.
#' Default is \code{maf.approach = "SNP"}.
#' @param maf.operator \code{maf.operator = "AND"} or default \code{maf.operator = "OR"}.
#' When filtering over LOCUS or SNP, do you want the local \code{"AND"}
#' global MAF to pass the thresholds, or ... you want the local \code{"OR"}
#' global MAF to pass the thresholds, to keep the marker?

#' @param marker.number (Integer or string of number or "all") Calculations with
#' fixed or subsample of your markers. Default= \code{"all"}.
#' e.g. To test 500, 1000, 2000 and all  the markers:
#' \code{marker.number = c(500, 1000, 2000, "all"}.
#' To use only 500 makers \code{marker.number = 500}.
#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").

#' @param sampling.method (character) Should the markers be randomly selected
#' \code{"random"} for a classic Leave-One-Out (LOO) assignment or
#' chosen based on ranked Fst \code{"ranked"}, used in a
#' Training-Holdout-Leave One Out (THL) assignment ?
#' @param THL (character, integer, proportion) For \code{sampling.method = "ranked"} only.
#' Default \code{1}, 1 individual sample is used as holdout. This individual is not
#' participating in the markers ranking. For each marker number,
#' the analysis will be repeated with all the indiviuals in the data set
#' (e.g. 500 individuals, 500 times 500, 1000, 2000 markers).
#' If a proportion is used e.g. \code{0.15},= 15% of individuals in each
#' populations are chosen randomly as holdout individuals.
#' With \code{THL = "all"} all individuals are used for ranking (not good) and
#' \code{iterations} argument below is set to \code{1} by default.
#' For the other THL values, you can create different holdout individuals lists
#' with the \code{iterations} argument below.
#' @param iterations With random marker selection the iterations argument =
#' the number of iterations to repeat marker resampling, default is \code{10}
#' With \code{marker.number = c(500, 1000)} and default iterations setting,
#' 500 markers will be randomly chosen 10 times and 1000 markers will be randomly
#' chosen 10 times. For the ranked method, using \code{THL = 1}, the analysis
#' will be repeated for each individuals in the data set for every
#' \code{marker.number} selected. With a proportion argument \code{THL = 0.15},
#' 15% of individuals in each populations are chosen randomly as holdout
#' individuals and this process is reapeated the number of times chosen by the
#' \code{iterations} value.


#' @param folder (optional) The name of the folder created in the working directory to save the files/results.
#' @param gsi_sim.filename (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. Default \code{gsi_sim_data.txt}.
#' The number of markers used will be appended to the name of the file.
#' @param keep.gsi.files (Boolean) Default \code{FALSE} The input and output gsi_sim files
#' will be deleted from the directory when finished processing.
#' With \code{TRUE}, remember to allocate a large chunk of the disk space for the analysis.
#' @param pop.levels (required) A character string with your populations ordered.
#' @param pop.labels (optional) A character string for your populations labels.
#' If you need to rename sampling sites in \code{pop.levels} or combined sites/pop
#' into a different names, here is the place.
#' @param pop.id.start The start of your population id
#' in the name of your individual sample.
#' @param pop.id.end The end of your population id
#' in the name of your individual sample.

#' @param pop.select (string) Conduct the assignment analysis on a
#' selected list of populations. Default = \code{NULL} for no selection and keep
#' all population.
#' e.g. \code{pop.select = "QUE"} to select QUE population samples.
#' \code{pop.select = c("QUE", "ONT")} to select QUE and ONT population samples.

#' @param subsample (Integer or Proportion) Default is no sumsampling, \code{subsample = NULL}.
#' With a proportion argument \code{subsample = 0.15}, 15 percent of individuals
#' in each populations are chosen randomly to represent the dataset.
#' With \code{subsample = 36}, 36 individuals in each populations are chosen
#' randomly to represent the dataset.
#' @param iterations.subsample (Integer) The number of iterations to repeat 
#' subsampling, default: \code{iterations.subsample = 1}.
#' With \code{subsample = 20} and \code{iterations.subsample = 10},
#' 20 individuals/populations will be randomly chosen 10 times.

#' @param baseline (optional) A character string with your baseline id.
#' From the \code{pop.id.start} and \code{pop.id.end} you isolate
#' the baseline and mixture group. Here you need to give the id for
#' baseline e.g. \code{c("QUE-ADU", "ONT-ADU")}.
#' @param mixture (optional) But required if bseline was selected. A character
#' string with your mixture id. e.g. \code{c("QUE-JUV", "ONT-JUV")}.

#' @param imputations Should a map-independent imputations of markers be
#' computed. Available choices are: (1) \code{FALSE} for no imputation.
#' (2) \code{"max"} to use the most frequent category for imputations.
#'  (3) \code{"rf"} using Random Forest algorithm. Default = \code{FALSE}.
#' @param imputations.group \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by populations. If you choose
#' global, turn the verbose to \code{TRUE}, to see progress.
#' Default = \code{"populations"}.
#' @param num.tree The number of trees to grow in Random Forest. Default is 100.
#' @param iteration.rf The number of iterations of missing data algorithm
#' in Random Forest. Default is 10.
#' @param split.number Non-negative integer value used to specify
#' random splitting in Random Forest. Default is 100.
#' @param verbose Logical. Should trace output be enabled on each iteration
#' in Random Forest ? Default is \code{FALSE}.
#' @param parallel.core (optional) The number of core for OpenMP shared-memory parallel
#' programming of Random Forest imputations. For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.
#' @details The imputations using Random Forest requires more time to compute
#' and can take several
#' minutes and hours depending on the size of the dataset and polymorphism of
#' the species used. e.g. with a low polymorphic taxa, and a data set
#' containing 30\% missing data, 5 000 haplotypes loci and 500 individuals
#' will require 15 min.
#' The Fst is based on Weir and Cockerham 1984 equations.
#' @return Depending on arguments selected, several files are written to the your
#' working directory or \code{folder}
#' The output in your global environment is a list. To view the assignment results
#' \code{$assignment} to view the ggplot2 figure \code{$plot.assignment}. 
#' See example below.

#' @note \code{GBS_assignment} assumes that the command line version of gsi_sim 
#' is properly installed and available on the command line, so it is executable from 
#' any directory (more info on how to do this, here 
#' \url{http://gbs-cloud-tutorial.readthedocs.org/en/latest/03_computer_setup.html?highlight=bash_profile#save-time}.
#' The easiest way is to put the binary, the \code{gsi_sim} executable,
#' in the folder \code{/usr/local/bin}. To compile gsi_sim, follow the 
#' instruction here: \url{https://github.com/eriqande/gsi_sim}.

#' @export
#' @rdname GBS_assignment
#' @import dplyr
#' @import foreach
#' @import parallel
#' @import doParallel
#' @import stringi
#' @importFrom purrr map
#' @importFrom purrr flatten

#' @examples
#' \dontrun{
#' assignment.treefrog <- GBS_assignment(
#' vcf.file = "batch_1.vcf",
#' whitelist.markers = "whitelist.vcf.txt",
#' snp.LD = NULL,
#' common.markers = TRUE,
#' marker.number = c(500, 5000, "all"),
#' sampling.method = "ranked",
#' THL = 0.3,
#' blacklist.id = "blacklist.id.lobster.tsv",
#' subsample = 25,
#' iterations.subsample = 10
#' gsi_sim.filename = "treefrog.txt",
#' keep.gsi.files = FALSE,
#' pop.levels = c("PAN", "COS")
#' pop.id.start = 5, pop.id.end = 7,
#' imputations = FALSE,
#' parallel.core = 12
#' )
#' 
#' Since the 'folder' argument is missing, it will be created automatically
#' inside your working directory.
#' 
#' To create a dataframe with the assignment results: 
#' assignment <- assignment.treefrog$assignment.
#' 
#' To plot the assignment using ggplot2 and facet 
#' (with subsample by current pop):
#' assignment.treefrog$plot.assignment + facet_grid(SUBSAMPLE~CURRENT).
#' 
#' To save the plot:
#' ggsave("assignment.treefrog.THL.subsample.pdf", height = 35, 
#' width = 60,dpi = 600, units = "cm", useDingbats = F)
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
#'  Regression and Classification (RF-SRC), R package version 1.6.1.
#' @references Ishwaran H. and Kogalur U.B. (2007). Random survival forests
#' for R. R News 7(2), 25-31.
#' @references Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S. (2008).
#' Random survival forests. Ann. Appl. Statist. 2(3), 841--860.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1"){
  utils::globalVariables(
    c('ID', '#CHROM', 'CHROM', 'FORMAT', 'INDIVIDUALS', 'FORMAT_ID', 'LOCUS',
      'POS', 'REF', 'ALT', 'POP_ID', 'READ_DEPTH', 'ALLELE_DEPTH', 'GL',
      'ERASE', 'GT', 'MARKERS', 'QQ', 'PQ', 'N', 'MAF_GLOBAL', 'MAF_LOCAL',
      'ALLELES', 'POP_ID', 'GT', 'INDIVIDUALS', 'MARKERS', 'POP_ID', 'nal',
      'ALLELES_GROUP', 'ALLELES', 'N_IND_GENE', 'P', 'N', 'nal_sq',
      'nal_sq_sum', 'nal_sq_sum_nt', 'npl', 'het', 'mho', 'mhom', 'dum',
      'dum1', 'SSG', 'ntal', 'SSP', 'ntalb', 'SSi', 'MSI', 'sigw', 'MSP',
      'siga', 'sigb', 'lsiga', 'lsigb', 'lsigw', 'FST', 'MARKERS',
      'MARKERS_ALLELES', 'ALLELES', 'POP_ID', 'INDIVIDUALS', 'filename',
      'ID', 'KEEPER', 'ASSIGN', 'OTHERS', 'CURRENT', 'INFERRED',
      'SECOND_BEST_POP', 'SCORE', 'SECOND_BEST_SCORE',
      'MARKER_NUMBER', 'MISSING_DATA', 'TOTAL', 'ASSIGNMENT_PERC',
      'MARKERS', 'CURRENT', 'INFERRED', 'MARKER_NUMBER', 'MISSING_DATA',
      'ITERATIONS', 'METHOD', 'TOTAL', 'MEAN_i', 'MEAN', 'ASSIGNMENT_PERC',
      'SE', 'MEDIAN', 'MIN', 'MAX', 'QUANTILE25', 'QUANTILE75', 'SE_MIN',
      'SE_MAX', '.', 'QUAL', 'FILTER', 'INFO', 'pb', 'SUBSAMPLE'
    )
  )
}

GBS_assignment <- function(vcf.file,
                           whitelist.markers = NULL,
                           blacklist.genotype = NULL,
                           snp.LD = NULL,
                           common.markers = NULL,
                           maf.local.threshold = NULL,
                           maf.global.threshold = NULL,
                           maf.pop.num.threshold = 1,
                           maf.approach = "SNP",
                           maf.operator = "OR",
                           marker.number = "all",
                           blacklist.id = NULL,
                           pop.levels,
                           pop.labels,
                           pop.id.start, 
                           pop.id.end,
                           pop.select = NULL,
                           subsample = NULL,
                           iterations.subsample = 1,
                           sampling.method,
                           THL = 1,
                           iterations = 10,
                           folder,
                           gsi_sim.filename = "gsi_sim_data.txt",
                           keep.gsi.files,
                           baseline = NULL,
                           mixture = NULL,
                           imputations = FALSE,
                           imputations.group = "populations",
                           num.tree = 100,
                           iteration.rf = 10,
                           split.number = 100,
                           verbose = FALSE,
                           parallel.core = NULL) {
  
  message("Assignment analysis using stackr and gsi_sim")
  
  # Checking for missing and/or default arguments ******************************
  if (missing(vcf.file)) stop("VCF file required")
  if (missing(whitelist.markers)) whitelist.markers <- NULL # no Whitelist
  if (missing(blacklist.genotype)) blacklist.genotype <- NULL # no genotype to erase
  if (missing(snp.LD)) snp.LD <- NULL
  if (missing(common.markers)) common.markers <- FALSE
  if (missing(maf.local.threshold)) maf.local.threshold <- NULL
  if (missing(maf.global.threshold)) maf.global.threshold <- NULL
  if (missing(maf.pop.num.threshold)) maf.pop.num.threshold <- 1
  if (missing(maf.approach)) maf.approach <- "SNP"
  if (missing(maf.operator)) maf.operator <- "OR"
  if (missing(marker.number)) marker.number <- "all"
  if (missing(blacklist.id)) blacklist.id <- NULL # No blacklist of ID
  if (missing(pop.levels)) stop("pop.levels required")
  if (missing(pop.labels)) pop.labels <- pop.levels # pop.labels
  if (missing(pop.id.start)) stop("pop.id.start required")
  if (missing(pop.id.end)) stop("pop.id.end required")
  if (missing(pop.select)) pop.select <- NULL
  if (missing(subsample)) subsample <- NULL
  if (missing(iterations.subsample)) iterations.subsample <- 1
  if (missing(sampling.method)) stop("Sampling method required")
  if (sampling.method == "ranked" & missing(THL)) THL <- 1 # THL
  if (missing(iterations)) iterations <- 10
  if (THL == "all") iterations <- 1
  if (THL == 1) iterations <- 1
  if (missing(gsi_sim.filename)) gsi_sim.filename <- "gsi_sim_data.txt"
  if (missing(keep.gsi.files)) keep.gsi.files <- FALSE
  if (missing(baseline)) baseline <- NULL # Baseline
  if (missing(mixture)) mixture <- NULL  # Mixture
  if (missing(imputations)) imputations <- FALSE
  if (missing(imputations.group)) imputations.group <- "populations"
  if (missing(num.tree)) num.tree <- 100
  if (missing(iteration.rf)) iteration.rf <- 10
  if (missing(split.number)) split.number <- 100
  if (missing(verbose)) verbose <- FALSE
  if (missing(parallel.core) | is.null(parallel.core)) parallel.core <- detectCores()-1
  
  # Create a folder based on filename to save the output files *****************
  if (missing(folder)){
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    if (imputations == "FALSE") {
      message("Imputations: not selected")
      directory <- stri_join(getwd(),"/", "assignment_analysis_", "method_", sampling.method, "_no_imputations_", file.date, "/", sep = "")
      dir.create(file.path(directory))
    } else {
      message("Imputations: selected")
      directory <- stri_join(getwd(),"/","assignment_analysis_", "method_", sampling.method, "_imputations_", imputations,"_", imputations.group, "_", file.date, "/", sep = "")
      dir.create(file.path(directory))
    }
    message(stri_join("Folder: ", directory))
    file.data <- NULL #unused object
  } else {
    directory <- stri_join(getwd(), "/", folder, "/", sep = "")
    dir.create(file.path(directory))
    message(stri_join("Folder: ", directory))
  }
  
  # Import/read VCF ************************************************************
  message("Importing the VCF...")
  vcf <- read_delim(
    vcf.file,
    delim = "\t",
    comment = "##",
    progress = interactive()
  ) %>%
    select(-c(QUAL, FILTER, INFO)) %>%
    rename(LOCUS = ID, CHROM = `#CHROM`) %>%
    mutate(
      CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1")
    )
  
  # Detect STACKS version ******************************************************
  if (stri_detect_fixed(vcf$FORMAT[1], "AD")) {
    stacks.version <- "new"
  } else{
    stacks.version <- "old"
  }
  vcf <- vcf %>% select(-FORMAT)
  
  # Whitelist of markers *******************************************************
  if (is.null(whitelist.markers)) { # no Whitelist
    message("No whitelist to apply to the VCF")
    vcf <- vcf
  } else { # with Whitelist of markers
    message("Filtering the VCF with the whitelist from your directory")
    whitelist.markers <- read_tsv(whitelist.markers, col_names = TRUE)
    columns.names.whitelist <- colnames(whitelist.markers)
    if ("CHROM" %in% columns.names.whitelist){
      whitelist.markers$CHROM <- as.character(whitelist.markers$CHROM)
    }
    vcf <- suppressWarnings(
      vcf %>%
        semi_join(whitelist.markers, by = columns.names.whitelist)
    )
  }
  
  # Tidying the VCF to make it easy to work on the data for conversion *********
  message("Making the VCF population wise")
  vcf <- vcf %>%
    tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(CHROM, LOCUS, POS, REF, ALT)) %>% # Gather individuals in 1 colummn
    mutate( # Make population ready
      POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
      POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = F), levels = unique(pop.labels), ordered =T),
      INDIVIDUALS =  as.character(INDIVIDUALS)
    )
  
  # Blacklist id ***************************************************************
  if (is.null(blacklist.id)) { # No blacklist of ID
    message("No individual blacklisted")
    vcf <- vcf
  } else { # With blacklist of ID
    message("Using the blacklisted id from the directory")
    blacklist.id <- read_tsv(blacklist.id, col_names = T)
    vcf <- suppressWarnings(
      vcf %>%
        anti_join(blacklist.id, by = "INDIVIDUALS") %>%
        mutate(POP_ID = droplevels(POP_ID))
    )
  }
  
  # Pop select *****************************************************************
  if (is.null(pop.select)){
    vcf <- vcf
  } else {
    message(stri_join(length(pop.select), "population(s) selected", sep = " "))
    vcf <- suppressWarnings(
      vcf %>%
        filter(POP_ID %in% pop.select)
    )
  }
  
  # Tidy VCF *******************************************************************
  message("Tidy vcf into factory for conversion into gsi_sim ...")
  if (stacks.version == "new"){ # with new version of stacks > v.1.29
    vcf <- vcf %>%
      tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"),
                      sep = ":", extra = "warn") %>%
      select(-c(READ_DEPTH, ALLELE_DEPTH, GL))
  } else { # stacks version prior to v.1.29 had no Allele Depth field...
    vcf <- vcf %>%
      tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "GL"),
                      sep = ":", extra = "warn") %>%
      select(-c(READ_DEPTH, GL))
  }
  
  # Blacklist genotypes ********************************************************
  if (is.null(blacklist.genotype)) { # no Whitelist
    message("No genotype to erase")
    vcf <- vcf
  } else {
    message("Erasing genotype with the blacklist")
    blacklist.genotype <- read_tsv(blacklist.genotype, col_names = TRUE)
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    if ("CHROM" %in% columns.names.blacklist.genotype){
      columns.names.blacklist.genotype$CHROM <- as.character(columns.names.blacklist.genotype$CHROM)
    }
    
    # control check to keep only whitelisted markers from the blacklist of genotypes
    if (!is.null(whitelist.markers)){
      message("Control check to keep only whitelisted markers 
              present in the blacklist of genotypes to erase.")
      # updating the whitelist of markers t have all columns that id markers
      whitelist.markers.ind <- vcf %>% select(CHROM, LOCUS, POS, INDIVIDUALS) %>% distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      # updating the blacklist.genotype
      blacklist.genotype <- suppressWarnings(semi_join(whitelist.markers.ind, blacklist.genotype, by = columns.names.blacklist.genotype))
    } else {
      blacklist.genotype <- blacklist.genotype
    }
    
    # control check to remove blacklisted individuals from the blacklist of genotypes
    if (!is.null(blacklist.id)){
      message("Control check to remove blacklisted individuals 
              present in the blacklist of genotypes to erase.")
      blacklist.genotype <- suppressWarnings(anti_join(blacklist.genotype, blacklist.id, by = "INDIVIDUALS"))
    } else {
      blacklist.genotype <- blacklist.genotype
    }
    
    # Add one column that will allow to include the blacklist in the dataset 
    # by x column(s) of markers
    blacklist.genotype <- mutate(.data = blacklist.genotype, ERASE = rep("erase", n()))
    
    vcf <- suppressWarnings(
      vcf %>%
        full_join(blacklist.genotype, by = c("CHROM", "LOCUS", "POS", "INDIVIDUALS")) %>%
        mutate(
          ERASE = stri_replace_na(str = ERASE, replacement = "ok"),
          GT = ifelse(ERASE == "erase", "./.", GT)
        ) %>% 
        select(-ERASE)
    )
  } # end erase genotypes
  
  # dump unused object
  blacklist.id <- NULL
  whitelist.markers <- NULL
  whitelist.markers.ind <- NULL
  
  # subsampling data ***********************************************************
  # Function:
  subsampling_data <- function(iterations.subsample, ...){
    # message(paste0("Creating data subsample: ", iterations.subsample))
    if (is.null(subsample)){
      subsample.select <- ind.pop.df %>% 
        mutate(SUBSAMPLE = rep(iterations.subsample, n()))
    } else{
      if (subsample > 1){ # integer
        subsample.select <- ind.pop.df %>%
          group_by(POP_ID) %>%
          sample_n(subsample, replace = FALSE) %>% # sampling individuals for each pop
          arrange(POP_ID, INDIVIDUALS) %>% 
          mutate(SUBSAMPLE = rep(iterations.subsample, n())) %>% 
          ungroup()
      }
      if (subsample < 1){ # proportion
        subsample.select <- ind.pop.df %>%
          group_by(POP_ID) %>%
          sample_frac(subsample, replace = FALSE) %>% # sampling individuals for each pop
          arrange(POP_ID, INDIVIDUALS) %>% 
          mutate(SUBSAMPLE = rep(iterations.subsample, n())) %>% 
          ungroup()
      }
    }
    return(subsample.select)
  } # end subsampling function
  
  # subsample <- 15 # test
  # subsample <- NULL # test
  # iterations.subsample <- 5 # test
  # iterations.subsample <- 1 # test
  
  # create the subsampling list
  ind.pop.df <- vcf %>% select(POP_ID, INDIVIDUALS) %>% distinct(POP_ID, INDIVIDUALS)
  subsample.list <- map(.x = 1:iterations.subsample, .f = subsampling_data, subsample = subsample)
  
  # keep track of subsampling individuals and write to directory
  if (is.null(subsample)){
    message("Subsampling: not selected")
  } else{
    message("Subsampling: selected")
    subsampling.individuals <- bind_rows(subsample.list)
    write_tsv(x = subsampling.individuals, path = paste0(directory, "subsampling.individuals.tsv"), col_names = TRUE, append = FALSE)
  } # end subsampling
  
  # unused objects
  subsampling.individuals <- NULL
  ind.pop.df <- NULL
  
  # assignment analysis ********************************************************
  # Function:
  assignment_function <- function(data, ...){
    # data <- subsample.list[[1]] # test
    subsampling.individuals <- data
    subsample.id <- unique(subsampling.individuals$SUBSAMPLE)
    
    if (!is.null(subsample)){
      message(paste("Analyzing subsample: ", subsample.id))
    }
    
    # Updating directories for subsampling
    if(is.null(subsample)){
      directory.subsample <- directory
    } else {
      directory.subsample <- paste0(directory, "subsample_", subsample.id, "/")
      dir.create(file.path(directory.subsample))
    }
    
    vcf <- vcf %>%
      semi_join(subsampling.individuals, by = c("POP_ID", "INDIVIDUALS"))
    
    # unused object
    data <- NULL
    subsampling.individuals <- NULL
    
    # LD control... keep only 1 SNP per haplotypes/reads (optional) ************
    if (is.null(snp.LD)){
      vcf <- vcf
      snp.LD <- NULL
    } else{
      message("Minimizing LD...")
      snp.locus <- vcf %>% select(LOCUS, POS) %>% distinct(POS)
      # Random selection
      if (snp.LD == "random"){
        snp.select <- snp.locus %>%
          group_by(LOCUS) %>%
          sample_n(size = 1, replace = FALSE)
        message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP randomly selected to keep 1 SNP per read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
      }
      
      # Fist SNP on the read
      if (snp.LD == "first"){
        snp.select <- snp.locus %>%
          group_by(LOCUS) %>%
          summarise(POS = min(POS))
        message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
      }
      
      # Last SNP on the read
      if (snp.LD == "last"){
        snp.select <- snp.locus %>%
          group_by(LOCUS) %>%
          summarise(POS = max(POS))
        message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
      }
      
      # filtering the VCF to minimize LD
      vcf <- vcf %>% semi_join(snp.select, by = c("LOCUS", "POS"))
      message("Filtering the tidy VCF to minimize LD by keeping only 1 SNP per short read/haplotype")
    } # end of snp.LD control
    
    # Unique markers id: combine CHROM, LOCUS and POS into MARKERS *************
    vcf <- vcf %>%
      mutate(
        POS = stri_pad_left(str = POS, width = 8, pad = "0"),
        LOCUS = stri_pad_left(str = LOCUS, width = 8, pad = "0")
      ) %>%
      arrange(CHROM, LOCUS, POS) %>%
      tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_")
    
    # Markers in common between all populations (optional) *********************
    if (common.markers == FALSE) {
      vcf <- vcf
    } else { # keep only markers present in all pop
      message("Using markers common in all populations")
      pop.number <- n_distinct(vcf$POP_ID)
      
      pop.filter <- vcf %>%
        filter(GT != "./.") %>%
        group_by(MARKERS) %>%
        filter(n_distinct(POP_ID) == pop.number) %>%
        arrange(MARKERS) %>%
        select(MARKERS) %>%
        distinct(MARKERS)
      
      message(stri_join("Number of original markers = ", n_distinct(vcf$MARKERS), 
                     "\n", "Number of markers present in all the populations = ", 
                     n_distinct(pop.filter$MARKERS), "\n", 
                     "Number of markers removed = ", 
                     n_distinct(vcf$MARKERS) - n_distinct(pop.filter$MARKERS))
      )
      vcf <- suppressWarnings(vcf %>% semi_join(pop.filter, by = "MARKERS"))
      pop.filter <- NULL # ununsed object
    } # end common markers
    
    # Minor Allele Frequency filter ********************************************
    if (is.null(maf.global.threshold) | is.null(maf.local.threshold)){ # no MAF
      vcf <- vcf
    } else { # with MAF
      message("Filtering the VCF with MAF")
      maf.local <- vcf %>%
        filter(GT != "./.") %>%
        group_by(MARKERS, POP_ID, REF, ALT) %>%
        summarise(
          N = as.numeric(n()),
          PQ = as.numeric(length(GT[GT == "1/0" | GT == "0/1"])),
          QQ = as.numeric(length(GT[GT == "1/1"]))
        ) %>%
        mutate(MAF_LOCAL = ((QQ * 2) + PQ) / (2 * N))
      
      maf.global <- maf.local %>%
        group_by(MARKERS) %>%
        summarise_each_(funs(sum), vars = c("N", "PQ", "QQ")) %>%
        mutate(MAF_GLOBAL = ((QQ * 2) + PQ) / (2 * N)) %>%
        select(MARKERS, MAF_GLOBAL)
      
      maf.data <- maf.global %>%
        left_join(maf.local, by = c("MARKERS")) %>%
        select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
      
      write_tsv(x = maf.data, path = paste0(directory.subsample,"maf.data.tsv"), col_names = TRUE, append = FALSE)
      message("The MAF table was written in your folder")
      
      # update the vcf with the maf info
      vcf <- full_join(vcf, maf.data, by = c("MARKERS", "POP_ID"))
      # vcf.bk <- vcf               # for test
      # vcf <- vcf.bk               # for test
      # maf.local.threshold <- 0.05 # for test
      # maf.global.threshold <- 0.1 # for test
      # maf.pop.num.threshold <- 1  # for test
      
      if (maf.approach == "haplotype"){
        vcf.maf <- tidyr::separate(data = vcf, col = MARKERS, into = c("CHROM", "LOCUS", "POS"), sep = "_", remove = FALSE, extra = "warn")
        
        if (maf.operator == "OR") {
          vcf.maf <- vcf %>%
            group_by(LOCUS, POP_ID) %>%
            summarise(
              MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
              MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
            ) %>%
            filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
            group_by(LOCUS) %>%
            tally() %>%
            filter(n >= maf.pop.num.threshold) %>%
            select(LOCUS) %>%
            left_join(vcf, by = "LOCUS") %>%
            arrange(LOCUS, POP_ID)
        } else { # AND operator between local and global maf
          vcf.maf <- vcf %>%
            group_by(LOCUS, POP_ID) %>%
            summarise(
              MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
              MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
            ) %>%
            filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
            group_by(LOCUS) %>%
            tally() %>%
            filter(n >= maf.pop.num.threshold) %>%
            select(LOCUS) %>%
            left_join(vcf, by = "LOCUS") %>%
            arrange(LOCUS, POP_ID)
        }
        vcf.maf <- vcf %>% select(-c(CHROM, LOCUS, POS))
      } else { # SNP approach
        if(maf.operator == "OR") {
          vcf.maf <- vcf %>%
            group_by(MARKERS, POP_ID) %>%
            summarise(
              MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
              MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
            ) %>%
            filter(MAF_LOCAL >= maf.local.threshold | MAF_GLOBAL >= maf.global.threshold) %>%
            group_by(MARKERS) %>%
            tally() %>%
            filter(n >= maf.pop.num.threshold) %>%
            select(MARKERS) %>%
            left_join(vcf, by = "MARKERS") %>%
            arrange(MARKERS, POP_ID)
        } else { # AND operator between local and global maf
          vcf.maf <- vcf %>%
            group_by(MARKERS, POP_ID) %>%
            summarise(
              MAF_GLOBAL = mean(MAF_GLOBAL, na.rm = TRUE),
              MAF_LOCAL = min(MAF_LOCAL, na.rm = TRUE)
            ) %>%
            filter(MAF_LOCAL >= maf.local.threshold & MAF_GLOBAL >= maf.global.threshold) %>%
            group_by(MARKERS) %>%
            tally() %>%
            filter(n >= maf.pop.num.threshold) %>%
            select(MARKERS) %>%
            left_join(vcf, by = "MARKERS") %>%
            arrange(MARKERS, POP_ID)
        }
      }
      message(stri_join("The number of MARKERS removed by the MAF filters = ", 
                     n_distinct(vcf$MARKERS)-n_distinct(vcf.maf$MARKERS), "\n", 
                     "The number of MARKERS before -> after the MAF filters: ", 
                     n_distinct(vcf$MARKERS)," -> ", n_distinct(vcf.maf$MARKERS), 
                     " MARKERS"))
      
      vcf <- vcf.maf %>% select(-c(MAF_LOCAL, MAF_GLOBAL))
      vcf.maf <- NULL # remove unused object
    } # end of MAF filters
    
    
    # Change the genotype coding  **********************************************
    # easier for the integration in downstream conversion to gsi_sim
    message("Recoding genotypes for gsi_sim")
    vcf <- vcf %>%
      mutate(
        REF= stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE), # replace nucleotide with numbers
        ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE),# replace nucleotide with numbers
        GT = ifelse(GT == "0/0", stri_join(REF, REF, sep = "_"),
                    ifelse(GT == "1/1",  stri_join(ALT, ALT, sep = "_"),
                           ifelse(GT == "0/1", stri_join(REF, ALT, sep = "_"),
                                  ifelse(GT == "1/0", stri_join(ALT, REF, sep = "_"), "0_0")
                           )
                    )
        )
      ) %>%
      arrange(MARKERS, POP_ID) %>%
      select(-c(REF, ALT))
    
    # more prep for the no imputation section
    gsim.prep <- vcf %>%
      tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
      tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID))
    
    # save.image("assignment.lobster.RData")  # test
    # save.image("assignment.lobster.subsample.RData")  # test
    # load("assignment.lobster.RData")  # test
    # load("assignment.lobster.subsample.RData")  # test
    
    # Imputations **************************************************************
    if (imputations != "FALSE"){
      
      vcf.prep <- vcf %>%
        mutate(
          GT = stri_replace_all_fixed(GT, pattern = "0_0", replacement = "NA", vectorize_all = FALSE),
          GT = replace(GT, which(GT == "NA"), NA)
        ) %>%
        group_by(INDIVIDUALS, POP_ID) %>% 
        tidyr::spread(data = ., key = MARKERS, value = GT) %>%
        ungroup() %>% 
        arrange(POP_ID, INDIVIDUALS)
      
      if (imputations == "rf") {
        
        # Parallel computations options
        options(rf.cores = parallel.core, mc.cores = parallel.core)
        
        # Start cluster registration backend
        cl <- parallel::makeCluster(parallel.core, methods = FALSE,outfile = "")
        
        # doSNOW::registerDoSNOW(cl)
        doParallel::registerDoParallel(cl)
        
        # imputations using Random Forest with the package randomForestSRC
        impute_markers_rf <- function(x){
          randomForestSRC::impute.rfsrc(data = x,
                                        ntree = num.tree,
                                        nodesize = 1,
                                        nsplit = split.number,
                                        nimpute = iteration.rf,
                                        do.trace = verbose)
        }
        
        # imputations by populations (default) or globally
        # default by pop
        if (imputations.group == "populations"){
          message("Imputations computed by populations, take a break...")
          df.split.pop <- split(x = vcf.prep, f = vcf.prep$POP_ID) # slip data frame by population
          pop.list <- names(df.split.pop) # list the pop
          imputed.dataset <-list() # create empty list
          # for (i in pop.list) {
          imputed.dataset <- foreach(i=pop.list, 
                                     .packages = c("plyr", "dplyr", "tidyr", 
                                                   "stringi", "readr", 
                                                   "randomForestSRC")
          ) %dopar% {
            sep.pop <- df.split.pop[[i]]
            sep.pop <- suppressWarnings(
              plyr::colwise(factor, exclude = NA)(sep.pop)
            )
            # message of progress for imputations by population
            message(paste("Completed imputations for pop ", i, sep = ""))
            imputed.dataset[[i]] <- impute_markers_rf(sep.pop)
          }
          # close parallel connection settings
          stopCluster(cl)
          message("Almost finished with the imputations...")
          vcf.imp <- suppressWarnings(as.data.frame(bind_rows(imputed.dataset)))
          
          # Second round of imputations: remove introduced NA if some pop don't have the markers by using
          # RF globally
          vcf.imp <- suppressWarnings(plyr::colwise(factor, exclude = NA)(vcf.imp)) # Make the columns factor
          vcf.imp <- impute_markers_rf(vcf.imp) # impute globally
          
          # dump unused objects
          df.split.pop <- NULL
          pop.list <- NULL
          sep.pop <- NULL
          imputed.dataset <- NULL
          vcf.prep <- NULL
          
        } else if (imputations.group == "global"){
          # Globally (not by pop_id)
          message("Imputations computed globally, take a break...")
          vcf.prep <- plyr::colwise(factor, exclude = NA)(vcf.prep)
          vcf.imp <- impute_markers_rf(vcf.prep)
          
          vcf.prep <- NULL # remove unused object
          
        }
        
      } else if (imputations == "max") {
        
        if (imputations.group == "populations"){
          message("Imputations computed by populations")
          
          vcf.imp <- suppressWarnings(
            vcf.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
              group_by(MARKERS, POP_ID) %>%
              mutate(
                GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                GT = replace(GT, which(GT == "NA"), NA)
              ) %>%
              # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
              # will take the global observed values by markers for those cases.
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              group_by(INDIVIDUALS, POP_ID) %>% 
              tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
          
          vcf.prep <- NULL # remove unused object
          
        } else if (imputations.group == "global"){
          # Globally (not by pop_id)
          message("Imputations computed globally")
          
          vcf.imp <- suppressWarnings(
            vcf.prep %>%
              tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
              group_by(MARKERS) %>%
              mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
              group_by(INDIVIDUALS, POP_ID) %>% 
              tidyr::spread(data = ., key = MARKERS, value = GT) %>%
              ungroup()
          )
          
          vcf.prep <- NULL # remove unused object
        }
      }
      
      # transform the imputed dataset into gsi_sim
      message("Imputed VCF into factory for conversion into gsi_sim...")
      gsi.prep.imp <- suppressWarnings(
        vcf.imp %>%
          tidyr::gather(key = MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>% # make tidy
          tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
          tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID)) # make tidy
      )
    } # End imputations
    
    # Sampling of markers ******************************************************
    # unique list of markers after all the filtering
    # if "all" is present in the list, change to the maximum number of markers
    unique.markers <- gsim.prep %>% 
      select(MARKERS) %>% 
      distinct(MARKERS) %>% 
      arrange(MARKERS)
    
    marker.number <- stri_replace_all_fixed(str = marker.number, pattern = "all", 
                                            replacement = nrow(unique.markers), 
                                            vectorize_all = TRUE)
    
    # Functions ******************************************************************
    # Fst function: Weir & Cockerham 1984
    fst_WC84 <- function(data, holdout.samples, ...){
      # data <- vcf # test
      # holdout.samples <- holdout$INDIVIDUALS # test
      
      pop.number <- n_distinct(data$POP_ID)
      
      if (is.null(holdout.samples)){ # use all the individuals
        data.genotyped <- data %>%
          filter(GT != "0_0")
      } else{ # with holdout set
        data.genotyped <- data %>%
          filter(GT != "0_0") %>% # remove missing genotypes
          filter(!INDIVIDUALS %in% holdout.samples) # remove supplementary individual before ranking markers with Fst
      }
      
      n.pop.locus <- data.genotyped %>%
        select(MARKERS, POP_ID) %>%
        group_by(MARKERS) %>%
        distinct(POP_ID) %>%
        tally %>%
        rename(npl = n)
      
      ind.count.locus <- data.genotyped %>%
        group_by(MARKERS) %>%
        tally
      
      ind.count.locus.pop <- data.genotyped %>%
        group_by(POP_ID, MARKERS) %>%
        tally %>%
        rename(nal = n) %>%
        mutate(
          nal_sq = nal^2,
          N_IND_GENE = nal*2
        )
      
      #common
      freq.al.locus <- data.genotyped %>%
        mutate(
          A1 = stri_sub(GT, 1, 1),
          A2 = stri_sub(GT, 3, 3)
        ) %>%
        select(-GT) %>%
        tidyr::gather(key = ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, POP_ID, MARKERS))
      
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
          select(MARKERS, ALLELES) %>%
          group_by(MARKERS, ALLELES) %>%
          distinct(MARKERS, ALLELES) %>%
          full_join(ind.count.locus, by = "MARKERS") %>%
          rename(ntal = n) %>%
          full_join(mean.n.pop.corrected.per.locus, by = "MARKERS")
      )
      
      fst.ranked <- suppressWarnings(
        data.genotyped %>%
          mutate(het = ifelse(stri_sub(GT, 1, 1) != stri_sub(GT, 3, 3), 1, 0)) %>%
          group_by(MARKERS, POP_ID) %>%
          summarise(mho = sum(het, na.rm = TRUE)) %>%  # = the number of heterozygote individuals per pop and markers
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
          ) %>%
          group_by(MARKERS) %>%
          summarise(
            lsiga = sum(siga, na.rm = TRUE),
            lsigb = sum(sigb, na.rm = TRUE),
            lsigw = sum(sigw, na.rm = TRUE),
            FST = round(lsiga/(lsiga + lsigb + lsigw), 6),
            FIS = round(lsigb/(lsigb + lsigw), 6)
          ) %>%
          arrange(desc(FST)) %>%
          select(MARKERS, FST) %>%
          mutate(RANKING = seq(from = 1, to = n()))
      )
      
      # select(MARKERS, FIS, FST)
      return(fst.ranked)
    } # end Fst function

    # Write the files
    write_gsi <- function (data, markers.names, imputations, filename, i, m, ...){
      
      data$POP_ID <- droplevels(x = data$POP_ID)
      n.individuals <- n_distinct(data$INDIVIDUALS)  # number of individuals
      pop <- data$POP_ID  # Create a vector with the population ordered by levels
      data <- suppressWarnings(data %>% select(-POP_ID))  # remove pop id
      gsi_sim.split <- split(data, pop)  # split gsi_sim by populations
      filename <- filename  # gsi_sim filename

      # filename modification based with or without imputations
      if (imputations == FALSE) {
        # message("Output...No imputation")
        filename <- stri_replace_all_fixed(filename,
                                           pattern = "txt",
                                           replacement = stri_join(
                                             i, m, 
                                             "no.imputation", "txt", sep = "."
                                           )
        )
      } else {
        # message("Output...With imputations")
        filename <- stri_replace_all_fixed(filename,
                                           pattern = "txt",
                                           replacement = stri_join(
                                             i, m, 
                                             "imputed", "txt", sep = "."
                                           )
        )
      }
      
      # directory <- getwd() # test
      # Line 1: number of individuals and the number of markers
      line1_gsi_sim <- as.data.frame(stri_join(n.individuals, m, sep = " "))
      write.table(line1_gsi_sim, file = paste0(directory.subsample, filename), col.names = FALSE, row.names = FALSE, quote = FALSE)
      
      # Markers names
      loci.table <- as.data.frame(markers.names)
      write_delim(x = loci.table, path = paste0(directory.subsample, filename), delim = "\n", append = TRUE, col_names = FALSE)
      
      # remaining lines, individuals and genotypes
      for (k in levels(pop)) {
        pop.line <- as.data.frame(stri_join("pop", k, sep = " "))
        write_delim(x = pop.line, path = paste0(directory.subsample, filename), delim = "\n", append = TRUE, col_names = FALSE)
        write_delim(x = gsi_sim.split[[k]], path = paste0(directory.subsample, filename), delim = " ", append = TRUE, col_names = FALSE)
      }
      # message(stri_join("Data file (no imputation):", filename, "\nWritten to the working directory:", directory, sep = " "))
      return(filename)
    } # end write gsi function
    
    # Assignment
    assignment_analysis <- function(data, select.markers, markers.names, missing.data, i, m, holdout, ...){
      # data <- gsim.prep #test
      # missing.data <- "no.imputation" #test
      data.select <- suppressWarnings(
        data %>%
          semi_join(select.markers, by = "MARKERS") %>%
          arrange(MARKERS) %>%  # make tidy
          tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELES, sep = "_") %>%
          arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>%
          # group_by(INDIVIDUALS, POP_ID) %>% 
          tidyr::spread(data = ., key = MARKERS_ALLELES, value = GT) %>%
          # ungroup() %>% 
          arrange(POP_ID, INDIVIDUALS)
      )
      
      if (is.null(mixture)) {
        # message("No baseline or mixture data")
        input <- write_gsi(data = data.select, markers.names = markers.names, imputations = imputations, filename = gsi_sim.filename, i = i, m = m)
      } else {
        # Baseline
        baseline.data <- suppressWarnings(
          data.select %>%
            filter(POP_ID %in% baseline) %>%
            arrange(POP_ID) %>%
            mutate(POP_ID = droplevels(POP_ID))
        )
        
        # gsi_sim baseline filename
        baseline.filename <- gsi_sim.filename
        marker.number.in.filename <- stri_join("baseline", i, m, "txt", sep = ".")
        baseline.filename <- stri_replace_all_fixed(baseline.filename, pattern = "txt",
                                                    replacement = marker.number.in.filename)
        
        # save file
        baseline.input <- write_gsi(data = baseline.data, markers.names = markers.names, imputations = imputations, filename = baseline.filename, i = i, m = m)
        # message(stri_join("Baseline data file:", filename, "\nWritten to the working directory:", directory, sep = " "))
        
        # Mixture
        mixture.data <- suppressWarnings(
          data.select %>%
            filter(POP_ID %in% mixture) %>%
            arrange(POP_ID) %>%
            mutate(POP_ID = droplevels(POP_ID))
        )
        
        # gsi_sim mixture filename
        mixture.filename <- gsi_sim.filename
        marker.number.in.filename <- stri_join("baseline", i, m, "txt", sep = ".")
        mixture.filename <- stri_replace_all_fixed(mixture.filename, pattern = "txt",
                                                   replacement = marker.number.in.filename)
        
        # save file
        mixture.input <- write_gsi(data = mixture.data, markers.names = markers.names, imputations = imputations, filename = mixture.filename, i = i, m = m)
        # message(stri_join("Mixture data file:", filename, "\nWritten to the working directory:", directory, sep = " "))
      } # end writing gsi files to disk
      
      # Run gsi_sim ------------------------------------------------------------
      if (is.null(mixture)) {
        input.gsi <- stri_join(directory.subsample,input)
        output.gsi <- stri_replace_all_fixed(input.gsi, pattern = "txt", replacement = "output.txt")
        setwd(directory.subsample)
        system(paste("gsi_sim -b", input.gsi, "--self-assign > ", output.gsi))
        } else{
        message("this option is under construction :)")
      }
      # Option remove the input file from directory to save space
      if (keep.gsi.files == FALSE){
        file.remove(input.gsi)
      }
      
      # Get Assignment results -------------------------------------------------
      # Keep track of the holdout individual
      if(sampling.method == "ranked"){
        if (THL == "all") {
          holdout.id <- NULL
        } else {
          holdout.id <- holdout$INDIVIDUALS
        }
      }
      
      # Number of markers
      n.locus <- m
      
      assignment <- suppressWarnings(
        read_delim(output.gsi, col_names = "ID", delim = "\t") %>%
          tidyr::separate(ID, c("KEEPER", "ASSIGN"), sep = ":/", extra = "warn") %>%
          filter(KEEPER == "SELF_ASSIGN_A_LA_GC_CSV") %>%
          tidyr::separate(ASSIGN, c("INDIVIDUALS", "ASSIGN"), sep = ";", extra = "merge") %>%
          tidyr::separate(ASSIGN, c("INFERRED", "OTHERS"), sep = ";", convert = TRUE, numerals = "no.loss", extra = "merge") %>%
          tidyr::separate(OTHERS, c("SCORE", "OTHERS"), sep = ";;", convert = TRUE, numerals = "no.loss", extra = "merge") %>%
          tidyr::separate(OTHERS, c("SECOND_BEST_POP", "OTHERS"), sep = ";", convert = TRUE, numerals = "no.loss", extra = "merge") %>%
          tidyr::separate(OTHERS, c("SECOND_BEST_SCORE", "OTHERS"), sep = ";;", convert = TRUE, numerals = "no.loss") %>%
          mutate(
            CURRENT = factor(stri_sub(INDIVIDUALS, pop.id.start, pop.id.end), levels = pop.levels, labels = pop.labels, ordered = T),
            CURRENT = droplevels(CURRENT),
            INFERRED = factor(INFERRED, levels = unique(pop.labels), ordered = T),
            INFERRED = droplevels(INFERRED),
            SECOND_BEST_POP = factor(SECOND_BEST_POP, levels = unique(pop.labels), ordered = T),
            SECOND_BEST_POP = droplevels(SECOND_BEST_POP),
            SCORE = round(SCORE, 2),
            SECOND_BEST_SCORE = round(SECOND_BEST_SCORE, 2),
            MARKER_NUMBER = as.numeric(rep(n.locus, n())),
            MISSING_DATA = rep(missing.data, n())
          ) %>%
          select(INDIVIDUALS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, MISSING_DATA) %>%
          arrange(CURRENT) %>%
          select(INDIVIDUALS, CURRENT, INFERRED, SCORE, MARKER_NUMBER, MISSING_DATA)
      )
      
      if(sampling.method == "ranked"){
        if (THL == "all") {
          assignment <- assignment
        } else{
          assignment <- filter(.data = assignment, INDIVIDUALS %in% holdout.id)
        }
      }
      
      if (keep.gsi.files == FALSE){
        file.remove(output.gsi)
      }
      
      # saving preliminary results
      if(sampling.method == "random"){
        assignment <- mutate(.data = assignment, ITERATIONS = rep(i, n()))
        #       write_tsv(x = assignment, path = paste0(directory,filename.ass.res), col_names = FALSE, append = TRUE) #create an empty file
      }
      
      if(sampling.method == "ranked"){
        # if (THL == 1 | THL == "all"){
        #         write_tsv(x = assignment, path = paste0(directory,filename.ass.res), col_names = FALSE, append = TRUE) #create an empty file
        #       }
        if (THL != 1 & THL != "all"){
          #         assignment <- assignment %>%
          #           group_by(CURRENT, INFERRED, MARKER_NUMBER, MISSING_DATA) %>%
          #           tally %>%
          #           group_by(CURRENT, MARKER_NUMBER) %>%
          #           mutate(TOTAL = sum(n)) %>%
          #           ungroup() %>%
          #           mutate(ASSIGNMENT_PERC = round(n/TOTAL*100, 0)) %>%
          #           select(CURRENT, INFERRED, MARKER_NUMBER, MISSING_DATA, ASSIGNMENT_PERC) %>%
          #           group_by(CURRENT, MARKER_NUMBER, MISSING_DATA) %>%
          #           tidyr::spread(data = ., key = INFERRED, value = ASSIGNMENT_PERC) %>%
          #           tidyr::gather(INFERRED, ASSIGNMENT_PERC, -c(CURRENT, MARKER_NUMBER, MISSING_DATA)) %>%
          #           mutate(ASSIGNMENT_PERC = as.numeric(stri_replace_na(ASSIGNMENT_PERC, replacement = 0))) %>%
          #           filter(as.character(CURRENT) == as.character(INFERRED)) %>%
          #           select(CURRENT, INFERRED, ASSIGNMENT_PERC, MARKER_NUMBER, MISSING_DATA) %>%
          #           mutate(ITERATIONS = rep(i, n()))
          # write_tsv(x = assignment, path = paste0(directory,filename.ass.res), col_names = FALSE, append = TRUE) #create an empty file
          assignment <- assignment %>%
            mutate(
              CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
              CURRENT = droplevels(CURRENT)
            ) %>% 
            group_by(CURRENT, MARKER_NUMBER, MISSING_DATA) %>%
            summarise(
              n = length(CURRENT[as.character(CURRENT) == as.character(INFERRED)]),
              TOTAL = length(CURRENT)
            ) %>%
            ungroup() %>% 
            mutate(
              ASSIGNMENT_PERC = round(n/TOTAL*100, 0),
              ITERATIONS = rep(i, n())
            ) %>% 
            select(-n, -TOTAL)
        }
      }
      return(assignment)
    } # end assignment analysis function
   
    # Random method ************************************************************
    if (sampling.method == "random"){
      message("Conducting Assignment analysis with markers selected randomly")
      # Number of times to repeat the sampling of markers
      iterations.list <- 1:iterations
      # iterations.list <- 1:200 # test
      
      # Plan A: use a list containing all the lists of marker combinations
      # Plan B: use nesting foreach loop (%:%)
      # Plan A is faster
      
      # Function: Random selection of marker function + iterations
      marker_selection <- function(iterations){
        m <- as.numeric(m)
        select.markers <- sample_n(tbl = unique.markers, size = m, replace = FALSE) %>%
          arrange(MARKERS) %>%
          mutate(
            ITERATIONS = rep(iterations, n()),
            MARKER_NUMBER = rep(m, n())
          )
      }
      markers.random.lists <- list()
      
      message("Making a list containing all the markers combinations")
      # Go through the function with the marker number selected
      for (m in marker.number){
        res <- purrr::map(.x = iterations.list, .f = marker_selection)
        markers.random.lists[[m]] <- res
      }
      markers.random.lists <- purrr::flatten(markers.random.lists)
      # test <- markers.random.selection.list[[101]]
      
      markers.random.lists.table <- as_data_frame(bind_rows(markers.random.lists))
      write_tsv(x = markers.random.lists.table, path = paste0(directory.subsample, "markers.random.tsv"), col_names = TRUE, append = FALSE)
      
      # Start cluster registration backend
      # parallel.core <- 8 # test
      cl <- parallel::makeCluster(parallel.core, methods = FALSE, outfile = "")
      # doSNOW::registerDoSNOW(cl)
      doParallel::registerDoParallel(cl)
      
      # Set seed for random sampling
      random.seed <- sample(x = 1:1000000, size = 1)
      # set.seed(random.seed)
      parallel::clusterSetRNGStream(cl = cl, iseed = random.seed)
      random.seed <- data.frame(RANDOM_SEED_NUMBER = random.seed)
      write_tsv(x = random.seed, path = paste0(directory.subsample, "random_seed_GBS_assignment.tsv"), col_names = TRUE, append = FALSE)
      
      mrl <- NULL
      res <- list()
      message("Starting parallel computations for the assignment analysis
First sign of progress may take some time
Progress can be monitored with activity in the folder...")
      
      # Progress Bar during parallel computations
      #     progress.max <- length(markers.random.lists)
      #     pb <- txtProgressBar(max = progress.max, title = "Assignment in progress", style = 3, width = 85)
      #     progress <- function(n) setTxtProgressBar(pb, n)
      #     opts <- list(progress = progress)
      
      # foreach
      #     assignment.res <- foreach(mrl=markers.random.lists, 
      #                               .options.snow=opts, 
      #                               .packages = c("dplyr", "tidyr", "stringi",
      #                                             "readr", "purrr")
      #     ) %dopar% {
      
      assignment.res <- foreach(
        mrl=markers.random.lists,
        .packages = c("dplyr", "tidyr", "stringi", "readr", "purrr"),
        .verbose = FALSE
      ) %dopar% {
        # mrl <- markers.random.lists[1] # test
        # mrl <- as_data_frame(purrr::flatten(mrl)) # test
        Sys.sleep(0.01)                             # for progress bar
        mrl <- data.frame(mrl)                      # marker random list
        i <- as.numeric(unique(mrl$ITERATIONS))     # iteration
        m <- as.numeric(unique(mrl$MARKER_NUMBER))  # number of marker selected
        select.markers <- mrl %>%                   # markers
          ungroup() %>% 
          select(MARKERS) %>% 
          arrange(MARKERS)
        
        # get the list of loci after filter
        markers.names <- unique(select.markers$MARKERS)
        assignment.no.imp <- assignment_analysis(data = gsim.prep, 
                                                 missing.data = "no.imputation", 
                                                 i  = i
        )
        
        # With imputations
        if (imputations != FALSE) {# with imputations
          if (imputations == "rf"){
            if (imputations.group == "populations"){
              missing.data <- "imputed RF populations"
            } else{
              missing.data <- "imputed RF global"
            }
          } else {
            if (imputations.group == "populations"){
              missing.data <- "imputed max populations"
            } else{
              missing.data <- "imputed max global"
            }
          }
          assignment.imp <- assignment_analysis(data = gsi.prep.imp, 
                                                missing.data = missing.data, 
                                                i = i
          )
        }
        
        #compile assignment results each marker number for the iteration
        if (imputations == FALSE) {
          assignment <- assignment.no.imp
        } else{
          assignment <- bind_rows(assignment.no.imp, assignment.imp)
        }
        assignment <- mutate(.data = assignment, ITERATIONS = rep(i, n()))
        return(assignment)
      } # End of iterations for both with and without imputations
      message("Summarizing the assignment analysis results")
      stopCluster(cl) # close parallel connection settings
      
      # Compiling the results
      assignment.res <- suppressWarnings(
        as_data_frame(bind_rows(assignment.res))%>%
          mutate(METHOD = rep("LOO", n()))
      )
      
      assignment.stats.pop <- suppressWarnings(
        assignment.res %>%
          group_by(CURRENT, INFERRED, MARKER_NUMBER, MISSING_DATA, ITERATIONS, METHOD) %>%
          tally %>%
          group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, ITERATIONS, METHOD) %>%
          mutate(TOTAL = sum(n)) %>%
          ungroup() %>%
          mutate(MEAN_i = round(n/TOTAL*100, 0)) %>%
          filter(as.character(CURRENT) == as.character(INFERRED)) %>%
          select(CURRENT, MEAN_i, MARKER_NUMBER, MISSING_DATA, ITERATIONS, METHOD) %>%
          mutate(
            CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = T),
            CURRENT = droplevels(CURRENT)
          ) %>%
          group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
          summarise(
            MEAN = round(mean(MEAN_i), 2),
            SE = round(sqrt(var(MEAN_i)/length(MEAN_i)), 2),
            MIN = round(min(MEAN_i), 2),
            MAX = round(max(MEAN_i), 2),
            MEDIAN = round(median(MEAN_i), 2),
            QUANTILE25 = round(quantile(MEAN_i, 0.25), 2),
            QUANTILE75 = round(quantile(MEAN_i, 0.75), 2)
          ) %>%
          arrange(CURRENT, MARKER_NUMBER)
      )
      
      pop.levels.assignment.stats.overall <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
      
      assignment.stats.overall <- assignment.stats.pop %>%
        group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
        rename(ASSIGNMENT_PERC = MEAN) %>%
        summarise(
          MEAN = round(mean(ASSIGNMENT_PERC), 2),
          SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
          MIN = round(min(ASSIGNMENT_PERC), 2),
          MAX = round(max(ASSIGNMENT_PERC), 2),
          MEDIAN = round(median(ASSIGNMENT_PERC), 2),
          QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
          QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2)
        ) %>%
        mutate(CURRENT = rep("OVERALL", n())) %>%
        arrange(CURRENT, MARKER_NUMBER)
      
      assignment.summary.stats <- suppressWarnings(
        bind_rows(assignment.stats.pop, assignment.stats.overall) %>%
          mutate(CURRENT = factor(CURRENT, levels = pop.levels.assignment.stats.overall, ordered = TRUE)) %>%
          arrange(CURRENT, MARKER_NUMBER) %>%
          mutate(
            SE_MIN = MEAN - SE,
            SE_MAX = MEAN + SE,
            ITERATIONS = rep(iterations, n())
          ) %>%
          select(CURRENT, MARKER_NUMBER, MEAN, MEDIAN, SE, MIN, MAX, QUANTILE25, QUANTILE75, SE_MIN, SE_MAX, METHOD, MISSING_DATA, ITERATIONS)
      )
      # Write the tables to directory
      # assignment results
      if (imputations == FALSE) {
        filename.assignment.res <- stri_join("assignment.res", "no.imputation", sampling.method, "tsv", sep = ".")
      } else{ # with imputations
        filename.assignment.res <- stri_join("assignment.res", "imputed", sampling.method, "tsv", sep = ".")
      }
      write_tsv(x = assignment.res, path = paste0(directory.subsample,filename.assignment.res), col_names = TRUE, append = FALSE)
      
      # assignment summary stats
      if (imputations == FALSE) {
        filename.assignment.sum <- stri_join("assignment.summary.stats", "no.imputation", sampling.method, "tsv", sep = ".")
      } else{ # with imputations
        filename.assignment.sum <- stri_join("assignment.summary.stats", "imputed", sampling.method, "tsv", sep = ".")
      }
      write_tsv(x = assignment.summary.stats, path = paste0(directory.subsample,filename.assignment.sum), col_names = TRUE, append = FALSE)
    } # end method random
    
    # Ranked method ************************************************************
    if (sampling.method == "ranked"){
      message("Conducting Assignment analysis with ranked markers")
      
      # List of all individuals
      ind.pop.df<- vcf %>% 
        select(POP_ID, INDIVIDUALS) %>% 
        distinct(POP_ID, INDIVIDUALS)
      
      # THL selection
      message("Using THL method, ranking Fst with training samples...")
      if (THL == 1){
        # Will go through the individuals in the list one by one.
        iterations.list <- ind.pop.df$INDIVIDUALS
        # Keep track of holdout individuals
        holdout.individuals <- ind.pop.df %>%
          mutate(ITERATIONS = stri_join("HOLDOUT", seq(1:n()), sep = "_"))
      } else if (THL == "all") { # no holdout for that one
        iterations.list <- iterations
        holdout.individuals <- NULL
        message("Warning: using all the individuals for ranking markers based on Fst\nNo holdout samples")
        message("Recommended reading: \nAnderson, E. C. (2010) Assessing the power of informative subsets of
                loci for population assignment: standard methods are upwardly biased.\nMolecular ecology resources 10, 4:701-710.")
      } else {
        # Create x (iterations) list of y (THL) proportion of individuals per pop.
        if (stri_detect_fixed(THL, ".") & THL < 1) {
          # iterations <- 5 # test
          # THL <- 0.4 # test
          holdout.individuals.list <- list()
          iterations.list <- 1:iterations
          for (x in 1:iterations){
            holdout.individuals <- ind.pop.df %>%
              group_by(POP_ID) %>%
              sample_frac(THL, replace = FALSE) %>%  # sampling fraction for each pop
              arrange(POP_ID, INDIVIDUALS) %>%
              ungroup() %>%
              select(INDIVIDUALS) %>%
              mutate(ITERATIONS = rep(x, n()))
            holdout.individuals.list[[x]] <- holdout.individuals
          }
          holdout.individuals <- as.data.frame(bind_rows(holdout.individuals.list))
        }
        
        # Create x (iterations) list of y (THL) individuals per pop.
        if (THL > 1) {
          holdout.individuals.list <- list()
          iterations.list <- 1:iterations
          for (x in 1:iterations){
            holdout.individuals <- ind.pop.df %>%
              group_by(POP_ID) %>%
              sample_n(THL, replace = FALSE) %>% # sampling individuals for each pop
              arrange(POP_ID, INDIVIDUALS) %>%
              ungroup() %>%
              select(INDIVIDUALS) %>%
              mutate(ITERATIONS = rep(x, n()))
            holdout.individuals.list[[x]] <- holdout.individuals
          }
          holdout.individuals <- as.data.frame(bind_rows(holdout.individuals.list))
        }
        message("Holdout samples saved in your folder")
      } # end tracking holdout individuals
      write_tsv(x = holdout.individuals, 
                path = paste0(directory.subsample,"holdout.individuals.tsv"), 
                col_names = TRUE, 
                append = FALSE
      )
      message("Holdout samples saved in your folder")
      
      # test
      #     # Preparing file for saving to directory preliminary results
      #     filename.ass.res <- "assignment.preliminiary.results.tsv"
      #     if (THL == 1 | THL == "all"){
      #       ass.res <- data_frame(INDIVIDUALS = character(0), 
      #                             CURRENT = character(0), 
      #                             INFERRED = character(0), 
      #                             SCORE = integer(0), 
      #                             MARKER_NUMBER = integer(0), 
      #                             MISSING_DATA = character(0)
      #       ) #create an empty dataframe
      #     } else{#THL != 1... > 1 number or < 1 for proportions
      #       ass.res <- data_frame(CURRENT = character(0),
      #                             INFERRED = character(0),
      #                             ASSIGNMENT_PERC = integer(0),
      #                             MARKER_NUMBER = integer(0),
      #                             ITERATIONS = integer(0),
      #                             MISSING_DATA = character(0)
      #       ) #create an empty dataframe
      #     }
      #     write_tsv(x = ass.res, 
      #               path = paste0(directory,filename.ass.res), 
      #               col_names = TRUE, 
      #               append = FALSE
      #     ) #create an empty file
      
      
      # Going through the loop of holdout individuals
      message("Starting parallel computations for the assignment analysis
First sign of progress may take some time
Progress can be monitored with activity in the folder...")
      
      # Progress Bar during parallel computations
      #     if (THL == 1){
      #       progress.max <- length(iterations.list)
      #     } else {
      #       progress.max <- iterations
      #     }
      #     pb <- txtProgressBar(max = progress.max, 
      #                          title = "Assignment in progress", 
      #                          style = 3, 
      #                          width = 85
      #     )
      #     progress <- function(n) setTxtProgressBar(pb, n)
      #     opts <- list(progress = progress)
      #     
      # Start cluster registration backend
      # cl <- parallel::makeCluster(parallel.core, methods = FALSE, outfile = "") # test
      # cl <- parallel::makeCluster(parallel.core)
      # cl <- parallel::makeCluster(parallel.core, methods = FALSE, outfile = paste0(directory, "parallel.computations.log")) # test
      # doSNOW::registerDoSNOW(cl)
      # doParallel::registerDoParallel(cl)
      
      # foreach
      # i <- NULL
#       assignment.res <- list()
      ## foreach with progress bar option using SNOW
      #     assignment.res <- foreach(
      #       i=iterations.list, .options.snow=opts, 
      #       .packages = c("dplyr", "tidyr", "stringi", "readr"),
      #       .verbose = FALSE
      #     ) %dopar% {
      
      #       assignment.res <- foreach(
      #         i=iterations.list, 
      #         .packages = c("dplyr", "tidyr", "stringi", "readr"),
      #         .verbose = FALSE
      #       ) %dopar% {
      
      assignment_ranking <- function(iterations.list, ...){
        # Sys.sleep(0.01)  # For progress bar
        # i <- "TRI_09" #test
        # i <- "CAR_01" #test
        # i <- 1
        # i <- 5
        i <- iterations.list
        
        # Ranking Fst with training dataset (keep holdout individuals out)
        message("Ranking markers based on Fst")
        if (THL == "all"){
          holdout <- NULL
          fst.ranked <- fst_WC84(data = vcf, holdout.samples = NULL)
          if (imputations != FALSE){
            fst.ranked.imp <- fst_WC84(data = vcf.imp, holdout.samples = NULL)
          }
        } else if (THL == 1) {
          holdout <- data.frame(INDIVIDUALS = i)
          fst.ranked <- fst_WC84(data = vcf, holdout.samples = holdout$INDIVIDUALS)
          if (imputations != FALSE){
            fst.ranked.imp <- fst_WC84(data = vcf.imp, holdout.samples = holdout$INDIVIDUALS)
          }
        } else {
          holdout <- data.frame(holdout.individuals.list[i])
          fst.ranked <- fst_WC84(data = vcf, holdout.samples = holdout$INDIVIDUALS)
          if (imputations != FALSE){
            fst.ranked.imp <- fst_WC84(data = vcf.imp, holdout.samples = holdout$INDIVIDUALS)
          }
        }
        
        # Saving Fst
#         if (THL != 1 & THL != "all") { # for THL != 1 (numbers and proportions)
#           i <- unique(holdout$ITERATIONS)
#         }
        
        fst.ranked.filename <- stri_join("fst.ranked_", i, ".tsv", sep = "") # No imputation
        write_tsv(x = fst.ranked, path = paste0(directory.subsample, fst.ranked.filename), 
                  col_names = TRUE, 
                  append = FALSE
        )
        
        if (imputations != FALSE){  # With imputations
          fst.ranked.filename.imp <- stri_join("fst.ranked_", i,
                                                "_imputed",
                                                ".tsv", 
                                                sep = ""
          ) 
          write_tsv(x = fst.ranked.imp, 
                    path = paste0(directory.subsample, fst.ranked.filename.imp), 
                    col_names = TRUE, 
                    append = FALSE
          )
        }
        
        # Markers numbers loop
        # Create empty lists to feed the results
        message("Going throught the marker.number")
        assignment.marker <- list()
        assignment_marker_loop <- function(m, ...){
          message("Marker number: ", m)
          # for (m in marker.number) {
          # message("Marker number: ", m)
          # m <- 200 # test
          # m <- 400 # test
          # rm(m) # test
          m <- as.numeric(m)
          # No imputation
#           select.markers <- fst.ranked %>%
#             filter(row_number() <= m) %>%
#             select(MARKERS)
          RANKING <- NULL
          select.markers <- filter(.data = fst.ranked, RANKING <= m) %>%
            select(MARKERS)
          
          # get the list of markers after filter
          markers.names <- unique(select.markers$MARKERS)
          
          # Assignment analysis without imputations
          assignment.no.imp <- assignment_analysis(data = gsim.prep,
                                                   select.markers = select.markers,
                                                   markers.names = markers.names,
                                                   missing.data = "no.imputation", 
                                                   i = i, 
                                                   m = m,
                                                   holdout = holdout
          )
          
          # With imputations
          if (imputations != FALSE) {  # with imputations
#             select.markers <- fst.ranked.imp %>%  # not the same in no imputation
#               filter(row_number() <= m) %>%
#               select(MARKERS)
            select.markers <- filter(.data = fst.ranked.imp, RANKING <= m) %>%
              select(MARKERS)
            
            # get the list of markers after filter
            markers.names <- unique(select.markers$MARKERS)  # not the same in no imputation
            
            if (imputations == "rf"){
              if (imputations.group == "populations"){
                missing.data <- "imputed RF populations"
              } else{
                missing.data <- "imputed RF global"
              }
            } else {
              if (imputations.group == "populations"){
                missing.data <- "imputed max populations"
              } else{
                missing.data <- "imputed max global"
              }
            }
            
            # Assignment analysis WITH imputations
            assignment.imp <- assignment_analysis(data = gsi.prep.imp,
                                                  select.markers = select.markers,
                                                  markers.names = markers.names,
                                                  missing.data = missing.data, 
                                                  i = i,
                                                  select.markers = select.markers,
                                                  m = m,
                                                  holdout = holdout
            )
          }
          
          #compile assignment results each marker number for the iteration
          if (imputations == FALSE) {# with imputations
            assignment <- assignment.no.imp
          } else{
            assignment <- bind_rows(assignment.no.imp, assignment.imp)
          }
          m <- as.character(m)
          assignment.marker[[m]] <- assignment
          return(assignment.marker)
        }  # End marker number loop for both with and without imputations
        
        assignment.marker <- map(
          .x = marker.number, 
          .f = assignment_marker_loop,
          fst.ranked = fst.ranked,
          i = i,
          vcf = vcf,
          gsim.prep = gsim.prep,
          pop.levels = pop.levels,
          pop.labels = pop.labels,
          pop.id.start =  pop.id.start,
          pop.id.end = pop.id.end,
          sampling.method = sampling.method,
          THL = THL,
          iterations = iterations,
          gsi_sim.filename = gsi_sim.filename,
          keep.gsi.files = keep.gsi.files,
          baseline = baseline,
          mixture = mixture,
          imputations = imputations,
          parallel.core = parallel.core
        )
        
        
        message("Summarizing the assignment analysis results by iterations and marker group")
        # if (THL == 1) {
#           assignment.res.summary <- suppressWarnings(
#             as_data_frame(bind_rows(assignment.marker)) %>%
#               mutate(METHOD = rep("THL", n()))
#           )
          
          assignment.res.summary <- suppressWarnings(
            bind_rows(purrr::flatten(assignment.marker)) %>%
              mutate(METHOD = rep("THL", n()))
          )
          
          res.filename <- stri_join("assignment_ind_", i, ".tsv", sep = "") # No imputation
          write_tsv(x = assignment.res.summary, path = paste0(directory.subsample, res.filename), 
                    col_names = TRUE, 
                    append = FALSE
          )
        # }
        return(assignment.res.summary)
      }  # End assignment ranking function
      # stopCluster(cl)  # close parallel connection settings
      
      # using mclapply
      assignment.res <- list()
      assignment.res <- parallel::mclapply(
        X = iterations.list, 
        FUN = assignment_ranking, 
        mc.preschedule = FALSE, 
        mc.silent = FALSE, 
        mc.cores = parallel.core,
        marker.number = marker.number
        )
#       , 
#         vcf = vcf,
#         gsi.prep = gsi.prep,
#         unique.markers = unique.markers,
#         ind.pop.df = ind.pop.df,
#         holdout.individuals = holdout.individuals,
#         marker.number = marker.number,
#         pop.levels = pop.levels,
#         pop.labels = pop.labels,
#         pop.id.start =  pop.id.start,
#         pop.id.end = pop.id.end,
#         sampling.method = sampling.method,
#         THL = THL,
#         iterations = iterations,
#         gsi_sim.filename = gsi_sim.filename,
#         keep.gsi.files = keep.gsi.files,
#         baseline = baseline,
#         mixture = mixture,
#         imputations = imputations,
#         imputations.group = imputations.group,
#         num.tree = num.tree,
#         iteration.rf = iteration.rf,
#         split.number = split.number,
#         verbose = verbose,
#         parallel.core = parallel.core
#       )
      
      # Compiling the results
      message("Compiling results")
      assignment.res.summary <- suppressWarnings(
       bind_rows(assignment.res) %>%
          mutate(METHOD = rep("THL", n()))
      )
      
      if (THL == 1 | THL == "all"){
        assignment.stats.pop <- assignment.res.summary %>%
          mutate(
            CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
            CURRENT = droplevels(CURRENT)
          ) %>% 
          group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
          summarise(
            n = length(CURRENT[as.character(CURRENT) == as.character(INFERRED)]),
            TOTAL = length(CURRENT)
          ) %>%
          ungroup() %>% 
          mutate(MEAN = round(n/TOTAL*100, 0)) %>% 
          select(-n, -TOTAL)
        #       assignment.stats.pop <- assignment.res.summary %>%
        #         group_by(CURRENT, INFERRED, MARKER_NUMBER, MISSING_DATA, METHOD) %>% 
        #         tally %>% 
        #         group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>% 
        #         mutate(TOTAL = sum(n)) %>% 
        #         ungroup() %>% 
        #         mutate(MEAN_i = round(n/TOTAL*100, 0)) %>% 
        #         filter(as.character(CURRENT) == as.character(INFERRED)) %>% 
        #         select(CURRENT, MEAN_i, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
        #         mutate(
        #           CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = T),
        #           CURRENT = droplevels(CURRENT)
        #         ) %>% 
        #         group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>% 
        #         summarise(
        #           MEAN = round(mean(MEAN_i), 2),
        #           SE = round(sqrt(var(MEAN_i)/length(MEAN_i)), 2),
        #           MIN = round(min(MEAN_i), 2),
        #           MAX = round(max(MEAN_i), 2),
        #           MEDIAN = round(median(MEAN_i), 2),
        #           QUANTILE25 = round(quantile(MEAN_i, 0.25), 2),
        #           QUANTILE75 = round(quantile(MEAN_i, 0.75), 2)
        #         ) %>%
        #         arrange(CURRENT, MARKER_NUMBER)
        
        pop.levels.assignment.stats.overall <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
        
        assignment.stats.overall <- assignment.stats.pop %>% 
          group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
          rename(ASSIGNMENT_PERC = MEAN) %>%
          summarise(
            MEAN = round(mean(ASSIGNMENT_PERC), 2),
            SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
            MIN = round(min(ASSIGNMENT_PERC), 2),
            MAX = round(max(ASSIGNMENT_PERC), 2),
            MEDIAN = round(median(ASSIGNMENT_PERC), 2),
            QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
            QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2)
          ) %>% 
          mutate(CURRENT = rep("OVERALL", n())) %>% 
          arrange(CURRENT, MARKER_NUMBER)
        
        assignment.summary.stats <- suppressWarnings(
          bind_rows(assignment.stats.pop, assignment.stats.overall) %>%
            mutate(CURRENT = factor(CURRENT, levels = pop.levels.assignment.stats.overall, ordered = TRUE)) %>%
            arrange(CURRENT, MARKER_NUMBER) %>%
            mutate(
              SE_MIN = MEAN - SE,
              SE_MAX = MEAN + SE
            )
        )
        
      } else {
        # THL != 1 or "all"
        # summary stats
        assignment.stats.pop <- assignment.res.summary %>%
          mutate(
            CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
            CURRENT = droplevels(CURRENT)
          ) %>%
          group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
          summarise(
            MEAN = round(mean(ASSIGNMENT_PERC), 2),
            SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
            MIN = round(min(ASSIGNMENT_PERC), 2),
            MAX = round(max(ASSIGNMENT_PERC), 2),
            MEDIAN = round(median(ASSIGNMENT_PERC), 2),
            QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
            QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2),
            SE_MIN = MEAN - SE,
            SE_MAX = MEAN + SE
          ) %>%
          arrange(CURRENT, MARKER_NUMBER)
        
        pop.levels.assignment.stats.overall <- c(levels(assignment.stats.pop$CURRENT), "OVERALL")
        
        assignment.stats.overall <- assignment.stats.pop %>%
          group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
          rename(ASSIGNMENT_PERC = MEAN) %>%
          summarise(
            MEAN = round(mean(ASSIGNMENT_PERC), 2),
            SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
            MIN = round(min(ASSIGNMENT_PERC), 2),
            MAX = round(max(ASSIGNMENT_PERC), 2),
            MEDIAN = round(median(ASSIGNMENT_PERC), 2),
            QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
            QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2),
            SE_MIN = MEAN - SE,
            SE_MAX = MEAN + SE
          ) %>%
          mutate(CURRENT = rep("OVERALL", n())) %>%
          arrange(CURRENT, MARKER_NUMBER)
        
        assignment.summary.stats <- suppressWarnings(
          bind_rows(assignment.stats.pop, assignment.stats.overall) %>%
            mutate(CURRENT = factor(CURRENT, levels = pop.levels.assignment.stats.overall, ordered = TRUE)) %>%
            arrange(CURRENT, MARKER_NUMBER)
        )
      } # end THL != 1
      
      # Write the tables to directory
      # assignment results
      if (imputations == FALSE) {
        filename.assignment.res <- stri_join("assignment.res", "no.imputation", sampling.method, "tsv", sep = ".")
      } else{ # with imputations
        filename.assignment.res <- stri_join("assignment.res", "imputed", sampling.method, "tsv", sep = ".")
      }
      write_tsv(x = assignment.res.summary, path = paste0(directory.subsample,filename.assignment.res), col_names = TRUE, append = FALSE)
      
      # assignment summary stats
      if (imputations == FALSE) {
        filename.assignment.sum <- stri_join("assignment.summary.stats", "no.imputation", sampling.method, "tsv", sep = ".")
      } else{ # with imputations
        filename.assignment.sum <- stri_join("assignment.summary.stats", "imputed", sampling.method, "tsv", sep = ".")
      }
      write_tsv(x = assignment.summary.stats, path = paste0(directory.subsample,filename.assignment.sum), col_names = TRUE, append = FALSE)
    } # end of ranked THL method
    
    # update the assignment with subsampling iterations id
    assignment.summary.stats <- assignment.summary.stats %>% 
      mutate(SUBSAMPLE = rep(subsample.id, n()))
    return(assignment.summary.stats)
} # end assignment_function
  
  res <- map(.x = subsample.list, .f = assignment_function,
             vcf = vcf,
             snp.LD = snp.LD,
             common.markers = common.markers,
             maf.local.threshold = maf.local.threshold,
             maf.global.threshold = maf.global.threshold,
             maf.pop.num.threshold = maf.pop.num.threshold,
             maf.approach = maf.approach,
             maf.operator = maf.operator,
             marker.number = marker.number,
             pop.levels = pop.levels,
             pop.labels = pop.labels,
             pop.id.start =  pop.id.start,
             pop.id.end = pop.id.end,
             sampling.method = sampling.method,
             THL = THL,
             iterations = iterations,
             gsi_sim.filename = gsi_sim.filename,
             keep.gsi.files = keep.gsi.files,
             baseline = baseline,
             mixture = mixture,
             imputations = imputations,
             imputations.group = imputations.group,
             num.tree = num.tree,
             iteration.rf = iteration.rf,
             split.number = split.number,
             verbose = verbose,
             parallel.core = parallel.core
  )
  res <- bind_rows(res)
  write_tsv(x = res, path = paste0(directory, "assignment.results.tsv"), col_names = TRUE, append = FALSE)
  
  # Summary of the subsampling iterations
  if(iterations.subsample > 1){
    res.pop <- res %>%
      filter(CURRENT != "OVERALL") %>%
      group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
      rename(ASSIGNMENT_PERC = MEAN) %>%
      summarise(
        MEAN = round(mean(ASSIGNMENT_PERC), 2),
        SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
        MIN = round(min(ASSIGNMENT_PERC), 2),
        MAX = round(max(ASSIGNMENT_PERC), 2),
        MEDIAN = round(median(ASSIGNMENT_PERC), 2),
        QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
        QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2)
        ) %>% 
      mutate(SUBSAMPLE = rep("OVERALL", n())) %>%
      arrange(CURRENT, MARKER_NUMBER)
    
    res.overall <- res.pop %>% 
      group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
      rename(ASSIGNMENT_PERC = MEAN) %>%
      summarise(
        MEAN = round(mean(ASSIGNMENT_PERC), 2),
        SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
        MIN = round(min(ASSIGNMENT_PERC), 2),
        MAX = round(max(ASSIGNMENT_PERC), 2),
        MEDIAN = round(median(ASSIGNMENT_PERC), 2),
        QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
        QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2)
      ) %>%
      mutate(
        CURRENT = rep("OVERALL", n()),
        SUBSAMPLE = rep("OVERALL", n())
      ) %>% 
      arrange(CURRENT, MARKER_NUMBER)
    
    res.pop.overall <- suppressWarnings(
      bind_rows(res.pop, res.overall) %>%
        mutate(CURRENT = factor(CURRENT, levels = levels(res.pop$CURRENT), ordered = TRUE)) %>%
        arrange(CURRENT, MARKER_NUMBER) %>%
        mutate(
          SE_MIN = MEAN - SE,
          SE_MAX = MEAN + SE
        ) %>%
        select(CURRENT, MARKER_NUMBER, MEAN, MEDIAN, SE, MIN, MAX, QUANTILE25, QUANTILE75, SE_MIN, SE_MAX, METHOD, MISSING_DATA, SUBSAMPLE)
    )

    res <- bind_rows(
      res %>% 
        mutate(SUBSAMPLE = as.character(SUBSAMPLE)), 
      res.pop.overall) %>% 
      mutate(SUBSAMPLE = factor(SUBSAMPLE, levels = c(1:iterations.subsample, "OVERALL"), ordered = TRUE)) %>% 
      arrange(CURRENT, MARKER_NUMBER, SUBSAMPLE)
    
    # unused objects
    res.pop.overall <- NULL
    res.overall <- NULL
    res.pop <- NULL
  } # End summary of the subsampling iterations
  
  # Assignment plot
  plot.assignment <- ggplot(res, aes(x = factor(MARKER_NUMBER), y = MEAN))+
    geom_point(size = 2, alpha = 0.5) +
    geom_errorbar(aes(ymin = SE_MIN, ymax = SE_MAX), width = 0.3) +
    scale_y_continuous(breaks = c(0, 10, 20 ,30, 40, 50, 60, 70, 80, 90, 100))+
    labs(x = "Marker number")+
    labs(y = "Assignment success (%)")+
    theme_bw()+
    theme(
      legend.position = "bottom",      
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"), 
      axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.x = element_text(size = 8, family = "Helvetica", face = "bold", angle = 90, hjust = 0.5, vjust = 0.5), 
      axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
      axis.text.y = element_text(size = 10, family = "Helvetica", face = "bold")
    )

  # results
  res.list <- list(assignment = res, plot.assignment = plot.assignment)
  return(res.list)
} # end GBS_assignment

