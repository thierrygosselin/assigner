# Write a gsi_sim file from STACKS VCF file

#' @name assignment_ngs
#' @title Assignment analysis of next-generation sequencing data (GBS/RADseq, 
#' SNP chip, etc) using gsi_sim and adegenet. 
#' @description \code{gsi_sim} is a tool for doing and simulating genetic stock
#' identification and developed by Eric C. Anderson.
#' The arguments in the \code{assignment_ngs} function were tailored for the
#' reality of GBS/RADseq data for assignment analysis while
#' maintaining a reproducible workflow.
#' Various input files is offered. Individuals, populations and
#' markers can be filtered and/or selected in several ways using blacklist,
#' whitelist and other arguments. Map-independent imputation of missing genotype
#' using Random Forest or the most frequent category is also available.
#' Markers can be randomly selected for a classic LOO (Leave-One-Out)
#' assignment or chosen based on ranked Fst for a thl
#' (Training, Holdout, Leave-one-out) assignment analysis.

#' @param data Options include the VCF (1) or an haplotype files (2) created in STACKS 
#' (\code{data = "batch_1.vcf"} and \code{data = "batch_1.haplotypes.tsv"}, 
#' respectively) or a data frame (3) with tab separating the 
#' genotypes in columns (\code{data = "data.assignment.tsv"}). 
#' The 1st column is the \code{POP_ID}, 2nd colum 
#' the \code{INDIVIDUALS} and the remaining columns are the markers IDs
#' containing genotypes in the format: 3 digits per allele
#' and no space between alleles (e.g. 235240 : allele1 = 235 and allele2 = 240).
#' Missing genotypes are coded \code{0} or \code{000000}. 
#' Note that the \code{POP_ID} column can be any hierarchical grouping. 
#' See the argument \code{strata} for other means of controlling grouping used 
#' in the assignment. The last option for data input is a PLINK file in 
#' \code{tped/tfam} format (e.g. \code{data =  "data.assignment.tped"}). 
#' The first 2 columns of the \code{tfam} file will be used for the 
#' \code{strata} argument below, unless a new one is provided. 
#' Columns 1, 3 and 4 of the \code{tped} are discarded. The remaining columns 
#' correspond to the genotype in the format \code{01/04} 
#' where \code{A = 01, C = 02, G = 03 and T = 04}. For \code{A/T} format, use 
#' PLINK or bash to convert.
#' Use VCFTOOLS \url{http://vcftools.sourceforge.net/} with \code{--plink-tped} 
#' to convert very large VCF file. For \code{.ped} file conversion to 
#' \code{.tped} use PLINK \url{http://pngu.mgh.harvard.edu/~purcell/plink/} 
#' with \code{--recode transpose}.
#' 
#' @param assignment.analysis Assignment analysis conducted with 
#' \code{assignment.analysis = "gsi_sim"} or 
#' \code{assignment.analysis = "adegenet"}.
#' 
#' @param whitelist.markers (optional) A whitelist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by chromosome and/or locus and/or by snp.
#' The whitelist is in the working directory (e.g. "whitelist.txt").
#' de novo CHROM column with 'un' need to be changed to 1. 
#' Default \code{NULL} for no whitelist of markers. In the VCF, the column ID is
#' the LOCUS identification.

#' @param monomorphic.out (optional) For PLINK file, should the monomorphic 
#' markers present in the dataset be filtered out ? 
#' Default: \code{monomorphic.out = TRUE}.

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

#' @param snp.ld (optional) For VCF file only. With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 3 options: \code{"random"} selection, \code{"first"} or
#' \code{"last"} SNP on the same short read/haplotype. 
#' Default: \code{snp.ld = NULL}.
#' Note that for other file type, use stackr package for haplotype file and 
#' create a whitelist, for plink and data frames, use PLINK linkage 
#' disequilibrium based SNP pruning option.
#' @param common.markers (optional) Logical. Default: \code{common.markers = TRUE}, 
#' will only keep markers in common (genotyped) between all the populations.


#' @param maf.thresholds (string, double, optional) String with 
#' local/populations and global/overall maf thresholds, respectively.
#' Default: \code{maf.thresholds = NULL}. 
#' e.g. \code{maf.thresholds = c(0.05, 0.1)} for a local maf threshold 
#' of 0.05 and a global threshold of 0.1. Available for VCF, PLINK and data frame 
#' files. Use stackr for haplotypes files.
#' @param maf.pop.num.threshold (integer, optional) When maf thresholds are used,
#' this argument is for the number of pop required to pass the maf thresholds
#' to keep the locus. Default: \code{maf.pop.num.threshold = 1}
#' @param maf.approach (character, optional). By \code{maf.approach = "SNP"} or 
#' by \code{maf.approach = "haplotype"}.
#' The function will consider the SNP or ID/LOCUS/haplotype/read MAF statistics 
#' to filter the markers.
#' Default is \code{maf.approach = "SNP"}. The \code{haplotype} approach is 
#' restricted to VCF file.
#' @param maf.operator (character, optional) \code{maf.operator = "AND"} or default \code{maf.operator = "OR"}.
#' When filtering over LOCUS or SNP, do you want the local \code{"AND"}
#' global MAF to pass the thresholds, or ... you want the local \code{"OR"}
#' global MAF to pass the thresholds, to keep the marker?

#' @param max.marker An optional integer useful to subsample marker number in 
#' large PLINK file. Default: \code{max.marker = NULL}. e.g. if the data set 
#' contains 200 000 markers and \code{max.marker = 10000} 10000 markers are
#' subsampled randomly from the 200000 markers. Use \code{whitelist.markers} to
#' keep specific markers.
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
#' Training-Holdout-Leave One Out (thl) assignment ?
#' @param thl (character, integer, proportion) For \code{sampling.method = "ranked"} only.
#' Default \code{1}, 1 individual sample is used as holdout. This individual is not
#' participating in the markers ranking. For each marker number,
#' the analysis will be repeated with all the indiviuals in the data set
#' (e.g. 500 individuals, 500 times 500, 1000, 2000 markers).
#' If a proportion is used e.g. \code{thl = 0.15}, 15 percent of individuals in each
#' populations are chosen randomly as holdout individuals.
#' With \code{thl = "all"} all individuals are used for ranking (not good) and
#' \code{iteration.method} argument below is set to \code{1} by default.
#' For the other thl values, you can create different holdout individuals lists
#' with the \code{iteration.method} argument below (bootstrap).
#' @param iteration.method With random marker selection the iterations argument =
#' the number of iterations to repeat marker resampling. 
#' Default: \code{iteration.method = 10}.
#' With \code{marker.number = c(500, 1000)} and default iterations setting,
#' 500 markers will be randomly chosen 10 times and 1000 markers will be randomly
#' chosen 10 times. For the ranked method, using \code{thl = 1}, the analysis
#' will be repeated for each individuals in the data set for every
#' \code{marker.number} selected. With a proportion argument \code{thl = 0.15},
#' 15 percent of individuals in each populations are chosen randomly as holdout
#' individuals and this process is reapeated the number of times chosen by the
#' \code{iteration.method} value.


#' @param folder (optional) The name of the folder created in the working directory to save the files/results.
#' @param gsi_sim.filename (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. Default \code{assignment_data.txt}.
#' The number of markers used will be appended to the name of the file.
#' @param keep.gsi.files (Boolean) Default \code{FALSE} The input and output gsi_sim files
#' will be deleted from the directory when finished processing.
#' With \code{TRUE}, remember to allocate a large chunk of the disk space for the analysis.
#' @param pop.levels (required) A character string with your populations ordered.
#' @param pop.labels (optional) A character string for your populations labels.
#' If you need to rename sampling sites in \code{pop.levels} or combined sites/pop
#' into a different names, here is the place.
#' @param pop.id.start The start of your population id
#' in the name of your individual sample. Your individuals are identified 
#' in this form : SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020,
#' then, \code{pop.id.start} = 5. If you didn't name your individuals
#' with the pop id in it, use the \code{strata} argument. 
#' @param pop.id.end The end of your population id
#' in the name of your individual sample. Your individuals are identified 
#' in this form : SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020,
#' then, \code{pop.id.end} = 7. If you didn't name your individuals
#' with the pop id in it, use the \code{strata} argument.
#' @param strata (optional) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. Default: \code{strata = NULL}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping.

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
#' @param iteration.subsample (Integer) The number of iterations to repeat 
#' subsampling, default: \code{iteration.subsample = 1}.
#' With \code{subsample = 20} and \code{iteration.subsample = 10},
#' 20 individuals/populations will be randomly chosen 10 times.


#' @param imputation.method Should a map-independent imputations of markers be
#' computed. Available choices are: (1) \code{FALSE} for no imputation.
#' (2) \code{"max"} to use the most frequent category for imputations.
#' (3) \code{"rf"} using Random Forest algorithm. 
#' Default: \code{imputation.method = FALSE}.
#' @param impute (character) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
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
#' @details You need to have either the \code{pop.id.start} and \code{pop.id.end}
#' or the \code{strata} argument, to identify your populations.
#' 
#' The imputations using Random Forest requires more time to compute
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

#' @note \code{assignment_ngs} assumes that the command line version of gsi_sim 
#' is properly installed into \code{file.path(system.file(package = "assigner"), "bin", "gsi_sim")}.
#' Things are set up so that it will try running gsi_sim, and if it does not find it, the 
#' program will throw an error and ask the user to run \code{\link{install_gsi_sim}}
#' which will do its best to put a usable copy of gsi_sim where it is needed.  To do 
#' so, you must be connected to the internet. If that doesn't work, you will
#' need to compile the program yourself, or get it yourself, and the manually copy
#' it to \code{file.path(system.file(package = "assigner"), "bin", "gsi_sim")}.
#' To compile gsi_sim, follow the 
#' instruction here: \url{https://github.com/eriqande/gsi_sim}.

#' @export
#' @rdname assignment_ngs
#' @import dplyr
#' @import parallel
#' @import stringi
#' @import adegenet
#' @importFrom purrr map
#' @importFrom purrr flatten
#' @importFrom purrr keep
#' @importFrom purrr discard
#' @importFrom data.table fread

#' @examples
#' \dontrun{
#' assignment.treefrog <- assignment_ngs(
#' data = "batch_1.vcf",
#' assignment.analysis = "gsi_sim",
#' whitelist.markers = "whitelist.vcf.txt",
#' snp.ld = NULL,
#' common.markers = TRUE,
#' marker.number = c(500, 5000, "all"),
#' sampling.method = "ranked",
#' thl = 0.3,
#' blacklist.id = "blacklist.id.treefrog.tsv",
#' subsample = 25,
#' iteration.subsample = 10
#' gsi_sim.filename = "treefrog.txt",
#' keep.gsi.files = FALSE,
#' pop.levels = c("PAN", "COS")
#' pop.id.start = 5, pop.id.end = 7,
#' imputation.method = FALSE,
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
#' @references Danecek P, Auton A, Abecasis G et al. (2011)
#' The variant call format and VCFtools.
#' Bioinformatics, 27, 2156-2158.
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, 
#' Bender D, et al. 
#' PLINK: a tool set for whole-genome association and population-based linkage 
#' analyses. 
#' American Journal of Human Genetics. 2007; 81: 559–575. doi:10.1086/519795
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

# required to pass the R CMD check and have 'no visible binding for global variable'
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c("ID", "#CHROM", "CHROM", "FORMAT", "INDIVIDUALS", "FORMAT_ID", "LOCUS",
      "POS", "REF", "ALT", "POP_ID", "READ_DEPTH", "ALLELE_DEPTH", "GL",
      "ERASE", "GT", "MARKERS", "QQ", "PQ", "N", "MAF_GLOBAL", "MAF_LOCAL",
      "ALLELES", "POP_ID", "GT", "INDIVIDUALS", "MARKERS", "POP_ID", "nal",
      "ALLELES_GROUP", "ALLELES", "N_IND_GENE", "P", "N", "nal_sq",
      "nal_sq_sum", "nal_sq_sum_nt", "npl", "het", "mho", "mhom", "dum",
      "dum1", "SSG", "ntal", "SSP", "ntalb", "SSi", "MSI", "sigw", "MSP",
      "siga", "sigb", "lsiga", "lsigb", "lsigw", "FST", "MARKERS",
      "MARKERS_ALLELES", "ALLELES", "POP_ID", "INDIVIDUALS", "filename",
      "ID", "KEEPER", "ASSIGN", "OTHERS", "CURRENT", "INFERRED",
      "SECOND_BEST_POP", "SCORE", "SECOND_BEST_SCORE", "NUMBER", "INDIVIDUALS_ALLELES",
      "MARKER_NUMBER", "MISSING_DATA", "TOTAL", "ASSIGNMENT_PERC",
      "MARKERS", "CURRENT", "INFERRED", "MISSING_DATA",
      "ITERATIONS", "METHOD", "TOTAL", "MEAN_i", "MEAN", "ASSIGNMENT_PERC",
      "SE", "MEDIAN", "MIN", "MAX", "QUANTILE25", "QUANTILE75", "SE_MIN",
      "SE_MAX", ".", "QUAL", "FILTER", "INFO", "pb", "SUBSAMPLE", "STRATA", 
      "sum.pop", "A1", "A2", "INDIVIDUALS_2", "Cnt", "Catalog ID", "GROUP",
      "COUNT", "MAX_COUNT_MARKERS", "hierarchy"
    )
  )
}

assignment_ngs <- function(data,
                           assignment.analysis,
                           whitelist.markers = NULL,
                           monomorphic.out = TRUE,
                           blacklist.genotype = NULL,
                           snp.ld = NULL,
                           common.markers = NULL,
                           maf.thresholds = NULL,
                           maf.pop.num.threshold = 1,
                           maf.approach = "SNP",
                           maf.operator = "OR",
                           max.marker = NULL,
                           marker.number = "all",
                           blacklist.id = NULL,
                           pop.levels,
                           pop.labels,
                           pop.id.start, 
                           pop.id.end,
                           strata = NULL,
                           pop.select = NULL,
                           subsample = NULL,
                           iteration.subsample = 1,
                           sampling.method,
                           thl = 1,
                           iteration.method = 10,
                           folder,
                           gsi_sim.filename = "assignment_data.txt",
                           keep.gsi.files,
                           imputation.method = FALSE,
                           impute = "genotypes",
                           imputations.group = "populations",
                           num.tree = 100,
                           iteration.rf = 10,
                           split.number = 100,
                           verbose = FALSE,
                           parallel.core = NULL,
                           ...) {
  
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  if (missing(assignment.analysis)) stop("assignment.analysis argument missing")
  if (assignment.analysis == "gsi_sim" & !gsi_sim_exists()){
    stop("Can't find the gsi_sim executable where it was expected at ", gsi_sim_binary_path(), ".  
         If you have internet access, you can install it
         from within R by invoking the function \"install_gsi_sim(fromSource = TRUE)\"")
  }
  if (assignment.analysis == "gsi_sim") message("Assignment analysis with gsi_sim")
  if (assignment.analysis == "adegenet") message("Assignment analysis with adegenet")
  if (missing(whitelist.markers)) whitelist.markers <- NULL # no Whitelist
  if (missing(monomorphic.out)) monomorphic.out <- TRUE # remove monomorphic
  if (missing(blacklist.genotype)) blacklist.genotype <- NULL # no genotype to erase
  if (missing(snp.ld)) snp.ld <- NULL
  if (missing(common.markers)) common.markers <- TRUE
  if (missing(maf.thresholds)) maf.thresholds <- NULL
  if (missing(maf.pop.num.threshold)) maf.pop.num.threshold <- 1
  if (missing(maf.approach)) maf.approach <- "SNP"
  if (missing(maf.operator)) maf.operator <- "OR"
  if (missing(max.marker)) max.marker <- NULL
  if (missing(marker.number)) marker.number <- "all"
  if (missing(blacklist.id)) blacklist.id <- NULL # No blacklist of ID
  if (missing(pop.levels)) stop("pop.levels required")
  if (missing(pop.labels)) pop.labels <- pop.levels # pop.labels
  if (missing(pop.id.start)) pop.id.start <- NULL
  if (missing(pop.id.end)) pop.id.end <- NULL
  if (missing(strata)) strata <- NULL
  if (missing(pop.select)) pop.select <- NULL
  if (missing(subsample)) subsample <- NULL
  if (missing(iteration.subsample)) iteration.subsample <- 1
  if (missing(sampling.method)) stop("Sampling method required")
  if (sampling.method == "ranked" & missing(thl)) thl <- 1 # thl
  if (missing(iteration.method)) iteration.method <- 10
  if (sampling.method  == "ranked" & thl == "all") iteration.method <- 1
  if (sampling.method == "ranked" & thl == 1) iteration.method <- 1
  if (missing(gsi_sim.filename)) gsi_sim.filename <- "assignment_data.txt"
  if (missing(keep.gsi.files)) keep.gsi.files <- FALSE
  if (missing(imputation.method)) imputation.method <- FALSE
  if (missing(imputations.group)) imputations.group <- "populations"
  if (imputation.method != FALSE & missing(impute)) stop("impute argument is necessary")
  if (imputation.method == FALSE & missing(impute)) impute <- NULL
  if (missing(num.tree)) num.tree <- 100
  if (missing(iteration.rf)) iteration.rf <- 10
  if (missing(split.number)) split.number <- 100
  if (missing(verbose)) verbose <- FALSE
  if (missing(parallel.core) | is.null(parallel.core)) parallel.core <- detectCores()-1
  if (missing(folder)) folder <- NULL
  
  # File type detection ********************************************************
  data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
  
  if (identical(data.type, "##fileformat=VCF") | stri_detect_fixed(str = data, pattern = ".vcf")) {
    data.type <- "vcf.file"
    message("File type: VCF")
  }
  
  if (stri_detect_fixed(str = data, pattern = ".tped")) {
    data.type <- "plink.file"
    message("File type: PLINK")
    if (!file.exists(stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
      stop("Missing tfam file with the same prefix as your tped")
    }
  } 
  
  if (stri_detect_fixed(str = data.type, pattern = "POP_ID") | stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS")) {
    data.type <- "df.file"
    message("File type: data frame of genotypes")
  }
  
  if (stri_detect_fixed(str = data.type, pattern = "Catalog")) {
    data.type <- "haplo.file"
    message("File type: haplotypes from stacks")
    if (is.null(blacklist.genotype)) {
      stop("blacklist.genotype file missing. 
              Use stackr's missing_genotypes function to create this blacklist")
    }
  }
  
  
  # Create a folder based on filename to save the output files *****************
  if (is.null(folder)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    if (imputation.method == "FALSE") {
      message("Map-imputation: no")
      directory <- stri_join(getwd(),"/", "assignment_analysis_", "method_", sampling.method, "_no_imputations_", file.date, "/", sep = "")
      dir.create(file.path(directory))
    } else {
      message("Map-imputation: yes")
      directory <- stri_join(getwd(),"/","assignment_analysis_", "method_", sampling.method, "_imputations_", imputation.method,"_", imputations.group, "_", file.date, "/", sep = "")
      dir.create(file.path(directory))
    }
    message(stri_join("Folder: ", directory))
    file.date <- NULL #unused object
  } else {
    directory <- stri_join(getwd(), "/", folder, "/", sep = "")
    dir.create(file.path(directory))
    message(stri_join("Folder: ", directory))
  }
  
  # Import whitelist of markers ************************************************
  if (is.null(whitelist.markers)) { # no Whitelist
    message("Whitelist of markers: no")
  } else { # with Whitelist of markers
    message("Whitelist of markers: yes")
    whitelist.markers <- read_tsv(whitelist.markers, col_names = TRUE)
    columns.names.whitelist <- colnames(whitelist.markers)
    if ("CHROM" %in% columns.names.whitelist) {
      whitelist.markers$CHROM <- as.character(whitelist.markers$CHROM)
    }
    if ("LOCUS" %in% columns.names.whitelist) {
      whitelist.markers$LOCUS <- as.character(whitelist.markers$LOCUS)
    }
    if ("POS" %in% columns.names.whitelist) {
      whitelist.markers$POS <- as.character(whitelist.markers$POS)
    }
  }
  
  if (data.type == "haplo.file") {
    whitelist.markers <- select(.data = whitelist.markers, LOCUS)
    columns.names.whitelist <- colnames(whitelist.markers)
  }
  
  # Import blacklist id ********************************************************
  if (is.null(blacklist.id)) { # No blacklist of ID
    message("Blacklisted individuals: no")
  } else { # With blacklist of ID
    message("Blacklisted individuals: yes")
    blacklist.id <- read_tsv(blacklist.id, col_names = TRUE)
  }
  
  # Import data ****************************************************************
  if (data.type == "vcf.file") { # VCF
    message("Importing the VCF...")
    
    if (is.null(strata) & is.null(pop.id.start) & is.null(pop.id.end)) {
      stop("pop.id.start and pop.id.end or strata arguments are required to 
           identify your populations with a VCF file")
    }
    
    input <- data.table::fread(
      input = data,
      sep = "\t",
      stringsAsFactors = FALSE, 
      header = TRUE,
      # Automatically filter with blacklist of id
      drop = c("QUAL", "FILTER", "INFO", blacklist.id$INDIVIDUALS),
      skip = "CHROM",
      showProgress = TRUE,
      verbose = FALSE
    ) %>% 
      as_data_frame() %>% 
      rename(LOCUS = ID, CHROM = `#CHROM`) %>%
      mutate(
        CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1"),
        POS = as.character(POS),
        LOCUS = as.character(LOCUS)
      )
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Detect STACKS version
    if (stri_detect_fixed(input$FORMAT[1], "AD")) {
      stacks.version <- "new"
    } else {
      stacks.version <- "old"
    }
    input <- input %>% select(-FORMAT)
    
    # Tidying the VCF to make it easy to work on the data for conversion
    message("Making the VCF population wise")
    input <- input %>%
      tidyr::gather(INDIVIDUALS, FORMAT_ID, -c(CHROM, LOCUS, POS, REF, ALT)) # Gather individuals in 1 colummn
    
    if (is.null(strata)){
      input <- input %>%
        mutate( # Make population ready
          POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
          POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered = TRUE),
          INDIVIDUALS =  as.character(INDIVIDUALS)
        )
    } else { # Make population ready with the strata provided
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA)
      
      input <- input %>%
        mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        left_join(strata.df, by = "INDIVIDUALS") %>% 
        mutate(POP_ID = factor(POP_ID, levels = unique(pop.labels), ordered =TRUE))
    }
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # Tidy VCF
    message("Tidying the vcf...")
    if (stacks.version == "new") { # with new version of stacks > v.1.29
      input <- input %>%
        tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "ALLELE_DEPTH", "GL"),
                        sep = ":", extra = "warn") %>%
        select(-c(READ_DEPTH, ALLELE_DEPTH, GL))
    } else { # stacks version prior to v.1.29 had no Allele Depth field...
      input <- input %>%
        tidyr::separate(FORMAT_ID, c("GT", "READ_DEPTH", "GL"),
                        sep = ":", extra = "warn") %>%
        select(-c(READ_DEPTH, GL))
    }
  } # End import VCF
  
  if (data.type == "plink.file") { # PLINK
    message("Importing the PLINK files...")
    strata.df <- data.table::fread(
      input = stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE),
      sep = " ", 
      header = FALSE, 
      stringsAsFactors = FALSE,
      verbose = FALSE,
      select = c(1,2),
      colClasses=list(character = c(1,2)),
      col.names = c("POP_ID", "INDIVIDUALS"),
      showProgress = TRUE, 
      data.table = FALSE)
    
    # remove "_" in individual name and replace with "-"
    strata.df$INDIVIDUALS <- stri_replace_all_fixed(str = strata.df$INDIVIDUALS, pattern = c("_", ":"), replacement = c("-", "-"), vectorize_all = FALSE)
    
    tped.header.prep <- strata.df %>% 
      select(INDIVIDUALS) %>%
      mutate(NUMBER = seq(1,n())) %>%
      mutate(ALLELE1 = rep("A1", n()), ALLELE2 = rep("A2", n())) %>%
      tidyr::gather(ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, NUMBER)) %>%
      arrange(NUMBER) %>% 
      select(-ALLELES_GROUP) %>% 
      tidyr::unite(INDIVIDUALS_ALLELES, c(INDIVIDUALS, ALLELES), sep = "_", remove = FALSE) %>% 
      arrange(NUMBER) %>% 
      mutate(NUMBER = seq(from = (1 + 4), to = n() + 4)) %>% 
      select(-ALLELES)
    
    tped.header.names <- c("LOCUS", tped.header.prep$INDIVIDUALS_ALLELES)
    tped.header.integer <- c(2, tped.header.prep$NUMBER)
    
    if (!is.null(blacklist.id)) { # using the blacklist of individuals
      whitelist.id <- tped.header.prep %>% 
        anti_join(blacklist.id, by = "INDIVIDUALS") %>% 
        arrange(NUMBER)
      tped.header.names <- c("LOCUS", whitelist.id$INDIVIDUALS_ALLELES)
      tped.header.integer <- c(2, whitelist.id$NUMBER)
    }
    
    input <- data.table::fread( # import PLINK
      input = data, 
      sep = " ", 
      header = FALSE, 
      stringsAsFactors = FALSE, 
      verbose = FALSE,
      select = tped.header.integer,
      col.names = tped.header.names,
      showProgress = TRUE,
      data.table = FALSE) %>% 
      as_data_frame() %>% 
      mutate(LOCUS = as.character(LOCUS))
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # To reduce the size of the dataset we subsample the markers with max.marker
    if (!is.null(max.marker)) {
      message("Using the max.marker to reduce the size of the dataset")
      input <- sample_n(tbl = input, size = max(as.numeric(max.marker)), replace = FALSE)
      
      max.marker.subsample.select <- input %>% 
        select(LOCUS) %>% 
        distinct(LOCUS) %>% 
        arrange(LOCUS)
      
      write_tsv(# save results
        x = max.marker.subsample.select, 
        path = paste0(directory,"max.marker.subsample.select.tsv")
      )
    }
    
    # Using the argument strata if provided to replace the current one
    if (!is.null(strata)) {
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA) %>% 
        mutate(INDIVIDUALS = stri_replace_all_fixed(str = INDIVIDUALS, 
                                                    pattern = c("_", ":"), 
                                                    replacement = c("-", "-"),
                                                    vectorize_all = FALSE)
        )
    }
    
    # Make tidy
    message("Tidying the PLINK file and integrating the tfam/strata file, for large dataset this may take several minutes...")
    input <- input %>% 
      tidyr::gather(key = INDIVIDUALS_ALLELES, value = GT, -LOCUS) %>%
      mutate(INDIVIDUALS = stri_replace_all_fixed(str = INDIVIDUALS_ALLELES, pattern = c("_A1", "_A2"), replacement = "", vectorize_all = FALSE)) %>% 
      left_join(strata.df, by = "INDIVIDUALS") %>% 
      mutate(
        POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered = TRUE),
        GT = stri_pad_left(str = GT, width = 3, pad = "0")
      )
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
    
    # removing untyped markers across all-pop
    remove.missing.gt <- input %>%
      select(LOCUS, GT) %>%
      filter(GT != "000")
    
    untyped.markers <- n_distinct(input$LOCUS) - n_distinct(remove.missing.gt$LOCUS)
    if (untyped.markers > 0) {
      message(paste0("Number of marker with 100 % missing genotypes: ", untyped.markers))
      input <- suppressWarnings(
        semi_join(input, 
                  remove.missing.gt %>% 
                    select(LOCUS) %>% 
                    distinct(LOCUS), 
                  by = "LOCUS")
      )
    }
    
    # Removing monomorphic markers
    if (monomorphic.out == TRUE) {
      message("Removing monomorphic markers...")
      mono.markers <- remove.missing.gt %>%
        group_by(LOCUS, GT) %>%
        distinct %>% 
        group_by(LOCUS) %>%
        tally %>% 
        filter(n == 1) %>% 
        select(LOCUS) %>% 
        arrange(LOCUS)
      
      # Remove the markers from the dataset
      input <- anti_join(input, mono.markers, by = "LOCUS")
      message(paste0("Number of monomorphic markers removed: ", n_distinct(mono.markers$LOCUS)))
    }
    # Unused objects
    tped.header.prep <- NULL
    tped.header.integer <- NULL
    tped.header.names <- NULL
    remove.missing.gt <- NULL
    mono.markers <- NULL
  } # End import PLINK
  
  if (data.type == "df.file") { # DATA FRAME OF GENOTYPES
    message("Importing the data frame")
    input <- data.table::fread(
      input = data, 
      sep = "\t", 
      header = TRUE, 
      stringsAsFactors = FALSE, 
      verbose = FALSE,
      showProgress = TRUE,
      data.table = FALSE
    ) %>% 
      as_data_frame() %>% 
      tidyr::gather(key = LOCUS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% 
      mutate(
        GT = as.character(GT),
        GT = stri_pad_left(str= GT, pad = "0", width = 6),
        INDIVIDUALS = stri_replace_all_fixed(str = INDIVIDUALS, 
                                                    pattern = c("_", ":"), 
                                                    replacement = c("-", "-"),
                                                    vectorize_all = FALSE)
        )
    
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      input <- suppressWarnings(anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    
    # Using the argument strata if provided to replace the current one
    if (is.null(strata)) {
      input <- input %>%
        mutate( # Make population ready
          POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered = TRUE),
          INDIVIDUALS =  as.character(INDIVIDUALS)
        )
      
      strata.df <- input %>% 
        select(INDIVIDUALS, POP_ID) %>% 
        distinct(INDIVIDUALS)
    } else {
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA) %>% 
        mutate(
          INDIVIDUALS = stri_replace_all_fixed(str = INDIVIDUALS, 
                                               pattern = c("_", ":"), 
                                               replacement = c("-", "-"),
                                               vectorize_all = FALSE)
        )
      
      input <- input %>%
        mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        select(-POP_ID) %>% 
        left_join(strata.df, by = "INDIVIDUALS") %>% 
        mutate(POP_ID = factor(POP_ID, levels = unique(pop.labels), ordered =TRUE))
    }
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
  } # End import data frame of genotypes
  
  if (data.type == "haplo.file") { # Haplotype file
    message("Importing the stacks haplotype file")
    input <- data.table::fread(
      input = data, 
      sep = "\t", 
      header = TRUE, 
      stringsAsFactors = FALSE, 
      verbose = FALSE,
      showProgress = TRUE,
      data.table = FALSE
    ) %>% 
      as_data_frame() %>% 
      select(-Cnt) %>% 
      rename(LOCUS = `Catalog ID`) %>%
      tidyr::gather(INDIVIDUALS, GT, -LOCUS) %>% 
      mutate(
        LOCUS = as.character(LOCUS),
        INDIVIDUALS = stri_replace_all_fixed(str = INDIVIDUALS, 
                                             pattern = c("_", ":"), 
                                             replacement = c("-", "-"),
                                             vectorize_all = FALSE)
      )
    
    # Filter with whitelist of markers
    if (!is.null(whitelist.markers)) {
      message("Filtering with whitelist of markers")
      input <- suppressWarnings(semi_join(input, whitelist.markers, by = columns.names.whitelist))
    }
    
    # Filter with blacklist of individuals
    if (!is.null(blacklist.id)) {
      message("Filtering with blacklist of individuals")
      input <- suppressWarnings(anti_join(input, blacklist.id, by = "INDIVIDUALS"))
    }
    
    # Using the argument strata if provided to replace the current one
    if (is.null(strata)){
      input <- input %>%
        mutate( # Make population ready
          POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
          POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered = TRUE),
          INDIVIDUALS =  as.character(INDIVIDUALS)
        )
    } else { # Make population ready with the strata provided
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA) %>% 
        mutate(
          INDIVIDUALS = stri_replace_all_fixed(str = INDIVIDUALS, 
                                               pattern = c("_", ":"), 
                                               replacement = c("-", "-"),
                                               vectorize_all = FALSE)
        )
      
      input <- input %>%
        mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
        left_join(strata.df, by = "INDIVIDUALS") %>% 
        mutate(POP_ID = factor(POP_ID, levels = unique(pop.labels), ordered =TRUE))
    }
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      input <- suppressWarnings(input %>% filter(POP_ID %in% pop.select))
    }
  } # End import haplotypes file
  
  # Blacklist genotypes ********************************************************
  if (is.null(blacklist.genotype)) { # no Whitelist
    message("Erasing genotype: no")
  } else {
    message("Erasing genotype: yes")
    blacklist.genotype <- read_tsv(blacklist.genotype, col_names = TRUE)
    blacklist.genotype <- suppressWarnings(
      plyr::colwise(as.character, exclude = NA)(blacklist.genotype)
    )
    columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    
    if ("CHROM" %in% columns.names.blacklist.genotype) {
      columns.names.blacklist.genotype$CHROM <- as.character(columns.names.blacklist.genotype$CHROM)
    }
    
    if (data.type == "haplo.file") {
      blacklist.genotype <- select(.data = blacklist.genotype, INDIVIDUALS, LOCUS)
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }
    
    # control check to keep only individuals in pop.select
    if (!is.null(pop.select)) {
      message("Control check to keep only individuals present in pop.select")
      # updating the blacklist.genotype
      if (is.null(strata)){
        blacklist.genotype <- suppressWarnings(
          blacklist.genotype  %>% 
            mutate( # Make population ready
              POP_ID = substr(INDIVIDUALS, pop.id.start, pop.id.end),
              POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = F), levels = unique(pop.labels), ordered =T),
              INDIVIDUALS =  as.character(INDIVIDUALS) 
            ) %>% 
            filter(POP_ID %in% pop.select) %>% 
            select(-POP_ID)
        )
      } else {
        blacklist.genotype <- suppressWarnings(
          blacklist.genotype %>%
            mutate(INDIVIDUALS =  as.character(INDIVIDUALS)) %>% 
            left_join(strata.df, by = "INDIVIDUALS") %>% 
            filter(POP_ID %in% pop.select) %>% 
            select(-POP_ID)
        )
      }
    }
    
    # control check to keep only whitelisted markers from the blacklist of genotypes
    if (!is.null(whitelist.markers)) {
      blacklist.genotype <- blacklist.genotype
      message("Control check to keep only whitelisted markers present in the blacklist of genotypes to erase.")
      # updating the whitelist of markers to have all columns that id markers
      if (data.type == "vcf.file"){
        whitelist.markers.ind <- input %>% select(CHROM, LOCUS, POS, INDIVIDUALS) %>% distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      } else {
        whitelist.markers.ind <- input %>% select(LOCUS, INDIVIDUALS) %>% distinct(LOCUS, INDIVIDUALS)
      }
      
      # updating the blacklist.genotype
      blacklist.genotype <- suppressWarnings(semi_join(whitelist.markers.ind, blacklist.genotype, by = columns.names.blacklist.genotype))
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }
    
    # control check to remove blacklisted individuals from the blacklist of genotypes
    if (!is.null(blacklist.id)) {
      message("Control check to remove blacklisted individuals present in the blacklist of genotypes to erase.")
      blacklist.genotype <- suppressWarnings(anti_join(blacklist.genotype, blacklist.id, by = "INDIVIDUALS"))
      columns.names.blacklist.genotype <- colnames(blacklist.genotype)
    }
    
    # Add one column that will allow to include the blacklist in the dataset 
    # by x column(s) of markers
    blacklist.genotype <- mutate(.data = blacklist.genotype, ERASE = rep("erase", n()))
    
    input <- suppressWarnings(
      input %>%
        full_join(blacklist.genotype, by = columns.names.blacklist.genotype) %>%
        mutate(ERASE = stri_replace_na(str = ERASE, replacement = "ok"))
    )
    
    if (data.type == "vcf.file") {
      input <- input %>% 
        mutate(GT = ifelse(ERASE == "erase", "./.", GT)) %>% 
        select(-ERASE)
    }
    
    if (data.type == "plink.file") {
      input <- input %>% 
        mutate(GT = ifelse(ERASE == "erase", "000", GT)) %>% 
        select(-ERASE)
    }
    
    if (data.type == "df.file") {
      input <- input %>% 
        mutate(GT = ifelse(ERASE == "erase", "000000", GT)) %>% 
        select(-ERASE)
    }
    
    if (data.type == "haplo.file") {
      input <- input %>% 
        mutate(GT = ifelse(ERASE == "erase", "-", GT)) %>% 
        select(-ERASE)
    }
    
  } # End erase genotypes
  
  # dump unused object
  blacklist.id <- NULL
  whitelist.markers <- NULL
  whitelist.markers.ind <- NULL
  blacklist.genotype <- NULL
  
  # subsampling data ***********************************************************
  # Function:
  subsampling_data <- function(iteration.subsample, ...) {
    # message(paste0("Creating data subsample: ", iteration.subsample))
    if (is.null(subsample)) {
      subsample.select <- ind.pop.df %>% 
        mutate(SUBSAMPLE = rep(iteration.subsample, n()))
    } else {
      if (subsample > 1) { # integer
        subsample.select <- ind.pop.df %>%
          group_by(POP_ID) %>%
          sample_n(subsample, replace = FALSE) %>% # sampling individuals for each pop
          arrange(POP_ID, INDIVIDUALS) %>% 
          mutate(SUBSAMPLE = rep(iteration.subsample, n())) %>% 
          ungroup()
      }
      if (subsample < 1) { # proportion
        subsample.select <- ind.pop.df %>%
          group_by(POP_ID) %>%
          sample_frac(subsample, replace = FALSE) %>% # sampling individuals for each pop
          arrange(POP_ID, INDIVIDUALS) %>% 
          mutate(SUBSAMPLE = rep(iteration.subsample, n())) %>% 
          ungroup()
      }
    }
    return(subsample.select)
  } # End subsampling function
  # create the subsampling list
  ind.pop.df <- input %>% select(POP_ID, INDIVIDUALS) %>% distinct(POP_ID, INDIVIDUALS)
  subsample.list <- map(.x = 1:iteration.subsample, .f = subsampling_data, subsample = subsample)
  
  # keep track of subsampling individuals and write to directory
  if (is.null(subsample)) {
    message("Subsampling: not selected")
  } else {
    message("Subsampling: selected")
    subsampling.individuals <- bind_rows(subsample.list)
    write_tsv(x = subsampling.individuals, path = paste0(directory, "subsampling.individuals.tsv"), col_names = TRUE, append = FALSE)
  } # End subsampling
  
  # unused objects
  subsampling.individuals <- NULL
  ind.pop.df <- NULL
  
  # assignment analysis ********************************************************
  # Function:
  assignment_function <- function(data, ...) {
    # data <- subsample.list[[1]] # test
    subsampling.individuals <- data
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
    input <- input %>%
      semi_join(subsampling.individuals, by = c("POP_ID", "INDIVIDUALS"))
    
    # unused object
    data <- NULL
    subsampling.individuals <- NULL
    
    # LD control... keep only 1 SNP per haplotypes/reads (optional) ************
    if (!is.null(snp.ld)) {
      if (data.type != "vcf.file") {
        stop("snp.ld is only available for VCF file, use stackr package for 
haplotype file and create a whitelist, for other file type, use 
             PLINK linkage disequilibrium based SNP pruning option")
      }
      message("Minimizing LD...")
      snp.locus <- input %>% select(LOCUS, POS) %>% distinct(POS)
      # Random selection
      if (snp.ld == "random") {
        snp.select <- snp.locus %>%
          group_by(LOCUS) %>%
          sample_n(size = 1, replace = FALSE)
        message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP randomly selected to keep 1 SNP per read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
      }
      
      # Fist SNP on the read
      if (snp.ld == "first") {
        snp.select <- snp.locus %>%
          group_by(LOCUS) %>%
          summarise(POS = min(POS))
        message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
      }
      
      # Last SNP on the read
      if (snp.ld == "last") {
        snp.select <- snp.locus %>%
          group_by(LOCUS) %>%
          summarise(POS = max(POS))
        message(stri_join("Number of original SNP = ", n_distinct(snp.locus$POS), "\n", "Number of SNP after keeping the first SNP on the read/haplotype = ", n_distinct(snp.select$POS), "\n", "Number of SNP removed = ", n_distinct(snp.locus$POS) - n_distinct(snp.select$POS)))
      }
      
      # filtering the VCF to minimize LD
      input <- input %>% semi_join(snp.select, by = c("LOCUS", "POS"))
      message("Filtering the tidy VCF to minimize LD by keeping only 1 SNP per short read/haplotype")
    } # End of snp.ld control
    
    # Unique markers id: for VCF combine CHROM, LOCUS and POS into MARKERS *****
    # For other type of file change LOCUS to MARKERS
    if (data.type != "vcf.file") {
      input <- input %>% rename(MARKERS = LOCUS)
    } else {
      input <- input %>%
        mutate(
          POS = stri_pad_left(str = POS, width = 8, pad = "0"),
          LOCUS = stri_pad_left(str = LOCUS, width = 8, pad = "0")
        ) %>%
        arrange(CHROM, LOCUS, POS) %>%
        tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "_")
    } # End Unique markers id
    
    # Markers in common between all populations (optional) *********************
    if (common.markers == TRUE) { # keep only markers present in all pop
      message("Using markers common in all populations:")
      pop.number <- n_distinct(input$POP_ID)
      
      if (data.type == "vcf.file") pop.filter <- input %>% filter(GT != "./.")
      if (data.type == "plink.file") pop.filter <- input %>% filter(GT != "000")
      if (data.type == "df.file") pop.filter <- input %>% filter(GT != "000000")
      if (data.type == "haplo.file") pop.filter <- input %>% filter(GT != "-")
      
      pop.filter <- pop.filter %>% 
        group_by(MARKERS) %>%
        filter(n_distinct(POP_ID) == pop.number) %>%
        arrange(MARKERS) %>%
        select(MARKERS) %>%
        distinct(MARKERS)
      
      message(stri_join("Number of original markers = ", n_distinct(input$MARKERS), 
                        "\n", "Number of markers present in all the populations = ", 
                        n_distinct(pop.filter$MARKERS), "\n", 
                        "Number of markers removed = ", 
                        n_distinct(input$MARKERS) - n_distinct(pop.filter$MARKERS))
      )
      input <- suppressWarnings(input %>% semi_join(pop.filter, by = "MARKERS"))
      pop.filter <- NULL # ununsed object
    } # End common markers
    
    # Minor Allele Frequency filter ********************************************
    # maf.thresholds <- c(0.05, 0.1) # test
    if (!is.null(maf.thresholds)) { # with MAF
      maf.local.threshold <- maf.thresholds[1]
      maf.global.threshold <- maf.thresholds[2]
      message("MAF filter: yes")
      
      if (data.type == "vcf.file") {
        maf.local <- input %>%
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
        
        maf.local <- NULL
        maf.global <- NULL
      } # end maf calculations with vcf
      
      if (data.type == "plink.file" | data.type == "df.file") {
        message("Calculating global and local MAF, this may take some time on large data set")
        
        # For data frame we split the alleles here to prep for MAF
        if (data.type == "df.file") { # for data frame of genotypes
          maf.data <- input %>% 
            tidyr::separate(data = ., col = GT, into = .(A1, A2), sep = 3, remove = TRUE) %>% 
            tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
            select(MARKERS, GT, POP_ID) %>% 
            filter(GT != "000")
        }
        
        if (data.type == "plink.file") { # For PLINK and common code below
          maf.data <- input %>%
            select(MARKERS, GT, POP_ID) %>% 
            filter(GT != "000")
        }
        
        maf.data <- maf.data %>%
          group_by(MARKERS, GT, POP_ID) %>%
          tally %>%
          arrange(MARKERS, GT) %>% 
          group_by(MARKERS, GT) %>%
          mutate(sum.pop = sum(n)) %>% 
          group_by(MARKERS) %>%
          mutate(MAF_GLOBAL = min(sum.pop)/sum(n)) %>% 
          group_by(MARKERS, POP_ID) %>%
          mutate(MAF_LOCAL = n/sum(n)) %>% 
          arrange(MARKERS, POP_ID, GT) %>% 
          group_by(MARKERS, POP_ID) %>% 
          filter(n == min(n)) %>% 
          distinct(MARKERS, POP_ID) %>% 
          select(MARKERS, POP_ID, MAF_LOCAL, MAF_GLOBAL)
      }# end maf calculations with PLINK or data frame of genotypes
      
      if (data.type == "haplo.file") {
        stop("MAF filtering is only available for haplotype file, use stackr
package and update your whitelist")
      }
      
      write_tsv(x = maf.data, 
                path = paste0(directory.subsample,"maf.data.tsv"), 
                col_names = TRUE, 
                append = FALSE
      )
      message("The MAF table was written in your folder")
      
      # # update the vcf with the maf info
      # input <- full_join(input, maf.data, by = c("MARKERS", "POP_ID"))
      if (maf.approach == "haplotype") {
        if (data.type != "vcf.file") {
          stop("The haplotype approach during MAF filtering is for VCF files only")
        }
        vcf.maf <- tidyr::separate(data = maf.data, 
                                   col = MARKERS, 
                                   into = c("CHROM", "LOCUS", "POS"), 
                                   sep = "_", 
                                   remove = FALSE, 
                                   extra = "warn"
        )
        
        if (maf.operator == "OR") {
          vcf.maf <- maf.data %>%
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
            left_join(input, by = "LOCUS") %>%
            arrange(LOCUS, POP_ID)
        } else { # AND operator between local and global maf
          vcf.maf <- maf.data %>%
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
            left_join(input, by = "LOCUS") %>%
            arrange(LOCUS, POP_ID)
        }
        vcf.maf <- vcf.maf %>% select(-c(CHROM, LOCUS, POS))
      } # end maf haplotype approach
      
      if (maf.approach == "SNP") { # SNP approach
        if (maf.operator == "OR") {
          vcf.maf <- maf.data %>%
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
            left_join(input, by = "MARKERS") %>%
            arrange(MARKERS, POP_ID)
        } else { # AND operator between local and global maf
          vcf.maf <- maf.data %>%
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
            left_join(input, by = "MARKERS") %>%
            arrange(MARKERS, POP_ID)
        }
      } # end maf snp approach
      
      
      message(stri_join("The number of MARKERS removed by the MAF filters = ", 
                        n_distinct(input$MARKERS)-n_distinct(vcf.maf$MARKERS), "\n", 
                        "The number of MARKERS before -> after the MAF filters: ", 
                        n_distinct(input$MARKERS)," -> ", n_distinct(vcf.maf$MARKERS), 
                        " MARKERS"))
      
      input <- vcf.maf
      
      # unused object
      vcf.maf <- NULL 
      maf.data <- NULL
    } # End of MAF filters
    
    # Change the genotype coding  **********************************************
    message("Recoding genotypes")
    # adegenet & gsi_sim
    
    if (data.type == "vcf.file" & assignment.analysis == "gsi_sim") { # for VCF input
      input <- input %>%
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
      
      gsi.prep <- input %>%
        tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
        tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID))
    }
    if (data.type == "vcf.file" & assignment.analysis == "adegenet") {
      genind.prep <- input %>%
        mutate(GT = stri_replace_all_fixed(str = GT, pattern = "/", replacement = ":", vectorize_all = FALSE)) %>%
        mutate(GT = stri_replace_all_fixed(str = GT, pattern = c("0:0", "1:1", "0:1", "1:0", ".:."), replacement = c("2_0", "0_2", "1_1", "1_1", "NA_NA"), vectorize_all = FALSE)) %>%
        select(-REF, -ALT) %>% 
        arrange(MARKERS, POP_ID) %>%
        tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_", extra = "drop", remove = TRUE) %>%
        tidyr::gather(key = ALLELES, value = COUNT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% # make tidy
        tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>%
        mutate(COUNT = replace(COUNT, which(COUNT == "NA"), NA)) %>% 
        group_by(POP_ID, INDIVIDUALS) %>%
        tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>%
        arrange(POP_ID, INDIVIDUALS)
    }
    
    if (data.type == "haplo.file") { # for haplotypes
      gsi.prep <- suppressWarnings(
        input %>%
          mutate(GT = stri_replace_all_fixed(GT, "-", "000/000", vectorize_all=F)) %>% 
          arrange(MARKERS) %>% 
          tidyr::separate(
            col = GT, into = c("A1", "A2"), 
            sep = "/", extra = "drop", remove = TRUE
          ) %>%
          mutate(
            A2 = stri_replace_na(str = A2, replacement = "000"),
            A2 = ifelse(A2 == "000", A1, A2)
          ) %>%
          tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
          ungroup()
      )
    }
    if (data.type == "haplo.file" & assignment.analysis == "gsi_sim") {
      input <- suppressWarnings(
        gsi.prep %>%
          filter(GT != "000") %>% 
          tidyr::spread(data = ., key = MARKERS, value = GT) %>% # this reintroduce the missing, but with NA
          ungroup() %>% 
          plyr::colwise(.fun = factor, exclude = NA)(.)
      )
      
      input <- suppressWarnings(
        input %>%
          ungroup() %>% 
          mutate_each(funs(as.integer), -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
          ungroup() %>% 
          tidyr::gather(data = ., key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
          mutate(
            GT = stri_replace_na(str = GT, replacement = "000"),
            GT = stri_pad_left(str = GT, width = 3, pad = "0")
          ) %>% 
          tidyr::spread(data = ., key = ALLELES, value = GT) %>% 
          tidyr::unite(data = ., GT, A1, A2, sep = "", remove = TRUE)
      )
      
      # for haplo.file we need to change back again the gsi.prep file
      gsi.prep <- input %>% 
        tidyr::separate(data = ., col = GT, into = .(A1, A2), sep = 3, remove = TRUE) %>% 
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) 
    }
    if (data.type == "df.file") { # For data frame of genotypes
      gsi.prep <- input %>% 
        tidyr::separate(data = ., col = GT, into = .(A1, A2), sep = 3, remove = TRUE) %>% 
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) 
    }
    if (data.type == "plink.file") { # for PLINK
      gsi.prep <- input %>% 
        tidyr::separate(col = INDIVIDUALS_ALLELES, into = c("INDIVIDUALS_2", "ALLELES"), sep = "_", remove = TRUE) %>% 
        select(-INDIVIDUALS_2)
    }
    if (data.type == "df.file" | data.type == "plink.file" | data.type == "haplo.file") {
      if (assignment.analysis == "adegenet" ) {
        genind.prep <- suppressWarnings(
          gsi.prep %>%
            filter(GT != "000") %>% 
            tidyr::spread(data = ., key = MARKERS, value = GT) %>% # this reintroduce the missing, but with NA
            ungroup() %>% 
            plyr::colwise(.fun = factor, exclude = NA)(.)
        )
        
        genind.prep <- suppressWarnings(
          genind.prep %>%
            ungroup() %>% 
            mutate_each(funs(as.integer), -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
            ungroup() %>% 
            tidyr::gather(data = ., key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
            mutate(GT = stri_replace_na(str = GT, replacement = "000")) %>%
            filter(GT != "000") %>%
            select(-ALLELES) %>%
            group_by(POP_ID, INDIVIDUALS, MARKERS, GT) %>% 
            tally %>%
            ungroup() %>%
            tidyr::unite(MARKERS_ALLELES, MARKERS, GT, sep = ":", remove = TRUE) %>%
            arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>% 
            group_by(POP_ID, INDIVIDUALS) %>% 
            tidyr::spread(data = ., key = MARKERS_ALLELES, value = n) %>%
            ungroup() %>%
            tidyr::gather(data = ., key = MARKERS_ALLELES, value = COUNT, -c(INDIVIDUALS, POP_ID)) %>% 
            tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = ":", remove = TRUE) %>% 
            mutate(COUNT = as.numeric(stri_replace_na(str = COUNT, replacement = "0"))) %>% 
            group_by(INDIVIDUALS, MARKERS) %>%
            mutate(MAX_COUNT_MARKERS = max(COUNT, na.rm = TRUE)) %>%
            ungroup() %>% 
            mutate(COUNT = ifelse(MAX_COUNT_MARKERS == 0, "erase", COUNT)) %>%
            select(-MAX_COUNT_MARKERS) %>% 
            mutate(COUNT = replace(COUNT, which(COUNT == "erase"), NA)) %>% 
            arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES) %>% 
            tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>%
            tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>% 
            arrange(POP_ID, INDIVIDUALS)
        )
      }
    }
    
    # gsi.prep <- gsi.prep %>% 
    #   arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES)
    
    # only adegenet
    if (assignment.analysis == "adegenet") {
      # genind arguments common to all data.type
      ind <- as.character(genind.prep$INDIVIDUALS)
      pop <- genind.prep$POP_ID
      genind.df <- genind.prep %>% ungroup() %>% 
        select(-c(INDIVIDUALS, POP_ID))
      rownames(genind.df) <- ind
      loc.names <- colnames(genind.df)
      strata <- genind.prep %>% ungroup() %>% select(INDIVIDUALS, POP_ID) %>% distinct(INDIVIDUALS, POP_ID)
      
      # genind constructor
      prevcall <- match.call()
      genind.object <- genind(tab = genind.df, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = NULL)
      
      # sum <- summary(genind.object) # test
      # sum$NA.perc # test
      
      ind <- NULL
      pop <- NULL
      genind.df <- NULL
      genind.prep <- NULL
    }
    
    # Imputations **************************************************************
    if (imputation.method != "FALSE") {
      message("Preparing the data for imputations")
      
      if (data.type == "vcf.file") { # for VCF input
        if (impute == "genotype") {
          input.prep <- input %>%
            select(-REF, -ALT) %>%
            mutate(GT = stri_replace_all_fixed(str = GT, pattern = "/", replacement = ":", vectorize_all = FALSE)) %>%
            mutate(
              GT = stri_replace_all_fixed(GT, pattern = ".:.", replacement = "NA", vectorize_all = FALSE),
              GT = replace(GT, which(GT == "NA"), NA)
            ) %>%
            group_by(INDIVIDUALS, POP_ID) %>% 
            tidyr::spread(data = ., key = MARKERS, value = GT) %>%
            ungroup() %>% 
            arrange(POP_ID, INDIVIDUALS)
        }
        
        if (impute == "allele") {
          input.prep <- input %>%
            select(-REF, -ALT) %>%
            mutate(GT = stri_replace_all_fixed(str = GT, pattern = "/", replacement = ":", vectorize_all = FALSE)) %>%
            mutate(
              GT = stri_replace_all_fixed(GT, pattern = ".:.", replacement = "NA:NA", vectorize_all = FALSE)
            ) %>% 
            tidyr::separate(col = GT, into = c("A1", "A2"), sep = ":") %>%  # separate the genotypes into alleles
            tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
            mutate(GT = replace(GT, which(GT == "NA"), NA)) %>%
            group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
            tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
            ungroup() %>% 
            arrange(POP_ID, INDIVIDUALS)
        }
        
      } # End VCF prep file for imputation
      
      if (data.type == "plink.file") {
        
        if (impute == "genotype"){
          input.prep <- gsi.prep %>% 
            group_by(MARKERS, INDIVIDUALS, POP_ID) %>% 
            tidyr::spread(data = ., key = ALLELES, value = GT) %>% 
            tidyr::unite(col = GT, A1, A2, sep = "", remove = TRUE) %>% 
            mutate(
              GT = stri_replace_all_fixed(GT, pattern = "000000", replacement = "NA", vectorize_all = FALSE),
              GT = replace(GT, which(GT == "NA"), NA)
            ) %>%
            group_by(INDIVIDUALS, POP_ID) %>% 
            tidyr::spread(data = ., key = MARKERS, value = GT) %>%
            ungroup() %>% 
            arrange(POP_ID, INDIVIDUALS)
        }
        
        if (impute == "allele"){
          input.prep <- gsi.prep %>%
            mutate(
              GT = stri_replace_all_fixed(GT, pattern = "000", replacement = "NA", vectorize_all = FALSE),
              GT = replace(GT, which(GT == "NA"), NA)
            ) %>%
            group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
            tidyr::spread(data = ., key = MARKERS, value = GT) %>%
            ungroup() %>% 
            arrange(POP_ID, INDIVIDUALS)
        }
        
        # glimpse(test)
      } # End PLINK prep file for imputation
      
      if (data.type == "df.file" | data.type == "haplo.file") { # for df input
        if (impute == "genotype") {
          input.prep <- input %>%
            mutate(
              GT = stri_replace_all_fixed(GT, pattern = "000000", replacement = "NA", vectorize_all = FALSE),
              GT = replace(GT, which(GT == "NA"), NA)
            ) %>%
            group_by(INDIVIDUALS, POP_ID) %>% 
            tidyr::spread(data = ., key = MARKERS, value = GT) %>%
            ungroup() %>% 
            arrange(POP_ID, INDIVIDUALS)
        }
        if (impute == "allele") {
          input.prep <- gsi.prep %>%
            mutate(
              GT = stri_replace_all_fixed(GT, pattern = "000", replacement = "NA", vectorize_all = FALSE),
              GT = replace(GT, which(GT == "NA"), NA)
            ) %>%
            group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
            tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
            ungroup() %>% 
            arrange(POP_ID, INDIVIDUALS)
        }
      } # End data frame prep file for imputation
      
      # keep stratification
      strata.df.impute <- input.prep %>% 
        select(INDIVIDUALS, POP_ID) %>% 
        distinct(INDIVIDUALS, POP_ID)
      
      # Imputation with Random Forest
      if (imputation.method == "rf") {
        # Parallel computations options
        options(rf.cores = parallel.core, mc.cores = parallel.core)
        
        # imputations using Random Forest with the package randomForestSRC
        impute_genotype_rf <- function(x) {
          randomForestSRC::impute.rfsrc(data = x,
                                        ntree = num.tree,
                                        nodesize = 1,
                                        nsplit = split.number,
                                        nimpute = iteration.rf,
                                        do.trace = verbose)
        } # End of imputation function
        
        # Random Forest by pop
        if (imputations.group == "populations") {
          message("Imputations computed by populations, take a break...")
          df.split.pop <- split(x = input.prep, f = input.prep$POP_ID) # slip data frame by population
          pop.list <- names(df.split.pop) # list the pop
          imputed.dataset <-list() # create empty list
          
          # Function to go through the populations
          impute_rf_pop <- function(pop.list, ...){
            sep.pop <- df.split.pop[[pop.list]]
            sep.pop <- suppressWarnings(
              plyr::colwise(factor, exclude = NA)(sep.pop)
            )
            # message of progress for imputations by population
            message(paste("Completed imputations for pop ", pop.list, sep = ""))
            # imputed.dataset[[i]] <- impute_markers_rf(sep.pop) # test with foreach
            imputed.dataset <- impute_genotype_rf(sep.pop)
            imputed.dataset <- suppressWarnings(
              plyr::colwise(as.character, exclude = NA)(imputed.dataset)
            )
            return(imputed.dataset)
          } # End impute_rf_pop
          
          input.imp <- list()
          input.imp <- parallel::mclapply(
            X = pop.list, 
            FUN = impute_rf_pop, 
            mc.preschedule = FALSE, 
            mc.silent = FALSE, 
            mc.cores = parallel.core
          )
          
          # Compiling the results
          message("Compiling imputations results")
          input.imp <- suppressWarnings(bind_rows(input.imp))
          
          # Second round of imputations (globally) to remove introduced NA 
          # In case that some pop don't have the markers
          input.imp <- suppressWarnings(plyr::colwise(factor, exclude = NA)(input.imp)) # Make the columns factor
          input.imp <- impute_genotype_rf(input.imp) # impute globally
          input.imp <- plyr::colwise(as.character, exclude = NA)(input.imp)
          
          # reintroduce the stratification
          input.imp <- suppressWarnings(
            left_join(strata.df.impute, 
                      input.imp %>% 
                        select(-POP_ID)
                      , by = "INDIVIDUALS") %>% 
              arrange(POP_ID, INDIVIDUALS) %>% 
              ungroup()
          )
          
          # dump unused objects
          df.split.pop <- NULL
          pop.list <- NULL
          sep.pop <- NULL
          imputed.dataset <- NULL
          input.prep <- NULL
          
        } # End imputation RF populations
        # Random Forest global
        if (imputations.group == "global") { # Globally (not by pop_id)
          message("Imputations computed globally, take a break...")
          input.prep <- plyr::colwise(factor, exclude = NA)(input.prep)
          input.imp <- impute_genotype_rf(input.prep)
          input.imp <- plyr::colwise(as.character, exclude = NA)(input.imp)
          
          # reintroduce the stratification
          input.imp <- suppressWarnings(
            left_join(strata.df.impute, 
                      input.imp %>% 
                        select(-POP_ID)
                      , by = "INDIVIDUALS") %>% 
              arrange(POP_ID, INDIVIDUALS) %>% 
              ungroup()
          )
          input.prep <- NULL # remove unused object
        } # End imputation RF global
        
        
        # data prep
        if (impute == "genotype") {
          input.imp <- suppressWarnings(
            input.imp %>%
              tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID))
          )
        }
        if (impute == "allele") {
          input.imp <- suppressWarnings(
            input.imp %>%
              tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES))
          )
        }
      } # End imputation RF
      
      # Imputation using the most common genotype
      if (imputation.method == "max") { # End imputation max
        if (imputations.group == "populations") {
          message("Imputations computed by populations")
          
          if (impute == "genotype"){
            input.imp <- suppressWarnings(
              input.prep %>%
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
                ungroup()
            )
          }
          
          if (impute == "allele"){
            input.imp <- suppressWarnings(
              input.prep %>%
                tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
                group_by(MARKERS, POP_ID) %>%
                mutate(
                  GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE)),
                  GT = replace(GT, which(GT == "NA"), NA)
                ) %>%
                # the next 2 steps are necessary to remove introduced NA if some pop don't have the markers
                # will take the global observed values by markers for those cases.
                group_by(MARKERS) %>%
                mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
                ungroup()
            )
          }
          
          input.prep <- NULL # remove unused object
          
        } # End imputation max populations 
        if (imputations.group == "global") {
          # Globally (not by pop_id)
          message("Imputations computed globally")
          if (impute == "genotype"){
            input.imp <- suppressWarnings(
              input.prep %>%
                tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
                group_by(MARKERS) %>%
                mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
                ungroup()
            )
          }
          
          if (impute == "allele"){
            input.imp <- suppressWarnings(
              input.prep %>%
                tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
                group_by(MARKERS) %>%
                mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
                ungroup()
            )
          }
          
          input.prep <- NULL # remove unused object
        } # End imputation max global 
      } # End imputations max
      
      # prepare the imputed dataset for gsi_sim or adegenet
      message("Preparing imputed data set...")
      if (assignment.analysis == "gsi_sim") {
        if (impute == "genotype") {
          if (data.type == "vcf.file") { # for VCF input
            gsi.prep.imp <- suppressWarnings(
              input.imp %>%
                # tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% # make tidy
                tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_", extra = "drop", remove = TRUE) %>%
                tidyr::gather(key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID))
            )
          }
          if (data.type == "plink.file" | data.type == "df.file" | data.type == "haplo.file") {
            gsi.prep.imp <- input.imp %>%
              # tidyr::gather(key = MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>% 
              tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3) %>%  # separate the genotypes into alleles
              tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID))
          }
        }
        if (impute == "allele") {
          gsi.prep.imp <- input.imp
        }
      }
      
      # adegenet
      if (assignment.analysis == "adegenet") {
        if (impute == "genotype") {
          if (data.type == "vcf.file") {
            genind.prep.imp <- input.imp %>% 
              # tidyr::gather(key = MARKERS, value = GT, -c(POP_ID, INDIVIDUALS, ALLELES)) %>% # make tidy
              tidyr::spread(data = ., key = ALLELES, value = GT) %>%
              tidyr::unite(GT, A1, A2, sep = ":", remove = TRUE) %>%
              mutate(GT = stri_replace_all_fixed(str = GT, pattern = c("0:0", "1:1", "0:1", "1:0"), replacement = c("2_0", "0_2", "1_1", "1_1"), vectorize_all = FALSE)) %>%
              arrange(MARKERS, POP_ID) %>%
              tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_", extra = "drop", remove = TRUE) %>%
              tidyr::gather(key = ALLELES, value = COUNT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% # make tidy
              tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>%
              group_by(POP_ID, INDIVIDUALS) %>%
              tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>%
              arrange(POP_ID, INDIVIDUALS)
          }
          if (data.type == "df.file" | data.type == "plink.file" | data.type == "haplo.file") {
            genind.prep.imp <- input.imp %>%
              # tidyr::gather(key = MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>% 
              tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3) %>%  # separate the genotypes into alleles
              tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID))
            
            genind.prep.imp <- suppressWarnings(
              genind.prep.imp %>%
                tidyr::spread(data = ., key = MARKERS, value = GT) %>% # this reintroduce the missing, but with NA
                ungroup() %>% 
                plyr::colwise(.fun = factor, exclude = NA)(.)
            )
            
            genind.prep.imp <- suppressWarnings(
              genind.prep.imp %>%
                ungroup() %>% 
                mutate_each(funs(as.integer), -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
                ungroup() %>% 
                tidyr::gather(data = ., key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
                select(-ALLELES) %>%
                group_by(POP_ID, INDIVIDUALS, MARKERS, GT) %>% 
                tally %>%
                ungroup() %>%
                tidyr::unite(MARKERS_ALLELES, MARKERS, GT, sep = ":", remove = TRUE) %>%
                arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>% 
                group_by(POP_ID, INDIVIDUALS) %>% 
                tidyr::spread(data = ., key = MARKERS_ALLELES, value = n) %>%
                ungroup() %>%
                tidyr::gather(data = ., key = MARKERS_ALLELES, value = COUNT, -c(INDIVIDUALS, POP_ID)) %>% 
                tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = ":", remove = TRUE) %>% 
                mutate(COUNT = as.numeric(stri_replace_na(str = COUNT, replacement = "0"))) %>% 
                ungroup() %>%
                arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES) %>% 
                tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>%
                tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>% 
                arrange(POP_ID, INDIVIDUALS)
            )
          }
        } # end impute genotype adegenet
        if (impute == "allele") {
          if (data.type == "vcf.file") {
            genind.prep.imp <- input.imp %>% 
              # tidyr::gather(key = MARKERS, value = GT, -c(POP_ID, INDIVIDUALS, ALLELES)) %>% # make tidy
              tidyr::spread(data = ., key = ALLELES, value = GT) %>%
              tidyr::unite(GT, A1, A2, sep = ":", remove = TRUE) %>%
              mutate(GT = stri_replace_all_fixed(str = GT, pattern = c("0:0", "1:1", "0:1", "1:0"), replacement = c("2_0", "0_2", "1_1", "1_1"), vectorize_all = FALSE)) %>%
              arrange(MARKERS, POP_ID) %>%
              tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_", extra = "drop", remove = TRUE) %>%
              tidyr::gather(key = ALLELES, value = COUNT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% # make tidy
              tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>%
              group_by(POP_ID, INDIVIDUALS) %>%
              tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>%
              arrange(POP_ID, INDIVIDUALS)
          }
          if (data.type == "df.file" | data.type == "plink.file" | data.type == "haplo.file") {
            genind.prep.imp <- suppressWarnings(
              input.imp %>%
                tidyr::spread(data = ., key = MARKERS, value = GT) %>% # this reintroduce the missing, but with NA
                ungroup() %>% 
                plyr::colwise(.fun = factor, exclude = NA)(.)
            )
            
            genind.prep.imp <- suppressWarnings(
              genind.prep.imp %>%
                ungroup() %>% 
                mutate_each(funs(as.integer), -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
                ungroup() %>% 
                tidyr::gather(data = ., key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
                select(-ALLELES) %>%
                group_by(POP_ID, INDIVIDUALS, MARKERS, GT) %>% 
                tally %>%
                ungroup() %>%
                tidyr::unite(MARKERS_ALLELES, MARKERS, GT, sep = ":", remove = TRUE) %>%
                arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>% 
                group_by(POP_ID, INDIVIDUALS) %>% 
                tidyr::spread(data = ., key = MARKERS_ALLELES, value = n) %>%
                ungroup() %>%
                tidyr::gather(data = ., key = MARKERS_ALLELES, value = COUNT, -c(INDIVIDUALS, POP_ID)) %>% 
                tidyr::separate(data = ., col = MARKERS_ALLELES, into = c("MARKERS", "ALLELES"), sep = ":", remove = TRUE) %>% 
                mutate(COUNT = as.numeric(stri_replace_na(str = COUNT, replacement = "0"))) %>% 
                ungroup() %>%
                arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES) %>% 
                tidyr::unite(MARKERS_ALLELES, MARKERS, ALLELES, sep = ".", remove = TRUE) %>%
                tidyr::spread(data = ., key = MARKERS_ALLELES, value = COUNT) %>% 
                arrange(POP_ID, INDIVIDUALS)
            )
          }
        } # end impute genotype adegenet
        
        # genind arguments common to all data.type
        genind.prep.imp <- genind.prep.imp %>% arrange(POP_ID, INDIVIDUALS)
        ind <- as.character(genind.prep.imp$INDIVIDUALS)
        pop <- genind.prep.imp$POP_ID
        genind.df <- genind.prep.imp %>%
          ungroup() %>% 
          select(-c(INDIVIDUALS, POP_ID))
        
        rownames(genind.df) <- ind
        loc.names <- colnames(genind.df)
        strata <- genind.prep.imp %>% 
          ungroup() %>% 
          select(INDIVIDUALS, POP_ID) %>% 
          distinct(INDIVIDUALS, POP_ID)
        
        # genind constructor
        prevcall <- match.call()
        genind.object.imp <- genind(tab = genind.df, pop = pop, prevcall = prevcall, ploidy = 2, type = "codom", strata = strata, hierarchy = NULL)
        
        # sum <- summary(genind.object.imp) # test
        # sum$NA.perc # test
        
        ind <- NULL
        pop <- NULL
        genind.df <- NULL
        # genind.prep <- NULL
        # genind.prep.imp <- NULL
        
      } # end adegenet
      
    } # End imputations
    
    # Sampling of markers ******************************************************
    # unique list of markers after all the filtering
    # if "all" is present in the list, change to the maximum number of markers
    unique.markers <- input %>% 
      select(MARKERS) %>% 
      distinct(MARKERS) %>% 
      arrange(MARKERS)
    
    
    marker.number <- stri_replace_all_fixed(str = marker.number, pattern = "all", 
                                            replacement = nrow(unique.markers), 
                                            vectorize_all = TRUE)
    marker.number <- as.numeric(marker.number)
    
    # In marker.number, remove marker group higher than the max number of markers
    removing.marker <- purrr::keep(.x = marker.number, .p = marker.number > nrow(unique.markers))
    
    if (length(removing.marker) > 0) {
      message.marker <- stri_c(removing.marker, collapse = ", ")
      message("Removing marker.number higher than the max number of markers: ", message.marker)
    }
    
    marker.number <- purrr::discard(.x = marker.number, .p = marker.number > nrow(unique.markers))
    
    
    # Functions ******************************************************************
    # Fst function: Weir & Cockerham 1984
    fst_WC84 <- function(data, holdout.samples, ...) {
      # data <- input # test
      # holdout.samples <- holdout$INDIVIDUALS # test
      
      # data.type <- data$GT[[1]] # VCF vs DF # test
      pop.number <- n_distinct(data$POP_ID)
      
      if (is.null(holdout.samples)) { # use all the individuals
        if (data.type == "vcf.file") { # for VCF input
          # if (stri_detect_fixed(data.type, "_")){ # for VCF file input # test
          data.genotyped <- data %>%
            filter(GT != "0_0")
        } else { # for df and haplotype files
          data.genotyped <- data %>%
            filter(GT != "000000") # Check for df and plink...
        }
      } else { # with holdout set
        if (data.type == "vcf.file") { # for VCF input
          # if (stri_detect_fixed(data.type, "_")){ # for VCF file input
          data.genotyped <- data %>%
            filter(GT != "0_0") %>% # remove missing genotypes
            filter(!INDIVIDUALS %in% holdout.samples) # remove supplementary individual before ranking markers with Fst
        } else { # for df and haplotype files
          data.genotyped <- data %>%
            filter(GT != "000000") %>% # remove missing genotypes
            filter(!INDIVIDUALS %in% holdout.samples) # remove supplementary individual before ranking markers with Fst
        }
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
      if (data.type == "vcf.file") { # for VCF input
        # if (stri_detect_fixed(data.type, "_")){ # for VCF file input
        freq.al.locus <- data.genotyped %>%
          mutate(A1 = stri_sub(GT, 1, 1), A2 = stri_sub(GT, 3, 3)) %>% 
          select(-GT) %>%
          tidyr::gather(key = ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, POP_ID, MARKERS))
        # tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # test takes longer
      } else { # for DF, haplo and plink file input
        freq.al.locus <- data.genotyped %>%
          tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3) %>%  # separate the genotypes into alleles
          tidyr::gather(key = ALLELES_GROUP, ALLELES, -c(INDIVIDUALS, POP_ID, MARKERS))
      }
      
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
      if (data.type == "vcf.file") { # for VCF input
        # if (stri_detect_fixed(data.type, "_")){ # for VCF file input
        data.genotyped.het <- data.genotyped %>%
          mutate(het = ifelse(stri_sub(GT, 1, 1) != stri_sub(GT, 3, 3), 1, 0))
      } else { # for DF file input
        data.genotyped.het <- data.genotyped %>%
          mutate(het = ifelse(stri_sub(GT, 1, 3) != stri_sub(GT, 4, 6), 1, 0))
      }
      
      fst.ranked <- suppressWarnings(
        data.genotyped.het %>%
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
          mutate(
            RANKING = seq(from = 1, to = n()),
            QUARTILE = ntile(FST,10)
          ) 
      )
      
      # select(MARKERS, FIS, FST)
      # remove unused objects
      data <- NULL
      data.genotyped <- NULL
      data.genotyped.het <- NULL
      n.pop.locus <- NULL
      ind.count.locus <- NULL
      ind.count.locus.pop <- NULL
      freq.al.locus <- NULL
      freq.al.locus.pop <- NULL
      freq.al.locus.global <- NULL
      mean.n.pop.corrected.per.locus <- NULL
      ncal <- NULL
      
      return(fst.ranked)
    } # End fst_WC84 function
    
    # Write the files
    write_gsi <- function (data, markers.names, filename, i, m, ...) {
      # data <- data.select
      # markers.names = markers.names
      # filename = filename
      # i = i
      # m = m
      # 
      data$POP_ID <- droplevels(x = data$POP_ID)
      n.individuals <- n_distinct(data$INDIVIDUALS)  # number of individuals
      pop <- data$POP_ID  # Create a vector with the population ordered by levels
      data <- suppressWarnings(data %>% select(-POP_ID))  # remove pop id
      gsi_sim.split <- split(data, pop)  # split gsi_sim by populations
      filename <- filename  # gsi_sim filename
      
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
      return(filename)
    } # End write_gsi function
    
    # Assignment with gsi_sim
      assignment_analysis <- function(data, select.markers, markers.names, missing.data, i, m, holdout, filename, ...) {
        # data <- gsi.prep #test
        # data <- gsi.prep.imp #test
        # data <- genind.prep #test
        # missing.data <- "no.imputation" #test
        data.select <- suppressWarnings(
          data %>%
            semi_join(select.markers, by = "MARKERS") %>%
            arrange(MARKERS) %>%  # make tidy
            tidyr::unite(col = MARKERS_ALLELES, MARKERS , ALLELES, sep = "_") %>%
            arrange(POP_ID, INDIVIDUALS, MARKERS_ALLELES) %>%
            tidyr::spread(data = ., key = MARKERS_ALLELES, value = GT) %>%
            arrange(POP_ID, INDIVIDUALS)
        )
        
        # Write gsi_sim input file to directory
        input.gsi <- write_gsi(data = data.select, markers.names = markers.names, filename = filename, i = i, m = m)
        
        # Run gsi_sim ------------------------------------------------------------
        input.gsi <- stri_join(directory.subsample, input.gsi)
        output.gsi <- stri_replace_all_fixed(input.gsi, pattern = "txt", replacement = "output.txt")
        setwd(directory.subsample)
        system(paste(gsi_sim_binary(), "-b", input.gsi, "--self-assign > ", output.gsi))
        
        # Option remove the input file from directory to save space
        if (keep.gsi.files == FALSE) file.remove(input.gsi)
        
        # Get Assignment results -------------------------------------------------
        # Keep track of the holdout individual
        if (sampling.method == "ranked") {
          if (thl == "all") {
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
            tidyr::separate(OTHERS, c("SECOND_BEST_SCORE", "OTHERS"), sep = ";;", convert = TRUE, numerals = "no.loss")
        )
        
        if (!is.null(pop.id.start)){
          assignment <- suppressWarnings(
            assignment %>%
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
                METHOD = rep(sampling.method, n()),
                MISSING_DATA = rep(missing.data, n())
              ) %>%
              select(INDIVIDUALS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
              arrange(CURRENT)
          )
        } else {# with strata.df info
          
          assignment <- suppressWarnings(
            assignment %>%
              mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>% 
              left_join(strata.df, by = "INDIVIDUALS") %>%
              rename(CURRENT = POP_ID) %>% 
              mutate(
                INFERRED = factor(INFERRED, levels = unique(pop.labels), ordered = TRUE),
                INFERRED = droplevels(INFERRED),
                SECOND_BEST_POP = factor(SECOND_BEST_POP, levels = unique(pop.labels), ordered = TRUE),
                SECOND_BEST_POP = droplevels(SECOND_BEST_POP),
                SCORE = round(SCORE, 2),
                SECOND_BEST_SCORE = round(SECOND_BEST_SCORE, 2),
                MARKER_NUMBER = as.numeric(rep(n.locus, n())),
                METHOD = rep(sampling.method, n()),
                MISSING_DATA = rep(missing.data, n())
              ) %>%
              select(INDIVIDUALS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
              arrange(CURRENT)
          )
        }
        if (sampling.method == "random") {
          assignment <- assignment %>% 
            mutate(ITERATIONS = rep(i, n())) %>% 
            select(INDIVIDUALS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA, ITERATIONS) %>%
            arrange(CURRENT)
        }
        if (sampling.method == "ranked") {
          if (thl == "all") {
            assignment <- assignment
          } else {
            assignment <- filter(.data = assignment, INDIVIDUALS %in% holdout.id)
          }
        }
        
        # Option remove the output file from directory to save space
        if (keep.gsi.files == FALSE) file.remove(output.gsi)
        
        # saving preliminary results
        if (sampling.method == "ranked") {
          
          assignment <- assignment %>%
            mutate(METHOD = rep(stri_join("ranked_thl_", thl) , n()))
          
          if (thl != 1 & thl != "all") {
            assignment <- assignment %>%
              mutate(
                CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
                CURRENT = droplevels(CURRENT),
                ITERATIONS = rep(i, n())
              )
          }
        }
        return(assignment)
      } # End assignment_analysis function
    
    # Assignment with adegenet
      assignment_analysis_adegenet <- function(data, select.markers, markers.names, missing.data, i, m, holdout, ...) {
        # data <- genind.object #test
        # missing.data <- "no.imputation" #test
        data.select <- data[loc = select.markers$MARKERS]
        
        # Run adegenet *********************************************************
        pop.data <- data.select@pop
        pop.data <- droplevels(pop.data)
        
        
        if (sampling.method == "random") {
          # Alpha-Score DAPC
          # When all the individuals are accounted for in the model construction
          dapc.best.optim.a.score <- optim.a.score(dapc(data.select, n.da = length(levels(pop.data)), n.pca = round((length(indNames(data.select))/3)-1, 0)), pop = pop.data, plot = FALSE)$best
          message(stri_paste("a-score optimisation for iteration:", i, sep = " ")) # message not working in parallel...
          
          # DAPC with all the data
          dapc.all <- dapc(data.select, n.da = length(levels(pop.data)), n.pca = dapc.best.optim.a.score, pop = pop.data)
          message(stri_paste("DAPC iteration:", i, sep = " "))
          message(stri_paste("DAPC marker group:", m, sep = " "))
        }
        
        if (sampling.method == "ranked") {
          
          # Alpha-Score DAPC training data
          training.data <- data.select[!indNames(data.select) %in% holdout$INDIVIDUALS] # training dataset
          pop.training <- training.data@pop
          pop.training <- droplevels(pop.training)
          
          dapc.best.optim.a.score <- optim.a.score(dapc(training.data, n.da = length(levels(pop.training)), n.pca = round(((length(indNames(training.data))/3)-1), 0)), pop = pop.training, plot = FALSE)$best
          message(stri_paste("a-score optimisation for iteration:", i, sep = " "))
          
          dapc.training <- dapc(training.data, n.da = length(levels(pop.training)), n.pca = dapc.best.optim.a.score, pop = pop.training)
          message(stri_paste("DAPC of training data set for iteration:", i, sep = " "))
          
          # DAPC holdout individuals
          holdout.data <- data.select[indNames(data.select) %in% holdout$INDIVIDUALS] # holdout dataset
          pop.holdout <- holdout.data@pop
          pop.holdout <- droplevels(pop.holdout)
          assignment.levels <- levels(pop.holdout) # for figure
          rev.assignment.levels <- rev(assignment.levels)  # for figure 
          
          dapc.predict.holdout <- predict.dapc(dapc.training, newdata = holdout.data)
          message(stri_paste("Assigning holdout data for iteration:", i, sep = " "))
        }
        
        
        # Get Assignment results -------------------------------------------------
        # Keep track of the holdout individual
        if (sampling.method == "ranked") {
          if (thl == "all") {
            holdout.id <- NULL
          } else {
            holdout.id <- holdout$INDIVIDUALS
          }
        }
        
        # Number of markers
        n.locus <- m
        
        if (sampling.method == "random") {
          assignment <- data_frame(ASSIGNMENT_PERC = summary(dapc.all)$assign.per.pop*100) %>% 
            bind_cols(data_frame(POP_ID = levels(pop.data))) %>%
            mutate(ASSIGNMENT_PERC = round(ASSIGNMENT_PERC, 2)) %>% 
            select(POP_ID, ASSIGNMENT_PERC)
        }        
        if (sampling.method == "ranked") {
          assignment <- data.frame(INDIVIDUALS = indNames(holdout.data), POP_ID = pop.holdout, ASSIGN = dapc.predict.holdout$assign, dapc.predict.holdout$posterior) %>% 
            rename(CURRENT = POP_ID, INFERRED = ASSIGN) %>%
            mutate(
              CURRENT = factor(CURRENT, levels = rev.assignment.levels, ordered = TRUE),
              INFERRED = factor(INFERRED, levels = assignment.levels, ordered = TRUE)
            )
        }
        
        assignment <- assignment %>% 
          mutate(
            METHOD = rep(sampling.method, n()),
            ITERATIONS = rep(i, n()),
            MARKER_NUMBER = as.numeric(rep(n.locus, n())),
            MISSING_DATA = rep(missing.data, n())
          )
        
        return(assignment)
      } # End assignment_analysis_adegenet function

    # Random method ************************************************************
    if (sampling.method == "random") {
      message("Conducting Assignment analysis with markers selected randomly")
      # Number of times to repeat the sampling of markers
      iterations.list <- 1:iteration.method
      # iterations.list <- 1:200 # test
      # Function: Random selection of marker function + iteration.method
      marker_selection <- function(iteration.method) {
        m <- as.numeric(m)
        select.markers <- sample_n(tbl = unique.markers, size = m, replace = FALSE) %>%
          arrange(MARKERS) %>%
          mutate(
            ITERATIONS = rep(iteration.method, n()),
            MARKER_NUMBER = rep(m, n())
          )
      }
      markers.random.lists <- list()
      
      message("Making a list containing all the markers combinations")
      # Go through the function with the marker number selected
      for (m in marker.number) {
        res <- purrr::map(.x = iterations.list, .f = marker_selection)
        markers.random.lists[[m]] <- res
      }
      markers.random.lists <- purrr::flatten(markers.random.lists)
      # test <- markers.random.selection.list[[101]]
      
      markers.random.lists.table <- as_data_frame(bind_rows(markers.random.lists))
      write_tsv(x = markers.random.lists.table, path = paste0(directory.subsample, "markers.random.tsv"), col_names = TRUE, append = FALSE)
      
      # Set seed for random sampling
      random.seed <- sample(x = 1:1000000, size = 1)
      # set.seed(random.seed)
      # parallel::clusterSetRNGStream(cl = cl, iseed = random.seed)
      random.seed <- data.frame(RANDOM_SEED_NUMBER = random.seed)
      write_tsv(x = random.seed, path = paste0(directory.subsample, "random_seed_assignment_ngs.tsv"), col_names = TRUE, append = FALSE)
      
      message("Starting parallel computations for the assignment analysis
First sign of progress may take some time
Progress can be monitored with activity in the folder...")
      mrl <- NULL
      holdout <- NULL
      assignment.random <- list()
      assignment_random <- function(markers.random.lists, ...) {
        mrl <- markers.random.lists
        # mrl <- markers.random.lists[1] # test
        mrl <- data.frame(mrl)                      # marker random list
        i <- as.numeric(unique(mrl$ITERATIONS))     # iteration
        m <- as.numeric(unique(mrl$MARKER_NUMBER))  # number of marker selected
        
        select.markers <- mrl %>%                   # markers
          ungroup() %>% 
          select(MARKERS) %>% 
          arrange(MARKERS)
        
        # get the list of loci after filter
        markers.names <- unique(select.markers$MARKERS)
        
        # Assignment analysis without imputations
        filename <- stri_replace_all_fixed(gsi_sim.filename,
                                           pattern = "txt",
                                           replacement = stri_join(
                                             i, m, 
                                             "no.imputation", "txt", sep = "."
                                           )
        )
        
        if (assignment.analysis == "gsi_sim") {
          assignment.no.imp <- assignment_analysis(data = gsi.prep,
                                                   select.markers = select.markers,
                                                   markers.names = markers.names,
                                                   missing.data = "no.imputation", 
                                                   i = i, 
                                                   m = m,
                                                   holdout = NULL,
                                                   filename = filename
          )
        }
        if (assignment.analysis == "adegenet") {
          assignment.no.imp <- assignment_analysis_adegenet(data = genind.object,
                                                            select.markers = select.markers,
                                                            markers.names = markers.names,
                                                            missing.data = "no.imputation", 
                                                            i = i, 
                                                            m = m,
                                                            holdout = NULL
          )
        }
        
        # unused objects
        filename <- NULL
        
        # With imputations
        if (imputation.method != FALSE) {# with imputations
          if (imputation.method == "rf") {
            if (imputations.group == "populations") {
              missing.data <- "imputed RF populations"
            } else {
              missing.data <- "imputed RF global"
            }
          } else {
            if (imputations.group == "populations") {
              missing.data <- "imputed max populations"
            } else {
              missing.data <- "imputed max global"
            }
          }
          # Assignment analysis WITH imputations
          filename <- stri_replace_all_fixed(gsi_sim.filename,
                                             pattern = "txt",
                                             replacement = stri_join(
                                               i, m, 
                                               "imputed", "txt", sep = "."
                                             )
          )
          if (assignment.analysis == "gsi_sim") {
            assignment.imp <- assignment_analysis(data = gsi.prep.imp,
                                                  select.markers = select.markers,
                                                  markers.names = markers.names,
                                                  missing.data = missing.data, 
                                                  i = i,
                                                  m = m,
                                                  holdout = holdout,
                                                  filename = filename
            )
          }
          if (assignment.analysis == "adegenet") {
            assignment.imp <- assignment_analysis_adegenet(data = genind.object.imp,
                                                           select.markers = select.markers,
                                                           markers.names = markers.names,
                                                           missing.data = missing.data, 
                                                           i = i,
                                                           m = m,
                                                           holdout = holdout
            )
          }
          
          # unused objects
          select.markers <- NULL
          markers.names <- NULL
        } # End with imputations
        
        #compile assignment results each marker number for the iteration
        if (imputation.method == FALSE) {
          assignment <- assignment.no.imp
          gsi.prep.imp <- NULL
          input.imp <- NULL
        } else {
          assignment <- bind_rows(assignment.no.imp, assignment.imp)
        }
        # assignment <- mutate(.data = assignment, ITERATIONS = rep(i, n()))
        return(assignment)
      } # End of iterations for both with and without imputations
      
      assignment.res <- NULL
      assignment.res <- parallel::mclapply(
        X = markers.random.lists, 
        FUN = assignment_random, 
        mc.preschedule = FALSE, 
        mc.silent = FALSE, 
        mc.cores = parallel.core
      )
      
      # Compiling the results
      message("Compiling results")
      if(assignment.analysis == "adegenet") {
        assignment.res <- suppressWarnings(
          bind_rows(assignment.res) %>% 
            rename(CURRENT = POP_ID) %>% 
            mutate(
              SUBSAMPLE = rep(subsample.id, n()),
              CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
              CURRENT = droplevels(CURRENT)
            ) %>% 
            arrange(CURRENT, MARKER_NUMBER, MISSING_DATA, ITERATIONS)
        )
      } else {
        assignment.res <- suppressWarnings(
          bind_rows(assignment.res) %>% 
            mutate(SUBSAMPLE = rep(subsample.id, n())) %>% 
            # arrange(POP_ID, INDIVIDUALS, MARKER_NUMBER, MISSING_DATA, ITERATIONS)
            arrange(CURRENT, INDIVIDUALS, MARKER_NUMBER, MISSING_DATA, ITERATIONS)
        )
      }
      
      # Write the tables to directory
      # assignment results
      if(assignment.analysis == "gsi_sim") {
        if (is.null(subsample)) {
          if (imputation.method == FALSE) {
            filename.assignment.res <- stri_join("assignment", sampling.method, "no.imputation", "results", "individuals", "iterations", "tsv", sep = ".")
          } else { # with imputations
            filename.assignment.res <- stri_join("assignment", sampling.method, "imputed", "results", "individuals", "iterations", "tsv", sep = ".")
          }
        } else {# with subsampling
          if (imputation.method == FALSE) {
            filename.assignment.res <- stri_join("assignment", sampling.method, "no.imputation", "results", "individuals","iterations", "subsample", subsample.id, "tsv", sep = ".")
          } else { # with imputations
            filename.assignment.res <- stri_join("assignment", sampling.method, "imputed", "results", "individuals", "iterations", "subsample", subsample.id, "tsv", sep = ".")
          }
        }
        write_tsv(x = assignment.res, path = paste0(directory.subsample,filename.assignment.res), col_names = TRUE, append = FALSE)
      } else { # with adegenet
        if (is.null(subsample)) {
          if (imputation.method == FALSE) {
            filename.assignment.res <- stri_join("assignment", sampling.method, "no.imputation", "results", "iterations", "tsv", sep = ".")
          } else { # with imputations
            filename.assignment.res <- stri_join("assignment", sampling.method, "imputed", "results", "iterations", "tsv", sep = ".")
          }
        } else {# with subsampling
          if (imputation.method == FALSE) {
            filename.assignment.res <- stri_join("assignment", sampling.method, "no.imputation", "results", "iterations", "subsample", subsample.id, "tsv", sep = ".")
          } else { # with imputations
            filename.assignment.res <- stri_join("assignment", sampling.method, "imputed", "results", "iterations", "subsample", subsample.id, "tsv", sep = ".")
          }
        }
        write_tsv(x = assignment.res, path = paste0(directory.subsample,filename.assignment.res), col_names = TRUE, append = FALSE)
      }
      
      if (assignment.analysis == "gsi_sim") {
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
      }
      
      if (assignment.analysis == "adegenet") {
        assignment.stats.pop <- suppressWarnings(
          assignment.res %>%
            group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
            summarise(
              MEAN = round(mean(ASSIGNMENT_PERC), 2),
              SE = round(sqrt(var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
              MIN = round(min(ASSIGNMENT_PERC), 2),
              MAX = round(max(ASSIGNMENT_PERC), 2),
              MEDIAN = round(median(ASSIGNMENT_PERC), 2),
              QUANTILE25 = round(quantile(ASSIGNMENT_PERC, 0.25), 2),
              QUANTILE75 = round(quantile(ASSIGNMENT_PERC, 0.75), 2)
            ) %>%
            ungroup %>% 
            mutate(
              CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
              CURRENT = droplevels(CURRENT)
            ) %>%
            arrange(CURRENT, MARKER_NUMBER)
        )
      }
      
      # Next step is common for gsi_sim and adegenet
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
            ITERATIONS = rep(iteration.method, n())
          ) %>%
          select(CURRENT, MARKER_NUMBER, MEAN, MEDIAN, SE, MIN, MAX, QUANTILE25, QUANTILE75, SE_MIN, SE_MAX, METHOD, MISSING_DATA, ITERATIONS)
      )
      
      
      # Write the tables to directory
      # assignment summary stats
      if (is.null(subsample)) {
        if (imputation.method == FALSE) {
          filename.assignment.sum <- stri_join("assignment", sampling.method, "no.imputation", "results", "summary.stats", "tsv", sep = ".")
        } else { # with imputations
          filename.assignment.sum <- stri_join("assignment", sampling.method, "imputed", "results", "summary.stats", "tsv", sep = ".")
        }
      } else {# with subsampling
        if (imputation.method == FALSE) {
          filename.assignment.sum <- stri_join("assignment", sampling.method, "no.imputation", "results", "summary.stats", "subsample", subsample.id, "tsv", sep = ".")
        } else { # with imputations
          filename.assignment.sum <- stri_join("assignment", sampling.method, "imputed", "results", "summary.stats", "subsample", subsample.id, "tsv", sep = ".")
        }
      }
      write_tsv(x = assignment.summary.stats, path = paste0(directory.subsample,filename.assignment.sum), col_names = TRUE, append = FALSE)
      
    } # End method random
    
    # Ranked method ************************************************************
    if (sampling.method == "ranked") {
      message("Conducting Assignment analysis with ranked markers")
      
      # List of all individuals
      ind.pop.df<- input %>% 
        ungroup %>% 
        select(POP_ID, INDIVIDUALS) %>% 
        distinct(POP_ID, INDIVIDUALS)
      
      # thl selection
      message("Using thl method, ranking Fst with training samples...")
      if (thl == 1) {
        # Will go through the individuals in the list one by one.
        iterations.list <- ind.pop.df$INDIVIDUALS
        # Keep track of holdout individuals
        holdout.individuals <- ind.pop.df %>%
          mutate(ITERATIONS = stri_join("HOLDOUT", seq(1:n()), sep = "_"))
      } else if (thl == "all") { # no holdout for that one
        iterations.list <- iteration.method
        holdout.individuals <- NULL
        message("Warning: using all the individuals for ranking markers based on Fst\nNo holdout samples")
        message("Recommended reading: \nAnderson, E. C. (2010) Assessing the power of informative subsets of
                loci for population assignment: standard methods are upwardly biased.\nMolecular ecology resources 10, 4:701-710.")
      } else {
        # Create x (iterations) list of y (thl) proportion of individuals per pop.
        if (stri_detect_fixed(thl, ".") & thl < 1) {
          # iteration.method <- 5 # test
          # thl <- 0.4 # test
          holdout.individuals.list <- list()
          iterations.list <- 1:iteration.method
          for (x in 1:iteration.method) {
            holdout.individuals <- ind.pop.df %>%
              group_by(POP_ID) %>%
              sample_frac(thl, replace = FALSE) %>%  # sampling fraction for each pop
              arrange(POP_ID, INDIVIDUALS) %>%
              ungroup() %>%
              select(INDIVIDUALS) %>%
              mutate(ITERATIONS = rep(x, n()))
            holdout.individuals.list[[x]] <- holdout.individuals
          }
          holdout.individuals <- as.data.frame(bind_rows(holdout.individuals.list))
        }
        
        # Create x (iterations) list of y (thl) individuals per pop.
        if (thl > 1) {
          holdout.individuals.list <- list()
          iterations.list <- 1:iteration.method
          for (x in 1:iteration.method) {
            holdout.individuals <- ind.pop.df %>%
              group_by(POP_ID) %>%
              sample_n(thl, replace = FALSE) %>% # sampling individuals for each pop
              arrange(POP_ID, INDIVIDUALS) %>%
              ungroup() %>%
              select(INDIVIDUALS) %>%
              mutate(ITERATIONS = rep(x, n()))
            holdout.individuals.list[[x]] <- holdout.individuals
          }
          holdout.individuals <- as.data.frame(bind_rows(holdout.individuals.list))
        }
      } # End tracking holdout individuals
      write_tsv(x = holdout.individuals, 
                path = paste0(directory.subsample,"holdout.individuals.tsv"), 
                col_names = TRUE, 
                append = FALSE
      )
      message("Holdout samples saved in your folder")
      
      # Going through the loop of holdout individuals
      message("Starting parallel computations for the assignment analysis
First sign of progress may take some time
Progress can be monitored with activity in the folder...")
      
      assignment_ranking <- function(iterations.list, ...) {
        # i <- "TRI_09" #test
        # i <- "CAR_01" #test
        # i <- 1
        # i <- 5
        i <- iterations.list
        
        # Ranking Fst with training dataset (keep holdout individuals out)
        message("Ranking markers based on Fst")
        if (thl == "all") {
          holdout <- NULL
          fst.ranked <- fst_WC84(data = input, holdout.samples = NULL)
          if (imputation.method != FALSE) {
            fst.ranked.imp <- fst_WC84(data = input.imp, holdout.samples = NULL)
          }
        } else if (thl == 1) {
          holdout <- data.frame(INDIVIDUALS = i)
          fst.ranked <- fst_WC84(data = input, holdout.samples = holdout$INDIVIDUALS)
          if (imputation.method != FALSE) {
            fst.ranked.imp <- fst_WC84(data = input.imp, holdout.samples = holdout$INDIVIDUALS)
          }
        } else { # thl proportion or > 1
          holdout <- data.frame(holdout.individuals.list[i])
          fst.ranked <- fst_WC84(data = input, holdout.samples = holdout$INDIVIDUALS)
          if (imputation.method != FALSE) {
            fst.ranked.imp <- fst_WC84(data = input.imp, holdout.samples = holdout$INDIVIDUALS)
          }
        }
        
        fst.ranked.filename <- stri_join("fst.ranked_", i, ".tsv", sep = "") # No imputation
        write_tsv(
          x = fst.ranked, 
          path = paste0(directory.subsample, fst.ranked.filename), 
          col_names = TRUE, 
          append = FALSE
        )
        
        if (imputation.method != FALSE) {  # With imputations
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
        
        # Markers numbers loop function
        message("Going throught the marker.number")
        assignment.marker <- list() # Create empty lists to feed the results
        assignment_marker_loop <- function(m, ...) {
          message("Marker number: ", m)
          # m <- 200 # test
          # m <- 400 # test
          m <- as.numeric(m)
          RANKING <- NULL
          select.markers <- filter(.data = fst.ranked, RANKING <= m) %>%
            select(MARKERS)
          
          # get the list of markers after filter
          markers.names <- unique(select.markers$MARKERS)
          
          # Assignment analysis without imputations
          filename <- stri_replace_all_fixed(gsi_sim.filename,
                                             pattern = "txt",
                                             replacement = stri_join(
                                               i, m, 
                                               "no.imputation", "txt", sep = "."
                                             )
          )
          if (assignment.analysis == "gsi_sim") {
            assignment.no.imp <- assignment_analysis(data = gsi.prep,
                                                     select.markers = select.markers,
                                                     markers.names = markers.names,
                                                     missing.data = "no.imputation", 
                                                     i = i, 
                                                     m = m,
                                                     holdout = holdout,
                                                     filename = filename
            )
          }
          
          if (assignment.analysis == "adegenet") {
            assignment.no.imp <- assignment_analysis_adegenet(data = genind.object,
                                                              select.markers = select.markers,
                                                              markers.names = markers.names,
                                                              missing.data = "no.imputation", 
                                                              i = i, 
                                                              m = m,
                                                              holdout = holdout
            )
          }
          # unused objects
          select.markers <- NULL
          markers.names <- NULL
          RANKING <- NULL
          filename <- NULL
          
          # With imputations
          if (imputation.method != FALSE) {  # with imputations
            
            select.markers <- filter(.data = fst.ranked.imp, RANKING <= m) %>%
              select(MARKERS)
            
            # get the list of markers after filter
            markers.names <- unique(select.markers$MARKERS)  # not the same in no imputation
            
            if (imputation.method == "rf") {
              if (imputations.group == "populations") {
                missing.data <- "imputed RF populations"
              } else {
                missing.data <- "imputed RF global"
              }
            } else {
              if (imputations.group == "populations") {
                missing.data <- "imputed max populations"
              } else {
                missing.data <- "imputed max global"
              }
            }
            
            # Assignment analysis WITH imputations
            filename <- stri_replace_all_fixed(gsi_sim.filename,
                                               pattern = "txt",
                                               replacement = stri_join(
                                                 i, m, 
                                                 "imputed", "txt", sep = "."
                                               )
            )
            if (assignment.analysis == "gsi_sim") {
              assignment.imp <- assignment_analysis(data = gsi.prep.imp,
                                                    select.markers = select.markers,
                                                    markers.names = markers.names,
                                                    missing.data = missing.data, 
                                                    i = i,
                                                    m = m,
                                                    holdout = holdout,
                                                    filename = filename
              )
            }
            if (assignment.analysis == "adegenet") {
              assignment.imp <- assignment_analysis_adegenet(data = gsi.prep.imp,
                                                             select.markers = select.markers,
                                                             markers.names = markers.names,
                                                             missing.data = missing.data, 
                                                             i = i,
                                                             m = m,
                                                             holdout = holdout
              )
            }
            
            # unused objects
            select.markers <- NULL
            markers.names <- NULL
            RANKING <- NULL
          } # End with imputations
          
          #compile assignment results each marker number for the iteration
          if (imputation.method == FALSE) {# with imputations
            assignment <- assignment.no.imp
            fst.ranked.imp <- NULL
            gsi.prep.imp <- NULL
            input.imp <- NULL
          } else {
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
          fst.ranked.imp = fst.ranked.imp,
          i = i,
          holdout = holdout,
          input = input,
          gsi.prep = gsi.prep,
          input.imp = input.imp,
          # vcf = vcf.imp, # was an error before, double check...
          gsi.prep.imp = gsi.prep.imp,
          pop.levels = pop.levels,
          pop.labels = pop.labels,
          pop.id.start =  pop.id.start,
          pop.id.end = pop.id.end,
          sampling.method = sampling.method,
          thl = thl,
          iteration.method = iteration.method,
          gsi_sim.filename = gsi_sim.filename,
          keep.gsi.files = keep.gsi.files,
          imputation.method = imputation.method,
          parallel.core = parallel.core
        )
        
        message("Summarizing the assignment analysis results by iterations and marker group")
        assignment.res.summary <- suppressWarnings(
          bind_rows(purrr::flatten(assignment.marker))
        )
        
        res.filename <- stri_join("assignment_", i, ".tsv", sep = "") # No imputation
        write_tsv(x = assignment.res.summary, path = paste0(directory.subsample, res.filename), 
                  col_names = TRUE, 
                  append = FALSE
        )
        return(assignment.res.summary)
      }  # End assignment ranking function
      
      # using mclapply
      assignment.res <- list()
      assignment.res <- parallel::mclapply(
        X = iterations.list, 
        FUN = assignment_ranking, 
        mc.preschedule = TRUE, 
        mc.silent = TRUE, 
        mc.cores = parallel.core,
        marker.number = marker.number
      )
      
      # Compiling the results
      message("Compiling results")
      assignment.res.summary <- suppressWarnings(
        bind_rows(assignment.res)%>% 
          mutate(SUBSAMPLE = rep(subsample.id, n())) %>% 
          # arrange(CURRENT, INDIVIDUALS, MARKER_NUMBER, MISSING_DATA, ITERATIONS)
          arrange(CURRENT, INDIVIDUALS, MARKER_NUMBER, MISSING_DATA)
      )
      
      # assignment results
      if (is.null(subsample)) {
        if (imputation.method == FALSE) {
          filename.assignment.res <- stri_join("assignment", sampling.method, "no.imputation", "results", "individuals", "iterations", "tsv", sep = ".")
        } else { # with imputations
          filename.assignment.res <- stri_join("assignment", sampling.method, "imputed", "results", "individuals", "iterations", "tsv", sep = ".")
        }
      } else {# with subsampling
        if (imputation.method == FALSE) {
          filename.assignment.res <- stri_join("assignment", sampling.method, "no.imputation", "results", "individuals","iterations", "subsample", subsample.id, "tsv", sep = ".")
        } else { # with imputations
          filename.assignment.res <- stri_join("assignment", sampling.method, "imputed", "results", "individuals", "iterations", "subsample", subsample.id, "tsv", sep = ".")
        }
      }
      write_tsv(x = assignment.res.summary, path = paste0(directory.subsample,filename.assignment.res), col_names = TRUE, append = FALSE)
      
      
      if (thl == 1 | thl == "all") {
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
        # thl != 1 or "all"
        # summary stats
        if (assignment.analysis == "adegenet") {
          assignment.res.summary.prep <- assignment.res.summary %>% 
            # assignment.res.summary <- assignment.res.summary %>% 
            group_by(CURRENT, INFERRED, ITERATIONS, MARKER_NUMBER, MISSING_DATA, METHOD) %>% 
            tally %>% 
            group_by(CURRENT) %>% 
            mutate(TOTAL = sum(n)) %>% 
            ungroup() %>% 
            mutate(ASSIGNMENT_PERC = round(n/TOTAL*100, 0)) %>% 
            filter(CURRENT == INFERRED) %>% 
            select(-n, -TOTAL)
        }
        
        if (assignment.analysis == "gsi_sim") {
          assignment.res.summary.prep <- assignment.res.summary %>% 
          group_by(CURRENT, MARKER_NUMBER, METHOD, MISSING_DATA, ITERATIONS) %>%
          summarise(
            n = length(CURRENT[as.character(CURRENT) == as.character(INFERRED)]),
            TOTAL = length(CURRENT)
          ) %>%
          ungroup() %>% 
          mutate(ASSIGNMENT_PERC = round(n/TOTAL*100, 0)) %>% 
          select(-n, -TOTAL)
        }
        
        if (is.null(subsample)) {
          if (imputation.method == FALSE) {
            filename.assignment.res.sum <- stri_join("assignment", sampling.method, "no.imputation", "results", "summary", "tsv", sep = ".")
          } else { # with imputations
            filename.assignment.res.sum <- stri_join("assignment", sampling.method, "imputed", "results", "summary", "tsv", sep = ".")
          }
        } else {# with subsampling
          if (imputation.method == FALSE) {
            filename.assignment.res.sum <- stri_join("assignment", sampling.method, "no.imputation", "results", "summary", "subsample", subsample.id, "tsv", sep = ".")
          } else { # with imputations
            filename.assignment.res.sum <- stri_join("assignment", sampling.method, "imputed", "results", "summary", "subsample", subsample.id, "tsv", sep = ".")
          }
        }
        write_tsv(x = assignment.res.summary.prep, path = paste0(directory.subsample,filename.assignment.res.sum), col_names = TRUE, append = FALSE)
        
        assignment.stats.pop <- assignment.res.summary.prep %>%
          # assignment.stats.pop <- assignment.res.summary %>%
          mutate(
            CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered = TRUE),
            CURRENT = droplevels(CURRENT)
          ) %>%
          group_by(CURRENT, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
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
          group_by(MARKER_NUMBER, METHOD, MISSING_DATA) %>%
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
      } # End thl != 1
      
      # Write the tables to directory
      # assignment summary stats
      if (is.null(subsample)) {
        if (imputation.method == FALSE) {
          filename.assignment.sum <- stri_join("assignment", sampling.method, "no.imputation", "results", "summary.stats", "tsv", sep = ".")
        } else { # with imputations
          filename.assignment.sum <- stri_join("assignment", sampling.method, "imputed", "results", "summary.stats", "tsv", sep = ".")
        }
      } else {# with subsampling
        if (imputation.method == FALSE) {
          filename.assignment.sum <- stri_join("assignment", sampling.method, "no.imputation", "results", "summary.stats", "subsample", subsample.id, "tsv", sep = ".")
        } else { # with imputations
          filename.assignment.sum <- stri_join("assignment", sampling.method, "imputed", "results", "summary.stats", "subsample", subsample.id, "tsv", sep = ".")
        }
      }
      write_tsv(x = assignment.summary.stats, path = paste0(directory.subsample,filename.assignment.sum), col_names = TRUE, append = FALSE)
    } # End of ranked thl method
    
    # update the assignment with subsampling iterations id
    assignment.summary.stats <- assignment.summary.stats %>% 
      mutate(SUBSAMPLE = rep(subsample.id, n()))
    return(assignment.summary.stats)
  } # End assignment_function
  
  res <- map(.x = subsample.list, .f = assignment_function,
             assignment.analysis = assignment.analysis,
             input = input,
             snp.ld = snp.ld,
             common.markers = common.markers,
             maf.thresholds = maf.thresholds,
             maf.pop.num.threshold = maf.pop.num.threshold,
             maf.approach = maf.approach,
             maf.operator = maf.operator,
             marker.number = marker.number,
             pop.levels = pop.levels,
             pop.labels = pop.labels,
             pop.id.start =  pop.id.start,
             pop.id.end = pop.id.end,
             sampling.method = sampling.method,
             thl = thl,
             iteration.method = iteration.method,
             gsi_sim.filename = gsi_sim.filename,
             directory = directory,
             keep.gsi.files = keep.gsi.files,
             imputation.method = imputation.method,
             impute = impute,
             imputations.group = imputations.group,
             num.tree = num.tree,
             iteration.rf = iteration.rf,
             split.number = split.number,
             verbose = verbose,
             parallel.core = parallel.core
  )
  res <- bind_rows(res)
  
  if (imputation.method == FALSE) {
    filename.res <- stri_join("assignment", sampling.method, "no.imputation.results.summary.stats.subsample", "tsv", sep = ".")
  } else { # with imputations
    filename.res <- stri_join("assignment", sampling.method, "imputed.results,summary.stats.subsample", "tsv", sep = ".")
  }
  write_tsv(x = res, path = paste0(directory,filename.res), col_names = TRUE, append = FALSE)
  
  # Summary of the subsampling iterations
  if (iteration.subsample > 1) {
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
      mutate(SUBSAMPLE = factor(SUBSAMPLE, levels = c(1:iteration.subsample, "OVERALL"), ordered = TRUE)) %>% 
      arrange(CURRENT, MARKER_NUMBER, SUBSAMPLE)
    
    if (imputation.method == FALSE) {
      filename.assignment.sum.subsample <- stri_join("assignment", sampling.method, "no.imputation.results.summary.stats.subsample.overall", "tsv", sep = ".")
    } else { # with imputations
      filename.assignment.sum.subsample <- stri_join("assignment", sampling.method, "imputed.results.summary.stats.subsample.overall", "tsv", sep = ".")
    }
    write_tsv(x = res, path = paste0(directory, filename.assignment.sum.subsample), col_names = TRUE, append = FALSE)
    
    
    
    # unused objects
    res.pop.overall <- NULL
    res.overall <- NULL
    res.pop <- NULL
  } # End summary of the subsampling iterations
  
  # Assignment plot
  if (imputation.method == FALSE) { # no imputation
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
        axis.text.x = element_text(size = 8, family = "Helvetica", face = "bold", angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = element_text(size = 10, family = "Helvetica", face = "bold")
      )
  } else { #with imputations
    
    plot.assignment <- ggplot(res, aes(x = factor(MARKER_NUMBER), y = MEAN))+
      geom_point(aes(colour = MISSING_DATA), size = 2, alpha = 0.8) +
      geom_errorbar(aes(ymin = SE_MIN, ymax = SE_MAX), width = 0.3) +
      scale_colour_manual(name = "Missing data", values = c("gray33", "dodgerblue"))+
      scale_y_continuous(breaks = c(0, 10, 20 ,30, 40, 50, 60, 70, 80, 90, 100))+
      labs(x = "Marker number")+
      labs(y = "Assignment success (%)")+
      theme_bw()+
      theme(
        legend.position = "bottom",      
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"), 
        axis.title.x = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.x = element_text(size = 8, family = "Helvetica", face = "bold", angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.y = element_text(size = 10, family = "Helvetica", face = "bold"), 
        axis.text.y = element_text(size = 10, family = "Helvetica", face = "bold")
      )
  } # end plot
  
  # results
  res.list <- list(assignment = res, plot.assignment = plot.assignment)
  return(res.list)
} # End assignment_ngs
