# Conduct mixture analysis with gsi_sim

#' @name assignment_mixture
#' @title Mixture/Baseline assignment analysis of massive parallel sequencing data (GBS/RADseq, 
#' SNP chip, etc) using \code{gsi_sim} and \code{\link[adegenet]{adegenet}}

#' @description \code{gsi_sim} is a tool for doing and simulating genetic stock
#' identification and developed by Eric C. Anderson.
#' The arguments in the \code{assignment_mixture} function were tailored for the
#' reality of GBS/RADseq data to assign mixture samples to baseline populations
#' while maintaining a reproducible workflow.
#' Various input files are offered. Individuals, populations and
#' markers can be filtered and/or selected in several ways using blacklist,
#' whitelist and other arguments. Map-independent imputation of missing genotype
#' using Random Forest or the most frequent category is also available.
#' Markers can be randomly selected for a classic LOO (Leave-One-Out)
#' assignment or chosen based on ranked Fst. For this, the baseline samples are
#' used for the training and the mixture samples as holdout. 
#' Classic Leave-one-out is then used to assign individual mixture samples.


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
#' @param mixture A text file mixture individual ID. The column header is 
#' \code{INDIVIDUALS} and the file is in the working directory (e.g. "mixture.txt").
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

#' @param monomorphic.out (optional) Logical. For PLINK file, should the monomorphic 
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
#' will only keep markers in common (genotyped) between all the baseline samples (populations).


#' @param maf.thresholds (string, double, optional) String with 
#' local/populations and global/overall Minor Allele Frequency (maf) thresholds, respectively.
#' Default: \code{maf.thresholds = NULL}. The maf is calculated on the baseline samples only.
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
#' chosen based on ranked Fst \code{"ranked"} using the baseline samples 
#' for the training and the mixture samples as holdout. 
#' Classic Leave-one-out is then used to assign individual mixture samples.
#' @param iteration.method With random marker selection the iterations argument =
#' the number of iterations to repeat marker resampling. 
#' Default: \code{iteration.method = 10}.
#' With \code{marker.number = c(500, 1000)} and default iterations setting,
#' 500 markers will be randomly chosen 10 times and 1000 markers will be randomly
#' chosen 10 times.

#' @param folder (optional) The name of the folder created in the working directory to save the files/results.
#' @param filename (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. Default \code{assignment_data.txt}.
#' The number of markers used will be appended to the name of the file.
#' @param keep.gsi.files (Logical) Default \code{FALSE} The input and output gsi_sim files
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
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID column. With
#' a PLINK file, the first 2 columns of the \code{tfam} are used. 
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
#' @param impute.mixture (Logical) Imputations of mixture samples.
#' Default: \code{impute.mixture = FALSE}. For no imputation. 
#' For \code{impute.mixture = TRUE} the imputations.group (see below)
#' for the mixture samples is automatically set to 
#' \code{imputations.group = "global"}. Warning: bias could be introduced by imputing
#' missing genotype in the mixture samples.
#' @param impute (character) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
#' @param imputations.group \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by populations. If you choose
#' global, turn the verbose to \code{TRUE}, to see progress.
#' Default: \code{imputations.group = "populations"}.
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
#' \code{$assignment}.

#' @note \code{assignment_mixture} assumes that the command line version of gsi_sim 
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
#' @rdname assignment_mixture
#' @import parallel
#' @import stringi
#' @import adegenet
#' @import dplyr
#' @importFrom stats var median quantile
#' @importFrom purrr map
#' @importFrom purrr flatten
#' @importFrom purrr keep
#' @importFrom purrr discard
#' @importFrom data.table fread

#' @examples
#' \dontrun{
#' # with adegenet DAPC for the assignment and sampling.method = "random":
#' assignment.treefrog <- assignment_mixture(
#' data = "batch_1.vcf",
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
#' pop.levels = c("PAN", "COS")
#' pop.id.start = 5, pop.id.end = 7,
#' imputation.method = FALSE,
#' parallel.core = 12
#' )
#' # with gsi_sim for the mixture assignment and sampling.method = "ranked"
#' # Here I also want to impute the genotypes of the data (baseline and mixture) 
#' # using random forest:
#' assignment.tuna <- assignment_mixture(
#' data = "data.frame.genotypes.tuna.tsv",
#' mixture = "cohort.tuna.tsv"
#' assignment.analysis = "gsi_sim",
#' common.markers = TRUE,
#' marker.number = c(100, 200, 300),
#' sampling.method = "ranked",
#' subsample = 25,
#' iteration.subsample = 5
#' filename = "tuna.txt",
#' keep.gsi.files = FALSE,
#' pop.levels = c("BAJ", "IND"),
#' imputation.method = "rf", 
#' impute.mixture = TRUE, 
#' impute = "genotype", 
#' imputations.group = "populations", 
#' num.tree = 100, 
#' iteration.rf = 10, 
#' split.number = 100, 
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

#' @seealso \code{gsi_sim} development page is available here: \url{https://github.com/eriqande/gsi_sim}

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
      "COUNT", "MAX_COUNT_MARKERS", "hierarchy", "ANALYSIS", "NUMBER_ITERATIONS",
      "TOTAL_ITERATIONS", "MEAN_ITERATIONS", "X1", "X2", "NUMBER_SUBSAMPLE",
      "TOTAL_SUBSAMPLE", "MEAN_SUBSAMPLE"
    )
  )
}

assignment_mixture <- function(data,
                               mixture,
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
                               iteration.method = 10,
                               folder,
                               filename = "assignment_data.txt",
                               keep.gsi.files,
                               imputation.method = FALSE,
                               impute.mixture = FALSE,
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
  if (missing(mixture)) stop("mixture file missing")
  if (missing(assignment.analysis)) stop("assignment.analysis argument missing")
  if (assignment.analysis == "gsi_sim" & !gsi_sim_exists()){
    stop("Can't find the gsi_sim executable where it was expected at ", gsi_sim_binary_path(), ".  
          If you have internet access, you can install it
          from within R by invoking the function \"install_gsi_sim(fromSource = TRUE)\"")
  }
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
  if (missing(iteration.method)) iteration.method <- 10
  if (missing(filename)) filename <- "assignment_data.txt"
  if (missing(keep.gsi.files)) keep.gsi.files <- FALSE
  if (missing(imputation.method)) imputation.method <- FALSE
  if (missing(impute.mixture)) impute.mixture <- FALSE
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
  # Base filename
  base.filename <- filename # filename will change from time to time in the function
  
  # Create a folder based on filename to save the output files *****************
  if (is.null(folder)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    if (imputation.method == FALSE) {
      message("Map-imputation: no")
      directory <- stri_join(getwd(),"/", "assignment_mixture_analysis", "_no_imputations_", file.date, "/", sep = "")
      dir.create(file.path(directory))
    } else {
      message("Map-imputation: yes")
      directory <- stri_join(getwd(),"/","assignment_mixture_analysis", "_imputations_", imputation.method,"_", imputations.group, "_", file.date, "/", sep = "")
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
        CHROM = stri_replace_all_fixed(CHROM, pattern = "un", replacement = "1")
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
    strata.df$INDIVIDUALS <- stri_replace_all_fixed(str = strata.df$INDIVIDUALS, pattern = "_", replacement = "-", vectorize_all = TRUE)
    
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
      data.table = FALSE)
    
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
        rename(POP_ID = STRATA)
    }
    
    # Make tidy
    message("Tidying the PLINK file and integrating the tfam/strata file, for large dataset this may take several minutes...")
    input <- input %>% 
      tidyr::gather(key = INDIVIDUALS_ALLELES, value = GT, -LOCUS) %>%
      mutate(INDIVIDUALS = stri_replace_all_fixed(str = INDIVIDUALS_ALLELES, pattern = c("_A1", "_A2"), replacement = "", vectorize_all = FALSE)) %>% 
      left_join(strata.df, by = "INDIVIDUALS") %>% 
      mutate(
        POP_ID = factor(POP_ID, levels = pop.levels, ordered =TRUE),
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
                    distinct(LOCUS), 
                  by = "LOCUS")
      )
    }
    
    # Removing monomorphic markers
    if (monomorphic.out == TRUE) {
      message("Removing monomorphic markers...")
      mono.markers <- remove.missing.gt %>%
        group_by(LOCUS, GT) %>%
        distinct(LOCUS, GT) %>% 
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
      colClasses = "character",
      verbose = FALSE,
      showProgress = TRUE,
      data.table = FALSE
    ) %>% 
      as_data_frame() %>% 
      tidyr::gather(key = LOCUS, value = GT, -c(INDIVIDUALS, POP_ID)) %>% 
      mutate(
        POP_ID = factor(stri_replace_all_fixed(POP_ID, pop.levels, pop.labels, vectorize_all = FALSE), levels = unique(pop.labels), ordered =T),
        GT = stri_pad_left(str = GT, width = 6, pad = "0")
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
      strata.df <- input %>% distinct(INDIVIDUALS, POP_ID)
    } else {
      strata.df <- read_tsv(file = strata, col_names = TRUE, col_types = "cc") %>% 
        rename(POP_ID = STRATA)
      input <- left_join(input, strata.df, by = "INDIVIDUALS")
    }
    
    # Pop select
    if (!is.null(pop.select)) {
      message(stri_join(length(pop.select), "population(s) selected", sep = " "))
      pop.select <- c(pop.select, "mixture") # we don't want to remove the mixture samples, yet.
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
      tidyr::gather(INDIVIDUALS, GT, -LOCUS)
    
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
  } # End import haplotypes file
  
  # Blacklist genotypes ********************************************************
  if (is.null(blacklist.genotype)) { # no Whitelist
    message("Erasing genotype: no")
  } else {
    message("Erasing genotype: yes")
    blacklist.genotype <- read_tsv(blacklist.genotype, col_names = TRUE)
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
        whitelist.markers.ind <- input %>% distinct(CHROM, LOCUS, POS, INDIVIDUALS)
      } else {
        whitelist.markers.ind <- input %>% distinct(LOCUS, INDIVIDUALS)
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
  
  # mixture data  **************************************************************
  mixture.df <- read_tsv(file = mixture, col_names = TRUE, col_types = "c")
  
  # subsampling data ***********************************************************
  # Function:
  subsampling_data <- function(iteration.subsample, ...) {
    # message(paste0("Creating data subsample: ", iteration.subsample))
    if (is.null(subsample)) {
      subsample.select <- ind.pop.df %>% 
        mutate(SUBSAMPLE = rep(iteration.subsample, n()))
    } else {
      # separate all the mixture samples
      mixture.select <- ind.pop.df %>% filter(INDIVIDUALS %in% mixture.df$INDIVIDUALS)
      
      # subsample the baseline
      if (subsample > 1) { # integer
        subsample.select <- ind.pop.df %>%
          filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
          group_by(POP_ID) %>%
          sample_n(subsample, replace = FALSE) %>% # sampling individuals for each pop
          arrange(POP_ID, INDIVIDUALS)
        
        # Join baseline and mixture back in 1 dataset
        subsample.select <- bind_rows(subsample.select, mixture.select) %>% 
          mutate(SUBSAMPLE = rep(iteration.subsample, n())) %>% 
          ungroup()
      }
      if (subsample < 1) { # proportion
        subsample.select <- ind.pop.df %>%
          filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
          group_by(POP_ID) %>%
          sample_frac(subsample, replace = FALSE) %>% # sampling individuals for each pop
          arrange(POP_ID, INDIVIDUALS)
        
        # Join baseline and mixture back in 1 dataset
        subsample.select <- bind_rows(subsample.select, mixture.select) %>% 
          mutate(SUBSAMPLE = rep(iteration.subsample, n())) %>% 
          ungroup()
      }
    }
    return(subsample.select)
  } # End subsampling function
  
  # create the subsampling list
  ind.pop.df <- input %>% distinct(POP_ID, INDIVIDUALS)
  subsample.list <- map(.x = 1:iteration.subsample, .f = subsampling_data, subsample = subsample)
  
  # keep track of subsampling individuals and write to directory
  if (!is.null(subsample)) {
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
      snp.locus <- input %>% select(LOCUS, POS) %>% distinct(POS, .keep_all = TRUE)
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
        tidyr::unite(MARKERS, c(CHROM, LOCUS, POS), sep = "__")
    } # End Unique markers id
    
    # Markers in common between all populations (optional) *********************
    # This need to be moved while doing the assignment
    if (common.markers == TRUE) { # keep only markers present in all pop
      message("Using markers common in all populations:")
      pop.number <- input %>%
        filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
        select(POP_ID) %>% 
        filter(POP_ID != "mixture")
      
      pop.number <- n_distinct(pop.number$POP_ID)
      
      if (data.type == "vcf.file") pop.filter <- input %>% filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% filter(GT != "./.")
      if (data.type == "plink.file") pop.filter <- input %>% filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% filter(GT != "000")
      if (data.type == "df.file") pop.filter <- input %>% filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% filter(GT != "000000")
      if (data.type == "haplo.file") pop.filter <- input %>% filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% filter(GT != "-")
      
      pop.filter <- pop.filter %>% 
        group_by(MARKERS) %>%
        filter(n_distinct(POP_ID) == pop.number) %>%
        arrange(MARKERS) %>%
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
          filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
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
            filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
            tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
            tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
            select(MARKERS, GT, POP_ID) %>% 
            filter(GT != "000")
        }
        
        if (data.type == "plink.file") { # For PLINK and common code below
          maf.data <- input %>%
            filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
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
          distinct(MARKERS, POP_ID, .keep_all = TRUE) %>% 
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
                                   sep = "__", 
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
      # input <- input %>%
      #   mutate(
      #     REF= stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE), # replace nucleotide with numbers
      #     ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("1", "2", "3", "4"), vectorize_all = FALSE),# replace nucleotide with numbers
      #     GT = ifelse(GT == "0/0", stri_join(REF, REF, sep = "_"),
      #                 ifelse(GT == "1/1",  stri_join(ALT, ALT, sep = "_"),
      #                        ifelse(GT == "0/1", stri_join(REF, ALT, sep = "_"),
      #                               ifelse(GT == "1/0", stri_join(ALT, REF, sep = "_"), "0_0")
      #                        )
      #                 )
      #     )
      #   ) %>%
      #   arrange(MARKERS, POP_ID) %>%
      #   select(-c(REF, ALT))
      # 
      # gsi.prep <- input %>%
      #   tidyr::separate(col = GT, into = c("A1", "A2"), sep = "_") %>%  # separate the genotypes into alleles
      #   tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID))
      
      input <- input %>%
        mutate(
          REF= stri_replace_all_fixed(str = REF, pattern = c("A", "C", "G", "T"), replacement = c("001", "002", "003", "004"), vectorize_all = FALSE), # replace nucleotide with numbers
          ALT = stri_replace_all_fixed(str = ALT, pattern = c("A", "C", "G", "T"), replacement = c("001", "002", "003", "004"), vectorize_all = FALSE),# replace nucleotide with numbers
          GT = ifelse(GT == "0/0", stri_join(REF, REF, sep = ""),
                      ifelse(GT == "1/1",  stri_join(ALT, ALT, sep = ""),
                             ifelse(GT == "0/1", stri_join(REF, ALT, sep = ""),
                                    ifelse(GT == "1/0", stri_join(ALT, REF, sep = ""), "000000")
                             )
                      )
          )
        ) %>%
        arrange(MARKERS, POP_ID) %>%
        select(-c(REF, ALT)) 
      
      gsi.prep <- input %>%
        tidyr::separate(data = ., col = GT, into = .(A1, A2), sep = 3, remove = TRUE) %>% 
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) 
      
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
        ungroup () %>%
        mutate(
          INDIVIDUALS = as.character(INDIVIDUALS),
          POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
          POP_ID = factor(POP_ID) # xvalDapc does accept pop as ordered factor
        ) %>% 
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
        tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) 
    }
    if (data.type == "df.file") { # For data frame of genotypes
      gsi.prep <- input %>% 
        tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
        tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) 
    }
    if (data.type == "plink.file") { # for PLINK
      message("Recoding genotypes for gsi_sim")
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
            ungroup () %>%
            mutate(
              INDIVIDUALS = as.character(INDIVIDUALS),
              POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
              POP_ID = factor(POP_ID) # xvalDapc does accept pop as ordered factor
            ) %>% 
            arrange(POP_ID, INDIVIDUALS)
        )
      }
    }
    
    # Conversion to adegenet genind object
    if (assignment.analysis == "adegenet") {
      # genind arguments common to all data.type
      genind.prep <- genind.prep %>% arrange(POP_ID, INDIVIDUALS)
      ind <- genind.prep$INDIVIDUALS
      pop <- genind.prep$POP_ID
      genind.df <- genind.prep %>% ungroup() %>% 
        select(-c(INDIVIDUALS, POP_ID))
      rownames(genind.df) <- ind
      loc.names <- colnames(genind.df)
      strata <- genind.prep %>% 
        ungroup() %>% 
        distinct(INDIVIDUALS, POP_ID) %>% 
        mutate(MIXTURE = ifelse(INDIVIDUALS %in% mixture.df$INDIVIDUALS, "mixture", "baseline"))
      
      
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
    
    # only gsi_sim
    # if (assignment.analysis == "gsi_sim"){
    #   gsi.prep <- gsi.prep %>% 
    #     arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES)
    # }
    
    # Imputations **************************************************************
    if (imputation.method != "FALSE") {
      message("Preparing the data for imputations")
      
      # if (data.type == "vcf.file") { # for VCF input
      #   if (impute == "genotype") {
      #     input.prep <- input %>%
      #       select(-REF, -ALT) %>%
      #       mutate(GT = stri_replace_all_fixed(str = GT, pattern = "/", replacement = ":", vectorize_all = FALSE)) %>%
      #       mutate(
      #         GT = stri_replace_all_fixed(GT, pattern = ".:.", replacement = "NA", vectorize_all = FALSE),
      #         GT = replace(GT, which(GT == "NA"), NA)
      #       ) %>%
      #       group_by(INDIVIDUALS, POP_ID) %>% 
      #       tidyr::spread(data = ., key = MARKERS, value = GT) %>%
      #       ungroup() %>% 
      #       arrange(POP_ID, INDIVIDUALS)
      #   }
      #   
      #   if (impute == "allele") {
      #     input.prep <- input %>%
      #       select(-REF, -ALT) %>%
      #       mutate(GT = stri_replace_all_fixed(str = GT, pattern = "/", replacement = ":", vectorize_all = FALSE)) %>%
      #       mutate(
      #         GT = stri_replace_all_fixed(GT, pattern = ".:.", replacement = "NA:NA", vectorize_all = FALSE)
      #       ) %>% 
      #       tidyr::separate(col = GT, into = c("A1", "A2"), sep = ":") %>%  # separate the genotypes into alleles
      #       tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
      #       mutate(GT = replace(GT, which(GT == "NA"), NA)) %>%
      #       group_by(INDIVIDUALS, POP_ID, ALLELES) %>% 
      #       tidyr::spread(data = ., key = MARKERS, value = GT) %>% 
      #       ungroup() %>% 
      #       arrange(POP_ID, INDIVIDUALS)
      #   }
      #   
      # } # End VCF prep file for imputation
      
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
      
      if (data.type == "df.file" | data.type == "haplo.file" | data.type == "vcf.file") { # for df input
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
      
      # If imputations.group is by populations, we remove mixture samples conduct the imputations on
      # the baseline by populations, but the imputation for the mixture samples is conducted globally, automatically. 
      input.prep.mixture <- input.prep %>% 
        filter(INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
        mutate(POP_ID = droplevels(POP_ID))
      
      input.prep.baseline <- input.prep %>% 
        filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
        mutate(POP_ID = droplevels(POP_ID))
      
      strata.df.impute <- input.prep %>% 
        distinct(INDIVIDUALS, POP_ID)
      
      input.prep <- NULL # remove unused object
      
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
          message("Imputations computed by populations for baseline samples, take a break...")
          df.split.pop <- split(x = input.prep.baseline, f = input.prep.baseline$POP_ID) # slip data frame by population
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
            mc.cleanup = TRUE,
            mc.cores = parallel.core
          )
          
          # Compiling the results
          input.imp <- suppressWarnings(bind_rows(input.imp))
          
          # Second round of imputations (globally) to remove introduced NA 
          # In case that some pop don't have the markers
          input.imp <- suppressWarnings(plyr::colwise(factor, exclude = NA)(input.imp)) # Make the columns factor
          input.imp <- impute_genotype_rf(input.imp) # impute globally
          input.imp <- plyr::colwise(as.character, exclude = NA)(input.imp)
          
          # combine the mixture (no imputation) + the imputed baseline
          input.imp <- suppressWarnings(
            bind_rows(input.imp, input.prep.mixture) %>% 
              select(-POP_ID) # remove the column POP_ID
          )
          
          if (impute.mixture == TRUE) {
            # impute globally the mixture samples
            message("Imputations computed globally for mixture samples")
            input.imp <- suppressWarnings(plyr::colwise(factor, exclude = NA)(input.imp)) # Make the columns factor
            input.imp <- impute_genotype_rf(input.imp) # impute globally
            input.imp <- plyr::colwise(as.character, exclude = NA)(input.imp)
            input.imp <- suppressWarnings(
              left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                arrange(POP_ID, INDIVIDUALS) %>% 
                ungroup()
            )
          }
          
          if (impute.mixture == FALSE) {
            input.imp <- suppressWarnings(
              left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                arrange(POP_ID, INDIVIDUALS) %>% 
                ungroup()
            )
          }
          
          # dump unused objects
          df.split.pop <- NULL
          pop.list <- NULL
          sep.pop <- NULL
          imputed.dataset <- NULL
          input.prep <- NULL
          input.prep.mixture <- NULL
          input.prep.baseline <- NULL
          
        } # End imputation RF populations
        # Random Forest global
        if (imputations.group == "global") { # Globally (not by pop_id)
          message("Imputations computed globally for baseline samples, take a break...")
          input.imp <- plyr::colwise(factor, exclude = NA)(input.prep.baseline)
          input.imp <- input.imp %>% select(-POP_ID) # remove the column POP_ID
          input.imp <- impute_genotype_rf(input.imp)
          input.imp <- plyr::colwise(as.character, exclude = NA)(input.imp)
          input.prep <- NULL # remove unused object
          
          # combine the mixture (no imputation) + the imputed baseline
          input.imp <- suppressWarnings(
            bind_rows(input.imp, 
                      input.prep.mixture %>% 
                        select(-POP_ID)# remove the column POP_ID
            )
          )
          
          if (impute.mixture == TRUE) {
            # impute globally the mixture samples
            message("Imputations computed globally for mixture samples")
            input.imp <- suppressWarnings(plyr::colwise(factor, exclude = NA)(input.imp)) # Make the columns factor
            input.imp <- impute_genotype_rf(input.imp) # impute globally
            input.imp <- plyr::colwise(as.character, exclude = NA)(input.imp)
            input.imp <- suppressWarnings(
              left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                arrange(POP_ID, INDIVIDUALS) %>% 
                ungroup()
            )
          }
          
          if (impute.mixture == FALSE) {
            input.imp <- suppressWarnings(
              left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                arrange(POP_ID, INDIVIDUALS) %>% 
                ungroup()
            )
          }
        } # End imputation RF global
        
        # data prep
        if (impute == "genotype") {
          input.imp <- suppressWarnings(
            input.imp %>%
              tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID))
          )
          
          # Replace the GT == NA for "000000"
          if (impute.mixture == FALSE) {
            input.imp <- input.imp %>% 
              mutate(
                GT = stri_replace_na(GT, replacement = "000000")
              )
          }
        }
        if (impute == "allele") {
          input.imp <- suppressWarnings(
            input.imp %>%
              tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES))
          )
          # Replace the GT == NA for "000"
          if (impute.mixture == FALSE) {
            input.imp <- input.imp %>% 
              mutate(
                GT = stri_replace_na(GT, replacement = "000")
              )
          }
        }
      } # End imputation RF
      
      # Imputation using the most common genotype
      if (imputation.method == "max") { # End imputation max
        if (imputations.group == "populations") {
          message("Imputations computed by populations for baseline samples")
          message(stri_paste("Imputing: ", impute))
          
          if (impute == "genotype"){
            input.imp <- suppressWarnings(
              input.prep.baseline %>%
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
                ungroup() %>% 
                select(-POP_ID)
            )
            # combine the mixture (no imputation) + the imputed baseline
            input.imp <- suppressWarnings(
              bind_rows(input.imp, 
                        input.prep.mixture %>% 
                          tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>% 
                          select(-POP_ID)
              )
            )
            
            if (impute.mixture == TRUE) {
              # impute globally the mixture samples
              message("Imputations computed globally for mixture samples")
              input.imp <- suppressWarnings(
                input.imp %>% 
                  group_by(MARKERS) %>%
                  mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>% 
                  ungroup()
              )
              
              input.imp <- suppressWarnings(
                left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                  arrange(POP_ID, INDIVIDUALS, MARKERS) %>% 
                  ungroup()
              )
            }
            
            if (impute.mixture == FALSE) {
              input.imp <- suppressWarnings(
                left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                  arrange(POP_ID, INDIVIDUALS, MARKERS) %>% 
                  mutate(
                    GT = stri_replace_na(GT, replacement = "000000")
                  ) %>% 
                  ungroup()
              )
            } 
          }
          
          if (impute == "allele"){
            input.imp <- suppressWarnings(
              input.prep.baseline %>%
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
                ungroup() %>% 
                select(-POP_ID)
            )
            
            # combine the mixture (no imputation) + the imputed baseline
            input.imp <- suppressWarnings(
              bind_rows(input.imp, 
                        input.prep.mixture %>% 
                          tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
                          select(-POP_ID)
              )
            )
            
            if (impute.mixture == TRUE) {
              # impute globally the mixture samples
              message("Imputations computed globally for mixture samples")
              input.imp <- suppressWarnings(
                input.imp %>% 
                  group_by(MARKERS) %>%
                  mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
                  ungroup()
              )
              
              input.imp <- suppressWarnings(
                left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                  arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES) %>% 
                  ungroup()
              )
            }
            
            if (impute.mixture == FALSE) {
              input.imp <- suppressWarnings(
                left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                  arrange(POP_ID, INDIVIDUALS, MARKERS) %>% 
                  mutate(
                    GT = stri_replace_na(GT, replacement = "000")
                  ) %>% 
                  ungroup()
              )
            } 
          }
          
          input.prep <- NULL # remove unused object
          input.prep.baseline <- NULL
          input.prep.mixture <- NULL
          
        } # End imputation max populations 
        if (imputations.group == "global") {
          # Globally (not by pop_id)
          message("Imputations computed globally for baseline samples")
          message(stri_paste("Imputing: ", impute))
          
          if (impute == "genotype"){
            input.imp <- suppressWarnings(
              input.prep.baseline %>%
                tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>%
                group_by(MARKERS) %>%
                mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
                ungroup() %>% 
                select(-POP_ID)
            )
            
            # combine the mixture (no imputation) + the imputed baseline
            input.imp <- suppressWarnings(
              bind_rows(input.imp, 
                        input.prep.mixture %>% 
                          tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID)) %>% 
                          select(-POP_ID)
              )
            )
            
            if (impute.mixture == TRUE) {
              # impute globally the mixture samples
              message("Imputations computed globally for mixture samples")
              input.imp <- suppressWarnings(
                input.imp %>% 
                  group_by(MARKERS) %>%
                  mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>% 
                  ungroup()
              )
              
              input.imp <- suppressWarnings(
                left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                  arrange(POP_ID, INDIVIDUALS, MARKERS) %>% 
                  ungroup()
              )
            }
            
            if (impute.mixture == FALSE) {
              input.imp <- suppressWarnings(
                left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                  arrange(POP_ID, INDIVIDUALS, MARKERS) %>% 
                  mutate(
                    GT = stri_replace_na(GT, replacement = "000000")
                  ) %>% 
                  ungroup()
              )
            } 
          }
          
          if (impute == "allele"){
            input.imp <- suppressWarnings(
              input.prep.baseline %>%
                tidyr::gather(MARKERS, GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>%
                group_by(MARKERS) %>%
                mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
                ungroup() %>% 
                select(-POP_ID)
            )
            
            # combine the mixture (no imputation) + the imputed baseline
            input.imp <- suppressWarnings(
              bind_rows(input.imp, 
                        input.prep.mixture %>% 
                          tidyr::gather(key = MARKERS, value = GT, -c(INDIVIDUALS, POP_ID, ALLELES)) %>% 
                          select(-POP_ID)
              )
            )
            if (impute.mixture == TRUE) {
              # impute globally the mixture samples
              message("Imputations computed globally for mixture samples")
              input.imp <- suppressWarnings(
                input.imp %>% 
                  group_by(MARKERS) %>%
                  mutate(GT = stri_replace_na(GT, replacement = max(GT, na.rm = TRUE))) %>%
                  ungroup()
              )
              
              input.imp <- suppressWarnings(
                left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                  arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES) %>% 
                  ungroup()
              )
            }
            if (impute.mixture == FALSE) {
              input.imp <- suppressWarnings(
                left_join(strata.df.impute, input.imp, by = "INDIVIDUALS") %>% 
                  arrange(POP_ID, INDIVIDUALS, MARKERS) %>% 
                  mutate(
                    GT = stri_replace_na(GT, replacement = "000")
                  ) %>% 
                  ungroup()
              )
            } 
          }
          
          input.prep <- NULL # remove unused object
        } # End imputation max global 
      } # End imputations max
      
      # test <- input.imp %>% filter(GT == "000000")
      # test <- input.imp %>% filter(GT == "000")
      # test <- input.imp %>% ungroup %>% filter(is.na(GT))
      
      # prepare the imputed dataset for gsi_sim or adegenet
      message("Preparing imputed data set for assignement analysis")
      if (assignment.analysis == "gsi_sim") {
        if (impute == "genotype") {
            gsi.prep.imp <- input.imp %>%
              tidyr::separate(col = GT, into = c("A1", "A2"), sep = 3) %>%  # separate the genotypes into alleles
              tidyr::gather(key = ALLELES, GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>% 
              select(POP_ID, INDIVIDUALS, MARKERS, ALLELES, GT) %>% 
              arrange(POP_ID, INDIVIDUALS, MARKERS, ALLELES)
        }
        if (impute == "allele") {
          gsi.prep.imp <- input.imp
        }
      }
      
      # adegenet
      if (assignment.analysis == "adegenet") {
        if (impute == "genotype") {
            genind.prep.imp <- input.imp %>%
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
                ungroup () %>%
                mutate(
                  INDIVIDUALS = as.character(INDIVIDUALS),
                  POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
                  POP_ID = factor(POP_ID) # xvalDapc does accept pop as ordered factor
                ) %>% 
                arrange(POP_ID, INDIVIDUALS)
            )
        } # end impute genotype adegenet
        if (impute == "allele") {
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
          distinct(INDIVIDUALS, POP_ID) %>% 
          mutate(MIXTURE = ifelse(INDIVIDUALS %in% mixture.df$INDIVIDUALS, "mixture", "baseline"))
        
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
    # Assignment with gsi_sim
    if (assignment.analysis == "gsi_sim") {
      assignment_analysis <- function(data, select.markers, markers.names, missing.data, i, m, filename, ...) {
        # data <- gsi.prep #test
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
        
        # Baseline
        baseline.data <- suppressWarnings(
          data.select %>%
            filter(!INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
            arrange(POP_ID) %>%
            mutate(POP_ID = droplevels(POP_ID))
        )
        # gsi_sim baseline filename
        baseline.filename <- filename
        baseline.filename <- stri_replace_all_fixed(baseline.filename, pattern = ".txt",
                                                    replacement = "_baseline.txt")
        
        # save input fileto directory
        baseline.input <- assigner::write_gsi_sim(data = baseline.data, markers.names = markers.names, directory.subsample = directory.subsample, filename = baseline.filename, i = i, m = m)
        
        # Mixture
        mixture.data <- suppressWarnings(
          data.select %>%
            filter(INDIVIDUALS %in% mixture.df$INDIVIDUALS) %>% 
            mutate(POP_ID = factor(rep("mixture", n())))
        )
        
        # gsi_sim mixture filename
        mixture.filename <- filename
        mixture.filename <- stri_replace_all_fixed(mixture.filename, pattern = ".txt",
                                                   replacement = "_mixture.txt")
        
        # save input fileto directory
        mixture.input <- assigner::write_gsi_sim(data = mixture.data, markers.names = markers.names, directory.subsample = directory.subsample, filename = mixture.filename, i = i, m = m)
        
        # Run gsi_sim ------------------------------------------------------------
        baseline.input.gsi <- stri_join(directory.subsample, baseline.input)
        mixture.input.gsi <- stri_join(directory.subsample, mixture.input)
        
        output.gsi <- stri_replace_all_fixed(mixture.input, pattern = "txt", replacement = "output.txt")
        output.gsi <- stri_join(directory.subsample, output.gsi)
        setwd(directory.subsample)
        system(paste("gsi_sim -b ", baseline.input.gsi, "-t ", mixture.input.gsi, " > ", output.gsi))
        
        # Option remove the input file from directory to save space
        if (keep.gsi.files == FALSE) {
          file.remove(baseline.input.gsi)
          file.remove(mixture.input.gsi)
        }
        
        # Get Assignment results -------------------------------------------------
        # Number of markers
        n.locus <- m
        
        assignment <- suppressWarnings(
          read_delim(output.gsi, col_names = "ID", delim = "\t") %>%
            tidyr::separate(ID, c("KEEPER", "ASSIGN"), sep = ":/", extra = "warn") %>%
            filter(KEEPER == "GMA_FULL_EM_INDIVS_CSV") %>%
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
                ANALYSIS = rep("mixture", n()),
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
              select(INDIVIDUALS, ANALYSIS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
              arrange(CURRENT)
          )
        } else {
          assignment <- suppressWarnings(
            assignment %>%
              mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>% 
              left_join(strata.df, by = "INDIVIDUALS") %>%
              rename(CURRENT = POP_ID) %>% 
              mutate(
                ANALYSIS = rep("mixture", n()),
                CURRENT = factor(CURRENT, levels = unique(pop.labels), ordered =TRUE),
                # CURRENT = factor(CURRENT, levels = pop.levels, labels = pop.labels, ordered = TRUE),
                CURRENT = droplevels(CURRENT),
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
              select(INDIVIDUALS, ANALYSIS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
              arrange(CURRENT)
          )
          
          if (sampling.method == "random") {
            assignment <- assignment %>% 
              mutate(ITERATIONS = rep(i, n())) %>% 
              select(INDIVIDUALS, ANALYSIS, CURRENT, INFERRED, SCORE, SECOND_BEST_POP, SECOND_BEST_SCORE, MARKER_NUMBER, METHOD, MISSING_DATA, ITERATIONS) %>%
              arrange(CURRENT)
          }
          
        }
        
        
        # Option remove the output file from directory to save space
        if (keep.gsi.files == FALSE) file.remove(output.gsi)
        
        return(assignment)
      } # End assignment_analysis function
    }
    
    # Assignment with adegenet
    if (assignment.analysis == "adegenet") {
      assignment_analysis_adegenet <- function(data, select.markers, markers.names, missing.data, i, m, holdout, ...) {
        # data <- genind.object #test
        # missing.data <- "no.imputation" #test
        data.select <- data[loc = select.markers$MARKERS]
        
        # Run adegenet *********************************************************
        pop.data <- data.select@pop
        pop.data <- droplevels(pop.data)
        
        # # Alpha-Score DAPC
        # # When all the individuals are accounted for in the model construction
        # dapc.best.optim.a.score <- optim.a.score(dapc(data.select, n.da = length(levels(pop.data)), n.pca = round((length(indNames(data.select))/3)-1, 0)), pop = pop.data, plot = FALSE)$best
        # message(stri_paste("a-score optimisation for iteration:", i, sep = " ")) # message not working in parallel...
        # 
        # # DAPC with all the data
        # dapc.all <- dapc(data.select, n.da = length(levels(pop.data)), n.pca = dapc.best.optim.a.score, pop = pop.data)
        # message(stri_paste("DAPC iteration:", i, sep = " "))
        # message(stri_paste("DAPC marker group:", m, sep = " "))
        
        # Alpha-Score DAPC training data
        training.data <- data.select[data.select@strata$MIXTURE == "baseline"]
        # training.data <- data.select[!indNames(data.select) %in% holdout$INDIVIDUALS] # training dataset
        # indNames(training.data)
        # training.data@strata
        pop.training <- training.data@pop
        pop.training <- droplevels(pop.training)
        
        dapc.best.optim.a.score <- optim.a.score(dapc(training.data, n.da = length(levels(pop.training)), n.pca = round(((length(indNames(training.data))/3)-1), 0)), pop = pop.training, plot = FALSE)$best
        message(stri_paste("a-score optimisation for iteration:", i, sep = " "))
        
        dapc.training <- dapc(training.data, n.da = length(levels(pop.training)), n.pca = dapc.best.optim.a.score, pop = pop.training)
        message(stri_paste("DAPC of training data set for iteration:", i, sep = " "))
        
        # DAPC holdout individuals
        holdout.data <- data.select[data.select@strata$MIXTURE == "mixture"]
        # indNames(holdout.data)
        # holdout.data@strata
        # holdout.data <- data.select[indNames(data.select) %in% holdout$INDIVIDUALS] # holdout dataset
        pop.holdout <- holdout.data@pop
        pop.holdout <- droplevels(pop.holdout)
        assignment.levels <- levels(pop.data)
        rev.assignment.levels <- rev(assignment.levels)
        
        dapc.predict.holdout <- predict.dapc(dapc.training, newdata = holdout.data)
        message(stri_paste("Assigning holdout data for iteration:", i, sep = " "))
        
        
        # Get Assignment results -----------------------------------------------
        
        # Number of markers
        n.locus <- m
        
        if (sampling.method == "ranked") {
          i <- "not available with sampling.method = ranked"
        }
        
        assignment <- data.frame(INDIVIDUALS = indNames(holdout.data), POP_ID = pop.holdout, ASSIGN = dapc.predict.holdout$assign, dapc.predict.holdout$posterior) %>% 
          rename(CURRENT = POP_ID, INFERRED = ASSIGN) %>%
          mutate(
            ANALYSIS = rep("mixture", n()),
            MARKER_NUMBER = as.numeric(rep(n.locus, n())),
            METHOD = rep(sampling.method, n()),
            MISSING_DATA = rep(missing.data, n()),
            SUBSAMPLE = rep(subsample.id, n()),
            CURRENT = factor(CURRENT, levels = rev.assignment.levels, ordered = TRUE),
            INFERRED = factor(INFERRED, levels = assignment.levels, ordered = TRUE),
            ITERATIONS = rep(i, n())
          )
        return(assignment)
      } # End assignment_analysis_adegenet function
    } # End assignment_function
    
    # Random method ************************************************************
    if (sampling.method == "random") {
      message("Conducting Assignment analysis with markers selected randomly")
      # Number of times to repeat the sampling of markers
      iterations.list <- 1:iteration.method
      
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
      write_tsv(x = random.seed, path = paste0(directory.subsample, "random_seed_assignment_mixture.tsv"), col_names = TRUE, append = FALSE)
      
      message("Starting parallel computations for the assignment analysis
First sign of progress may take some time
Progress can be monitored with activity in the folder...")
      mrl <- NULL
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
        filename <- stri_replace_all_fixed(base.filename,
                                           pattern = ".txt",
                                           replacement = stri_join(
                                             "", "iteration", i, "markers", m, 
                                             "no_imputation.txt", sep = "_"
                                           )
        )
        if (assignment.analysis == "gsi_sim") {
          assignment.no.imp <- assignment_analysis(data = gsi.prep,
                                                   select.markers = select.markers,
                                                   markers.names = markers.names,
                                                   missing.data = "no.imputation", 
                                                   i = i, 
                                                   m = m,
                                                   filename = filename
          )
        }
        if (assignment.analysis == "adegenet") {
          assignment.no.imp <- assignment_analysis_adegenet(data = genind.object,
                                                            select.markers = select.markers,
                                                            markers.names = markers.names,
                                                            missing.data = "no.imputation", 
                                                            i = i, 
                                                            m = m
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
          filename <- stri_replace_all_fixed(base.filename,
                                             pattern = ".txt",
                                             replacement = stri_join(
                                               "", "iteration", i, "markers", m,
                                               "imputed.txt", sep = "_"
                                             )
          )
          
          if (assignment.analysis == "gsi_sim") {
            assignment.imp <- assignment_analysis(data = gsi.prep.imp,
                                                  select.markers = select.markers,
                                                  markers.names = markers.names,
                                                  missing.data = missing.data, 
                                                  i = i,
                                                  m = m,
                                                  filename = filename
            )
          }
          if (assignment.analysis == "adegenet") {
            assignment.imp <- assignment_analysis_adegenet(data = genind.object.imp,
                                                           select.markers = select.markers,
                                                           markers.names = markers.names,
                                                           missing.data = missing.data, 
                                                           i = i,
                                                           m = m
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
        mc.cleanup = TRUE,
        mc.cores = parallel.core
      )
      
      # Compiling the results
      message("Compiling results")
      assignment.res <- suppressWarnings(
        bind_rows(assignment.res) %>% 
          mutate(SUBSAMPLE = rep(subsample.id, n())) %>% 
          arrange(INDIVIDUALS, MARKER_NUMBER, MISSING_DATA, ITERATIONS)
      )
      
      # Write to the directory assignment results
      if (imputation.method == FALSE) {
        filename.assignment.res <- stri_join("assignment.mixture", "no.imputation", sampling.method, "tsv", sep = ".")
      } else { # with imputations
        filename.assignment.res <- stri_join("assignment.mixture", "imputed", sampling.method, "tsv", sep = ".")
      }
      write_tsv(x = assignment.res, path = paste0(directory.subsample, filename.assignment.res), col_names = TRUE, append = FALSE)
      
      if (assignment.analysis == "gsi_sim") {
        assignment.mixture.summary.stats <- assignment.res %>% 
          group_by(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, INFERRED, SUBSAMPLE) %>%
          summarise(
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
          mutate(
            TOTAL_ITERATIONS = rep(iteration.method, n())
          ) %>% 
          select(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE, INFERRED, NUMBER_ITERATIONS, TOTAL_ITERATIONS, MEAN_ITERATIONS, MEAN, SE, MIN, MAX, MEDIAN, QUANTILE25, QUANTILE75) %>% 
          arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE)
      }
      
      if (assignment.analysis == "adegenet") {
        assignment.mixture.summary.stats <- suppressWarnings(
          assignment.res %>%
            ungroup() %>%
            mutate(CURRENT = factor(CURRENT)) %>% 
            group_by(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE) %>%
            summarise(
              NUMBER_ITERATIONS = length(ITERATIONS),
              MEAN_ITERATIONS = round((NUMBER_ITERATIONS/iteration.method)*100, 2)
            ) %>%
            ungroup() %>%
            arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE, CURRENT, INFERRED)
        )
      }
      
      # Next step is common for gsi_sim and adegenet
      # Write the tables to directory
      # assignment summary stats
      if (imputation.method == FALSE) {
        filename.assignment.sum <- stri_join("assignment.mixture.summary.results", "no.imputation", sampling.method, "tsv", sep = ".")
      } else { # with imputations
        filename.assignment.sum <- stri_join("assignment.mixture.summary.results", "imputed", sampling.method, "tsv", sep = ".")
      }
      write_tsv(x = assignment.mixture.summary.stats, path = paste0(directory.subsample,filename.assignment.sum), col_names = TRUE, append = FALSE)
    } # End method random
    
    # Ranked method ************************************************************
    if (sampling.method == "ranked") {
      message("Conducting Assignment analysis with ranked markers")
      
      # List of all individuals
      ind.pop.df<- input %>% 
        ungroup %>% 
        distinct(POP_ID, INDIVIDUALS)
      
      message("Using thl method, ranking Fst with training samples...")
      holdout.individuals <- mixture.df
      
      write_tsv(x = holdout.individuals, 
                path = paste0(directory.subsample,"holdout.individuals.tsv"), 
                col_names = TRUE, 
                append = FALSE
      )
      message("Holdout samples = mixture samples: saved in your folder")
      
      # Going through the loop of holdout individuals
      message("Starting parallel computations for the assignment analysis
First sign of progress may take some time
Progress can be monitored with activity in the folder...")
      
      # assignment_ranking <- function(iterations.list, ...) {
      
      # Ranking Fst with training dataset (keep holdout individuals out)
      message("Ranking markers based on Fst with training samples")
      fst.ranked <- assigner::fst_WC84(
        data = input,
        holdout.samples = holdout.individuals$INDIVIDUALS
        )$fst.ranked
      
      write_tsv(
        x = fst.ranked, 
        path = paste0(directory.subsample, "fst_ranked.tsv"), 
        col_names = TRUE, 
        append = FALSE
      )
      if (imputation.method != FALSE) {
        fst.ranked.imp <- assigner::fst_WC84(
          data = input.imp, 
          holdout.samples = holdout.individuals$INDIVIDUALS
          )$fst.ranked
        write_tsv(x = fst.ranked.imp, 
                  path = paste0(directory.subsample, "fst_ranked_imputed.tsv"), 
                  col_names = TRUE, 
                  append = FALSE
        )
      }
      
      # Markers numbers loop function
      message("Going throught the marker.number")
      # assignment.marker <- list() # Create empty lists to feed the results
      
      i <- NULL
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
        filename <- stri_replace_all_fixed(base.filename,
                                           pattern = ".txt",
                                           replacement = stri_join(
                                             "", "markers", m, 
                                             "no_imputation.txt", sep = "_"
                                           )
        )
        
        
        if (assignment.analysis == "gsi_sim") {
          assignment.no.imp <- assignment_analysis(data = gsi.prep,
                                                   select.markers = select.markers,
                                                   markers.names = markers.names,
                                                   missing.data = "no.imputation",
                                                   i = NULL,
                                                   m = m,
                                                   filename = filename
          )
        }
        
        if (assignment.analysis == "adegenet") {
          assignment.no.imp <- assignment_analysis_adegenet(data = genind.object,
                                                            select.markers = select.markers,
                                                            markers.names = markers.names,
                                                            missing.data = "no.imputation", 
                                                            i = NULL, 
                                                            m = m
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
          filename <- stri_replace_all_fixed(base.filename,
                                             pattern = ".txt",
                                             replacement = stri_join(
                                               "", "markers", m, 
                                               "imputed.txt", sep = "_"
                                             )
          )
          
          if (assignment.analysis == "gsi_sim") {
            assignment.imp <- assignment_analysis(data = gsi.prep.imp,
                                                  select.markers = select.markers,
                                                  markers.names = markers.names,
                                                  missing.data = missing.data, 
                                                  i = NULL,
                                                  m = m,
                                                  filename = filename
            )
          }
          if (assignment.analysis == "adegenet") {
            assignment.imp <- assignment_analysis_adegenet(data = genind.object.imp,
                                                           select.markers = select.markers,
                                                           markers.names = markers.names,
                                                           missing.data = missing.data, 
                                                           i = NULL,
                                                           m = m
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
        return(assignment)
      }  # End marker number loop for both with and without imputations
      
      # using mclapply
      assignment.res <- list()
      assignment.res <- parallel::mclapply(
        X = marker.number, 
        FUN = assignment_marker_loop,
        mc.preschedule = FALSE, 
        mc.silent = FALSE, 
        mc.cleanup = TRUE,
        mc.cores = parallel.core,
        fst.ranked = fst.ranked,
        fst.ranked.imp = fst.ranked.imp,
        i = NULL,
        input = input,
        gsi.prep = gsi.prep,
        input.imp = input.imp,
        gsi.prep.imp = gsi.prep.imp,
        pop.levels = pop.levels,
        pop.labels = pop.labels,
        pop.id.start =  pop.id.start,
        pop.id.end = pop.id.end,
        sampling.method = sampling.method,
        iteration.method = iteration.method,
        filename = filename,
        keep.gsi.files = keep.gsi.files,
        imputation.method = imputation.method,
        parallel.core = parallel.core
      )
      
      # Compiling the results
      message("Compiling results")
      assignment.res <- suppressWarnings(
        bind_rows(assignment.res) %>% 
          mutate(SUBSAMPLE = rep(subsample.id, n())) %>% 
          arrange(INDIVIDUALS, MARKER_NUMBER, MISSING_DATA)
      )
      
      # Write to the directory assignment results
      if (imputation.method == FALSE) {
        filename.assignment.res <- stri_join("assignment.mixture", "no.imputation", sampling.method, "tsv", sep = ".")
      } else { # with imputations
        filename.assignment.res <- stri_join("assignment.mixture", "imputed", sampling.method, "tsv", sep = ".")
      }
      write_tsv(x = assignment.res, path = paste0(directory.subsample, filename.assignment.res), col_names = TRUE, append = FALSE)
    } # End of ranked thl method
    
    return(assignment.res)
  } # End assignment_function
  
  res <- map(.x = subsample.list, .f = assignment_function,
             assignment.analysis = assignment.analysis,
             mixture.df = mixture.df,
             strata.df = strata.df,
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
             filename = filename,
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
  
  if (!is.null(subsample)){
    write_tsv(x = res, path = paste0(directory, "assignment.mixture.results.tsv"), col_names = TRUE, append = FALSE)
  }
  
  # Summary of the subsampling iterations
  if (sampling.method == "random") {
    if (assignment.analysis == "gsi_sim") {
      assignment.mixture.summary.subsample <- res %>% 
        group_by(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, INFERRED, SUBSAMPLE) %>%
        summarise(
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
        mutate(
          TOTAL_ITERATIONS = rep(iteration.method, n())
        ) %>% 
        select(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE, INFERRED, NUMBER_ITERATIONS, TOTAL_ITERATIONS, MEAN_ITERATIONS, MEAN, SE, MIN, MAX, MEDIAN, QUANTILE25, QUANTILE75) %>% 
        arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE)
    }
    if (assignment.analysis == "adegenet") {
      assignment.mixture.summary.subsample <- res %>% 
        select(-X1, -X2) %>% 
        ungroup() %>%
        mutate(CURRENT = factor(CURRENT)) %>% 
        group_by(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE) %>%
        summarise(
          NUMBER_ITERATIONS = length(ITERATIONS),
          MEAN_ITERATIONS = round((NUMBER_ITERATIONS/iteration.method)*100, 2)
        ) %>%
        ungroup() %>%
        arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, SUBSAMPLE, CURRENT, INFERRED) %>% 
        group_by(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, MISSING_DATA) %>%
        summarise(
          MEAN_SUBSAMPLE = round(mean(MEAN_ITERATIONS), 2),
          SE = round(sqrt(stats::var(MEAN_ITERATIONS)/length(MEAN_ITERATIONS)), 2),
          MIN = round(min(MEAN_ITERATIONS), 2),
          MAX = round(max(MEAN_ITERATIONS), 2),
          MEDIAN = round(stats::median(MEAN_ITERATIONS), 2),
          QUANTILE25 = round(stats::quantile(MEAN_ITERATIONS, 0.25), 2),
          QUANTILE75 = round(stats::quantile(MEAN_ITERATIONS, 0.75), 2)
        ) %>%
        ungroup() %>%
        arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, MISSING_DATA, CURRENT, INFERRED)
    }
  } # end random
  if (sampling.method == "ranked") {
    if (assignment.analysis == "gsi_sim") {
      assignment.mixture.summary.subsample <- res %>% 
        group_by(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA, INFERRED) %>%
        summarise(
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
        mutate(
          TOTAL_SUBSAMPLE = rep(iteration.subsample, n())
        ) %>% 
        select(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA, INFERRED, NUMBER_SUBSAMPLE, TOTAL_SUBSAMPLE, MEAN_SUBSAMPLE, MEAN, SE, MIN, MAX, MEDIAN, QUANTILE25, QUANTILE75) %>% 
        arrange(INDIVIDUALS, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA)
    }
    if (assignment.analysis == "adegenet") {
      assignment.mixture.summary.subsample <- res %>%
        select(-X1, -X2, -ITERATIONS) %>%
        ungroup() %>%
        mutate(CURRENT = factor(CURRENT)) %>%
        group_by(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
        summarise(
          NUMBER_SUBSAMPLE = length(SUBSAMPLE),
          MEAN_SUBSAMPLE = round((NUMBER_SUBSAMPLE/iteration.subsample)*100, 2)) %>% 
        ungroup() %>% 
        mutate(TOTAL_SUBSAMPLE = rep(iteration.subsample, n())) %>% 
        select(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA, NUMBER_SUBSAMPLE, TOTAL_SUBSAMPLE, MEAN_SUBSAMPLE) %>% 
        arrange(INDIVIDUALS, CURRENT, INFERRED, ANALYSIS, MARKER_NUMBER, METHOD, MISSING_DATA)
    }
  } # end ranked
  
  
  
  # assignment summary results
  if (imputation.method == FALSE) {
    filename.assignment.sum <- stri_join("assignment.mixture.summary.results", "no.imputation", sampling.method, "tsv", sep = ".")
  } else { # with imputations
    filename.assignment.sum <- stri_join("assignment.mixture.summary.results", "imputed", sampling.method, "tsv", sep = ".")
  }
  write_tsv(x = assignment.mixture.summary.subsample, path = paste0(directory,filename.assignment.sum), col_names = TRUE, append = FALSE)
  
  # results
  res.list <- list(assignment = res, assignment.mixture.summary.results = assignment.mixture.summary.subsample)
  return(res.list)
} # End assignment_mixture

