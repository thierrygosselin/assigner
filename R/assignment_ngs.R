# Assignment analysis of massive parallel sequencing data
#' @name assignment_ngs

#' @title Assignment analysis of massive parallel sequencing data (GBS/RADseq, 
#' SNP chip, etc) using \href{https://github.com/eriqande/gsi_sim}{gsi_sim} and \code{\link[adegenet]{adegenet}}. 

#' @description
#' The arguments in the \code{assignment_ngs} function were tailored for the
#' reality of GBS/RADseq data for assignment analysis while
#' maintaining a reproducible workflow. The assignment analysis can be conducted
#' using \href{https://github.com/eriqande/gsi_sim}{gsi_sim}, a tool 
#' for doing and simulating genetic stock identification and 
#' developed by Eric C. Anderson, or 
#' \href{https://github.com/thibautjombart/adegenet}{adegenet}, 
#' a R package developed by Thibaul Jombart.
#' Various input files are offered. Individuals, populations and
#' markers can be filtered and/or selected in several ways using blacklist,
#' whitelist and other arguments. Map-independent imputation of missing genotype
#' using Random Forest or the most frequent category is also available.
#' Markers can be randomly selected for a classic LOO (Leave-One-Out)
#' assignment or chosen based on ranked Fst for a thl
#' (Training, Holdout, Leave-one-out) assignment analysis.

#' @param data 7 options: vcf, plink, stacks haplotype file, genind, genepop, 
#' and a data frame in wide or long/tidy format. \emph{See details}.
#' The function uses 
#' \href{https://github.com/thierrygosselin/stackr}{stackr} 
#' \code{\link[stackr]{read_long_tidy_wide}} and 
#' \code{\link[stackr]{tidy_genomic_data}}.

#' @param assignment.analysis Assignment analysis conducted with 
#' \code{assignment.analysis = "gsi_sim"} or 
#' \code{assignment.analysis = "adegenet"}.

#' @param sampling.method (character) Should the markers be randomly selected
#' \code{sampling.method == "random"} for a classic Leave-One-Out (LOO) assignment or
#' chosen based on ranked Fst \code{sampling.method == "ranked"}, used in a
#' Training-Holdout-Leave One Out (thl) assignment ?

#' @param adegenet.dapc.opt (optional, character) \strong{Argument available only when 
#' using:
#' \code{assignment.analysis = "adegenet"} with
#' \code{sampling.method == "random"}}.
#' 
#' The assignment using dapc can use the optimized alpha score 
#' \code{adegenet.dapc.opt == "optim.a.score"} or 
#' cross-validation \code{adegenet.dapc.opt == "xval"}
#' for stability of group membership probabilities. 
#' For fine tuning the trade-off between power of discrimination and over-fitting.
#' See \pkg{adegenet} documentation for more details. 
#' 
#' \code{adegenet.dapc.opt == "xval"} doesn't work with missing data, so it's 
#' only available with \strong{imputed data} (i.e. imputation.method == "rf" or "max").
#' With non imputed data or the default: \code{adegenet.dapc.opt == "optim.a.score"}.


#' @param adegenet.n.rep (optional, integer) 
#' When \code{adegenet.dapc.opt == "xval"}, the number of replicates to be 
#' carried out at each level of PC retention. 
#' Default: \code{adegenet.n.rep = 30}.
#' See \pkg{adegenet} documentation for more details.
#' @param adegenet.training (optional, numeric) 
#' When \code{adegenet.dapc.opt == "xval"}, the proportion of data (individuals) 
#' to be used for the training set.
#' Default: \code{adegenet.training = 0.9}, if all groups have >= 10 members. 
#' Otherwise, training.set scales automatically to the largest proportion 
#' that still ensures all groups will be present in both training 
#' and validation sets.
#' See \pkg{adegenet} documentation for more details.

#' @param thl (character, integer, proportion) For \code{sampling.method = "ranked"} only.
#' Default \code{thl = 1}, 1 individual sample is used as holdout. This individual is not
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

#' @param subsample (Integer or Proportion) Default is no subsampling, \code{subsample = NULL}.
#' With a proportion argument \code{subsample = 0.15}, 15 percent of individuals
#' in each populations are chosen randomly to represent the dataset.
#' With \code{subsample = 36}, 36 individuals in each populations are chosen
#' randomly to represent the dataset.

#' @param iteration.subsample (Integer) The number of iterations to repeat 
#' subsampling, default: \code{iteration.subsample = 1}.
#' With \code{subsample = 20} and \code{iteration.subsample = 10},
#' 20 individuals/populations will be randomly chosen 10 times.

#' @param marker.number (Integer or string of number or "all") Calculations with
#' fixed or subsample of your markers. Default= \code{"all"}.
#' e.g. To test 500, 1000, 2000 and all  the markers:
#' \code{marker.number = c(500, 1000, 2000, "all")}.
#' To use only 500 makers \code{marker.number = 500}.


#' @param blacklist.id (optional) A blacklist with individual ID and
#' a column header 'INDIVIDUALS'. The blacklist is in the working directory
#' (e.g. "blacklist.txt").
#' Default: \code{blacklist.id = NULL}.

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
#' with 'un' need to be changed to 1. 
#' Default: \code{blacklist.genotype = NULL} for no blacklist of 
#' genotypes to erase.

#' @param whitelist.markers (optional) A whitelist containing CHROM (character
#' or integer) and/or LOCUS (integer) and/or
#' POS (integer) columns header. To filter by chromosome and/or locus and/or by snp.
#' The whitelist is in the working directory (e.g. "whitelist.txt").
#' de novo CHROM column with 'un' need to be changed to 1. 
#' In the VCF, the column ID is the LOCUS identification.
#' Default \code{whitelist.markers = NULL} for no whitelist of markers.

#' @param monomorphic.out (optional) Should the monomorphic 
#' markers present in the dataset be filtered out ? 
#' Default: \code{monomorphic.out = TRUE}.

#' @param snp.ld (optional) \strong{For VCF file only}. 
#' SNP short distance linkage disequilibrium pruning. With anonymous markers from
#' RADseq/GBS de novo discovery, you can minimize linkage disequilibrium (LD) by
#' choosing among these 3 options: \code{"random"} selection, \code{"first"} or
#' \code{"last"} SNP on the same short read/haplotype. For long distance linkage
#' disequilibrium pruning, see details below.
#' Default: \code{snp.ld = NULL}.

#' @param common.markers (optional) Logical. Default: \code{common.markers = TRUE}, 
#' will only keep markers in common (genotyped) between all the populations.

#' @param maf.thresholds (string, double, optional) String with 
#' local/populations and global/overall maf thresholds, respectively.
#' e.g. \code{maf.thresholds = c(0.05, 0.1)} for a local maf threshold 
#' of 0.05 and a global threshold of 0.1. Available for VCF, PLINK and data frame 
#' files. Use stackr for haplotypes files and use the whitelist of markers.
#' Default: \code{maf.thresholds = NULL}. 

#' @param maf.pop.num.threshold (integer, optional) When maf thresholds are used,
#' this argument is for the number of pop required to pass the maf thresholds
#' to keep the locus. Default: \code{maf.pop.num.threshold = 1}

#' @param maf.approach (character, optional). By \code{maf.approach = "SNP"} or 
#' by \code{maf.approach = "haplotype"}.
#' The function will consider the SNP or ID/LOCUS/haplotype/read MAF statistics 
#' to filter the markers.
#' Default is \code{maf.approach = "SNP"}. The \code{haplotype} approach is 
#' restricted to VCF file.

#' @param maf.operator (character, optional) \code{maf.operator = "AND"} or 
#' default \code{maf.operator = "OR"}.
#' When filtering over LOCUS or SNP, do you want the local \code{"AND"}
#' global MAF to pass the thresholds, or ... you want the local \code{"OR"}
#' global MAF to pass the thresholds, to keep the marker?

#' @param max.marker An optional integer useful to subsample marker number in 
#' large PLINK file. e.g. if the data set 
#' contains 200 000 markers and \code{max.marker = 10000} 10000 markers are
#' subsampled randomly from the 200000 markers. Use \code{whitelist.markers} to
#' keep specific markers.
#' Default: \code{max.marker = NULL}.

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
#' Default: \code{pop.labels = NULL}. When pop.levels is not null and pop.labels
#' is not specified. pop.labels = pop.levels.
#' If you find this too complicated, there is also the
#' \code{strata} argument that can do the same thing, see below.

#' @param strata (optional for data frame and PLINK files, 
#' required for VCF and haplotypes files) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping. 
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' Default: \code{strata = NULL}.

#' @param pop.select (string, optional) Selected list of populations for 
#' the analysis. e.g. \code{pop.select = c("QUE", "ONT")} to select \code{QUE}
#'and \code{ONT} population samples (out of 20 pops).
# Default: \code{pop.select = NULL} 


#' @param folder (optional) The name of the folder created in the working 
#' directory to save the files/results. Default: \code{folder = NULL} will create
#' the folder for you with data and time.

#' @param filename (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. 
#' Several info will be appended to the name of the file.
#' Default \code{assignment_data.txt}.

#' @param keep.gsi.files (logical, optional) With the default, 
#' the intermediate input and output gsi_sim files will be deleted from the 
#' directory when finished processing. I you decide to keep the files, with 
#' \code{keep.gsi.files =TRUE}, remember to allocate a large chunk of the disk 
#' space for the analysis.
#' Default: \code{keep.gsi.files = FALSE} 

#' @param imputation.method (character, optional) 
#' Methods available for map-independent imputations of missing genotype: 
#' (1) \code{"max"} to use the most frequent category for imputations.
#' (2) \code{"rf"} using Random Forest algorithm. 
#' Default: no imputation \code{imputation.method = NULL}.

#' @param impute (character, optional) Imputation on missing genotype 
#' \code{impute = "genotype"} or alleles \code{impute = "allele"}.
#' Default: \code{"genotype"}.

#' @param imputations.group (character, optional) \code{"global"} or \code{"populations"}.
#' Should the imputations be computed globally or by populations. If you choose
#' global, turn the verbose to \code{TRUE}, to see progress.
#' Default = \code{"populations"}.

#' @param num.tree (integer, optional) The number of trees to grow in Random Forest. 
#' Default: \code{num.tree = 100}.

#' @param iteration.rf (integer, optional) The number of iterations of missing data algorithm
#' in Random Forest. 
#' Default: \code{iteration.rf = 10}.

#' @param split.number (integer, optional) Non-negative integer value used to specify
#' random splitting in Random Forest. 
#' Default: \code{split.number = 100}.

#' @param verbose (logical, optional) Should trace output be enabled on each iteration
#' in Random Forest ? 
#' Default: \code{verbose = FALSE}.

#' @param parallel.core (optional) The number of core for OpenMP shared-memory parallel
#' programming of Random Forest imputations. For more info on how to install the
#' OpenMP version see \code{\link[randomForestSRC]{randomForestSRC-package}}.
#' If not selected \code{detectCores()-1} is used as default.

#' @details 
#' \strong{Input files:}
#' \enumerate{
#' \item VCF file (e.g. \code{data = "batch_1.vcf"}). 
#' To make the VCF population ready, you need the \code{strata} argument.
#' 
#' \item haplotype file created in STACKS (e.g. \code{data = "batch_1.haplotypes.tsv"}).
#' To make the haplotype file population ready, you need the \code{strata} argument.
#' 
#' \item Data frame
#' To discriminate the long from the wide format, 
#' the function \pkg{stackr} \code{\link[stackr]{read_long_tidy_wide}} searches 
#' for "MARKERS" in column names (TRUE = long format).
#' The data frame is tab delimitted.

#' \strong{Wide format:}
#' The wide format cannot store metadata info.
#' The wide format starts with these 2 id columns: 
#' \code{INDIVIDUALS}, \code{POP_ID} (that refers to any grouping of individuals), 
#' the remaining columns are the markers in separate columns storing genotypes.
#' 
#' \strong{Long/Tidy format:}
#' The long format is considered to be a tidy data frame and can store metadata info. 
#' (e.g. from a VCF see \pkg{stackr} \code{\link[stackr]{tidy_genomic_data}}). A minimum of 4 columns
#' are required in the long format: \code{INDIVIDUALS}, \code{POP_ID}, 
#' \code{MARKERS} and \code{GENOTYPE or GT}. The rest are considered metata info.
#' 
#' \strong{2 genotypes formats are available:}
#' 6 characters no separator: e.g. \code{001002 of 111333} (for heterozygote individual).
#' 6 characters WITH separator: e.g. \code{001/002 of 111/333} (for heterozygote individual).
#' The separator can be any of these: \code{"/", ":", "_", "-", "."}.
#' 
#' \emph{How to get a tidy data frame ?}
#' \pkg{stackr} \code{\link[stackr]{tidy_genomic_data}} can transform 6 genomic data formats 
#' in a tidy data frame.
#' 
#' \item PLINK file in 
#' \code{tped/tfam} format (e.g. \code{data =  "data.assignment.tped"}). 
#' The first 2 columns of the \code{tfam} file will be used for the 
#' \code{strata} argument below, unless a new one is provided. 
#' Columns 1, 3 and 4 of the \code{tped} are discarded. The remaining columns 
#' correspond to the genotype in the format \code{01/04} 
#' where \code{A = 01, C = 02, G = 03 and T = 04}. For \code{A/T} format, use 
#' PLINK or bash to convert.
#' Use \href{http://vcftools.sourceforge.net/}{VCFTOOLS} with \code{--plink-tped} 
#' to convert very large VCF file. For \code{.ped} file conversion to 
#' \code{.tped} use \href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK} 
#' with \code{--recode transpose},
#' 
#' \item \code{\link[adegenet]{genind}} object from \code{\link[adegenet]{adegenet}}.
#' 
#' \item genepop data file (e.g. \code{data = kiwi_data.gen}). Here, the function can only use
#' alleles encoded with 3 digits.
#' }

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

#' @note \code{assignment_ngs} assumes that the command line version of 
#' \href{https://github.com/eriqande/gsi_sim}{gsi_sim} 
#' is properly installed into \code{file.path(system.file(package = "assigner"), "bin", "gsi_sim")}.
#' Things are set up so that it will try running gsi_sim, and if it does not find it, the 
#' program will throw an error and ask the user to run \code{\link{install_gsi_sim}}
#' which will do its best to put a usable copy of gsi_sim where it is needed.  To do 
#' so, you must be connected to the internet. If that doesn't work, you will
#' need to compile the program yourself, or get it yourself, and the manually copy
#' it to \code{file.path(system.file(package = "assigner"), "bin", "gsi_sim")}.
#' To compile \href{https://github.com/eriqande/gsi_sim}{gsi_sim}, follow the 
#' instruction here: \url{https://github.com/eriqande/gsi_sim}.

#' @export
#' @rdname assignment_ngs
#' @import parallel
#' @import stringi
#' @import adegenet
#' @import dplyr
#' @import stackr
#' @importFrom stats var median quantile
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
#' filename = "treefrog.txt",
#' keep.gsi.files = FALSE, 
#' strata = "strata.treefrog.tsv",
#' pop.levels = c("PAN", "COS")
#' imputation.method = NULL,
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
#' To view the full range of y values = Assignment success(%): 
#' assignment.treefrog$plot.assignment + 
#' facet_grid(SUBSAMPLE~CURRENT) + 
#' scale_y_continuous(limits = c(0,100)) 
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
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795
#' @references Jombart T, Devillard S, Balloux F. 
#' Discriminant analysis of principal components: a new method for the analysis 
#' of genetically structured populations. 
#' BMC Genet. 2010:11: 94. doi:10.1186/1471-2156-11-94
#' @references Jombart T, Ahmed I. adegenet 1.3-1: new tools for the analysis 
#' of genome-wide SNP data. 
#' Bioinformatics. 2011:27: 3070–3071. doi:10.1093/bioinformatics/btr521

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
      "COUNT", "MAX_COUNT_MARKERS", "hierarchy", "GT_VCF"
    )
  )
}

assignment_ngs <- function(
  data,
  assignment.analysis,
  sampling.method,
  adegenet.dapc.opt = "optim.a.score",
  adegenet.n.rep = 30,
  adegenet.training = 0.9,
  thl = 1,
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
  strata = NULL,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  imputation.method = NULL,
  impute = "genotype",
  imputations.group = "populations",
  num.tree = 100,
  iteration.rf = 10,
  split.number = 100,
  verbose = FALSE,
  folder = NULL,
  filename = "assignment_data.txt",
  keep.gsi.files = FALSE,
  parallel.core = detectCores()-1,
  ...) {
  
  # Checking for missing and/or default arguments ******************************
  if (missing(data)) stop("Input file missing")
  if (missing(assignment.analysis)) stop("assignment.analysis argument missing")
  if (missing(sampling.method)) stop("sampling.method argument missing")
  if (assignment.analysis == "gsi_sim" & !gsi_sim_exists()){
    stop("Can't find the gsi_sim executable where it was expected at ", gsi_sim_binary_path(), ".  
         If you have internet access, you can install it
         from within R by invoking the function \"install_gsi_sim(fromSource = TRUE)\"")
  }
  if (assignment.analysis == "gsi_sim") message("Assignment analysis with gsi_sim")
  if (assignment.analysis == "adegenet") message("Assignment analysis with adegenet")
  if (!is.null(pop.levels) & is.null(pop.labels)) pop.labels <- pop.levels
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  
  # Create a folder based on filename to save the output files *****************
  if (is.null(folder)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    if (is.null(imputation.method)) {
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
  
  # File type detection ********************************************************
  if(adegenet::is.genind(data)){
    data.type <- "genind.file"
    message("File type: genind object")
  } else {
    data.type <- readChar(con = data, nchars = 16L, useBytes = TRUE)
    if (identical(data.type, "##fileformat=VCF") | stri_detect_fixed(str = data, pattern = ".vcf")) {
      data.type <- "vcf.file"
      # message("File type: VCF")
    }
    if (stri_detect_fixed(str = data, pattern = ".tped")) {
      data.type <- "plink.file"
      # message("File type: PLINK")
      if (!file.exists(stri_replace_all_fixed(str = data, pattern = ".tped", replacement = ".tfam", vectorize_all = FALSE))) {
        stop("Missing tfam file with the same prefix as your tped")
      }
    } 
    if (stri_detect_fixed(str = data.type, pattern = "POP_ID") | stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") | stri_detect_fixed(str = data.type, pattern = "MARKERS")) {
      data.type <- "df.file"
      # message("File type: data frame of genotypes")
    }
    if (stri_detect_fixed(str = data.type, pattern = "Catalog")) {
      # data.type <- "haplo.file"
      message("File type: haplotypes from stacks")
      if (is.null(blacklist.genotype)) {
        stop("blacklist.genotype file missing. 
             Use stackr's missing_genotypes function to create this blacklist")
      }
    }
    if (stri_detect_fixed(str = data, pattern = ".gen")) {
      # data.type <- "genepop.file"
      message("File type: genepop")
    } 
    
  } # end file type detection
  
  # Strata argument required for VCF and haplotypes files **********************
  if (data.type == "haplo.file" | data.type == "vcf.file") {
    if (is.null(strata)) stop("strata argument is required")
  }
  
  # Import input ***************************************************************
  input <- stackr::tidy_genomic_data(
    data = data, 
    vcf.metadata = FALSE,
    blacklist.id = blacklist.id, 
    blacklist.genotype = blacklist.genotype, 
    whitelist.markers = whitelist.markers, 
    monomorphic.out = monomorphic.out, 
    max.marker = max.marker,
    # snp.ld = NULL, 
    common.markers = FALSE, 
    # maf.thresholds = NULL, 
    # maf.pop.num.threshold = 1, 
    # maf.approach = "snp", 
    # maf.operator = "or",
    strata = strata, 
    pop.levels = pop.levels, 
    pop.labels = pop.labels, 
    pop.select = pop.select,
    filename = NULL
  )
  
  # create a strata.df
  strata.df <- input %>% 
    select(INDIVIDUALS, POP_ID) %>% 
    distinct(INDIVIDUALS)
  strata <- strata.df
  pop.levels <- levels(input$POP_ID)
  pop.labels <- pop.levels
  
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
    write_tsv(
      x = subsampling.individuals, 
      path = paste0(directory, "subsampling.individuals.tsv"), 
      col_names = TRUE, 
      append = FALSE
    )
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
    
    # # Unique markers id *********************************************************
    # if (data.type != "vcf.file") {
    #   input <- input %>% rename(MARKERS = LOCUS)
    # }# End Unique markers id
    
    # Markers in common between all populations (optional) *********************
    if (common.markers) { # keep only markers present in all pop
      message("Using markers common in all populations:")
      pop.number <- n_distinct(input$POP_ID)
      
      pop.filter <- input %>% filter(GT != "000000")
      
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
          filter(GT_VCF != "./.") %>%
          group_by(MARKERS, POP_ID, REF, ALT) %>%
          summarise(
            N = as.numeric(n()),
            PQ = as.numeric(length(GT_VCF[GT_VCF == "1/0" | GT_VCF == "0/1"])),
            QQ = as.numeric(length(GT_VCF[GT_VCF == "1/1"]))
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
        
        # We split the alleles here to prep for MAF
        maf.data <- input %>%
          select(MARKERS,POP_ID, INDIVIDUALS, GT) %>%
          tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
          tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
          filter(GT != "000")
        
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
                path = "maf.data.tsv",
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
    
    # Adegenet  ****************************************************************
    
    if (assignment.analysis == "adegenet" ) {
      message("Preparing adegenet genind object")
      genind.prep <- suppressWarnings(
        input %>%
          tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
          tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID))  %>%
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
          mutate(
            INDIVIDUALS = as.character(INDIVIDUALS),
            POP_ID = as.character(POP_ID), # required to be able to do xvalDapc with adegenet.
            POP_ID = factor(POP_ID) # xvalDapc does accept pop as ordered factor
          ) %>% 
          arrange(POP_ID, INDIVIDUALS)
      )
      
      # genind arguments common to all data.type
      ind <- genind.prep$INDIVIDUALS
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
      loc.names <- NULL
      strata <- NULL
      prevcall <- NULL
      genind.prep <- NULL
    }
    
    # Imputations **************************************************************
    if (!is.null(imputation.method)) {
      
      input.imp <- stackr::stackr_imputations_module(
        data = input, 
        imputation.method = imputation.method, 
        impute = impute, 
        imputations.group = imputations.group, 
        num.tree = num.tree, 
        iteration.rf = iteration.rf, 
        split.number = split.number, 
        verbose = verbose, 
        parallel.core = parallel.core, 
        filename = "dataset.imputed.tsv"
      )
      
      # prepare the imputed dataset for gsi_sim or adegenet
      message("Preparing imputed data set for assignement analysis")
      
      # adegenet
      if (assignment.analysis == "adegenet") {
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
        
        # genind arguments common to all data.type
        ind <- genind.prep.imp$INDIVIDUALS
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
        loc.names <- NULL
        strata <- NULL
        prevcall <- NULL
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
    
    
    marker.number <- as.numeric(
      stri_replace_all_fixed(
        str = marker.number, pattern = "all", replacement = nrow(unique.markers), 
        vectorize_all = TRUE
      )
    )
    
    # In marker.number, remove marker group higher than the max number of markers
    removing.marker <- purrr::keep(.x = marker.number, .p = marker.number > nrow(unique.markers))
    
    if (length(removing.marker) > 0) {
      message(
        "Removing marker.number higher than the max number of markers: ", 
        stri_c(removing.marker, collapse = ", ")
      )
    }
    marker.number <- purrr::discard(.x = marker.number, .p = marker.number > nrow(unique.markers))
    
    # # remove marker number below 10
    # removing.low.marker.number <- purrr::keep(.x = marker.number, .p = marker.number < 10)
    # if (length(removing.low.marker.number) > 0) {
    #   message(
    #     "Removing marker.number lower than the min (10 markers): ", 
    #     stri_c(removing.low.marker.number, collapse = ", ")
    #   )
    # }
    # marker.number <- purrr::discard(.x = marker.number, .p = marker.number < 10)
    # 
    # Functions ******************************************************************
    # Assignment with gsi_sim
    assignment_analysis <- function(
      data, select.markers, 
      markers.names, missing.data, i, m, holdout, filename, ...) {
      # data <- input #test
      # data <- genind.prep #test
      # data <- genind.object.imp # test
      # missing.data <- "no.imputation" #test
      
      data.select <- suppressWarnings(
        data %>%
          semi_join(select.markers, by = "MARKERS") %>%
          arrange(POP_ID, INDIVIDUALS, MARKERS)
      )
      
      # Write gsi_sim input file to directory
      input.gsi <- assigner::write_gsi_sim(
        data = data.select, 
        pop.levels = pop.levels, 
        pop.labels = pop.labels, 
        strata = NULL, 
        filename = filename
      )
      
      # Run gsi_sim ------------------------------------------------------------
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
      
      assignment <- suppressWarnings(
        assignment %>%
          mutate(INDIVIDUALS = as.character(INDIVIDUALS)) %>% 
          left_join(strata.df, by = "INDIVIDUALS") %>%
          rename(CURRENT = POP_ID) %>% 
          mutate(
            INFERRED = factor(INFERRED, levels = pop.levels, ordered = TRUE),
            INFERRED = droplevels(INFERRED),
            SECOND_BEST_POP = factor(SECOND_BEST_POP, levels = pop.levels, ordered = TRUE),
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
              CURRENT = factor(CURRENT, levels = pop.levels, ordered = TRUE),
              CURRENT = droplevels(CURRENT),
              ITERATIONS = rep(i, n())
            )
        }
      }
      return(assignment)
    } # End assignment_analysis function
    
    # Assignment with adegenet
    assignment_analysis_adegenet <- function(
      data, select.markers, 
      adegenet.dapc.opt, adegenet.n.rep, adegenet.training, 
      parallel.core, markers.names, missing.data, i, m, holdout, ...) {
      # data <- genind.object #test
      # missing.data <- "no.imputation" #test
      data.select <- data[loc = select.markers$MARKERS]
      
      # Run adegenet *********************************************************
      pop.data <- data.select@pop
      pop.data <- droplevels(pop.data)
      
      
      if (sampling.method == "random") {
        # DAPC optimized a-score 
        if (adegenet.dapc.opt == "optim.a.score") {
          dapc.best.optim.a.score <- optim.a.score(
            dapc(data.select, 
                 n.da = length(levels(pop.data)), 
                 n.pca = round((length(indNames(data.select))/3)-1, 0)), 
            pop = pop.data, 
            plot = FALSE
          )$best
          message(stri_paste("a-score optimisation for iteration:", i, sep = " ")) # message not working in parallel...
          
          # DAPC
          dapc.assignment <- dapc(data.select, n.da = length(levels(pop.data)), n.pca = dapc.best.optim.a.score, pop = pop.data)
          message(stri_paste("DAPC iteration:", i, sep = " "))
          message(stri_paste("DAPC marker group:", m, sep = " "))
        }
        
        # DAPC with Cross-Validation
        if (adegenet.dapc.opt == "xval") {
          dapc.assignment <- xvalDapc(
            x = data.select@tab, 
            grp = pop(data.select),
            n.da = length(levels(pop.data)),
            n.pca.max = round((length(indNames(data.select))/3)-1, 0), 
            n.rep = adegenet.n.rep , 
            training.set = adegenet.training, 
            result = "groupMean", 
            center = TRUE, 
            scale = FALSE, 
            xval.plot = FALSE, 
            parallel = "multicore", 
            ncpus = parallel.core
          )$DAPC
          
          message(stri_paste("DAPC iteration:", i, sep = " "))
          message(stri_paste("DAPC marker group:", m, sep = " "))
        }
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
        assignment <- data_frame(ASSIGNMENT_PERC = summary(dapc.assignment)$assign.per.pop*100) %>% 
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
      # iterations.list <- 1:10 # test
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
      
      message(
        "Starting parallel computations for the assignment analysis.
First sign of progress may take some time.
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
        filename.no.imp <- stri_paste(directory.subsample, filename, sep = "")
        filename.no.imp <- stri_replace_all_fixed(
          filename.no.imp,
          pattern = "txt",
          replacement = stri_join(
            i, m, 
            "no.imputation", "txt", sep = "."
          )
        )
        
        if (assignment.analysis == "gsi_sim") {
          assignment.no.imp <- assignment_analysis(
            data = input,
            select.markers = select.markers,
            markers.names = markers.names,
            missing.data = "no.imputation", 
            i = i, 
            m = m,
            holdout = NULL,
            filename = filename.no.imp
          )
        }
        if (assignment.analysis == "adegenet") {
          assignment.no.imp <- assignment_analysis_adegenet(
            data = genind.object,
            select.markers = select.markers,
            adegenet.dapc.opt = "optim.a.score",
            adegenet.n.rep = adegenet.n.rep, 
            adegenet.training = adegenet.training,
            markers.names = markers.names,
            missing.data = "no.imputation", 
            i = i, 
            m = m,
            holdout = NULL
          )
        }
        
        # With imputations
        if (!is.null(imputation.method)) {# with imputations
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
          filename.imp <- stri_paste(directory.subsample, filename, sep = "")
          filename.imp <- stri_replace_all_fixed(
            filename.imp,
            pattern = "txt",
            replacement = stri_join(
              i, m, 
              "imputed", "txt", sep = "."
            )
          )
          if (assignment.analysis == "gsi_sim") {
            assignment.imp <- assignment_analysis(
              data = input.imp,
              select.markers = select.markers,
              markers.names = markers.names,
              missing.data = missing.data, 
              i = i,
              m = m,
              holdout = holdout,
              filename = filename.imp
            )
          }
          if (assignment.analysis == "adegenet") {
            assignment.imp <- assignment_analysis_adegenet(
              data = genind.object.imp,
              select.markers = select.markers,
              adegenet.dapc.opt = adegenet.dapc.opt, 
              adegenet.n.rep = adegenet.n.rep, 
              adegenet.training = adegenet.training, 
              parallel.core = parallel.core,
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
        
        # Compile assignment results each marker number for the iteration
        if (is.null(imputation.method)) {
          assignment <- assignment.no.imp
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
              CURRENT = factor(CURRENT, levels = pop.levels, ordered = TRUE),
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
          if (is.null(imputation.method)) {
            filename.assignment.res <- stri_join(
              "assignment", sampling.method, "no.imputation", "results", 
              "individuals", "iterations", "tsv", sep = "."
            )
          } else { # with imputations
            filename.assignment.res <- stri_join(
              "assignment", sampling.method, "imputed", "results", "individuals", 
              "iterations", "tsv", sep = "."
            )
          }
        } else {# with subsampling
          if (is.null(imputation.method)) {
            filename.assignment.res <- stri_join(
              "assignment", sampling.method, "no.imputation", "results", "individuals",
              "iterations", "subsample", subsample.id, "tsv", sep = "."
            )
          } else { # with imputations
            filename.assignment.res <- stri_join(
              "assignment", sampling.method, "imputed", "results", "individuals", 
              "iterations", "subsample", subsample.id, "tsv", sep = "."
            )
          }
        }
        write_tsv(x = assignment.res, 
                  path = paste0(directory.subsample,filename.assignment.res), 
                  col_names = TRUE, 
                  append = FALSE
        )
      } else { # with adegenet
        if (is.null(subsample)) {
          if (is.null(imputation.method)) {
            filename.assignment.res <- stri_join(
              "assignment", sampling.method, "no.imputation", "results", 
              "iterations", "tsv", sep = "."
            )
          } else { # with imputations
            filename.assignment.res <- stri_join(
              "assignment", sampling.method, "imputed", "results", "iterations", 
              "tsv", sep = "."
            )
          }
        } else {# with subsampling
          if (is.null(imputation.method)) {
            filename.assignment.res <- stri_join(
              "assignment", sampling.method, "no.imputation", "results", 
              "iterations", "subsample", subsample.id, "tsv", sep = "."
            )
          } else { # with imputations
            filename.assignment.res <- stri_join(
              "assignment", sampling.method, "imputed", "results", "iterations",
              "subsample", subsample.id, "tsv", sep = "."
            )
          }
        }
        write_tsv(x = assignment.res, 
                  path = paste0(directory.subsample,filename.assignment.res), 
                  col_names = TRUE, append = FALSE
        )
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
              CURRENT = factor(CURRENT, levels = pop.levels, ordered = T),
              CURRENT = droplevels(CURRENT)
            ) %>%
            group_by(CURRENT, MARKER_NUMBER, MISSING_DATA, METHOD) %>%
            summarise(
              MEAN = round(mean(MEAN_i), 2),
              SE = round(sqrt(stats::var(MEAN_i)/length(MEAN_i)), 2),
              MIN = round(min(MEAN_i), 2),
              MAX = round(max(MEAN_i), 2),
              MEDIAN = round(stats::median(MEAN_i), 2),
              QUANTILE25 = round(stats::quantile(MEAN_i, 0.25), 2),
              QUANTILE75 = round(stats::quantile(MEAN_i, 0.75), 2)
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
              SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
              MIN = round(min(ASSIGNMENT_PERC), 2),
              MAX = round(max(ASSIGNMENT_PERC), 2),
              MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
              QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
              QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
            ) %>%
            ungroup %>% 
            mutate(
              CURRENT = factor(CURRENT, levels = pop.levels, ordered = TRUE),
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
          SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
          MIN = round(min(ASSIGNMENT_PERC), 2),
          MAX = round(max(ASSIGNMENT_PERC), 2),
          MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
          QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
          QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
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
        if (is.null(imputation.method)) {
          filename.assignment.sum <- stri_join(
            "assignment", sampling.method, "no.imputation", "results", 
            "summary.stats", "tsv", sep = "."
          )
        } else { # with imputations
          filename.assignment.sum <- stri_join(
            "assignment", sampling.method, "imputed", "results", "summary.stats",
            "tsv", sep = "."
          )
        }
      } else {# with subsampling
        if (is.null(imputation.method)) {
          filename.assignment.sum <- stri_join(
            "assignment", sampling.method, "no.imputation", "results", 
            "summary.stats", "subsample", subsample.id, "tsv", sep = "."
          )
        } else { # with imputations
          filename.assignment.sum <- stri_join(
            "assignment", sampling.method, "imputed", "results", "summary.stats",
            "subsample", subsample.id, "tsv", sep = "."
          )
        }
      }
      write_tsv(x = assignment.summary.stats, 
                path = paste0(directory.subsample,filename.assignment.sum), 
                col_names = TRUE, append = FALSE
      )
      
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
        # i <- 3
        i <- iterations.list
        
        # Ranking Fst with training dataset (keep holdout individuals out)
        message("Ranking markers based on Fst")
        # THL = "all" ---------
        if (thl == "all") {
          holdout <- NULL
          fst.ranked <- assigner::fst_WC84(
            data = input, 
            pop.levels = pop.levels, pop.labels = pop.labels, strata = NULL, 
            holdout.samples = NULL
          )$fst.ranked
          if (!is.null(imputation.method)) {
            fst.ranked.imp <- assigner::fst_WC84(
              data = input.imp, 
              pop.levels = pop.levels, pop.labels = pop.labels, strata = NULL, 
              holdout.samples = NULL
            )$fst.ranked
          }
          ## THL = 1------------
        } else if (thl == 1) {
          holdout <- data.frame(INDIVIDUALS = i)
          fst.ranked <- assigner::fst_WC84(
            data = input, 
            pop.levels = pop.levels, pop.labels = pop.labels, strata = NULL, 
            holdout.samples = holdout$INDIVIDUALS
          )$fst.ranked
          if (!is.null(imputation.method)) {
            fst.ranked.imp <- assigner::fst_WC84(
              data = input.imp, 
              pop.levels = pop.levels, pop.labels = pop.labels, strata = NULL, 
              holdout.samples = holdout$INDIVIDUALS
            )$fst.ranked
          }
        } else { # thl proportion or > 1 ----------
          holdout <- data.frame(holdout.individuals.list[i])
          fst.ranked <- assigner::fst_WC84(
            data = input, 
            pop.levels = pop.levels, pop.labels = pop.labels, strata = NULL, 
            holdout.samples = holdout$INDIVIDUALS
          )$fst.ranked
          if (!is.null(imputation.method)) {
            fst.ranked.imp <- assigner::fst_WC84(
              data = input.imp, 
              pop.levels = pop.levels, pop.labels = pop.labels, strata = NULL, 
              holdout.samples = holdout$INDIVIDUALS
            )$fst.ranked
          }
        }
        
        fst.ranked.filename <- stri_join("fst.ranked_", i, ".tsv", sep = "") # No imputation
        write_tsv(
          x = fst.ranked, 
          path = paste0(directory.subsample, fst.ranked.filename), 
          col_names = TRUE, 
          append = FALSE
        )
        
        if (!is.null(imputation.method)) {  # With imputations
          fst.ranked.filename.imp <- stri_join(
            "fst.ranked_", i, "_imputed",".tsv", sep = ""
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
          # m <- 100 # test
          # m <- 150 # test
          # m <- 200 # test
          # m <- 400 # test
          m <- as.numeric(m)
          RANKING <- NULL
          select.markers <- filter(.data = fst.ranked, RANKING <= m) %>%
            select(MARKERS)
          
          # get the list of markers after filter
          markers.names <- unique(select.markers$MARKERS)
          
          # Assignment analysis without imputations
          filename.no.imp <- stri_paste(directory.subsample, filename, sep = "")
          filename.no.imp <- stri_replace_all_fixed(
            filename.no.imp, pattern = "txt",
            replacement = stri_join(
              i, m, "no.imputation", "txt", sep = "."
            )
          )
          if (assignment.analysis == "gsi_sim") {
            assignment.no.imp <- assignment_analysis(
              data = input,
              select.markers = select.markers,
              markers.names = markers.names,
              missing.data = "no.imputation", 
              i = i, 
              m = m,
              holdout = holdout,
              filename = filename.no.imp
            )
          }
          
          if (assignment.analysis == "adegenet") {
            assignment.no.imp <- assignment_analysis_adegenet(
              data = genind.object,
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
          
          # With imputations
          if (!is.null(imputation.method)) {  # with imputations
            
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
            filename.imp <- stri_paste(directory.subsample, filename, sep = "")
            filename.imp <- stri_replace_all_fixed(
              filename.imp, pattern = "txt", replacement = stri_join(
                i, m, "imputed", "txt", sep = "."
              )
            )
            if (assignment.analysis == "gsi_sim") {
              assignment.imp <- assignment_analysis(
                data = input.imp,
                select.markers = select.markers,
                markers.names = markers.names,
                missing.data = missing.data, 
                i = i,
                m = m,
                holdout = holdout,
                filename = filename.imp
              )
            }
            if (assignment.analysis == "adegenet") {
              assignment.imp <- assignment_analysis_adegenet(
                data = genind.object.imp,
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
          if (is.null(imputation.method)) {# with imputations
            assignment <- assignment.no.imp
            fst.ranked.imp <- NULL
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
          input.imp = input.imp,
          # vcf = vcf.imp, # was an error before, double check...
          pop.levels = pop.levels,
          pop.labels = pop.labels,
          sampling.method = sampling.method,
          thl = thl,
          iteration.method = iteration.method,
          filename = filename,
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
        if (is.null(imputation.method)) {
          filename.assignment.res <- stri_join("assignment", sampling.method, "no.imputation", "results", "individuals", "iterations", "tsv", sep = ".")
        } else { # with imputations
          filename.assignment.res <- stri_join("assignment", sampling.method, "imputed", "results", "individuals", "iterations", "tsv", sep = ".")
        }
      } else {# with subsampling
        if (is.null(imputation.method)) {
          filename.assignment.res <- stri_join("assignment", sampling.method, "no.imputation", "results", "individuals","iterations", "subsample", subsample.id, "tsv", sep = ".")
        } else { # with imputations
          filename.assignment.res <- stri_join("assignment", sampling.method, "imputed", "results", "individuals", "iterations", "subsample", subsample.id, "tsv", sep = ".")
        }
      }
      write_tsv(x = assignment.res.summary, path = paste0(directory.subsample,filename.assignment.res), col_names = TRUE, append = FALSE)
      
      
      if (thl == 1 | thl == "all") {
        assignment.stats.pop <- assignment.res.summary %>%
          mutate(
            CURRENT = factor(CURRENT, levels = pop.levels, ordered = TRUE),
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
            SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
            MIN = round(min(ASSIGNMENT_PERC), 2),
            MAX = round(max(ASSIGNMENT_PERC), 2),
            MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
            QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
            QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
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
        # if (assignment.analysis == "adegenet") {
        #   assignment.res.summary.prep <- assignment.res.summary %>% 
        #     # assignment.res.summary <- assignment.res.summary %>% 
        #     group_by(CURRENT, INFERRED, ITERATIONS, MARKER_NUMBER, MISSING_DATA, METHOD) %>% 
        #     tally %>% 
        #     group_by(CURRENT) %>% 
        #     mutate(TOTAL = sum(n)) %>% 
        #     ungroup() %>% 
        #     mutate(ASSIGNMENT_PERC = round(n/TOTAL*100, 0)) %>% 
        #     filter(CURRENT == INFERRED) %>% 
        #     select(-n, -TOTAL)
        # }
        
        # if (assignment.analysis == "gsi_sim") {
        assignment.res.summary.prep <- assignment.res.summary %>% 
          group_by(CURRENT, MARKER_NUMBER, METHOD, MISSING_DATA, ITERATIONS) %>%
          summarise(
            n = length(CURRENT[as.character(CURRENT) == as.character(INFERRED)]),
            TOTAL = length(CURRENT)
          ) %>%
          ungroup() %>% 
          mutate(ASSIGNMENT_PERC = round(n/TOTAL*100, 0)) %>% 
          select(-n, -TOTAL)
        # }
        
        if (is.null(subsample)) {
          if (is.null(imputation.method)) {
            filename.assignment.res.sum <- stri_join("assignment", sampling.method, "no.imputation", "results", "summary", "tsv", sep = ".")
          } else { # with imputations
            filename.assignment.res.sum <- stri_join("assignment", sampling.method, "imputed", "results", "summary", "tsv", sep = ".")
          }
        } else {# with subsampling
          if (is.null(imputation.method)) {
            filename.assignment.res.sum <- stri_join("assignment", sampling.method, "no.imputation", "results", "summary", "subsample", subsample.id, "tsv", sep = ".")
          } else { # with imputations
            filename.assignment.res.sum <- stri_join("assignment", sampling.method, "imputed", "results", "summary", "subsample", subsample.id, "tsv", sep = ".")
          }
        }
        write_tsv(x = assignment.res.summary.prep, path = paste0(directory.subsample,filename.assignment.res.sum), col_names = TRUE, append = FALSE)
        
        assignment.stats.pop <- assignment.res.summary.prep %>%
          # assignment.stats.pop <- assignment.res.summary %>%
          mutate(
            CURRENT = factor(CURRENT, levels = pop.levels, ordered = TRUE),
            CURRENT = droplevels(CURRENT)
          ) %>%
          group_by(CURRENT, MARKER_NUMBER, METHOD, MISSING_DATA) %>%
          summarise(
            MEAN = round(mean(ASSIGNMENT_PERC), 2),
            SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
            MIN = round(min(ASSIGNMENT_PERC), 2),
            MAX = round(max(ASSIGNMENT_PERC), 2),
            MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
            QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
            QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2),
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
            SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
            MIN = round(min(ASSIGNMENT_PERC), 2),
            MAX = round(max(ASSIGNMENT_PERC), 2),
            MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
            QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
            QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2),
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
        if (is.null(imputation.method)) {
          filename.assignment.sum <- stri_join("assignment", sampling.method, "no.imputation", "results", "summary.stats", "tsv", sep = ".")
        } else { # with imputations
          filename.assignment.sum <- stri_join("assignment", sampling.method, "imputed", "results", "summary.stats", "tsv", sep = ".")
        }
      } else {# with subsampling
        if (is.null(imputation.method)) {
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
             sampling.method = sampling.method,
             thl = thl,
             iteration.method = iteration.method,
             filename = filename,
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
  
  if (is.null(imputation.method)) {
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
        SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
        MIN = round(min(ASSIGNMENT_PERC), 2),
        MAX = round(max(ASSIGNMENT_PERC), 2),
        MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
        QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
        QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
      ) %>% 
      mutate(SUBSAMPLE = rep("OVERALL", n())) %>%
      arrange(CURRENT, MARKER_NUMBER)
    
    res.overall <- res.pop %>% 
      group_by(MARKER_NUMBER, MISSING_DATA, METHOD) %>%
      rename(ASSIGNMENT_PERC = MEAN) %>%
      summarise(
        MEAN = round(mean(ASSIGNMENT_PERC), 2),
        SE = round(sqrt(stats::var(ASSIGNMENT_PERC)/length(ASSIGNMENT_PERC)), 2),
        MIN = round(min(ASSIGNMENT_PERC), 2),
        MAX = round(max(ASSIGNMENT_PERC), 2),
        MEDIAN = round(stats::median(ASSIGNMENT_PERC), 2),
        QUANTILE25 = round(stats::quantile(ASSIGNMENT_PERC, 0.25), 2),
        QUANTILE75 = round(stats::quantile(ASSIGNMENT_PERC, 0.75), 2)
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
    
    if (is.null(imputation.method)) {
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
  if (is.null(imputation.method)) { # no imputation
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
