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

#' @param data 6 options: vcf (to make vcf population ready, see details below),
#' plink, stacks haplotype file, genind, genepop, 
#' and a data frame in wide format. \emph{See details}.

#' @param strata (optional for data frame and PLINK files, 
#' required for VCF and haplotypes files) A tab delimited file with 2 columns with header:
#' \code{INDIVIDUALS} and \code{STRATA}. With a 
#' data frame of genotypes the strata is the INDIVIDUALS and POP_ID columns, with
#' PLINK files, the \code{tfam} first 2 columns are used. 
#' If a \code{strata} file is specified, the strata file will have
#' precedence. The \code{STRATA} column can be any hierarchical grouping. 
#' To create a strata file see \code{\link[stackr]{individuals2strata}}.
#' If you have already run 
#' \href{http://catchenlab.life.illinois.edu/stacks/}{stacks} on your data, 
#' the strata file is similar to a stacks `population map file`, make sure you 
#' have the required column names  (\code{INDIVIDUALS} and \code{STRATA}).
#' Default: \code{strata = NULL}.

#' @param mixture (optional) A file in the working directory (e.g. "mixture.txt")
#' with mixture individual ID. The column header is \code{INDIVIDUALS}.
#' Default: \code{mixture = NULL} but see details on the different ways to 
#' identify your unknown or mixture samples.
#' 
#' @param assignment.analysis Assignment analysis conducted with 
#' \code{assignment.analysis = "gsi_sim"} or 
#' \code{assignment.analysis = "adegenet"}.

#' @inheritParams stackr::tidy_genomic_data 

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

#' @param folder (optional) The name of the folder created in the working 
#' directory to save the files/results. Default: \code{folder = NULL}.

#' @param filename (optional) The name of the file written to the directory.
#' Use the extension ".txt" at the end. 
#' Default \code{filename = assignment_data.txt}.
#' The number of markers used will be appended to the name of the file.

#' @param keep.gsi.files (Logical) With default: \code{keep.gsi.files = FALSE}, 
#' the input and output gsi_sim files
#' will be deleted from the directory when finished processing.
#' With \code{keep.gsi.files = TRUE}, remember to allocate a large chunk of the disk space for the analysis.

#' @param subsample (Integer or Proportion) Default is no sumsampling, 
#' \code{subsample = NULL}.
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
#' e.g. To test 500, 1000, 2000 and all the markers:
#' \code{marker.number = c(500, 1000, 2000, "all")}.
#' To use only 500 makers \code{marker.number = 500}.


#' @inheritParams stackr::stackr_imputations_module 

#' @param impute.mixture (Logical) Imputations of mixture samples.
#' Default: \code{impute.mixture = FALSE}. For no imputation. 
#' For \code{impute.mixture = TRUE} the imputations.group (see below)
#' for the mixture samples is automatically set to 
#' \code{imputations.group = "global"}. Warning: bias could be introduced by imputing
#' missing genotype in the mixture samples.

#' @details
#' 
#' \strong{2 options to identify your unknown or mixture samples, using:}
#' \itemize{
#'   \item \strong{strata} argument with "mixture" or "unknown" 
#'   instead of the grouping id (e.g. populations) for the other samples.
#'    \item \strong{mixture} argument using a file in the working directory 
#'    (e.g. "mixture.txt") with mixture individual ID. 
#'    The column header is \code{INDIVIDUALS}.
#'   }
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
#' which will do its best to put a usable copy of gsi_sim where it is needed. To do 
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
#' strata = "strata.tsv",
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
#' imputation.method = NULL,
#' parallel.core = 12
#' )
#' # with gsi_sim for the mixture assignment and sampling.method = "ranked"
#' # Here I also want to impute the genotypes of the data (baseline and mixture) 
#' # using random forest:
#' assignment.tuna <- assignment_mixture(
#' data = "data.frame.genotypes.tuna.tsv",
#' mixture = "cohort.tuna.tsv",
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
#' Regression and Classification (RF-SRC), R package version 1.6.1.
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


assignment_mixture <- function(
  data,
  strata = NULL,
  mixture = NULL,
  assignment.analysis,
  sampling.method,
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
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  imputation.method = NULL,
  impute.mixture = FALSE,
  impute = "genotype",
  imputations.group = "populations",
  num.tree = 100,
  iteration.rf = 10,
  split.number = 100,
  verbose = FALSE,
  folder = NULL,
  filename = "assignment_data.txt",
  keep.gsi.files = FALSE,
  parallel.core = detectCores()-1
) {
  cat("#######################################################################\n")
  cat("#################### assigner::assignment_mixture #####################\n")
  cat("#######################################################################\n")
  
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
  
  # POP_ID in gsi_sim does not like spaces, we need to remove space in everything touching POP_ID...
  # pop.levels, pop.labels, pop.select, strata, etc
  if (!is.null(pop.levels) & is.null(pop.labels)) {
    pop.levels <- stri_replace_all_fixed(pop.levels, pattern = " ", replacement = "_", vectorize_all = FALSE)
    pop.labels <- pop.levels
  }
  if (!is.null(pop.labels)) {
    pop.labels <- stri_replace_all_fixed(pop.labels, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  if (!is.null(pop.labels) & is.null(pop.levels)) stop("pop.levels is required if you use pop.labels")
  if (!is.null(pop.select)) {
    pop.select <- stri_replace_all_fixed(pop.select, pattern = " ", replacement = "_", vectorize_all = FALSE)
  }
  
  # File type detection----------------------------------------------------------
  if(is.genind(data)){
    data.type <- "genind.file"
    # message("File type: genind object")
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
    if (stri_detect_fixed(str = data.type, pattern = "POP_ID") | stri_detect_fixed(str = data.type, pattern = "INDIVIDUALS") | stri_detect_fixed(str = data.type, pattern = "MARKERS")| stri_detect_fixed(str = data.type, pattern = "LOCUS")) {
      data.type <- "df.file"
      # message("File type: data frame of genotypes")
    }
    
    
    if (stri_detect_fixed(str = data.type, pattern = "Catalog")) {
      data.type <- "haplo.file"
      # message("File type: haplotypes from stacks")
      # if (is.null(blacklist.genotype)) {
      #   stop("blacklist.genotype file missing. 
      #        Use stackr's missing_genotypes function to create this blacklist")
      # }
    }
    if (stri_detect_fixed(str = data, pattern = ".gen")) {
      data.type <- "genepop.file"
      # message("File type: genepop")
    } 
  } # end file type detection
  
  if(data.type == "haplo.file") {
    message("With stacks haplotype file the maf.approach is automatically set to: haplotype")
    maf.approach <- "SNP"
    # confusing, but because the haplotpe file doesn't have snp info, only locus info
    # it's treated as markers/snp info and filtered the same way as the approach by SNP.
    # but it's really by haplotype
  }
  
  if (maf.approach == "haplotype") {
    if (data.type != "vcf.file" | data.type != "haplo.file") {
      stop("The haplotype approach during MAF filtering is for VCF and
           stacks haplotypes file, only. Use the snp approach for the other file types")
    }
  }
  
  # Strata argument required for VCF and haplotypes files-----------------------
  if (data.type == "haplo.file" | data.type == "vcf.file") {
    if (is.null(strata)) stop("strata argument is required")
  }
  
  # Base filename
  base.filename <- filename # filename will change from time to time in the function
  
  # Create a folder based on filename to save the output files *****************
  if (is.null(folder)) {
    # Get date and time to have unique filenaming
    file.date <- stri_replace_all_fixed(Sys.time(), pattern = " EDT", replacement = "")
    file.date <- stri_replace_all_fixed(file.date, pattern = c("-", " ", ":"), replacement = c("", "@", ""), vectorize_all = FALSE)
    file.date <- stri_sub(file.date, from = 1, to = 13)
    
    if (is.null(imputation.method)) {
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
  
  # change "unknown" to "mixture" for simplicity of pop_id below ---------------
  input$POP_ID <- stri_replace_all_fixed(
    str = input$POP_ID, 
    pattern = "unknnown", 
    replacement = "mixture", 
    vectorize_all = FALSE
  )
  
  # mixture data  **************************************************************
  if ("mixture" %in% unique(input$POP_ID)) {
    mixture.df <- input %>% 
      filter(POP_ID == "mixture") %>% 
      distinct(INDIVIDUALS)
  }
  
  if (!is.null(mixture)) {
    mixture.df <- read_tsv(file = mixture, col_names = TRUE, col_types = "c")
    input <- mutate(.data = input, POP_ID = if_else(INDIVIDUALS %in% mixture.df$INDIVIDUALS, "mixture", as.character(POP_ID)))
  }
  
  # subsampling data ***********************************************************
  # Function:
  subsampling_data <- function(iteration.subsample, ...) {
    if (is.null(subsample)) {
      subsample.select <- ind.pop.df %>% 
        mutate(SUBSAMPLE = rep(iteration.subsample, n()))
    } else {
      # separate all the mixture samples
      mixture.select <- ind.pop.df %>% filter(POP_ID == "mixture")
      
      # subsample the baseline
      if (subsample > 1) { # integer
        subsample.select <- ind.pop.df %>%
          filter(POP_ID != "mixture") %>% 
          group_by(POP_ID) %>%
          sample_n(subsample, replace = FALSE) %>% # sampling individuals for each pop
          arrange(POP_ID, INDIVIDUALS)
      }
      if (subsample < 1) { # proportion
        subsample.select <- ind.pop.df %>%
          filter(POP_ID != "mixture") %>% 
          group_by(POP_ID) %>%
          sample_frac(subsample, replace = FALSE) %>% # sampling individuals for each pop
          arrange(POP_ID, INDIVIDUALS)
      }
      # Join baseline and mixture back in 1 dataset
      subsample.select <- bind_rows(subsample.select, mixture.select) %>% 
        mutate(SUBSAMPLE = rep(iteration.subsample, n())) %>% 
        ungroup()
    }
    return(subsample.select)
  } # End subsampling function
  
  # create the subsampling list
  ind.pop.df <- input %>% distinct(POP_ID, INDIVIDUALS)
  subsample.list <- map(.x = 1:iteration.subsample, .f = subsampling_data, subsample = subsample)
  
  # keep track of subsampling individuals and write to directory
  if (is.null(subsample)) {
    message("Subsampling: not selected")
  } else {
    message("Subsampling: selected")
    subsampling.individuals <- bind_rows(subsample.list)
    write_tsv(
      x = subsampling.individuals, 
      path = paste0(directory, "subsampling_individuals.tsv"), 
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
    
    # LD control... keep only 1 SNP per haplotypes/reads (optional) ------------
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
    
    
    # Markers in common between all populations (optional) ---------------------
    # This need to be moved while doing the assignment
    if (common.markers) { # keep only markers present in all pop
      message("Using markers common in all baseline populations:")
      pop.number <- n_distinct(input$POP_ID[input$POP_ID != "mixture"])
      
      pop.filter <- input %>% 
        filter(POP_ID != "mixture") %>% 
        filter(GT != "000000") %>% 
        group_by(MARKERS) %>%
        filter(n_distinct(POP_ID) == pop.number) %>%
        arrange(MARKERS) %>%
        ungroup() %>% 
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
          filter(POP_ID != "mixture") %>% 
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
            filter(POP_ID != "mixture") %>% 
            tidyr::separate(data = ., col = GT, into = c("A1", "A2"), sep = 3, remove = TRUE) %>% 
            tidyr::gather(data = ., key = ALLELES, value = GT, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
            select(MARKERS, GT, POP_ID) %>% 
            filter(GT != "000")
        }
        
        if (data.type == "plink.file") { # For PLINK and common code below
          maf.data <- input %>%
            filter(POP_ID != "mixture") %>% 
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

    # Keep a strata df ---------------------------------------------------------
    strata.df <- distinct(.data = input, INDIVIDUALS, POP_ID)
    
    # Adegenet no imputations --------------------------------------------------
    if (assignment.analysis == "adegenet" ) {
      genind.object <- stackr::write_genind(data = input)
    }
    
    input.baseline <- filter(.data = input, POP_ID != "mixture")
    input.mixture <- filter(.data = input, POP_ID == "mixture")
    
    # Imputations --------------------------------------------------------------
    if (!is.null(imputation.method)) {
      message("Preparing the data for imputations")
      
      input.baseline <- filter(.data = input, POP_ID != "mixture")
      input.mixture <- filter(.data = input, POP_ID == "mixture")
      # strata.df.impute <- distinct(.data = input, INDIVIDUALS, POP_ID)
      
      # imputation for the mixture samples, if selected, is always conducted globally. 
      message("Imputations of baseline samples ...")
      input.baseline.imp <- stackr::stackr_imputations_module(
        data = input.baseline, 
        imputation.method = imputation.method, 
        impute = impute, 
        imputations.group = imputations.group, 
        num.tree = num.tree, 
        iteration.rf = iteration.rf, 
        split.number = split.number, 
        verbose = verbose, 
        parallel.core = parallel.core, 
        filename = NULL
      )
      
      # combine the mixture (no imputation) + the imputed baseline
      input.imp <- suppressWarnings(
        bind_rows(input.baseline.imp, input.mixture)
      )
      
      if (impute.mixture) {
        # impute globally the mixture samples
        message("Imputations computed globally for mixture samples:")
        input.imp <- stackr::stackr_imputations_module(
          data = input.imp, 
          imputation.method = imputation.method, 
          impute = impute, 
          imputations.group = "global", 
          num.tree = num.tree, 
          iteration.rf = iteration.rf, 
          split.number = split.number, 
          verbose = verbose, 
          parallel.core = parallel.core, 
          filename = NULL
        )
              } # End impute.mixture
      
      input.baseline.imp <- filter(.data = input.imp, POP_ID != "mixture")
      input.mixture.imp <- filter(.data = input.imp, POP_ID == "mixture")
      
      # test <- input.imp %>% filter(GT == "000000")
      # prep. adegenet
      if (assignment.analysis == "adegenet") {
        genind.object.imp <- stackr::write_genind(data = input.imp)
      } # end adegenet
    } # end imputations

  # Sampling of markers---------------------------------------------------------
  # unique list of markers after all the filtering
  # if "all" is present in the list, change to the maximum number of markers
  unique.markers <- input %>% 
    distinct(MARKERS) %>% 
    arrange(MARKERS)
  
  marker.number <- stri_replace_all_fixed(
    str = marker.number, 
    pattern = "all", 
    replacement = nrow(unique.markers), 
    vectorize_all = TRUE
  )
  marker.number <- as.numeric(marker.number)
  
  # In marker.number, remove marker group higher than the max number of markers
  removing.marker <- purrr::keep(.x = marker.number, .p = marker.number > nrow(unique.markers))
  
  if (length(removing.marker) > 0) {
    message(
      "Removing marker.number higher than the max number of markers: ", 
      stri_c(removing.marker, collapse = ", ")
    )
  }
  marker.number <- purrr::discard(.x = marker.number, .p = marker.number > nrow(unique.markers))
  
  
  # Functions ------------------------------------------------------------------
  # Assignment with gsi_sim
  assignment_analysis <- function(data, select.markers, markers.names, missing.data, i, m, filename, ...) {
    # data <- input #test
    # missing.data <- "no.imputation" #test
    data.select <- suppressWarnings(
      data %>%
        semi_join(select.markers, by = "MARKERS") %>%
        arrange(POP_ID, INDIVIDUALS, MARKERS)
    )
    
    input.baseline <- filter(.data = data, POP_ID != "mixture")
    input.mixture <- filter(.data = data, POP_ID == "mixture")
    
    # Baseline filename
    baseline.filename <- stri_replace_all_fixed(
      filename, 
      pattern = ".txt",
      replacement = "_baseline.txt", 
      vectorize_all = FALSE
      )
    baseline.input.gsi <- stri_join(directory.subsample, baseline.filename)
    
    # save input file to directory
    baseline.input <- assigner::write_gsi_sim(
      data = input.baseline, 
      markers.names = markers.names, 
      # directory.subsample = directory.subsample,
      # i = i,
      # m = m,
      filename = baseline.input.gsi
      )
    
    # Mixture filename
    mixture.filename <- stri_replace_all_fixed(
      filename, 
      pattern = ".txt",
      replacement = "_mixture.txt",
      vectorize_all = FALSE
    )
    mixture.input.gsi <- stri_join(directory.subsample, mixture.filename)
    
    # save input file to directory
    mixture.input <- assigner::write_gsi_sim(
      data = input.mixture, 
      markers.names = markers.names, 
      # directory.subsample = directory.subsample, 
      # i = i, 
      # m = m,
      filename = mixture.input.gsi
      )
    
    # Run gsi_sim ------------------------------------------------------------
    output.gsi <- stri_replace_all_fixed(mixture.input, pattern = "txt", replacement = "output.txt")
    # output.gsi <- stri_join(directory.subsample, output.gsi)
    setwd(directory.subsample)
    system(paste("gsi_sim -b ", baseline.input.gsi, "-t ", mixture.input.gsi, " > ", output.gsi))
    
    # Option remove the input file from directory to save space
    if (!keep.gsi.files) {
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
    
    
    # Option remove the output file from directory to save space
    if (!keep.gsi.files) file.remove(output.gsi)
    
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
    training.data <- data.select[data.select@strata$POP_ID != "mixture"]
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
    holdout.data <- data.select[data.select@strata$POP_ID == "mixture"]
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
      filename <- stri_replace_all_fixed(
        base.filename,
        pattern = ".txt",
        replacement = stri_join(
          "", "iteration", i, "markers", m, 
          "no_imputation.txt", sep = "_"
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
          filename = filename
        )
      }
      if (assignment.analysis == "adegenet") {
        assignment.no.imp <- assignment_analysis_adegenet(
          data = genind.object,
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
        filename <- stri_replace_all_fixed(
          base.filename,
          pattern = ".txt",
          replacement = stri_join(
            "", "iteration", i, "markers", m,
            "imputed.txt", sep = "_"
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
            filename = filename
          )
        }
        if (assignment.analysis == "adegenet") {
          assignment.imp <- assignment_analysis_adegenet(
            data = genind.object.imp,
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
    if (is.null(imputation.method)) {
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
    if (is.null(imputation.method)) {
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
    if (!is.null(imputation.method)) {
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
      filename <- stri_replace_all_fixed(
        base.filename,
        pattern = ".txt",
        replacement = stri_join(
          "", "markers", m, 
          "no_imputation.txt", sep = "_"
        )
      )
      
      
      if (assignment.analysis == "gsi_sim") {
        assignment.no.imp <- assignment_analysis(
          data = input,
          select.markers = select.markers,
          markers.names = markers.names,
          missing.data = "no.imputation",
          i = NULL,
          m = m,
          filename = filename
        )
      }
      
      if (assignment.analysis == "adegenet") {
        assignment.no.imp <- assignment_analysis_adegenet(
          data = genind.object,
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
        filename <- stri_replace_all_fixed(
          base.filename,
          pattern = ".txt",
          replacement = stri_join(
            "", "markers", m, 
            "imputed.txt", sep = "_"
          )
        )
        
        if (assignment.analysis == "gsi_sim") {
          assignment.imp <- assignment_analysis(
            data = input.imp,
            select.markers = select.markers,
            markers.names = markers.names,
            missing.data = missing.data, 
            i = NULL,
            m = m,
            filename = filename
          )
        }
        if (assignment.analysis == "adegenet") {
          assignment.imp <- assignment_analysis_adegenet(
            data = genind.object.imp,
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
      if (is.null(imputation.method)) {# with imputations
        assignment <- assignment.no.imp
        fst.ranked.imp <- NULL
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
      input.imp = input.imp,
      pop.levels = pop.levels,
      pop.labels = pop.labels,
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
    if (is.null(imputation.method)) {
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
if (is.null(imputation.method)) {
  filename.assignment.sum <- stri_join("assignment.mixture.summary.results", "no.imputation", sampling.method, "tsv", sep = ".")
} else { # with imputations
  filename.assignment.sum <- stri_join("assignment.mixture.summary.results", "imputed", sampling.method, "tsv", sep = ".")
}
write_tsv(x = assignment.mixture.summary.subsample, path = paste0(directory,filename.assignment.sum), col_names = TRUE, append = FALSE)

# results
res.list <- list(assignment = res, assignment.mixture.summary.results = assignment.mixture.summary.subsample)
cat("############################## completed ##############################\n")
return(res.list)
} # End assignment_mixture

