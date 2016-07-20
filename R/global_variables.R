# remove NOTE about no visible binding for global variable during R CMD check --
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
      "COUNT", "MAX_COUNT_MARKERS", "hierarchy", "GT_VCF", "ANALYSIS",
      "MEAN_ITERATIONS", "MEAN_SUBSAMPLE", "NUMBER_ITERATIONS",
      "NUMBER_SUBSAMPLE", "TOTAL_ITERATIONS", "TOTAL_SUBSAMPLE", "X1", "X2",
      "strata.df")
  )
}
