# Run from this directory. Reads completed experiments; never runs simulations.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, dir.exists(args[[1]]))
root <- normalizePath(args[[1]])
files <- c(coverage = "coverage-sensitivity-20260917/summary.rds",
  calls = "gl-gt-comparison-20260917/summary.rds",
  technical = "technical-structure-20260918/report-data.rds",
  oracle = "weak-mixture-diagnostic-20260918/diagnostic-summary.rds")
paths <- file.path(root, files)
stopifnot(all(file.exists(paths)))
data <- lapply(paths, readRDS)
names(data) <- names(files)
data$provenance <- data.frame(source = files, md5 = unname(tools::md5sum(paths)))
dir.create("data", showWarnings = FALSE)
saveRDS(data, "data/review-results.rds")
sources <- c("stock_likelihood_mixture.R", "stock_mixture_blocks.R",
             "stock_mixture_gds.R", "stock_mixture_validation.R",
             "assign_likelihood_dapc.R")
source.paths <- file.path("../../../R", sources)
stopifnot(all(file.exists(source.paths)))
dir.create("source-snapshot", showWarnings = FALSE)
stopifnot(all(file.copy(source.paths, file.path("source-snapshot", sources), overwrite = TRUE)))
saveRDS(data.frame(file = sources, md5 = unname(tools::md5sum(source.paths))),
        "data/source-checksums.rds")
writeLines(capture.output(sessionInfo()), "data/preparation-session.txt")
