# Reproducible stockR hard-genotype comparison----------------------------------
# Run from assigner's source directory. No source datasets are modified.
# Outputs are printed; all intermediate VCF/GDS files are temporary.
library(stockR)
source("R/assign_likelihood_dapc.R")
source("R/stock_likelihood_mixture.R")
source("R/stock_mixture_blocks.R")
source("R/stock_mixture_gds.R")
source("R/stock_mixture_validation.R")

run_stockr_validation <- function() {
  genotype <- get("myData", envir = asNamespace("stockR"))
  truth <- attr(genotype, "grps")
  groups <- get("stocks", envir = asNamespace("stockR"))$init.grps
  reference <- stockR::stockSTRUCTURE(genotype, K = 3, start.grps = groups,
    control = list(method = "EM", quiet = TRUE, EM.maxit = 500,
                   EM.eps = 1e-9, EM.tau.step = 1))
  g <- reference$data$SNPdata
  tau <- stockR:::get.init.taus(groups, magic.num = reference$control$EM.tau.init.max, K = 3)
  initial <- stockR:::calcCondMeans(reference$data, tau, K = 3)
  rownames(initial) <- as.character(seq_len(nrow(g)))
  colnames(initial) <- paste0("STOCK_", 1:3)
  folder <- tempfile("stockr_validation_"); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  vcf <- file.path(folder, "example.vcf")
  header <- c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000000>",
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Near-certain hard-genotype likelihood">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", paste0("s", seq_len(ncol(g)))), collapse = "\t"))
  records <- vapply(seq_len(nrow(g)), function(l) {
    value <- c("0,999,999", "999,0,999", "999,999,0")[g[l, ] + 1L]
    value[is.na(value)] <- "."
    paste(c("1", l, paste0("v", l), "A", "C", ".", "PASS", ".", "PL", value),
          collapse = "\t")
  }, character(1))
  writeLines(c(header, records), vcf)
  path <- file.path(folder, "example.gds")
  SeqArray::seqVCF2GDS(vcf, path, verbose = FALSE)
  elapsed <- system.time(fit <- .sm_fit_gds(path, initial,
    proportion = unname(stockR:::calcPis(tau)), n.starts = 1, block.snps = 500,
    max.iter = 500, tolerance = 1e-10, epsilon = 1 / (10 * ncol(g))))
  # Compare partitions independently of arbitrary component numbering.
  ours <- max.col(fit$membership)
  theirs <- max.col(reference$postProbs)
  pair <- function(z) outer(z, z, "==")
  effective <- (!is.na(g)) %*% fit$membership
  identified <- effective > 1e-10
  frequency.difference <- abs(fit$frequency -
      pmin(pmax(reference$condMeans, fit$epsilon), 1 - fit$epsilon))
  report <- list(stockR.version = as.character(packageVersion("stockR")),
    samples = ncol(g), input.markers = nrow(genotype), retained.markers = nrow(g),
    converged = fit$converged, iterations = fit$iterations,
    elapsed.seconds = unname(elapsed["elapsed"]),
    stockR.loglik = reference$margLogL, GLPL.loglik = fit$loglik,
    max.membership.difference = max(abs(fit$membership - reference$postProbs)),
    uninformative.stock.marker.cells = sum(!identified),
    max.frequency.difference.after.floor = max(frequency.difference[identified]),
    max.unidentified.frequency.difference = max(c(0, frequency.difference[!identified])),
    same.stockR.partition = identical(pair(ours), pair(theirs)),
    same.generating.partition = all(pair(ours) == pair(truth)))
  print(report)
  stopifnot(report$converged, report$same.stockR.partition,
            report$max.membership.difference < 1e-5,
            report$max.frequency.difference.after.floor < 1e-5)
  if ("--extended" %in% commandArgs(trailingOnly = TRUE)) {
    # Modest subset for workflow verification, not a coverage study. The full
    # 5000-marker matched-initialization comparison above remains separate.
    gds <- SeqArray::seqOpen(path)
    on.exit(SeqArray::seqClose(gds), add = TRUE, after = FALSE)
    SeqArray::seqSetFilter(gds, variant.id = 1:300, verbose = FALSE)
    subset.fit <- .sm_fit_gds(gds, initial[1:300, ], n.starts = 3, seed = 51)
    cat("300-marker fitted partition matches generating groups:",
        all(pair(max.col(subset.fit$membership)) == pair(truth)), "\n")
    bootstrap <- .sm_bootstrap_gds(gds, subset.fit, B = 8, n.starts = 2, seed = 53)
    print(bootstrap$diagnostics)
    folds <- stats::setNames(rep(1:3, length.out = ncol(g)), paste0("s", seq_len(ncol(g))))
    cv <- .sm_cv_gds(gds, K = 1:3, folds = folds, n.starts = 3,
                     max.iter = 200, seed = 54)
    print(cv$summary)
    print(cv$diagnostics)
    cat("Best predictive K in the tested grid:", cv$best.predictive.K, "\n")
    if ("--sensitivity" %in% commandArgs(trailingOnly = TRUE)) {
      for (bound in c(.001, .01)) {
        alternate <- .sm_cv_gds(gds, K = 1:3, folds = folds, n.starts = 3,
          max.iter = 200, seed = 54, epsilon = bound)
        cat("Frequency-bound sensitivity, epsilon =", bound, "\n")
        print(alternate$summary)
      }
    }
  }
  invisible(report)
}
run_stockr_validation()
