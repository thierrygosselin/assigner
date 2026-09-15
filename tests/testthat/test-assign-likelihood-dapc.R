likelihood_gds <- function(folder, partial = FALSE, field = "PL") {
  samples <- paste0("s", seq_len(16))
  variants <- paste0("v", seq_len(12))
  group <- rep(c("A", "B"), each = 8)
  header <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=100000>",
    if (field == "PL")
      '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">'
    else
      '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10 genotype likelihoods">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", samples), collapse = "\t")
  )
  records <- character(length(variants))
  for (l in seq_along(variants)) {
    values <- ifelse(group == "A", "0,40,80", "80,40,0")
    if (l %% 3L == 0L) values[c(2, 10)] <- "20,0,20"
    if (field == "GL") {
      values <- vapply(strsplit(values, ",", fixed = TRUE), function(x) {
        paste(format(-as.numeric(x) / 10, scientific = FALSE, trim = TRUE),
              collapse = ",")
      }, character(1))
    }
    if (l == 1L) values[16] <- "."
    if (partial && l == 2L) values[1] <- "0,.,20"
    records[l] <- paste(c("chr1", l * 100, variants[l], "A", "G", ".",
                           "PASS", ".", field, values), collapse = "\t")
  }
  vcf <- file.path(folder, "likelihoods.vcf")
  gds <- file.path(folder, "likelihoods.gds")
  writeLines(c(header, records), vcf)
  SeqArray::seqVCF2GDS(vcf, gds, ignore.chr.prefix = "", verbose = FALSE)
  list(path = gds, samples = samples,
       groups = stats::setNames(group, samples))
}

test_that("assign_likelihood_dapc fits PL data and projects held-out samples", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  fixture <- likelihood_gds(folder)
  training <- fixture$samples[c(1:6, 9:14)]
  fit <- assign_likelihood_dapc(
    fixture$path,
    groups = fixture$groups,
    likelihood = "PL",
    n.components = 1,
    training.samples = training,
    max.iter = 30,
    tolerance = 1e-4,
    verbose = FALSE
  )

  expect_s3_class(fit, "assign_likelihood_dapc")
  expect_identical(fit$likelihood, "PL")
  expect_equal(dim(fit$scores), c(16, 1))
  expect_equal(dim(fit$score.covariance), c(1, 1, 16))
  expect_true(all(is.finite(fit$scores)))
  expect_true(all(is.finite(fit$score.covariance)))
  expect_equal(unname(rowSums(fit$probability)), rep(1, 16), tolerance = 1e-8)
  expect_identical(fit$assignments$TRAINING, fixture$samples %in% training)
  expect_gte(mean(fit$assignments$ASSIGNED_GROUP == fixture$groups), 0.9)
  expect_equal(length(fit$diagnostics$rmse), fit$diagnostics$iterations)
  expect_true(all(fit$eigenvalues >= 0))
  expect_equal(fit$genomic.covariance,
               t(fit$genomic.covariance), tolerance = 1e-10)
})

test_that("assign_likelihood_dapc accepts GL and STRATA metadata", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  fixture <- likelihood_gds(folder, field = "GL")
  strata <- data.frame(
    INDIVIDUALS = fixture$samples,
    STRATA = unname(fixture$groups),
    PLATE = rep(c("P1", "P2"), 8)
  )
  fit <- assign_likelihood_dapc(fixture$path, groups = strata,
                         n.components = 1, max.iter = 30,
                         tolerance = 1e-4, verbose = FALSE)

  expect_identical(fit$likelihood, "GL")
  expect_identical(colnames(fit$scores), "PC1")
  expect_gte(mean(fit$assignments$ASSIGNED_GROUP == strata$STRATA), 0.9)
})

test_that("assign_likelihood_dapc validates likelihood vectors and training labels", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  fixture <- likelihood_gds(folder)
  bad.groups <- fixture$groups
  bad.groups[fixture$samples[1]] <- NA
  expect_error(
    assign_likelihood_dapc(fixture$path, bad.groups, n.components = 1,
                    training.samples = fixture$samples, max.iter = 1,
                    verbose = FALSE),
    "Every training sample"
  )

  other <- tempfile(); dir.create(other)
  on.exit(unlink(other, recursive = TRUE), add = TRUE)
  partial <- likelihood_gds(other, partial = TRUE)
  expect_error(
    assign_likelihood_dapc(partial$path, partial$groups, n.components = 1,
                    max.iter = 1, verbose = FALSE),
    "Partially missing PL"
  )
})

test_that("allele-frequency EM recovers ALT frequency from decisive calls", {
  dosage <- c(0, 1, 2, 2)
  log.likelihood <- array(-1000, dim = c(length(dosage), 3, 1))
  for (i in seq_along(dosage)) {
    log.likelihood[i, dosage[i] + 1L, 1] <- 0
  }

  observed <- .ld_allele_frequency(log.likelihood)
  expect_equal(observed, sum(dosage) / (2 * length(dosage)),
               tolerance = 1e-7)
})

test_that("allele-frequency EM agrees with direct likelihood optimization", {
  # Independent numerical check of the observed-data likelihood for a locus.
  likelihood <- rbind(
    c(0.90, 0.09, 0.01),
    c(0.15, 0.75, 0.10),
    c(0.02, 0.18, 0.80),
    c(0.25, 0.50, 0.25),
    c(0.70, 0.25, 0.05)
  )
  log.likelihood <- array(log(likelihood),
                          dim = c(nrow(likelihood), 3, 1))
  em <- .ld_allele_frequency(
    log.likelihood,
    max.iter = 1000,
    tolerance = 1e-10
  )

  negative.log.likelihood <- function(q) {
    prior <- c((1 - q)^2, 2 * q * (1 - q), q^2)
    -sum(log(drop(likelihood %*% prior)))
  }
  direct <- stats::optimize(
    negative.log.likelihood,
    interval = c(1e-5, 1 - 1e-5),
    tol = 1e-12
  )$minimum

  expect_equal(em, direct, tolerance = 1e-5)
})
