# Actual VCF-to-GDS fixtures, independent of DAPC test helpers--------------------
sm_gds_fixture <- function(folder, field = "PL", partial = FALSE) {
  header <- c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    if (field == "PL") '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihood">'
    else '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Likelihood">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "s1", "s2", "s3", "s4"), collapse = "\t"))
  records <- vapply(1:3, function(l) {
    values <- if (field == "PL") c("0,20,60", "60,20,0", "20,0,20", "0,0,0")
      else c("0,-2,-6", "-6,-2,0", "-2,0,-2", "0,0,0")
    if (l == 2L) values[4] <- "."
    if (partial && l == 1L) values[1] <- "0,.,20"
    paste(c("1", l, paste0("v", l), "A", "C", ".", "PASS", ".", field, values),
          collapse = "\t")
  }, character(1))
  vcf <- file.path(folder, paste0(field, ".vcf"))
  path <- file.path(folder, paste0(field, ".gds"))
  writeLines(c(header, records), vcf)
  SeqArray::seqVCF2GDS(vcf, path, verbose = FALSE)
  path
}

test_that("regularized GDS fits retain objective and training-only prior settings", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  path <- sm_gds_fixture(folder)
  initial <- matrix(.5, 3, 1, dimnames = list(as.character(1:3), "A"))
  weights <- setNames(c(1,0,0,0), paste0("s",1:4))
  pooled <- .sm_fit_gds(path, initial, sample.weights=weights,
    regularization="pooled", prior.strength=2, n.starts=1)
  symmetric <- .sm_fit_gds(path, initial, sample.weights=weights,
    regularization="symmetric", prior.strength=2, n.starts=1)
  expect_true(all(pooled$prior.center < .01))
  expect_equal(symmetric$prior.center, rep(.5,3))
  expect_true(all(symmetric$frequency > pooled$frequency))
  expect_true(all(diff(symmetric$objective.history) >= -1e-7))
  expect_equal(symmetric$loglik, sum(symmetric$sample.loglik*weights))
  expect_equal(symmetric$regularization, "symmetric")
  expect_true(all(is.na(symmetric$likelihood.information.audit$MEMBERSHIP_CORRELATION)))
  expect_gte(symmetric$near.best.starts,1)
  expect_error(.sm_fit_gds(path,initial,prior.strength=1), "prior.strength")
  expect_error(.sm_fit_gds(path,initial,regularization="pooled"), "prior.strength")
  expect_error(.sm_fit_gds(path,initial,init.max.snps=0), "init.max.snps")
})

test_that("real GL and PL GDS files agree and metadata does not drive fitting", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  pl <- sm_gds_fixture(folder)
  gl <- sm_gds_fixture(folder, "GL")
  frequency <- matrix(c(rep(.3, 3), rep(.7, 3)), 3, 2,
                      dimnames = list(as.character(1:3), c("A", "B")))
  metadata <- data.frame(INDIVIDUALS = paste0("s", 4:1), STRATA = c("X", "X", "Y", "Y"))
  a <- .sm_fit_gds(pl, frequency, strata = metadata)
  b <- .sm_fit_gds(gl, frequency[3:1, ], strata = transform(metadata, STRATA = "ALL"))
  expect_equal(a$membership, b$membership)
  expect_equal(a$frequency, b$frequency)
  expect_identical(a$likelihood, "PL")
  expect_identical(b$likelihood, "GL")
  expect_identical(a$strata$INDIVIDUALS, paste0("s", 1:4))
  expect_equal(unname(a$membership[4, ]), unname(a$proportion))
  expect_true(is.na(a$assignments$ASSIGNED_STOCK[4]))
  expect_equal(a$assignments$INFORMATIVE_MARKERS[4], 0)
  expect_error(.sm_fit_gds(pl, frequency, likelihood = "GL"), "does not contain")
  expect_error(.sm_fit_gds(pl, frequency, max.genotypes = 1), "RAM guard")
})

test_that("borrowed GDS filters survive success and error", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  path <- sm_gds_fixture(folder)
  gds <- SeqArray::seqOpen(path)
  on.exit(SeqArray::seqClose(gds), add = TRUE, after = FALSE)
  SeqArray::seqSetFilter(gds, sample.id = c("s1", "s3", "s4"),
                         variant.id = c(1, 3), verbose = FALSE)
  before <- SeqArray::seqGetFilter(gds)
  frequency <- matrix(.4, 2, 1, dimnames = list(c("3", "1"), "A"))
  fit <- .sm_fit_gds(gds, frequency)
  expect_identical(SeqArray::seqGetFilter(gds), before)
  expect_identical(fit$sample.ids, c("s1", "s3", "s4"))
  expect_equal(fit$variant.ids, c(1, 3))
  expect_identical(rownames(fit$frequency), c("1", "3"))
  expect_error(.sm_fit_gds(gds, unname(frequency)), "variant.id")
  expect_identical(SeqArray::seqGetFilter(gds), before)
  # A borrowed file remains readable after both calls.
  expect_identical(SeqArray::seqGetData(gds, "sample.id"), fit$sample.ids)
})

test_that("malformed GDS likelihoods fail and release owned handles", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  path <- sm_gds_fixture(folder, partial = TRUE)
  frequency <- matrix(.4, 3, 1, dimnames = list(as.character(1:3), "A"))
  expect_error(.sm_fit_gds(path, frequency), "partially missing")
  gds <- SeqArray::seqOpen(path, readonly = FALSE)
  expect_s4_class(gds, "SeqVarGDSClass")
  SeqArray::seqClose(gds)
})

test_that("stored genometranslator metadata is aligned to active samples", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  path <- sm_gds_fixture(folder)
  gds <- SeqArray::seqOpen(path, readonly = FALSE)
  on.exit(SeqArray::seqClose(gds), add = TRUE, after = FALSE)
  root <- gdsfmt::addfolder.gdsn(gds, "genometranslator")
  meta <- gdsfmt::addfolder.gdsn(root, "individuals.meta")
  gdsfmt::add.gdsn(meta, "INDIVIDUALS", paste0("s", 4:1))
  gdsfmt::add.gdsn(meta, "STRATA", c("north", "south", "east", "west"))
  SeqArray::seqSetFilter(gds, sample.id = c("s1", "s2", "s4"), verbose = FALSE)
  frequency <- matrix(.4, 3, 1, dimnames = list(as.character(1:3), "A"))
  fit <- .sm_fit_gds(gds, frequency)
  expect_identical(fit$strata$STRATA, c("west", "east", "north"))
  expect_identical(fit$strata$INDIVIDUALS, c("s1", "s2", "s4"))
  expect_error(.sm_fit_gds(gds, frequency,
    strata = data.frame(INDIVIDUALS = "s1")), "cover active")
})

test_that("GDS block sizes preserve fits, filters and random state", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  path <- sm_gds_fixture(folder)
  frequency <- matrix(c(rep(.3, 3), rep(.7, 3)), 3, 2,
                      dimnames = list(as.character(1:3), c("A", "B")))
  gds <- SeqArray::seqOpen(path)
  on.exit(SeqArray::seqClose(gds), add = TRUE, after = FALSE)
  before <- SeqArray::seqGetFilter(gds)
  set.seed(812)
  rng <- .Random.seed; kind <- RNGkind()
  a <- .sm_fit_gds(gds, frequency, n.starts = 3, seed = 37, block.snps = 1)
  expect_identical(.Random.seed, rng)
  expect_identical(RNGkind(), kind)
  b <- .sm_fit_gds(gds, frequency, n.starts = 3, seed = 37, block.snps = 2)
  expect_equal(a$frequency, b$frequency, tolerance = 1e-9)
  expect_equal(a$membership, b$membership, tolerance = 1e-9)
  expect_equal(a$start.diagnostics, b$start.diagnostics, tolerance = 1e-9)
  expect_identical(SeqArray::seqGetFilter(gds), before)
  expect_equal(nrow(a$start.diagnostics), 3)
  expect_equal(sum(a$start.diagnostics$SELECTED), 1)
  expect_equal(a$loglik, max(a$start.diagnostics$LOGLIK[a$start.diagnostics$CONVERGED]))
  # Total data exceed this limit, but each one-marker block fits.
  capped <- .sm_fit_gds(gds, frequency, n.starts = 1, max.genotypes = 4)
  expect_equal(capped$block.snps, 1)
  expect_null(capped$informative)
  expect_equal(unname(capped$sample.information), c(3, 3, 3, 0))
  # Direct parity with the in-memory reference, using the same start.
  raw <- SeqArray::seqGetData(gds, "annotation/format/PL")
  x <- array(as.numeric(raw$data), c(4, 3, 3))
  full <- .sm_fit(x, frequency, encoding = "PL")
  expect_equal(unname(capped$frequency), unname(full$frequency), tolerance = 1e-9)
  expect_equal(unname(capped$membership), unname(full$membership), tolerance = 1e-9)
  expect_equal(capped$history, full$history, tolerance = 1e-9)
})

test_that("block controls reject invalid values before fitting", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  path <- sm_gds_fixture(folder)
  frequency <- matrix(.4, 3, 1, dimnames = list(as.character(1:3), "A"))
  expect_error(.sm_fit_gds(path, frequency, block.snps = 0), "positive integer")
  expect_error(.sm_fit_gds(path, frequency, n.starts = 1.5), "positive integer")
  expect_error(.sm_fit_gds(path, frequency, seed = NA), "seed")
  expect_warning(.sm_fit_gds(path, frequency, n.starts = 1, max.iter = 1), "converge")
})
