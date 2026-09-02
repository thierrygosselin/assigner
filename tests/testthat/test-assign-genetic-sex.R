test_that("validated Y/W panels assign sex and restore active GDS filters", {
  skip_if_not_installed("SeqArray")

  sample.id <- c(paste0("F", 1:6), paste0("M", 1:6), "ambiguous")
  y1 <- c(rep("./.", 6), rep("0/1", 6), "0/1")
  y2 <- c(rep("./.", 6), rep("0/1", 6), "./.")
  w1 <- c(rep("0/1", 6), rep("./.", 6), "0/1")
  w2 <- c(rep("0/1", 6), rep("./.", 6), "./.")
  background <- rep("0/1", length(sample.id))
  records <- list(y1, y2, w1, w2, background)
  vcf <- tempfile(fileext = ".vcf")
  gds.file <- tempfile(fileext = ".gds")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=1000>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
      "INFO", "FORMAT", sample.id), collapse = "\t"),
    vapply(seq_along(records), function(i) paste(c(
      "1", i * 100, paste0("marker", i), "A", "G", ".", "PASS", ".",
      "GT", records[[i]]
    ), collapse = "\t"), character(1))
  ), vcf)
  SeqArray::seqVCF2GDS(vcf, gds.file, verbose = FALSE)
  gds <- SeqArray::seqOpen(gds.file)
  on.exit(SeqArray::seqClose(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, variant.id = 1:4, verbose = FALSE)
  before <- SeqArray::seqGetFilter(gds)

  panel <- tibble::tibble(
    VARIANT_ID = 1:4,
    MARKERS = paste0("marker", 1:4),
    ASSIGNMENT_DIRECTION = rep(c("Y-like", "W-like"), each = 2),
    PRESENCE_SOURCE = "genotype_call"
  )
  metadata <- tibble::tibble(
    INDIVIDUALS = sample.id,
    SEX = c(rep("F", 6), rep("M", 6), "U")
  )
  result <- assign_genetic_sex(
    data = gds,
    panel = panel,
    metadata = metadata,
    sex.column = "SEX",
    min.panel.markers = 2,
    min.score = 0.8,
    min.margin = 0.2,
    save.results = FALSE,
    verbose = FALSE
  )
  after <- SeqArray::seqGetFilter(gds)

  expect_s3_class(result, "assign_genetic_sex")
  expect_true(result$active_selection_restored)
  expect_identical(before$sample.sel, after$sample.sel)
  expect_identical(before$variant.sel, after$variant.sel)
  expect_equal(result$assignment$INFERRED_SEX[1:6], rep("F", 6))
  expect_equal(result$assignment$INFERRED_SEX[7:12], rep("M", 6))
  expect_equal(result$assignment$INFERRED_SEX[13], "U")
  expect_true(all(result$assignment$CONCORDANT[1:12]))
  expect_equal(nrow(result$marker_scores), length(sample.id) * 4)
})

test_that("panel requirements and sex labels are validated", {
  expect_equal(
    .ags_normalise_sex(c("female", "F", "MALE", "m", "unknown", NA)),
    c("F", "F", "M", "M", "U", "U")
  )
  expect_error(.ags_probability(2, "min.score"), "zero to one")
  expect_error(.ags_count(1.5, "markers", 1L), "whole number")
})
