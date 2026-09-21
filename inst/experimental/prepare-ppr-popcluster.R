# Prepare matched, training-only PopCluster inputs; never run the executable.
# Rscript prepare-ppr-popcluster.R EXISTING_MARKER_BENCHMARK NEW_OUTPUT
ppr_encode <- function(g) {
  if (!is.matrix(g) || any(!is.na(g) & !(g %in% 0:2)))
    stop("Expected a matrix containing only 0, 1, 2 or NA.")
  g[is.na(g)] <- 3L
  apply(g, 1L, paste0, collapse = "")
}
ppr_decode <- function(lines) {
  if (!length(lines) || length(unique(nchar(lines))) != 1L ||
      any(grepl("[^0123]", lines))) stop("Invalid SNP encoding.")
  g <- do.call(rbind, lapply(strsplit(lines, "", fixed = TRUE), as.integer))
  g[g == 3L] <- NA_integer_
  g
}
ppr_prepare <- function(input, output) {
  input <- normalizePath(input, mustWork = TRUE)
  if (file.exists(output)) stop("Use a new output directory.")
  manifest <- readRDS(file.path(input, "manifest.rds"))
  required <- c("ID", "CONDITION", "REPLICATE", "MARKERS", "DRIFT", "TRUE_K")
  stopifnot(all(required %in% names(manifest)), !anyDuplicated(manifest$ID))
  dir.create(output, recursive = TRUE)
  records <- vector("list", nrow(manifest))
  for (j in seq_len(nrow(manifest))) {
    row <- manifest[j, ]
    source <- file.path(input, paste(row$CONDITION, row$REPLICATE, sep = "_"), "truth.rds")
    d <- readRDS(source)
    g <- d$train[, seq_len(row$MARKERS), drop = FALSE]
    # Match run-stockr-native.R exactly. Existing benchmark is complete GT.
    stopifnot(!anyNA(g))
    keep <- which(apply(g, 2L, function(x) length(unique(x)) > 1L))
    g <- g[, keep, drop = FALSE]
    row$TRAIN_N <- nrow(g)
    row$RETAINED_L <- ncol(g)
    row$PPR_NOMINAL <- nrow(g) * row$MARKERS * row$DRIFT^2
    row$PPR_RETAINED <- nrow(g) * ncol(g) * row$DRIFT^2
    row$SOURCE_MD5 <- unname(tools::md5sum(source))
    row$EXPORTED <- row$TRUE_K == 3L && row$MARKERS == 10000L
    if (row$EXPORTED) {
      dest <- file.path(output, row$ID); dir.create(dest)
      ids <- sprintf("sample%04d", seq_len(nrow(g)))
      encoded <- ppr_encode(g)
      stopifnot(identical(unname(ppr_decode(encoded)), unname(g)))
      writeLines(as.vector(rbind(ids, encoded)), file.path(dest, "train.dat"))
      # Manual v1.7: 24 ordered settings. Mixture only, no population labels,
      # no scaling, pooled-frequency prior. Three starts = pilot, not full K scan.
      settings <- c(nrow(g), ncol(g), 1, 0, 930000+j, "train.dat", "mixture",
                    0, 3, 3, 3, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2)
      stopifnot(length(settings) == 24L)
      writeLines(as.character(settings), file.path(dest, "mixture.PcPjt"))
      # Evaluation labels are deliberately not written into the program input.
      saveRDS(list(sample.ids = ids, original.ids = rownames(g),
                   marker.indices = keep, truth = d$group,
                   source = source, source.md5 = row$SOURCE_MD5),
              file.path(dest, "evaluation-only.rds"))
    }
    records[[j]] <- row
  }
  result <- do.call(rbind, records)
  saveRDS(list(manifest = result, session = sessionInfo(),
    note = "DRIFT is the generating parameter, not measured WC84 FST; training N only."),
    file.path(output, "audit.rds"))
  print(result[, c("ID", "TRAIN_N", "RETAINED_L", "PPR_RETAINED", "EXPORTED")], row.names = FALSE)
  invisible(result)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  stopifnot(length(args) == 2L)
  ppr_prepare(args[1], args[2])
}
