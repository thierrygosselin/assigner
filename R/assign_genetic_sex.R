#' Assign genetic sex from a validated Y/W marker panel
#'
#' Classify individuals from the presence and absence pattern of validated
#' Y-like or W-like markers. The function consumes the assignment panel written
#' by [radr::sexy_markers()] and a GDS. It never changes recorded sex and returns
#' `U` when evidence is insufficient, weak, or contradictory.
#'
#' @section Scoring:
#' A Y-like marker supports male when present and female when absent. A W-like
#' marker supports female when present and male when absent. Female and male
#' scores are the proportions of panel markers agreeing with each pattern.
#' Assignment requires at least `min.panel.markers`, a winning score of at least
#' `min.score`, and a female-male score difference of at least `min.margin`.
#'
#' @section Validation:
#' A panel should be selected in discovery samples and evaluated in independent
#' samples, families, batches, or populations. Reusing discovery samples can
#' substantially overestimate accuracy. Marker absence is especially sensitive
#' to coverage and batch effects. Recorded sex is used only for concordance
#' reporting and never contributes to classification.
#'
#' @param data A GDS filepath or an open `SeqVarGDSClass` object.
#' @param panel A data frame or TSV written as `sex_assignment_panel.tsv` by
#'   [radr::sexy_markers()]. Required columns are `VARIANT_ID`, `MARKERS`, and
#'   `ASSIGNMENT_DIRECTION`.
#' @param metadata Optional data frame or TSV containing `INDIVIDUALS` and,
#'   optionally, `sex.column`. If `NULL`, individual metadata are read from GDS.
#' @param sex.column Recorded-sex column used only for concordance reporting.
#' @param min.panel.markers Minimum number of matched panel markers required.
#' @param min.score Minimum winning agreement proportion.
#' @param min.margin Minimum absolute difference between female and male scores.
#' @param path.folder Parent directory for the numbered result folder.
#' @param save.results Write assignment, marker-score, panel, and summary TSVs.
#' @param verbose Display progress and summary messages.
#'
#' @return An `assign_genetic_sex` object containing `assignment`,
#'   `marker_scores`, `panel`, `summary`, and `path.folder`.
#' @export
#' @examples
#' \dontrun{
#' result <- assigner::assign_genetic_sex(
#'   data = "validation.gds",
#'   panel = "sex_assignment_panel.tsv",
#'   metadata = "validation_metadata.tsv",
#'   sex.column = "SEX"
#' )
#' result$assignment
#' }
assign_genetic_sex <- function(
    data,
    panel,
    metadata = NULL,
    sex.column = "STRATA",
    min.panel.markers = 2L,
    min.score = 0.8,
    min.margin = 0.2,
    path.folder = NULL,
    save.results = TRUE,
    verbose = TRUE
) {
  force(data)
  .start <- tgbase::startup(
    package = "assigner", f.name = "assign_genetic_sex", verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)

  min.panel.markers <- .ags_count(
    min.panel.markers, "min.panel.markers", 1L
  )
  .ags_probability(min.score, "min.score")
  .ags_probability(min.margin, "min.margin")
  .ags_flag(save.results, "save.results")
  .ags_flag(verbose, "verbose")
  if (!is.character(sex.column) || length(sex.column) != 1L ||
      is.na(sex.column) || !nzchar(sex.column)) {
    rlang::abort("`sex.column` must be one non-empty column name.")
  }

  panel <- .ags_read_table(panel, "panel")
  required <- c("VARIANT_ID", "MARKERS", "ASSIGNMENT_DIRECTION")
  missing.columns <- setdiff(required, names(panel))
  if (length(missing.columns)) {
    rlang::abort(paste0(
      "The panel is missing: ", paste(missing.columns, collapse = ", "), "."
    ))
  }
  panel$ASSIGNMENT_DIRECTION <- as.character(panel$ASSIGNMENT_DIRECTION)
  if (!nrow(panel)) rlang::abort("The assignment panel contains no markers.")
  if (any(!panel$ASSIGNMENT_DIRECTION %in% c("Y-like", "W-like"))) {
    rlang::abort("`ASSIGNMENT_DIRECTION` must contain only Y-like or W-like.")
  }
  if (anyNA(panel$VARIANT_ID) || anyDuplicated(panel$VARIANT_ID)) {
    rlang::abort("Panel `VARIANT_ID` values must be unique and non-missing.")
  }

  opened.here <- FALSE
  if (inherits(data, "SeqVarGDSClass")) {
    gds <- data
  } else if (is.character(data) && length(data) == 1L && !is.na(data) &&
             file.exists(data) && grepl("\\.gds$", data, ignore.case = TRUE)) {
    gds <- SeqArray::seqOpen(data)
    opened.here <- TRUE
  } else {
    rlang::abort(
      "`data` must be a GDS filepath or an open SeqVarGDSClass object."
    )
  }
  on.exit({
    if (opened.here) try(SeqArray::seqClose(gds), silent = TRUE)
  }, add = TRUE)
  selection.before <- SeqArray::seqGetFilter(gds)
  SeqArray::seqFilterPush(gds)
  filter.pushed <- TRUE
  on.exit({
    if (filter.pushed) try(SeqArray::seqFilterPop(gds), silent = TRUE)
  }, add = TRUE)

  sample.id <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  active.variant.id <- SeqArray::seqGetData(gds, "variant.id")
  match.index <- match(as.character(panel$VARIANT_ID),
    as.character(active.variant.id))
  missing.panel <- is.na(match.index)
  if (any(missing.panel) && verbose) {
    .ags_message(
      sum(missing.panel), " panel marker(s) were absent from the active GDS."
    )
  }
  panel <- panel[!missing.panel, , drop = FALSE]
  match.index <- match.index[!missing.panel]
  if (nrow(panel) < min.panel.markers) {
    rlang::abort(paste0(
      "Only ", nrow(panel), " panel marker(s) matched the active GDS; ",
      min.panel.markers, " are required."
    ))
  }
  gds.variant.id <- active.variant.id[match.index]
  SeqArray::seqSetFilter(gds, variant.id = gds.variant.id, verbose = FALSE)

  presence.source <- if ("PRESENCE_SOURCE" %in% names(panel)) {
    unique(panel$PRESENCE_SOURCE)
  } else {
    "genotype_call"
  }
  if (length(presence.source) != 1L ||
      !presence.source %in% c("genotype_call", "read_depth")) {
    rlang::abort("The panel must use one supported `PRESENCE_SOURCE`.")
  }
  dosage <- .ags_get_dosage(gds, sample.id)
  if (presence.source == "read_depth") {
    threshold <- unique(panel$COVERAGE_THRESHOLD)
    threshold <- threshold[is.finite(threshold)]
    if (length(threshold) != 1L) {
      rlang::abort("A read-depth panel requires one `COVERAGE_THRESHOLD`.")
    }
    depth <- .ags_get_depth(gds, gds.variant.id, sample.id)
    if (is.null(depth)) {
      rlang::abort(
        "The panel requires read depth, but compatible GDS depth was not found."
      )
    }
    present <- is.finite(depth) & depth >= threshold
  } else {
    present <- is.finite(dosage)
  }
  colnames(present) <- as.character(panel$MARKERS)
  rownames(present) <- sample.id

  y.like <- panel$ASSIGNMENT_DIRECTION == "Y-like"
  expected.male <- matrix(
    rep(y.like, each = length(sample.id)), nrow = length(sample.id)
  )
  observed <- present
  male.agreement <- observed == expected.male
  female.agreement <- observed != expected.male
  male.score <- unname(rowMeans(male.agreement))
  female.score <- unname(rowMeans(female.agreement))
  winning.score <- pmax(female.score, male.score)
  score.margin <- abs(female.score - male.score)
  inferred <- ifelse(male.score > female.score, "M",
    ifelse(female.score > male.score, "F", "U"))
  sufficient <- nrow(panel) >= min.panel.markers &
    winning.score >= min.score & score.margin >= min.margin
  inferred[!sufficient] <- "U"

  sample.metadata <- .ags_metadata(metadata, gds, sample.id)
  recorded <- rep(NA_character_, length(sample.id))
  if (!is.null(sample.metadata) && sex.column %in% names(sample.metadata)) {
    sample.metadata <- sample.metadata[
      match(sample.id, sample.metadata$INDIVIDUALS), , drop = FALSE
    ]
    recorded <- .ags_normalise_sex(sample.metadata[[sex.column]])
  }
  assignment <- tibble::tibble(
    INDIVIDUALS = sample.id,
    RECORDED_SEX = recorded,
    INFERRED_SEX = inferred,
    FEMALE_SCORE = female.score,
    MALE_SCORE = male.score,
    WINNING_SCORE = winning.score,
    SCORE_MARGIN = score.margin,
    N_PANEL_MARKERS = nrow(panel),
    N_PRESENT_MARKERS = rowSums(present),
    STATUS = ifelse(inferred == "U", "unresolved", "assigned"),
    CONCORDANT = ifelse(
      recorded %in% c("F", "M") & inferred %in% c("F", "M"),
      recorded == inferred, NA
    )
  )
  present.vector <- as.vector(present)
  expected.male.vector <- rep(y.like, each = length(sample.id))
  marker.scores <- tibble::tibble(
    INDIVIDUALS = rep(sample.id, times = nrow(panel)),
    VARIANT_ID = rep(panel$VARIANT_ID, each = length(sample.id)),
    MARKERS = rep(panel$MARKERS, each = length(sample.id)),
    ASSIGNMENT_DIRECTION = rep(
      panel$ASSIGNMENT_DIRECTION, each = length(sample.id)
    ),
    PRESENT = present.vector,
    SUPPORTS_SEX = ifelse(present.vector == expected.male.vector, "M", "F")
  )
  summary <- tibble::tibble(
    metric = c(
      "samples", "assigned", "unresolved", "recorded_known",
      "recorded_concordant", "panel_markers", "y_like", "w_like"
    ),
    value = c(
      length(sample.id), sum(inferred != "U"), sum(inferred == "U"),
      sum(recorded %in% c("F", "M")), sum(assignment$CONCORDANT, na.rm = TRUE),
      nrow(panel), sum(y.like), sum(!y.like)
    )
  )

  if (is.null(path.folder)) path.folder <- getwd()
  dir.create(path.folder, recursive = TRUE, showWarnings = FALSE)
  path.folder <- normalizePath(path.folder, mustWork = TRUE)
  result.folder <- file.path(
    path.folder, paste0("assign_genetic_sex_", file.date)
  )
  if (dir.exists(result.folder)) {
    suffix <- 2L
    while (dir.exists(paste0(result.folder, "_", suffix))) suffix <- suffix + 1L
    result.folder <- paste0(result.folder, "_", suffix)
  }
  if (save.results) {
    dir.create(result.folder, recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(
      assignment, file.path(result.folder, "genetic_sex_assignments.tsv")
    )
    readr::write_tsv(
      marker.scores, file.path(result.folder, "genetic_sex_marker_scores.tsv")
    )
    readr::write_tsv(
      panel, file.path(result.folder, "sex_assignment_panel_used.tsv")
    )
    readr::write_tsv(
      summary, file.path(result.folder, "genetic_sex_summary.tsv")
    )
  } else {
    result.folder <- NA_character_
  }

  SeqArray::seqFilterPop(gds)
  filter.pushed <- FALSE
  selection.after <- SeqArray::seqGetFilter(gds)
  restored <- identical(selection.before$sample.sel, selection.after$sample.sel) &&
    identical(selection.before$variant.sel, selection.after$variant.sel)
  if (!restored) rlang::abort("The active GDS selection was not restored.")
  if (opened.here) {
    SeqArray::seqClose(gds)
    opened.here <- FALSE
  }
  if (verbose) {
    .ags_message("Panel markers used: ", nrow(panel), ".")
    .ags_message(
      "Assigned: ", sum(inferred != "U"), "; unresolved: ",
      sum(inferred == "U"), "."
    )
    if (any(!is.na(assignment$CONCORDANT))) {
      .ags_message(
        "Recorded-sex concordance: ",
        sum(assignment$CONCORDANT, na.rm = TRUE), " of ",
        sum(!is.na(assignment$CONCORDANT)), "."
      )
    }
    if (save.results) .ags_message("Results written to: ", result.folder)
  }
  result <- list(
    assignment = assignment,
    marker_scores = marker.scores,
    panel = panel,
    summary = summary,
    active_selection_restored = restored,
    path.folder = result.folder,
    settings = list(
      min.panel.markers = min.panel.markers,
      min.score = min.score,
      min.margin = min.margin,
      presence.source = presence.source
    )
  )
  class(result) <- "assign_genetic_sex"
  result
}

#' @export
print.assign_genetic_sex <- function(x, ...) {
  cat("Genetic-sex assignment\n")
  cat("  Samples: ", nrow(x$assignment), "\n", sep = "")
  cat("  Assigned: ", sum(x$assignment$INFERRED_SEX != "U"), "\n", sep = "")
  cat("  Unresolved: ", sum(x$assignment$INFERRED_SEX == "U"), "\n", sep = "")
  cat("  Panel markers: ", nrow(x$panel), "\n", sep = "")
  invisible(x)
}

.ags_read_table <- function(x, name) {
  if (is.data.frame(x)) return(tibble::as_tibble(x))
  if (is.character(x) && length(x) == 1L && !is.na(x) && file.exists(x)) {
    return(readr::read_tsv(x, show_col_types = FALSE, progress = FALSE))
  }
  rlang::abort(paste0("`", name, "` must be a data frame or TSV filepath."))
}

.ags_metadata <- function(metadata, gds, sample.id) {
  if (is.null(metadata)) {
    metadata <- tryCatch(
      genometranslator::extract_individuals_metadata(
        gds = gds, whitelist = TRUE
      ), error = function(error) NULL
    )
  } else {
    metadata <- .ags_read_table(metadata, "metadata")
  }
  if (is.null(metadata) || !"INDIVIDUALS" %in% names(metadata)) return(NULL)
  metadata$INDIVIDUALS <- as.character(metadata$INDIVIDUALS)
  metadata[metadata$INDIVIDUALS %in% sample.id, , drop = FALSE]
}

.ags_normalise_sex <- function(x) {
  value <- toupper(trimws(as.character(x)))
  out <- rep("U", length(value))
  out[value %in% c("F", "FEMALE")] <- "F"
  out[value %in% c("M", "MALE")] <- "M"
  out
}

.ags_get_dosage <- function(gds, sample.id) {
  dosage <- as.matrix(SeqArray::seqGetData(gds, "$dosage_alt"))
  if (nrow(dosage) == length(sample.id)) return(dosage)
  if (ncol(dosage) == length(sample.id)) return(t(dosage))
  rlang::abort("GDS dosage dimensions do not match active sample IDs.")
}

.ags_get_depth <- function(gds, variant.id, sample.id) {
  normalise <- function(x) {
    if (!is.matrix(x)) return(NULL)
    if (nrow(x) == length(sample.id)) return(x)
    if (ncol(x) == length(sample.id)) return(t(x))
    NULL
  }
  depth <- tryCatch(
    normalise(SeqArray::seqGetData(gds, "annotation/format/DP")),
    error = function(error) NULL
  )
  if (!is.null(depth)) return(depth)
  .ags_embedded_depth(gds, variant.id, sample.id)
}

.ags_embedded_depth <- function(gds, variant.id, sample.id) {
  metadata.node <- NULL
  for (path in c("genometranslator/genotypes.meta", "radiator/genotypes.meta")) {
    node <- gdsfmt::index.gdsn(gds$root, path = path, silent = TRUE)
    if (!is.null(node)) {
      metadata.node <- node
      break
    }
  }
  if (is.null(metadata.node) ||
      !"READ_DEPTH" %in% gdsfmt::ls.gdsn(metadata.node)) return(NULL)
  full.sample <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(gds$root, "sample.id")
  ))
  full.variant <- gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(gds$root, "variant.id")
  )
  marker.node <- gdsfmt::index.gdsn(metadata.node, "M_SEQ", silent = TRUE)
  depth.node <- gdsfmt::index.gdsn(metadata.node, "READ_DEPTH", silent = TRUE)
  if (is.null(marker.node) || is.null(depth.node)) return(NULL)
  marker <- gdsfmt::read.gdsn(marker.node)
  depth <- gdsfmt::read.gdsn(depth.node)
  expected <- length(full.sample) * length(full.variant)
  if (length(marker) != expected || length(depth) != expected) return(NULL)
  marker.matrix <- matrix(marker, nrow = length(full.sample))
  if (!all(marker.matrix[1L, ] == full.variant)) return(NULL)
  depth.matrix <- matrix(depth, nrow = length(full.sample))
  s <- match(sample.id, full.sample)
  v <- match(variant.id, full.variant)
  if (anyNA(s) || anyNA(v)) return(NULL)
  depth.matrix[s, v, drop = FALSE]
}

.ags_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    rlang::abort(paste0("`", name, "` must be one number from zero to one."))
  }
}

.ags_count <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != as.integer(x) || x < minimum) {
    rlang::abort(paste0(
      "`", name, "` must be a whole number of at least ", minimum, "."
    ))
  }
  as.integer(x)
}

.ags_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    rlang::abort(paste0("`", name, "` must be TRUE or FALSE."))
  }
}

.ags_message <- function(...) {
  text <- paste0(..., collapse = "")
  message(paste(strwrap(text, width = 80L, exdent = 2L), collapse = "\n"))
}
