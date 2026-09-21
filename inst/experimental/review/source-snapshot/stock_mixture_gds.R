#' Experimental GDS interface for a genotype-likelihood stock mixture
#'
#' Fit the internal blockwise mixture prototype using active GDS samples
#' and variants. This is a development interface, not a validated stock
#' delineation workflow. Scott Foster's review is required before release.
#'
#' @param data An open SeqArray GDS object or an existing GDS path. Paths are
#'   opened read-only and closed on exit. Borrowed handles remain open; active
#'   filters are preserved on success and error.
#' @param frequency Starting ALT frequencies, a numeric marker-by-stock matrix.
#'   Row names must be the active GDS `variant.id` values, not VCF marker names.
#'   Rows are matched by ID; extra, missing and duplicate IDs are rejected.
#'   Columns define arbitrary component labels, not biological stock identities.
#' @param proportion Optional starting mixture proportions in frequency-column
#'   order. If named, names must match component labels and are reordered.
#' @param likelihood `"auto"` prefers PL, then GL. It never falls back to GT
#'   or GP. Select GL explicitly when both fields exist and PL is unsuitable.
#' @param strata Optional metadata data frame containing unique `INDIVIDUALS`.
#'   Otherwise read genometranslator metadata when present. Extra metadata rows
#'   are allowed; all active sample IDs must be covered. Metadata including
#'   `STRATA` is returned for diagnostics only, never used in fitting.
#' @param max.iter,tolerance,epsilon Controls passed to the internal EM engine.
#' @param block.snps Maximum markers read per block (default 1000).
#' @param n.starts Number of sequential starts (default 5). First uses supplied
#'   frequencies; remaining starts follow `start.strategy`. No STRATA constraints
#'   are used. Seven starts cover one mixed cycle; additional cycles can split
#'   other lower-K components. No fixed number guarantees the global optimum.
#' @param start.strategy `"random"` preserves random soft starts. Opt-in `"mixed"`
#'   cycles through PCA/k-means starts using 2, 5 and 10 PCs (rank permitting), a
#'   random soft start, a split of a fitted K-1 component, and a logit-frequency
#'   perturbation of the best solution. Split seeds use five random lower-K starts.
#'   Failed starts are reported, not silently replaced. For K=1 the split slot
#'   uses a random start. K-means uses 20 restarts, not simulated annealing.
#' @param init.max.snps Maximum evenly spaced active markers used for PCA starts.
#'   Training dosages are temporarily completed with training pooled expectations
#'   only for initialization. Training-only weighted PCA uses at most 10 axes.
#'   The sample-marker seed matrix is capped by `max.genotypes`; SVD workspace
#'   adds memory. Mean completion can still reflect structured missingness.
#' @param regularization `"none"` retains maximum likelihood. `"symmetric"` uses
#'   a fixed prior center of 0.5; `"pooled"` uses training-only pooled ML frequencies.
#'   The center is fixed before fitting each K; this is empirical regularization,
#'   not a jointly estimated hierarchical model or integrated predictive posterior.
#' @param prior.strength Nonnegative pseudo-allele count. Must be zero for `"none"`
#'   and positive otherwise. For center q and strength s the Beta prior has
#'   parameters 1+s*q and 1+s*(1-q). Choose strength independently of final test
#'   data; excessive shrinkage may obscure private alleles or weak differentiation.
#' @param seed Non-negative integer seed. Caller RNG state and kind are restored.
#' @param sample.weights Optional non-negative sampling weights named by active
#'   sample ID. Zero-weight samples are projected but do not estimate parameters.
#'   Weights multiply individual log mixture likelihoods, not posterior odds.
#' @param max.genotypes Maximum sample-marker cells in one block (default five
#'   million). Blocks shrink automatically to meet this bound. This is not a
#'   total RAM limit: frequency and membership matrices remain in memory.
#'   Vectorized updates also use several temporary arrays proportional to the
#'   block size, processing one component at a time rather than all components.
#'
#' @return Internal fit with named membership/frequency matrices, sample and
#'   variant IDs, aligned metadata, assignments, selected field, call and
#'   package version. Uninformative samples retain prior membership but have
#'   no assigned stock. Assignments with exactly tied probabilities are also NA.
#'   `frequency.support` contains membership-weighted informative sample counts;
#'   `frequency.estimable` flags counts above 0.0000000001. Unsupported frequency
#'   values are optimizer placeholders, not identified allele-frequency estimates.
#'   `objective` and `objective.history` include the log-prior penalty, whereas
#'   `loglik`, `history` and `sample.loglik` contain data likelihood only. Starts
#'   are ranked by objective among converged fits. Prior parameters and seed
#'   settings are returned. Penalized objectives must not be compared across K.
#'   `likelihood.information.audit` reports training membership associations with the
#'   proportion of flat/missing likelihood vectors (not sequencing depth or a
#'   literal missing-GT rate). Correlation can reflect biological or technical
#'   structure; it is a diagnostic, not a contamination test.
#' @details Reporting includes assignment probability margin, normalized entropy
#'   and status. These are conditional diagnostics, not calibrated confidence;
#'   K=1 ambiguity metrics are NA. Start diagnostics include label-invariant
#'   partition ARI against the selected start and comparable training sample counts,
#'   excluding ties. Partition storage is bounded by `max.genotypes`; when exceeded,
#'   comparison is explicitly not assessed. No confidence cutoff is imposed.
#'
#' Requires diploid biallelic likelihoods in REF/REF, REF/ALT, ALT/ALT
#'   order. Three likelihoods alone do not prove that the caller used a diploid
#'   model; verify the original calling settings. GL is log10 likelihood, PL is
#'   non-negative integer Phred likelihood. Fully missing and flat vectors are
#'   uninformative, partial vectors are rejected. All-uninformative markers stop
#'   the analysis rather than silently changing the marker set.
#'
#' No genotype calls are made and no GDS nodes or files are modified. This
#'   version streams likelihood blocks and selects the highest-objective
#'   converged start. Normalization and EM updates are vectorized within each
#'   bounded block; the complete likelihood array is never loaded by this interface.
#'   If none converge, it returns the best unfinished fit with
#'   warnings. Per-start diagnostics include failures, convergence and likelihood
#'   differences; equal likelihoods do not prove identical partitions. Convergence does not
#'   establish a global optimum. Membership is not ancestry, and genetic
#'   components do not establish demographic independence. Record the Git
#'   commit as well as the package version when using this development code.
#'
#' The bundled stockR hard-genotype example validates a matched-start limiting
#'   case, not low-coverage likelihood calibration or reliable K selection.
#'   Inspect frequency support and repeat starts: small components, rare alleles,
#'   linked markers, relatives, admixture and batch-specific missingness can
#'   undermine interpretation. Frequencies constrained by `epsilon` can strongly
#'   affect prediction of held-out individuals. The bound is not an estimated
#'   sequencing-error rate or a shrinkage prior. Bootstrap intervals condition
#'   on K; their coverage remains to be established.
#' @seealso [.sm_bootstrap_gds()], [.sm_cv_gds()]
#' @md
#' @keywords internal
.sm_fit_gds <- function(data, frequency, proportion = NULL,
                        likelihood = c("auto", "PL", "GL"), strata = NULL,
                        max.iter = 500L, tolerance = 1e-8, epsilon = 1e-8,
                        max.genotypes = 5e6, block.snps = 1000L,
                        n.starts = 5L, seed = 1L, sample.weights = NULL,
                        start.strategy = c("random", "mixed"), init.max.snps = 1000L,
                        regularization = c("none", "symmetric", "pooled"),
                        prior.strength = 0) {
  call <- match.call()
  start.strategy <- match.arg(start.strategy)
  regularization <- match.arg(regularization)
  if (length(prior.strength) != 1L || !is.numeric(prior.strength) ||
      !is.finite(prior.strength) || prior.strength < 0 ||
      (regularization == "none" && prior.strength != 0) ||
      (regularization != "none" && prior.strength == 0))
    stop("Use zero prior.strength for none, or a positive strength for regularization.")
  if (length(init.max.snps) != 1L || !is.numeric(init.max.snps) ||
      !is.finite(init.max.snps) || init.max.snps < 1 || init.max.snps != floor(init.max.snps))
    stop("init.max.snps must be a positive integer.")
  likelihood <- match.arg(likelihood)
  context <- .ld_open_gds(data)
  on.exit(.ld_close_gds(context), add = TRUE)
  gds <- context$gds

  # Preserve GDS order, not the caller's metadata or starting-matrix order-------
  ids <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  variants <- SeqArray::seqGetData(gds, "variant.id")
  keys <- as.character(variants)
  if (!length(ids) || !length(keys)) stop("No active samples or variants.")
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      anyNA(keys) || anyDuplicated(keys)) stop("GDS IDs must be unique and non-missing.")
  if (is.null(sample.weights)) sample.weights <- stats::setNames(rep(1, length(ids)), ids)
  if (!is.numeric(sample.weights) || is.null(names(sample.weights)) ||
      anyDuplicated(names(sample.weights)) || !setequal(names(sample.weights), ids) ||
      any(!is.finite(sample.weights)) || any(sample.weights < 0) ||
      !is.finite(sum(sample.weights)) || sum(sample.weights) <= 0)
    stop("sample.weights must uniquely name every active sample and be non-negative.")
  sample.weights <- sample.weights[ids]
  if (length(max.genotypes) != 1L || !is.numeric(max.genotypes) ||
      !is.finite(max.genotypes) || max.genotypes < 1) stop("Invalid max.genotypes.")
  for (name in c("block.snps", "n.starts", "max.iter")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max)
      stop(name, " must be a positive integer.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) stop("Invalid seed.")
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon <= 0 || epsilon >= .5 || !is.numeric(tolerance) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0)
    stop("Invalid epsilon or tolerance.")
  block.snps <- min(block.snps, floor(max.genotypes / length(ids)), length(keys))
  if (block.snps < 1) stop("One marker exceeds the max.genotypes RAM guard.")
  alleles <- as.character(SeqArray::seqGetData(gds, "allele"))
  if (anyNA(alleles) || any(lengths(strsplit(alleles, ",", fixed = TRUE)) != 2L)) {
    stop("Only biallelic variants are supported.")
  }
  if (!is.matrix(frequency) || !is.numeric(frequency) ||
      is.null(rownames(frequency)) || anyDuplicated(rownames(frequency)) ||
      !setequal(rownames(frequency), keys) || ncol(frequency) < 1L) {
    stop("frequency requires one named row per active variant.id and at least one stock.")
  }
  frequency <- frequency[match(keys, rownames(frequency)), , drop = FALSE]
  if (any(!is.finite(frequency)) || any(frequency < 0 | frequency > 1))
    stop("Starting frequencies must lie in [0,1].")
  frequency <- pmin(pmax(frequency, epsilon), 1 - epsilon)
  stocks <- colnames(frequency)
  if (is.null(stocks)) stocks <- paste0("STOCK_", seq_len(ncol(frequency)))
  if (anyNA(stocks) || any(!nzchar(stocks)) || anyDuplicated(stocks)) {
    stop("Component names must be unique and non-empty.")
  }
  if (!is.null(proportion) && !is.null(names(proportion))) {
    if (anyDuplicated(names(proportion)) || !setequal(names(proportion), stocks)) {
      stop("Named proportions must match component names.")
    }
    proportion <- proportion[stocks]
  }
  if (is.null(proportion)) proportion <- rep(1 / length(stocks), length(stocks))
  if (!is.numeric(proportion) || length(proportion) != length(stocks) ||
      any(!is.finite(proportion)) || any(proportion <= 0)) stop("Invalid proportions.")
  proportion <- proportion / sum(proportion)
  if (any(proportion == 0) || any(!is.finite(proportion))) stop("Invalid proportions.")

  # Metadata is an audit layer; it never constrains inferred membership---------
  if (is.null(strata)) {
    node <- gdsfmt::index.gdsn(gds, "genometranslator/individuals.meta", silent = TRUE)
    if (!is.null(node)) {
      strata <- genometranslator::extract_individuals_metadata(gds, whitelist = FALSE)
    }
  }
  if (is.null(strata)) strata <- tibble::tibble(INDIVIDUALS = ids)
  if (!is.data.frame(strata) || !"INDIVIDUALS" %in% names(strata)) {
    stop("strata must contain INDIVIDUALS.")
  }
  metadata.ids <- as.character(strata$INDIVIDUALS)
  if (anyNA(metadata.ids) || anyDuplicated(metadata.ids) ||
      !all(ids %in% metadata.ids)) stop("Metadata must uniquely cover active samples.")
  strata <- tibble::as_tibble(strata[match(ids, metadata.ids), , drop = FALSE])

  # Read genuine likelihoods without completing or hard-calling genotypes-------
  field <- .ld_choose_likelihood(gds, likelihood)
  blocks <- split(seq_along(keys), ceiling(seq_along(keys) / block.snps))
  read.block <- function(index) {
    # Restrict only variants: retain caller's active sample selection. The
    # outer filter stack restores the original selection on every exit path.
    SeqArray::seqSetFilter(gds, variant.id = variants[index], verbose = FALSE)
    raw <- SeqArray::seqGetData(gds, paste0("annotation/format/", field))
    if (!is.list(raw) || !all(c("length", "data") %in% names(raw)) ||
        length(raw$length) != length(index) || anyNA(raw$length) ||
        any(raw$length != 3L) || !is.matrix(raw$data) ||
        !identical(dim(raw$data), c(length(ids), 3L * length(index)))) {
      stop("Expected exactly three diploid likelihoods per active sample and variant.")
    }
    x <- array(as.numeric(raw$data), c(length(ids), 3L, length(index)))
    # SeqArray missing floating-point GL values use a NaN sentinel.
    x[is.nan(x)] <- NA_real_
    .sm_prepare(x, field)
  }
  # Preflight once, before trying starts: malformed data are not optimization
  # failures. Store counts, not a full sample-by-marker information matrix.
  information <- numeric(length(ids))
  marker.information <- integer(length(keys))
  training.information <- numeric(length(keys))
  for (index in blocks) {
    block <- read.block(index)
    information <- information + rowSums(block$informative)
    marker.information[index] <- colSums(block$informative)
    training.information[index] <- colSums(block$informative * sample.weights)
  }
  if (any(marker.information == 0L)) stop("Uninformative marker(s): ",
    paste(head(keys[marker.information == 0L], 10), collapse = ", "))
  if (any(training.information == 0)) stop("No training information for marker(s): ",
    paste(head(keys[training.information == 0], 10), collapse = ", "))
  prior.center <- rep(.5, length(keys))
  if (regularization == "pooled") {
    pooled <- .sm_block_fit(read.block, blocks, matrix(.5, length(keys), 1), 1,
      information, max.iter, tolerance, epsilon, sample.weights)
    if (!pooled$converged) stop("Training pooled frequency estimate did not converge.")
    prior.center <- pooled$frequency[, 1]
  }
  fit <- .sm_multistart(read.block, blocks, frequency, proportion, information,
    n.starts, seed, max.iter, tolerance, epsilon, sample.weights,
    start.strategy, init.max.snps, max.genotypes, prior.center, prior.strength)
  fit$regularization <- regularization
  # An absent stock-marker combination has no estimable frequency even if the
  # optimizer retained a finite starting value. Expose support rather than
  # presenting those values as evidence for fixed/private alleles.
  support <- matrix(0, length(keys), length(stocks), dimnames = list(keys, stocks))
  for (index in blocks) {
    block <- read.block(index)
    support[index, ] <- crossprod(block$informative, fit$membership * sample.weights)
  }
  fit$frequency.support <- support
  fit$frequency.estimable <- support > 1e-10
  fit$sample.weights <- sample.weights
  names(fit$sample.loglik) <- ids
  fit$sample.information <- stats::setNames(information, ids)
  uninformative <- 1-information/length(keys)
  training <- sample.weights > 0 & information > 0
  fit$likelihood.information.audit <- tibble::tibble(
    STOCK = stocks,
    MEAN_UNINFORMATIVE_PROPORTION = vapply(seq_along(stocks), function(j)
      stats::weighted.mean(uninformative, fit$membership[, j]*sample.weights), numeric(1)),
    MEMBERSHIP_CORRELATION = vapply(seq_along(stocks), function(j) {
      x <- fit$membership[training,j]; y <- uninformative[training]
      if (length(x)<2L || stats::sd(x)==0 || stats::sd(y)==0) NA_real_
      else stats::cor(x,y)
    }, numeric(1)))
  fit$marker.information <- stats::setNames(marker.information, keys)
  fit$block.snps <- block.snps
  dimnames(fit$membership) <- list(ids, stocks)
  dimnames(fit$frequency) <- list(keys, stocks)
  names(fit$proportion) <- stocks
  best <- apply(fit$membership, 1L, max)
  tied <- rowSums(fit$membership == best) > 1L
  assigned <- stocks[max.col(fit$membership, ties.method = "first")]
  assigned[information == 0L | tied] <- NA_character_
  fit$assignments <- tibble::tibble(
    INDIVIDUALS = ids, ASSIGNED_STOCK = assigned, PROBABILITY = best,
    INFORMATIVE_MARKERS = unname(information)
  )
  # Descriptive ambiguity, not calibrated accuracy or an invented cutoff.
  runner.up <- if (length(stocks) > 1L)
    apply(fit$membership, 1, function(p) sort(p, decreasing = TRUE)[2]) else rep(NA_real_, length(ids))
  entropy <- if (length(stocks) > 1L) -rowSums(ifelse(fit$membership > 0,
    fit$membership * log(fit$membership), 0)) / log(length(stocks)) else rep(NA_real_, length(ids))
  fit$assignments$PROBABILITY_MARGIN <- best - runner.up
  fit$assignments$NORMALIZED_ENTROPY <- entropy
  fit$assignments$ASSIGNMENT_STATUS <- ifelse(information == 0L, "no_information",
    ifelse(!fit$converged, "unfinished_fit", ifelse(length(stocks) == 1L,
      "single_component_not_evidence_of_panmixia", ifelse(tied, "tied", "conditional_assignment"))))
  fit$sample.ids <- ids
  fit$variant.ids <- variants
  fit$strata <- strata
  fit$likelihood <- field
  fit$call <- call
  fit$package.version <- as.character(utils::packageVersion("assigner"))
  class(fit) <- c("stock_likelihood_mixture", "list")
  fit
}
