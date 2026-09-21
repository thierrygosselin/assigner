# Reproducible resampling without changing the caller's RNG---------------------
.sm_seed_scope <- function(seed, code) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) || seed > .Machine$integer.max) stop("Invalid seed.")
  kind <- RNGkind()
  had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had) assign(".Random.seed", old, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  force(code)
}

# Minimum-cost bijection by dynamic programming, not greedy label matching------
# Restrict K to 12: exact assignment uses O(K * 2^K) memory/time here.
.sm_align <- function(reference, candidate) {
  if (!identical(dim(reference), dim(candidate)) || any(!is.finite(reference)) ||
      any(!is.finite(candidate))) stop("Invalid frequency matrices for alignment.")
  k <- ncol(reference)
  if (k > 12L) stop("Exact alignment currently supports at most 12 components.")
  cost <- matrix(0, k, k)
  for (i in seq_len(k)) for (j in seq_len(k))
    cost[i, j] <- mean((reference[, i] - candidate[, j])^2)
  count <- 2^k
  dp <- rep(Inf, count); dp[1] <- 0
  choice <- integer(count)
  for (mask in 0:(count - 2L)) {
    used <- as.logical(intToBits(mask)[seq_len(k)])
    row <- sum(used) + 1L
    for (j in which(!used)) {
      next.mask <- bitwOr(mask, bitwShiftL(1L, j - 1L))
      value <- dp[mask + 1L] + cost[row, j]
      if (value < dp[next.mask + 1L]) {
        dp[next.mask + 1L] <- value
        choice[next.mask + 1L] <- j
      }
    }
  }
  order <- integer(k); mask <- count - 1L
  for (i in k:1L) {
    order[i] <- choice[mask + 1L]
    mask <- bitwXor(mask, bitwShiftL(1L, order[i] - 1L))
  }
  list(order = order, mean.squared.distance = dp[count] / k)
}

# Capture warnings in audit tables, never discard failed refits silently--------
.sm_audit_fit <- function(code) {
  messages <- character()
  result <- withCallingHandlers(tryCatch(force(code), error = identity),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(fit = result, warnings = paste(unique(messages), collapse = "; "))
}

# Fixed-K whole-individual bootstrap-------------------------------------------
# Internal development interface; see inst/experimental/stock-mixture-review.md.
# Individuals are the independent resampling units. All their linked markers
# stay together. This is an ordinary nonparametric bootstrap, NOT stockBOOT's
# Bayesian bootstrap, and does not calibrate a misspecified independent-SNP
# model. No new biological populations or unseen alleles can be generated.
#' Experimental fixed-K individual bootstrap for a likelihood mixture
#'
#' Refit an ordinary nonparametric individual bootstrap, align component labels,
#' and summarize percentile intervals. This is not stockR's Bayesian bootstrap.
#' Interval coverage has not yet been established for low-coverage genomic data.
#' @param data GDS path or open SeqArray handle with the reference fit's exact
#'   active samples and variants, in the same order. Input filters are preserved.
#' @param fit Converged result from `.sm_fit_gds()` fitted with unit sample weights.
#' @param B Positive integer number of bootstrap replicates; default 200.
#' @param level Central percentile interval level, strictly between zero and one.
#' @param seed Non-negative integer seed; caller RNG state is restored.
#' @param frequency.variants Character GDS variant IDs for frequency intervals.
#'   Empty by default to avoid retaining large frequency arrays. Requested cells
#'   must have informative membership-weighted support.
#' @param n.starts Number of starts per replicate. First starts from the reference
#'   fit; others use its stored start strategy. Regularization settings are also
#'   inherited; pooled centers are re-estimated within each bootstrap replicate.
#'   This does not guarantee global optima.
#' @param max.iter,tolerance EM iteration limit and relative likelihood tolerance.
#' @param block.snps Maximum variants per likelihood block.
#' @param max.values Maximum number of retained bootstrap numeric values.
#' @return List with `intervals`, replicate `diagnostics`, IDs, seeds and settings.
#'   Intervals are withheld (`NULL`) if any replicate is unusable. Diagnostics
#'   distinguish convergence from usability and retain warnings and errors.
#' @details Samples are resampled as whole individuals, retaining their marker
#'   vectors. Conditional membership probabilities are not exponentiated by
#'   bootstrap weights. Results condition on K and the chosen marker panel.
#'   Related individuals require a different resampling design. Marker linkage,
#'   unseen alleles, component merging and uncertain K are not corrected here.
#'   Exact frequency-based label alignment supports K <= 12; a successful match
#'   does not prove that components remain biologically comparable.
#' @seealso [.sm_fit_gds()], [.sm_cv_gds()]
#' @md
#' @keywords internal
.sm_bootstrap_gds <- function(data, fit, B = 200L, level = .95, seed = 1L,
                              frequency.variants = character(), n.starts = 3L,
                              max.iter = 500L, tolerance = 1e-8,
                              block.snps = 1000L, max.values = 5e6) {
  if (!inherits(fit, "stock_likelihood_mixture") || !isTRUE(fit$converged))
    stop("Supply a converged reference fit.")
  if (any(fit$sample.weights != 1)) stop("Reference fit must use unit sampling weights.")
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B < 2 || B != floor(B))
    stop("B must be an integer >= 2.")
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must lie strictly between zero and one.")
  context <- .ld_open_gds(data)
  on.exit(.ld_close_gds(context), add = TRUE)
  ids <- as.character(SeqArray::seqGetData(context$gds, "sample.id"))
  variants <- SeqArray::seqGetData(context$gds, "variant.id")
  if (!identical(ids, fit$sample.ids) || !identical(variants, fit$variant.ids))
    stop("Active sample and variant IDs/order must match the reference fit.")
  if (anyDuplicated(frequency.variants) || !all(frequency.variants %in% rownames(fit$frequency)))
    stop("Unknown or duplicated frequency.variants.")
  rows <- match(frequency.variants, rownames(fit$frequency))
  if (length(rows) && any(!fit$frequency.estimable[rows, , drop = FALSE]))
    stop("Requested reference frequencies include unsupported stock-marker combinations.")
  n <- length(ids); k <- ncol(fit$frequency)
  if (k > 12) stop("Bootstrap alignment supports at most 12 components.")
  values <- as.double(B) * k * (n + 1 + length(rows))
  if (length(max.values) != 1L || !is.finite(max.values) || max.values < 1 || values > max.values)
    stop("Bootstrap draw storage exceeds max.values.")
  .sm_seed_scope(seed, {
    stocks <- colnames(fit$frequency)
    replicates <- paste0("BOOT_", seq_len(B))
    membership <- array(NA_real_, c(n, k, B), dimnames = list(ids, stocks, replicates))
    proportion <- matrix(NA_real_, B, k, dimnames = list(replicates, stocks))
    frequency <- array(NA_real_, c(length(rows), k, B),
                       dimnames = list(frequency.variants, stocks, replicates))
    audit <- vector("list", B)
    # Precompute only small replicate seeds. The GDS fit restores our RNG.
    seeds <- sample.int(.Machine$integer.max, B)
    for (b in seq_len(B)) {
      weights <- stats::setNames(tabulate(sample.int(n, n, replace = TRUE), n), ids)
      result <- .sm_audit_fit(.sm_fit_gds(context$gds, fit$frequency,
        proportion = fit$proportion, likelihood = fit$likelihood,
        sample.weights = weights, n.starts = n.starts, seed = seeds[b],
        max.iter = max.iter, tolerance = tolerance, epsilon = fit$epsilon,
        block.snps = block.snps, start.strategy = fit$start.strategy,
        init.max.snps = fit$init.max.snps,
        max.genotypes = fit$init.max.cells,
        regularization = fit$regularization, prior.strength = fit$prior.strength))
      current <- result$fit
      ok <- !inherits(current, "error") && isTRUE(current$converged)
      converged <- ok
      reason <- if (inherits(current, "error")) conditionMessage(current) else NA_character_
      distance <- NA_real_
      if (ok) {
        alignment <- .sm_align(fit$frequency, current$frequency)
        order <- alignment$order
        distance <- alignment$mean.squared.distance
        membership[, , b] <- current$membership[, order, drop = FALSE]
        proportion[b, ] <- current$proportion[order]
        if (length(rows)) frequency[, , b] <- current$frequency[rows, order, drop = FALSE]
        if (length(rows) && any(!current$frequency.estimable[rows, order, drop = FALSE])) {
          ok <- FALSE
          reason <- "Requested frequencies unsupported in this replicate."
        }
      }
      audit[[b]] <- tibble::tibble(REPLICATE = b, CONVERGED = converged, USABLE = ok,
        ALIGNMENT_DISTANCE = distance, WARNINGS = result$warnings,
        ERROR = reason)
    }
    diagnostics <- dplyr::bind_rows(audit)
    # Failed replicates can be non-random. Do not hide them by reporting
    # intervals conditional on only successful optimization runs.
    intervals <- NULL
    if (all(diagnostics$USABLE)) {
      probabilities <- c((1 - level) / 2, .5, 1 - (1 - level) / 2)
      intervals <- list(
        membership = apply(membership, c(1, 2), stats::quantile, probs = probabilities),
        proportion = apply(proportion, 2, stats::quantile, probs = probabilities),
        frequency = if (length(rows)) apply(frequency, c(1, 2), stats::quantile,
                                             probs = probabilities) else NULL)
    } else warning("Incomplete bootstrap: intervals withheld; inspect diagnostics and rerun.")
    list(intervals = intervals, diagnostics = diagnostics, B = B, level = level,
         seed = seed, sample.ids = ids, stocks = colnames(fit$frequency),
         frequency.variants = frequency.variants,
         method = "whole-individual nonparametric bootstrap, conditional on K and markers")
  })
}

# Compare K by held-out individual predictive log likelihood--------------------
# Fold labels are independent sampling units chosen by the analyst; relatives
# and repeated samples should be placed in the same fold. This never treats
# STRATA as known stocks. Initialization is data-independent (0.5), and only
# training individuals update frequencies, proportions and additional starts.
#' Experimental held-out predictive comparison of likelihood mixtures
#'
#' Fit each K using training samples and score held-out marginal likelihoods.
#' This diagnostic is not validated as an estimator of biological stock count:
#' the stockR three-group example exposes sensitivity to frequency bounds and
#' initialization. Do not interpret the returned winner as demographic evidence.
#' @param data GDS path or open SeqArray handle. Active filters are preserved.
#' @param K Unique positive integers to compare; must include 1.
#' @param folds Fold labels named by every active sample ID, with at least two
#'   nonempty folds. Keep related or repeated samples in the same fold. Each
#'   training set must contain more individuals than the fitted K.
#' @param n.starts Starts per fold and K; at least two when K exceeds one.
#'   The first start is symmetric and additional starts break that symmetry.
#' @param seed Non-negative integer seed; caller RNG state is restored.
#' @param max.iter,tolerance EM iteration limit and relative likelihood tolerance.
#' @param epsilon Allele-frequency bound strictly between zero and one half.
#'   This is not an estimated prior. Rare held-out alleles make predictive scores
#'   sensitive to this bound, particularly for small training components.
#' @param block.snps Maximum variants per likelihood block.
#' @param likelihood `"auto"` prefers PL then GL; neither GT nor GP is substituted.
#' @param start.strategy,init.max.snps Initialization controls; see `.sm_fit_gds()`.
#' @param regularization,prior.strength Fixed regularization specification; see
#'   `.sm_fit_gds()`. Pooled centers are estimated separately in each training
#'   fold. Held-out scores exclude prior penalties. Hyperparameter selection on
#'   these same scores requires a separate outer evaluation to avoid optimism.
#' @return List containing per-K `summary`, per-fold `diagnostics`, per-sample
#'   `heldout` scores, `best.predictive.K`, boundary flag, folds and settings.
#'   Any incomplete K comparison withholds the winner. Exact ties select smaller K.
#' @details `fold.support` retains all tied fold winners. `summary` includes
#'   `GAP_FROM_BEST`; `support.flags` reports incomplete comparisons, differing
#'   fold winners, grid boundaries and fitting warnings. These are descriptive:
#'   overlapping training folds do not supply independent confidence estimates.
#'   K=1 means no supported split in this comparison, not proven panmixia. Absence
#'   of flags does not validate structure. Across-fold membership agreement is
#'   not calculated. No universal support or assignment cutoff is imposed.
#'
#' Held-out samples have zero fitting weight. Frequencies and mixture
#'   proportions are estimated from training samples only. Scores integrate over
#'   stock membership, rather than selecting the best component. Uninformative
#'   held-out individuals are excluded and counts are reported. Upstream marker
#'   selection and quality control are not repeated within folds. Compare scores
#'   only for the same observations; normalized likelihood constants cancel in
#'   that comparison, not between different datasets. Convergence does not prove
#'   optimal initialization, calibrated probabilities or correct model assumptions.
#' @seealso [.sm_fit_gds()], [.sm_bootstrap_gds()]
#' @md
#' @keywords internal
.sm_cv_gds <- function(data, K = 1:4, folds, n.starts = 5L, seed = 1L,
                       max.iter = 500L, tolerance = 1e-8, epsilon = 1e-8,
                       block.snps = 1000L, likelihood = c("auto", "PL", "GL"),
                       start.strategy = c("random", "mixed"), init.max.snps = 1000L,
                       regularization = c("none", "symmetric", "pooled"), prior.strength = 0) {
  likelihood <- match.arg(likelihood)
  start.strategy <- match.arg(start.strategy)
  regularization <- match.arg(regularization)
  if (!is.numeric(n.starts) || length(n.starts) != 1L || !is.finite(n.starts) ||
      n.starts < 1 || n.starts != floor(n.starts)) stop("n.starts must be a positive integer.")
  if (!is.numeric(K) || !length(K) || any(!is.finite(K)) || any(K < 1 | K != floor(K)) ||
      anyDuplicated(K) || !1 %in% K) stop("K must contain unique positive integers including 1.")
  context <- .ld_open_gds(data)
  on.exit(.ld_close_gds(context), add = TRUE)
  ids <- as.character(SeqArray::seqGetData(context$gds, "sample.id"))
  variants <- as.character(SeqArray::seqGetData(context$gds, "variant.id"))
  if (is.null(names(folds)) || anyDuplicated(names(folds)) || !setequal(names(folds), ids) ||
      anyNA(folds)) stop("folds must uniquely name every active sample.")
  folds <- as.character(folds[ids])
  if (any(!nzchar(folds)) || length(unique(folds)) < 2) stop("At least two nonempty folds required.")
  if (n.starts < 2 && any(K > 1)) stop("K > 1 requires at least two starts to break symmetry.")
  .sm_seed_scope(seed, {
    audit <- list(); samples <- list(); z <- 0L
    for (fold in unique(folds)) for (k in sort(K)) {
      z <- z + 1L
      train <- folds != fold
      if (sum(train) <= k) stop("Each training fold must have more individuals than K.")
      initial <- matrix(.5, length(variants), k,
                        dimnames = list(variants, paste0("STOCK_", seq_len(k))))
      result <- .sm_audit_fit(.sm_fit_gds(context$gds, initial,
        likelihood = likelihood,
        sample.weights = stats::setNames(as.numeric(train), ids),
        n.starts = n.starts, seed = sample.int(.Machine$integer.max, 1),
        max.iter = max.iter, tolerance = tolerance, epsilon = epsilon,
        block.snps = block.snps, start.strategy = start.strategy, init.max.snps = init.max.snps,
        regularization = regularization, prior.strength = prior.strength))
      current <- result$fit
      ok <- !inherits(current, "error") && isTRUE(current$converged)
      score <- NA_real_; count <- 0L
      minimum.component <- NA_real_; boundary <- NA_real_
      if (ok) {
        minimum.component <- min(colSums(current$membership[train, , drop = FALSE]))
        boundary <- mean(current$frequency <= epsilon | current$frequency >= 1 - epsilon)
        test <- !train & current$sample.information > 0
        count <- sum(test)
        if (count) {
          score <- sum(current$sample.loglik[test])
          samples[[z]] <- tibble::tibble(K = k, FOLD = fold,
            INDIVIDUALS = ids[test], LOGLIK = unname(current$sample.loglik[test]))
        } else ok <- FALSE
      }
      audit[[z]] <- tibble::tibble(K = k, FOLD = fold, CONVERGED = ok,
        N_TEST = count, TEST_LOGLIK = score, MIN_COMPONENT_N = minimum.component,
        FRACTION_BOUNDARY = boundary, WARNINGS = result$warnings,
        ERROR = if (inherits(current, "error")) conditionMessage(current) else NA_character_)
    }
    diagnostics <- dplyr::bind_rows(audit)
    summary <- diagnostics |>
      dplyr::group_by(K) |>
      dplyr::summarise(COMPLETE = all(CONVERGED),
        N_TEST = sum(N_TEST), TOTAL_LOGLIK = if (all(CONVERGED)) sum(TEST_LOGLIK) else NA_real_,
        .groups = "drop") |>
      dplyr::mutate(MEAN_LOGLIK = TOTAL_LOGLIK / N_TEST)
    best <- NA_integer_
    if (all(summary$COMPLETE)) {
      best <- summary$K[which.max(summary$MEAN_LOGLIK)]
    } else warning("Incomplete K comparison: no best K reported; inspect fold diagnostics.")
    # Fold scores share training observations: gaps and votes are descriptive,
    # not independent replicates, confidence intervals or posterior K weights.
    fold.support <- diagnostics |>
      dplyr::group_by(FOLD) |>
      dplyr::mutate(FOLD_COMPLETE = all(CONVERGED & is.finite(TEST_LOGLIK) & N_TEST > 0),
        MEAN_LOGLIK = TEST_LOGLIK / N_TEST,
        FOLD_WINNER = if (all(CONVERGED & is.finite(TEST_LOGLIK) & N_TEST > 0))
          MEAN_LOGLIK == max(MEAN_LOGLIK) else rep(NA, dplyr::n())) |>
      dplyr::ungroup()
    summary$GAP_FROM_BEST <- if (is.na(best)) NA_real_ else
      max(summary$MEAN_LOGLIK) - summary$MEAN_LOGLIK
    flags <- character()
    if (is.na(best)) flags <- c(flags, "incomplete_comparison")
    if (!is.na(best) && best == 1L) flags <- c(flags, "no_supported_split_not_proof_of_panmixia")
    if (!is.na(best) && best == max(K)) flags <- c(flags, "winner_at_grid_boundary")
    if (!is.na(best) && any(fold.support$FOLD_WINNER & fold.support$K != best, na.rm = TRUE))
      flags <- c(flags, "fold_winners_disagree_or_tie")
    if (any(nzchar(diagnostics$WARNINGS))) flags <- c(flags, "fit_warnings_recorded")
    list(summary = summary, diagnostics = diagnostics, fold.support = fold.support,
         support.flags = flags,
         support.interpretation = "Descriptive support only; no calibrated K confidence or demographic-stock claim",
         heldout = dplyr::bind_rows(samples), best.predictive.K = best,
         best.at.upper.boundary = !is.na(best) && best == max(K),
         folds = stats::setNames(folds, ids), seed = seed,
         epsilon = epsilon, likelihood = likelihood, start.strategy = start.strategy,
         regularization = regularization, prior.strength = prior.strength,
         interpretation = "Predictive comparison for this K grid, not a biological stock count")
  })
}
