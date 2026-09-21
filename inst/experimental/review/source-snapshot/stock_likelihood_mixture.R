# Experimental genotype-likelihood stock mixture-------------------------------
# Statistical reference: Foster SD, Feutry P, Grewe PM, Berry O, Hui FKC,
# Davies CR (2018), Molecular Ecology Resources 18:1310-1325.
# https://doi.org/10.1111/1755-0998.12920
# Inspected stockR 1.0.76: calcCondProbs, calcCondMeans, calcPis, calcTaus
# (Scott D. Foster; GPL >= 2). This is a new likelihood extension, not a
# verbatim copy. Scott's statistical/code review is required before release.
# Internal prototype: no automatic K selection, bootstrap,
# sample-group constraints or production-scale memory guarantees yet.

# Validate and normalize likelihoods-------------------------------------------
.sm_prepare <- function(x, encoding = c("GL", "PL", "log")) {
  encoding <- match.arg(encoding)
  if (!is.numeric(x) || length(dim(x)) != 3L || dim(x)[2L] != 3L ||
      any(dim(x) == 0L)) stop("Expected sample x 3 genotypes x marker array.")
  if (any(is.nan(x))) stop("NaN likelihoods are not supported.")
  if (encoding == "PL" && any(x < 0 | x != floor(x), na.rm = TRUE)) {
    stop("PL must contain non-negative whole numbers.")
  }
  if (encoding == "PL") x <- -x * log(10) / 10
  if (encoding == "GL") x <- x * log(10)
  # Vectorize within the caller's bounded marker block. Keep matrix dimensions
  # explicitly, including one-sample and one-marker blocks.
  n <- dim(x)[1L]; m <- dim(x)[3L]
  a <- matrix(x[, 1L, ], n, m)
  b <- matrix(x[, 2L, ], n, m)
  c <- matrix(x[, 3L, ], n, m)
  missing <- is.na(a) & is.na(b) & is.na(c)
  if (any((is.na(a) | is.na(b) | is.na(c)) & !missing) ||
      any(x == Inf, na.rm = TRUE) ||
      any(a == -Inf & b == -Inf & c == -Inf, na.rm = TRUE)) {
    stop("Malformed or partially missing likelihood vector.")
  }
  a[missing] <- b[missing] <- c[missing] <- 0
  peak <- pmax(a, b, c)
  a <- a - peak; b <- b - peak; c <- c - peak
  x[, 1L, ] <- a; x[, 2L, ] <- b; x[, 3L, ] <- c
  list(log.likelihood = x, informative = a != 0 | b != 0 | c != 0)
}

# E step: membership and stock-conditional expected allele counts---------------
.sm_estep <- function(x, frequency, proportion) {
  n <- dim(x)[1L]; loci <- dim(x)[3L]; k <- ncol(frequency)
  conditional <- matrix(0, n, k)
  dosage <- array(0, c(n, loci, k))
  for (j in seq_len(k)) for (l in seq_len(loci)) {
    prior <- stats::dbinom(0:2, 2, frequency[l, j], log = TRUE)
    for (i in seq_len(n)) {
      joint <- x[i, , l] + prior
      peak <- max(joint)
      normalizer <- peak + log(sum(exp(joint - peak)))
      # A constant likelihood is exactly neutral, not a weak genotype call.
      conditional[i, j] <- conditional[i, j] +
        if (all(x[i, , l] == 0)) 0 else normalizer
      dosage[i, l, j] <- sum((0:2) * exp(joint - normalizer))
    }
  }
  joint <- sweep(conditional, 2L, log(proportion), "+")
  peak <- apply(joint, 1L, max)
  normalizer <- peak + log(rowSums(exp(joint - peak)))
  list(membership = exp(joint - normalizer), dosage = dosage,
       loglik = sum(normalizer), conditional = conditional)
}

# One deterministic EM fit from explicit starting values-----------------------
# Frequencies are constrained to [epsilon, 1-epsilon] throughout both steps.
# Unlike stockR's likelihood-only floor, this maximizes a consistently bounded
# objective. All-missing/flat vectors are excluded from parameter updates.
# Entirely uninformative samples receive the fitted mixture prior but do not
# estimate it; entirely uninformative markers are rejected as unidentifiable.
.sm_fit <- function(x, frequency, proportion = NULL, encoding = "GL",
                    max.iter = 500L, tolerance = 1e-8, epsilon = 1e-8) {
  prepared <- .sm_prepare(x, encoding)
  x <- prepared$log.likelihood
  observed <- prepared$informative
  if (any(colSums(observed) == 0L)) stop("Uninformative marker: no frequency information.")
  if (!is.matrix(frequency) || !is.numeric(frequency) ||
      nrow(frequency) != dim(x)[3L] || ncol(frequency) < 1L ||
      any(!is.finite(frequency)) || any(frequency < 0 | frequency > 1)) {
    stop("frequency must be a marker x stock matrix in [0,1].")
  }
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0 || epsilon >= 0.5 ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
      length(max.iter) != 1L || !is.finite(max.iter) || max.iter < 1 || max.iter != floor(max.iter)) {
    stop("Invalid epsilon, tolerance or max.iter.")
  }
  k <- ncol(frequency)
  if (is.null(proportion)) proportion <- rep(1 / k, k)
  if (length(proportion) != k || any(!is.finite(proportion)) || any(proportion <= 0)) {
    stop("Provide one positive mixture proportion per stock.")
  }
  proportion <- proportion / sum(proportion)
  frequency <- pmin(pmax(frequency, epsilon), 1 - epsilon)
  informative.samples <- rowSums(observed) > 0L
  state <- .sm_estep(x, frequency, proportion)
  history <- state$loglik
  converged <- FALSE
  for (iteration in seq_len(max.iter)) {
    tau <- state$membership
    proportion <- colMeans(tau[informative.samples, , drop = FALSE])
    if (any(proportion <= 0)) stop("Empty component; use a different initialization.")
    for (j in seq_len(k)) for (l in seq_len(nrow(frequency))) {
      weights <- tau[, j] * observed[, l]
      if (sum(weights) > 0) {
        frequency[l, j] <- sum(weights * state$dosage[, l, j]) / (2 * sum(weights))
      }
    }
    frequency <- pmin(pmax(frequency, epsilon), 1 - epsilon)
    next.state <- .sm_estep(x, frequency, proportion)
    change <- next.state$loglik - state$loglik
    if (change < -1e-7 * (1 + abs(state$loglik))) stop("EM likelihood decreased.")
    history <- c(history, next.state$loglik)
    state <- next.state
    if (abs(change) <= tolerance * (1 + abs(state$loglik))) {
      converged <- TRUE
      break
    }
  }
  list(frequency = frequency, proportion = proportion,
       membership = state$membership, loglik = state$loglik,
       history = history, converged = converged, iterations = iteration,
       informative = observed, epsilon = epsilon)
}
