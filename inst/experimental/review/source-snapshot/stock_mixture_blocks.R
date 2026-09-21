# Blockwise likelihood calculations--------------------------------------------
# Vectorize over one bounded marker block and one component at a time. Never
# allocate the full sample x marker x stock dosage array of the dense reference.
.sm_marker <- function(x, p) {
  joint <- sweep(x, 2L, stats::dbinom(0:2, 2, p, log = TRUE), "+")
  peak <- apply(joint, 1L, max)
  normalizer <- peak + log(rowSums(exp(joint - peak)))
  posterior <- exp(joint - normalizer)
  list(loglik = ifelse(rowSums(abs(x)) == 0, 0, normalizer),
       dosage = posterior[, 2L] + 2 * posterior[, 3L])
}

.sm_block_moments <- function(block, p) {
  n <- dim(block$log.likelihood)[1L]; m <- length(p)
  prior <- vapply(p, function(q) stats::dbinom(0:2, 2, q, log = TRUE), numeric(3))
  a <- sweep(matrix(block$log.likelihood[, 1L, ], n, m), 2L, prior[1L, ], "+")
  b <- sweep(matrix(block$log.likelihood[, 2L, ], n, m), 2L, prior[2L, ], "+")
  c <- sweep(matrix(block$log.likelihood[, 3L, ], n, m), 2L, prior[3L, ], "+")
  peak <- pmax(a, b, c)
  normalizer <- peak + log(exp(a - peak) + exp(b - peak) + exp(c - peak))
  dosage <- exp(b - normalizer) + 2 * exp(c - normalizer)
  normalizer[!block$informative] <- 0
  list(loglik = normalizer, dosage = dosage)
}

.sm_block_estep <- function(read.block, blocks, frequency, proportion, n,
                            sample.weights = rep(1, n)) {
  conditional <- matrix(0, n, ncol(frequency))
  for (index in blocks) {
    block <- read.block(index)
    for (j in seq_len(ncol(frequency))) {
      moments <- .sm_block_moments(block, frequency[index, j])
      conditional[, j] <- conditional[, j] + rowSums(moments$loglik)
    }
  }
  joint <- sweep(conditional, 2L, log(proportion), "+")
  peak <- apply(joint, 1L, max)
  normalizer <- peak + log(rowSums(exp(joint - peak)))
  list(membership = exp(joint - normalizer),
       loglik = sum(normalizer * sample.weights), sample.loglik = normalizer)
}

.sm_block_fit <- function(read.block, blocks, frequency, proportion, information,
                          max.iter, tolerance, epsilon,
                          sample.weights = rep(1, length(information)),
                          prior.center = rep(.5, nrow(frequency)), prior.strength = 0) {
  n <- length(information)
  penalty <- function(p) if (prior.strength == 0) 0 else
    prior.strength * sum(prior.center * log(p) + (1-prior.center) * log1p(-p))
  state <- .sm_block_estep(read.block, blocks, frequency, proportion, n, sample.weights)
  history <- state$loglik
  objective <- state$loglik + penalty(frequency)
  objective.history <- objective
  converged <- FALSE
  for (iteration in seq_len(max.iter)) {
    tau <- state$membership
    fit.weights <- sample.weights * (information > 0L)
    proportion <- colSums(tau * fit.weights) / sum(fit.weights)
    if (any(proportion <= 0)) stop("Empty component; use a different initialization.")
    # Re-read likelihoods to compute expected counts using OLD frequencies and
    # full-genome responsibilities. Updating membership per block is not EM.
    for (index in blocks) {
      block <- read.block(index)
      for (j in seq_len(ncol(frequency))) {
        moments <- .sm_block_moments(block, frequency[index, j])
        weights <- block$informative * (tau[, j] * sample.weights)
        counts <- colSums(weights)
        valid <- counts > 0 | prior.strength > 0
        update <- (colSums(weights * moments$dosage) + prior.strength * prior.center[index]) /
          (2 * counts + prior.strength)
        frequency[index[valid], j] <- pmin(1 - epsilon, pmax(epsilon, update[valid]))
      }
    }
    next.state <- .sm_block_estep(read.block, blocks, frequency, proportion, n, sample.weights)
    next.objective <- next.state$loglik + penalty(frequency)
    change <- next.objective - objective
    if (change < -1e-7 * (1 + abs(objective))) stop("EM objective decreased.")
    history <- c(history, next.state$loglik)
    objective.history <- c(objective.history, next.objective)
    objective <- next.objective
    state <- next.state
    if (abs(change) <= tolerance * (1 + abs(objective))) {
      converged <- TRUE
      break
    }
  }
  list(frequency = frequency, proportion = proportion,
       membership = state$membership, loglik = state$loglik, history = history,
       converged = converged, iterations = iteration, epsilon = epsilon,
       sample.loglik = state$sample.loglik, objective = objective,
       objective.history = objective.history, prior.center = prior.center,
       prior.strength = prior.strength)
}

# Training-only, bounded PCA seed space----------------------------------------
# Dosages and temporary completion are for initialization ONLY. Neither these
# completed values nor held-out samples enter the final likelihood estimator.
.sm_start_space <- function(read.block, m, information, sample.weights,
                            max.snps = 1000L, max.cells = 5e6) {
  train <- which(sample.weights > 0 & information > 0)
  count <- min(m, max.snps, floor(max.cells / length(train)))
  if (length(train) < 2L || count < 1L) stop("Insufficient initialization space.")
  index <- unique(as.integer(round(seq(1, m, length.out = count))))
  dosage <- matrix(0, length(train), length(index))
  for (l in seq_along(index)) {
    block <- read.block(index[l])
    x <- matrix(block$log.likelihood[train, , 1], length(train), 3)
    w <- sample.weights[train] * block$informative[train, 1]
    if (sum(w) == 0) stop("No training information for PCA seed marker.")
    p <- .5
    for (iteration in 1:1000) {
      moments <- .sm_marker(x, p)
      next.p <- min(1-1e-8, max(1e-8, sum(w*moments$dosage)/(2*sum(w))))
      if (abs(next.p-p) < 1e-7) break
      p <- next.p
    }
    if (iteration == 1000) stop("Pooled initialization frequency did not converge.")
    value <- .sm_marker(x, next.p)$dosage
    value[!block$informative[train, 1]] <- 2*next.p
    dosage[, l] <- value - stats::weighted.mean(value, sample.weights[train])
  }
  # Weighted covariance orientation; project unique training samples afterwards.
  rank <- min(10L, nrow(dosage)-1L, ncol(dosage))
  sv <- svd(dosage * sqrt(sample.weights[train]), nu = 0L, nv = rank)
  keep <- which(sv$d[seq_len(rank)] > max(sv$d)*1e-10)
  if (!length(keep)) stop("No variable PCA initialization coordinates.")
  list(train = train, scores = dosage %*% sv$v[, keep, drop = FALSE], markers = index)
}

.sm_partition_start <- function(read.block, blocks, base, information,
                                sample.weights, labels, epsilon) {
  n <- length(information); k <- ncol(base)
  weights <- matrix(.01/k, n, k)
  present <- which(!is.na(labels))
  weights[cbind(present, labels[present])] <- 1-.01+.01/k
  frequency <- base
  for (index in blocks) {
    block <- read.block(index)
    for (l in seq_along(index)) {
      dosage <- .sm_marker(matrix(block$log.likelihood[, , l], n, 3L), .5)$dosage
      for (j in seq_len(k)) {
        w <- weights[, j]*block$informative[, l]*sample.weights
        frequency[index[l], j] <- sum(w*dosage)/(2*sum(w))
      }
    }
  }
  list(frequency = pmin(pmax(frequency, epsilon), 1-epsilon),
       proportion = colSums(weights*(sample.weights*(information>0))) /
         sum(sample.weights*(information>0)))
}

# Independent soft partitions for additional starting values--------------------
# These weights initialize frequencies only. They are not fixed memberships,
# ancestry proportions or STRATA labels. No hard genotype calls are made.
.sm_block_start <- function(read.block, blocks, base, information, epsilon,
                            sample.weights = rep(1, length(information))) {
  n <- length(information); k <- ncol(base)
  weights <- matrix(stats::rexp(n * k), n, k)
  weights <- weights / rowSums(weights)
  start <- base
  for (index in blocks) {
    block <- read.block(index)
    for (l in seq_along(index)) {
      dosage <- .sm_marker(matrix(block$log.likelihood[, , l], n, 3L),
                           mean(base[index[l], ]))$dosage
      for (j in seq_len(k)) {
        w <- weights[, j] * block$informative[, l] * sample.weights
        start[index[l], j] <- sum(w * dosage) / (2 * sum(w))
      }
    }
  }
  list(frequency = pmin(pmax(start, epsilon), 1 - epsilon),
       proportion = colSums(weights * (sample.weights * (information > 0))) /
         sum(sample.weights * (information > 0)))
}

# Run sequentially: retain only the best fit, not a large fit per start----------
.sm_partition_agreement <- function(a, b) {
  keep <- !is.na(a) & !is.na(b)
  a <- a[keep]; b <- b[keep]; n <- length(a)
  if (n < 2L) return(c(ARI = NA_real_, N = n))
  choose2 <- function(x) x * (x - 1) / 2
  tab <- table(a, b)
  expected <- sum(choose2(rowSums(tab))) * sum(choose2(colSums(tab))) / choose2(n)
  upper <- (sum(choose2(rowSums(tab))) + sum(choose2(colSums(tab)))) / 2
  ari <- if (upper == expected) as.numeric(sum(choose2(tab)) == upper) else
    (sum(choose2(tab)) - expected) / (upper - expected)
  c(ARI = ari, N = n)
}

.sm_multistart <- function(read.block, blocks, frequency, proportion, information,
                           n.starts, seed, max.iter, tolerance, epsilon,
                           sample.weights = rep(1, length(information)),
                           start.strategy = "random", init.max.snps = 1000L,
                           init.max.cells = 5e6, prior.center = rep(.5, nrow(frequency)),
                           prior.strength = 0) {
  old.kind <- RNGkind()
  had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old.kind))
    if (had.seed) assign(".Random.seed", old.seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  rows <- vector("list", n.starts)
  # Retain compact partitions only, with the existing cell guard. Do not retain
  # all frequency/membership matrices merely to report start agreement.
  compare.rows <- which(information > 0 & sample.weights > 0)
  retain.partitions <- as.double(length(compare.rows)) * n.starts <= init.max.cells
  partitions <- vector("list", n.starts)
  best <- NULL; selected <- NA_integer_
  space <- NULL; lower <- NULL
  for (s in seq_len(n.starts)) {
    method <- if (s == 1L) "supplied" else if (start.strategy == "random") "random" else
      c("pca2", "pca5", "pca10", "random", "split", "perturb")[(s-2L) %% 6L + 1L]
    if (method == "split" && ncol(frequency) == 1L) method <- "random"
    attempt <- tryCatch({
      start <- if (s == 1L) list(frequency = frequency, proportion = proportion)
        else .sm_block_start(read.block, blocks, frequency, information, epsilon, sample.weights)
      k <- ncol(frequency)
      if (startsWith(method, "pca") || (method == "split" && k > 1L)) {
        if (is.null(space)) space <- .sm_start_space(read.block, nrow(frequency),
          information, sample.weights, init.max.snps, init.max.cells)
        scores <- space$scores; train <- space$train
        labels <- rep(NA_integer_, length(information))
        if (startsWith(method, "pca")) {
          pcs <- min(ncol(scores), as.integer(sub("pca", "", method)))
          labels[train] <- if (k == 1L) 1L else stats::kmeans(
            scores[, seq_len(pcs), drop = FALSE], k, nstart = 20L, iter.max = 100L)$cluster
        } else {
          if (is.null(lower)) lower <- .sm_multistart(read.block, blocks,
            frequency[, seq_len(k-1L), drop = FALSE], rep(1/(k-1L), k-1L),
            information, 5L, seed, max.iter, tolerance, epsilon, sample.weights,
            prior.center = prior.center, prior.strength = prior.strength)
          if (!lower$converged) stop("Lower-K split seed did not converge.")
          labels[train] <- max.col(lower$membership[train, , drop = FALSE], ties.method = "first")
          # Cycle through each lower-K component across successive split starts.
          target <- ((s-6L) %/% 6L) %% (k-1L) + 1L
          split.rows <- which(labels[train] == target)
          if (length(split.rows) < 2) stop("Split component contains fewer than two samples.")
          split <- stats::kmeans(scores[split.rows, , drop = FALSE], 2L, nstart = 20L)$cluster
          labels[train[split.rows[split == 2L]]] <- k
        }
        start <- .sm_partition_start(read.block, blocks, frequency, information,
          sample.weights, labels, epsilon)
      } else if (method == "perturb" && !is.null(best)) {
        start <- list(frequency = pmin(pmax(stats::plogis(stats::qlogis(best$frequency) +
          matrix(stats::rnorm(length(frequency), sd = .5), nrow(frequency))), epsilon), 1-epsilon),
          proportion = best$proportion)
      }
      .sm_block_fit(read.block, blocks, start$frequency, start$proportion,
                     information, max.iter, tolerance, epsilon, sample.weights,
                     prior.center, prior.strength)
    }, error = identity)
    if (inherits(attempt, "error")) {
      rows[[s]] <- tibble::tibble(START = s, METHOD = method, LOGLIK = NA_real_, OBJECTIVE = NA_real_,
        CONVERGED = FALSE, ITERATIONS = NA_integer_, ERROR = conditionMessage(attempt))
      next
    }
    rows[[s]] <- tibble::tibble(START = s, METHOD = method, LOGLIK = attempt$loglik, OBJECTIVE = attempt$objective,
      CONVERGED = attempt$converged, ITERATIONS = attempt$iterations, ERROR = NA_character_)
    if (retain.partitions) {
      tau <- attempt$membership[compare.rows, , drop = FALSE]
      label <- max.col(tau, ties.method = "first")
      label[rowSums(tau == apply(tau, 1, max)) > 1L] <- NA_integer_
      partitions[[s]] <- label
    }
    # Prefer converged fits; if none converge, return the best unfinished fit
    # with an explicit warning at the interface rather than hiding it.
    if (is.null(best) || (attempt$converged && !best$converged) ||
        (identical(attempt$converged, best$converged) && attempt$objective > best$objective)) {
      best <- attempt
      selected <- s
    }
  }
  diagnostics <- dplyr::bind_rows(rows)
  if (is.null(best)) stop("All starts failed: ", paste(unique(diagnostics$ERROR), collapse = "; "))
  diagnostics$SELECTED <- diagnostics$START == selected
  diagnostics$LOGLIK_MINUS_SELECTED <- diagnostics$LOGLIK - best$loglik
  diagnostics$OBJECTIVE_MINUS_SELECTED <- diagnostics$OBJECTIVE - best$objective
  agreement <- t(vapply(seq_len(n.starts), function(s) {
    if (!retain.partitions || is.null(partitions[[s]])) return(c(ARI = NA_real_, N = NA_real_))
    .sm_partition_agreement(partitions[[selected]], partitions[[s]])
  }, numeric(2)))
  diagnostics$PARTITION_ARI_VS_SELECTED <- agreement[, 1]
  diagnostics$N_COMPARABLE <- agreement[, 2]
  best$start.partition.assessment <- if (retain.partitions)
    "Label-invariant training partition agreement; ties excluded; inspect N_COMPARABLE" else
    "Not assessed: partition storage exceeds init.max.cells"
  if (any(!is.na(diagnostics$ERROR))) warning("One or more starts failed; inspect start.diagnostics.")
  if (!best$converged) {
    warning("No starts converged; returning an unfinished fit. Do not interpret it as final.")
  } else if (any(!diagnostics$CONVERGED & is.na(diagnostics$ERROR))) {
    warning("One or more starts did not converge; inspect start.diagnostics.")
  }
  best$start.diagnostics <- diagnostics
  best$selected.start <- selected
  best$seed <- seed
  best$start.strategy <- start.strategy
  best$init.max.snps <- init.max.snps
  best$init.max.cells <- init.max.cells
  best$near.best.starts <- sum(diagnostics$CONVERGED &
    is.finite(diagnostics$OBJECTIVE) &
    abs(diagnostics$OBJECTIVE-best$objective) <= tolerance*(1+abs(best$objective)))
  best
}
