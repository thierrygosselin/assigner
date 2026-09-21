# Diagnose initialization separately from held-out allele-frequency prediction.
# Run from the source directory; stockR's internal simulated example is required.
# Oracle labels are a diagnostic ONLY, never an unsupervised initialization.
source("R/stock_likelihood_mixture.R")
source("R/stock_mixture_blocks.R")
g <- get("myData", asNamespace("stockR"))[1:300, ]
truth <- attr(get("myData", asNamespace("stockR")), "grps")
print(table(truth))
n <- ncol(g); m <- nrow(g)
x <- array(0, c(n, 3L, m))
for (l in seq_len(m)) for (i in seq_len(n)) {
  if (!is.na(g[l, i])) x[i, , l] <- ifelse(0:2 == g[l, i], 0, -Inf)
}
prepared <- .sm_prepare(x, "log")
reader <- function(index) list(log.likelihood = prepared$log.likelihood[, , index, drop = FALSE],
                              informative = prepared$informative[, index, drop = FALSE])
blocks <- list(seq_len(m))
information <- rowSums(prepared$informative)
folds <- rep(1:3, length.out = n)
# Count private/rare alleles missed in each group's training observations.
for (fold in 1:3) for (group in sort(unique(truth))) {
  training <- g[, folds != fold & truth == group, drop = FALSE]
  testing <- g[, folds == fold & truth == group, drop = FALSE]
  p <- rowMeans(training, na.rm = TRUE) / 2
  absent <- sum((p == 0 & rowSums(testing > 0, na.rm = TRUE) > 0) |
                (p == 1 & rowSums(testing < 2, na.rm = TRUE) > 0), na.rm = TRUE)
  cat("Fold", fold, "group", group, "training n", ncol(training),
      "held-out n", ncol(testing), "markers with unseen allele", absent, "\n")
}
rows <- list(); z <- 0L
for (bound in c(1e-8, .001, .01)) for (fold in 1:3) for (k in 1:3) {
  train <- folds != fold
  for (method in c("random", "pca", if (k == 3L) "oracle")) {
    set.seed(5400 + fold * 10 + k)
    initial <- matrix(.5, m, k)
    if (method != "random") {
      # Training-only centering and temporary mean completion for starts only.
      # The actual fit still uses the full genotype likelihood objective.
      training <- t(g[, train, drop = FALSE])
      means <- colMeans(training, na.rm = TRUE)
      for (l in seq_len(m)) training[is.na(training[, l]), l] <- means[l]
      labels <- if (method == "oracle") as.integer(factor(truth[train])) else if (k == 1L)
        rep(1L, sum(train)) else stats::kmeans(stats::prcomp(training)$x[, 1:5],
                                             centers = k, nstart = 25)$cluster
      for (j in seq_len(k)) initial[, j] <- colMeans(training[labels == j, , drop = FALSE]) / 2
    }
    fit <- .sm_multistart(reader, blocks, pmin(pmax(initial, bound), 1-bound),
      rep(1/k, k), information, if (method == "random") 10L else 1L,
      54L + fold, 200L, 1e-8, bound, as.numeric(train))
    z <- z + 1L
    partition <- max.col(fit$membership[train, , drop = FALSE])
    rows[[z]] <- data.frame(epsilon = bound, fold, K = k, method,
      converged = fit$converged, training = fit$loglik,
      heldout = sum(fit$sample.loglik[!train]),
      truth.match = all(outer(partition, partition, "==") ==
                         outer(truth[train], truth[train], "==")))
    cat("Completed", z, "of 63 comparisons\n")
  }
}
result <- do.call(rbind, rows)
print(result, row.names = FALSE)
print(aggregate(cbind(training, heldout, truth.match) ~ epsilon + K + method,
                result, sum), row.names = FALSE)
stopifnot(all(result$converged))
