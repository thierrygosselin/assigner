#' Likelihood-aware discriminant analysis
#'
#' Estimate principal components directly from genotype likelihoods stored in
#' GDS, then discriminate groups in the retained PC space.
#'
#' @param data An open `SeqArray` GDS object or the path to a GDS file.
#' @param groups A grouping vector named with sample identifiers, a data frame
#'   containing `INDIVIDUALS` and `group.column`, or `NULL` to use individual
#'   metadata stored in the GDS. See **Strata and group labels**.
#' @param group.column Name of the strata column whose values define the groups
#'   to discriminate. The default is `"STRATA"`.
#' @param likelihood Genotype-likelihood field. `"auto"` prefers `PL`, then
#'   `GL`.
#' @param n.components Number of likelihood-aware principal components.
#' @param training.samples Optional sample identifiers used to estimate the
#'   PCA and group distributions. By default, samples with known group
#'   labels are used. Other samples are projected without refitting marker
#'   parameters.
#' @param group.prior Either `"equal"`, `"empirical"`, or a named numeric
#'   vector containing one positive prior probability per group.
#' @param max.iter Maximum number of alternating likelihood iterations.
#' @param tolerance Convergence tolerance for the reconstruction root mean
#'   squared error.
#' @param minimum.maf Minimum global minor-allele frequency estimated from the
#'   training likelihoods. Markers below this threshold are excluded before PCA.
#' @param covariance.ridge Positive diagonal regularizer used for group
#'   covariance matrices.
#' @param verbose Display progress messages wrapped to 80 characters.
#'
#' @return An object of class `assign_likelihood_dapc`. It contains sample scores,
#'   score covariance matrices, marker parameters, group parameters,
#'   assignment probabilities, likelihood diagnostics, sample and marker IDs,
#'   and the matched call.
#'
#' @details
#' This function does not convert uncertain genotypes into hard calls. It uses
#' an iterative genotype-likelihood PCA related to PCAngsd, followed by linear
#' discriminant analysis of the retained principal components.
#' For sample \eqn{i}, marker \eqn{l}, and genotype \eqn{g}, it maximizes a
#' posterior model of the form
#'
#' \deqn{P(D_{il} | z_i) = sum_g P(D_{il} | G_{il}=g)
#' P(G_{il}=g | z_i).}
#'
#' The genotype prior is diploid Hardy-Weinberg with an individual allele
#' frequency estimated iteratively from a truncated singular-value
#' decomposition. Posterior genotype dosages update those frequencies and the
#' final covariance diagonal uses posterior squared deviations, rather than the
#' square of a posterior mean. Eigendecomposition of this uncertainty-corrected
#' covariance produces the principal components. Posterior genotype variance is
#' propagated approximately into a covariance matrix for each sample's scores.
#'
#' `PL` and `GL` are likelihoods conditional on the sequencing data, not
#' posterior genotype probabilities. `PL` is converted with
#' \eqn{10^{-PL/10}} and `GL` with \eqn{10^{GL}} before combining it with the
#' model-derived genotype prior. Likelihood vectors are normalized for
#' numerical stability without changing their relative information.
#'
#' PCA is fitted without using group labels. Group labels enter only after PCA,
#' when group means and a pooled within-group
#' covariance are estimated, following linear discriminant analysis. Assignment
#' integrates each sample's approximate PC-score covariance into the group
#' density. Returned probabilities are model-based assignments, not genomic
#' ancestry proportions.
#'
#' If every labelled sample is used for training, assignments for those samples
#' are in-sample descriptions and are not measures of predictive accuracy.
#' Proper validation requires externally defined training and test samples, or
#' repeated refitting within resampling folds. This function never selects
#' dimensions or tuning parameters using test labels.
#'
#' The initial implementation is intentionally restricted to diploid,
#' biallelic variants having exactly three complete `PL` or `GL` values. Fully
#' missing likelihood vectors are allowed and contribute no likelihood.
#' Partially missing or malformed vectors are rejected. Multiallelic, haploid,
#' polyploid, phased, pooled, and gVCF reference-block models are not currently
#' supported.
#'
#' Technical effects in depth, library preparation, plate, lane, mapping, or
#' variant calling can structure genotype likelihoods. A likelihood-aware model
#' preserves uncertainty but does not remove batch effects. Inspect likelihood
#' entropy, missing likelihoods, depth, and inferred coordinates against all
#' relevant technical and biological metadata.
#'
#' @section Strata and group labels:
#' In the assigner ecosystem, **strata** means the complete sample metadata table:
#' one row per individual and potentially many biological or technical columns,
#' such as population, sampling site, sex, year, plate, library, or sequencing
#' lane. A **group label** is the value from one selected strata column that is
#' used as the response in the discriminant step.
#'
#' With `groups = NULL`, the function reads the individual metadata stored in
#' the GDS and uses `group.column = "STRATA"` by default. Supplying a metadata
#' data frame uses the same rule. A named vector is a shortcut containing the
#' selected group labels directly; its names must be the GDS sample IDs.
#'
#' Selecting a column as the grouping variable does not adjust for the other
#' strata columns. For example, choosing `STRATA` as `group.column` makes its
#' biological strata the classification target; `PLATE` and `LIBRARY` remain audit
#' variables and are not silently included as covariates. Associations between
#' inferred components, assignments, and those technical variables must be
#' examined explicitly.
#'
#' @section Lifecycle:
#' `assign_likelihood_dapc()` is experimental. It is available for use, but its
#' interface, numerical implementation, and returned components may change as
#' the method is evaluated on additional datasets. Analyses should record the
#' assigner package version and, when using a development build, the exact Git
#' commit. See `vignette("life_cycle", package = "assigner")`.
#'
#' @section Choosing an assignment method:
#' [assign_individuals()] evaluates the probability of an individual's data
#' directly under allele frequencies estimated for each candidate source. It is
#' the preferred generative assignment method when source populations are
#' predefined and its assumptions are defensible. `assign_likelihood_dapc()`
#' instead provides likelihood-aware genomic coordinates followed by
#' discriminant classification. It is useful when visualization and
#' dimension-reduced discrimination are part of the question. Its assignment
#' probabilities are classifier outputs and are not interchangeable with the
#' population likelihoods returned by [assign_individuals()].
#'
#' @references
#' Patterson, N., Price, A. L., & Reich, D. (2006). Population structure and
#' eigenanalysis. *PLoS Genetics*, 2, e190.
#'
#' Jombart, T., Devillard, S., & Balloux, F. (2010). Discriminant analysis of
#' principal components: a new method for the analysis of genetically
#' structured populations. *BMC Genetics*, 11, 94.
#'
#' Meisner, J., & Albrechtsen, A. (2018). Inferring population structure and
#' admixture proportions in low-depth NGS data. *Genetics*, 210, 719-731.
#'
#' Meisner, J., & Albrechtsen, A. (2019). Testing for Hardy-Weinberg equilibrium
#' in structured populations using genotype or low-depth next generation
#' sequencing data. *Molecular Ecology Resources*, 19, 1144-1152.
#'
#' @examples
#' \dontrun{
#' fit <- assign_likelihood_dapc(
#'   data = "genotypes.gds",
#'   group.column = "STRATA",
#'   likelihood = "PL",
#'   n.components = 3
#' )
#' fit$assignments
#' }
#'
#' @author Thierry Gosselin \email{thierry.gosselin@@csiro.au}
#' @export
assign_likelihood_dapc <- function(
    data,
    groups = NULL,
    group.column = "STRATA",
    likelihood = c("auto", "PL", "GL"),
    n.components = 2L,
    training.samples = NULL,
    group.prior = c("equal", "empirical"),
    max.iter = 100L,
    tolerance = 1e-5,
    minimum.maf = 0.01,
    covariance.ridge = 1e-6,
    verbose = TRUE
) {
  likelihood <- match.arg(likelihood)
  .ld_check_count(n.components, "n.components", 1L)
  .ld_check_count(max.iter, "max.iter", 1L)
  .ld_check_number(tolerance, "tolerance", lower = 0, open = TRUE)
  .ld_check_number(minimum.maf, "minimum.maf", lower = 0)
  if (minimum.maf >= 0.5) rlang::abort("`minimum.maf` must be smaller than 0.5.")
  .ld_check_number(covariance.ridge, "covariance.ridge", lower = 0,
                   open = TRUE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    rlang::abort("`verbose` must be TRUE or FALSE.")
  }
  context <- .ld_open_gds(data)
  on.exit(.ld_close_gds(context), add = TRUE)
  gds <- context$gds
  sample.ids <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  variant.ids <- SeqArray::seqGetData(gds, "variant.id")
  if (length(sample.ids) < 3L) rlang::abort("At least three samples are required.")
  if (anyDuplicated(sample.ids)) rlang::abort("Active GDS sample IDs must be unique.")

  alleles <- as.character(SeqArray::seqGetData(gds, "allele"))
  if (length(alleles) != length(variant.ids) ||
      any(lengths(strsplit(alleles, ",", fixed = TRUE)) != 2L)) {
    rlang::abort("All active variants must be biallelic.")
  }
  field <- .ld_choose_likelihood(gds, likelihood)
  log.likelihood <- .ld_read_likelihood(gds, field, length(sample.ids),
                                        length(variant.ids))
  observed <- apply(log.likelihood, c(1L, 3L), function(x) all(is.finite(x)))
  sample.information <- rowSums(observed)
  marker.information <- colSums(observed)
  if (any(sample.information == 0L)) {
    rlang::abort(paste0("Samples without observed genotype likelihoods: ",
      paste(sample.ids[sample.information == 0L], collapse = ", "), "."))
  }
  if (any(marker.information == 0L)) {
    rlang::abort("Every active marker must have at least one likelihood vector.")
  }

  group <- .ld_groups(gds, groups, group.column, sample.ids)
  if (is.null(training.samples)) {
    training <- !is.na(group)
  } else {
    if (!is.character(training.samples) || anyNA(training.samples) ||
        anyDuplicated(training.samples) ||
        !all(training.samples %in% sample.ids)) {
      rlang::abort("`training.samples` must be unique active GDS sample IDs.")
    }
    training <- sample.ids %in% training.samples
  }
  if (sum(training) <= n.components + 1L) {
    rlang::abort("Training requires more samples than `n.components + 1`.")
  }
  if (anyNA(group[training])) {
    rlang::abort("Every training sample must have a group label.")
  }
  group.training <- droplevels(factor(group[training]))
  if (nlevels(group.training) < 2L) {
    rlang::abort("At least two training groups are required.")
  }
  if (n.components >= min(sum(training), length(variant.ids))) {
    rlang::abort("`n.components` must be smaller than both training samples and markers.")
  }

  # Estimate the ALT-allele frequency without hard-calling genotypes. The
  # symmetric interval below turns the ALT frequency into an implicit MAF
  # filter while preserving the REF/ALT orientation used by the GDS.
  frequency <- .ld_allele_frequency(
    log.likelihood[training, , , drop = FALSE]
  )
  keep.variant <- frequency >= minimum.maf & frequency <= 1 - minimum.maf
  if (sum(keep.variant) <= n.components) {
    rlang::abort("Too few markers remain after the likelihood-based MAF filter.")
  }
  excluded.variant.ids <- variant.ids[!keep.variant]
  variant.ids <- variant.ids[keep.variant]
  alleles <- alleles[keep.variant]
  frequency <- frequency[keep.variant]
  log.likelihood <- log.likelihood[, , keep.variant, drop = FALSE]
  marker.information <- marker.information[keep.variant]

  .ld_message(verbose, "Fitting ", n.components, " likelihood-aware PC(s) ",
              "from ", sum(training), " training samples and ",
              length(variant.ids), " biallelic variants using ", field, ".")
  model <- .ld_fit_pca(log.likelihood[training, , , drop = FALSE],
                       frequency, n.components, max.iter, tolerance, verbose)

  scores <- matrix(NA_real_, nrow = length(sample.ids), ncol = n.components,
                   dimnames = list(sample.ids, paste0("PC", seq_len(n.components))))
  score.covariance <- array(NA_real_,
    dim = c(n.components, n.components, length(sample.ids)),
    dimnames = list(colnames(scores), colnames(scores), sample.ids))
  scores[training, ] <- model$scores
  score.covariance[, , training] <- model$score.covariance
  if (any(!training)) {
    .ld_message(verbose, "Projecting ", sum(!training),
                " non-training sample(s) without refitting marker parameters.")
    projected <- .ld_project_pca(log.likelihood[!training, , , drop = FALSE],
                                 model, max.iter, tolerance)
    scores[!training, ] <- projected$scores
    score.covariance[, , !training] <- projected$score.covariance
  }

  priors <- .ld_group_prior(group.training, group.prior)
  classifier <- .ld_fit_groups(scores[training, , drop = FALSE], group.training,
                               covariance.ridge)
  probability <- .ld_assign(scores, score.covariance, classifier, priors)
  assignments <- tibble::tibble(
    INDIVIDUALS = sample.ids,
    OBSERVED_GROUP = as.character(group),
    TRAINING = training,
    ASSIGNED_GROUP = colnames(probability)[max.col(probability, ties.method = "first")],
    MAX_PROBABILITY = apply(probability, 1L, max)
  )
  probability.table <- tibble::as_tibble(probability, rownames = "INDIVIDUALS")

  result <- list(
    call = match.call(),
    likelihood = field,
    sample.ids = sample.ids,
    variant.ids = variant.ids,
    excluded.variant.ids = excluded.variant.ids,
    alleles = alleles,
    training = training,
    scores = scores,
    score.covariance = score.covariance,
    allele.frequency = stats::setNames(model$allele.frequency, variant.ids),
    eigenvalues = model$eigenvalues,
    genomic.covariance = model$genomic.covariance,
    marker.loadings = structure(model$loadings,
      dimnames = list(as.character(variant.ids), colnames(scores))),
    group.parameters = classifier,
    group.prior = priors,
    probability = probability,
    assignments = assignments,
    probability.table = probability.table,
    diagnostics = list(
      converged = model$converged,
      iterations = model$iterations,
      rmse = model$rmse,
      sample.observed.likelihoods = stats::setNames(sample.information, sample.ids),
      marker.observed.likelihoods = stats::setNames(marker.information, variant.ids)
    )
  )
  class(result) <- "assign_likelihood_dapc"
  if (!model$converged) {
    warning("Likelihood-aware PCA reached `max.iter` without convergence.")
  }
  .ld_message(verbose, "Likelihood-aware grouping completed. Training assignments ",
              "are descriptive unless assessed by external validation.")
  result
}

.ld_open_gds <- function(data) {
  owned <- !inherits(data, "SeqVarGDSClass")
  if (owned && (!is.character(data) || length(data) != 1L || is.na(data) ||
      !file.exists(data))) {
    rlang::abort("`data` must be an open SeqArray GDS or an existing GDS path.")
  }
  gds <- if (owned) SeqArray::seqOpen(data, readonly = TRUE) else data
  SeqArray::seqFilterPush(gds)
  list(gds = gds, owned = owned)
}

.ld_close_gds <- function(context) {
  SeqArray::seqFilterPop(context$gds)
  if (context$owned) SeqArray::seqClose(context$gds)
}

.ld_choose_likelihood <- function(gds, requested) {
  present <- vapply(c("PL", "GL"), function(field) {
    !is.null(gdsfmt::index.gdsn(gds, paste0("annotation/format/", field),
                               silent = TRUE))
  }, logical(1))
  if (requested == "auto") {
    if (!any(present)) rlang::abort("The GDS contains neither PL nor GL.")
    return(names(present)[which(present)[1L]])
  }
  if (!present[[requested]]) rlang::abort(paste0("The GDS does not contain ", requested, "."))
  requested
}

.ld_read_likelihood <- function(gds, field, n.sample, n.variant) {
  raw <- SeqArray::seqGetData(gds, paste0("annotation/format/", field))
  if (!is.list(raw) || !all(c("length", "data") %in% names(raw)) ||
      length(raw$length) != n.variant || any(raw$length != 3L)) {
    rlang::abort(paste0(field,
      " must contain exactly three values for every active biallelic variant."))
  }
  value <- raw$data
  if (!is.matrix(value) || !identical(dim(value), c(n.sample, 3L * n.variant))) {
    rlang::abort(paste0("Unexpected SeqArray layout for ", field, "."))
  }
  value <- array(as.numeric(value), dim = c(n.sample, 3L, n.variant))
  missing <- apply(value, c(1L, 3L), function(x) all(is.na(x)))
  partial <- apply(value, c(1L, 3L), function(x) anyNA(x) && !all(is.na(x)))
  if (any(partial)) rlang::abort(paste0("Partially missing ", field,
                                       " vectors are not supported."))
  if (field == "PL") {
    if (any(value[!is.na(value)] < 0) ||
        any(value[!is.na(value)] != floor(value[!is.na(value)]))) {
      rlang::abort("PL values must be non-negative whole numbers.")
    }
    value <- -value * log(10) / 10
  } else {
    if (any(!is.finite(value[!is.na(value)]))) {
      rlang::abort("GL values must be finite log10 likelihoods.")
    }
    value <- value * log(10)
  }
  for (i in seq_len(n.sample)) for (l in seq_len(n.variant)) {
    if (!missing[i, l]) value[i, , l] <- value[i, , l] - max(value[i, , l])
  }
  value
}

.ld_groups <- function(gds, groups, column, sample.ids) {
  if (!is.character(column) || length(column) != 1L || is.na(column) ||
      !nzchar(column)) rlang::abort("`group.column` must name one metadata column.")
  if (is.null(groups)) {
    groups <- genometranslator::extract_individuals_metadata(gds, whitelist = FALSE)
  }
  if (is.atomic(groups) && !is.data.frame(groups)) {
    if (is.null(names(groups)) || anyDuplicated(names(groups)) ||
        !all(sample.ids %in% names(groups))) {
      rlang::abort("A grouping vector must be named with all active sample IDs.")
    }
    answer <- as.character(groups[sample.ids])
  } else {
    if (!is.data.frame(groups) ||
        !all(c("INDIVIDUALS", column) %in% names(groups))) {
      rlang::abort("Grouping metadata requires INDIVIDUALS and `group.column`.")
    }
    ids <- as.character(groups$INDIVIDUALS)
    if (anyNA(ids) || anyDuplicated(ids)) {
      rlang::abort("Grouping metadata sample IDs must be unique and non-missing.")
    }
    rows <- match(sample.ids, ids)
    if (anyNA(rows)) rlang::abort("Grouping metadata must cover all active samples.")
    answer <- as.character(groups[[column]][rows])
  }
  answer[!nzchar(trimws(answer))] <- NA_character_
  factor(answer, levels = unique(answer[!is.na(answer)]))
}

.ld_logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(x - m)))
}

.ld_log_prior <- function(p) {
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  cbind(2 * log1p(-p), log(2) + log(p) + log1p(-p), 2 * log(p))
}

.ld_posterior_moments <- function(log.likelihood, probability) {
  n <- dim(log.likelihood)[1L]
  loci <- dim(log.likelihood)[3L]
  mean <- matrix(NA_real_, n, loci)
  variance <- matrix(NA_real_, n, loci)

  # These loops operate on many very small three-genotype vectors. Keeping the
  # calculation in base R avoids repeatedly creating tibbles in a hot path.
  for (l in seq_len(loci)) for (i in seq_len(n)) {
    ll <- log.likelihood[i, , l]
    if (all(is.na(ll))) ll <- rep(0, 3L)
    prior <- .ld_log_prior(probability[i, l])
    posterior <- exp(ll + prior - .ld_logsumexp(ll + prior))
    mean[i, l] <- sum(0:2 * posterior)
    variance[i, l] <- sum((0:2)^2 * posterior) - mean[i, l]^2
  }
  list(mean = mean, variance = variance)
}

.ld_allele_frequency <- function(log.likelihood, max.iter = 100L,
                                 tolerance = 1e-7) {
  loci <- dim(log.likelihood)[3L]
  frequency <- rep(0.25, loci)

  # EM update for a diploid biallelic locus under Hardy-Weinberg proportions:
  # calculate posterior ALT dosages from the likelihoods and current allele
  # frequency, then divide their sum by twice the informative sample count.
  for (iteration in seq_len(max.iter)) {
    previous <- frequency
    for (l in seq_len(loci)) {
      observed <- which(apply(log.likelihood[, , l, drop = FALSE], 1L,
                              function(x) all(is.finite(x))))
      expected.alt.dosage <- numeric(length(observed))
      for (j in seq_along(observed)) {
        i <- observed[j]
        prior <- .ld_log_prior(frequency[l])
        posterior <- exp(log.likelihood[i, , l] + prior -
                         .ld_logsumexp(log.likelihood[i, , l] + prior))
        expected.alt.dosage[j] <- sum(0:2 * posterior)
      }
      frequency[l] <- mean(expected.alt.dosage) / 2
    }

    # A small numerical boundary prevents log(0); it is not a biological
    # pseudocount. Near-fixed loci are subsequently governed by minimum.maf.
    frequency <- pmin(pmax(frequency, 1e-5), 1 - 1e-5)
    if (sqrt(mean((frequency - previous)^2)) < tolerance) break
  }
  frequency
}

.ld_truncated_svd <- function(centered, components) {
  decomposition <- base::svd(centered, nu = components, nv = components)
  list(
    reconstruction = decomposition$u[, seq_len(components), drop = FALSE] %*%
      (diag(decomposition$d[seq_len(components)], components) %*%
       t(decomposition$v[, seq_len(components), drop = FALSE])),
    marker.axes = decomposition$u[, seq_len(components), drop = FALSE],
    sample.axes = decomposition$v[, seq_len(components), drop = FALSE],
    singular.values = decomposition$d[seq_len(components)]
  )
}

.ld_fit_pca <- function(log.likelihood, frequency, components, max.iter, tolerance,
                        verbose) {
  n <- dim(log.likelihood)[1L]
  loci <- dim(log.likelihood)[3L]
  population.probability <- matrix(frequency, n, loci, byrow = TRUE)
  moments <- .ld_posterior_moments(log.likelihood, population.probability)
  centered <- sweep(moments$mean, 2L, 2 * frequency, "-")
  centered[is.na(centered)] <- 0
  unscaled <- .ld_truncated_svd(t(centered), components)
  individual.frequency <- pmin(pmax(
    sweep(unscaled$reconstruction, 1L, 2 * frequency, "+") / 2,
    1e-4), 1 - 1e-4)
  history <- numeric(max.iter)
  converged <- FALSE
  for (iteration in seq_len(max.iter)) {
    previous <- unscaled$reconstruction
    moments <- .ld_posterior_moments(log.likelihood,
                                     t(individual.frequency))
    centered <- sweep(moments$mean, 2L, 2 * frequency, "-")
    centered[is.na(centered)] <- 0
    unscaled <- .ld_truncated_svd(t(centered), components)
    individual.frequency <- pmin(pmax(
      sweep(unscaled$reconstruction, 1L, 2 * frequency, "+") / 2,
      1e-4), 1 - 1e-4)
    history[iteration] <- sqrt(mean((unscaled$reconstruction - previous)^2))
    .ld_message(verbose, "Iteration ", iteration,
                ": PC loading RMSE = ",
                format(history[iteration], digits = 7), ".")
    if (history[iteration] < tolerance) {
      converged <- TRUE
      break
    }
  }

  moments <- .ld_posterior_moments(log.likelihood,
                                   t(individual.frequency))
  scale <- sqrt(2 * frequency * (1 - frequency))
  standardized <- sweep(
    sweep(moments$mean, 2L, 2 * frequency, "-"), 2L, scale, "/")
  standardized[is.na(standardized)] <- 0
  genomic.covariance <- tcrossprod(standardized) / loci
  diagonal <- rowSums(sweep(
    moments$variance + sweep(moments$mean, 2L, 2 * frequency, "-")^2,
    2L, scale^2, "/"), na.rm = TRUE) / loci
  diag(genomic.covariance) <- diagonal
  eig <- eigen(genomic.covariance, symmetric = TRUE)
  values <- pmax(eig$values[seq_len(components)], 0)
  vectors <- eig$vectors[, seq_len(components), drop = FALSE]
  scores <- sweep(vectors, 2L, sqrt(values * loci), "*")
  colnames(scores) <- paste0("PC", seq_len(components))
  marker.axes <- crossprod(standardized, vectors)
  denominator <- sqrt(pmax(values * loci, .Machine$double.eps))
  marker.axes <- sweep(marker.axes, 2L, denominator, "/")
  colnames(marker.axes) <- colnames(scores)
  score.covariance <- .ld_score_covariance(moments$variance, marker.axes,
                                           scale)
  list(
    allele.frequency = frequency,
    individual.frequency = individual.frequency,
    unscaled.marker.axes = unscaled$marker.axes,
    scores = scores,
    score.covariance = score.covariance,
    loadings = marker.axes,
    eigenvalues = values,
    genomic.covariance = genomic.covariance,
    rmse = history[seq_len(iteration)],
    iterations = iteration,
    converged = converged
  )
}

.ld_score_covariance <- function(genotype.variance, marker.axes, scale) {
  n <- nrow(genotype.variance)
  components <- ncol(marker.axes)
  answer <- array(0, c(components, components, n))
  for (i in seq_len(n)) {
    weight <- genotype.variance[i, ] / scale^2
    weight[!is.finite(weight)] <- 0
    answer[, , i] <- crossprod(marker.axes,
      sweep(marker.axes, 1L, weight, "*"))
  }
  answer
}

.ld_project_pca <- function(log.likelihood, model, max.iter, tolerance) {
  n <- dim(log.likelihood)[1L]
  loci <- dim(log.likelihood)[3L]
  frequency <- model$allele.frequency
  individual.frequency <- matrix(frequency, loci, n)
  coefficients <- matrix(0, ncol(model$unscaled.marker.axes), n)
  for (iteration in seq_len(max.iter)) {
    previous <- coefficients
    moments <- .ld_posterior_moments(log.likelihood,
                                     t(individual.frequency))
    centered <- sweep(moments$mean, 2L, 2 * frequency, "-")
    centered[is.na(centered)] <- 0
    coefficients <- crossprod(model$unscaled.marker.axes, t(centered))
    reconstruction <- model$unscaled.marker.axes %*% coefficients
    individual.frequency <- pmin(pmax(
      sweep(reconstruction, 1L, 2 * frequency, "+") / 2,
      1e-4), 1 - 1e-4)
    if (sqrt(mean((coefficients - previous)^2)) < tolerance) break
  }
  moments <- .ld_posterior_moments(log.likelihood,
                                   t(individual.frequency))
  scale <- sqrt(2 * frequency * (1 - frequency))
  standardized <- sweep(
    sweep(moments$mean, 2L, 2 * frequency, "-"), 2L, scale, "/")
  standardized[is.na(standardized)] <- 0
  scores <- standardized %*% model$loadings
  colnames(scores) <- colnames(model$scores)
  list(scores = scores,
       score.covariance = .ld_score_covariance(moments$variance,
                                               model$loadings, scale))
}

.ld_fit_groups <- function(scores, group, covariance.ridge) {
  levels <- levels(group); components <- ncol(scores)
  means <- matrix(NA_real_, length(levels), components,
                  dimnames = list(levels, colnames(scores)))
  for (k in seq_along(levels)) {
    values <- scores[group == levels[k], , drop = FALSE]
    means[k, ] <- colMeans(values)
  }
  residual <- scores - means[as.integer(group), , drop = FALSE]
  denominator <- nrow(scores) - length(levels)
  if (denominator <= 0L) {
    rlang::abort("Training requires residual degrees of freedom within groups.")
  }
  covariance <- crossprod(residual) / denominator +
    diag(covariance.ridge, components)
  list(mean = means, covariance = covariance)
}

.ld_group_prior <- function(group, prior) {
  levels <- levels(group)
  if (is.character(prior)) {
    prior <- match.arg(prior, c("equal", "empirical"))
    value <- if (prior == "equal") rep(1 / length(levels), length(levels)) else
      as.numeric(table(group)) / length(group)
    return(stats::setNames(value, levels))
  }
  if (!is.numeric(prior) || is.null(names(prior)) || any(!is.finite(prior)) ||
      any(prior <= 0) || !setequal(names(prior), levels)) {
    rlang::abort("A numeric `group.prior` must name one positive value per group.")
  }
  prior <- prior[levels] / sum(prior[levels])
  prior
}

.ld_assign <- function(scores, score.covariance, classifier, prior) {
  groups <- rownames(classifier$mean)
  log.probability <- matrix(NA_real_, nrow(scores), length(groups),
                            dimnames = list(rownames(scores), groups))
  for (i in seq_len(nrow(scores))) for (k in seq_along(groups)) {
    covariance <- classifier$covariance + score.covariance[, , i]
    difference <- scores[i, ] - classifier$mean[k, ]
    root <- chol(covariance)
    quadratic <- sum(backsolve(root, difference, transpose = TRUE)^2)
    log.det <- 2 * sum(log(diag(root)))
    log.probability[i, k] <- log(prior[k]) - 0.5 *
      (length(difference) * log(2 * pi) + log.det + quadratic)
  }
  probability <- log.probability
  for (i in seq_len(nrow(probability))) {
    normalizer <- .ld_logsumexp(probability[i, ])
    probability[i, ] <- exp(probability[i, ] - normalizer)
  }
  probability
}

.ld_check_count <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < minimum ||
      x != floor(x)) rlang::abort(paste0("`", name, "` must be a whole number >= ", minimum, "."))
}

.ld_check_number <- function(x, name, lower, open = FALSE) {
  invalid <- !is.numeric(x) || length(x) != 1L || !is.finite(x) ||
    if (open) x <= lower else x < lower
  if (invalid) rlang::abort(paste0("`", name, "` must be a finite number ",
    if (open) "> " else ">= ", lower, "."))
}

.ld_message <- function(verbose, ...) {
  if (!verbose) return(invisible(NULL))
  message(paste(strwrap(paste0(...), width = 80L), collapse = "\n"))
  invisible(NULL)
}
