#' Hierarchical AMOVA for incomplete genomic data
#'
#' Computes a locus-wise analysis of molecular variance (AMOVA). Missing
#' genotypes are handled independently at each locus and all sample sizes,
#' degrees of freedom, and unequal-sample-size coefficients are recalculated
#' from the observations that actually contribute to that locus. This follows
#' the locus-wise strategy used by Stacks rather than deleting every individual
#' with at least one missing genotype.
#'
#' @param data A GDS file path or open GDS object accepted by
#'   [genometranslator::read_genome()]. A long genomic data frame is also
#'   accepted for programmatic use and testing.
#' @param strata Optional strata file or data frame accepted by
#'   [genometranslator::read_strata()]. It must contain `INDIVIDUALS` and the
#'   columns named in `hierarchy`. When supplied, these hierarchy columns take
#'   precedence over metadata stored in the GDS.
#' @param hierarchy Character vector naming nested grouping columns, ordered
#'   from highest to lowest level (for example `c("REGION", "POP_ID")`).
#' @param individual,marker Column names containing individual and marker IDs.
#' @param value Optional genotype-value column. With `value = NULL`, the
#'   function derives diploid alternate-allele dosage from `ALT_DOSAGE`, numeric
#'   `GT`, or six-digit calibrated `GT`. Supply a column name explicitly for
#'   haplotypes, nucleotide sequences, or custom distances.
#' @param distance Molecular distance. `"identity"` is the 0/1 allelic or
#'   haplotypic distance used for a Meirmans-style standardized statistic;
#'   `"nucleotide"` is the proportion of nucleotide differences between
#'   haplotype strings; `"euclidean"` uses squared dosage differences; and
#'   `"manhattan"` uses absolute dosage differences. A function may also be
#'   supplied; it must return a squared-distance matrix.
#' @param missing Missing-data strategy. `"locuswise"` uses all called
#'   observations at each locus. `"filter"` first applies marker-level call-rate
#'   thresholds and then works locus-wise. `"complete"` uses only individuals
#'   called at every retained marker. No imputation is performed.
#' @param min_call_rate Minimum called proportion required in every lowest-level
#'   group when `missing = "filter"`.
#' @param min_groups Minimum number of represented lowest-level groups required
#'   at a locus.
#' @param min_individuals Minimum called individuals required per represented
#'   lowest-level group.
#' @param euclidean How to handle a non-Euclidean squared-distance matrix.
#'   `"check"` records the result, `"lingoes"` applies a Lingoes correction,
#'   and `"none"` skips validation.
#' @param standardized Calculate a Meirmans-style standardized statistic. This
#'   is available for `distance = "identity"`, where the theoretical maximum is
#'   defined. It is deliberately not obtained by dividing distances by their
#'   observed sample maximum.
#' @param permutations Number of hierarchy-aware randomizations. Zero disables
#'   permutation tests.
#' @param seed Optional random seed used for permutations.
#' @param tolerance Numerical tolerance for Euclidean checks.
#'
#' @return An object of class `assigner_amova` containing global statistics,
#'   locus variance components, missing-data diagnostics, and permutation
#'   results when requested.
#' @export
#' @examples
#' \dontrun{
#' # Hierarchy stored in GDS sample metadata:
#' fit <- amova_genomic(
#'   data = "genomic_data.gds",
#'   hierarchy = c("REGION", "POP_ID"),
#'   distance = "euclidean",
#'   missing = "locuswise",
#'   standardized = FALSE
#' )
#'
#' # Hierarchy stored separately:
#' fit <- amova_genomic(
#'   data = "genomic_data.gds",
#'   strata = "amova_strata.tsv",
#'   hierarchy = c("REGION", "POP_ID")
#' )
#' }
#'
#' @references
#' Excoffier L, Smouse PE, Quattro JM (1992). Analysis of molecular
#' variance inferred from metric distances among DNA haplotypes. Genetics 131,
#' 479-491.
#'
#' Meirmans PG (2006). Using the AMOVA framework to estimate a standardized
#' genetic differentiation measure. Evolution 60, 2399-2402.
#' @md
amova_genomic <- function(
  data,
  hierarchy,
  strata = NULL,
  individual = "INDIVIDUALS",
  marker = "MARKERS",
  value = NULL,
  distance = c("euclidean", "identity", "nucleotide", "manhattan"),
  missing = c("locuswise", "filter", "complete"),
  min_call_rate = 0.8,
  min_groups = 2L,
  min_individuals = 2L,
  euclidean = c("check", "lingoes", "none"),
  standardized = identical(distance, "identity"),
  permutations = 0L,
  seed = NULL,
  tolerance = sqrt(.Machine$double.eps)
) {
  if (!is.character(hierarchy) || !length(hierarchy)) {
    stop("`hierarchy` must name at least one grouping column.", call. = FALSE)
  }
  if (is.data.frame(data)) {
    data <- as.data.frame(data)
  } else {
    data <- genometranslator::read_genome(data = data, import.metadata = TRUE)
    if (!is.data.frame(data)) {
      stop("`genometranslator::read_genome()` did not return a genomic data frame.", call. = FALSE)
    }
    missing_hierarchy <- setdiff(hierarchy, names(data))
    if (length(hierarchy) == 1L && length(missing_hierarchy) == 1L &&
        "STRATA" %in% names(data)) {
      data[[hierarchy]] <- data$STRATA
    }
  }
  if (!is.null(strata)) {
    strata_data <- genometranslator::read_strata(strata = strata)$strata
    strata_needed <- unique(c(individual, hierarchy))
    strata_absent <- setdiff(strata_needed, names(strata_data))
    if (length(strata_absent)) {
      stop("Missing strata columns: ", paste(strata_absent, collapse = ", "), call. = FALSE)
    }
    data <- data[setdiff(names(data), hierarchy)]
    data <- dplyr::left_join(
      data,
      dplyr::distinct(strata_data[strata_needed]),
      by = individual
    )
  }
  prepared <- .amova_prepare_value(data, marker = marker, value = value,
                                   distance = distance)
  data <- prepared$data
  value <- prepared$value
  needed <- unique(c(individual, marker, value, hierarchy))
  absent <- setdiff(needed, names(data))
  if (length(absent)) {
    stop("Missing columns: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  missing <- match.arg(missing)
  euclidean <- match.arg(euclidean)
  if (is.character(distance)) distance <- match.arg(distance)
  if (!is.function(distance) && standardized && distance != "identity") {
    stop("`standardized = TRUE` currently requires `distance = \"identity\"`; ",
         "a sample-normalized distance is not a theoretical maximum.", call. = FALSE)
  }
  permutations <- as.integer(permutations)
  if (is.na(permutations) || permutations < 0L) {
    stop("`permutations` must be a non-negative integer.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  x <- data[needed]
  names(x)[match(c(individual, marker, value), names(x))] <-
    c(".individual", ".marker", ".value")
  if (anyDuplicated(x[c(".individual", ".marker")])) {
    stop("Each individual-marker combination must occur once.", call. = FALSE)
  }
  if (anyNA(x$.individual) || anyNA(x$.marker)) {
    stop("Individual and marker IDs cannot be missing.", call. = FALSE)
  }
  if (anyNA(x[hierarchy])) {
    stop("Hierarchy columns cannot contain missing values.", call. = FALSE)
  }
  .amova_validate_nested(x, hierarchy, ".individual")

  retained <- rep(TRUE, nrow(x))
  marker_audit <- .amova_marker_audit(x, hierarchy[length(hierarchy)])
  marker_audit$retained <- marker_audit$groups_present >= min_groups &
    marker_audit$minimum_called >= min_individuals
  if (missing == "filter") {
    marker_audit$retained <- marker_audit$retained &
      marker_audit$minimum_call_rate >= min_call_rate
  }
  keep_markers <- marker_audit$.marker[marker_audit$retained]
  x <- x[x$.marker %in% keep_markers, , drop = FALSE]
  if (!nrow(x)) stop("No markers satisfy the missing-data requirements.", call. = FALSE)

  if (missing == "complete") {
    called <- !is.na(x$.value)
    counts <- tapply(called, x$.individual, sum)
    complete_ids <- names(counts)[counts == length(unique(x$.marker))]
    x <- x[x$.individual %in% complete_ids, , drop = FALSE]
    if (!nrow(x)) stop("No complete individuals remain.", call. = FALSE)
  }

  loci <- split(x, x$.marker, drop = TRUE)
  fits <- lapply(loci, function(z) {
    z <- z[!is.na(z$.value), , drop = FALSE]
    low_counts <- table(z[[hierarchy[length(hierarchy)]]])
    z <- z[z[[hierarchy[length(hierarchy)]]] %in%
             names(low_counts)[low_counts >= min_individuals], , drop = FALSE]
    if (length(unique(z[[hierarchy[length(hierarchy)]]])) < min_groups) return(NULL)
    gr <- .amova_nested_factors(z, hierarchy)
    if (any(vapply(gr, nlevels, integer(1)) < 2L)) return(NULL)
    d2 <- .amova_distance(z$.value, distance)
    euc <- .amova_euclidean(d2, tolerance)
    correction <- 0
    if (!euc$euclidean && euclidean == "lingoes") {
      correction <- -euc$minimum_eigenvalue
      d2[row(d2) != col(d2)] <- d2[row(d2) != col(d2)] + 2 * correction
      euc <- .amova_euclidean(d2, tolerance)
    }
    fit <- .amova_fit(d2, gr)
    fit$euclidean <- euc$euclidean
    fit$minimum_eigenvalue <- euc$minimum_eigenvalue
    fit$correction <- correction
    fit$n <- nrow(z)
    fit$groups <- nlevels(gr[[length(gr)]])
    fit$marker <- as.character(z$.marker[1])
    fit$gr <- gr
    fit$raw_gr <- z[hierarchy]
    fit$d2 <- d2
    fit
  })
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("No locus was estimable after filtering.", call. = FALSE)

  component_names <- names(fits[[1]]$sigma)
  same_components <- vapply(fits, function(z) identical(names(z$sigma), component_names), logical(1))
  if (!all(same_components)) {
    stop("The represented hierarchy is inconsistent among retained loci.", call. = FALSE)
  }
  components <- do.call(rbind, lapply(fits, function(z) z$sigma))
  rownames(components) <- vapply(fits, `[[`, character(1), "marker")
  global_sigma <- colSums(components)
  phi <- .amova_phi(global_sigma)

  standardized_result <- NULL
  if (standardized) {
    max_components <- do.call(rbind, lapply(fits, function(z) {
      max_d2 <- matrix(1, nrow(z$d2), ncol(z$d2))
      diag(max_d2) <- 0
      .amova_fit_max_between(z$d2, max_d2, z$gr)$sigma
    }))
    sigma_max <- colSums(max_components)
    phi_max <- .amova_phi(sigma_max)
    standardized_result <- rep(NA_real_, length(phi))
    names(standardized_result) <- names(phi)
    standardized_result["PHI_ST"] <- phi["PHI_ST"] / phi_max["PHI_ST"]
    standardized_result[!is.finite(standardized_result)] <- NA_real_
  }

  permutation_result <- NULL
  if (permutations > 0L) {
    observed <- global_sigma
    permuted <- array(NA_real_, c(permutations, length(observed)),
                      dimnames = list(NULL, names(observed)))
    for (b in seq_len(permutations)) {
      for (j in seq_along(observed)) {
        permuted[b, j] <- sum(vapply(fits, function(z) {
          gp <- .amova_permute_component(z$raw_gr, j)
          .amova_fit(z$d2, gp)$sigma[j]
        }, numeric(1)))
      }
    }
    p <- (colSums(sweep(permuted, 2, observed, `>=`), na.rm = TRUE) + 1) /
      (permutations + 1)
    p[length(p)] <- NA_real_
    permutation_result <- list(components = permuted, p_value = p,
                               iterations = permutations, seed = seed)
  }

  locus_table <- data.frame(
    MARKERS = rownames(components),
    components,
    N = vapply(fits, `[[`, numeric(1), "n"),
    GROUPS = vapply(fits, `[[`, numeric(1), "groups"),
    EUCLIDEAN = vapply(fits, `[[`, logical(1), "euclidean"),
    MIN_EIGENVALUE = vapply(fits, `[[`, numeric(1), "minimum_eigenvalue"),
    CORRECTION = vapply(fits, `[[`, numeric(1), "correction"),
    row.names = NULL,
    check.names = FALSE
  )
  global <- data.frame(
    statistic = names(phi),
    estimate = unname(phi),
    standardized = if (is.null(standardized_result)) NA_real_ else unname(standardized_result),
    row.names = NULL
  )
  structure(
    list(
      call = match.call(), global = global, variance_components = global_sigma,
      per_locus = locus_table, marker_audit = marker_audit,
      permutation = permutation_result,
      settings = list(distance = if (is.function(distance)) "custom" else distance,
                      missing = missing, hierarchy = hierarchy,
                      value = value, euclidean = euclidean,
                      standardized = standardized)
    ),
    class = "assigner_amova"
  )
}

.amova_prepare_value <- function(data, marker, value, distance) {
  if (!is.null(value)) return(list(data = data, value = value))
  if (is.function(distance) && is.null(value)) {
    stop("A custom distance requires an explicit `value` column.", call. = FALSE)
  }
  if ("ALT_DOSAGE" %in% names(data)) {
    data$.AMOVA_DOSAGE <- as.integer(data$ALT_DOSAGE)
    return(list(data = data, value = ".AMOVA_DOSAGE"))
  }
  if (!"GT" %in% names(data)) {
    stop("No genotype column found. Supply `value`, `ALT_DOSAGE`, or `GT`.", call. = FALSE)
  }
  gt <- as.character(data$GT)
  if (all(is.na(gt) | gt %in% c("0", "1", "2"))) {
    data$.AMOVA_DOSAGE <- as.integer(gt)
    return(list(data = data, value = ".AMOVA_DOSAGE"))
  }
  if (!all(is.na(gt) | grepl("^[0-9]{6}$", gt))) {
    stop(paste0(
      "Automatic AMOVA dosage requires `ALT_DOSAGE`, numeric 0/1/2 `GT`, ",
      "or six-digit calibrated `GT`; otherwise supply `value`."
    ), call. = FALSE)
  }
  a1 <- substr(gt, 1L, 3L)
  a2 <- substr(gt, 4L, 6L)
  missing_gt <- is.na(gt) | gt == "000000"
  a1[missing_gt] <- NA_character_
  a2[missing_gt] <- NA_character_
  alleles <- data.frame(
    marker = rep(as.character(data[[marker]]), 2L),
    allele = c(a1, a2), stringsAsFactors = FALSE
  )
  alleles <- alleles[!is.na(alleles$allele), , drop = FALSE]
  allele_list <- split(alleles$allele, alleles$marker)
  n_alleles <- vapply(allele_list, function(z) length(unique(z)), integer(1))
  if (any(n_alleles > 2L)) {
    stop("Automatic dosage currently requires biallelic markers.", call. = FALSE)
  }
  ref <- vapply(allele_list, function(z) sort(unique(z))[1L], character(1))
  marker_ref <- unname(ref[as.character(data[[marker]])])
  dosage <- as.integer(a1 != marker_ref) + as.integer(a2 != marker_ref)
  dosage[missing_gt] <- NA_integer_
  data$.AMOVA_DOSAGE <- dosage
  list(data = data, value = ".AMOVA_DOSAGE")
}

#' @export
print.assigner_amova <- function(x, ...) {
  cat("Hierarchical genomic AMOVA\n")
  cat("  loci:", nrow(x$per_locus), "\n")
  cat("  missing-data strategy:", x$settings$missing, "\n")
  cat("  distance:", x$settings$distance, "\n\n")
  print(x$global, row.names = FALSE)
  invisible(x)
}

.amova_validate_nested <- function(x, hierarchy, individual) {
  id_h <- unique(x[c(individual, hierarchy)])
  if (anyDuplicated(id_h[[individual]])) {
    stop("Each individual must have one hierarchy assignment.", call. = FALSE)
  }
  if (length(hierarchy) > 1L) {
    for (j in 2:length(hierarchy)) {
      link <- unique(id_h[c(hierarchy[j - 1L], hierarchy[j])])
      if (any(table(link[[hierarchy[j]]]) > 1L)) {
        stop("`", hierarchy[j], "` is not nested within `", hierarchy[j - 1L], "`.",
             call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

.amova_marker_audit <- function(x, lowest) {
  ids <- unique(x[c(".individual", lowest)])
  totals <- data.frame(group = names(table(ids[[lowest]])),
                       total = as.numeric(table(ids[[lowest]])))
  called <- x[!is.na(x$.value), c(".marker", ".individual", lowest), drop = FALSE]
  names(called)[names(called) == lowest] <- "group"
  tab <- stats::aggregate(called$.individual,
                          called[c(".marker", "group")],
                          function(v) length(unique(v)))
  names(tab)[3] <- "called"
  tab <- merge(tab, totals, by = "group", all.x = TRUE, sort = FALSE)
  tab$rate <- tab$called / tab$total
  markers <- unique(as.character(x$.marker))
  do.call(rbind, lapply(markers, function(m) {
    z <- tab[tab$.marker == m, , drop = FALSE]
    data.frame(.marker = m, groups_present = nrow(z),
               minimum_called = if (nrow(z)) min(z$called) else 0,
               minimum_call_rate = if (nrow(z)) min(z$rate) else 0)
  }))
}

.amova_nested_factors <- function(z, hierarchy) {
  out <- vector("list", length(hierarchy))
  names(out) <- hierarchy
  for (j in seq_along(hierarchy)) {
    out[[j]] <- interaction(z[hierarchy[seq_len(j)]], drop = TRUE,
                            lex.order = TRUE, sep = "/")
  }
  as.data.frame(out, stringsAsFactors = TRUE)
}

.amova_distance <- function(value, method) {
  if (is.function(method)) {
    ans <- method(value)
    if (inherits(ans, "dist")) ans <- as.matrix(ans)
    if (!is.matrix(ans) || nrow(ans) != length(value) || ncol(ans) != length(value)) {
      stop("A custom distance function must return an n by n squared-distance matrix.", call. = FALSE)
    }
    return(ans)
  }
  if (method == "nucleotide") {
    value <- as.character(value)
    widths <- nchar(value)
    if (length(unique(widths)) != 1L) stop("Haplotype strings must have equal length.", call. = FALSE)
    chars <- strsplit(value, "", fixed = TRUE)
    ans <- outer(seq_along(chars), seq_along(chars), Vectorize(function(i, j) {
      mean(chars[[i]] != chars[[j]])
    }))
    return(ans)
  }
  if (method == "identity") return(outer(value, value, `!=`) * 1)
  if (!is.numeric(value)) stop("Dosage distances require numeric genotype values.", call. = FALSE)
  delta <- outer(value, value, `-`)
  if (method == "euclidean") delta^2 else abs(delta)
}

.amova_euclidean <- function(d2, tolerance) {
  n <- nrow(d2)
  H <- diag(n) - matrix(1 / n, n, n)
  eig <- eigen(-0.5 * H %*% d2 %*% H, symmetric = TRUE, only.values = TRUE)$values
  cutoff <- tolerance * max(1, max(abs(eig)))
  list(euclidean = min(eig) >= -cutoff, minimum_eigenvalue = min(eig))
}

.amova_fit <- function(d2, gr) {
  n <- nrow(d2)
  nlv <- ncol(gr)
  N <- lapply(gr, tabulate)
  ssd <- numeric(nlv + 2L)
  ssd[nlv + 2L] <- sum(d2) / (2 * n)
  for (i in seq_len(nlv)) {
    p <- gr[[i]]
    denom <- N[[i]][p]
    ssd[i + 1L] <- sum((d2 / (2 * denom))[outer(p, p, `==`)])
  }
  if (nlv > 1L) for (i in 2:nlv) ssd[i] <- ssd[i] - ssd[i + 1L]
  ssd[1] <- ssd[nlv + 2L] - sum(ssd[-(nlv + 2L)])
  df <- numeric(nlv + 2L)
  df[seq_len(nlv)] <- lengths(N)
  df[nlv + 1L] <- n
  for (i in (nlv + 1L):2L) df[i] <- df[i] - df[i - 1L]
  df[1] <- df[1] - 1L
  df[nlv + 2L] <- n - 1L
  if (any(df[-length(df)] <= 0)) stop("Insufficient replication in the represented hierarchy.", call. = FALSE)
  msd <- ssd / df
  coef <- .amova_ncoef(gr, N, n)
  sigma <- .amova_varcomp(msd, nlv, coef, names(gr))
  list(ssd = ssd, df = df, msd = msd, coefficient = coef, sigma = sigma)
}

.amova_ncoef <- function(gr, N, n) {
  nlv <- ncol(gr); nig <- N[[nlv]]; npop <- length(nig)
  if (nlv == 1L) return((n - sum(nig^2) / n) / (npop - 1))
  if (nlv == 2L) {
    g <- gr[[1]][match(seq_len(npop), as.integer(gr[[2]]))]
    bysum <- function(v, f) unlist(lapply(split(v, f), sum))
    A <- sum(bysum(nig^2, g) / bysum(nig, g))
    G <- nlevels(gr[[1]])
    c((n - A) / sum(tabulate(g) - 1),
      (A - sum(nig^2) / n) / (G - 1),
      (n - sum(bysum(nig, g)^2 / n)) / (G - 1))
  } else {
    out <- numeric(nlv + 1L)
    out[nlv] <- (n - sum(nig^2) / n) / (npop - 1)
    out[nlv + 1L] <- 1
    for (i in seq_len(nlv - 1L)) {
      group <- gr[[i]]
      g <- group[match(seq_len(npop), as.integer(gr[[i + 1L]]))]
      sizes <- unlist(lapply(split(nig, g), sum))
      out[i] <- (n - sum(sizes^2) / sum(sizes)) / (nlevels(group) - 1)
    }
    out
  }
}

.amova_varcomp <- function(msd, nlv, coef, level_names) {
  if (nlv == 1L) ans <- c((msd[1] - msd[2]) / coef, msd[2]) else {
    ans <- numeric(nlv + 1L)
    if (nlv == 2L) {
      ans[3] <- msd[3]
      ans[2] <- (msd[2] - ans[3]) / coef[1]
      ans[1] <- (msd[1] - ans[3] - coef[2] * ans[2]) / coef[3]
    } else {
      ans[nlv + 1L] <- msd[nlv + 1L]
      for (i in nlv:1L) {
        sel <- i:(nlv + 1L)
        ans[i] <- (msd[i] - sum(coef[sel] * ans[sel])) / coef[i]
      }
    }
  }
  names(ans) <- c(level_names, "Within")
  ans
}

.amova_phi <- function(sigma) {
  total <- sum(sigma)
  nlv <- length(sigma) - 1L
  if (!is.finite(total) || total == 0) {
    nm <- if (nlv == 1L) "PHI_ST" else if (nlv == 2L) c("PHI_ST", "PHI_CT", "PHI_SC") else paste0("PHI_", names(sigma)[-length(sigma)])
    return(stats::setNames(rep(NA_real_, length(nm)), nm))
  }
  if (nlv == 1L) return(stats::setNames(unname(sigma[1] / total), "PHI_ST"))
  if (nlv == 2L) {
    return(stats::setNames(
      c(sum(sigma[1:2]) / total,
        unname(sigma[1] / total),
        unname(sigma[2] / sum(sigma[2:3]))),
      c("PHI_ST", "PHI_CT", "PHI_SC")
    ))
  }
  cumulative <- rev(cumsum(rev(sigma[-length(sigma)]))) / total
  stats::setNames(cumulative, paste0("PHI_", names(sigma)[-length(sigma)]))
}

.amova_fit_max_between <- function(observed, maximum, gr) {
  low <- gr[[ncol(gr)]]
  use <- observed
  use[outer(low, low, `!=`)] <- maximum[outer(low, low, `!=`)]
  .amova_fit(use, gr)
}

.amova_permute_component <- function(raw_gr, component) {
  nlv <- ncol(raw_gr)
  out <- raw_gr
  if (component > nlv) return(.amova_nested_factors(out, names(out)))
  if (nlv == 1L) {
    out[[1]] <- sample(out[[1]])
    return(.amova_nested_factors(out, names(out)))
  }
  if (component == nlv) {
    parent <- interaction(out[seq_len(nlv - 1L)], drop = TRUE, lex.order = TRUE)
    out[[nlv]] <- unsplit(lapply(split(out[[nlv]], parent), sample), parent)
    return(.amova_nested_factors(out, names(out)))
  }
  child <- interaction(out[seq_len(component + 1L)], drop = TRUE, lex.order = TRUE)
  unit <- !duplicated(child)
  unit_child <- child[unit]
  labels <- out[[component]][unit]
  if (component == 1L) {
    labels <- sample(labels)
  } else {
    parent <- interaction(out[seq_len(component - 1L)], drop = TRUE, lex.order = TRUE)
    unit_parent <- parent[unit]
    labels <- unsplit(lapply(split(labels, unit_parent), sample), unit_parent)
  }
  map <- stats::setNames(as.character(labels), as.character(unit_child))
  out[[component]] <- unname(map[as.character(child)])
  .amova_nested_factors(out, names(out))
}
