# Machine-learning assignment engines -----------------------------------------

assign_individuals_ml <- function(data, strata = NULL,
                                  engine = c("elastic_net", "xgboost", "tabpfn"),
                                  folds = 5L, random.seed = NULL, verbose = TRUE) {
  engine <- match.arg(engine)
  if (is.null(random.seed)) random.seed <- sample.int(1e6, 1L)
  set.seed(random.seed)
  x <- genometranslator::read_genome(data, import.metadata = TRUE)
  if (!is.null(strata)) {
    s <- genometranslator::read_strata(strata)$strata
    x <- dplyr::select(x, -tidyselect::any_of("STRATA")) |>
      dplyr::left_join(dplyr::select(s, INDIVIDUALS, STRATA), by = "INDIVIDUALS")
  }
  if (!all(c("INDIVIDUALS", "MARKERS", "STRATA") %in% names(x)))
    rlang::abort("Data must contain INDIVIDUALS, MARKERS, and STRATA.")
  x <- assignment_dosage(x)
  wide <- x |> dplyr::select(INDIVIDUALS, STRATA, MARKERS, DOSAGE) |>
    tidyr::pivot_wider(names_from = MARKERS, values_from = DOSAGE)
  meta <- dplyr::select(wide, INDIVIDUALS, STRATA)
  mat <- as.matrix(dplyr::select(wide, -INDIVIDUALS, -STRATA))
  storage.mode(mat) <- "double"
  ref <- which(!is.na(meta$STRATA) & nzchar(as.character(meta$STRATA)))
  if (!length(ref)) rlang::abort("No reference individuals have STRATA.")
  classes <- sort(unique(as.character(meta$STRATA[ref])))
  y <- match(as.character(meta$STRATA), classes) - 1L
  folds <- max(2L, min(as.integer(folds), min(table(y[ref]))))
  fold.id <- integer(nrow(mat))
  for (z in split(ref, y[ref])) fold.id[z] <- sample(rep(seq_len(folds), length.out = length(z)))
  probabilities <- matrix(NA_real_, nrow(mat), length(classes), dimnames = list(NULL, classes))
  for (k in seq_len(folds)) {
    test <- ref[fold.id[ref] == k]; train <- setdiff(ref, test)
    probabilities[test, ] <- assignment_ml_fit(engine, mat[train,,drop=FALSE], y[train],
                                                 mat[test,,drop=FALSE], classes, random.seed)
  }
  unknown <- setdiff(seq_len(nrow(mat)), ref)
  if (length(unknown)) probabilities[unknown, ] <- assignment_ml_fit(
    engine, mat[ref,,drop=FALSE], y[ref], mat[unknown,,drop=FALSE], classes, random.seed)
  scores <- tibble::as_tibble(probabilities) |>
    dplyr::mutate(INDIVIDUALS = meta$INDIVIDUALS, CURRENT = meta$STRATA) |>
    tidyr::pivot_longer(tidyselect::all_of(classes), names_to = "CANDIDATE",
                        values_to = "PROBABILITY") |>
    dplyr::select(INDIVIDUALS, CURRENT, CANDIDATE, PROBABILITY)
  assignment <- scores |> dplyr::group_by(INDIVIDUALS) |>
    dplyr::slice_max(PROBABILITY, n = 1L, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(INDIVIDUALS, CURRENT, INFERRED = CANDIDATE,
                     PROBABILITY, CORRECT = dplyr::if_else(
                       is.na(CURRENT), NA, as.character(CURRENT) == INFERRED))
  missingness <- tibble::tibble(
    INDIVIDUALS = meta$INDIVIDUALS,
    N_MISSING = rowSums(is.na(mat)),
    PROP_MISSING = rowMeans(is.na(mat)))
  list(engine = engine, assignment = assignment, probabilities = scores,
       missingness = missingness, random.seed = random.seed, folds = folds)
}

assignment_dosage <- function(x) {
  if ("ALT_DOSAGE" %in% names(x)) return(dplyr::mutate(x, DOSAGE = as.integer(ALT_DOSAGE)))
  g <- as.character(x$GT)
  if (all(is.na(g) | g %in% c("0","1","2"))) return(dplyr::mutate(x, DOSAGE = as.integer(GT)))
  if (!all(is.na(g) | grepl("^[0-9]{6}$", g)))
    rlang::abort("ML engines require ALT_DOSAGE, numeric GT, or six-digit GT.")
  x$A1 <- substr(g,1,3); x$A2 <- substr(g,4,6)
  x$A1[g == "000000" | is.na(g)] <- x$A2[g == "000000" | is.na(g)] <- NA
  a <- x |> dplyr::select(MARKERS,A1,A2) |>
    tidyr::pivot_longer(c(A1,A2), values_to="ALLELE") |>
    dplyr::filter(!is.na(ALLELE)) |> dplyr::group_by(MARKERS) |>
    dplyr::summarise(REF_ALLELE=sort(unique(ALLELE))[1], N=dplyr::n_distinct(ALLELE), .groups="drop")
  if (any(a$N > 2)) rlang::abort("ML engines currently require biallelic markers.")
  x |> dplyr::left_join(a, by="MARKERS") |>
    dplyr::mutate(DOSAGE=dplyr::if_else(is.na(A1)|is.na(A2), NA_integer_,
      as.integer(A1 != REF_ALLELE)+as.integer(A2 != REF_ALLELE))) |>
    dplyr::select(-A1,-A2,-REF_ALLELE,-N)
}

assignment_ml_fit <- function(engine, train, y, test, classes, seed) {
  if (engine == "elastic_net") {
    tgbase::check_package("glmnet")
    means <- colMeans(train, na.rm=TRUE); means[!is.finite(means)] <- 0
    train <- sweep(train,2,means,FUN=function(z,m) ifelse(is.na(z),m,z))
    test <- sweep(test,2,means,FUN=function(z,m) ifelse(is.na(z),m,z))
    fit <- glmnet::cv.glmnet(train, factor(y), family="multinomial", alpha=0.5,
                             nfolds=max(3,min(5,min(table(y)))))
    p <- stats::predict(fit, test, s="lambda.min", type="response")[,,1]
    return(p[,as.character(seq_along(classes)-1L),drop=FALSE])
  }
  if (engine == "xgboost") {
    tgbase::check_package("xgboost")
    fit <- xgboost::xgb.train(
      params=list(objective="multi:softprob", num_class=length(classes),
                  eval_metric="mlogloss", max_depth=4, eta=0.05,
                  subsample=0.8, colsample_bytree=0.8, seed=seed),
      data=xgboost::xgb.DMatrix(train,label=y,missing=NA), nrounds=200, verbose=0)
    return(matrix(stats::predict(fit,xgboost::xgb.DMatrix(test,missing=NA)),
                  ncol=length(classes),byrow=TRUE))
  }
  tgbase::check_package("reticulate")
  tabpfn <- reticulate::import("tabpfn", delay_load=FALSE)
  fit <- tabpfn$TabPFNClassifier()
  fit$fit(train, y)
  as.matrix(fit$predict_proba(test))
}
