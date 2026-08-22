#' Legacy name for assignment evaluation
#'
#' `assignment_ngs()` is the former name of [evaluate_assignment()]. New code
#' should use the more descriptive name.
#'
#' @inheritParams evaluate_assignment
#' @param ... Additional arguments passed to [evaluate_assignment()].
#' @return The value returned by [evaluate_assignment()].
#' @keywords internal
#' @export
assignment_ngs <- function(
  data,
  strata = NULL,
  pop.levels = NULL,
  assignment.analysis = c("gsi_sim", "adegenet"),
  markers.sampling = c("ranked", "random"),
  marker.number = "all",
  thl = 1,
  iteration.method = 10,
  subsample = NULL,
  iteration.subsample = 1,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1,
  ...
) {
  .Deprecated("evaluate_assignment", package = "assigner")
  evaluate_assignment(
    data = data,
    strata = strata,
    pop.levels = pop.levels,
    assignment.analysis = assignment.analysis,
    markers.sampling = markers.sampling,
    marker.number = marker.number,
    thl = thl,
    iteration.method = iteration.method,
    subsample = subsample,
    iteration.subsample = iteration.subsample,
    verbose = verbose,
    parallel.core = parallel.core,
    ...
  )
}
