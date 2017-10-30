#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL



.onUnload <- function(libpath) {
  library.dynam.unload("grur", libpath)
}
