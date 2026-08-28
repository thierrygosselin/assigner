#' Plot a population tree from pairwise FST
#'
#' Build an UPGMA or neighbour-joining tree from a square pairwise FST matrix,
#' optionally estimate bootstrap support, and draw the result with `ggtree`.
#'
#' @param data A numeric, square, symmetric pairwise FST matrix with matching,
#'   non-empty row and column names.
#' @param tree.method Clustering method: `"upgma"` or `"nj"`.
#' @param bootstrap Number of bootstrap replicates. Use `NULL` or a value below
#'   three to skip bootstrapping.
#' @param scale Logical; display a tree scale bar.
#' @param ladderize Logical; ladderize the displayed tree.
#' @param xlim Optional numeric vector of length two passed to
#'   [ggplot2::xlim()]. Use `NULL` to let ggplot2 choose the limits.
#' @param plot.name Output filename stem. When `NULL`, files are not written.
#' @param tree.width,tree.height Output dimensions in centimetres.
#'
#' @return A `ggtree` plot. The underlying `phylo` object and bootstrap support
#'   are available in the plot attributes `tree` and `bootstrap`.
#'
#' @examples
#' \dontrun{
#' fit <- fst_WC84(data)
#' plot_fst_tree(
#'   data = fit$pairwise.fst.full.matrix,
#'   plot.name = "population_tree",
#'   bootstrap = 1000L
#' )
#' }
#'
#' @references
#' Yu G (2020). Using ggtree to visualize data on tree-like structures.
#' *Current Protocols in Bioinformatics*, 69, e96.
#'
#' @export
plot_fst_tree <- function(
    data,
    tree.method = c("upgma", "nj"),
    bootstrap = 1000L,
    scale = TRUE,
    ladderize = FALSE,
    xlim = NULL,
    plot.name = "tree",
    tree.width = 15,
    tree.height = 15
) {
  tree.method <- match.arg(tree.method)
  if (!requireNamespace("ape", quietly = TRUE)) {
    rlang::abort("Package `ape` is required; install it with `install.packages(\"ape\")`.")
  }
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    rlang::abort(paste0(
      "Bioconductor package `ggtree` is required; install it with ",
      "`BiocManager::install(\"ggtree\")`."
    ))
  }

  data <- as.matrix(data)
  if (!is.numeric(data) || nrow(data) != ncol(data) || nrow(data) < 2L) {
    rlang::abort("`data` must be a numeric square matrix with at least two populations.")
  }
  if (is.null(rownames(data)) || is.null(colnames(data)) ||
      !identical(rownames(data), colnames(data))) {
    rlang::abort("`data` must have matching row and column names.")
  }
  if (anyNA(data) || !isTRUE(all.equal(data, t(data)))) {
    rlang::abort("`data` must be complete and symmetric.")
  }

  tree_fun <- switch(
    tree.method,
    upgma = function(x) ape::as.phylo(stats::hclust(stats::as.dist(x), method = "average")),
    nj = function(x) ape::nj(stats::as.dist(x))
  )
  tree <- tree_fun(data)
  bootstrap.value <- NULL

  if (!is.null(bootstrap) && bootstrap >= 3L) {
    bootstrap <- as.integer(bootstrap)[1L]
    bootstrap.value <- ape::boot.phylo(
      phy = tree,
      x = data,
      FUN = tree_fun,
      block = 1,
      B = bootstrap,
      trees = identical(tree.method, "upgma"),
      rooted = identical(tree.method, "upgma"),
      mc.cores = 1L
    )
    if (is.numeric(bootstrap.value)) {
      bootstrap.value <- round(100 * bootstrap.value / bootstrap)
      tree$node.label <- bootstrap.value
    } else {
      bootstrap.value <- NULL
      message("Not enough groupings to estimate bootstrap support.")
    }
  }

  tree.figure <- ggtree::ggtree(tree, ladderize = ladderize) +
    ggtree::geom_tiplab(size = 3, hjust = -0.05)
  if (isTRUE(scale)) tree.figure <- tree.figure + ggtree::geom_treescale()
  if (!is.null(xlim)) tree.figure <- tree.figure + ggplot2::xlim(xlim)

  attr(tree.figure, "tree") <- tree
  attr(tree.figure, "bootstrap") <- bootstrap.value

  if (!is.null(plot.name)) {
    ggplot2::ggsave(
      filename = paste0(plot.name, ".pdf"),
      plot = tree.figure,
      width = tree.width,
      height = tree.height,
      dpi = 600,
      units = "cm",
      useDingbats = FALSE
    )
    ggplot2::ggsave(
      filename = paste0(plot.name, ".png"),
      plot = tree.figure,
      width = tree.width,
      height = tree.height,
      dpi = 300,
      units = "cm"
    )
  }

  tree.figure
}
