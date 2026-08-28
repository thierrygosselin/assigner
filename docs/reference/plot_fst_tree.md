# Plot a population tree from pairwise FST

Build an UPGMA or neighbour-joining tree from a square pairwise FST
matrix, optionally estimate bootstrap support, and draw the result with
\`ggtree\`.

## Usage

``` r
plot_fst_tree(
  data,
  tree.method = c("upgma", "nj"),
  bootstrap = 1000L,
  scale = TRUE,
  ladderize = FALSE,
  xlim = NULL,
  plot.name = "tree",
  tree.width = 15,
  tree.height = 15
)
```

## Arguments

- data:

  A numeric, square, symmetric pairwise FST matrix with matching,
  non-empty row and column names.

- tree.method:

  Clustering method: \`"upgma"\` or \`"nj"\`.

- bootstrap:

  Number of bootstrap replicates. Use \`NULL\` or a value below three to
  skip bootstrapping.

- scale:

  Logical; display a tree scale bar.

- ladderize:

  Logical; ladderize the displayed tree.

- xlim:

  Optional numeric vector of length two passed to \[ggplot2::xlim()\].
  Use \`NULL\` to let ggplot2 choose the limits.

- plot.name:

  Output filename stem. When \`NULL\`, files are not written.

- tree.width, tree.height:

  Output dimensions in centimetres.

## Value

A \`ggtree\` plot. The underlying \`phylo\` object and bootstrap support
are available in the plot attributes \`tree\` and \`bootstrap\`.

## References

Yu G (2020). Using ggtree to visualize data on tree-like structures.
\*Current Protocols in Bioinformatics\*, 69, e96.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fst_WC84(data)
plot_fst_tree(
  data = fit$pairwise.fst.full.matrix,
  plot.name = "population_tree",
  bootstrap = 1000L
)
} # }
```
