# Assignment plot of genotype likelihood (Dlr).

The function generate a figure similar to Paetkau's et al. (2004) Fig 6.

## Usage

``` r
plot_dlr(
  data,
  pop.pairwise,
  dlr = NULL,
  x.dlr = NULL,
  y.dlr = NULL,
  fst = NULL,
  x.fst = NULL,
  y.fst = NULL,
  plot.dlr = TRUE,
  plot.width = NULL,
  plot.height = NULL,
  plot.dpi = NULL
)
```

## Arguments

- data:

  The assigment object from
  [`dlr`](https://thierrygosselin.github.io/assigner/reference/dlr.md)
  output.

- pop.pairwise:

  List of pairwise populations to generate the Dlr plot.

- dlr:

  (optional) Character string with Dlr value.

- x.dlr:

  (optional) Position to the x-axis of the Dlr value.

- y.dlr:

  (optional) Position to the y-axis of the Dlr value.

- fst:

  (optional) Character string with Fst value.

- x.fst:

  (optional) Position to the x-axis of the Fst value.

- y.fst:

  (optional) Position to the y-axis of the Fst value.

- plot.dlr:

  (optional) By default will write the plot in the working directory.

- plot.width:

  (optional) Width in cm of the figure.

- plot.height:

  (optional) height in cm of the figure.

- plot.dpi:

  (optional) Number of dpi for the figure (e.g 600).

## Value

A list with the assignment table and the assignment plot.

## References

Paetkau D, Slade R, Burden M, Estoup A (2004) Genetic assignment methods
for the direct, real-time estimation of migration rate: a
simulation-based exploration of accuracy and power. Molecular Ecology,
13, 55-65.

Meirmans PG, Van Tienderen PH (2004) genotype and genodive: two programs
for the analysis of genetic diversity of asexual organisms. Molecular
Ecology Notes, 4, 792-794.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
