# Genotype likelihood ratio distance (Dlr)

The function computes Paetkau's et al. (1997) genotype likelihood ratio
distance (Dlr).

## Usage

``` r
dlr(
  data,
  strata = NULL,
  plots = FALSE,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```

## Arguments

- data:

  A result returned by \[assign_individuals()\], its complete
  \`likelihoods\` table, or a GenoDive assignment output file.

- strata:

  A tab delimited file with 2 columns with header: `INDIVIDUALS` and
  `STRATA`. The `STRATA` column is used here as the populations id of
  your sample.

- plots:

  (optional) Generate Dlr plots for all the pairwise populations in the
  dataset. The plots are `ggplot2` objects that can be modified with the
  proper `ggplot2` syntax. Default: `plots = FALSE`.

- filename:

  (optional) Name of the file prefix for the matrix and the table
  written in the working directory.

- parallel.core:

  (optional) The number of core for parallel computation. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

## Value

A list with 5 objects:

1.  the assignment results (\$assignment),

2.  the dlr pairwise table (\$dlr.table),

3.  the lower diagonal dlr distance matrix (\$dlr.dist),

4.  a data.frame with the dlr distance mirrored (\$dlr.matrix),

5.  the list of dlr plots (\$dlr.plots)

## References

Paetkau D, Slade R, Burden M, Estoup A (2004) Genetic assignment methods
for the direct, real-time estimation of migration rate: a
simulation-based exploration of accuracy and power. Molecular Ecology,
13, 55-65.

Paetkau D, Waits LP, Clarkson PL, Craighead L, Strobeck C (1997) An
empirical evaluation of genetic distance statistics using microsatellite
data from bear (Ursidae) populations. Genetics, 147, 1943-1957.

Meirmans PG, Van Tienderen PH (2004) genotype and genodive: two programs
for the analysis of genetic diversity of asexual organisms. Molecular
Ecology Notes, 4, 792-794.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
dlr <- assigner::dlr(
data = "assignment.gdv", strata = "my.strata.tsv", plots = TRUE)

# to get the plots list:
plot.list <- dlr$dlr.plots
# access and isolate in different object a plot with $
} # }
```
