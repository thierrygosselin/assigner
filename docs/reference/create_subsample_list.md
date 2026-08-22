# Create a list of subsampled individuals

Generate a list of subsampled individuals from strata for replicated
subsampling. Optionally writes the output to file.

## Usage

``` r
create_subsample_list(
  subsample,
  strata.bk,
  iteration.subsample,
  path.folder = NULL,
  verbose = FALSE
)
```

## Arguments

- subsample:

  Integer or \`"min"\`: number of individuals to subsample per stratum,
  or \`"min"\` to use the minimum stratum size. Default: `NULL`.

- strata.bk:

  A tibble of the original strata information, with `STRATA_SEQ` as a
  grouping column.

- iteration.subsample:

  Number of subsampling replicates. Default: `1L`.

- path.folder:

  Path to the directory where subsample individuals file will be
  written. Default: `NULL`.

- verbose:

  (logical) Show messages during execution. Default: `FALSE`.

## Value

A list of tibbles, each containing subsampled individuals for one
iteration.
