# fst_nei_subsample

Function that link all with subsampling

## Usage

``` r
fst_nei_subsample(
  x,
  input,
  pop.levels = NULL,
  pop.labels = NULL,
  strata = NULL,
  holdout.samples = NULL,
  pairwise = FALSE,
  ci = FALSE,
  iteration.ci = 100,
  quantiles.ci = c(0.025, 0.975),
  digits = 9,
  subsample = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE
)
```
