# fst_subsample

Function that link all with subsampling

## Usage

``` r
fst_subsample(
  x,
  data,
  snprelate = FALSE,
  strata = NULL,
  holdout.samples = NULL,
  pairwise = FALSE,
  ci = FALSE,
  iteration.ci = 100,
  quantiles.ci = c(0.025, 0.975),
  digits = 9,
  subsample = NULL,
  path.folder = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  forking = FALSE,
  ...
)
```
