# assigner parallel function

Updating assigner to use future

## Usage

``` r
assigner_future(
  .x,
  .f,
  flat.future = c("int", "chr", "dfr", "dfc", "walk", "drop"),
  split.vec = FALSE,
  split.with = NULL,
  split.chunks = 4L,
  parallel.core = parallel::detectCores() - 1,
  forking = FALSE,
  ...
)
```
