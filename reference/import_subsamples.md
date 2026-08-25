# Import individual's assignment results of different subsample folder.

This function will import all the individual's assignment results of
different subsample folder into R.

## Usage

``` r
import_subsamples(dir.path, imputations = FALSE)
```

## Arguments

- dir.path:

  The path to the directory containing the subsample folders.

- imputations:

  (logical) Was the data imputed or not. Default: `imputations = FALSE`

## Value

A data frame of all the individual's assignment, with iterations and
subsample.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
subsamples.data <- import_subsamples(
dir.path = "assignment_analysis_method_random_imputations_rf_populations",
imputations = TRUE
)
} # }
```
