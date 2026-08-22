# Import the fst ranking from all the subsample runs inside an assignment folder.

This function will import all the fst ranking from all the subsample
runs inside an assignment folder.

## Usage

``` r
import_subsamples_fst(dir.path)
```

## Arguments

- dir.path:

  The path to the directory containing the subsample folders.

## Value

A data frame of all the Fst and ranking.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
subsamples.data <- import_subsamples_fst(
dir.path = "assignment_analysis_method_ranked_no_imputations_20160425@2321"
)
} # }
```
