# Write Temporary Data to Parquet

This function filters the input data based on subsample information,
removes missing genotypes, and optionally removes holdout individuals.
It then writes the resulting data to a Parquet file.

## Usage

``` r
write_temp_data_fst(
  subsample.list = NULL,
  data,
  holdout.samples = NULL,
  strata = NULL,
  file.date = NULL,
  path.folder = NULL,
  verbose = TRUE
)
```

## Arguments

- subsample.list:

  A data frame or tibble containing information about the subsamples. It
  should include a column \`SUBSAMPLE\` used to identify each subsample.

- data:

  A data frame or tibble containing the dataset that will be filtered
  based on the subsample and holdout samples.

- holdout.samples:

  A vector of IDs or individual identifiers to be excluded from the
  dataset (optional).

- strata:

  A data frame or tibble containing information about the strata, used
  to identify and filter holdout samples.

- file.date:

  A character string representing the date (or timestamp) to be included
  in the filename for the Parquet file.

- path.folder:

  The folder path where the resulting Parquet file will be saved.

- verbose:

  A logical indicating whether to display progress messages. Defaults to
  \`TRUE\`.

## Value

A character string representing the path to the saved Parquet file.

## Details

The function filters the \`data\` by the \`ID_SEQ\` values found in the
\`subsample.list\`. Missing genotypes (denoted by \`GT == "000000"\`)
are removed. If \`holdout.samples\` are provided, those individuals are
removed from the \`data\` based on matching IDs in the \`strata\` data
frame. The filtered data is then written to a Parquet file with the
specified filename and path.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage
write_temp_data_fst(
  subsample.list = subsample_info,
  data = my_data
)
} # }
```
