# Plot the distribution of FST values

Creates a histogram of FST values from a tibble or data frame.

## Usage

``` r
plot_fst_distribution(data, binwidth = 0.01)
```

## Arguments

- data:

  A tibble or data frame containing a numeric \`FST\` column.

- binwidth:

  A numeric value specifying the histogram bin width. Default is
  \`0.01\`.

## Value

A \`ggplot2\` object showing the histogram of FST values.
