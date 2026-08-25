# Format mean and range of numeric vector as string

Produces a string like "0.123 \[0.100 - 0.150\]" for a numeric vector.

## Usage

``` r
format_mean_range(
  x,
  meanx = TRUE,
  min.max = TRUE,
  digits = 3,
  scientific = FALSE
)
```

## Arguments

- x:

  A numeric vector.

- meanx:

  Logical; whether to include the mean value in the output. Defaults to
  \`TRUE\`.

- min.max:

  Logical; whether to include the minimum and maximum values in the
  output. Defaults to \`TRUE\`.

- digits:

  Number of digits \*after the decimal point\*.

- scientific:

  Logical; whether to use scientific notation.

## Value

A character string with formatted mean, min, and max values.
