# Internal Matrix Utilities for Pairwise Distance Data

These internal utility functions convert long-format pairwise population
data into symmetric matrices, typically used for visualising pairwise
FST or genetic distance estimates.

## Usage

``` r
make_upper_matrix(data, value.col, fill.value = "")

make_full_symmetric_matrix(upper.matrix, diagonal.value = "0")
```

## Arguments

- data:

  A data frame or tibble with \`POP1\`, \`POP2\`, and a value column.

- value.col:

  The unquoted column to use as matrix values.

- fill.value:

  The value to fill missing cells (default: \`""\`).

- upper.matrix:

  A square matrix (upper triangle filled).

- diagonal.value:

  The value to assign to the diagonal (default: \`"0"\`).

## Value

A matrix. The first function returns an upper-triangular matrix, and the
second a full symmetric matrix with diagonal.

## Details

\`make_upper_matrix()\` converts long-form pairwise data into a wide
upper-triangular matrix. \`make_full_symmetric_matrix()\` takes that
upper matrix and completes the symmetric form.
