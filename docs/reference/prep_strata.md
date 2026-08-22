# Prepare stratification data

Read and join strata information to the main dataset, standardising
population levels.

## Usage

``` r
prep_strata(
  data,
  strata,
  pop.levels = NULL,
  blacklist.id = NULL,
  verbose = FALSE
)
```

## Arguments

- data:

  A tibble containing genetic data.

- strata:

  Path to the strata file or a tibble containing strata information.

- pop.levels:

  (optional) A character vector specifying the population levels to use
  for factor levels of \`STRATA\`. Default: `NULL`.

- blacklist.id:

  (optional) A vector of individual IDs to blacklist. Default: `NULL`.

- verbose:

  (logical) Show progress messages. Default: `FALSE`.

## Value

A tibble with updated \`STRATA\` column and joined strata info.
