# Minimize arguments list for parallel execution

Removes unused or large objects from a list of arguments. Designed to
reduce memory overhead when passing arguments to parallel workers.

## Usage

``` r
minimize_args_list(args, keep, verbose = TRUE)
```

## Arguments

- args:

  A named list of arguments (e.g., args.for.fst)

- keep:

  A character vector of argument names to keep

- verbose:

  Show what is being dropped (default: TRUE)

## Value

A slimmed-down list with only relevant arguments
