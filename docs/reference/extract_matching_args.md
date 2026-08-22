# Extract arguments matching a function's formals

Extracts named arguments from an environment (e.g., the calling
environment) that match the formal arguments of a target function.

Useful when forwarding arguments from a parent function to a subfunction
without manually repeating each parameter.

## Usage

``` r
extract_matching_args(from.env, to.fn, .evaluate = TRUE, .exclude = NULL)
```

## Arguments

- from.env:

  The environment from which to extract arguments (typically
  \`environment()\` or \`rlang::caller_env()\`).

- to.fn:

  The target function whose formal arguments will be matched.

- .evaluate:

  Logical, whether to force evaluation (default: TRUE).

- .exclude:

  Optional character vector of argument names to exclude.

## Value

A named list of matched and evaluated arguments (excluding
\`.exclude\`).

## Examples

``` r
if (FALSE) { # \dontrun{
args.for.fst <- extract_matching_args(from.env = environment(), to.fn = assigner::compute_fst)
result <- rlang::exec(assigner::compute_fst, x = my_data, !!!args.for.fst)
} # }
```
