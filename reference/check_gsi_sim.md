# Check the gsi_sim installation

\`gsi_sim\`, developed by assigner co-author Eric C. Anderson, is an
external command-line program and is not installed by the R package.
This function searches the active \`PATH\` and validates that the
executable can be started.

## Usage

``` r
check_gsi_sim(error = TRUE, verbose = TRUE)
```

## Arguments

- error:

  Logical. Stop with installation guidance when \`gsi_sim\` cannot be
  used. Default: `error = TRUE`.

- verbose:

  Logical. Display the detected executable. Default: `verbose = TRUE`.

## Value

Invisibly returns the executable path, or an empty string when it is
unavailable and `error = FALSE`.

## References

Anderson EC (2010). Assessing the power of informative subsets of loci
for population assignment: standard methods are upwardly biased.
Molecular Ecology Resources, 10, 701-710.
