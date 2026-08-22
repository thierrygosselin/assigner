# Monte Carlo assessment of an assignment reference baseline

Repeatedly remove known reference individuals, treat them as an
artificial mixture of unknown origin, and assign them against the
remaining reference baseline. This evaluates assignment performance
without allowing test individuals to contribute to reference allele
frequencies.

## Usage

``` r
assess_assignment_mc(
  data,
  strata = NULL,
  repetitions = 50L,
  mixture.size = 100L,
  min.remaining = 5L,
  random.seed = NULL,
  ...,
  verbose = TRUE
)
```

## Arguments

- data:

  Genomic data accepted by \[assign_individuals()\].

- strata:

  Optional strata metadata accepted by
  \[genometranslator::read_strata()\]. Default: `strata = NULL`.

- repetitions:

  Number of Monte Carlo replicates. Default: `repetitions = 50`.

- mixture.size:

  Total number of reference individuals held out in every replicate.
  Default: `mixture.size = 100`.

- min.remaining:

  Minimum individuals retained in every source population. Default:
  `min.remaining = 5`.

- random.seed:

  Optional simulation seed. Default: `random.seed = NULL`.

- ...:

  Additional arguments passed to \[assign_individuals()\].

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

## Value

A list with replicate-level assignments, accuracy summaries, and the
identities held out in every replicate.

## Details

Holdout numbers are allocated approximately in proportion to the
reference population sizes, while retaining at least \`min.remaining\`
individuals in every population. Marker selection and other
data-dependent preprocessing should be performed independently within
each training baseline when they form part of the inferential workflow.

## References

Moran BM, Anderson EC (2019). Bayesian inference from the conditional
genetic stock identification model. Canadian Journal of Fisheries and
Aquatic Sciences, 76(4), 551-560.
[doi:10.1139/cjfas-2018-0016](https://doi.org/10.1139/cjfas-2018-0016) .
