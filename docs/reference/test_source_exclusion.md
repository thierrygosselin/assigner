# Monte Carlo exclusion test for candidate source populations

Test whether an individual's multilocus genotype is unusually improbable
under each candidate source population. This complements the relative
ranking returned by \[assign_individuals()\]: an individual can rank
first for a source while still being incompatible with every sampled
source.

## Usage

``` r
test_source_exclusion(
  result,
  data = NULL,
  simulations = 10000L,
  alpha = 0.01,
  frequency.floor = NULL,
  random.seed = NULL,
  verbose = TRUE
)
```

## Arguments

- result:

  Result returned by \[assign_individuals()\] using the native
  called-genotype likelihood engine.

- data:

  Optional genomic data used for the assignment. Supplying it lets
  simulations reproduce the exact non-missing marker set of each
  individual. Without it, the function matches only the individual's
  number of scored markers by sampling from markers available in the
  candidate source. Default: `data = NULL`.

- simulations:

  Number of simulated genotypes per candidate source. Default:
  `simulations = 10000`.

- alpha:

  Lower-tail probability used for exclusion. Default: `alpha = 0.01`.

- frequency.floor:

  Minimum allele frequency used to avoid zero probabilities. With
  \`NULL\`, use \`1 / (2 \* N + 1)\` from the reference count for each
  source-marker combination. Default: `frequency.floor = NULL`.

- random.seed:

  Optional simulation seed. Default: `random.seed = NULL`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

## Value

A list containing one row per individual and candidate source with the
observed likelihood, Monte Carlo lower-tail probability, exclusion call,
and simulation threshold; an individual summary; and simulation
settings.

## Details

For each source and marker, genotypes are generated under Hardy-Weinberg
proportions using the reference alternate-allele frequency. The same
likelihood model used for assignment scores each simulated genotype. The
Monte Carlo probability is \`(b + 1) / (B + 1)\`, where \`b\` is the
number of simulated likelihoods no greater than the observed likelihood
and \`B\` is \`simulations\`.

This is a source-exclusion diagnostic, not proof of origin or migration.
Results depend on representative reference sampling, allele-frequency
estimates, Hardy-Weinberg assumptions, genotyping error, missingness,
linkage disequilibrium, family structure, and whether the true source
was sampled. Linked markers make the product likelihood overconfident.
Repeat the test on linkage-pruned data or biologically defined marker
blocks. Candidate inversion regions should be analysed explicitly
because suppressed recombination can contribute many correlated markers.

The current simulator is for called diploid genotypes. GL/PL assignment
needs coverage-aware read or genotype-likelihood simulation and is not
silently approximated here.

## References

Cornuet JM, Piry S, Luikart G, Estoup A, Solignac M (1999). New methods
employing multilocus genotypes to select or exclude populations as
origins of individuals. Genetics, 153, 1989-2000.

Paetkau D, Slade R, Burden M, Estoup A (2004). Genetic assignment
methods for the direct, real-time estimation of migration rate: a
simulation-based exploration of accuracy and power. Molecular Ecology,
13, 55-65.
[doi:10.1046/j.1365-294X.2004.02008.x](https://doi.org/10.1046/j.1365-294X.2004.02008.x)
.

Manel S, Gaggiotti OE, Waples RS (2005). Assignment methods: matching
biological questions with appropriate techniques. Trends in Ecology &
Evolution, 20, 136-142.
[doi:10.1016/j.tree.2004.12.004](https://doi.org/10.1016/j.tree.2004.12.004)
.
