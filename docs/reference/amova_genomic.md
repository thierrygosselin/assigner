# Hierarchical AMOVA for incomplete genomic data

Computes a locus-wise analysis of molecular variance (AMOVA). Missing
genotypes are handled independently at each locus and all sample sizes,
degrees of freedom, and unequal-sample-size coefficients are
recalculated from the observations that actually contribute to that
locus. This follows the locus-wise strategy used by Stacks rather than
deleting every individual with at least one missing genotype.

## Usage

``` r
amova_genomic(
  data,
  hierarchy,
  strata = NULL,
  individual = "INDIVIDUALS",
  marker = "MARKERS",
  value = NULL,
  distance = c("euclidean", "identity", "nucleotide", "manhattan"),
  missing = c("locuswise", "filter", "complete"),
  min.call.rate = 0.8,
  min.groups = 2L,
  min.individuals = 2L,
  euclidean = c("check", "lingoes", "none"),
  standardized = identical(distance, "identity"),
  permutations = 0L,
  seed = NULL,
  alpha = 0.05,
  resampling = c("none", "locus", "block"),
  bootstrap = 0L,
  confidence = 0.95,
  block = NULL,
  chromosome = NULL,
  position = NULL,
  block.size = NULL,
  population.jackknife = FALSE,
  tolerance = sqrt(.Machine$double.eps)
)
```

## Arguments

- data:

  A GDS file path or open GDS object accepted by
  [`genometranslator::read_genome()`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html).
  A long genomic data frame is also accepted for programmatic use and
  testing.

- hierarchy:

  Character vector naming nested grouping columns, ordered from highest
  to lowest level (for example `c("REGION", "STRATA")`).

- strata:

  Optional strata file or data frame accepted by
  [`genometranslator::read_strata()`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.html).
  It must contain `INDIVIDUALS` and the columns named in `hierarchy`.
  When supplied, these hierarchy columns take precedence over metadata
  stored in the GDS.

- individual, marker:

  Column names containing individual and marker IDs.

- value:

  Optional genotype-value column. With `value = NULL`, the function
  derives diploid alternate-allele dosage from `ALT_DOSAGE`, numeric
  `GT`, or six-digit calibrated `GT`. Supply a column name explicitly
  for haplotypes, nucleotide sequences, or custom distances.

- distance:

  Molecular distance. `"identity"` is the 0/1 allelic or haplotypic
  distance used for a Meirmans-style standardized statistic;
  `"nucleotide"` is the proportion of nucleotide differences between
  haplotype strings; `"euclidean"` uses squared dosage differences; and
  `"manhattan"` uses absolute dosage differences. A function may also be
  supplied; it must return a squared-distance matrix.

- missing:

  Missing-data strategy. `"locuswise"` uses all called observations at
  each locus. `"filter"` first applies marker-level call-rate thresholds
  and then works locus-wise. `"complete"` uses only individuals called
  at every retained marker. No imputation is performed.

- min.call.rate:

  Minimum called proportion required in every lowest-level group when
  `missing = "filter"`.

- min.groups:

  Minimum number of represented lowest-level groups required at a locus.

- min.individuals:

  Minimum called individuals required per represented lowest-level
  group.

- euclidean:

  How to handle a non-Euclidean squared-distance matrix. `"check"`
  records the result, `"lingoes"` applies a Lingoes correction, and
  `"none"` skips validation.

- standardized:

  Calculate a Meirmans-style standardized statistic. This is available
  for `distance = "identity"`, where the theoretical maximum is defined.
  It is deliberately not obtained by dividing distances by their
  observed sample maximum.

- permutations:

  Number of hierarchy-aware randomizations. Zero disables permutation
  tests. A single hierarchy randomization is shared by all loci in an
  iteration, so loci are not treated as independent population-level
  replicates.

- seed:

  Optional random seed used for permutations.

- alpha:

  Significance level used to flag hierarchy tests whose finite
  permutation space cannot attain the requested level.

- resampling:

  Marker uncertainty method: `"none"`, `"locus"`, or `"block"`.

- bootstrap:

  Number of marker-bootstrap replicates.

- confidence:

  Confidence level for percentile bootstrap intervals.

- block:

  Optional column containing genomic-block identifiers.

- chromosome, position:

  Columns used with `block.size` to construct physical genomic blocks
  when `block` is not supplied.

- block.size:

  Positive physical block width.

- population.jackknife:

  Calculate leave-one-lowest-level-population-out sensitivity estimates.

- tolerance:

  Numerical tolerance for Euclidean checks.

## Value

An object of class `assigner_amova` containing global statistics, locus
variance components, missing-data diagnostics, and permutation results
when requested.

## Details

`assigner` makes the finite-permutation limitation explicit and
auditable, especially for genomic datasets where thousands of loci can
give a misleading impression of high replication. A single hierarchy
randomization is applied consistently to every locus in an iteration.
Loci therefore cannot masquerade as additional population-level
replicates. The permutation design table reports the exchangeable unit,
number of unique allocations, theoretical minimum P-value, and whether
`alpha` is attainable for each component.

A permutation P-value addresses compatibility with a specified null
model; it does not describe the precision or biological magnitude of a
Phi statistic. Locus and genomic-block intervals describe
marker-sampling uncertainty; the population jackknife describes
sensitivity to particular sampled populations. Neither is a substitute
for population replication in a higher-level component, and the returned
uncertainty report labels that distinction explicitly.

## References

Excoffier L, Smouse PE, Quattro JM (1992). Analysis of molecular
variance inferred from metric distances among DNA haplotypes. Genetics
131, 479-491.

Meirmans PG (2006). Using the AMOVA framework to estimate a standardized
genetic differentiation measure. Evolution 60, 2399-2402.

Fitzpatrick BM (2009). Power and sample size for nested analysis of
molecular variance. Molecular Ecology 18, 3961-3966.

## See also

The package vignette [AMOVA for incomplete genomic
data](https://thierrygosselin.github.io/assigner/doc/amova_incomplete_genomic_data.md),
available with
[`vignette("amova_incomplete_genomic_data", package = "assigner")`](https://thierrygosselin.github.io/assigner/articles/amova_incomplete_genomic_data.md),
provides the theoretical background, implementation comparisons,
missing-data guidance, finite-permutation audit, and interpretation of
genomic AMOVA results.

## Examples

``` r
if (FALSE) { # \dontrun{
# Hierarchy stored in GDS sample metadata:
fit <- amova_genomic(
  data = "genomic_data.gds",
  hierarchy = c("REGION", "STRATA"),
  distance = "euclidean",
  missing = "locuswise",
  standardized = FALSE
)

# Hierarchy stored separately:
fit <- amova_genomic(
  data = "genomic_data.gds",
  strata = "amova_strata.tsv",
  hierarchy = c("REGION", "STRATA")
)
} # }
```
