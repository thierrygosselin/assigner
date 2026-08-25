# AMOVA for incomplete genomic data

## Purpose

Analysis of molecular variance (AMOVA) partitions molecular variation
among nested biological levels using a matrix of squared distances
([Excoffier et al. 1992](#ref-Excoffier1992)). It is a general framework
rather than a single genetic distance or a single statistic.

Several mature programs and R packages already implement AMOVA. The goal
of
[`assigner::amova_genomic()`](https://thierrygosselin.github.io/assigner/reference/amova_genomic.md)
is not to replace or discredit them. Its narrower goal is to make the
observations contributing to AMOVA explicit when a genomic matrix is
incomplete: which individuals and populations contributed at each locus,
which loci were retained, whether a distance was Euclidean, and how the
global variance components were assembled.

This vignette explains the theory, compares existing implementations,
and documents the choices made in `assigner`.

## AMOVA in brief

### From observations to squared distances

For $`n`$ molecular observations, let $`\Delta=(\delta_{ij}^2)`$ be the
matrix of squared distances. The total sum of squared deviations can be
written as

``` math
SS_T = \frac{1}{2n}\sum_{i=1}^{n}\sum_{j=1}^{n}\delta_{ij}^2.
```

Equivalent sums are calculated within populations and, when present,
within successively higher groups. Differences between these sums
provide the sums of squares assigned to each hierarchical level.
Dividing by the corresponding degrees of freedom gives mean squares.
Method-of-moments equations then estimate variance components while
accounting for unequal sample sizes ([Excoffier et al.
1992](#ref-Excoffier1992)).

For a hierarchy `REGION/STRATA`, the variance components are
conventionally written as

- $`\sigma_a^2`$: among regions;
- $`\sigma_b^2`$: among populations within regions;
- $`\sigma_c^2`$: within populations.

The associated statistics are

``` math
\Phi_{ST}=\frac{\sigma_a^2+\sigma_b^2}
                 {\sigma_a^2+\sigma_b^2+\sigma_c^2},
```

``` math
\Phi_{CT}=\frac{\sigma_a^2}
                 {\sigma_a^2+\sigma_b^2+\sigma_c^2},
\qquad
\Phi_{SC}=\frac{\sigma_b^2}{\sigma_b^2+\sigma_c^2}.
```

The notation describes correlations of molecular differences at
specified hierarchical levels. It should not be confused automatically
with every estimator called $`F_{ST}`$. Weir and Cockerham’s estimator,
for example, is derived from allelic correlations and sampling theory
rather than from an arbitrary molecular-distance matrix ([Weir and
Cockerham 1984](#ref-WeirCockerham1984)).

### A nested-looking table is not necessarily an AMOVA hierarchy

AMOVA requires more than columns that can be arranged as
`REGION/STRATA`. The lower units must be biologically meaningful
replicates for the higher-level inference. Capture locality, life stage,
migratory route, inferred cluster, and reproductive population are not
automatically interchangeable.

Two recent highly migratory fish studies illustrate the distinction.
Mikles et al. combined whole-genome data with independently observed
movements of Atlantic bluefin tuna from Gulf of Mexico and Mediterranean
spawning grounds ([Mikles et al. 2026](#ref-Mikles2026)). Adults and
larvae form nested-looking categories, but they are life stages rather
than replicated populations within each spawning region. With only two
primary spawning populations, higher-level AMOVA has essentially no
population replication. Their movement, ordination, differentiation, and
demographic analyses addressed different and appropriate questions.

Chevrier et al. sampled swordfish at many Indian Ocean localities and
recovered a north-south signal associated strongly with a possible
chromosomal inversion ([Chevrier et al. 2024](#ref-Chevrier2024)).
Localities could form the lowest AMOVA level, but capture sites may mix
reproductive origins, and testing a grouping discovered from the same
SNPs is post hoc rather than independent confirmation ([Meirmans
2015](#ref-Meirmans2015)). A useful AMOVA would be descriptive or use
regions specified before examining the genotypes, and would compare
linkage-pruned neutral loci, the inversion, and the genome with that
region removed.

Before fitting a hierarchy, ask whether the grouping was specified
independently, whether lower units are biological replicates, whether
enough units occur in every group, and whether one linked region could
dominate the component.
[`amova_genomic()`](https://thierrygosselin.github.io/assigner/reference/amova_genomic.md)
can audit permutation support and linkage-aware marker uncertainty, but
it cannot decide whether a metadata column represents a valid biological
population. That decision belongs to the sampling design.

### The distance is part of the biological model

The original AMOVA paper emphasizes that alternative distance matrices
encode alternative assumptions about molecular evolution ([Excoffier et
al. 1992](#ref-Excoffier1992)). Examples include:

- identity distance: zero for identical states and one otherwise;
- nucleotide distance: the number or proportion of nucleotide
  substitutions;
- allele-sharing distances for codominant genotypes;
- Euclidean or Manhattan distances between allele dosages;
- user-defined distances appropriate to a marker system.

Changing the distance can change the estimand. Consequently, a function
should record the distance used rather than label all analyses simply as
“AMOVA”.

Classical AMOVA uses a squared Euclidean distance matrix. A
non-Euclidean dissimilarity can generate negative eigenvalues and
variance components that do not have the usual Euclidean decomposition.
`assigner` therefore records an eigenvalue diagnostic and can apply a
Lingoes correction. A correction changes the geometry and must be
reported; it is not merely a numerical detail.

### Negative variance components

Sampling estimates of variance components can be negative. A negative
estimate does not necessarily indicate a programming error, and
replacing it silently with zero changes the estimator.
[`amova_genomic()`](https://thierrygosselin.github.io/assigner/reference/amova_genomic.md)
retains negative components. Biological interpretation should consider
sampling uncertainty, the distance, and the sampling design.

## Standardized differentiation is a second calculation

$`\Phi_{ST}`$ may have a low attainable maximum when within-population
diversity is high. Meirmans proposed standardizing the statistic by its
maximum possible value conditional on the observed within-population
component ([Meirmans 2006](#ref-Meirmans2006)):

``` math
\Phi'_{ST}=\frac{\Phi_{ST}}{\Phi_{ST(max)}}.
```

The maximum must be defined by the distance model. It is not, in
general, the largest distance observed in a sample. Dividing every locus
by its observed maximum makes the result depend on sample composition
and missingness.

Following the logic used by GenoDive and Stacks, `assigner` currently
computes the standardized statistic only with identity distance. The
observed within-population distances are retained, whereas
between-population distances are replaced by their theoretical maximum.
Other distances can be used for $`\Phi_{ST}`$, but they do not receive
an automatic $`\Phi'_{ST}`$ unless a defensible theoretical maximum is
supplied.

This is why $`\Phi_{ST}`$ and $`\Phi'_{ST}`$ need not use the same
distance matrix. Reporting both can be useful because ordinary and
standardized differentiation answer related but distinct questions
([Meirmans and Hedrick 2011](#ref-MeirmansHedrick2011)).

## Why incomplete genomic data need special treatment

With hundreds or thousands of RADseq loci, requiring every individual to
be called at every locus can discard most of the dataset. At the other
extreme, replacing missing values with zero treats absence of
information as a measured allele count. A single pairwise distance
calculated from a different set of loci for every pair can also be
difficult to audit and may not be Euclidean.

Missingness is additionally capable of being informative. Coverage,
restriction site polymorphism, mapping quality, population ancestry, and
laboratory batch can all affect whether a genotype is observed. No
distance formula can make such missingness ignorable by itself.

`assigner` therefore treats missing-data policy as part of the analysis:

1.  **Locus-wise** (default): remove uncalled observations only for the
    current locus, then recompute the represented hierarchy, sample
    sizes, degrees of freedom, and unequal-sample-size coefficients.
2.  **Filter**: apply explicit call-rate and representation thresholds,
    then use the locus-wise calculation.
3.  **Complete**: retain only individuals called at all retained loci.
    This is provided as a diagnostic or for genuinely dense datasets.

No imputation is performed. Imputation belongs in `grur`, where multiple
completed datasets and their uncertainty can be generated explicitly. A
later `assigner` interface can combine AMOVA estimates across those
datasets rather than hiding imputation inside the variance calculation.

## Existing implementations

The following comparison describes differences in scope and interface.
It is not a ranking of scientific quality.

### `ade4`

[`ade4::amova()`](https://adeverse.github.io/ade4/reference/amova.html)
is a direct and general implementation based on haplotype-count tables,
an optional distance object, and optional nested structures. It verifies
that a supplied distance is Euclidean and
[`ade4::randtest()`](https://adeverse.github.io/ade4/reference/randtest.html)
provides randomization tests. `ade4` is a foundational
multivariate-analysis package developed by an expert group that includes
Stéphane Dray, Anne-Béatrice Dufour, Jean Thioulouse, Daniel Chessel,
Thibaut Jombart, and contributors including Pierre Legendre ([Dray and
Dufour 2007](#ref-DrayDufour2007)).

Its AMOVA interface expects the user to provide valid sample counts and
a complete distance object. That is a sensible separation of
responsibilities, but it means that locus-dependent missingness must be
resolved before calling the function. It does not itself maintain an
audit of different contributing individuals and populations for
thousands of genomic loci.

### `pegas`

`pegas::amova()` uses a convenient formula interface such as
`distance ~ region/population`. It accepts a distance matrix, calculates
hierarchical variance components, and implements level-specific
permutations ([Paradis 2010](#ref-Paradis2010)). Its source uses the
exact Excoffier coefficients for one- and two-level hierarchies and
documented approximate formulae for deeper hierarchies.

As with `ade4`, the distance matrix is the analytical input. Missing
values in that matrix are warned about, but constructing a defensible
distance matrix from sparse genomic observations remains the user’s
responsibility.

### `poppr`

`poppr::poppr.amova()` provides an important bridge from `adegenet`
genetic objects to either the `ade4` or `pegas` AMOVA engines. It
supports clonal data, haplotype splitting, mixed reproductive systems,
missing-data preprocessing, and corrections for non-Euclidean distances
([Kamvar et al. 2014](#ref-Kamvar2014), [2015](#ref-Kamvar2015)).

For `genind` objects, missing data are handled before AMOVA through
options such as removing loci or genotypes or replacing missing allele
frequencies. For `genlight`, its bitwise distance can scale pairwise
comparisons for missing observations. These are useful general policies,
but they differ from rebuilding the AMOVA sampling coefficients
separately from the observations present at every locus. `poppr` also
has a broader and different purpose: population genetics with particular
strength for clonal and partially clonal organisms.

### GenoDive

GenoDive implements hierarchical AMOVA, ordinary and standardized
differentiation, numerous distances, filtering, and simple procedures
for filling missing genotypes ([Meirmans 2020](#ref-Meirmans2020)). The
published description lists three procedures: randomly drawing alleles
from overall allele frequencies, randomly drawing alleles from
population allele frequencies, and restoring polyploid dosage from
incomplete genotypes and population allele frequencies. These are more
specific than replacing every missing value by a column mean, but they
are single-completion, frequency-based approaches. They should not be
described as a modern multiple-imputation framework designed to
propagate uncertainty from the high and structured missingness often
observed in RADseq datasets. GenoDive’s standardized AMOVA is closely
connected to Meirmans’ derivation of $`\Phi'_{ST}`$([Meirmans
2006](#ref-Meirmans2006)), making it an important reference
implementation and a useful independent comparison for `assigner`.

GenoDive is a standalone application rather than an R package with an
exposed R implementation. The published papers document its
capabilities, but not every internal computational decision can be
verified from the material distributed with `assigner`. Comparisons
should therefore use controlled input datasets and document software
versions and settings rather than infer undocumented details.

### Stacks

The `populations` program in Stacks was designed around RADseq data
([Rochette et al. 2019](#ref-RochetteCatchen2019)). Its source
illustrates a useful genomic strategy: observed haplotypes are tabulated
independently at each locus, populations without valid observations are
excluded at that locus, and the AMOVA sample-size coefficients are
recalculated from the represented haplotypes. It supports minimum sample
proportions per population and minimum numbers of represented
populations.

Stacks calculates haplotype $`\Phi_{ST}`$ using nucleotide-substitution
distances. For standardized differentiation, it separately uses identity
distances and constructs a maximum between-population distance. These
choices strongly influenced the `assigner` design.

Stacks is nevertheless a complete RADseq processing pipeline with its
own data model and output conventions. It is not intended as a general R
AMOVA engine for arbitrary tidy genomic inputs or user-defined
distances.

## What is different in `assigner`?

[`assigner::amova_genomic()`](https://thierrygosselin.github.io/assigner/reference/amova_genomic.md)
combines established AMOVA calculations with an explicit
incomplete-genomic-data contract:

- long genomic data are analyzed locus by locus;
- the observations and hierarchical units contributing to every locus
  are recalculated rather than inferred from the full input strata;
- global statistics are ratios of **summed locus variance components**,
  not arithmetic means of locus $`\Phi`$-statistics;
- marker call rates, represented groups, sample sizes, retained loci,
  eigenvalues, and geometric corrections are returned for audit;
- distance methods have explicit names and allow a custom
  squared-distance function;
- ordinary $`\Phi`$-statistics and standardized $`\Phi'_{ST}`$ are kept
  as separate calculations;
- permutation occurs at the exchangeable unit appropriate to the tested
  hierarchy level;
- imputation is intentionally outside the function.

The contribution is therefore not a new definition of AMOVA. It is a
reproducible way to connect AMOVA to incomplete next-generation genotype
data while exposing decisions that are otherwise commonly made during
preprocessing.

## Running `amova_genomic()`

The public input follows the rest of `assigner`: `data` may be a GDS
file path or an open GDS object handled by
[`genometranslator::read_genome()`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html).
Diploid alternate-allele dosage is detected automatically from standard
GDS genotype fields. A long data frame remains supported for
programmatic use and examples. Hierarchy columns are ordered from
highest to lowest and may come from GDS metadata or the separate
`strata` argument.

``` r

set.seed(2026)

example_amova <- expand.grid(
  INDIVIDUALS = sprintf("ind_%02d", 1:24),
  MARKERS = paste0("loc_", 1:12),
  stringsAsFactors = FALSE
)

individual_strata <- data.frame(
  INDIVIDUALS = sprintf("ind_%02d", 1:24),
  REGION = rep(c("north", "south"), each = 12),
  POP_ID = rep(c("n1", "n2", "s1", "s2"), each = 6)
)

example_amova <- merge(example_amova, individual_strata, by = "INDIVIDUALS")
population_effect <- c(n1 = 0.1, n2 = 0.6, s1 = 1.2, s2 = 1.7)
example_amova$GT_BIN <- rbinom(
  nrow(example_amova), 2,
  plogis(-1 + population_effect[example_amova$POP_ID])
)
example_amova$GT_BIN[sample(nrow(example_amova), 30)] <- NA

head(example_amova)
#>   INDIVIDUALS MARKERS REGION POP_ID GT_BIN
#> 1      ind_01   loc_1  north     n1      1
#> 2      ind_01   loc_6  north     n1     NA
#> 3      ind_01  loc_11  north     n1      0
#> 4      ind_01   loc_5  north     n1      0
#> 5      ind_01  loc_10  north     n1      1
#> 6      ind_01   loc_8  north     n1      0
```

### Locus-wise AMOVA

``` r

fit <- amova_genomic(
  data = example_amova,
  hierarchy = c("REGION", "POP_ID"),
  value = "GT_BIN",
  distance = "euclidean",
  missing = "locuswise",
  min.groups = 2,
  min.individuals = 2,
  standardized = FALSE
)

fit
#> Hierarchical genomic AMOVA
#>   loci: 12 
#>   missing-data strategy: locuswise 
#>   distance: euclidean 
#> 
#>  statistic   estimate standardized
#>     PHI_ST 0.32412536           NA
#>     PHI_CT 0.30824170           NA
#>     PHI_SC 0.02296129           NA
head(fit$per_locus)
#>   MARKERS     REGION      POP_ID    Within  N GROUPS EUCLIDEAN MIN_EIGENVALUE
#> 1   loc_1 -0.1828342  0.26878189 0.4739583 20      4      TRUE  -1.246316e-15
#> 2  loc_10  0.4545455 -0.05000000 0.3000000 23      4      TRUE  -2.039205e-15
#> 3  loc_11  0.3264038 -0.05487714 0.4188889 19      4      TRUE  -2.809005e-15
#> 4  loc_12  0.5063272  0.03518519 0.2981481 22      4      TRUE  -2.556524e-15
#> 5   loc_2  0.2358076 -0.06604938 0.6814815 22      4      TRUE  -1.323316e-15
#> 6   loc_3  0.1589945 -0.06736720 0.3598039 21      4      TRUE  -1.693783e-15
#>   CORRECTION
#> 1          0
#> 2          0
#> 3          0
#> 4          0
#> 5          0
#> 6          0
head(fit$marker_audit)
#>   .marker groups_present minimum_called minimum_call_rate retained
#> 1   loc_1              4              4         0.6666667     TRUE
#> 2   loc_6              4              5         0.8333333     TRUE
#> 3  loc_11              4              4         0.6666667     TRUE
#> 4   loc_5              4              5         0.8333333     TRUE
#> 5  loc_10              4              5         0.8333333     TRUE
#> 6   loc_8              4              5         0.8333333     TRUE
```

### Filtering before locus-wise calculation

``` r

fit_filtered <- amova_genomic(
  data = example_amova,
  hierarchy = c("REGION", "POP_ID"),
  value = "GT_BIN",
  distance = "euclidean",
  missing = "filter",
  min.call.rate = 0.75,
  min.groups = 4,
  min.individuals = 2,
  standardized = FALSE
)

fit_filtered$global
#>   statistic     estimate standardized
#> 1    PHI_ST  0.340159897           NA
#> 2    PHI_CT  0.342975060           NA
#> 3    PHI_SC -0.004284713           NA
```

### Standardized differentiation

For a Meirmans-style standardized statistic, use an identity state and
identity distance. The `standardized` column is populated for `PHI_ST`;
it is not silently generalized to hierarchy-specific standardized
statistics whose maxima have not been implemented.

``` r

fit_identity <- amova_genomic(
  data = haplotype_data,
  hierarchy = c("REGION", "POP_ID"),
  value = "HAPLOTYPE",
  distance = "identity",
  missing = "locuswise",
  standardized = TRUE
)
```

### Hierarchical permutations

``` r

fit_permuted <- amova_genomic(
  data = example_amova,
  hierarchy = c("REGION", "POP_ID"),
  value = "GT_BIN",
  distance = "euclidean",
  missing = "locuswise",
  standardized = FALSE,
  permutations = 999,
  seed = 2026
)

fit_permuted$permutation$p_value
fit_permuted$permutation$design
```

#### The finite permutation space is part of the result

Fitzpatrick’s key observation is combinatorial, not a competing
estimator of AMOVA ([Fitzpatrick 2009](#ref-Fitzpatrick2009)). To test a
group-level component, populations are the exchangeable replicates. More
individuals or loci can improve estimation of within-population genetic
variation, but they do not create more independent ways to allocate the
sampled populations among groups.

For group sizes $`k_1,\ldots,k_g`$,
[`amova_genomic()`](https://thierrygosselin.github.io/assigner/reference/amova_genomic.md)
calculates the number of unique allocations from the multinomial
coefficient and removes duplicate partitions caused by groups of equal
size. For two groups containing two populations each, there are only
three unique partitions. Consequently the smallest exact probability is
$`1/3`$, and a group-level test cannot attain $`\alpha=0.05`$,
regardless of how many SNPs or individuals were sampled.

The `design` table reports, for every tested component:

- the unit that is exchangeable under its null hypothesis;
- the number of those units in the sampling design;
- the number of unique allocations;
- the theoretical minimum attainable $`P`$-value; and
- whether the requested `alpha` can be attained.

The function warns before interpretation when a test cannot attain
`alpha`. This is intentionally different from merely reporting the Monte
Carlo resolution $`1/(B+1)`$ for `B` sampled permutations ([Phipson and
Smyth 2010](#ref-PhipsonSmyth2010)). The returned `p_value` is not
allowed below the exact finite-space bound; `monte_carlo_p` is retained
for auditing the random approximation.

#### Genomic loci do not multiply the high-level replication

With incomplete genomic data, each locus can contain a different subset
of individuals. Nevertheless, one hierarchy randomization is drawn per
component and iteration and then applied consistently to every retained
locus. Permuting the hierarchy independently at every locus would
incorrectly turn loci into additional population-level replication and
could make a weak sampling design look much more informative than it is.

Fitzpatrick’s conclusion has been carried into later sampling-design
guidance ([Werth 2011](#ref-Werth2011)): for a hierarchical group test,
sampling additional populations per group is more useful than adding
individuals to a few populations. The limitation is also a standard
property of restricted randomization tests. It is therefore appropriate
to implement it as an audit and validity safeguard, not to present it as
a new form of AMOVA.

The AMOVA implementations reviewed for this vignette provide
user-selected randomization counts, but their ordinary result objects do
not generally expose the design-specific number of unique high-level
allocations or warn that the chosen significance level is unattainable.
This does **not** mean their AMOVA estimators are wrong, nor that their
authors rejected Fitzpatrick’s argument. The distinction in `assigner`
is that the limitation is calculated and reported alongside the
component test, which is particularly important when thousands of loci
can give a misleading impression of replication.

#### P-values are not effect-size uncertainty

Permutation tests and confidence intervals answer different questions. A
permutation $`P`$-value asks how compatible the observed component is
with a particular null randomization. It does not show how precisely a
$`\Phi`$-statistic was estimated, and a small value does not imply
biologically important differentiation. Conversely, a scientifically
meaningful estimate can remain uncertain or have an unattainable
permutation threshold when few populations were sampled.

Confidence intervals for $`F`$- and $`\Phi`$-statistics are consequently
useful: they display uncertainty around the estimated magnitude and
allow readers to judge whether zero and biologically negligible effects
remain plausible. They should be reported alongside the estimate rather
than treated merely as another route to a significance label.

For hierarchical genomic AMOVA, however, the bootstrap design is part of
the estimand:

- resampling individuals within populations represents uncertainty from
  the sampled individuals, but cannot create population-level
  replication;
- resampling populations within their parent groups can address a
  higher-level component when populations are the intended sampling
  units;
- resampling single SNPs as if independent is anticonservative under
  linkage; a locus bootstrap is appropriate only for effectively
  independent loci, while a genomic block bootstrap requires genomic
  positions and defensible blocks;
- uncertainty caused by imputation or non-random missingness is a
  separate layer and should eventually be propagated from `grur`
  completed datasets.

[`amova_genomic()`](https://thierrygosselin.github.io/assigner/reference/amova_genomic.md)
therefore offers explicitly labelled locus and genomic-block percentile
intervals. Locus resampling assumes effectively independent loci; block
resampling is preferred when linkage can be represented by an existing
block column or by chromosome, position, and a defensible physical block
width.

``` r

fit_uncertainty <- amova_genomic(
  data = gds,
  hierarchy = c("REGION", "STRATA"),
  resampling = "block",
  chromosome = "CHROM", position = "POSITION", block.size = 1e6,
  bootstrap = 999, confidence = 0.95,
  population.jackknife = TRUE, seed = 2026
)
fit_uncertainty$uncertainty$report
fit_uncertainty$uncertainty$marker$intervals
fit_uncertainty$uncertainty$population.jackknife
```

The population jackknife refits AMOVA after omitting each lowest-level
population. It measures leverage, not population-sampling uncertainty.
Marker intervals likewise do not represent population replication.

## Validation strategy

Before treating this implementation as a reference estimator, it should
be validated in layers:

1.  **Complete balanced data:** compare sums of squares, mean squares,
    variance components, and $`\Phi`$-statistics with `ade4`, `pegas`,
    GenoDive, and Stacks where their input models overlap.
2.  **Complete unbalanced data:** verify unequal-sample-size
    coefficients against worked calculations and independent software.
3.  **Incomplete data under MCAR:** introduce known call rates and
    compare bias, variance, coverage, and locus retention.
4.  **Population- and genotype-associated missingness:** demonstrate
    where a locus-wise analysis remains biased because the missingness
    mechanism is not ignorable.
5.  **Distance sensitivity:** compare identity, nucleotide,
    allele-sharing, and dosage distances without assuming they estimate
    the same biological quantity.
6.  **Permutation calibration:** verify type-I error and power for each
    hierarchy level, with special attention to the number of
    exchangeable populations and groups.
7.  **Performance:** profile realistic RADseq datasets before moving
    numerical kernels to C++. Statistical orchestration should remain
    visible in R.

## Recommended reading

The following sequence gives a useful path through the subject:

1.  Excoffier, Smouse, and Quattro for the original distance-based AMOVA
    framework ([Excoffier et al. 1992](#ref-Excoffier1992)).
2.  Weir and Cockerham for sampling-based (F)-statistics and why
    estimators with similar labels need not be identical ([Weir and
    Cockerham 1984](#ref-WeirCockerham1984)).
3.  Meirmans for standardized AMOVA differentiation ([Meirmans
    2006](#ref-Meirmans2006)).
4.  Meirmans and Hedrick for the interpretation of $`F_{ST}`$,
    standardized (F’\_{ST}), and related measures ([Meirmans and Hedrick
    2011](#ref-MeirmansHedrick2011)).
5.  Fitzpatrick for the limits of hierarchical permutation power
    ([Fitzpatrick 2009](#ref-Fitzpatrick2009)).
6.  The `ade4`, `pegas`, and `poppr` papers for their respective R
    ecosystems ([Dray and Dufour 2007](#ref-DrayDufour2007); [Paradis
    2010](#ref-Paradis2010); [Kamvar et al. 2014](#ref-Kamvar2014),
    [2015](#ref-Kamvar2015)).
7.  The GenoDive and Stacks 2 papers for practical population-genomic
    workflows ([Meirmans 2020](#ref-Meirmans2020); [Rochette et al.
    2019](#ref-RochetteCatchen2019)).

## Reproducibility checklist

An AMOVA report should state at least:

- software and version;
- hierarchy and sampling unit at every level;
- molecular distance and whether it was squared;
- missing-data and filtering policy;
- number of retained loci and samples per level;
- Euclidean diagnostic and any correction;
- whether global values were obtained from summed components or averaged
  locus statistics;
- definition of the maximum used for a standardized statistic;
- permutation unit, restrictions, seed, and number of permutations;
- whether negative variance components were retained or truncated.

Chevrier, Thomas, Dominique A. Cowart, Anne-Elise Nieblas, et al. 2024.
“Population Structure of the Swordfish, Xiphias Gladius, Across the
Indian Ocean Using Next-Generation Sequencing.” *ICES Journal of Marine
Science* 82 (5): fsae179. <https://doi.org/10.1093/icesjms/fsae179>.

Dray, Stéphane, and Anne-Béatrice Dufour. 2007. “The Ade4 Package:
Implementing the Duality Diagram for Ecologists.” *Journal of
Statistical Software* 22 (4): 1–20.
<https://doi.org/10.18637/jss.v022.i04>.

Excoffier, Laurent, Peter E. Smouse, and Joseph M. Quattro. 1992.
“Analysis of Molecular Variance Inferred from Metric Distances Among DNA
Haplotypes: Application to Human Mitochondrial DNA Restriction Data.”
*Genetics* 131: 479–91. <https://doi.org/10.1093/genetics/131.2.479>.

Fitzpatrick, Benjamin M. 2009. “Power and Sample Size for Nested
Analysis of Molecular Variance.” *Molecular Ecology* 18 (19): 3961–66.
<https://doi.org/10.1111/j.1365-294X.2009.04314.x>.

Kamvar, Zhian N., Jonah C. Brooks, and Niklaus J. Grünwald. 2015. “Novel
r Tools for Analysis of Genome-Wide Population Genetic Data with
Emphasis on Clonality.” *Frontiers in Genetics* 6: 208.
<https://doi.org/10.3389/fgene.2015.00208>.

Kamvar, Zhian N., Javier F. Tabima, and Niklaus J. Grünwald. 2014.
“Poppr: An r Package for Genetic Analysis of Populations with Clonal,
Partially Clonal, and/or Sexual Reproduction.” *PeerJ* 2: e281.
<https://doi.org/10.7717/peerj.281>.

Meirmans, Patrick G. 2006. “Using the AMOVA Framework to Estimate a
Standardized Genetic Differentiation Measure.” *Evolution* 60 (11):
2399–402. <https://doi.org/10.1111/j.0014-3820.2006.tb01874.x>.

Meirmans, Patrick G. 2015. “Seven Common Mistakes in Population Genetics
and How to Avoid Them.” *Molecular Ecology* 24 (13): 3223–31.
<https://doi.org/10.1111/mec.13243>.

Meirmans, Patrick G. 2020. “GenoDive Version 3.0: Easy-to-Use Software
for the Analysis of Genetic Data of Diploids and Polyploids.” *Molecular
Ecology Resources* 20 (4): 1126–31.
<https://doi.org/10.1111/1755-0998.13145>.

Meirmans, Patrick G., and Philip W. Hedrick. 2011. “Assessing Population
Structure: FST and Related Measures.” *Molecular Ecology Resources* 11
(1): 5–18. <https://doi.org/10.1111/j.1755-0998.2010.02927.x>.

Mikles, Chloe S., Camille M. L. S. Pagniello, Eyal Bigal, et al. 2026.
“Adaptive Genomic Divergence Parallels Migratory Behavior in Atlantic
Bluefin Tuna.” *Current Biology* 36: 2518–36.
<https://doi.org/10.1016/j.cub.2026.04.006>.

Paradis, Emmanuel. 2010. “Pegas: An r Package for Population Genetics
with an Integrated-Modular Approach.” *Bioinformatics* 26 (3): 419–20.
<https://doi.org/10.1093/bioinformatics/btp696>.

Phipson, Belinda, and Gordon K. Smyth. 2010. “Permutation p-Values
Should Never Be Zero: Calculating Exact p-Values When Permutations Are
Randomly Drawn.” *Statistical Applications in Genetics and Molecular
Biology* 9 (1): Article 39. <https://doi.org/10.2202/1544-6115.1585>.

Rochette, Nicolas C., Angel G. Rivera-Colón, and Julian M. Catchen.
2019. “Stacks 2: Analytical Methods for Paired-End Sequencing Improve
RADseq-Based Population Genomics.” *Molecular Ecology* 28 (21): 4737–54.
<https://doi.org/10.1111/mec.15253>.

Weir, Bruce S., and C. Clark Cockerham. 1984. “Estimating f-Statistics
for the Analysis of Population Structure.” *Evolution* 38 (6): 1358–70.
<https://doi.org/10.2307/2408641>.

Werth, Silke. 2011. “Optimal Sample Sizes and Allelic Diversity in
Studies of the Genetic Variability of Mycobiont and Photobiont
Populations.” *The Lichenologist* 43 (1): 73–81.
<https://doi.org/10.1017/S0024282910000563>.
