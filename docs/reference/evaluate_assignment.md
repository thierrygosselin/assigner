# Population assignment analysis for genomic data

Evaluate population-assignment performance using random or FST-ranked
marker subsets, with optional individual subsampling and
cross-validation. Assignment assumptions are listed in the section
below.

- **Input file:** various file format are supported (see `data` argument
  below).

- **Cross-Validations:** Markers can be randomly selected for a classic
  LOO (Leave-One-Out) assignment or chosen based on ranked Fst for a thl
  (Training, Holdout, Leave-one-out) assignment analysis.

- **Assignment analysis:** conducted in
  [gsi_sim](https://github.com/eriqande/gsi_sim), developed by assigner
  co-author Eric C. Anderson for genetic stock identification and
  simulation, or [adegenet](https://github.com/thibautjombart/adegenet),
  an R package developed by Thibaul Jombart.

- **Parallel:** The assignment can be conduncted on multiple CPUs. The R
  GUI is unstable with this functions, I recommend using
  [RStudio](https://www.rstudio.com/products/rstudio/download/).

## Usage

``` r
evaluate_assignment(
  data,
  strata = NULL,
  pop.levels = NULL,
  assignment.analysis = c("gsi_sim", "adegenet"),
  markers.sampling = c("ranked", "random"),
  marker.number = "all",
  thl = 1,
  iteration.method = 10,
  subsample = NULL,
  iteration.subsample = 1,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1,
  ...
)
```

## Arguments

- data:

  Several input format are accepted. assigner uses genometranslator
  [`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
  to import the data. See function documentation for more details.

- strata:

  Optional strata data or filename passed to readers that support it.
  Default: `strata = NULL`.

- pop.levels:

  (optional, string) This refers to the levels in a factor. In this
  case, the id of the pop. Use this argument to have the pop ordered
  your way instead of the default alphabetical or numerical order. e.g.
  `pop.levels = c("QUE", "ONT", "ALB")` instead of the default
  `pop.levels = c("ALB", "ONT", "QUE")`. White spaces in population
  names are replaced by underscore. Default: `pop.levels = NULL`.

- assignment.analysis:

  (character) Assignment analysis conducted with
  `assignment.analysis = "gsi_sim"` or
  `assignment.analysis = "adegenet"`. See **Details** section below for
  installing [gsi_sim](https://github.com/eriqande/gsi_sim).

- markers.sampling:

  (character) 2 options for markers selection:

  1.  `markers.sampling == "random"` the subset of markers are selected
      randomly, this is the classic *Leave-One-Out* (LOO) assignment.

  2.  `markers.sampling == "ranked"` the subset of markers are first
      ranked based on an overall *decreasing* Fst values. The Fst is
      computed with
      [`fst_WC84`](https://thierrygosselin.github.io/assigner/reference/fst_WC84.md)
      function, that uses a fast implementation of Weir and Cockerham
      1984 Fst/Theta equations. This selection method is used during
      *Training-Holdout-Leave One Out* (thl) assignment. How many
      markers are selected is controlled with argument `thl`.

- marker.number:

  (Integer or string of number or "all") The assignment analysis can use
  all your markers (default) or a subset of your markers. e.g. To test
  500, 1000, 2000 and all the markers:
  `marker.number = c(500, 1000, 2000, "all")`. To use only 500 makers
  `marker.number = 500`. How those markers are sampled is determined
  with the argument `markers.sampling`, next. Default =
  `marker.number = "all"`.

- thl:

  (character, integer, proportion) For `markers.sampling = "ranked"`
  only. Several options are available:

  1.  `thl = 1`: Only 1 individual sample is used as holdout. This
      individual is not participating in the markers ranking. For each
      marker number, the analysis will be repeated with all the
      indiviuals in the data set (e.g. 500 individuals, 500 times 500,
      1000, 2000 markers). This is the default.

  2.  `proportion`: e.g. `thl = 0.15`, 15 percent of the individuals, in
      each strata, are chosen randomly as holdout individuals.

  3.  `thl = "all"`: all individuals are used for ranking (not good) and
      the argument `iteration.method = 1` is set by default. This is for
      testing only.

  Different lists of holdout individuals are generated when the argument
  `iteration.method` (bootstrap) is used.

- iteration.method:

  (integer) With **random markers selection** method, the iterations
  argument = the number of iterations to repeat marker resampling.
  Default: `iteration.method = 10`.

  With `marker.number = c(500, 1000)` and default iterations setting,
  500 markers will be randomly chosen 10 times and 1000 markers will be
  randomly chosen 10 times.

  **Notes:** If all the markers are used, with `marker.number = "all"`
  or in a series of marker number groupings
  `marker.number = c(200, 500, "all")`, the number of iteration is
  automatically set to 1. The remaining groupings are unaffected.

  With **ranked makers selection** method, using `thl = 1`, the analysis
  will be repeated for each individuals in the data set for every
  `marker.number` selected. With a proportion argument `thl = 0.15`, 15
  percent of individuals in each populations are chosen randomly as
  holdout individuals and this process is reapeated the number of times
  chosen by the `iteration.method` value.

- subsample:

  (Integer or Character, optional) This argument subsamples individuals
  to evaluate and reduce the influence of unequal reference-population
  sizes. Unequal sizes can make allele-frequency estimates more precise
  in large strata and bias assignment toward them. With
  `subsample = 36`, 36 individuals in each populations are chosen
  randomly to represent the dataset. This integer as to be smaller than
  the population with min sample size, if higher, the minimum sample
  size found in the data will be used as default. In doubt, use
  `subsample = "min"`, the function will use the smallest population
  sample size found in the data. The number of times this process is
  repeated is controlled by the argument `iteration.subsample`. Default:
  `subsample = NULL` (no subsampling). This is best used as an
  exploratory sensitivity analysis: repeat subsampling and compare
  assignment results with the full dataset. For low-coverage GL/PL data,
  equal numbers do not necessarily provide equal information; use the
  effective sample sizes returned by
  [`assign_individuals`](https://thierrygosselin.github.io/assigner/reference/assign_individuals.md)
  to guide population or coverage equalization.

- iteration.subsample:

  (Integer) The number of iterations to repeat subsampling of
  individuals. With `subsample = 20` and `iteration.subsample = 10`, 20
  individuals/populations will be randomly chosen 10 times. Default:
  `iteration.subsample = 1`.

- verbose:

  Logical. Display progress messages. For GDS input, the current number
  of samples and markers and a summary of active filters are displayed.
  Default: `verbose = TRUE`.

- parallel.core:

  Number of processor cores passed to readers that support parallel
  processing. Default: `parallel.core = parallel::detectCores() - 1`.

- ...:

  (optional) To pass further argument for fine-tuning the function (see
  advanced section below).

## Value

Depending on arguments selected, several folders and files are written:

1.  `assigner_evaluate_assignment_args_date@time.tsv`: For
    reproducibility, the function call, arguments and values used along
    the default arguments.

2.  `assignment.plot.pdf`: The assignment figure.

3.  `assignment.results.summary.stats.tsv`: tibble of summarized
    assignment statistics.

4.  `assignment.ranked.results.summary.stats.all.subsamples.tsv`: When
    subsampling is used, this file contains the raw results of all
    subsample before summarizing.

**THL: Training, Holdout, Leave-one-out**:

Intermediate files are written, you can inspect specific iterations
and/or subsample:

1.  `assignment.ranked.results.iterations.raw.tsv`: tibble with raw
    intermediate results of all iterations.

2.  `assignment.ranked.results.iterations.summary.tsv`: tibble with
    intermediate summary over iterations.

3.  `holdout.individuals.tsv`: tibble with the holdout individuals and
    associated iteration and random seed number.

**LOO: Leave-one-out**:

Intermediate files are written, you can inspect specific iterations
and/or subsample:

1.  `assignment.random.results.iterations.raw.tsv`: tibble with raw
    intermediate results of all iterations.

2.  `markers.random.tsv`: tibble with the random markers selected for
    each iteration with associated random seed number.

**Folders**

The names have the different iterations *i* starting with
`assignment_`*i* contains:

- `assignment_`*i*`.tsv`: the assignment result, for the iteration.

- `fst.ranked_`*i*`.tsv`: for THL method, the ranked Fst per markers,
  for the iteration.

- `gsi_sim_seeds`: the `gsi_sim` random seeds when this program is used,
  for the iteration.

The output in your global environment is a list. To view the assignment
results `$assignment` to view the ggplot2 figure `$assignment.plot`. See
example below.

## Marker selection and high-grading bias

Anderson (2010) demonstrated that assignment accuracy is biased upward
when the same individuals are used both to select an informative marker
panel and to evaluate that panel. Ordinary leave-one-out assignment does
not remove this high-grading bias because the focal individual has
already influenced marker ranking. With \`markers.sampling = "ranked"\`,
use a genuine Training-Holdout-Leave-one-out design. Holdout individuals
must not influence FST ranking, the number of markers selected,
thresholds, parameter tuning, or other analytical choices. Use \`thl =
"all"\` only for development or descriptive exploration, not to estimate
future assignment accuracy.

## Advance mode

Ideally, forget about this section. For advance users, through
*dots-dots-dots ...* you can pass several arguments for fine-tuning the
function:

1.  `adegenet.dapc.opt` (optional, character) Argument available only
    when using: `assignment.analysis = "adegenet"` with
    `markers.sampling == "random"`.

    The assignment using dapc can use the optimized alpha score
    `adegenet.dapc.opt == "optim.a.score"` or cross-validation
    `adegenet.dapc.opt == "xval"` for stability of group membership
    probabilities. For fine tuning the trade-off between power of
    discrimination and over-fitting. See adegenet documentation for more
    details. `adegenet.dapc.opt == "xval"` doesn't work with missing
    data, so it's only available with full dataset or **imputed
    dataset**. With non imputed data or the default:
    `adegenet.dapc.opt == "optim.a.score"`.

2.  `adegenet.n.rep`: (optional, integer) When
    `adegenet.dapc.opt == "xval"`, the number of replicates to be
    carried out at each level of PC retention. Default:
    `adegenet.n.rep = 30`. See adegenet documentation for more details.

3.  `adegenet.training`: (optional, numeric) When
    `adegenet.dapc.opt == "xval"`, the proportion of data (individuals)
    to be used for the training set. Default: `adegenet.training = 0.9`,
    if all groups have \>= 10 members. Otherwise, training.set scales
    automatically to the largest proportion that still ensures all
    groups will be present in both training and validation sets. See
    adegenet documentation for more details.

4.  `folder`: (optional) The name of the folder created in the working
    directory to save the files/results. Default: `folder = NULL` will
    create the folder for you with data and time.

5.  `filename`: (optional) The name of the file written to the
    directory. Use the extension ".txt" at the end. Several info will be
    appended to the name of the file. Default `assignment_data.txt`.

6.  `keep.gsi.files`: (logical, optional) With the default, the
    intermediate input and output gsi_sim files will be deleted from the
    directory when finished processing. I you decide to keep the files,
    with `keep.gsi.files = TRUE`, remember to allocate a large chunk of
    the disk space for the analysis. Default: `keep.gsi.files = FALSE`

7.  `random.seed`: (integer, optional) For reproducibility, set an
    integer that will be used inside function that requires randomness.
    With default, a random number is generated and printed in the
    appropriate output. Default: `random.seed = NULL`.

8.  `full.y.range`: (logical, optional) By default the y axis print min
    and max values are determied automatically from the data. It might
    be more usefull to see the full range of the scale from 0 to 100, in
    this case use `full.y.range = TRUE`. This can also be modified after
    running the analysis. See example below. Default:
    `full.y.range = FALSE`.

## Data quality and assignment assumptions

Assignment accuracy is meaningful only when it reflects biological
differentiation rather than technical or sampling structure. Manel et
al. (2005) emphasize matching the method and its assumptions to the
biological question. Before running or interpreting this evaluation,
examine:

- **Candidate sources and strata.** Candidate populations must be
  defined coherently and sampled representatively. An unsampled true
  source, pooled substructure, admixture, or an unsuitable choice of
  strata can produce a forced but misleading assignment.

- **Linkage disequilibrium.** Highly linked markers provide correlated
  evidence and can cause a genomic region to be counted many times. Thin
  or prune markers, work with justified linkage blocks, and compare
  results across alternative panels. Any marker ranking or pruning used
  to evaluate performance must be learned from training data only.

- **Hardy-Weinberg departure.** Strong departure can indicate genotyping
  error, paralogs, inbreeding, family structure, admixture, or
  incorrectly combined populations. Diagnose the cause rather than
  applying an automatic Hardy-Weinberg filter.

- **Genotyping error.** Allelic dropout, low depth, duplicated or
  paralogous loci, inconsistent allele coding, and platform-specific
  errors can create false population signals. Use technical replicates
  and inspect depth, allele balance, reproducibility, and error
  patterns.

- **Missingness and batch effects.** Missingness, sequencing depth,
  extraction batch, library, lane, and genotype-processing history
  should not be confounded with strata. A random holdout can preserve
  this confounding and report excellent but non-transferable accuracy.
  When possible, evaluate on a separate batch, collection, site, or time
  period.

- **Duplicates and relatedness.** Duplicate samples or close relatives
  shared between training and holdout data leak information and inflate
  performance. Split related groups together when family structure is
  present.

- **Rare variants and sample size.** Low minor-allele counts are
  unstable when LOO or THL removes individuals, and small or uneven
  source samples yield imprecise allele frequencies and asymmetric
  power. Examine sensitivity to marker and sample composition.

- **Independent validation.** Report confusion among every source, not
  only overall success. Repeat the analysis across random seeds and
  marker panels, and prefer external validation when the intended
  samples differ from the reference data in space, time, laboratory
  batch, or sequencing technology.

This function does not silently perform general genomic quality control.
Read and standardize data with
[genometranslator](https://thierrygosselin.github.io/genometranslator/),
then inspect and filter it with
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow before assignment analysis.

## Life cycle

Map-independent imputation of missing genotype is avaible in my other R
package called [grur](https://github.com/thierrygosselin/grur).

Use [grur](https://github.com/thierrygosselin/grur) to :

1.  **Visualize your missing data:** before imputing your genotypes,
    visualize your missing data. Several visual tools are available
    inside [grur](https://github.com/thierrygosselin/grur) to help you
    decide the best strategy after.

2.  **Optimize:** use [grur](https://github.com/thierrygosselin/grur)
    imputation module and other functions to optimize the imputations of
    your dataset. You need to test arguments. Failing to conduct tests
    and adjust imputations arguments will **generate artifacts** and/or
    **exacerbate bias**. Using defaults is not optional here...

3.  **genome_translator:** use the output argument inside
    [grur](https://github.com/thierrygosselin/grur) imputation module to
    generate the required formats for assigner (e.g. a tidy dataset)

**Deprecated arguments:**

- `sampling.method`: renamed `markers.sampling`.

## Dependencies

The `adegenet` engine is supplied by the CRAN package of the same name.
The alternative [gsi_sim](https://github.com/eriqande/gsi_sim) engine,
developed by assigner co-author Eric C. Anderson, is an external
command-line program and is not installed by assigner. Follow the
complete source-installation instructions in the package get-started
vignette and validate the executable with
[`check_gsi_sim`](https://thierrygosselin.github.io/assigner/reference/check_gsi_sim.md)
before selecting `assignment.analysis = "gsi_sim"`.

## References

Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008) An
improved method for predicting the accuracy of genetic stock
identification. Canadian Journal of Fisheries and Aquatic Sciences 65,
7:1475-1486.

Anderson EC (2010). Assessing the power of informative subsets of loci
for population assignment: standard methods are upwardly biased.
Molecular Ecology Resources, 10(4), 701-710.
[doi:10.1111/j.1755-0998.2010.02846.x](https://doi.org/10.1111/j.1755-0998.2010.02846.x)
.

Weir BS, Cockerham CC (1984) Estimating F-Statistics for the Analysis of
Population Structure. Evolution, 38, 1358-1370.

Jombart T, Devillard S, Balloux F. Discriminant analysis of principal
components: a new method for the analysis of genetically structured
populations. BMC Genet. 2010:11: 94. doi:10.1186/1471-2156-11-94

Jombart T, Ahmed I. adegenet 1.3-1: new tools for the analysis of
genome-wide SNP data. Bioinformatics. 2011:27: 3070-3071.
doi:10.1093/bioinformatics/btr521

Manel S, Gaggiotti OE, Waples RS (2005). Assignment methods: matching
biological questions with appropriate techniques. Trends in Ecology &
Evolution, 20(3), 136-142. doi:10.1016/j.tree.2004.12.004.

Paetkau D, Slade R, Burden M, Estoup A (2004). Genetic assignment
methods for the direct, real-time estimation of migration rate: a
simulation-based exploration of accuracy and power. Molecular Ecology,
13, 55-65. doi:10.1046/j.1365-294X.2003.02008.x.

Moran BM, Anderson EC (2019). Bayesian inference from the conditional
genetic stock identification model. Canadian Journal of Fisheries and
Aquatic Sciences, 76(4), 551-560. doi:10.1139/cjfas-2018-0016.

## See also

\[assign_individuals()\],
[gsi_sim](https://github.com/eriqande/gsi_sim), and
[rubias](https://github.com/eriqande/rubias)

\`gsi_sim\`, developed by assigner co-author Eric C. Anderson, is used
here to evaluate individual assignment and marker-panel performance.
\`rubias\`, developed by Eric C. Anderson and Ben Moran, addresses the
related but distinct conditional GSI problem of jointly estimating
mixed-stock proportions and individual posterior assignments. See
\[assign_individuals()\] for a fuller comparison.

## Author

Thierry Gosselin <thierrygosselin@icloud.com> and Eric C. Anderson

## Examples

``` r
if (FALSE) { # \dontrun{
assignment.treefrog <- evaluate_assignment(
 data = "batch_1.vcf",
 strata = "strata.treefrog.tsv",
 assignment.analysis = "gsi_sim",
 marker.number = c(500, 5000, "all"),
 markers.sampling = "ranked",
 thl = 0.3
 )
 # To create a dataframe with the assignment results:
   assignment <- assignment.treefrog$assignment
 # To plot the assignment using ggplot2 and facet
   fig <- assignment.treefrog$assignment.plot
 # To view the full range of y values = Assignment success(%):
   fig + ggplot2::scale_y_continuous(limits = c(0,100))

 # Or use the ... argument: full.y.range = TRUE
} # }
```
