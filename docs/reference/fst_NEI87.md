# A fast implementation of Nei's 1987 Fst (overall and paiwise estimates)

The function calculates for diploid genomes the classic Nei's Gst
(1987), Nei's G'st (prime), that comes with a correction for the bias
that stems from sampling a limited number of populations. Also
calculated is Jost's D, an index of population differentiation that is
independent of the amount of within-population diversity (Jost, 2008).
Both overall and pairwise Fst can be estimated with confidence intervals
based on bootstrap of markers (resampling with replacement). The
function should give identical results *at the 4th decimal* when tested
against `genet.dist` in `hierfstat` and with the Fst computed in
`Calculate Distances` or
[GenoDive](http://www.bentleydrummer.nl/software/software/GenoDive.md).
The fastest computation is still
[GenoDive](http://www.bentleydrummer.nl/software/software/GenoDive.md),
but here, the R solution computes confidence intervals and it's very
fast. The computations takes advantage of tidyverse packages and
parallel. The impact of unbalanced design on estimates can be tested by
using the subsample argument.

*Special concerns for genome-wide estimate and filtering bias*

During computation, the function first starts by keeping only the
polymorphic markers in common between the populations. Keep this in mind
when filtering your markers to use this function characteristic
strategically to get better genome-wide estimate. This is even more
important when your project involves more than 2 populations that
evolved more by neutral processes (e.g. genetic drift) than by natural
selection (see the vignette for more details).

## Usage

``` r
fst_NEI87(
  data,
  pop.levels = NULL,
  pop.labels = NULL,
  strata = NULL,
  holdout.samples = NULL,
  pairwise = FALSE,
  ci = FALSE,
  iteration.ci = 100,
  quantiles.ci = c(0.025, 0.975),
  subsample = NULL,
  iteration.subsample = 1,
  digits = 9,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A file in the working directory or object in the global environment in
  wide or long (tidy) formats. To import, the function uses internally
  [genometranslator](https://github.com/thierrygosselin/genometranslator)
  [`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html).
  See details for more info.

  *How to get a tidy data frame ?*
  [genometranslator](https://github.com/thierrygosselin/genometranslator)
  [`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
  can import supported genomic data formats in a tidy data frame (VCF,
  PLINK, genind, genlight, gtypes, genepop, stacks haplotype file,
  hierfstat, ...). You can also use this function to filter your dataset
  using whitelist of markers, blacklist of individuals and genotypes.

- pop.levels:

  (optional, string) This refers to the levels in a factor. In this
  case, the id of the pop. Use this argument to have the pop ordered
  your way instead of the default alphabetical or numerical order. e.g.
  `pop.levels = c("QUE", "ONT", "ALB")` instead of the default
  `pop.levels = c("ALB", "ONT", "QUE")`. Default: `pop.levels = NULL`.

- pop.labels:

  (optional, string) Use this argument to rename/relabel your pop or
  combine your pop. e.g. To combine `"QUE"` and `"ONT"` into a new pop
  called `"NEW"`: (1) First, define the levels for your pop with
  `pop.levels` argument: `pop.levels = c("QUE", "ONT", "ALB")`. (2)
  then, use `pop.labels` argument:
  `pop.levels = c("NEW", "NEW", "ALB")`.#' To rename `"QUE"` to `"TAS"`:
  `pop.labels = c("TAS", "ONT", "ALB")`. Default: `pop.labels = NULL`.
  If you find this too complicated, there is also the `strata` argument
  that can do the same thing, see below.

- strata:

  (optional) A tab delimited file with 2 columns with header:
  `INDIVIDUALS` and `STRATA`. If a `strata` file is specified, the
  strata file will have precedence over any grouping found input file
  (`data`). The `STRATA` column can be any hierarchical grouping.
  Default: `strata = NULL`.

- holdout.samples:

  (optional) Samples that don't participate in the Fst computation
  (supplementary). Data frame with one column `INDIVIDUALS`. Default:
  `holdout.samples = NULL`.

- pairwise:

  (logical, optional) With `pairwise = TRUE`, the pairwise NEI87 Fst is
  calculated between populations. Default: `pairwise = FALSE`.

- ci:

  (logical, optional) Compute bootstrapped confidence intervals.
  Default: `ci = FALSE`.

- iteration.ci:

  (integer, optional) The number of iterations for the boostraps
  (resampling with replacement of markers). Default:
  `iteration.ci = 100`.

- quantiles.ci:

  (double, optional) The quantiles for the bootstrapped confidence
  intervals. Default: `quantiles.ci = c(0.025,0.975)`.

- subsample:

  (Integer or character) With `subsample = 36`, 36 individuals in each
  populations are chosen randomly to represent the dataset. With
  `subsample = "min"`, the minimum number of individual/population found
  in the data is used automatically. Default is no subsampling,
  `subsample = NULL`.

- iteration.subsample:

  (Integer) The number of iterations to repeat subsampling. With
  `subsample = 20` and `iteration.subsample = 10`, 20
  individuals/populations will be randomly chosen 10 times. Default:
  `iteration.subsample = 1`.

- digits:

  (optional, integer) The number of decimal places to be used in
  results. Default: `digits = 9`.

- parallel.core:

  (optional) The number of core for parallel computation of pairwise
  Fst. If not selected `detectCores() - 1` is used as default.

- verbose:

  (logical, optional) `verbose = TRUE` to be chatty during execution.
  Default: `verbose = FALSE`.

- ...:

  other parameters passed to the function.

## Value

The function returns a list with several objects. When sumsample is
selected the objects end with `.subsample`.

- `$subsampling.individuals`: the combinations of individuals and
  subsamples,

- `$fst.markers`: Nei's Gst, Nei's G'st and Jost's D by markers,

- `$fst.ranked`: Nei's Gst, Nei's G'st and Jost's D by markers ranked by
  Nei's Gst,

- `$fst.overall`: the overall Nei's Gst, Nei's G'st and Jost's D by
  markers with confidence intervals.

- `$fis.markers`: the fis by markers,

- `$fis.overall`: the mean fis overall markers with confidence intervals
  and the number of markers,

- `$fst.plot`: the histogram of the overall G'st per markers,

- `$pairwise.fst`: pairwise Nei's Gst, Nei's G'st and Jost's D in
  long/tidy data frame and the number of markers ,

- `$pairwise.fst.upper.matrix`: the pairwise fst prime in a upper
  triangle matrix,

- `$pairwise.fst.full.matrix`: the pairwise fst prime matrix (duplicated
  upper and lower triangle),

- `$pairwise.fst.ci.matrix`: matrix with pairwise fst prime in the upper
  triangle and the confidence intervals in the lower triangle.

- when subsample is selected `$pairwise.fst.subsample.mean` is a summary
  of all pairwise comparisons subsample. The mean is calculated accross
  summary statistics.

## Details

**Input data:**

To discriminate the long from the wide format, the function
genometranslator
[`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
searches for `MARKERS or LOCUS` in column names (TRUE = long format).
The data frame is tab delimitted. **Wide format:** The wide format
cannot store metadata info. The wide format starts with these 2 id
columns: `INDIVIDUALS`, `STRATA` (that refers to any grouping of
individuals), the remaining columns are the markers in separate columns
storing genotypes.

**Long/Tidy format:** The long format is considered to be a tidy data
frame and can store metadata info. (e.g. from a VCF; see
genometranslator
[`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)).
A minimum of 4 columns are required in the long format: `INDIVIDUALS`,
`STRATA`, `MARKERS or LOCUS` and `GENOTYPE or GT`. The rest are
considered metata info.

**2 genotypes formats are available:** 6 characters no separator: e.g.
`001002 of 111333` (for heterozygote individual). 6 characters WITH
separator: e.g. `001/002 of 111/333` (for heterozygote individual). The
separator can be any of these: `"/", ":", "_", "-", "."`.

*How to get a tidy data frame ?* genometranslator
[`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
can import supported genomic data formats in a tidy data frame.

## Note

**Negative Fst** are technical artifact of the computation (see Roesti
el al. 2012) and are automatically replaced with zero inside this
function.

**Why no p-values ?**

There is no null hypothesis testing with *P*-values. Confidence
intervals provided with the *F*-statistics enables more reliable
conclusions about the biological trends in the data.

## References

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press

Meirmans PG, Van Tienderen PH (2004) genotype and genodive: two programs
for the analysis of genetic diversity of asexual organisms. Molecular
Ecology Notes, 4, 792-794.

Roesti M, Salzburger W, Berner D. (2012) Uninformative polymorphisms
bias genome scans for signatures of selection. BMC Evol Biol., 12:94.
doi:10.1111/j.1365-294X.2012.05509.x

Jost L. (2008) G(ST) and its relatives do not measure differentiation.
Molecular Ecology. 17: 4015-4026.

## See also

From
[GenoDive](http://www.bentleydrummer.nl/software/software/GenoDive.md)
manual: *'In general, rather than to test differentiation between all
pairs of populations, it is adviseable to perform an overall test of
population differentiation, possibly using a hierarchical population
structure, (see AMOVA)'*

To compute an AMOVA, use
[GenoDive](http://www.bentleydrummer.nl/software/software/GenoDive.md)

[hierfstat](https://github.com/jgx65/hierfstat/)

Link for
[GenoDive](http://www.bentleydrummer.nl/software/software/GenoDive.md)

For Fisher's exact test and p-values per markers see `mmod` `diff_test`.

[`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
to import supported genomic data format in tidy data frames.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
 require(assigner)
 wombat.fst.pairwise <- fst_NEI87(
 data = "wombat.filtered.tidy.tsv",
 sep = "/",
 pop.levels = c("ATL", "MLE", "BIS", "PMO", "SOL", "TAS", "ECU"),
 holdout.samples = NULL,
 pairwise = TRUE,
 ci = TRUE,
 iteration.ci = 10000,
 quantiles.ci = c(0.025,0.975),
 parallel.core = 8,
 verbose = TRUE
 )
 #To get the overall Fst estimate:
 wombat.fst.pairwise$fst.overall
 #To get the Fst plot:
 wombat.fst.pairwise$fst.plots
 #To get the pairwise Fst values with confidence intervals in a data frame:
 wombat.fst.pairwise$pairwise.fst
 #To get the pairwise fst and ci matrix in a data frame:
 # rename, data frame, put rownames in column
 pairwise.fst.ci.df <- data.frame(pairwise.fst.ci.matrix) %>% add_rownames("POP")
} # }
```
