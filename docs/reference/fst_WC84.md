# A fast implementation of Weir and Cockerham (1984) Fst/Theta (overall and paiwise estimates)

The function calculates Weir and Cockerham (1984) Fst for diploid
genomes. Both overall and pairwise Fst can be estimated with confidence
intervals based on bootstrap of markers (resampling with replacement).
The function gives identical results *at the 9th decimal* when tested
against `genet.dist` in `hierfstat`. Using the argument
`snprelate = TRUE` will compute the Fst with
[SNPRelate](https://github.com/zhengxwen/SNPRelate). This implementation
gives slightly upward bias values but provided the fastest computations
I know, but it doesn't compute confidence intervals, for now. For an R
implementation, `fst_WC84` is very fast. The computations takes
advantage of dplyr, tidyr, purrr, parallel and SNPRelate. The impact of
unbalanced design on estimates can be tested by using the subsample
argument (see advance mode section).

*Special concerns for genome-wide estimate and filtering bias*

During computation, the function first starts by keeping only the
polymorphic markers in common between the populations. Keep this in mind
when filtering your markers to use this function characteristic
strategically to get better genome-wide estimate. This is even more
important when your project involves more than 2 populations that
evolved more by neutral processes (e.g. genetic drift) than by natural
selection (see the vignette for more details).

*Fst along linkage groups*

For reference-guided genomes with informative `CHROM` and `POS`
metadata, the function automatically reports WC84 Fst for each linkage
group and in non-overlapping genomic windows. The windowed estimates are
displayed in a chromosome-concatenated Manhattan-style plot. With
exactly two strata, these are pairwise estimates; with more strata, they
are overall estimates among all included strata. Broad, spatially
coherent regions of elevated differentiation can help identify regions
for follow-up, including candidate inversions. An Fst pattern alone does
not establish an inversion. Candidate regions should be evaluated with
complementary evidence such as absolute divergence, nucleotide
diversity, Tajima's D, linkage disequilibrium, local PCA, read mapping,
or structural-variant analyses. De novo datasets using synthetic
`CHROM_1` metadata are not plotted automatically.

The `"kernel"` option follows the Gaussian smoothing philosophy used by
the STACKS `populations` program: estimates are evaluated along the
reference genome, marker contributions decline with distance, and the
kernel is truncated at three standard deviations on either side. Here,
the weights are applied to the WC84 variance components before Fst is
calculated.

## Usage

``` r
fst_WC84(
  data,
  snprelate = FALSE,
  strata = NULL,
  pop.levels = NULL,
  pairwise = FALSE,
  ci = FALSE,
  iteration.ci = 100,
  quantiles.ci = c(0.025, 0.975),
  heatmap.fst = FALSE,
  linkage.group = NULL,
  window.size = 25000L,
  window.method = c("fixed", "sliding", "kernel"),
  window.step = NULL,
  kernel.sigma = 150000,
  outlier.sd = 4,
  digits = 9,
  filename = "fst_WC84",
  parallel.core = parallel::detectCores() - 2,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html).
  You can also use this function to filter your dataset using whitelist
  of markers, blacklist of individuals and genotypes.

- snprelate:

  (optional, logical) Use
  [SNPRelate](https://github.com/zhengxwen/SNPRelate) to compute the
  Fst. It's the fastest computation I've seen so far!

  However, testing with different RADseq datasets as shown several
  upward bias with
  [`SNPRelate::snpgdsFst`](https://rdrr.io/pkg/SNPRelate/man/snpgdsFst.html)
  (last version tested was v.1.16.0). I compared the results with
  assigner, hierfstat and strataG (results available upon request). The
  SNPRelate author as not given me good reason to belive the issue is
  fully resolved, consequently, the option is no longer available, until
  further notice. Default: `snprelate = FALSE`

- strata:

  (optional, data frame) A tab delimited file with 2 columns with
  header: `INDIVIDUALS` and `STRATA`. If a `strata` file is specified,
  the strata file will have precedence over any grouping found data file
  (`data`). The `STRATA` column can be any hierarchical grouping.
  Default: `strata = NULL`.

- pop.levels:

  (optional, string) This refers to the levels in a factor. In this
  case, the id of the pop. Use this argument to have the pop ordered
  your way instead of the default alphabetical or numerical order. e.g.
  `pop.levels = c("QUE", "ONT", "ALB")` instead of the default
  `pop.levels = c("ALB", "ONT", "QUE")`. Default: `pop.levels = NULL`.

- pairwise:

  (optional, logical) With `pairwise = TRUE`, the pairwise WC84 Fst is
  calculated between populations. Default: `pairwise = FALSE`.

- ci:

  (optional, logical) Compute bootstrapped confidence intervals.
  Default: `ci = FALSE`.

- iteration.ci:

  (optional, integer) The number of iterations for the boostraps
  (resampling with replacement of markers). Default:
  `iteration.ci = 100`.

- quantiles.ci:

  (optional, double) The quantiles for the bootstrapped confidence
  intervals. Default: `quantiles.ci = c(0.025,0.975)`.

- heatmap.fst:

  (logical) Generate a heatmap with the Fst values in lower matrix and
  CI in the upper matrix. The heatmap can also be generated separately
  after the Fst analysis using the separate function:
  [`heatmap_fst`](https://thierrygosselin.github.io/assigner/reference/heatmap_fst.md).
  Default: `heatmap.fst = FALSE`.

- linkage.group:

  (optional, logical) Generate chromosome or linkage-group summaries and
  a genome-position plot. With `linkage.group = NULL`, the feature is
  activated automatically when informative `CHROM` and numeric `POS`
  metadata indicate a reference-guided genome. Use `TRUE` to force it
  for a one-linkage-group genome, or `FALSE` to disable it. Default:
  `linkage.group = NULL`.

- window.size:

  (optional, integer) Width of each fixed, non-overlapping genomic
  window, measured in base pairs (bp). WC84 variance components from
  markers inside each window are recombined to estimate window-level
  Fst. For example, `window.size = 25000L` creates adjacent 25 kb
  windows: 1-25,000 bp, 25,001-50,000 bp, and so forth. This argument
  does not currently define a sliding or overlapping window. Default:
  `window.size = 25000L` (25 kb).

- window.method:

  (optional, character) Genomic window method. `"fixed"` uses adjacent,
  non-overlapping windows and is the default. `"sliding"` uses
  overlapping rectangular windows of `window.size` bp, advanced by
  `window.step` bp. `"kernel"` evaluates Gaussian-weighted WC84
  components on a regular grid; markers within three standard deviations
  of each grid position contribute to the estimate. Default:
  `window.method = "fixed"`.

- window.step:

  (optional, integer) Distance in base pairs between consecutive window
  centres for `window.method = "sliding"` or `"kernel"`. With `NULL`,
  the step is 20 percent of `window.size` for sliding windows and
  `window.size` for kernel smoothing. It is ignored for fixed windows.
  Default: `window.step = NULL`.

- kernel.sigma:

  (optional, numeric) Standard deviation of the Gaussian kernel in base
  pairs when `window.method = "kernel"`. Markers farther than
  `3 * kernel.sigma` from a grid position receive zero weight. Following
  the STACKS default, `kernel.sigma = 150000` uses markers up to 450 kb
  on either side. Default: `kernel.sigma = 150000`.

- outlier.sd:

  (optional, numeric) Number of standard deviations above the mean
  window Fst used for the exploratory horizontal outlier threshold. Use
  `NULL` to omit the threshold. Default: `outlier.sd = 4`.

- digits:

  (optional, integer) The number of decimal places to be used in
  results. Default: `digits = 9`.

- filename:

  (optional, character) Give filename prefix for the output directory,
  this will trigger saving results. Default: `filename = "fst_WC84"`.

- parallel.core:

  (optional, integer) The number of core for parallel computation of
  pairwise Fst. See also the advance mode section below. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  (optional, logical) `verbose = TRUE` to be chatty during execution.
  Default: `verbose = FALSE`.

- ...:

  other parameters passed to the function.

## Value

The function returns a list with several objects. When sumsample is
selected the objects end with `.subsample`.

- `$subsampling.individuals`: the combinations of individuals and
  subsamples,

- `$sigma.loc`: the variance components per locus, with (`lsiga`: among
  populations, `lsigb`: among individuals within populations, `lsigw`:
  within individuals)

- `$fst.markers`: the fst by markers,

- `$fst.ranked`: the fst ranked,

- `$fst.overall`: the mean fst overall markers and the number of markers

- `$fis.markers`: the fis by markers,

- `$fis.overall`: the mean fis overall markers and the number of
  markers,

- `$fst.plot`: the histogram of the overall Fst per markers,

- `$fst.linkage.groups`: descriptive per-linkage-group summaries,

- `$fst.windows`: WC84 estimates in non-overlapping genomic windows,

- `$fst.genome.plot`: windowed WC84 Fst in a chromosome-concatenated
  genome plot,

- `$pairwise.fst`: the pairwise fst in long/tidy data frame and the
  number of markers ,

- `$pairwise.fst.upper.matrix`: the pairwise fst in a upper triangle
  matrix,

- `$pairwise.fst.full.matrix`: the pairwise fst matrix (duplicated upper
  and lower triangle),

- `$pairwise.fst.ci.matrix`: matrix with pairwise fst in the upper
  triangle and the confidence intervals in the lower triangle.

- when subsample is selected `$pairwise.fst.subsample.mean` is a summary
  of all pairwise comparisons subsample. The mean is calculated accross
  summary statistics.

## Note

**Negative Fst** are technical artifact of the computation (see Roesti
el al. 2012) and are automatically replaced with zero inside this
function.

**Why no p-values ?**

There is no null hypothesis testing with *P*-values and its rarely if
ever the appropriate model in population genomics, despite its
popularity with molecular ecologists interested in population
differentiation.

"The important scientific question is the real magnitude of the
differentiation, not the smallness of the *P* value" -Lou Jost

Confidence intervals provided with the *F*-statistics enables more
reliable conclusions about the biological trends in the data. A
confidence interval describes the statistical uncertainty of the
*F*-statistics estimate. If the confidence interval include zero, then
the null hypothesis cannot be rejected. If the confidence interval does
not include zero, the null hypothesis can be rejected and you can also
have an appreciation of the real magnitude of the statistical
differentiation whether its large or small.

## Advance mode

*dots-dots-dots ...* allows to pass several arguments for fine-tuning
the function:

1.  `filter.monomorphic` (logical, optional) By default monomorphic
    markers present in the dataset are removed (and it should stay that
    way...). Default: `filter.monomorphic = TRUE`.

2.  `holdout.samples` (optional, data frame) Samples that don't
    participate in the Fst computation (supplementary). Data frame with
    one column `INDIVIDUALS`. This argument is used inside assignment
    analysis. Default: `holdout.samples = NULL`.

3.  `subsample` (Integer or character) With `subsample = 36`, 36
    individuals in each populations are chosen randomly to represent the
    dataset. With `subsample = "min"`, the minimum number of
    individual/population found in the data is used automatically.
    Default is no subsampling, `subsample = NULL`.

4.  `iteration.subsample` (Integer) The number of iterations to repeat
    subsampling. With `subsample = 20` and `iteration.subsample = 10`,
    20 individuals/populations will be randomly chosen 10 times.
    Default: `iteration.subsample = 1`.

5.  `calibrate.alleles` (logical) Un-calibrated alleles can bias
    estimate and by default the function expect that the REF/ALT alleles
    are calibrated. Using `calibrate.alleles = TRUE`, can take a bit
    more time. Default: `calibrate.alleles = FALSE`.

## References

Excoffier L, Smouse PE, Quattro JM. Analysis of molecular variance
inferred from metric distances among DNA haplotypes: application to
human mitochondrial DNA restriction data. Genetics. 1992;131: 479-491.

Jost L, D vs. GST: Response to Heller and Siegismund (2009) and Ryman
and Leimar (2009). Molecular Ecology. 2009; 18:10 2088-2091.
https://doi.org/10.1111/j.1365-294x.2009.04186.x

Meirmans PG, Van Tienderen PH (2004) genotype and genodive: two programs
for the analysis of genetic diversity of asexual organisms. Molecular
Ecology Notes, 4, 792-794.

Michalakis Y, Excoffier L. A generic estimation of population
subdivision using distances between alleles with special reference for
microsatellite loci. Genetics. 1996;142: 1061-1064.

Weir BS, Cockerham CC (1984) Estimating F-Statistics for the Analysis of
Population Structure. Evolution, 38, 1358-1370.

Roesti M, Salzburger W, Berner D. (2012) Uninformative polymorphisms
bias genome scans for signatures of selection. BMC Evolutionary Biology,
12, 94.
[doi:10.1186/1471-2148-12-94](https://doi.org/10.1186/1471-2148-12-94) .

Catchen J, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013). Stacks:
an analysis tool set for population genomics. Molecular Ecology, 22,
3124-3140. [doi:10.1111/mec.12354](https://doi.org/10.1111/mec.12354) .

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS. A
high-performance computing toolset for relatedness and principal
component analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
doi:10.1093/bioinformatics/bts606

## See also

From
[GenoDive](http://www.bentleydrummer.nl/software/software/GenoDive.md)
manual: *'In general, rather than to test differentiation between all
pairs of populations, it is advisable to perform an overall test of
population differentiation, possibly using a hierarchical population
structure, (see AMOVA)'*. To compute an AMOVA, use
[GenoDive](http://www.bentleydrummer.nl/software/software/GenoDive.md)
or `Phi_st_Meirmans` in `mmod`.

[hierfstat](https://github.com/jgx65/hierfstat/)

For Fisher's exact test and p-values per markers see `mmod` `diff_test`.

**Vignette for this function:** [how to do the pairwise and overall Fst
with confidence intervals and build the phylogenetic
tree](https://www.dropbox.com/s/tiq4yenzmgzc2f5/fst_confidence_intervals.html?dl=0)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
wombat.fst.pairwise <- fst_WC84(
    data = "wombat.filtered.tidy.tsv",
    pop.levels = c("ATL", "MLE", "BIS", "PMO", "SOL", "TAS", "ECU"),
    pairwise = TRUE,
    ci = TRUE,
    iteration.ci = 10000,
    quantiles.ci = c(0.025,0.975),
    parallel.core = 8,
    verbose = TRUE,
    filename = "wombat",
    heatmap.fst = TRUE
)

# To get the overall Fst estimate:
wombat.fst.pairwise$fst.overall

# To get the Fst plot:
wombat.fst.pairwise$fst.plot

#To get the pairwise Fst values with confidence intervals in a data frame:
df <- wombat.fst.pairwise$pairwise.fst
} # }
```
