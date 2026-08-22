
# assigner <a href="https://thierrygosselin.github.io/assigner/"><img src="man/figures/logo.png" align="right" height="160" alt="assigner logo" /></a>

<!-- badges: start -->

[![Project status:
active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Package
version](https://img.shields.io/badge/package%20version-0.7.1-orange.svg)](https://github.com/thierrygosselin/assigner)
[![R
version](https://img.shields.io/badge/R-%3E%3D%204.4.0-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/assigner.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/assigner)
<!-- badges: end -->

`assigner` helps answer a practical question: **how well can genomic
data assign individuals to candidate source populations?**

Use `assigner` to:

- assign individuals with a native population-genetic likelihood method
  written entirely in R, using called genotypes (`GT`) or genotype
  likelihoods (`GL` and `PL`);
- retain the likelihood of every individual in every candidate
  population, not only the best and second-best assignments;
- compare likelihood assignment with elastic net, XGBoost, experimental
  TabPFN, `adegenet`, or Eric C. Anderson’s
  [`gsi_sim`](https://github.com/eriqande/gsi_sim) workflows;
- evaluate assignment honestly with holdout samples, leave-one-out
  estimation, and repeated marker panels;
- avoid the high-grading bias described by Anderson (2010) by keeping
  marker selection inside the training data and holdout individuals
  untouched;
- compare random marker panels with panels ranked by population
  differentiation;
- calculate and visualize Dlr as a measure of assignment resolution
  between population pairs;
- quantify effective reference sample size for GL/PL data and identify
  potentially important information imbalances among source populations;
- calculate an optional coverage-aware diagnostic for individuals that
  may originate from an unsampled source;
- estimate global and pairwise FST with confidence intervals; and
- keep analyses reproducible with recorded arguments, random seeds,
  intermediate results, and dated output folders.

Ordinary leave-one-out assignment does not remove the high-grading bias
created when the same individuals were already used to select apparently
informative loci. See [Anderson
(2010)](https://doi.org/10.1111/j.1755-0998.2010.02846.x) and the
get-started vignette for the Training-Holdout-Leave-one-out design.

`assigner` accepts genomic data prepared with
[genometranslator](https://thierrygosselin.github.io/genometranslator/)
and expects genomic quality control and filtering to have been completed
first, for example with [radr](https://thierrygosselin.github.io/radr/).

The package distinguishes assignment to predefined reference populations
from clustering, mixture estimation, and statistical migrant detection.
The [get-started
vignette](https://thierrygosselin.github.io/assigner/articles/get_started.html#theory-and-interpretation)
explains these distinctions, leave-one-out assignment, unsampled
alleles, missing data, validation, and the foundational references
behind the methods.

## Installation

``` r
install.packages("remotes")
remotes::install_github("thierrygosselin/assigner")
```

[`gsi_sim`](https://github.com/eriqande/gsi_sim), developed by assigner
co-author Eric C. Anderson, is optional and is not installed by the R
package. Complete installation instructions are in the [get-started
vignette](https://thierrygosselin.github.io/assigner/articles/get_started.html).
After installing it, verify the setup with:

``` r
assigner::check_gsi_sim()
```

## Start here

- [Get
  started](https://thierrygosselin.github.io/assigner/articles/get_started.html)
- [Function
  reference](https://thierrygosselin.github.io/assigner/reference/index.html)
- [FST confidence
  intervals](https://thierrygosselin.github.io/assigner/articles/fst_confidence_intervals.html)
- [Development
  roadmap](https://github.com/thierrygosselin/assigner/blob/main/ROADMAP.md)
- [Issue tracker](https://github.com/thierrygosselin/assigner/issues)

Record the package version with `packageVersion("assigner")` and obtain
the recommended citation with `citation("assigner")`.
