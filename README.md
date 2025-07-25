
# assigner <a href='http://thierrygosselin.github.io/assigner/'><img src='man/figures/logo.png' align="right" height="160" /></a>

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://tidyverse.org/lifecycle/#maturing)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.1.0-6666ff.svg)](https://cran.r-project.org/)
[![packageversion](https://img.shields.io/badge/Package%20version-0.7.0-orange.svg)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2025--07--01-brightgreen.svg)](/commits/master)
[![R-CMD-check](https://github.com/thierrygosselin/assigner/workflows/R-CMD-check/badge.svg)](https://github.com/thierrygosselin/assigner/actions)
[![macOS](https://github.com/thierrygosselin/assigner/workflows/macOS-latest%20(release)/badge.svg)](https://github.com/thierrygosselin/assigner/actions?workflow=macOS-latest%20(release))
[![Linux](https://github.com/thierrygosselin/assigner/workflows/check-linux//badge.svg)](https://github.com/thierrygosselin/assigner/actions?workflow=check-linux)
[![Windows](https://github.com/thierrygosselin/assigner/workflows/R-CMD-check/badge.svg)](https://github.com/thierrygosselin/assigner/actions?workflow=check-windows)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/assigner.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/assigner)
<!-- badges: end -->

The name **assigner** \|əˈsʌɪn\| is rooted in the Latin word
*assignare*. Its first use in French dates back to XIIIe.

For the logo, I was inspired by the [Northern Atlantic
Octupus](https://500px.com/photo/3885192/nature-a-reflexion-by-thierry-gosselin).
I was fortunate to spend a lot of times with one during my PhD. These
incredible creatures have 8 arms and [thousands of
suckers](https://500px.com/photo/3885141/stuck-on-me-by-thierry-gosselin)
that they can control independently. Octopus are really the best
multitaskers. The logo was designed by the artist [Claude
Thivierge](http://www.claudethivierge.com).

Genomic datasets produced by next-generation sequencing techniques that
reduce the size of the genome (e.g. genotype-by-sequencing (GBS) and
restriction-site-associated DNA sequencing (RADseq)) have a huge number
of markers that hold great potential and promises for assignment
analysis. After hitting the bioinformatic wall with the different
workflow, you’ll likely end up with several folders containing whitelist
and blacklist of markers and individuals, data sets with various *de
novo* and/or filtering parameters and … missing data. This reality of
GBS/RADseq data is quite hard on GUI software traditionally used for
population assignment analysis. The end results are usually poor data
exploration, constrained by time, and poor reproducibility.

**assigner** was tailored to make it easy to conduct population
assignment analysis using GBS/RADseq data within R. Additionally,
combining the use of tools like [R
Notebook](http://rmarkdown.rstudio.com/r_notebooks.html),
[RStudio](https://www.rstudio.com) and [GitHub](https://github.com) will
make effortless documenting your workflows and pipelines.

The **keywords** here to remember:

- 3 differents algorithms implemented: frequentist, likelihood and
  machine learning
- cross-validation techniques: classic Leave-One-Out (LOO) and Training,
  Holdout, Leave-one-out (THL) with marker selection
- resampling/bootstrap/subsampling
- fast Fst WC84 implementation)
- ggplot2-based plotting!
- <https://thierrygosselin.github.io/assigner/>

## Installation

To try out the dev version of **assigner**:

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("thierrygosselin/assigner")
library(assigner)
```

If you plan on using `gsi_sim` inside assigner, you need an additional
step:

With UNIX

``` r
assigner::install_gsi_sim(fromSource = TRUE)
```

With PC

``` r
assigner::install_gsi_sim()
```

- web site and additional info:
  <https://thierrygosselin.github.io/assigner/>
- [Computer setup - installation -
  troubleshooting](http://thierrygosselin.github.io/assigner/articles/rad_genomics_computer_setup.html)
- [assigner’s
  assumptions](http://thierrygosselin.github.io/assigner/reference/assignment_ngs.html#assumptions)
- [assigner’s
  features](http://thierrygosselin.github.io/assigner/FEATURES.html)
- [Function’s
  documentation](http://thierrygosselin.github.io/assigner/reference/index.html)
- [Vignettes](http://thierrygosselin.github.io/assigner/articles/index.html)
- How to cite assigner: inside R type `citation("assigner")`

## [Life cycle](https://thierrygosselin.github.io/assigner/articles/life_cycle.html)

assigner is maturing, but in order to make the package better, changes
are inevitable. Experimental functions will change, argument names will
change. Your codes and workflows might break from time to time until
**assigner is stable**. Consequently, depending on your tolerance to
change, assigner might not be for you.

- Philosophy, major changes and deprecated functions/arguments are
  documented in life cycle section of functions.
- The latest changes are documented
  ([here](https://thierrygosselin.github.io/assigner/articles/life_cycle.html))
  and in [changelog, versions, new features and bug
  history](https://thierrygosselin.github.io/assigner/news/index.html)
- [issues](https://github.com/thierrygosselin/assigner/issues/new/choose)
  and
  [contributions](https://github.com/thierrygosselin/assigner/issues/new/choose)
