
# assigner <img src="docs/logo.png" align="right" alt="" />

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://tidyverse.org/lifecycle/#maturing)
[![Travis-CI Build
Status](https://travis-ci.org/thierrygosselin/assigner.svg?branch=master)](https://travis-ci.org/thierrygosselin/assigner)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/assigner?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/assigner)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/assigner)](http://cran.r-project.org/package=assigner)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/assigner.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/assigner)
[![packageversion](https://img.shields.io/badge/Package%20version-0.5.5-orange.svg)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2019--04--28-brightgreen.svg)](/commits/master)

-----

The name **assigner** |əˈsʌɪn| is rooted in the latin word *assignare*.
It’s first use in french dates back to XIIIe.

Genomic datasets produced by next-generation sequencing techniques that
reduce the size of the genome (e.g. genotype-by-sequencing (GBS) and
restriction-site-associated DNA sequencing (RADseq)) have a huge numbers
of markers that hold great potential and promises for assignment
analysis. After hitting the bioinformatic wall with the different
workflows you’ll likely end up with several folders containing whitelist
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

The **keywords** here to remember: 3 differents algorithms implemented
with frequentist, likelihood and the latest machine learning methods,
marker selection (with a fast Fst WC84 implementation), cross-validation
techniques (classic Leave-One-Out and Training, Holdout, Leave-one-out),
resampling/bootstrap/subsampling, ggplot2-based plotting\!

## Installation

To try out the dev version of **assigner**:

``` r
if (!require("pak")) install.packages("pak")
pak::pkg_install("thierrygosselin/assigner")
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
  - [Computer setup and
    troubleshooting](http://thierrygosselin.github.io/assigner/articles/rad_genomics_computer_setup.html)
  - [Function’s
    documentation](http://thierrygosselin.github.io/assigner/reference/index.html)
  - [Assigner’s
    features](http://thierrygosselin.github.io/assigner/FEATURES.html)
  - [Vignettes](http://thierrygosselin.github.io/assigner/articles/index.html)
  - How to cite assigner: inside R type `citation("assigner")`
