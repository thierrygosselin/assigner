---
output: github_document
---

# assigner <a href="https://thierrygosselin.github.io/assigner/"><img src="man/figures/logo.png" align="right" height="160" alt="assigner logo" /></a>



<!-- badges: start -->
[![Project status: active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Package version](https://img.shields.io/badge/package%20version-0.7.1-orange.svg)](https://github.com/thierrygosselin/assigner)
[![R version](https://img.shields.io/badge/R-%3E%3D%204.4.0-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/14548/thierrygosselin/assigner.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/assigner)
<!-- badges: end -->

`assigner` provides reproducible population-assignment workflows and fast FST
estimators. It accepts genomic data prepared with
[genometranslator](https://thierrygosselin.github.io/genometranslator/) and
assumes quality control and filtering were completed first, for example with
[radr](https://thierrygosselin.github.io/radr/).

`evaluate_assignment()` supports random marker sampling and FST-ranked
Training–Holdout–Leave-one-out analyses using `adegenet` or the external
command-line program `gsi_sim`.

`assign_individuals()` provides a native likelihood-based assignment engine
implemented entirely in R. It retains the likelihood of every individual under
every candidate stratum, and its output can be passed directly to `dlr()`.

## Installation

```r
install.packages("remotes")
remotes::install_github("thierrygosselin/assigner")
```

`gsi_sim` is optional and is not installed by the R package. Complete
installation instructions are in the
[get-started vignette](https://thierrygosselin.github.io/assigner/articles/get_started.html).
After installing it, verify the setup with:

```r
assigner::check_gsi_sim()
```

## Start here

- [Get started](https://thierrygosselin.github.io/assigner/articles/get_started.html)
- [Function reference](https://thierrygosselin.github.io/assigner/reference/index.html)
- [FST confidence intervals](https://thierrygosselin.github.io/assigner/articles/fst_confidence_intervals.html)
- [Issue tracker](https://github.com/thierrygosselin/assigner/issues)

Record the package version with `packageVersion("assigner")` and obtain the
recommended citation with `citation("assigner")`.
