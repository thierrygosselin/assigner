[![Travis-CI Build Status](https://travis-ci.org/thierrygosselin/assigner.svg?branch=master)](https://travis-ci.org/thierrygosselin/assigner) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/assigner?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/assigner) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/assigner)](http://cran.r-project.org/package=assigner) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![DOI](https://zenodo.org/badge/14548/thierrygosselin/assigner.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/assigner)

[![packageversion](https://img.shields.io/badge/Package%20version-0.3.6-orange.svg)](commits/master) [![Last-changedate](https://img.shields.io/badge/last%20change-2016--11--14-brightgreen.svg)](/commits/master)

------------------------------------------------------------------------

assigner
========

This is the development page of the **assigner** package for the R software. The name **assigner** |əˈsʌɪn| is rooted in the latin word *assignare*. It's first use in french dates back to XIIIe.

Next-generation sequencing techniques that reduce the size of the genome (e.g. genotype-by-sequencing (GBS) and restriction-site-associated DNA sequencing (RADseq)) produce huge numbers of markers that hold great potential and promises for assignment analysis. After hitting the bioinformatic wall with the different workflows you'll likely end up with several folders containing whitelist and blacklist of markers and individuals, data sets with various *de novo* and/or filtering parameters and missing data. This reality of GBS/RADseq data is quite hard on GUI software traditionally used for assignment analysis. The end results is usually poor data exploration, constrained by time, and poor reproducibility.

**assigner** was tailored to make it easy to conduct assignment analysis using GBS/RADseq data within R. Additionally, combining the use of tools like [R Notebook](http://rmarkdown.rstudio.com/r_notebooks.html), [RStudio](https://www.rstudio.com) and [GitHub](https://github.com) will make effortless documenting your workflows and pipelines.

**assigner features:**

-   **Conduct assignment analysis** and **mixture analysis** using [gsi\_sim](https://github.com/eriqande/gsi_sim), a tool developed by Eric C. Anderson (see Anderson et al. 2008 and Anderson 2010) or [adegenet](https://github.com/thibautjombart/adegenet), a R package developed by Thibaul Jombart, to conduct the assignment analysis.
-   The input file are:
    1.  **VCF** (Danecek et al. 2011),
    2.  **PLINK** tped/tfam format (Purcell et al. 2007), usefull for very large files (&gt; 50 000 markers),
    3.  **genind** and **genlight** objets from [adegenet](https://github.com/thibautjombart/adegenet) (Jombart et al. 2010; Jombart and Ahmed, 2011)
    4.  **genepop** for traditionalists... (Raymond and Rousset 1995; Rousset 2008)
    5.  **haplotypes data frame file** produced by [STACKS](http://catchenlab.life.illinois.edu/stacks/) (Catchen et al. 2011, 2013). (e.g. *batch\_x.haplotypes.tsv*)
    6.  a data frame of genotypes in wide or long/tidy format
-   Markers can be randomly selected for a **classic LOO (Leave-One-Out) assignment** or chosen based on **ranked Fst** (Weir & Cockerham, 1984) for a **THL (Training, Holdout, Leave-one-out) assignment analysis** (reviewed in Anderson 2010)
-   Use `iteration.method` and/or `iteration.subsample` arguments to resample markers or individuals to get statistics!
-   assigner provides a fast implementation of Weir and Cockerham (1984) Fst/Theta and Nei's fst (1987). Both **overall** and **pairwise Fst** can be estimated with **confidence intervals** based on bootstrap of markers (resampling with replacement).
-   **Map-independent imputation** of missing genotype or alleles using **Random Forest** or the most frequent category is also available to test the impact of missing data on assignment analysis
-   **Filters:**
-   Individuals, populations and markers can be **filtered** and/or selected in several ways using **blacklist, whitelist** and other arguments
-   Genotypes of poor quality (e.g. in coverage, **genotype likelihood** or sequencing errors) can be erased prior to imputations or assignment analysis with the use of a `blacklist.genotype` argument
-   The impact of the **minor allele frequency**, MAF, (local and global) can also be easily explored with custom thresholds
-   Use the `whitelist.markers` argument to explore the impact of filters used in other software
-   Import and summarise the assignment results from [GenoDive](http://www.bentleydrummer.nl/software/software/GenoDive.html) (Meirmans and Van Tienderen, 2004)
-   Using [GenoDive](http://www.bentleydrummer.nl/software/software/GenoDive.html) output, compute the **genotype likelihood ratio distance metric (Dlr)** (Paetkau's et al. 1997, 2004)
-   `ggplot2`-based plotting to view assignment results and create publication-ready figures
-   Fast computations with optimized codes to run in parallel!

Installation
------------

To try out the dev version of **assigner**, follow the 3 steps below:

**Step 1:** Install or load the package **devtools**

``` r
if (!require("devtools")) install.packages("devtools") # to install
```

**Step 2:** Install **assigner** and install [gsi\_sim](https://github.com/eriqande/gsi_sim) from source

``` r
devtools::install_github("thierrygosselin/assigner", build_vignettes = TRUE)  # to install WITH vignettes
assigner::install_gsi_sim(fromSource = TRUE) # to install gsi_sim from source
```

**Step 3 (optional): Parallel computing** Install an OpenMP enabled [randomForestSRC](http://www.ccs.miami.edu/~hishwaran/rfsrc.html) package to do imputation in parallel. Follow the steps in this [vignette](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_imputations_parallel.Rmd). You don't need to do this when updating **assigner**.

\*\*<Notes:**>

-   **Problems during installation:** see this [vignette](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)
-   I recommend using [RStudio](https://www.rstudio.com/products/rstudio/download/) to run **assigner**. The R GUI is unstable with functions using parallel.
-   **Windows users**: 1. Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/). 2. To have *assigner* run in parallel, use [parallelsugar](https://github.com/nathanvan/parallelsugar). Easy to install and use ([instructions](https://github.com/nathanvan/parallelsugar#installation)). **Dependencies:**
-   **Imports:** adegenet, data.table, ggplot2, lazyeval, parallel, purrr, randomForestSRC, readr, stringi, stringr, tidyr, utils, plyr, dplyr (&gt;= 0.5.0), stackr (&gt;= 0.2.9)
-   **Suggests:** devtools, knitr, rmarkdown
-   **Remotes:** github::thierrygosselin/stackr

A quick way to install/load required packages and start using my packages (copy/paste the whole block):

``` r
if (!require("pacman")) install.packages("pacman")
library("pacman")
pacman::p_load(devtools, tidyverse, data.table, parallel, lazyeval, randomForestSRC)
pacman::p_load(devtools, tidyverse, data.table, parallel, lazyeval, randomForestSRC)
devtools::install_github("thierrygosselin/stackr", build_vignettes = TRUE)
library("stackr")
devtools::install_github("thierrygosselin/assigner", build_vignettes = TRUE)
library("assigner")
assigner::install_gsi_sim(fromSource = TRUE) 
}
```

Parallel computation and Cloud Computing
----------------------------------------

Most of the function in **assigner** were designed to be as fast as possible. Using computer with 16GB RAM is recommended. With more CPU and Memory comes faster computation time. If you decide to keep intermediate files during assignment analysis, you will need a large external drive (disk space is cheap). Solid State Drive and thunderbolt cables will provide fast input/output.

If disk space and computer power is an issue, cloud computing with [Google Cloud Compute Engine](https://cloud.google.com/compute/) and [Amazon Elastic Cloud Compute](https://aws.amazon.com/ec2/) is cheap and can be used easily.

A tutorial and pipeline along an Amazon Machine Image (AMI) are available in our [tutorial-workflow](http://gbs-cloud-tutorial.readthedocs.org/en/latest/).

The AMI is preloaded with **gsi\_sim** and the required **R packages**. Following a few steps: [link](http://gbs-cloud-tutorial.readthedocs.org/en/latest/10_use_rstudio.html), you can have [RStudio server](https://www.rstudio.com/) running and used through your web browser!

The Amazon image can be imported into Google Cloud Compute Engine to start a new compute engine virtual machine: [link](https://cloud.google.com/compute/docs/creating-custom-image#import_an_ami_image).

New features
------------

Change log, version, new features and bug history now lives in the [NEWS.md file](https://github.com/thierrygosselin/assigner/blob/master/NEWS.md)

**v.0.3.6**

-   bug fix `assignment_ngs` during imputation, the imputation module could not recognise that REF/ALT allele are not necessary or usefull for assignment analysis. \*enhancement to `assignment_ngs` and `assignment_mixture` so that when `marker.number` include `"all"` the `iteration.method` is set automatically to `1` when conducting the assignment with all the markers. Iterations at this point is useless.
-   random seed number is now stored in the appropriate files.
-   `assignment_mixture`: with `assignment.analysis = "gsi_sim` the unknown/mixture samples are compared with baseline populations using common markers between the pair. Now, the tables include the number of markers used. The summary provides the mean number of markers. This number will change each time randomness is used.

**v.0.3.5**

-   bug fix in population not recognise properly

**v.0.3.4**

-   `fst_NEI87`: very fast function that can compute: the overall and pairwise Nei's (1987) fst and f'st (prime). Bootstrap resampling of markers is avalaible to build Confidence Intervals. The estimates are available as a data frame and a matrix with upper diagonal filled with Fst values and lower diagonal filled with the confidence intervals. Jost's D is also given ;)

**v.0.3.3**

-   `fst_WC84`: bug fix, the function was not properly configured for multi-allelic markers (e.g. microsatellite, and haplotype format from STACKS). Thanks to Craig McDougall for catching this.

For previous news: [NEWS.md file](https://github.com/thierrygosselin/assigner/blob/master/NEWS.md)

Examples and Vignettes
----------------------

-   **From a browser:**
    -   [installation problems](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)
    -   [parallel computing during imputations](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_imputations_parallel.Rmd)
    -   [Pairwise and overall Fst with confidence intervals + building phylogenetic tree](https://github.com/thierrygosselin/assigner/blob/master/vignettes/fst_confidence_intervals.Rmd)
    -   Why use assigner to compute Fst of W&C (1984) in R ? [(Rmd Notebook)](https://github.com/thierrygosselin/assigner/blob/master/vignettes/fst_comparisons.Rmd) and [(html vignette)](https://www.dropbox.com/s/fl4s30enh6nnxiy/fst_comparisons.nb.html?dl=0)
    -   How to run several whitelists of markers in `assigner::assignment_ngs`[(html vignette)](https://www.dropbox.com/s/btt2oqx2qec601x/assignment.whitelists.Rmd?dl=0)
    -   **github** : [link](https://github.com/thierrygosselin/assigner/tree/master/vignettes)
-   **inside R**

    ``` r
    # The install must have the build_vignettes option:
    install_github("thierrygosselin/assigner", build_vignettes = TRUE)  # to install WITH vignettes
    ```

``` r
browseVignettes("assigner") # To browse vignettes
# For Pairwise and overall Fst with confidence intervals + building phylogenetic tree:
vignette("fst_confidence_intervals")
# To get top markers for assignment:
vignette("top_markers_assignment")
# To visualize individual's assignment results:
vignette("visualize_individual_assignment_results")
```

Vignettes are in development, check periodically for updates.

Roadmap of future developments:
-------------------------------

-   The ability to provide the ranking of markers based on other statistics, something else than the Fst (Weir and Cockerham, 1984) currently used in the function. e.g. Informativeness for Assignment Measure (In, Rosenberg et al. 2003), the Absolute Allele Frequency Differences (delta, δ, Shriver et al., 1997).
-   Provide ranking from other software: e.g. Toolbox for Ranking and Evaluation of SNPs [TRES](http://mlkd.csd.auth.gr/bio/tres/), [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) and [OutFLANK](https://github.com/whitlock/OutFLANK).
-   Would be very cool to use genotype likelihood information to get more accurate assignment.
-   Use Shiny and ggvis when subplots or facets are available
-   CRAN
-   ...suggestions ?

Contributions:
--------------

This package has been developed in the open, and it wouldn’t be nearly as good without your contributions. There are a number of ways you can help me make this package even better:

-   If you don’t understand something, please let me know.
-   Your feedback on what is confusing or hard to understand is valuable.
-   If you spot a typo, feel free to edit the underlying page and send a pull request.

New to pull request on github ? The process is very easy:

-   Click the edit this page on the sidebar.
-   Make the changes using github’s in-page editor and save.
-   Submit a pull request and include a brief description of your changes.
-   “Fixing typos” is perfectly adequate.

Citation:
---------

To get the citation for **assigner**, inside R:

``` r
citation("assigner")
```

References
----------

Anderson EC, Waples RS, Kalinowski ST (2008) An improved method for predicting the accuracy of genetic stock identification. Canadian Journal of Fisheries and Aquatic Sciences, 65, 1475–1486.

Anderson EC (2010) Assessing the power of informative subsets of loci for population assignment: standard methods are upwardly biased. Molecular Ecology Resources, 10, 701–710.

Catchen JM, Amores A, Hohenlohe PA et al. (2011) Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. G3, 1, 171–182.

Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) Stacks: an analysis tool set for population genomics. Molecular Ecology, 22, 3124–3140.

Danecek P, Auton A, Abecasis G et al. (2011) The variant call format and VCFtools. Bioinformatics, 27, 2156–2158.

Foll M, Gaggiotti O (2008) A Genome-Scan Method to Identify Selected Loci Appropriate for Both Dominant and Codominant Markers: A Bayesian Perspective. Genetics, 180, 977–993.

Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival, Regression and Classification (RF-SRC), R package version 1.6.1.

Jombart T (2008) adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24, 1403–1405.

Jombart T, Ahmed I (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics, 27, 3070–3071.

Kavakiotis I, Triantafyllidis A, Ntelidou D et al. (2015) TRES: Identification of Discriminatory and Informative SNPs from Population Genomic Data. Journal of Heredity, 106, 672–676.

Meirmans PG, Van Tienderen PH (2004) genotype and genodive: two programs for the analysis of genetic diversity of asexual organisms. Molecular Ecology Notes, 4, 792-794.

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press.

Paetkau D, Slade R, Burden M, Estoup A (2004) Genetic assignment methods for the direct, real-time estimation of migration rate: a simulation-based exploration of accuracy and power. Molecular Ecology, 13, 55-65.

Paetkau D, Waits LP, Clarkson PL, Craighead L, Strobeck C (1997) An empirical evaluation of genetic distance statistics using microsatellite data from bear (Ursidae) populations. Genetics, 147, 1943-1957.

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. American Journal of Human Genetics. 2007; 81: 559–575. <doi:10.1086/519795>

Raymond M. & Rousset F, (1995). GENEPOP (version 1.2): population genetics software for exact tests and ecumenicism. J. Heredity, 86:248-249

Rosenberg NA, Li LM, Ward R, Pritchard JK (2003) Informativeness of genetic markers for inference of ancestry. American Journal of Human Genetics, 73, 1402–1422.

Rousset F. genepop'007: a complete re-implementation of the genepop software for Windows and Linux. Molecular Ecology Resources. 2008, 8: 103-106. <doi:10.1111/j.1471-8286.2007.01931.x>

Shriver MD, Smith MW, Jin L et al. (1997) Ethnic-affiliation estimation by use of population-specific DNA markers. American Journal of Human Genetics, 60, 957.

Weir BS, Cockerham CC (1984) Estimating F-Statistics for the Analysis of Population Structure. Evolution, 38, 1358–1370.

Whitlock MC, Lotterhos KE (2015) Reliable Detection of Loci Responsible for Local Adaptation: Inference of a Null Model through Trimming the Distribution of FST\*. The American Naturalist, S000–S000.
