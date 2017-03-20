
<!-- README.md is generated from README.Rmd. Please edit that file -->
assigner<img src="assigner_logo.png" align="right"/>
====================================================

[![Travis-CI Build Status](https://travis-ci.org/thierrygosselin/assigner.svg?branch=master)](https://travis-ci.org/thierrygosselin/assigner) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/thierrygosselin/assigner?branch=master&svg=true)](https://ci.appveyor.com/project/thierrygosselin/assigner) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/assigner)](http://cran.r-project.org/package=assigner) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![DOI](https://zenodo.org/badge/14548/thierrygosselin/assigner.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/assigner)

[![packageversion](https://img.shields.io/badge/Package%20version-0.4.2-orange.svg)](commits/master) [![Last-changedate](https://img.shields.io/badge/last%20change-2017--03--20-brightgreen.svg)](/commits/master)

------------------------------------------------------------------------

This is the development page of the **assigner** package for the R software. The name **assigner** |əˈsʌɪn| is rooted in the latin word *assignare*. It's first use in french dates back to XIIIe.

Next-generation sequencing techniques that reduce the size of the genome (e.g. genotype-by-sequencing (GBS) and restriction-site-associated DNA sequencing (RADseq)) produce huge numbers of markers that hold great potential and promises for assignment analysis. After hitting the bioinformatic wall with the different workflows you'll likely end up with several folders containing whitelist and blacklist of markers and individuals, data sets with various *de novo* and/or filtering parameters and missing data. This reality of GBS/RADseq data is quite hard on GUI software traditionally used for assignment analysis. The end results is usually poor data exploration, constrained by time, and poor reproducibility.

**assigner** was tailored to make it easy to conduct assignment analysis using GBS/RADseq data within R. Additionally, combining the use of tools like [R Notebook](http://rmarkdown.rstudio.com/r_notebooks.html), [RStudio](https://www.rstudio.com) and [GitHub](https://github.com) will make effortless documenting your workflows and pipelines.

### Features:

<table style="width:100%;">
<colgroup>
<col width="26%" />
<col width="73%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Caracteristics</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><strong>Analysis</strong></td>
<td align="left"><a href="https://github.com/eriqande/gsi_sim">gsi_sim</a> (Anderson et al., 2008; Anderson, 2010) or <a href="https://github.com/thibautjombart/adegenet">adegenet</a> (Jombart, 2008; Jombart and Ahmed, 2011)<br><br>Conducted on: samples of <strong>known origin (population)</strong> or <strong>unknown/mixture</strong></td>
</tr>
<tr class="even">
<td align="left"><strong>Input files</strong></td>
<td align="left"><a href="https://samtools.github.io/hts-specs/">VCF</a> (Danecek et al., 2011), <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#tr">PLINK tped/tfam</a> (Purcell et al., 2007), <a href="https://github.com/thibautjombart/adegenet">adegenet genind and genlight</a> (Jombart et al., 2010; Jombart and Ahmed, 2011), <a href="https://github.com/EricArcher/strataG">strataG gtypes</a>, <a href="http://genepop.curtin.edu.au">Genepop</a> (Raymond and Rousset, 1995; Rousset, 2008), <a href="http://catchenlab.life.illinois.edu/stacks/">STACKS haplotype file</a> (Catchen et al., 2011, 2013), dataframes of genotypes in wide or long/tidy format</td>
</tr>
<tr class="odd">
<td align="left"><strong>Markers selections</strong></td>
<td align="left">Random or ranked based on W&amp;C's Fst (Weir and Cockerham, 1984) or Nei's Gst (Nei, 1987)</td>
</tr>
<tr class="even">
<td align="left"><strong>Cross-Validations</strong></td>
<td align="left"><strong>classic LOO</strong> (Leave-One-Out) or <strong>THL</strong> (Training, Holdout, Leave-one-out) (see Anderson, 2010)</td>
</tr>
<tr class="odd">
<td align="left"><strong>Resampling/Bootstrap/Subsampling</strong></td>
<td align="left">markers and individuals</td>
</tr>
<tr class="even">
<td align="left"><strong>Imputations</strong></td>
<td align="left"><strong>Map-independent</strong> imputations of missing genotypes.<br>Using <strong>Random Forest</strong> or the most frequent category.<br> Imputations can be conducted <strong>overall samples</strong> or <strong>by populations</strong>.</td>
</tr>
<tr class="odd">
<td align="left"><strong>Filters</strong></td>
<td align="left">Whitelist and blacklist arguments provide the way to select/filter(include/exclude):<br>genotypes, markers, individuals and populations.<br><br>See the dedicated filtering arguments on the next table below.</td>
</tr>
<tr class="even">
<td align="left"><strong><a href="http://ggplot2.org">ggplot2</a>-based plotting</strong></td>
<td align="left">View assignment results and create publication-ready figures</td>
</tr>
<tr class="odd">
<td align="left"><strong>Parallel</strong></td>
<td align="left">Codes designed and optimized for fast computations running imputations, iterations, etc. in parallel. Works with all OS: Linux, Mac and PC.</td>
</tr>
</tbody>
</table>

#### Specialized filtering arguments tailored for RADseq data and assignment analysis:

<table style="width:76%;">
<colgroup>
<col width="23%" />
<col width="52%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Filters</th>
<th align="left">Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><strong><code>pop.select</code></strong></td>
<td align="left">Keep selected populations</td>
</tr>
<tr class="even">
<td align="left"><strong><code>blacklist.id</code></strong></td>
<td align="left">Discard individuals</td>
</tr>
<tr class="odd">
<td align="left"><strong><code>whitelist.markers</code></strong></td>
<td align="left">Keep favorite markers</td>
</tr>
<tr class="even">
<td align="left"><strong><code>blacklist.genotype</code></strong></td>
<td align="left">Genotypes of poor quality (e.g. in <strong>coverage</strong>, <strong>genotype likelihood</strong> or <strong>sequencing errors</strong>) can be <strong>erased</strong> prior to imputations or assignment analysis</td>
</tr>
<tr class="odd">
<td align="left"><strong><code>maf.thresholds</code></strong></td>
<td align="left">Test different thresholds of the minor allele frequency (local and global MAF)</td>
</tr>
<tr class="even">
<td align="left"><strong><code>common.markers</code></strong></td>
<td align="left">Keep only markers in common between populations (optimized for pairwise analysis)</td>
</tr>
<tr class="odd">
<td align="left"><strong><code>monomorphic.out</code></strong></td>
<td align="left">Discard monomorphic markers</td>
</tr>
<tr class="even">
<td align="left"><strong><code>snp.ld</code></strong></td>
<td align="left">Prune markers based on short-distance linkage disequilibrium</td>
</tr>
<tr class="odd">
<td align="left"><strong><code>marker.number</code></strong></td>
<td align="left">Test different groupings of marker numbers automatically (selected randomly or based on the ranking method described in the table above)</td>
</tr>
<tr class="even">
<td align="left"><strong><code>More filters?</code></strong></td>
<td align="left">For more info on how to filter RADseq data, see <a href="https://github.com/thierrygosselin/stackr">stackr</a></td>
</tr>
</tbody>
</table>

#### Independent modules:

-   **Fst**: assigner provides **one of the fastest** implementation of Weir and Cockerham (1984) Fst/Theta and Nei's fst (1987). Both **overall** and **pairwise Fst** can be estimated with **confidence intervals** based on bootstrap of markers (resampling with replacement). See the \[vignettes for this below\].(<https://github.com/thierrygosselin/assigner#examples-and-vignettes>)
-   **AMOVA:** Analysis of Molecular Variance that run fast and works with missing values.
-   **Dlr:** Compute the genotype likelihood ratio distance metric (Paetkau's et al. 1997, 2004).

Installation
------------

To try out the dev version of **assigner**, follow the 3 steps below:

**Step 1:** Install or load the package **devtools** and **bioconductor** (for SNPRelate)

``` r
if (!require("devtools")) install.packages("devtools") # to install
source("https://bioconductor.org/biocLite.R") # for bioconductor dependencies
biocLite() # for bioconductor dependencies
```

**Step 2:** Install **assigner** and install [gsi\_sim](https://github.com/eriqande/gsi_sim) from source

``` r
devtools::install_github("thierrygosselin/assigner", build_vignettes = TRUE)  # to install WITH vignettes
assigner::install_gsi_sim(fromSource = TRUE) # to install gsi_sim from source
```

**Step 3 (optional): Parallel computing** Install an OpenMP enabled [randomForestSRC](http://www.ccs.miami.edu/~hishwaran/rfsrc.html) and [data.table](https://github.com/Rdatatable/data.table) packages (e.g. to do imputations in parallel). Follow the steps in this [vignette](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_imputations_parallel.Rmd). You don't need to do this when updating **assigner**.

**Notes**

-   **Problems during installation:** see this [vignette](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)
-   I recommend using [RStudio](https://www.rstudio.com/products/rstudio/download/) to run **assigner**. The R GUI is unstable with functions using parallel.
-   **Windows users**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

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

**v.0.4.2**

-   until `SNPRelate` [bias issue](https://github.com/zhengxwen/SNPRelate/issues/21) is resolved the option is unavailable
-   better use of `pbmcapply` for Windows
-   **imputations** is being reworked and will be buggy until the next update. The codes are being completely re-written and arguments will change (for the better).

**v.0.4.1**

-   debug code to work in parallel with Windows
-   code cleaning to prep for CRAN

**v.0.4.0**

-   `assignment_ngs` and `assignment_mixture` code cleaning to prep for CRAN and make them easier to debug.

**v.0.3.9**

-   I'm pleased to announce that `assigner` now works in parallel with **Windows**
-   bug fix introduce in last commit in `write_gsi_sim` where the file was not created properly from an internal module.

For previous news: [NEWS.md file](https://github.com/thierrygosselin/assigner/blob/master/NEWS.md)

Examples and Vignettes
----------------------

-   **From a browser:**
    -   [installation problems](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_installation_problems.Rmd)
    -   [parallel computing during imputations](https://github.com/thierrygosselin/stackr/blob/master/vignettes/vignette_imputations_parallel.Rmd)
    -   [Pairwise and overall Fst with confidence intervals + building phylogenetic tree](https://github.com/thierrygosselin/assigner/blob/master/vignettes/fst_confidence_intervals.Rmd)
    -   Why use assigner to compute Fst of W&C (1984) in R ? [(Rmd Notebook)](https://www.dropbox.com/s/djxvudh2amx2ayw/fst_comparisons.Rmd?dl=0) and [(html vignette)](https://www.dropbox.com/s/fl4s30enh6nnxiy/fst_comparisons.nb.html?dl=0)
    -   How to run several whitelists of markers in `assigner::assignment_ngs`[(html vignette)](https://www.dropbox.com/s/btt2oqx2qec601x/assignment.whitelists.Rmd?dl=0)
    -   **github** : [link](https://github.com/thierrygosselin/assigner/tree/master/vignettes)
-   **inside R** During installation of the package you need to select the buil vignette option: `install_github("thierrygosselin/assigner", build_vignettes = TRUE)`

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

Goudet J.hierfstat, a package for r to compute and test hierarchical F-statistics. Molecular Ecology Notes.2005; 5: 184-186. <doi:10.1111/j.1471-8286.2004.00828.x>

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

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS. A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics. 2012;28: 3326-3328. <doi:10.1093/bioinformatics/bts606>
