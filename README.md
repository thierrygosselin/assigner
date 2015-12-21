# assigner
[![Travis-CI Build Status](https://travis-ci.org/thierrygosselin/assigner.svg?branch=master)](https://travis-ci.org/thierrygosselin/assigner) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/assigner)](http://cran.r-project.org/package=assigner) [![DOI](https://zenodo.org/badge/14548/thierrygosselin/assigner.svg)](https://zenodo.org/badge/latestdoi/14548/thierrygosselin/assigner)



Next-generation sequencing techniques that reduce the size of the genome (e.g. genotype-by-sequencing (GBS) and
restriction-site-associated DNA sequencing (RADseq)) produce huge numbers of markers that hold great potential and promises for assignment analysis. After hitting the bioinformatic wall with the different workflows you'll likely end up with several folders containing whitelist and blacklist of markers and individuals, data sets with various de novo and/or filtering parameters and missing data. This reality of GBS/RAD data is quite hard on GUI software traditionally used for assignment analysis. The end results is usually poor data exploration, constrained by time, and poor reproducibility.

**assigner** was tailored to make it easy to conduct assignment analysis using GBS/RAD data within R. Additionally, combining the use of tools like [RStudio] (https://www.rstudio.com) and [GitHub] (https://github.com) will make effortless documenting your workflows and pipelines.


This is the development page of the **assigner** package for the R software.

**Use assigner to:**

* conduct assignment analysis using [gsi_sim] (https://github.com/eriqande/gsi_sim), a tool developed 
by Eric C. Anderson or [adegenet] (https://github.com/thibautjombart/adegenet), a R package
developed by Thibaul Jombart
* The input file is a the VCF (*batch_x.vcf*) or haplotype (*batch_1.haplotypes.tsv*) file, produced by [STACKS] (http://catchenlab.life.illinois.edu/stacks/)
* Individuals, populations and markers can be filtered and/or selected in several ways using blacklist, 
whitelist and other arguments
* Map-independent imputation of missing genotype using Random Forest or the most frequent category is also available to test the impact of missing data on assignment analysis
* Genotypes of poor quality (e.g. in coverage, genotype likelihood or sequencing errors) can be erased prior to imputations or assignment analysis with the use of a `blacklist.genotype` argument.
* Markers can be randomly selected for a classic LOO (Leave-One-Out) assignment or 
chosen based on ranked Fst (Weir & Cockerham, 1984) for a THL (Training, Holdout, Leave-one-out) assignment analysis
* use `iterations` argument to resample markers or individuals to get statistics!
* The impact of the minor allele frequency (local and global) can also be easily explored with custom thresholds
* `ggplot2`-based plotting to view results of assignment and create publication-ready figures


## Installation
You can try out the dev version of **assigner**. Follow the 2 steps below:

Step 1 You will need the package *devtools*
```r
install.packages("devtools") # to install
library(devtools) # to load
```

Step 2 Install **assigner**:
```r
install_github("thierrygosselin/assigner") # to install
library(assigner) # to load
```

**Note:**

Sometimes you'll get warnings while installing dependencies required for **assigner** or other R packages.
```r
Warning: cannot remove prior installation of package ‘stringi’
```

To solve this problem: 

Option 1. Delete the problematic packages manually and reinstall. On MAC computers, in the **Finder**, use the shortcut **cmd+shift+g**, or in the menu bar : **GO -> Go to Folder**, copy/paste the text below:
```r
/Library/Frameworks/R.framework/Resources/library
#Delete the problematic packages.
```

Option 2. If you know your way around the terminal and understand the consequences of using **sudo rm -R** command, here something faster to remove problematic packages:
```r
sudo rm -R /Library/Frameworks/R.framework/Resources/library/package_name
# Changing 'package_name' to the problematic package.
# Reinstall the package.
```

**Dependencies**

Here the list of packages that **assigner** is depending on:
```r
dplyr, reshape2, ggplot2, readr, stringi, tidyr, purrr, lazyeval, adegenet, randomForestSRC, foreach, parallel, doParallel
```
If you don't have them, no worries, it's intalled automatically during **assigner** installation. If you have them, it's your job to update them, because i'm usually using the latest versions...


## New

**v.0.1.0**


## Roadmap of future developments:

* The ability to provide the ranking of markers based on something else than Fst provided by the function. For example, Toolbox for Ranking and Evaluation of SNPs [TRES](http://mlkd.csd.auth.gr/bio/tres/), [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) and [OutFLANK](https://github.com/whitlock/OutFLANK).
* Documentation and vignettes
* Workflow tutorial
* Use Shiny and ggvis when subplots or facets are available
* CRAN
* ...suggestions ?

## Contributions:

This package has been developed in the open, and it wouldn’t be nearly as good without your contributions. There are a number of ways you can help me make this package even better:  
* If you don’t understand something, please let me know. 
* Your feedback on what is confusing or hard to understand is valuable. 
* If you spot a typo, feel free to edit the underlying page and send a pull request.

New to pull request on github ? The process is very easy:  
* Click the edit this page on the sidebar.
* Make the changes using github’s in-page editor and save.
* Submit a pull request and include a brief description of your changes. 
* “Fixing typos” is perfectly adequate.


## Example 
Under construction


Vignettes are in development, check periodically for updates.
