# Assigner features

## Features

| Caracteristics | Description |
|:-------------------|:--------------------------------------------------------|
|**Analysis**| [gsi_sim](https://github.com/eriqande/gsi_sim) (Anderson et al., 2008; Anderson, 2010) or [adegenet](https://github.com/thibautjombart/adegenet) (Jombart, 2008; Jombart and Ahmed, 2011)<br><br>Conducted on: samples of **known origin (population)** or **unknown/mixture**|
|**Input files**| assigner relies on `radiator::tidy_genomic_data` for imports |
|**Markers selections**|Random or ranked based on W&C's Fst (Weir and Cockerham, 1984) or Nei's Gst (Nei, 1987)|
|**Cross-Validations**|Avoid high-grading bias, use **classic LOO** (Leave-One-Out) or **THL** (Training, Holdout, Leave-one-out) (see Anderson, 2010)|
|**Resampling/Bootstrap/Subsampling**|markers and individuals|
|**Filters**| `strata`: use the strata to **whitelist samples** and/or **select populations** `??radiator::read_strata`<br>`whitelist.markers`: kleep favorite markers<br>`marker.number`: test different groupings of marker numbers automatically (selected randomly or based on the ranking method described in the table above)<br><br>For more info on how to filter RADseq data, see [radiator](https://github.com/thierrygosselin/radiator).|
|**[ggplot2](http://ggplot2.org)-based plotting**|View assignment results and create publication-ready figures|
|**Parallel**|Codes designed and optimized for fast computations running imputations, iterations, etc. in parallel. Works with all OS: Linux, Mac and PC.|


## Independent modules

* **Fst**: assigner provides **one of the fastest** implementation of Weir and Cockerham (1984) Fst/Theta and Nei's fst (1987). Both **overall** and **pairwise Fst** can be estimated with **confidence intervals** based on bootstrap of markers (resampling with replacement).
* **Dlr:** Compute the genotype likelihood ratio distance metric (Paetkau's et al. 1997, 2004).


## Parallel computation and Cloud Computing

Most of the function in **assigner** were designed to be as fast as possible. 
Using computer with 16GB RAM is recommended. 
With more CPU and Memory comes faster computation time. 
If you decide to keep intermediate files during assignment analysis, 
you will need a large external drive (disk space is cheap). 
Solid State Drive and thunderbolt cables will provide fast input/output. 

If disk space and computer power is an issue, cloud computing with [DigitalOcean](https://www.r-bloggers.com/a-new-image-on-digitalocean-to-start-using-rstudio-server-without-waiting-more-than-2-minutes-2/), [Google Cloud Compute Engine](https://cloud.google.com/compute/) or 
[Amazon Elastic Cloud Compute](https://aws.amazon.com/ec2/) is cheap and can be used easily. 

A tutorial and pipeline along an Amazon Machine Image (AMI) are available 
in our [tutorial-workflow](http://gbs-cloud-tutorial.readthedocs.org/en/latest/). 

The AMI is preloaded with **gsi_sim** and the required **R packages**. 
Following a [few steps](http://gbs-cloud-tutorial.readthedocs.org/en/latest/10_use_rstudio.html), 
you can have [RStudio server](https://www.rstudio.com/) running and used through your web browser!

The Amazon image can be [imported into Google Cloud Compute Engine to start a new compute engine virtual machine](https://cloud.google.com/compute/docs/creating-custom-image#import_an_ami_image). 

## References

Anderson EC, Waples RS, Kalinowski ST (2008) An improved method for predicting the accuracy of genetic stock identification. Canadian Journal of Fisheries and Aquatic Sciences, 65, 1475–1486.

Anderson EC (2010) Assessing the power of informative subsets of loci for population assignment: standard methods are upwardly biased. Molecular Ecology Resources, 10, 701–710.

Catchen JM, Amores A, Hohenlohe PA et al. (2011) Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences. G3, 1, 171–182.

Catchen JM, Hohenlohe PA, Bassham S, Amores A, Cresko WA (2013) Stacks: an analysis tool set for population genomics. Molecular Ecology, 22, 3124–3140.

Danecek P, Auton A, Abecasis G et al. (2011) The variant call format and VCFtools. Bioinformatics, 27, 2156–2158.

Foll M, Gaggiotti O (2008) A Genome-Scan Method to Identify Selected Loci Appropriate for Both Dominant and Codominant Markers: A Bayesian Perspective. Genetics, 180, 977–993.

Goudet J.hierfstat, a package for r to compute and test hierarchical F-statistics.
Molecular Ecology Notes.2005; 5: 184-186. doi:10.1111/j.1471-8286.2004.00828.x


Ishwaran H. and Kogalur U.B. (2015). Random Forests for Survival, Regression and Classification (RF-SRC), R package version 1.6.1.

Jombart T (2008) adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24, 1403–1405.

Jombart T, Ahmed I (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics, 27, 3070–3071.

Kavakiotis I, Triantafyllidis A, Ntelidou D et al. (2015) TRES: Identification of Discriminatory and Informative SNPs from Population Genomic Data. Journal of Heredity, 106, 672–676.

Meirmans PG, Van Tienderen PH (2004) genotype and genodive: two programs for the analysis of genetic diversity of asexual organisms. Molecular Ecology Notes, 4, 792-794.

Nei M. (1987) Molecular Evolutionary Genetics. Columbia University Press.

Paetkau D, Slade R, Burden M, Estoup A (2004) Genetic assignment methods for the direct, real-time estimation of migration rate: a simulation-based exploration of accuracy and power. Molecular Ecology, 13, 55-65.

Paetkau D, Waits LP, Clarkson PL, Craighead L, Strobeck C (1997) An empirical evaluation of genetic distance statistics using microsatellite data from bear (Ursidae) populations. Genetics, 147, 1943-1957.

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, et al. 
PLINK: a tool set for whole-genome association and population-based linkage analyses. 
American Journal of Human Genetics. 2007; 81: 559–575. doi:10.1086/519795

Raymond M. & Rousset F, (1995). GENEPOP (version 1.2): population genetics 
software for exact tests and ecumenicism. J. Heredity, 86:248-249

Rosenberg NA, Li LM, Ward R, Pritchard JK (2003) Informativeness of genetic markers for inference of ancestry. American Journal of Human Genetics, 73, 1402–1422.

Rousset F. genepop'007: a complete re-implementation of the genepop software
for Windows and Linux. Molecular Ecology Resources. 2008, 8: 103-106. 
doi:10.1111/j.1471-8286.2007.01931.x

Shriver MD, Smith MW, Jin L et al. (1997) Ethnic-affiliation estimation by use of population-specific DNA markers. American Journal of Human Genetics, 60, 957.

Weir BS, Cockerham CC (1984) Estimating F-Statistics for the Analysis of Population Structure. Evolution, 38, 1358–1370.

Whitlock MC, Lotterhos KE (2015) Reliable Detection of Loci Responsible for Local Adaptation: Inference of a Null Model through Trimming the Distribution of FST*. The American Naturalist, S000–S000.

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS. A high-performance computing toolset for relatedness and principal component analysis of SNP data. Bioinformatics. 2012;28: 3326-3328. doi:10.1093/bioinformatics/bts606
