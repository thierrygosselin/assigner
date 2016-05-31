# assigner v.0.2.8

* While changing some lines with `tidyr::spread` and `tidyr::gather` for `data.table::dcast.data.table` and `data.table::melt.data.table` to make the code faster, I forgot to split genotype into alleles for `gsi_sim`.
* please update both **stackr** and **assigner**
* the build error from Travis will be fixed soon. 
It should not affect the package "experience"" in any way.

# assigner v.0.2.7

* you need to update [stackr] (https://github.com/thierrygosselin/stackr) to v.0.2.7
to appreciate this new version of assigner.

* updated `assignment_ngs` with the separate stackr modules to simplify the function.

* new data file available for `assignment_ngs`: `genepop` and `genind` object.

* `assignment_ngs` now accept any vcf input file! i.e. it’s no longer limited to stacks vcf. 

* new arguments in `assignment_ngs`. The assignment using dapc can now use 
the optimized alpha score `adegenet.dapc.opt == "optim.a.score"` or 
the cross-validation `adegenet.dapc.opt == "xval"`. This is useful for fine tuning 
the trade-off between power of discrimination and over-fitting 
(for stability of group membership probabilities).
Cross validation with `adegenet.dapc.opt == "xval"` doesn't work with 
missing data, so it's only available with *imputed data* 
(i.e. `imputation.method == "rf" or "max"`). 
With non imputed data or the default: the optimized alpha-score is used 
(`adegenet.dapc.opt == "optim.a.score"`).
When using `adegenet.dapc.opt == "xval"`, 2 new arguments are available: 
(1) `adegenet.n.rep` and (2) `adegenet.training`. See documentation for details.

* removed arguments in `assignment_ngs`. Removed the `pop.id.start` and `pop.id.end`
arguments that were confusing people. For those used to these arguments, 
they are now recycled in the new function `individuals2strata`
in [stackr] (https://github.com/thierrygosselin/stackr). 
The strata file created by this function can be used with the `strata` argument in
`assignment_ngs`.

* 2 modified arguments in `assignment_ngs`: (1) `gsi_sim.filename` is now `filename`; and 
(2) if you didn't use the imputation argument, replace `imputation.method = FALSE`
to `imputation.method = NULL` or leave the argument missing.

* simplified sections of codes in `assignment_ngs` that dealt with `strata`, 
`pop.levels` and `pop.labels`.

* new function: `write_gsi_sim`. Write a gsi_sim file from a data frame (wide or long/tidy). 
Used internally in [assigner] (https://github.com/thierrygosselin/assigner)
and might be of interest for users.

# assigner v.0.2.6

* Added a `NEWS.md` file to track changes to the package.
* `fst_WC84` is now a separate and very fast function that can compute: the overall and pairwise Weir and Cockerham 1984 Theta/Fst. Bootstrap resampling of markers is avalaible to build Confidence Intervals (For Louis Bernatchez and his students;). The estimates are available as a data frame and a matrix with upper diagonal filled with Fst values and lower diagonal filled with the confidence intervals.

# assigner v.0.2.5
* cleaner code for strata section
* bug fix restricted to `assignment_ngs` + `assignment.analysis = "adegenet"` + `sampling.method = "ranked"`.
A line at the beginning of a gsi_sim code section was deleted makig the assignment with adegenet go through that chunk of code and causing 100% assignment! if (assignment.analysis = "gsi_sim") {code} prevent this problem...

# assigner v.0.2.4
* bug fixed using adegenet that was introduced in v.0.2.3
* introducing a new function `import_subsamples_fst` to import the fst ranking 
results from all the subsample runs inside an assignment folder.

# assigner v.0.2.3
* bug fixed in the compilation results section when no pop.id.start and end are used.

# assigner v.0.2.1
* updated the function `assignment_mixture` with `sampling.method = "ranked"` and 
`assignment.analysis = "adegenet"`.

# assigner v.0.2.0
* new function: `assignment_mixture` for mixture analysis.

# assigner v.0.1.9
* Simplified gsi_sim install

# assigner v.0.1.8
* You can now opt between [gsi_sim] (https://github.com/eriqande/gsi_sim) or [adegenet] (https://github.com/thibautjombart/adegenet), a R package developed by Thibaul Jombart, to conduct the assignment analysis

# assigner v.0.1.7
* New input file: Re-introduced the haplotype data frame file from stacks.
* Argument name change: `imputations` is now `impute.method`.
* New argument: `impute` with 2 options: `impute = "genotype"` or `impute = "allele"`.

# assigner v.0.1.6
* Input file argument is now `data` and covers the three types of files the 
function can use: VCF file, PLINK tped/tfam or data frame of genotypes file.
* Huge number of markers (> 50 000 markers) can now be imported in PLINK 
tped/tfam format. The first 2 columns of the `tfam` file will be used for the 
`strata` argument, unless a new one is provided. Columns 1, 3 and 4 of the 
`tped` are discarded. The remaining columns correspond to the genotype in the 
format `01/04` where `A = 01, C = 02, G = 03 and T = 04`. For `A/T` format, use 
PLINK or bash to convert. Use [VCFTOOLS] (http://vcftools.sourceforge.net/) with 
`--plink-tped` to convert very large VCF file. For `.ped` file conversion to 
`.tped` use [PLINK] (http://pngu.mgh.harvard.edu/~purcell/plink/) with
`--recode transpose`.

# assigner v.0.1.5
* bug fix in `method = "random"` and `imputation`

# assigner v.0.1.4
* Changed function name, from `GBS_assignment` to `assignment_ngs`. Stands for
assignment with next-generation sequencing data.
* New argument `df.file` if you don't have a VCF file. See documentation.
* New argument `strata` if you don't have population id or other metadata info 
in the individual name.  See documentation.

# assigner v.0.1.3
* Changed arguments `THL` to `thl` and `snp.LD` to `snp.ld` to follow convention.
* `iterations.subsample` changed to `iteration.subsample`.
* `iterations` changed to `iteration.method` to avoid confusion with other iteration arguments.
* Removed `baseline` and `mixture` arguments from the function `GBS_assignment`.
These options will be re-introduce later in a separate function.
* Using `marker.number` higher than the number of markers in the data set was causing
problems. This could arise when using arguments that removed markers from the dataset
(e.g. `snp.ld`, `common.markers`, and `maf` filters).

# assigner v.0.1.2
* new version to update with gsi_sim new install instruction for Linux and Mac.
After re-installing assigner package, follow the instruction to re-install 
the new [gsi_sim] (https://github.com/eriqande/gsi_sim). 
And delete the old binary 'gsisim' in the /usr/local/bin folder
with the following Terminal command: `sudo rm /usr/local/bin/gsisim`