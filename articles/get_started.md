# Get started with assigner

## Package roles

`assigner` evaluates population-assignment performance from genomic
data. Use
[genometranslator](https://thierrygosselin.github.io/genometranslator/)
to read genomic formats and
[radr](https://thierrygosselin.github.io/radr/) to explore and filter
them. `assigner` does not silently perform general genomic quality
control.

## Install assigner

``` r

install.packages("remotes")
remotes::install_github("thierrygosselin/assigner")
library(assigner)
```

The `adegenet` engine is an R dependency. The alternative `gsi_sim`
engine, developed by assigner co-author Eric C. Anderson, is an external
command-line program and is not bundled with `assigner`.

## Install gsi_sim

`gsi_sim`, developed by Eric C. Anderson, is not currently available
through Bioconda. Install it from its [source
repository](https://github.com/eriqande/gsi_sim) inside the shared
`genomics` Conda environment.

### macOS

``` bash
conda activate genomics
conda install -c conda-forge git clang make
cd ~/Downloads
git clone --recurse-submodules https://github.com/eriqande/gsi_sim.git
cd gsi_sim
./Compile_gsi_sim.sh
mkdir -p "$CONDA_PREFIX/bin"
cp gsi_sim-Darwin "$CONDA_PREFIX/bin/gsi_sim"
chmod +x "$CONDA_PREFIX/bin/gsi_sim"
gsi_sim --help
```

### Linux

``` bash
conda activate genomics
conda install -c conda-forge git gcc make
cd ~/Downloads
git clone --recurse-submodules https://github.com/eriqande/gsi_sim.git
cd gsi_sim
./Compile_gsi_sim.sh
mkdir -p "$CONDA_PREFIX/bin"
cp gsi_sim-Linux "$CONDA_PREFIX/bin/gsi_sim"
chmod +x "$CONDA_PREFIX/bin/gsi_sim"
gsi_sim --help
```

Start R from an activated `genomics` environment, then verify discovery:

``` r

assigner::check_gsi_sim()
```

The upstream repository includes an old Windows binary. The modern
`assigner` workflow expects an executable named `gsi_sim` on `PATH`;
Windows users should use a compatible Unix-like environment and confirm
it with
[`check_gsi_sim()`](https://thierrygosselin.github.io/assigner/reference/check_gsi_sim.md).

## Prepare data and strata

Start with maintained sample metadata. `INDIVIDUALS` must match the
genomic data, and `STRATA` identifies candidate source populations.
Import source data with
[`genometranslator::read_genome()`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html),
perform quality control with `radr`, and pass the prepared GDS or tidy
genomic data to `assigner`.

``` r

data("data_assigner_sim_01", package = "assigner")
genome <- data_assigner_sim_01
```

## Theory and interpretation

Population assignment is a **classification** problem when candidate
source populations are defined in advance and represented by reference
samples. For each individual, the likelihood method evaluates how
probable its multilocus genotype is under the allele frequencies
estimated in every candidate source. The source with the greatest
likelihood is the assignment, but all source likelihoods remain
scientifically useful and are retained by `assigner`.

This framework relies on more than the numerical maximum. Reference
samples should represent the candidate sources, allele-frequency
estimates should be reliable, markers should carry information that
distinguishes sources, and the relevant sources should have been
sampled. Manel, Gaggiotti, and Waples (2005) emphasize starting with the
biological question because individual assignment, population exclusion,
migrant detection, mixture analysis, and unsupervised clustering answer
different questions and require different assumptions.

### Reference individuals and leave-one-out estimation

When an individual belongs to the reference data, its genotype can
influence the allele frequencies used to evaluate it. This creates a
favourable bias toward its recorded source. With `leave.one.out = TRUE`,
`assigner` removes the focal individual’s allele contribution from its
home stratum before calculating that likelihood. Unknown individuals,
whose `STRATA` value is missing, do not contribute to any reference
frequency.

Reference population size matters. Small samples produce noisy
allele-frequency estimates and can exaggerate assignment confidence.
GenoDive uses more than five individuals as a minimum eligibility rule
for source populations; this is a useful warning threshold, not a
universal guarantee of adequate power. Paetkau et al. (2004) found that
substantially larger samples could be needed for well-calibrated migrant
tests. Inspect the number of called reference genotypes per stratum and
marker rather than relying only on the total sample count.

### Alleles absent from a reference sample

An allele carried by an individual can be absent from a finite reference
sample even when it exists at low frequency in the population. Using an
observed frequency of zero would force the multilocus likelihood to
zero. The native engine therefore applies a lower allele-frequency
bound. Its default is adaptive, `1 / (2N + 1)`, using the number of
called reference individuals at that marker. To reproduce the fixed
baseline described by Paetkau et al. (2004) and used by GenoDive, use:

``` r

native_fixed_floor <- assigner::assign_individuals(
  data = genome,
  frequency.floor = 0.005
)
```

The floor prevents numerical impossibility; it does not supply evidence
that a candidate source was adequately sampled.

### Unequal source populations and effective sample size

Unequal reference populations have always been an important assignment
concern: larger populations provide more stable allele-frequency
estimates and can gain an artificial advantage over small, poorly
represented populations. One of assigner’s particularly useful
exploratory tools is repeated population equalization:

``` r

balanced <- assigner::evaluate_assignment(
  data = genome,
  assignment.analysis = "gsi_sim",
  markers.sampling = "random",
  subsample = "min",
  iteration.subsample = 20
)
```

This repeatedly draws the same number of reference individuals from
every population. Compare the distribution of results with the full
analysis; do not keep only the most favourable replicate.

For called genotypes, equalizing population counts is a clear starting
point. For GL/PL data, equal counts can still contain unequal
information because coverage and genotype uncertainty differ.
[`assign_individuals()`](https://thierrygosselin.github.io/assigner/reference/assign_individuals.md)
therefore returns `effective.sample.size`, based on observed Fisher
information. DeSaix et al. (2024) show that equalizing this effective
size, rather than sample count alone, can materially improve
low-coverage assignment.

The result also contains `effective.sample.size.markers` and
`individual.effective.sample.size`. The first reveals markers or
populations with little effective information despite nominal sample
size; the second helps identify reference individuals whose low coverage
contributes much less than a fully observed diploid sample. With the
adaptive default, GL/PL frequency bounds use this marker-specific
effective size when it is available.

### Marker selection and high-grading bias

Anderson (2010) demonstrated a common but easily overlooked problem. If
the same individuals are used first to identify apparently informative
loci and then to evaluate the selected panel, predicted assignment
accuracy is systematically too high. He called this high-grading bias.
The problem becomes especially important when many candidate loci are
screened, reference samples are small, or candidate populations are
weakly differentiated.

Ordinary leave-one-out assignment is not sufficient in this situation.
Removing the focal individual while estimating allele frequencies does
not undo its earlier influence on marker ranking. If marker selection
used that individual, information has already crossed from evaluation
into training.

Use one of the following designs:

- **Training and holdout:** select markers and estimate reference
  frequencies using training individuals, then assess the chosen panel
  using untouched holdout individuals.
- **Training-Holdout-Leave-one-out:** use training individuals alone to
  select the marker panel. After marker selection, use all available
  reference individuals to improve allele-frequency estimates, but
  exclude each holdout individual’s own genotype when evaluating it.
- **Independent validation:** evaluate the completed panel on new
  known-origin samples, ideally including a separate sequencing or
  processing batch.

The holdout data must remain untouched during marker ranking, parameter
tuning, threshold selection, and decisions about the form of the
analysis. Repeatedly examining holdout performance and changing the
panel in response turns the holdout into training data. A final
independent dataset is then required.

This warning applies beyond likelihood assignment. FST ranking,
machine-learning feature selection, LD pruning thresholds, missingness
filters, and decisions about the number of loci must all be performed
inside the training portion when their performance is being evaluated.

### Assignment is not a migrant test

[`assign_individuals()`](https://thierrygosselin.github.io/assigner/reference/assign_individuals.md)
ranks the supplied candidate sources. It does not run the Monte Carlo
exclusion procedure described by Cornuet et al. (1999) and Paetkau et
al. (2004), and it does not estimate a migration rate. A best assignment
can still be poor in absolute terms, and the true source can be absent
from the reference data. Migrant detection requires an explicit
statistic, a calibrated null distribution, a stated error threshold, and
a design appropriate to whether all possible sources were sampled.

Paetkau et al. (2004) distinguish the likelihood in the sampled home
population from a likelihood ratio comparing home with the best
alternative. The ratio can have greater power when all sources are
sampled; the home likelihood is safer when possible sources are missing.
Their Dlr measure is a population-pair summary useful for assessing
assignment resolution, not an individual posterior probability.

## Native likelihood assignment

The native assignment workflow provides:

- explicit `"GT"`, `"GL"`, and `"PL"` genotype modes;
- EM allele-frequency estimation for GL/PL reference data;
- integration over all three diploid genotypes instead of forcing
  uncertain low-coverage data into hard calls;
- leave-one-out reference-frequency estimation;
- the complete likelihood of every individual in every candidate
  population;
- Fisher-information effective sample sizes for diagnosing unequal
  reference information; and
- an optional coverage-aware diagnostic for a possible unsampled source.

The native engine estimates allele frequencies in every reference
stratum and calculates the genotype log-likelihood of every individual
under every candidate stratum. It retains the complete likelihood table,
not only the best and second-best assignments.

``` r

native <- assigner::assign_individuals(
  data = genome,
  leave.one.out = TRUE
)

native$assignment
native$likelihoods
```

The native engine can use called genotypes or genotype likelihoods. The
GL/PL implementation follows the population-assignment framework of
DeSaix et al. (2024), accounting for uncertainty in both focal genotypes
and reference-population allele frequencies. Do not substitute the most
likely low-depth genotype and interpret it as certain.

Use explicit field-standard genotype modes:

``` r

gl_assignment <- assigner::assign_individuals(
  data = low_coverage_genome,
  genotype.method = "GL"
)

gl_assignment$effective.sample.size
gl_assignment$effective.sample.size.markers
gl_assignment$individual.effective.sample.size
```

`"GL"` accepts log10-scaled likelihoods or normalized likelihoods in
`GL_HOM_REF`, `GL_HET`, and `GL_HOM_ALT`. `"PL"` accepts Phred-scaled
values in the analogous `PL_` columns. There is intentionally no
automatic mode: the analysis records an explicit choice between `"GT"`,
`"GL"`, and `"PL"`.

When allele depths are available, request the exploratory
unsampled-source diagnostic:

``` r

gl_assignment <- assigner::assign_individuals(
  data = low_coverage_genome,
  genotype.method = "GL",
  unsampled.source = TRUE
)
```

This requires `ALLELE_REF_DEPTH` and `ALLELE_ALT_DEPTH`. A strongly
negative `UNSAMPLED_Z` indicates that the individual’s data are less
compatible with its assigned source than expected from reference
individuals with comparable read depth. The default flag at `-3` is
exploratory, not universal; calibrate it with the reference distribution
and known controls.

If the same individuals can be downsampled to resemble the coverage of
the samples that will ultimately be assigned, keep the full data in
`data` for reference-frequency estimation and supply the lower-coverage
likelihoods only for evaluation:

``` r

coverage_matched <- assigner::assign_individuals(
  data = full_reference_gl,
  evaluation.data = downsampled_reference_gl,
  genotype.method = "GL",
  leave.one.out = TRUE
)
```

This asks a realistic question: how well would lower-coverage
individuals be assigned by the stronger reference baseline? The two
inputs must describe the same individual-marker combinations. Generate
the downsampled GLs independently from reads or likelihood-aware tools;
do not manufacture uncertainty by merely adding noise to called
genotypes.

To validate the exploratory unsampled-source statistic against scores
generated independently with WGSassign:

``` r

comparison <- assigner::compare_unsampled_scores(
  gl_assignment,
  wgsassign_scores,
  id.column = "INDIVIDUALS",
  score.column = "WGSASSIGN_Z"
)

comparison$pearson
comparison$scores
```

The comparison function joins and summarizes existing results; it
neither runs nor reproduces WGSassign.

For reference individuals, leave-one-out estimation avoids using the
focal individual to estimate the allele frequencies against which it is
evaluated. The complete likelihood table can be used directly to
calculate Dlr:

``` r

native_dlr <- assigner::dlr(native)
native_dlr$dlr.table
native_dlr$dlr.matrix
```

### Missing-data patterns and marker blocks

For called genotypes, the full likelihood table contains `Z_LIKELIHOOD`,
an expected-likelihood standardization based on the exact markers
observed for each individual. It avoids imputing genotypes and makes
likelihood compatibility more comparable when missing-data patterns
differ. It does not repair batch effects or make non-random missingness
harmless.

Use biological blocks whenever chromosome or linkage-group information
is available:

``` r

blocked <- assigner::assign_individuals(
  genome,
  marker.blocks = marker_metadata |>
    dplyr::transmute(MARKERS, MARKER_BLOCK = CHROM)
)

blocked$marker.block.likelihoods
```

For an exploratory partition when genomic positions are unavailable, an
integer creates interleaved marker blocks. Compare whether the same
candidate repeatedly wins across blocks. These blocks are sensitivity
diagnostics, not independent replicates and not a substitute for linkage
information.

### Collections and broader reporting units

When sample metadata contains a broader grouping such as `REPUNIT`,
request both levels in the same analysis:

``` r

two_levels <- assigner::assign_individuals(
  genome,
  strata = sample_metadata,
  reporting.unit = "REPUNIT"
)

two_levels$assignment
two_levels$reporting.unit$assignment
```

The collection-level result is preserved. This matters when several
collections are difficult to distinguish individually but jointly form a
biologically useful reporting unit.

### Monte Carlo reference-baseline assessment

Leave-one-out evaluates one known individual at a time. A complementary
Monte Carlo assessment repeatedly removes several reference individuals,
treats them as an artificial mixture of unknown origin, and assigns them
against the remaining baseline:

``` r

mc <- assigner::assess_assignment_mc(
  genome,
  repetitions = 50,
  mixture.size = 100,
  min.remaining = 5,
  random.seed = 20260822
)

mc$replicate.summary
mc$stratum.summary
```

The function retains at least `min.remaining` individuals in every
source and records every holdout identity. If marker selection or
filtering is being evaluated, repeat it using only each replicate’s
retained training baseline; otherwise information can leak from the
artificial mixture into the model.

## Machine-learning assignment

Optional engines return the probability of every individual under every
candidate stratum. Install the R engines with:

``` r

install.packages(c("glmnet", "xgboost"))
```

``` r

elastic <- assigner::assign_individuals(
  genome, engine = "elastic_net", folds = 5, random.seed = 20260822
)

boosted <- assigner::assign_individuals(
  genome, engine = "xgboost", folds = 5, random.seed = 20260822
)
```

TabPFN is an experimental Python engine. In Terminal:

``` bash
conda activate genomics
python -m pip install tabpfn
```

Install `reticulate` in R, start R from the activated environment, and
run:

``` r

install.packages("reticulate")
tabpfn <- assigner::assign_individuals(
  genome, engine = "tabpfn", folds = 5, random.seed = 20260822
)
```

Missingness requires special attention for every engine. Elastic net
uses training-fold mean imputation, XGBoost learns missing-value
routing, and TabPFN uses native preprocessing with a missingness
indicator. Consequently, apparent assignment power can reflect
sequencing batch or processing history. Always inspect
`result$missingness`, keep preprocessing inside validation folds, and
when possible validate on an independent sequencing batch.

Likelihood and machine-learning scores should not be interpreted
identically. The likelihood engine encodes a population-genetic model
and reports genotype log-likelihoods. The machine-learning engines
estimate predictive class probabilities under their fitted models.
Compare methods with honest holdout or nested cross-validation, preserve
every preprocessing step within the training fold, and inspect confusion
among all candidate strata rather than accuracy alone.

## Explore the dataset through more than one lens

Trying the same dataset with several well-established programs should be
a natural part of population-assignment work. It is scientific curiosity
and a good research habit, not an attempt to select whichever result
looks best. Different programs emphasize different questions,
assumptions, priors, validation procedures, and treatments of missing or
uncertain genotypes. Agreement can increase confidence; disagreement is
often more informative because it identifies an assumption, population,
marker panel, or individual that deserves closer inspection.

Useful complementary perspectives include:

- [GenoDive](https://www.patrickmeirmans.com/software/genodive/),
  developed by Patrick G. Meirmans, for classical likelihood assignment,
  exclusion and migrant-testing workflows in an accessible graphical
  environment;
- [`gsi_sim`](https://github.com/eriqande/gsi_sim), developed by Eric C.
  Anderson, for assignment power, marker panels, holdouts and
  simulations;
- [`rubias`](https://github.com/eriqande/rubias), developed by Eric C.
  Anderson and Ben Moran, when mixture proportions and conditional
  genetic stock identification are the main estimands;
- [`WGSassign`](https://github.com/mgdesaix/WGSassign), developed by
  Matthew G. DeSaix and collaborators, for low-coverage whole-genome
  data represented by genotype likelihoods; and
- [`adegenet`](https://cran.r-project.org/package=adegenet), developed
  by Thibaut Jombart and collaborators, for multivariate and predictive
  analyses such as DAPC.

These programs are not interchangeable, so compare like with like: use
the same reference individuals, candidate sources, marker set and
validation split, and record every frequency correction, prior,
missing-data rule and random seed. Do not treat a majority vote among
programs as proof. Instead, ask why a result changes and whether the
difference is biological, statistical or technical. Whenever possible,
finish with independent samples or a sequencing batch that was not used
for marker selection, model fitting or parameter tuning.

## How assigner differs from rubias

[`rubias`](https://github.com/eriqande/rubias) was developed by Eric C.
Anderson and Ben Moran and implements Bayesian inference for the
conditional genetic stock identification model. The methods are related,
but their primary questions differ:

- [`assigner::assign_individuals()`](https://thierrygosselin.github.io/assigner/reference/assign_individuals.md)
  asks which supplied reference population best explains each
  individual’s genetic data and retains all source likelihoods;
- [`assigner::evaluate_assignment()`](https://thierrygosselin.github.io/assigner/reference/evaluate_assignment.md)
  asks how assignment performance changes across marker panels,
  holdouts, and population subsamples; and
- `rubias` is particularly designed to estimate the proportions
  contributed by reporting units or populations to a mixture, jointly
  with individual posterior assignments.

Because rubias performs a mixture analysis, inferred mixture proportions
can inform individual posterior membership.
[`assign_individuals()`](https://thierrygosselin.github.io/assigner/reference/assign_individuals.md)
does not estimate mixture proportions and evaluates individuals directly
against the reference strata. Neither approach is universally
preferable: use the method that matches the biological estimand.

## Training-Holdout-Leave-one-out with gsi_sim

This workflow uses [`gsi_sim`](https://github.com/eriqande/gsi_sim),
developed by assigner co-author Eric C. Anderson. Eric also developed
assigner’s original `gsi_sim` wrapper and contributed to the package’s
early development.

Holdout samples do not contribute to marker ranking or training, which
avoids evaluating assignment on the same samples used to select markers.

``` r

result <- assigner::evaluate_assignment(
  data = genome,
  assignment.analysis = "gsi_sim",
  markers.sampling = "ranked",
  marker.number = c(100, 250, "all"),
  thl = 0.20,
  iteration.method = 5,
  random.seed = 20260822
)
```

The dated output folder contains the argument log, holdout selections,
raw results, summaries, and plots. Retain it with the analysis.

## R-based assignment with adegenet

``` r

result_adegenet <- assigner::evaluate_assignment(
  data = genome,
  assignment.analysis = "adegenet",
  markers.sampling = "random",
  marker.number = c(100, "all"),
  iteration.method = 5,
  random.seed = 20260822
)
```

## Estimate FST

``` r

fst <- assigner::fst_WC84(
  data = genome,
  pairwise = TRUE,
  ci = TRUE,
  iteration.ci = 100,
  verbose = TRUE
)
```

See the [FST
vignette](https://thierrygosselin.github.io/assigner/articles/fst_confidence_intervals.md)
for interpretation and additional examples.

## References

- Cornuet J-M, Piry S, Luikart G, Estoup A, Solignac M (1999). New
  methods employing multilocus genotypes to select or exclude
  populations as origins of individuals. *Genetics* 153: 1989-2000.
- Anderson EC (2010). Assessing the power of informative subsets of loci
  for population assignment: standard methods are upwardly biased.
  *Molecular Ecology Resources* 10(4): 701-710.
  <doi:10.1111/j.1755-0998.2010.02846.x>.
- DeSaix MG, Rodriguez MD, Ruegg KC, Anderson EC (2024). Population
  assignment from genotype likelihoods for low-coverage whole-genome
  sequencing data. *Methods in Ecology and Evolution* 15: 493-510.
  <doi:10.1111/2041-210X.14286>.
- Manel S, Gaggiotti OE, Waples RS (2005). Assignment methods: matching
  biological questions with appropriate techniques. *Trends in Ecology &
  Evolution* 20(3): 136-142. <doi:10.1016/j.tree.2004.12.004>.
- Moran BM, Anderson EC (2019). Bayesian inference from the conditional
  genetic stock identification model. *Canadian Journal of Fisheries and
  Aquatic Sciences* 76(4): 551-560. <doi:10.1139/cjfas-2018-0016>.
- Paetkau D, Calvert W, Stirling I, Strobeck C (1995). Microsatellite
  analysis of population structure in Canadian polar bears. *Molecular
  Ecology* 4: 347-354.
- Paetkau D, Slade R, Burden M, Estoup A (2004). Genetic assignment
  methods for the direct, real-time estimation of migration rate: a
  simulation-based exploration of accuracy and power. *Molecular
  Ecology* 13: 55-65. <doi:10.1046/j.1365-294X.2003.02008.x>.
