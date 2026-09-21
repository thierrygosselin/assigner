# GL/PL stock mixture: validation and remaining limits

Larger-marker comparison: [stockR benchmark results](stockr-marker-results.md),
with the [paired design](stockr-marker-benchmark.md). Known-K recovery and
unknown-K selection are reported separately.

These are internal development functions. Scott Foster's review remains a
requirement before public release. Numerical success does not establish
calibrated uncertainty, a true stock count, or demographic independence.

The separate [controlled grur benchmark](grur-structure-benchmark.md) records
known-truth simulations varying differentiation, sample sizes, marker count,
missingness and sequencing depth. It also documents why conditional membership
probabilities must not be interpreted as calibrated assignment confidence.
The subsequent [truth-initialization diagnostic](truth-initialization-diagnostic.md)
tests whether privileged starting frequencies resolve the weak-structure results.

## stockR comparison

Run `Rscript inst/experimental/validate-stockr-example.R` from assigner's source
directory. Add `--extended` for the small bootstrap and cross-validation run.
The script uses stockR 1.0.76's internal `myData` example (100 individuals,
5000 SNPs, three generating groups), not an empirical tuna dataset. It does
not modify stockR, its example data, or the installed assigner package.

Known genotypes are encoded as PL 0/999/999 (permuted appropriately); missing
genotypes remain missing. This is a near-certain-genotype integration test,
not a reconstruction of read-level uncertainty. Both algorithms use the same
initial grouping-derived frequencies, and the comparison accounts for
stockR's frequency floor of 1/(10*nFish).

The full-data run converged in two iterations from that matched initialization.
Membership probabilities and supported allele frequencies matched stockR
to the reported numerical precision. Both partitions matched the generating
groups. Three stock-marker combinations had zero effective support; their
frequencies differed (maximum 0.577486) without changing the likelihood.
Those cells must not be treated as estimated stock-specific frequencies.
`frequency.support` and `frequency.estimable` now expose this distinction.
These comparisons do not show that arbitrary starts find the global optimum.

## Sampling weights

The engine now maximizes the weighted observed-data objective

$$
\ell_w = \sum_i w_i \log\left[\sum_k \pi_k
   \prod_l \sum_g L_{ilg}\Pr(g\mid p_{kl})\right].
$$

Weights multiply each individual's marginal log likelihood. They do NOT
exponentiate stock membership odds. Zero-weight individuals have no role in
frequency, mixture-proportion or random-start estimation, but are evaluated
under the fitted parameters. This supports holdouts and counted bootstrap
resamples without duplicating GDS samples. A marker without informative
positive-weight training observations stops that refit, rather than borrowing
information from held-out samples or silently changing the marker set.

## Uncertainty at a specified K

```r
uncertainty <- assigner:::.sm_bootstrap_gds(
  data = gds, fit = fit, B = 200, seed = 42,
  n.starts = 5, frequency.variants = c("1", "2")
)
uncertainty$diagnostics
uncertainty$intervals
```

- Whole individuals are resampled with replacement, represented as integer
  sampling weights. Their complete genotype-likelihood profiles stay together.
- This is an ordinary nonparametric individual bootstrap. It is NOT Scott's
  stockBOOT Bayesian bootstrap; its weighting convention is explicitly defined
  above and must be reviewed on its own merits.
- Each replicate refits the mixture, then aligns component labels to the
  reference using a minimum-cost one-to-one frequency match (exact assignment,
  currently K <= 12). Inspect alignment distances: matching labels does not
  establish that a replicate recovered the same biological components.
- Percentile intervals summarize membership probabilities, mixture proportions
  and explicitly selected variant frequencies. They are conditional on K,
  the observed markers, model assumptions and the sampling design.
  Mixture proportions describe this sample; they are not stock abundance
  fractions unless the sampling design justifies that interpretation.
- Missing/failed/nonconverged replicates are retained in diagnostics. If any
  replicate is unusable, intervals are withheld rather than silently computed
  from successful replicates only. Requested frequencies without adequate
  numerical support also make intervals unavailable.
- Individuals must be independent sampling units. Do not apply this routine
  unchanged to family clusters or repeated observations. It does not currently
  support grouped/stratified survey bootstrap weights.
- Linked markers stay together within each resampled individual, but the model
  still treats markers as independent. This bootstrap does not guarantee
  calibration of an overconfident composite likelihood. Marker-block sensitivity
  and replicated coverage experiments remain required.
- Private/fixed alleles and unsupported combinations require care: resampling
  cannot invent unobserved alleles. Narrow boundary intervals are not proof of
  certainty. The epsilon floor is numerical, not a biological prior.
- Membership draws use O(B*n*K) storage; frequency draws are stored only for
  requested variants. max.values limits draw allocation. Whole fit objects
  are not retained for all replicates.

The extended stockR run uses only eight replicates on the first 300 markers
to check execution. This is not an adequate interval-coverage experiment or
a recommended bootstrap replicate count for inference.

All eight initial replicate fits converged. The 300-marker reference fit also
recovered the generating partition; this does not ensure that predictive K
selection will recover the same grouping (see below).

## Comparing K without in-sample likelihood selection

```r
# Supply named fold IDs; keep relatives/repeated observations in the same fold.
comparison <- assigner:::.sm_cv_gds(
  data = gds, K = 1:4, folds = fold.ids,
  n.starts = 5, seed = 42, likelihood = "PL"
)
comparison$summary
comparison$diagnostics
comparison$heldout
```

Every fold fits frequencies and proportions using training individuals only.
Starting frequencies are data-independent (0.5), with additional randomized
training-only starts required when K > 1. No reference full-data solution,
known population labels, or STRATA constraints initialize the folds.
The score for each held-out individual is its marginal mixture log likelihood;
its group is integrated out, not chosen to maximize its test score.

K=1 is mandatory. All K values use the same fold assignments and active
marker set. Reported scores omit data-only likelihood normalization constants,
which cancel between K values on the same observations. Do not compare scores
between different marker sets, likelihood encodings with different scaling,
or unrelated datasets. Fully uninformative test individuals are excluded and
the counts are reported. Marker selection/QC performed before this routine is
not cross-validated; any data-driven selection must be accounted for separately.

The returned best.predictive.K maximizes mean held-out log likelihood across
the supplied grid, with the smaller K winning an exact tie. Any failed fold/K
withholds the overall choice. Inspect all warnings, likelihood differences,
fold composition and best.at.upper.boundary. A boundary optimum means the grid
has not bracketed an optimum. Repeat defensible fold allocations; one partition
and a tiny score improvement are not strong evidence for a particular K.

This is predictive model comparison, not a validated automatic biological
stock-count estimator. Linked regions, batch effects, relatives, poor calling
models and admixture can all influence it. Results can be sensitive to the
frequency bound and sparse/private alleles. There is no formal calibrated
significance threshold or confidence set for K in this prototype.

### A deliberately retained failed stock-count validation

On the first 300 markers of the bundled stockR example, with three folds,
three starts and epsilon = 0.00000001, all nine fold/K fits converged but
held-out scoring preferred K=1 rather than the generating K=3:

| K | Approximate total held-out log likelihood |
|---|---:|
| 1 | -13891 |
| 2 | -14165 |
| 3 | -14959 |

This is a failure to recover the generating stock count, not proof that K=1
is biologically correct. The result is retained rather than tuning it away.
Sparse stock-specific allele-frequency estimates, boundary regularization and
optimization require investigation. The script's optional --sensitivity flag
repeats the same folds and seeds at epsilon 0.001 and 0.01; these are diagnostic
comparisons, not a post-hoc rule for choosing the bound that recovers truth.

The completed sensitivity runs also converged in every fold. At epsilon 0.001,
the approximate total scores for K=1/2/3 were -13536/-13030/-13290; at 0.01,
they were -13481/-12827/-12979. Both preferred K=2, still not the generating
K=3. This establishes sensitivity to boundary treatment, not its sole cause.
The default K=3 folds had minimum effective component sizes around 17–22,
whereas the generating groups are highly unbalanced; local optima and
initialization deserve particular investigation. A converged fit is not a
guarantee of the correct partition or global likelihood maximum. K selection
is therefore NOT ready for inferential recommendations.

## Current verification scope

104 focused expectations pass across the four stock-mixture test files, with
no unexpected warnings or failures. They include weighted holdout exclusion,
exact label permutation, reproducible bootstrap intervals, preservation of
GDS filters/RNG, and withholding conclusions for incomplete refits. Two small
synthetic fixtures select the expected K=1 and K=2, respectively. Those easy
fixtures do not override the failed K=3 stockR example above.

## What remains to validate

- Repeated simulations measuring interval coverage and false/missed structure,
  rather than merely recovering one easy generating partition.
- Low coverage, structured missingness, unequal populations, rare/private
  alleles, relatedness, linkage/inversions, batch effects and admixture.
- Alternative starts and label ambiguity, plus frequency-bound sensitivity.
- Realistic runtime and memory across bootstrap replicates and K grids.
- Scott's review of the likelihood, weighting, uncertainty and K interpretation.
# Follow-up: separating optimization from prediction (2026-09-16)

Reproduce from the source root with
`Rscript inst/experimental/investigate-mixture-k.R`.
This diagnostic uses 300 markers from stockR 1.0.76's internal simulated data,
100 individuals, three deterministic held-out folds, and the same blockwise
EM engine. Near-certain likelihoods encode the supplied genotypes; this is not
a low-depth simulation. Generating group sizes are **60, 6 and 34**.

Three initialization strategies were examined separately:

- Random: symmetric supplied start plus nine random soft starts; select the
  best converged training likelihood.
- PCA: temporary training-only mean completion, five training PCs, k-means
  with 25 starts, then one likelihood EM fit. No held-out values initialize it.
- Oracle: known **training** generating labels initialize frequencies, then
  unconstrained EM. This is a diagnostic, not an available unsupervised method.

All 63 fits converged. For K = 3, oracle initialization recovered the generating
training partition in all nine fold/bound combinations. Neither random nor PCA
initialization recovered it exactly in any of those nine combinations. Oracle
training likelihood exceeded both alternatives in every combination checked.
These results establish an optimization problem, but not that the oracle is
the global maximum. PCA initialization alone is not a demonstrated solution.

Summed held-out log likelihoods (larger is better):

| Frequency bound | K=1 | K=2, random | K=3, random | K=3, PCA | K=3, oracle |
|---|---:|---:|---:|---:|---:|
| 0.00000001 | -13891.26 | -14289.86 | -15004.18 | -14838.06 | -14486.00 |
| 0.001 | -13536.17 | -13049.52 | -13342.83 | -13228.06 | -13112.74 |
| 0.01 | -13481.19 | -12833.40 | -12996.75 | -12950.07 | -12879.40 |

Even the correct-partition oracle does not beat the competing smaller-K score.
Thus the missed generating K cannot be attributed solely to local optima.
The smallest group has four training and two held-out individuals in each fold.
Among its 300 markers, respectively **26, 33 and 26** contain a held-out allele
absent from that group's training observations. Boundary frequency estimates
penalize these observations, explaining substantial sensitivity to epsilon.
This identifies a concrete mechanism, not proof that it explains every score
difference. Finite-sample predictive K need not equal generating K.

## Consequences for the implementation

No frequency default or core estimator was changed to force K = 3. Production
starts remain random soft starts; PCA above is a diagnostic comparator only.
The held-out helper remains experimental and must not be advertised as validated
stock-count estimation. The next methodological decision is a coherent
frequency-regularization/predictive model (rather than arbitrary clipping),
reviewed with Scott, plus initialization improvements assessed on independent
replicate simulations. Hyperparameters must not be selected using the final
held-out evaluation or knowledge of its generating labels.

Other limitations needing validation remain: weak differentiation, unequal
sampling, GL/PL calibration at low depth, population-structured missingness,
private alleles, related individuals, linkage, admixture, bootstrap coverage,
and ambiguous label alignment. The three unsupported stock-marker cells in the
full example remain explicitly flagged; positive support is not a guarantee
of an accurate estimate. The current 104 focused checks pass, but are software
checks, not substitutes for these statistical validation experiments.

## Implemented extensions and pilot results (2026-09-16)

The preceding results are retained as the historical baseline. New opt-in
controls provide mixed initialization and coherent Beta MAP regularization;
defaults remain unchanged. This is still internal research code awaiting Scott's
review, not a validated automatic stock-count procedure.

`start.strategy = "mixed"` cycles through 2/5/10-PC k-means starts, random soft
starts, lower-K component splits and perturbations of the best frequencies.
The PCA matrix contains only positive-weight informative training individuals,
uses bounded marker subsets and temporary pooled completion, and never replaces
the original likelihoods. Lower-K splits cycle through components across starts.
All starts are fully optimized; there is no short-run screening or annealing.
Starts that fail are retained in the fitting diagnostics.

`regularization = "symmetric"` uses q=0.5, while `"pooled"` estimates q from
training likelihoods. For strength s the prior is Beta(1+s*q, 1+s*(1-q)).
EM maximizes data log likelihood plus the log-prior penalty; it does not assume
that unpenalized likelihood must increase. Held-out scores exclude the penalty.
Pooled q is fixed within a fit and re-estimated within each training fold or
bootstrap replicate. This is not integrated Bayesian prediction or a fully
estimated hierarchical model. Strength must not be tuned on the final test set.

### Original stockR example, 300 markers

Thirteen mixed starts, three folds, and a pre-specified pooled strength of 2:

| Model | K=1 held-out sum | K=2 held-out sum | K=3 held-out sum | Winner |
|---|---:|---:|---:|---:|
| No regularization | -13891.26 | -14319.83 | -14841.44 | 1 |
| Pooled Beta MAP, s=2 | -13891.26 | -13339.11 | -13459.73 | 2 |

All 18 selected fits converged. Mixed starts improve K=3 training likelihood
over the earlier ten-random-start run, but still do not recover the exact known
training partition. Regularization improves held-out prediction but does not
make K=3 win. These results explicitly remain a failed generating-K benchmark.

### Independent simulation ladder: 60 markers

Run `Rscript inst/experimental/validate-mixture-ladder.R`.
Seven scenarios each have two independent simulation seeds, 60 individuals,
two simulation-stratified folds, K=1:4, seven mixed starts, and strengths 0 and 2.
The regularized comparisons use pooled centers. The generating labels define
the simulation folds and evaluation only, never the fitted starts or parameters.
All 224 selected fits converged. Three runs reported failed individual starts;
these warnings are retained, and alternative starts supplied converged fits.

| Scenario | Generating K | Winner without regularization (two seeds) | Winner with pooled s=2 |
|---|---:|---|---|
| One-population control | 1 | 1, 1 | 1, 1 |
| Two-population control | 2 | 2, 2 | 2, 2 |
| Balanced groups (20/20/20) | 3 | 2, 3 | 3, 3 |
| Unequal groups (36/4/20) | 3 | 2, 2 | 3, 3 |
| Weak differentiation | 3 | 1, 1 | 1, 1 |
| Low depth (mean 2 reads; correct read-error model) | 3 | 3, 3 | 3, 3 |
| Group-structured missingness (5%, 25%, 55%) | 3 | 1, 1 | 3, 3 |

At the generating K, both estimators recover the exact training partition in
all four fold/seed combinations of each strongly differentiated scenario,
including the unbalanced and low-depth scenarios. Neither does so in the weak
scenario. This separates optimization/partition recovery from predictive K.
Two seeds are far too few to estimate reliability, power or error rates.

### Larger-marker comparison: 180 markers

`Rscript inst/experimental/validate-mixture-ladder.R --markers=180` repeats the
same seven scenario designs with two seeds and 180 simulated markers. All 224
selected fits converged; two runs reported failed individual starts. These are
fresh simulated marker panels, not nested subsets of the 60-marker data.

The regularized winners for both seeds remain 1 for the one-population control,
2 for the two-population control, 3 for balanced, unbalanced, low-depth and
structured-missingness scenarios, and 1 for weak differentiation. Without
regularization, winners remain 2/3 for balanced groups and 2/2 for unequal groups;
structured missingness gives 1/2. Adding markers does not by itself resolve
small-group frequency prediction. Across both marker counts this is 448 selected
fits, but still only two seeds per scenario/count: not a power or coverage study.

### Missingness-only negative control

`validate-mixture-null-missingness.R` simulates one population, 60 individuals,
60 markers and three technical groups with missingness probabilities 5%, 30%
and 65%. With mixed starts and pooled s=2, the held-out sums for K=1:4 are
-2444.756, -2510.845, -2549.568 and -2581.929: K=1 wins. This single pilot does
not prove immunity to technical clustering, mapping bias or allele dropout.

### Scope of validation

The uncertain-genotype MAP update also agrees with independent one-dimensional
numerical optimization in a unit test. Other checks cover RNG preservation,
held-out exclusion, GDS conversion/filter preservation, objective monotonicity,
memory guards and propagation through bootstrap/CV. Private-allele loss under
shrinkage, LD/relatedness, misspecified GLs, real low-coverage datasets, Bayesian
frequency integration and interval coverage still need dedicated studies.

The focused software suite now has **141 passing expectations**, with zero
failures, warnings or skips. The vignette renders and all three internal help
files parse successfully. These checks do not replace the validation limits
listed above. No package installation, export, commit or push was performed.
