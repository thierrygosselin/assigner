# Controlled population-structure benchmark

## Purpose and scope

This benchmark separates known population structure from genotype uncertainty and
missingness. It is an experimental validation of the GL mixture, not evidence
that genetic clusters are demographically independent management stocks.

The generator is internal to grur (`R/simulation_benchmark.R`). It uses independent,
biallelic, diploid Hardy–Weinberg loci. An ancestral frequency q is drawn uniformly
between 0.05 and 0.95. Each population frequency is drawn independently from
Beta(q(1-F)/F, (1-q)(1-F)/F); F = 0 uses q directly. This is a controlled
Balding–Nichols construction, not a migration or demographic simulation.
The generating drift parameter is not asserted to equal realized sample FST.

Complete-truth and explicitly masked-truth global and pairwise WC84 estimates
are calculated using hierfstat. An equal-population frequency variance ratio
is recorded separately as TRUE_FREQUENCY_GST; it is not labelled WC84.
MASKED_TRUTH_WC84 does not measure low-depth genotype-calling error.

## Dataset design

Two independent replicates are generated for each of 12 scenarios:

- Three drift levels: 0.001, 0.01, 0.1; 30 individuals per group and 100 markers.
- Sample size: 10 or 100 individuals per group.
- Marker count: 500 markers.
- Unequal samples: 60/34/6.
- Independent random missingness: 30%.
- Group-associated missingness: 5%/25%/50%.
- Poisson sequencing depth with mean 3 and allele observation error 0.01.
- One-population null and strongly differentiated two-population control.

Unless stated otherwise, scenarios use three groups of 30, 100 markers and
drift 0.01. Scenarios use different truth seeds: comparisons are not paired
subsamples of the same underlying population realization. Two replicates are
a smoke-test scale, not enough to estimate power or false-positive rates.

The missingness and low-depth scenarios use drift 0.01, for which clean-data
recovery was already poor in the initial pilot. They cannot isolate a loss of
recovery caused by missingness. A follow-up must use paired observations of the
same strongly structured truth, alongside paired null controls, before claiming
missingness robustness or estimating its effect.

## Observation and truth separation

Truth genotypes, frequencies and STRATA remain in separate RDS files. GL-only
GDS files contain the observed likelihoods, not truth group metadata. Explicitly
masked cells and zero-depth cells have equal likelihoods for all genotypes.
Finite-depth likelihoods come from binomial alternate-read counts conditional
on genotype and sequencing error. Complete-data GLs are near-certain synthetic
likelihoods. This model omits mapping bias, allelic imbalance and relatedness.
Synthetic chromosome/position fields are indexing placeholders, not linkage.

## Fitting protocol

Run from the assigner repository:

```sh
Rscript inst/experimental/run-grur-structure-benchmark.R /path/to/new/output
```

The destination must not exist. No inputs are overwritten. The default generates
all 24 datasets but fits only replicate 1 of the three clean differentiation
levels and the K = 1/2 controls. `--extended` fits every generated dataset.

To evaluate the remaining datasets without regenerating inputs or overwriting
the pilot, use the continuation script from the same repository:

```sh
Rscript inst/experimental/continue-grur-structure-benchmark.R ORIGINAL_OUTPUT NEW_RESULTS 4
Rscript inst/experimental/summarize-grur-structure-benchmark.R ORIGINAL_OUTPUT NEW_RESULTS NEW_SUMMARY
```

Continuation uses the same settings and manifest-indexed seeds as the full
runner. Up to four independent processes each open their own GDS handles.
The original five results are reused, not re-estimated. Each additional dataset
has its own progress log, fitted objects and result checkpoint. Errors remain
explicit in the combined results. Source checksums, session information and
input GDS checksums are recorded; GDS checksums are verified after completion.
The summary retains fold warnings/errors, individual-start diagnostics and
membership–missingness audits, plus a plot of within-dataset predictive scores.

Each selected dataset receives a known-K diagnostic fit and two-fold predictive
cross-validation across K = 1:5. Fits use seven mixed starting points, pooled
regularization of strength 2 and a maximum of 300 iterations. Truth labels are
used only for evaluation, not initialization or likelihood fitting. The known-K
fit deliberately uses the generating K; the cross-validation candidates do not.
Alternating folds distribute each generated group across the training and test
sets. This is not a leave-population-out prediction exercise.

Adjusted Rand index (ARI) evaluates label-invariant partition recovery at known
K, excluding unassigned samples whose count is separately reported. Predictive
K is a distinct diagnostic. Convergence is not proof of a global optimum.
The chosen prior and limited starts are experimental settings, not calibrated
defaults or a guarantee of recovering K.

## Files and reproducibility

The initial run is stored in:

`/Users/gos040/.codex/.chatgpt-projects/g-p-6a8255ffbc948191b3592f138fe06cf5/artifacts/grur-structure-20260916`

This contains scenario definitions, seeds, session information, manifest,
truth and observation audit RDS files, GL GDS files, fitted objects and CV
objects. The final `pilot-results.txt` reports only completed fitted scenarios.
Generated-but-unfitted datasets must not be presented as validation results.

## Initial pilot results (16 September 2026)

All five selected known-K fits converged and assigned all 90 samples. All
cross-validation summaries were complete.

| Generating groups | Measured truth WC84 | Known-K ARI | Predictive K |
| --- | --- | --- | --- |
| 3, drift 0.001 | 0.00132 | 0.013 | 1 |
| 3, drift 0.01 | 0.01673 | 0.074 | 1 |
| 3, drift 0.1 | 0.09464 | 0.934 | 3 |
| 1, null | Not applicable | 1.000 | 1 |
| 2, drift 0.1 | 0.12086 | 1.000 | 2 |

These results show recovery of strong structure and no extra groups in this
single null replicate. Weak structure was not recovered with 100 markers under
these settings. This does not distinguish insufficient information from
optimization or regularization limitations; more markers, repeated simulations
and initialization sensitivity are needed. The other 19 datasets were not
fitted in the initial run; their completed continuation is reported below.

## Completed 24-dataset benchmark (16 September 2026)

All 24 known-K fits and all 240 selected cross-validation fold fits converged.
All K comparisons were complete, with zero fold warnings and zero fold errors.
The 168 starts retained for the known-K fits also all converged; this does not
mean that they agreed or reached global optima. Input GDS checksums were unchanged.
The five focused mixture regression-test files passed without failures.

| Scenario | Sample counts | Markers | Selected K, replicates 1 / 2 | Known-K ARI, replicates 1 / 2 |
| --- | --- | --- | --- | --- |
| Drift 0.001 | 30/30/30 | 100 | 1 / 1 | 0.013 / 0.003 |
| Drift 0.01 | 30/30/30 | 100 | 1 / 1 | 0.074 / 0.006 |
| Drift 0.1 | 30/30/30 | 100 | 3 / 3 | 0.934 / 1.000 |
| Small samples | 10/10/10 | 100 | 1 / 1 | -0.045 / 0.148 |
| Large samples | 100/100/100 | 100 | 1 / 1 | 0.071 / 0.045 |
| More markers | 30/30/30 | 500 | 1 / 1 | 0.149 / 0.123 |
| Unequal samples | 60/34/6 | 100 | 1 / 1 | 0.114 / -0.001 |
| Random missingness | 30/30/30 | 100 | 1 / 1 | 0.021 / 0.017 |
| Structured missingness | 30/30/30 | 100 | 1 / 1 | 0.020 / 0.004 |
| Mean depth 3 | 30/30/30 | 100 | 1 / 1 | -0.003 / 0.091 |
| One-population null | 90 | 100 | 1 / 1 | 1.000 / 1.000 |
| Strong two-population control | 45/45 | 100 | 2 / 2 | 1.000 / 1.000 |

Rows without an explicit drift level use 0.01, except the one-population null
and two-population control. The null's known-K ARI of 1 is trivial at K = 1;
its useful check is that cross-validation did not select additional components.

Strong structure was recovered in both replicates. The weaker settings were
not resolved under this design, including increasing samples OR markers. This
is not a universal detection threshold, proof that structure is absent, or a
demonstration that optimization and regularization are optimal. Missingness and
depth results require the paired controls described above for causal comparisons.
Two replicates per scenario do not support precise power or false-positive rates.

Additional fits and raw results:

`/Users/gos040/.codex/.chatgpt-projects/g-p-6a8255ffbc948191b3592f138fe06cf5/artifacts/grur-structure-extended-20260916`

Combined tables, per-fold/start audits, conditional-membership audit and
`predictive-K.png`:

`/Users/gos040/.codex/.chatgpt-projects/g-p-6a8255ffbc948191b3592f138fe06cf5/artifacts/grur-structure-summary-20260916`

The two 300-individual datasets took approximately 24.7 and 20.2 minutes each
for their known-K fit plus full K comparison while other jobs also ran.
These are batch timings, not isolated hardware benchmarks.

## Conditional membership is not calibrated confidence

The continuation's `more_markers_rep1` known-K fit has mean maximum membership
approximately 0.9999995, but its ARI against truth is only 0.149. These are
in-sample memberships conditional on a fitted three-component model and its
estimated parameters. They do not integrate parameter or K uncertainty, and
high values do not establish recovery of the generating populations.

The report retains this conditional-membership audit alongside ARI and held-out
K scores. Do not present a maximum membership of 0.95 as demonstrated 95%
assignment accuracy. Calibration requires independent labelled test samples and
replicated experiments; estimating training membership alone cannot provide it.
The larger-sample runs also reach different objectives from different starts,
so initialization remains a separate validation concern.

## Demographic simulation branch

grur::simulate_rad() now has sample.size independent of ne. A scalar is shared
across populations; a vector supplies per-population sample sizes for a single
num.pops setting. NULL retains legacy sampling. Forward simulations initialize
at ne, then sample surviving individuals without replacement; insufficient
survivors cause an explicit error. Tests cover settings and subsampling.

The current R environment lacks rmetasim and no fastsimcoal executable was
located. Consequently this run does not validate real fastsimcoal/rmetasim
execution. The controlled GL generator does not require these backends.

## Next validation stages

The [truth-initialization diagnostic](truth-initialization-diagnostic.md) is now
complete. It separates blind initialization, truth-started estimation and
fixed-known-frequency references, and documents why improved starts alone do
not resolve the weak-differentiation predictive-K problem.

Performance also needs attention before a large factorial grid. The 300-sample
runs are much slower than the smaller cases. A read-only diagnostic compared
the existing `.sm_marker()` row maximum (`apply(..., 1, max)`) with a three-column
`pmax()` calculation. For 1,000 calls on a synthetic 300-by-3 likelihood matrix,
elapsed times were 0.161 versus 0.045 seconds and outputs matched `all.equal()`.
This is a small kernel timing, not an end-to-end speedup claim or validation of
all edge cases. No mixture-engine code was changed during the benchmark.
Profile full fits and add edge-case numerical-equivalence tests before adopting
this or other optimizations.

The staged datasets are now fitted. Next, test independent-sample assignment
calibration and paired strong-structure/null missingness controls. Increase
independent replicates before estimating error rates. Use paired marker panels
to isolate marker-count effects. Add rare/private alleles, linkage, family structure,
asymmetric demographics and technical biases before making applied claims.
The ancestral-frequency bounds exclude very rare ancestral variants here.
Scott Foster's final methodological and code review remains necessary.
