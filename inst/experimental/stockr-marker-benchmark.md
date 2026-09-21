# Paired stockR / assigner marker-count benchmark

## Purpose

Compare group recovery at known K separately from selection of K. Simulation
labels supply evaluation truth, never initialization or model-selection criteria.
This is a controlled clean-data benchmark, not validation of management stocks.

## Locked design

- Marker panels: 1,000, 3,000, 5,000, 10,000, nested within each realization.
- Training and separate evaluation cohorts: 300 individuals each, balanced groups.
- Three-population drift parameters: 0.001, 0.005, 0.01; a two-population control
  at 0.01; and a one-population null. Drift parameters are not measured sample FST.
- Two independent population-frequency realizations per condition: 40 panels.
  The four nested panels are correlated, not four independent replicates.
- Complete diploid biallelic independent loci, ancestral frequencies uniform
  between 0.05 and 0.95, HWE within each generating group.
- stockR receives genotype counts; the GL engine receives equivalent GDS GL
  vectors (0 for the observed genotype, -80 log10 for alternatives).
- Seven identical blind starting partitions: PCA/k-means with 2, 5, 10 PCs,
  then four random partitions. PCA uses at most 1,000 evenly spaced markers and
  training samples only; k-means uses 20 restarts. K = 1 needs one analytic/EM start.
- Initial soft membership follows stockR's default softness 2/(K+1), and both
  engines start from the corresponding allele frequencies and proportions.
- stockR 1.0.76, ordinary EM, maximum 300 iterations, parameter tolerance 1e-8.
  GL engine uses objective tolerance 1e-8 and frequency bound 1e-8. These stopping
  rules and bounds differ; equal numeric tolerances do not make them equivalent.
- stockR fits are selected by converged finite training likelihood; GL fits by
  converged objective within a regularization setting. Failed starts remain visible.

## Stages

### Known K

Supply generating K to both methods. Compare stockR with unpenalized GL and
pooled-regularized GL at strength 100. This isolates the unpenalized implementation
comparison and the effect of shrinkage. Strength 100 is a predeclared candidate,
not an established optimum. It is NOT selected using test performance.

Measure independent-cohort ARI, assigned counts, mean maximum membership,
predictive likelihood, fit time and failed starts. ARI is not percentage correct;
conditional membership is not proven calibrated. Fit times exclude shared PCA
initialization and data-generation time, but include method-specific fitting work.

### Unknown K

A separate `cv` mode evaluates K = 1:4 using five-fold predictive validation.
stockR retains its unpenalized model. The GL method jointly selects K and pooled
strength from 2, 20, 100. Choices are locked before final test projection. This
is a common external predictive criterion, NOT a native stockR model-selection
routine. stockR's returned AIC/BIC are retained in fit objects, but are not
presented as equivalent to predictive CV or used to select these models.

One deterministic, balanced fold allocation per dataset; no repeated folds yet.
An incomplete grid withholds selection. Exact/near score ties within 1e-8 favor
smaller K, then weaker regularization. A selected K = 4 would be a search-boundary
result, not confirmation that four groups are present.

## Important fairness and implementation qualifications

- Training-invariant markers are removed for BOTH methods, because stockR removes
  them internally. Held-out data never determine that mask. Record retained counts.
- stockR's likelihood frequency floor depends on training sample count and differs
  for K = 1. Prediction uses the training bound, not the evaluation sample count.
- The shared-start protocol is not stockR's complete default workflow and does
  not evaluate its simulated/deterministic annealing alternatives. It is also
  different from assigner's earlier mixed split/perturb starts. Comparisons with
  earlier experiments therefore are not purely marker-count contrasts.
- In this hard-call limit only, the GL wrapper uses pooled genotype means for
  the prior center and shared start frequencies. Uncertain-GL experiments must
  estimate those from likelihoods; this wrapper must not be used for them.
- stockR returns no explicit convergence flag. Its iteration trace is used to
  reject iteration-cap exits; the returned fits and traces remain available.
- The independent cohort shares the generating frequencies. Two realizations
  are exploratory replication, not a precise estimate of statistical power.

## GDS scaling and tests

`benchmark-streaming-engine.R` reads bounded GL blocks and retains information
counts rather than a full sample-by-marker mask. Vectorized normalization and EM
match the existing calculation in regression checks. The block guard defaults
to 100,000 sample-marker cells, with 250 markers per block at N = 300 (75,000 cells).
This is NOT a total RAM bound: frequency/membership matrices, initialization and
the benchmark's paired hard-call matrices also consume memory. No per-method
peak-RAM comparison is claimed from block sizes.

Tests cover GL/PL/log normalization, flat and missing vectors, malformed inputs,
different block sizes, weighted fits, subset ordering, read guards, handle closure,
unchanged input hashes, reproducible shared starts, stockR likelihood projection
and K = 1 agreement. Production package defaults and installed packages are unchanged.

## Run / resume

From the assigner repository:

```sh
Rscript inst/experimental/test-streaming-engine.R
Rscript inst/experimental/test-stockr-paired.R
Rscript inst/experimental/run-stockr-marker-benchmark.R NEW_OUTPUT fixed 2
Rscript inst/experimental/run-stockr-marker-benchmark.R NEW_OUTPUT cv 2
```

The optional fourth argument is a comma-separated list of job row numbers.
Completed jobs are reused. Source hashes and stockR version must match the saved
audit when resuming. Failures are saved, not silently dropped or retried. Use a
fresh output directory for changed source or configuration. The fixed and CV
stages have separate results and progress logs.

Local output: project workspace `artifacts/stockr-marker-benchmark-20260916`.
Execution status and measured results will be recorded separately after each stage.
