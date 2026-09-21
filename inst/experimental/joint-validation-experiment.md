# Joint sample/marker expansion and predictive regularization

## Prespecified first-stage design — 16 September 2026

This experiment follows the [truth-start diagnostic](truth-initialization-diagnostic.md).
It uses no true-frequency starts. All starts are blind seven-start mixed fits.

- Conditions: three equally sampled populations with generating drift 0.01,
  and a one-population null. One new independent frequency realization per
  condition, seeds 83211 and 83212. This is exploratory, not a power study.
- Nested panels: 90 training individuals with 100 markers, and 300 training
  individuals with 500 markers. The smaller panel is an exact sample/marker
  subset of the larger panel. Both dimensions change together; this tests joint
  expansion, not their separate causal effects.
- An independent cohort of 300 individuals is generated under the same
  population frequencies. Its first 100 markers form the small-panel test data.
  Its genotypes are never used for fitting, pooled prior centers, initialization
  or hyperparameter selection. GDS sample IDs are disjoint from training IDs.
- Complete, near-certain GL observations isolate frequency-estimation and
  prediction issues. This is not a low-depth or missingness benchmark.
- Candidate K = 1:4; pooled regularization strengths 2, 20 and 100. Strength 2
  is the original baseline. No grid changes are made after viewing outcomes.
- Two-fold and five-fold training cross-validation, with five-fold predeclared
  as primary. There is one deterministic fold allocation for each; repeated
  random fold allocations are not yet tested. Round-robin folds balance the
  generating groups through sample order; labels do not enter the estimator.
- Each approach selects K and strength jointly by mean held-out, unpenalized
  log likelihood. Scores within 1e-8 per individual tie; prefer smaller K, then
  weaker regularization. An incomplete grid withholds selection.
- Also retain the best K at fixed strength 2 for each fold policy. All choices
  are written to a locked selection file BEFORE opening the evaluation GDS.
- Refit selected models on the full training cohort and predict the independent
  cohort. Report scores relative to a fitted K = 1 reference on the SAME test
  observations, test ARI, assigned counts and conditional membership values.
  Do not choose a preferred method or strength retrospectively from test ARI.
- Selected K = 4 is flagged as the search boundary. K = 3 is not forced.

An independent cohort conditional on the same frequencies measures prediction
to that population realization, not generalization across different evolutionary
histories. More independent simulation realizations are needed for the latter.

## Computation and verification

The experimental accelerator in `benchmark-vectorized-engine.R` vectorizes the
same EM equations over markers within each component. Normalized GL blocks are
cached only for small benchmark datasets (maximum one million sample-marker
cells). This is not a production-scale RAM guarantee or a new statistical model.
Production source functions and installed packages remain unchanged. An isolated
benchmark process substitutes the tested block update while retaining the
existing PCA/random/split/perturb initialization machinery.

Equivalence checks passed for one-row/one-marker edge cases, impossible genotype
states, flat likelihoods, uncertain GLs, missingness, zero-weight heldouts,
strengths 0/2/20/100, weighted EM and seven-start fits against the GDS implementation.
The design test checks incomplete-grid handling, deterministic ties, exact nested
panels and sample separation. No evaluation probabilities are claimed calibrated
merely because those numerical checks pass.

All five existing focused mixture regression-test files also passed with the
vectorized block update substituted in an isolated R session, including dense
versus blockwise calculations, GDS fixtures, initialization and validation tests.

The four datasets require 336 selected CV fold fits, each using seven starts,
plus selected full-training refits. Candidate-start failures and warnings remain
in saved diagnostics. Convergence does not establish a global optimum.

## Reproduction

Run from the assigner repository:

```sh
Rscript inst/experimental/test-benchmark-vectorized-engine.R
Rscript inst/experimental/run-joint-validation.R NEW_OUTPUT
Rscript inst/experimental/test-joint-validation-design.R NEW_OUTPUT
```

The output directory must not already exist. This run is stored at:

`/Users/gos040/.codex/.chatgpt-projects/g-p-6a8255ffbc948191b3592f138fe06cf5/artifacts/grur-joint-validation-20260916`

It contains locked design, source hashes, manifests, truth, separate train/test
GDS inputs, per-dataset progress logs, fold allocations, checkpointed start and
fold diagnostics, locked selections, final fits and evaluation predictions.
Input checksums are checked again after completion.

## Interpretation limits

The experiment can indicate whether joint expansion and stronger shrinkage help
under this controlled model. It cannot establish a universal minimum sample
size, validate real stock counts, demonstrate demographic independence, or
calibrate assignment probabilities across species and sequencing conditions.
Scott Foster's review and independent-method comparisons remain required.
