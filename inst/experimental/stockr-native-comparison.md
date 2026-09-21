# Native stockR workflow comparison

## Completed results

All 40 panels completed: 400 native calls, no recorded warnings/errors or
trace-classified iteration-cap exits. Single-call and best-of-three protocols
selected the same K in every panel.

### Selecting K using native AIC and BIC

**Marker panels: 1,000, 3,000, 5,000 and 10,000 SNPs, analysed separately.**
Each generating condition had two frequency realizations, with 300 training fish
and 300 independent test fish per panel. Candidate K ranged from 1 to 4;
AIC/BIC selection used training data only. Marker panels were nested within each
realization, not additional independent simulation replicates.

| Generating condition | True K | SNP panels tested | Native AIC K | Native BIC K |
| --- | ---: | --- | ---: | ---: |
| Drift F = 0.001 | 3 | 1,000 / 3,000 / 5,000 / 10,000 | 1 | 1 |
| Drift F = 0.005 | 3 | 1,000 / 3,000 / 5,000 / 10,000 | 1 | 1 |
| Drift F = 0.010 | 3 | 1,000 / 3,000 / 5,000 / 10,000 | 3 | 1 |
| Drift F = 0.010 | 2 | 1,000 / 3,000 / 5,000 / 10,000 | 2 | 1 |
| No population structure | 1 | 1,000 / 3,000 / 5,000 / 10,000 | 1 | 1 |

The displayed selections were identical across all four marker panels and both
frequency realizations, under both single-call and best-of-three native fitting
protocols. Drift F is the simulation setting, not a measured WC84 FST value.

At drift 0.005, supplying K = 3 gives mean independent-test ARI of 0.465,
0.970, 1.000 and 1.000 at 1,000/3,000/5,000/10,000 markers respectively
(best-of-three native fits). Thus the 10,000-marker result is a selection
limitation under these criteria, not failure to find the groups given K.
At drift 0.001, known-K recovery remains near chance even with 10,000 markers.

The six-panel direct comparison with existing GL CV results is now available.
At drift 0.005/10,000 markers, native stockR AIC and BIC choose K = 1 whereas
regularized GL predictive CV chooses K = 3 with test ARI = 1. At drift 0.01,
native AIC and GL CV both choose K = 3; native BIC chooses K = 1. These compare
different selection workflows and do not isolate GL encoding or establish
general superiority. Only one realization per condition has paired GL CV.

Native regression checks passed for actual defaults, returned IC formulas,
training likelihood reproduction and sample-order-preserving projection.
No package defaults, installed packages or public mixture APIs were changed.

Detailed reproducible output is in the local project workspace under
`artifacts/stockr-native-20260916-complete/`: `native-results.md`,
`paired-native-vs-GL-CV.md`, `summary.rds`, source audit and per-panel fits.
The initial interrupted run is retained separately in `stockr-native-20260916`;
it is not included in these results.

## Summary of K-selection approaches explored

These are different workflows, not interchangeable versions of native stockR.
Drift F below denotes the simulation setting. Recovery with K supplied and
selection of K are separate outcomes.

| Approach | Implementation and scope | Main observation |
| --- | --- | --- |
| Native AIC | stockR 1.0.76; all 40 panels; K = 1–4 | Selected the generating K at drift 0.01, but K = 1 at 0.005 and 0.001. |
| Native BIC | Same native fits and panels | Selected K = 1 throughout, including structured panels with excellent known-K recovery. |
| Five-fold predictive CV | External held-out scoring applied to stockR with shared blind starts; six panels, one realization per displayed condition | Selected K = 1 at drift 0.005 with 10,000 SNPs; selected K = 3 at 0.01 with 1,000 and 10,000 SNPs. |
| Predictive CV with pooled regularization | Our near-certain-GL prototype; same six-panel comparison; jointly tuned K and strength | Selected K = 3 at drift 0.005 with 10,000 SNPs, but not at 0.001. This tests shrinkage as well as selection, not a benefit of GL encoding alone. |
| Repeated-fold diagnostics and one-standard-error heuristic | Later GL-prototype experiments, not native stockR benchmarks | Did not establish a calibrated solution to low-coverage null splitting. Fold agreement alone did not certify the correct K. |

The first four rows used complete simulated genotypes (near-certain likelihoods
for the GL prototype). The later GL work examined sequencing uncertainty and
different designs; its failures must not be attributed to native stockR.
Additional starts and tighter convergence were optimization checks, not separate
K-selection criteria. Neither fixed-K bootstrap membership uncertainty nor
known-K assignment accuracy validates a K-selection rule.

See [shared-start predictive comparison](stockr-marker-results.md) and
[later prototype evidence review](stock-mixture-evidence-review.md).

## Three options for further testing

These are research proposals, not demonstrated fixes or changes to stockR.

1. **Calibrate evidence for an extra group by simulation.** Define a predictive
   gain statistic, simulate under the simpler fitted model, and rerun the entire
   fitting/selection procedure to estimate how often noise produces that gain.
   Start with K = 1 versus a split, then investigate sequential comparisons.
   Match sample size, marker count, depth and missingness; separately stress-test
   linkage and misspecification. Account for searching across candidate K values.
   Validate any threshold on new null and structured datasets. A nonsignificant
   gain means insufficient evidence, not proof of panmixia. This targets false
   splits; it may also reduce sensitivity to weak real structure.

2. **Test regularized prediction in a controlled, nested comparison.** Use the
   same genotype representation, starts, folds and candidate K grid, comparing
   unregularized and pooled-shrinkage fits. Tune strength and K within training
   data, then assess the locked workflow on untouched outer data. Begin with
   complete GT before extending to GL/PL. This isolates the promising shrinkage
   effect at drift 0.005 without claiming encoding superiority. Retain native
   AIC/BIC as baselines; do not insert a shrinkage fit into the ordinary parameter
   count formula without methodological justification. The existing six-panel
   pilot motivates this test but is not independent confirmation.

3. **Add resampling-based support and an explicit unresolved outcome.** Measure
   whether assignments of a common held-out cohort recur across training-sample
   resamples and independent marker blocks, using label-invariant agreement.
   Assess support alongside predictive gain, rather than forcing a single K.
   Include technical-batch and missingness sensitivity checks: a stable artifact
   can still be wrong. Calibrate any support thresholds on independent simulations
   and report both incorrect selections and abstentions. This extends beyond the
   start agreement and fold-winner diagnostics already implemented.

Suggested order: establish null calibration, perform the matched shrinkage
comparison, then evaluate whether stability-based abstention adds value. Across
all three, retain 1,000/3,000/5,000/10,000-marker panels, report realized WC84 FST,
use fresh frequency realizations, and score K recovery separately from assignment
accuracy. No tested result supports a universal failure threshold at FST < 0.02.

## Why this is a separate comparison

The earlier shared-start EM/predictive-CV comparison did not evaluate the
recommended stockR initialization or its returned information criteria.
Do not call a CV preference for K = 1 a failure of native stockR K selection.

## Source audit: stockR 1.0.76

Inspected the installed `stockSTRUCTURE` help, entire introductory vignette,
`stockSTRUCTURE`, `set.control`, `get.start.grps`, `getReturnObject`, `K1special`
and `getIC`. Public reference:
https://cran.r-universe.dev/stockR/doc/manual.html

- Recommended default: `Kmeans.EM`, `start.grps = NULL`, 25 k-means starts,
  up to 100 PCs calculated using all retained markers.
- Native source uses uncentered, unscaled PCA; we do not replace it.
- Actual default convergence settings: tolerance 0.00001, maximum 100 EM
  iterations, minimum 5, initial tau step 0.5. The last two differ from the help
  text's 1 and 1. We use installed behavior and retain each resolved control list.
- Returned AIC = -2 logL + 2p; BIC = -2 logL + log(N) p, where
  p = K * retained markers + K - 1 and N is the number of fish.
  Each returned criterion is verified numerically, not redefined.
- The vignette supplies K and recommends bootstrapping assignment uncertainty;
  it does not establish one universally preferred automatic K-selection rule.
  Accordingly we report AIC and BIC separately, not a single invented native rule.
- Native fitting returns no explicit convergence flag. Iteration-cap exits are
  flagged using the audited trace convention and excluded from selection.

## Locked native run

Same 40 nested panels, 300 training and 300 independent test fish as the earlier
benchmark: four marker counts, five conditions, two frequency realizations.
This is exploratory replication, not a precise power estimate. Reused datasets
are not a fresh confirmation set for the previously investigated GL method.

For each panel, fit K = 1–4 with unmodified native defaults except quiet output.
For K > 1, run three independent seeded native calls. Report both:

1. First call only: native single-run behavior.
2. Best converged likelihood of three calls per K: explicitly enhanced repeat
   protocol, NOT a package default. Each call already uses 25 k-means starts.

K = 1 uses its native analytic solution. No population labels, true allele
frequencies or test observations enter initialization or selection. Training-only
constant-marker removal is shared across candidate K. Each fish is its own
sampling group; no known breeding-group information is supplied.

Lock AIC/BIC choices before test projection. An incomplete converged K grid
withholds selection. Exact IC ties prefer the smaller K. Report known-K recovery
separately and retain warnings, errors, controls, likelihoods and timings.
Saved audit includes exact runner/helper sources and session information.

## Interpretation boundaries

Increasing the number of markers also increases the number of separately
estimated allele-frequency parameters. Moving from K = 1 to K = 3 adds
2M + 2 parameters for M retained markers. At M = 10,000 and N = 300, AIC
therefore requires a log-likelihood improvement greater than 20,002, and BIC
greater than approximately 57,043.5. These are consequences of the native
formula, not implementation errors. More markers need not automatically change
the preferred K when the gain per marker remains below the penalty per marker.
This explains a possible separation between accurate conditional classification
and conservative K selection; it does not prove the criterion is universally
inappropriate or that regularized predictive CV is universally better.

Compare native AIC/BIC with GL predictive selection only as different workflows,
not an isolated estimator comparison. The matched-method benchmark remains the
appropriate comparison for that second question. GL selection currently covers
six panels only; a full paired workflow power comparison requires completing the
remaining GL grid plus new, independent simulation replicates after freezing
the GL procedure. Annealing alternatives, bootstrap uncertainty, missingness,
uneven groups and uncertain GLs remain distinct extensions, not tested here.

## Reproduction

From the assigner repository:

```sh
Rscript inst/experimental/run-stockr-native.R EXISTING_MARKER_BENCHMARK NEW_OUTPUT
```

The runner uses two workers, refuses existing output directories, and records
per-panel progress and fit objects. It does not modify package defaults or data.
