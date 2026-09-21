# stockR comparison: larger marker panels

**Scope:** this page reports the shared-start comparison, not native stockR
initialization/AIC/BIC. See [native workflow comparison](stockr-native-comparison.md)
for the subsequent native audit and experiment. A K = 1 choice under the shared
CV protocol must not be described as a failure of stockR's native workflow.

## Completed known-K stage

All 40 panels completed: four marker counts, five generating conditions and two
independent frequency realizations. All 120 model fits converged; no starting
points failed. K was supplied in this stage, not selected.

Each panel uses 300 training and 300 separate evaluation individuals. The methods
are stockR EM, an unpenalized near-certain-GL mixture, and a pooled-regularized GL
mixture with fixed strength 100. Seven shared blind starts are used for K > 1.

Mean independent-test ARI for stockR across the two three-population realizations:

| Generating drift | 1,000 markers | 3,000 markers | 5,000 markers | 10,000 markers |
| --- | ---: | ---: | ---: | ---: |
| 0.001 | 0.002 | -0.0005 | 0.005 | 0.017 |
| 0.005 | 0.533 | 0.885 | 0.995 | 1.000 |
| 0.01 | 0.980 | 1.000 | 1.000 | 1.000 |

The unpenalized GL results match these partition scores except at drift 0.005,
5,000 markers, where its mean ARI is 0.990 rather than 0.995. This supports
agreement in the hard-genotype limit, not identical behavior under every bound
or stopping rule. Strength 100 makes small, nonuniform changes in group recovery.

High conditional membership remains potentially misleading: at drift 0.001 and
10,000 markers, stockR and unpenalized GL have mean maximum membership 0.956 and
0.944 across the two replicates despite ARI near zero. Strength 100 lowers these
values to about 0.916 but does not restore correct grouping. Membership certainty
conditional on a fitted model is not demonstrated assignment accuracy.

The two-population control at drift 0.01 has mean ARI 0.980 at 1,000 markers and
1.000 at the larger panels in all three models. The null supplied with K = 1
is necessarily trivial and is NOT evidence against false-positive K selection.

ARI measures agreement adjusted for chance, not percent correctly classified.
Generating drift is not identical to realized WC84 FST. Measured WC84 for the
10,000-marker three-group panels was 0.000958/0.000989, 0.005018/0.005022, and
0.009968/0.009938 for the 0.001, 0.005 and 0.01 targets, respectively.
The per-panel values are saved in `measured-WC84.rds`. These results do not
establish universal marker requirements or performance under uncertain GLs.

## Exploratory all-marker initialization check

The main starts use at most 1,000 markers for PCA. After observing poor recovery
at drift 0.001, we repeated its two 10,000-marker cases using ALL training markers
for the three PCA starts, while preserving the four random starts. This is a
post-primary exploratory analysis, not fresh held-out validation of an optimized
initialization method. The test samples remain independent of training, but their
earlier results motivated the sensitivity check.

| Model | Replicate 1 ARI | Replicate 2 ARI |
| --- | ---: | ---: |
| stockR | 0.0682 | 0.0376 |
| Unpenalized GL | 0.0682 | 0.0376 |
| GL strength 100 | 0.0678 | 0.0733 |

All six fits converged without failed starts. Using more markers for initialization
alone did not recover the weak groups. This does not establish global optimality
or an information-theoretic impossibility. More independent samples, other
optimization approaches and stronger independent controls remain worth testing.

### Privileged known-frequency diagnostic

Using the true generating frequencies at drift 0.001 and 10,000 markers gives
test ARI = 1 in BOTH replicates. Mean test log-likelihood gains over the fixed
true pooled one-population model are 5.2311 and 5.3305 per individual. Both
references use known frequencies; this avoids comparing a known-frequency model
against an estimated one and attributing all improvement to structure.

This demonstrates signal when the generating parameters are known. It suggests
learning/estimation or optimization limits in the blind fits, not a complete
absence of information. It does not show that a blind estimator can recover
those parameters from these training cohorts. Truth was used only in this
separately labelled diagnostic, never for the benchmark fits or CV selection.
Outputs: `artifacts/stockr-marker-oracle-20260916`.

## Unknown-K stage

The six-panel five-fold predictive-selection pilot is complete: 480 CV candidate
fits and 2,640 starts, with no failed starts or recorded fit warnings. Candidate
K was 1–4; GL jointly selected pooled strength from 2, 20 and 100. Both methods
used the same external predictive criterion, not native stockR AIC/BIC selection.
Choices were locked before projecting the independent test cohort.

| Generating condition | Markers | stockR selected K | GL selected K (strength) | stockR test ARI | GL test ARI |
| --- | ---: | ---: | ---: | ---: | ---: |
| K = 3, drift 0.001 | 10,000 | 1 | 1 (2) | 0 | 0 |
| K = 3, drift 0.005 | 10,000 | 1 | 3 (100) | 0 | 1 |
| K = 3, drift 0.01 | 1,000 | 3 | 3 (100) | 0.980 | 0.990 |
| K = 3, drift 0.01 | 10,000 | 3 | 3 (100) | 1 | 1 |
| K = 1, null | 1,000 | 1 | 1 (2) | 1 | 1 |
| K = 1, null | 10,000 | 1 | 1 (2) | 1 | 1 |

Both null controls select K = 1. Their perfect ARI is trivial once K = 1 is
selected; two null panels do not establish a false-positive rate. These are
exploratory results with one realization per displayed condition. The remaining
34 of 40 selection panels have NOT run; their performance must not be inferred
from known-K recovery. Strength 100 is the upper grid edge, not an established
optimum or recommended default. Replication, uncertain GLs, missingness, LD and
unbalanced samples remain to be tested.

An important completed contrast is drift 0.005 at 10,000 markers: the external
predictive criterion selects K = 1 for stockR but K = 3 with GL strength 100.
The selected regularized model has test ARI = 1. Its CV score improves by 4.939
log-likelihood units per individual over K = 1, whereas stockR's K = 3 score is
15.616 units below K = 1. Both methods recover the groups when K = 3 is supplied.

This is evidence for the benefit of regularization in this realization, NOT
evidence that GL encoding itself outperforms called genotypes: these GLs are
near-certain, and the model families differ in shrinkage. It is also not a claim
about native stockR AIC/BIC selection or every stockR optimization method.

## Reproducibility and limitations

See [design and implementation](stockr-marker-benchmark.md). Source snapshots,
hashes, seeds, complete fit objects, per-start diagnostics, predictions and GDS
checksums are saved with the run. stockR is version 1.0.76. Only ordinary EM with
shared starts was tested, not every stockR initialization/annealing option.

The streaming reader's maximum observed likelihood block was 75,000 sample-marker
cells (300 x 250), including at 10,000 markers. This is not total peak RAM.
The comparison still holds paired hard-call matrices for the simulation and
initialization; generic GL workflows must not derive starts or priors from truth.

Plots and detailed results are in the local project workspace:

- `artifacts/stockr-marker-summary-20260916-v2/known-K-recovery.png`
- `artifacts/stockr-marker-summary-20260916-v2/fixed-results.txt`
- `artifacts/stockr-marker-benchmark-20260916/`
- `artifacts/stockr-all-marker-starts-20260916/`
- `artifacts/stockr-selection-summary-20260916/selection-results.txt`
- `artifacts/stockr-selection-summary-20260916/predictive-K-selection.png`

The vectorized normalization and EM block updates have now been integrated into
the existing internal implementation. All six focused regression-test files pass,
including scalar-reference moments and weighted/block-size invariance checks.
The integrated engine agrees with preserved pre-change code at strengths 0/2/100.
The actual GDS interface also reproduces the saved 10,000-marker fit within
0.000001 for frequencies and membership, without changing the GDS checksum.

The benchmark was launched before that integration and retains its exact source
snapshot in `source-snapshot.rds`; `engine-before-scaling.rds` preserves the old
core for comparison. Source-hash guards intentionally require a new run directory
after code changes. Running jobs retained their originally loaded helpers. A
subsequent reporting-only correction accounts for stockR's early-convergence
trace truncation: saved stockR start iteration counts can be one too small, but
iteration-cap classification, selected fits and numerical results are unaffected.
The helper PL error text was also aligned with the existing validation wording.
Future paired benchmark runs call the integrated `.sm_block_fit` and
`.sm_block_estep` directly; the matched-start/projection regression test passes
with that implementation.

Package arguments, defaults and installed packages remain unchanged.
Scott Foster's review is still required before public release.
