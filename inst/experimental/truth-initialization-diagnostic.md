# Truth-initialization diagnostic — 16 September 2026

## Main finding

Better initialization helps recover the simulated groups, particularly with
500 markers, but does not solve the current predictive-K problem. The blind
solutions already have higher fitted objectives than the truth-started solutions
(apart from a numerical tie in one strong-structure replicate). More starts
aimed solely at maximizing that objective therefore cannot be assumed to recover
the generating groups.

This is not evidence of a software bug or that maximum likelihood is inherently
invalid. It is evidence that finite-sample fit, recovery of generating groups,
and prediction under estimated parameters differ in this experiment.

## What was run

Ten existing GDS datasets: two replicates each of generating drift 0.001, 0.01,
0.1, 300 individuals at drift 0.01, and 500 markers at drift 0.01. Other settings
are 90 individuals, equal population sizes and 100 markers. These are clean,
near-certain GL inputs, not a low-depth or missingness test.

Three comparisons use the same observations:

1. The saved blind, seven-start K = 3 fit.
2. A single K = 3 EM fit initialized at the true generating frequencies and
   proportions. Frequencies are then estimated freely. The pooled regularizer
   (strength 2), numerical bounds, tolerance and 300-iteration limit are unchanged.
3. A privileged reference that keeps true frequencies and proportions fixed,
   integrating over unknown individual group membership. It never fits parameters.

The truth start was also refitted with each of the original two training folds.
Held-out individuals have zero training weight; pooled prior centers are estimated
from training observations only. However, initialization uses population
frequencies unavailable in real analyses. These are privileged diagnostic
predictions, NOT ordinary cross-validation or a deployable K-selection method.
True labels are used for ARI evaluation, not to constrain fitted membership.

All ten full-data truth-start fits and all 20 fold fits converged without warnings
or errors. Source checksums matched the original benchmark engine. Input GDS
checksums remained unchanged. Tests cover label-invariant ARI, ties, singleton
behavior, and equality of frozen K = 1/K = 3 scores when component frequencies
are identical. No production mixture code, defaults or installed package changed.

## Recovery of the generating partition

ARI is label-invariant; 1 is perfect agreement and values near zero indicate
chance-level agreement. Ranges below cover two replicates, not confidence intervals.

| Design | Blind fit ARI | Truth-start fit ARI | Fixed true-model ARI |
| --- | --- | --- | --- |
| Drift 0.001; 90 samples, 100 markers | 0.003–0.013 | 0.003–0.008 | -0.006–-0.003 |
| Drift 0.01; 90 samples, 100 markers | 0.006–0.074 | 0.211–0.304 | 0.433–0.561 |
| Drift 0.1; 90 samples, 100 markers | 0.934–1.000 | 0.966–1.000 | 0.966–1.000 |
| Drift 0.01; 300 samples, 100 markers | 0.045–0.071 | 0.163–0.177 | 0.394–0.460 |
| Drift 0.01; 90 samples, 500 markers | 0.123–0.149 | 0.901–1.000 | 0.934–1.000 |

The fixed true model is a reference, not a guaranteed upper bound on realized
ARI. Individual classification accuracy and evidence for population number are
also different quantities; ambiguous individual assignments alone do not prove
that K cannot be inferred.

## Prediction: separate signal from parameter-estimation costs

It would be misleading to compare only true K = 3 frequencies with an estimated
K = 1 frequency vector. Knowing the frequencies can improve scores even without
population structure. We therefore added a frozen K = 1 reference using the
true mixture-weighted pooled allele frequencies under the single-population HWE
model. Both frozen references use known frequencies, not estimated ones.

Scores below are mean log-likelihood differences per individual for the same
observations. Positive favors K = 3. Ranges again represent two replicates.

| Design | Truth-started fitted K3 minus fitted K1 | Frozen true K3 minus frozen pooled K1 |
| --- | --- | --- |
| Drift 0.001; 90 samples, 100 markers | -2.063 to -1.709 | -0.00605 to 0.00012 |
| Drift 0.01; 90 samples, 100 markers | -2.017 to -1.330 | 0.170 to 0.298 |
| Drift 0.1; 90 samples, 100 markers | 3.362 to 3.444 | 5.067 to 6.262 |
| Drift 0.01; 300 samples, 100 markers | -0.826 to -0.628 | 0.199 to 0.221 |
| Drift 0.01; 90 samples, 500 markers | -9.572 to -8.929 | 1.894 to 2.122 |

Truth initialization generally improved K = 3 prediction over the blind K = 3
fits, but K = 1 still predicted better in every weaker-signal design. The
strong-structure results remained favorable to K = 3. This comparison adds a
privileged K = 3 start; it does not rerun a new complete K = 1:5 search.

For 500 markers, the known-frequency model separates groups and out-scores the
known pooled reference, while finite-training-sample estimates lose that
predictive advantage. This points to estimation variance/overfitting as a major
limitation alongside initialization. It does not establish an optimal prior,
a universal sample-size threshold, or impossibility of blind K recovery.

The apparent certainty is also informative: in the weakest 100-marker case,
the true fixed model has mean maximum membership about 0.46–0.47. After fitting,
the truth-started model reports about 0.99 despite chance-level recovery. These
in-sample conditional memberships must not be advertised as calibrated accuracy.

## Recommended next experiment

1. Generate nested sample AND marker subsets from the same population-frequency
   realization, increasing both dimensions together. Use separate training and
   final evaluation replicates to keep design choices away from final assessment.
2. Compare two-fold with repeated five-fold prediction: retain more individuals
   per group for parameter estimation. Report score differences and consistency,
   not merely the winning K. Use identical fold assignments across methods.
3. Evaluate a prespecified range of shrinkage strengths and a hierarchical
   frequency model as candidates. The current strength 2 is not established as
   adequate. Stronger shrinkage could reduce estimation noise but could also
   remove weak structure; choose settings within training data, not from the
   final held-out or truth-label evaluations.
4. Benchmark an independent likelihood-aware method on matched data and folds.
   Profile and validate performance improvements before multiplying long fits.

Do not select a method merely because it returns the known simulated K = 3.
Include one-/two-population controls, private alleles and later model
misspecification tests. Scott Foster's methodological review remains necessary.

## Reproduction and outputs

From the assigner repository:

```sh
Rscript inst/experimental/test-truth-diagnostic-helpers.R
Rscript inst/experimental/diagnose-truth-starts.R ORIGINAL_OUTPUT EXTENDED_RESULTS NEW_DIAGNOSTIC 4
Rscript inst/experimental/summarize-truth-starts.R ORIGINAL_OUTPUT NEW_DIAGNOSTIC NEW_SUMMARY
```

Existing destinations cause errors; original inputs/results are read-only.

Raw fits, per-fold records and audit:

`/Users/gos040/.codex/.chatgpt-projects/g-p-6a8255ffbc948191b3592f138fe06cf5/artifacts/grur-truth-diagnostic-20260916`

Combined summary and per-individual frozen-reference scores:

`/Users/gos040/.codex/.chatgpt-projects/g-p-6a8255ffbc948191b3592f138fe06cf5/artifacts/grur-truth-summary-20260916`
