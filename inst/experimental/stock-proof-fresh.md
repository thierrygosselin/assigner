# Fresh-replicate proof of concept

## Completed findings

All 16 panels completed. At generating drift 0.005, recovery of the true K = 3
across three fresh frequency realizations was:

| Markers | Fixed-strength GL predictive selection | Native stockR AIC | Native stockR BIC |
| --- | ---: | ---: | ---: |
| 1,000 | 1 / 3 | 0 / 3 | 0 / 3 |
| 3,000 | 2 / 3 | 0 / 3 | 0 / 3 |
| 5,000 | 2 / 3 | 0 / 3 | 0 / 3 |
| 10,000 | 2 / 3 | 0 / 3 | 0 / 3 |

At 10,000 markers, the GL workflow selected K = 3, 2, 3, with independent-test
ARI 1, 0.570, 1. Native AIC selected K = 1, 2, 2 and native BIC selected K = 1
in all three cases. Supplying the generating K to native stockR instead gave
perfect test recovery in all three 5,000- and 10,000-marker panels. This separates
good conditional fitting from imperfect K selection.

Both fresh null controls selected K = 1 in all workflows. Both weak controls
also selected K = 1: differentiation near 0.001 remains unresolved here.
These two null realizations do not estimate a general false-positive rate.

Measured WC84 FST at 10,000 markers was 0.005009, 0.004993 and 0.005075 for the
intermediate realizations, and 0.000970 and 0.000942 for the weak controls.
The result is not explained by accidentally simulating a much higher FST.

All 1,760 CV starts and 64 final GL starts converged. All 160 native EM fits
passed the trace-based cap check. One native K = 4 start for intermediate
replicate 3 at 1,000 markers emitted a k-means ten-iteration convergence warning;
its subsequent EM converged, it was retained in the audit, and neither native
criterion selected K = 4. No warning was silently removed from the results.

**Proof-of-concept conclusion:** fixed pooled regularization with predictive
selection recovered weak structure in fresh datasets more often than the native
AIC/BIC workflows tested here, but still underestimated K in some replicates.
This supports further development, not a universal success rate or a benefit
attributable to GL encoding alone. No settings were retuned during this run.

Full individual results, known-K diagnostics and fit audits are in the output
folder's `proof-results.md` and `summary.rds`; measured values are retained in
`measured-WC84.rds`. The generation/fitting source snapshot is in
`design-and-audit.rds`. The native adapter and hard-call-limit projection tests
passed, as did the per-panel GL score and unchanged-input checks.

## Question

Can a fixed pooled-regularized mixture with predictive K selection recover weak
structure reproducibly? This is a proof of concept, not a claim of universal
superiority over stockR or proof of demographically independent management stocks.

## Locked before fitting

- Three fresh, independent allele-frequency realizations at generating drift
  0.005, K = 3, each with nested 1,000/3,000/5,000/10,000-marker panels.
- Two fresh weak controls (drift 0.001, K = 3) and two fresh null controls
  (K = 1), each with 10,000 markers. Total: 16 panels from seven independent
  frequency realizations. Nested panels are NOT independent replicates.
- 300 training fish and 300 separately generated test fish per realization.
  Balanced populations, independent biallelic HWE loci, complete genotypes.
- The same existing grur simulation engine; new seeds beginning at 1920101.
  Generating drift is not exact measured sample FST.
- Existing GDS GL fitting engine; near-certain GLs equivalent to called genotypes.
  No uncertain-GL, missingness or LD benefit is claimed.
- GL regularization strength fixed at 100 from the earlier exploratory pilot.
  Fivefold predictive validation selects only K = 1–4. Seven blind starts for
  K > 1, the same 2/5/10-PC and random-start recipe as the previous benchmark.
  No modification in response to results. Equal folds, training-only marker
  masks, scores averaged across folds; score ties within 0.00000001 prefer smaller K.
- Native stockR Kmeans.EM defaults, three independent native calls per K > 1,
  best converged likelihood per K. Native AIC and BIC reported separately.
  This best-of-three protocol is an explicit enhancement, not the package default.
- Choices saved before test projection. Labels/true frequencies never enter
  fitting or selection. Incomplete candidate grids withhold selection.
- Native convergence is trace-classified, as in the preceding audited comparison.
  GL convergence uses the existing engine flag. Failed starts remain recorded.
- Hard-call binomial projection is used for speed only in this near-certain-GL
  experiment; full-fit scores are checked against the GDS GL engine to tolerance
  0.00000001. This projection must not be reused for uncertain GL experiments.

## Outputs and interpretation

Retain design, source snapshot, session information, seeds, input checksums,
per-fold scores/start diagnostics, native fits, selected GL fits and independent
test predictions. Report the number of successful K recoveries rather than
presenting three replicates as a precise power estimate. Test ARI measures
chance-adjusted partition agreement, not percent correct assignments.

This compares practical workflows with different selection criteria. It does
not isolate likelihood encoding, initialization or regularization individually;
the earlier matched-method benchmark addresses a different question.

Reproduce from assigner:

```sh
Rscript inst/experimental/run-stock-proof.R NEW_OUTPUT 4
```

Output: local project workspace `artifacts/stock-proof-fresh-20260916`.
Production defaults and installed packages are unchanged.
