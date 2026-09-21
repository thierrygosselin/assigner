# Sample size and regularization: controlled follow-up

## Completed results

All 15 fresh panels completed. Both fixed100 and tuned GL workflows selected
K = 3 in all nine structured panels. Therefore this fresh grid does NOT show an
incremental K-recovery benefit from tuning: fixed100 already succeeded here.
It does show a paired sample-size benefit for assignment accuracy:

| Fish per population | Total training N | Mean test ARI, fixed100 | Mean test ARI, tuned | Native AIC correct K |
| --- | ---: | ---: | ---: | ---: |
| 100 | 300 | 0.9440 | 0.9343 | 0 / 3 |
| 200 | 600 | 0.9933 | 0.9933 | 3 / 3 |
| 300 | 900 | 1.0000 | 1.0000 | 3 / 3 |

All 8,250 CV starts, 138 final GL starts and 150 native fits passed their
convergence checks; no native or CV-initialization warnings were recorded.
All final GDS likelihood comparisons and input checksum checks passed.
The full audit and individual rows are in `strength-size-results.md` and
`summary.rds` in the verified output folder.

Tuning selected strength 200 in all nine structured panels and improved their
independent-test mean log likelihood by 0.161–0.983 units per fish relative to
fixed100. Improved predictive density is not necessarily better classification:
at N = 300, two tuned partitions had slightly lower ARI. At larger N the
partitions agreed. Native BIC selected K = 1 in all structured panels.

All six nested null panels selected K = 1 in every workflow. These represent
two independent null realizations, not six independent controls. At K = 1 all
tested strengths give the same pooled prediction; strength 20 wins by the
declared tie rule, not evidence of a biologically preferable null strength.

### Separate revisit of the earlier missed realization

This is a post-hoc development diagnostic, not part of the fresh validation.
Original folds and initial partitions were retained, and the strength100 CV
scores reproduced the previous GDS results to within 0.000000000119.

| Markers | Previous fixed100 K | Jointly selected K / strength | Independent-test ARI |
| --- | ---: | --- | ---: |
| 3,000 | 1 | 3 / 500 | 0.9601 |
| 10,000 | 2 | 3 / 500 | 1.0000 |

These two panels share one previously identified difficult realization. Stronger
shrinkage helped it, but 500 is the upper grid boundary; it is not established as
optimal or as a new default. The grid was unchanged after seeing the results.
The diagnostic is reproducible with `revisit-strength-failures.R` and is retained
under local workspace `artifacts/strength-failure-revisit-20260916`.

**Practical conclusion:** prioritize adequate population sampling and allow
training-only strength selection rather than assuming one fixed value is best.
Keep predictive scoring separate from classification accuracy. Improved
initialization and repeated CV remain separate, untested changes in this follow-up.
No production defaults or public APIs were changed.

## Locked design

Question: which improves predictive K recovery near FST 0.005: more training
individuals, selecting pooled regularization strength, or both?

- 3,000 markers held fixed, rather than confounding sample and marker counts.
- Three new K = 3 frequency realizations at drift 0.005; two new K = 1 nulls.
- Structured training groups contain 100, 200 or 300 fish each (N = 300/600/900).
  Nulls use the same TOTAL sample sizes; they contain one population.
- Training subsets are nested within each realization. The same independent
  300-fish test cohort is used across its sample-size comparisons.
  Fifteen panels are five independent realizations, NOT fifteen replicates.
- Existing grur generator, independent complete HWE loci, seeds 2480001–2480005.
- Near-certain GLs stored and read from GDS, not uncertain-likelihood validation.
- Fivefold predictive validation selects K = 1–4 and strength from
  20, 50, 100, 200, 500. Also retain fixed strength 100 with K selected on exactly
  the same folds. No test-sample scores or generating labels select parameters.
- Training-only invariant-marker masks and pooled frequencies. All strengths
  share the SAME seven initial partitions for each fold/K. The existing
  PCA 2/5/10-PC plus four-random-start recipe is unchanged.
- Equal folds; mean held-out log likelihood; ties within 0.00000001 select
  smaller K, then weaker strength. Incomplete grids withhold selection.
- Native stockR Kmeans.EM with three independent native calls per K; best
  converged likelihood and returned AIC/BIC provide the comparator at every N.
- No repeated-CV, initialization or adaptive-shrinkage change in this experiment.
  Those are separate follow-ups, not simultaneous modifications.

## Validated computational shortcut

The new `benchmark-hard-limit.R` helper is ONLY an experimental accelerator for
complete vectors with one GL = 0 and the alternatives GL = -80. It reads and
validates these from GDS in bounded blocks, then stores a bounded call matrix
(maximum here 900 x 3,000). It rejects missing/uncertain vectors. It is not an
imputation method, production replacement or extension to arbitrary GDS inputs.

In this limit the conditional posterior genotype equals its observed count to
numerical precision. Binomial log likelihood and weighted counts can therefore
be calculated using matrix products. The prior, updates, starting probabilities,
bounds, objective and stopping criterion are unchanged.

`test-hard-limit.R` compares frequencies, proportions, membership, log likelihood
and penalized objective with the actual streaming GDS EM engine across K = 1–3
and strengths 0/20/50/100/200/500. Comparisons pass within 0.0000001; uncertain and
missing inputs are rejected. Every final experimental fit is additionally
checked against the streaming GDS likelihood evaluator at tolerance 0.00000001.
Input checksums are verified after each panel. No production defaults changed.

## Interpretation and outputs

The fresh fixed100 versus tuned comparison isolates the selection of shrinkage
within a fixed fitting recipe. Nested sample sizes make paired comparisons
possible, but the small number of independent realizations is not a power study.
Selection at strength 500 is an edge-of-grid result, not proof that 500 is optimal.
Nulls check obvious overclustering, not a precise false-positive rate.

Records include source snapshots/hashes, session info, seeds, folds, all CV scores
and start diagnostics, selected parameters locked before test access, native
fits, final predictions and warnings. GDS source and test data remain unchanged.

From the assigner repository:

```sh
Rscript inst/experimental/test-hard-limit.R
Rscript inst/experimental/run-strength-size.R NEW_OUTPUT 4
```

Output: local project workspace `artifacts/strength-size-20260916-verified`.
The initial interrupted output `strength-size-20260916` is retained for audit,
not included in results. Its final equality check rejected sample-name attributes
despite a maximum numerical discrepancy of approximately 0.0000000000255 in the
investigated panel. The corrected check verifies sample IDs/order explicitly and
then compares unnamed numeric scores. The same seeds/settings were rerun; this
was a validation-adapter correction, not model tuning.
The runner can reuse completed panels only when source hashes and stockR version
match. Saved failures require inspection; they are not silently retried.
