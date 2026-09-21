# Evidence and interpretation review — 18 September 2026

For Scott Foster's methodological/code review. Experimental work, not endorsement
or a released stock-delineation method. This is a discrete mixture, not admixture;
genetic groups do not alone establish demographically independent management stocks.

## Consolidated findings

- Coverage pilot: 3,000 independent SNPs, 150 training fish (50/group for K=3),
  300 independent test fish, two frequency realizations per condition. Drift F
  is a simulator parameter, not measured WC84 FST. There is no universal 10× cutoff.
- At F=0.005, mean test ARI with K supplied rose from 0.569 at 3× GL to 0.841
  at 5× and 0.936 at 20×. Matching maximum-likelihood calls scored 0.307,
  0.773 and 0.941. High-coverage differences are small; this is not proof of GT superiority.
- Strict normalized-likelihood filtering discarded substantial low-coverage
  information. This was a controlled calling rule, not production-caller GQ.
- Automatic K selection can split null data at low coverage. A separate six-null
  1× experiment split five realizations under the tested maximum-CV-score rule.
  These counts are not general false-positive estimates. One-SE selection,
  tighter convergence and additional starts did not provide a validated solution.
- At F=0.001, true genotypes still gave chance-level blind assignment. But frozen
  generating frequencies gave mean independent-test ARI 0.806. Correct training
  labels plus estimated frequencies gave only 0.146 with strength 500. Finite
  training-frequency estimation is an important limitation, not proof of absent signal.
- True-frequency initialization did not close that gap. At F=0.001, the training
  objective preferred blind solutions over better-generalizing truth-start fits.
- Conditional beta-binomial prediction and a joint mean-field variational pilot
  did not rescue weak-group recovery. They used true GT, K in {1,3}, training-only
  tuning and the same two null/weak/moderate realizations. They do not resolve GL
  false splits or constitute independent confirmation. Neither is adopted.

## Reporting decisions now implemented

Keep the numerical predictive K winner separate from descriptive support. Expose
fold winners including ties, score gaps, boundary/incomplete-comparison flags,
and fit warnings. Report label-invariant start partition agreement with comparable
sample counts and explicit memory-limited non-assessment. Give per-individual
margin, entropy and status, without a scientifically unvalidated confidence cutoff.

No report field proves the chosen K, calibrates individual assignment accuracy,
or establishes management-stock independence. Across-fold membership stability,
repeated-fold stability and calibration remain distinct from fold K-score agreement.

## Still required before release

1. Independent null and structured replication, especially under GL, uneven depth,
   population-structured missingness, linkage, relatedness and unequal sampling.
2. Evidence-backed rules for operational abstention; do not turn descriptive
   metrics into arbitrary hard thresholds or treat K=1 as a biological conclusion.
3. Review likelihood conventions, regularization, selection objective, uncertainty
   calibration and public API/documentation with Scott.

## Reproducibility record

Full RDS outputs and runnable scripts remain in the development workspace under
`artifacts/`: `coverage-sensitivity-20260917`, `gl-gt-comparison-20260917`,
`independent-null-20260917`, `weak-mixture-diagnostic-20260918`,
`frequency-prediction-20260918`, and `joint-vb-20260918`. Each contains audit/settings
records and panel outputs. These large local artifacts are not bundled with the
package; provide them separately for review. Historical source hashes describe
the benchmark versions, not subsequent reporting-only changes.
