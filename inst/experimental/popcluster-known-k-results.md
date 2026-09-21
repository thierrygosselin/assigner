# PopCluster versus native stockR: matched known-K pilot

Completed 21 September 2026. The initial six-panel pilot and eight additional
weak-differentiation replicates all completed successfully.
Here, K is the number of populations. This is a training-recovery comparison,
not automatic K selection, prediction, evaluation of genotype likelihoods
(GL) or their phred-scaled representation (PL), or validation of a universal
Patterson–Price–Reich (PPR) threshold.

## Design and checks

- 300 training individuals: 100 per true population; 10,000 independent
  single-nucleotide polymorphisms (SNPs).
- True K = 3 supplied to both methods; two independent simulated datasets per
  generating differentiation parameter. These parameters are not measured
  population differentiation (FST) using the Weir–Cockerham (1984) estimator.
- Exact same training genotypes and marker order as the existing native stockR
  benchmark; all 10,000 markers survived its training-variability filter.
- PopCluster 1.7.0.0, build 20250715: Model=1, Scaling=0, AllSNP=1,
  DataForm=1, no population labels, pooled-frequency prior, three replicate runs.
- Native stockR results reused from the completed benchmark, three independent
  native initialization calls. Highest training likelihood selects the run for
  each program. Truth is used only for evaluation, never run selection.
- Verified sample IDs/order, missing counts, SNP counts, and every individual's
  heterozygous genotype count against the original matrices.
- The adjusted Rand index (ARI) measures agreement with true groups, adjusted
  for chance. It is invariant to arbitrary cluster numbering; 1 means perfect
  recovery and values near 0 indicate chance-level agreement.
  stockR training classifications were obtained from its fitted mixture parameters;
  no posterior ties occurred.

## Initial six-panel results

In the table, F is the generating differentiation parameter, N the total number
of training individuals, and L the retained independent marker count.

| Generating F | Replicate | PPR S = N L F² | PopCluster training ARI | stockR training ARI |
|---|---|---|---|---|
| 0.001 | 1 | 3 | 0.063873 | 0.029431 |
| 0.001 | 2 | 3 | 0.311616 | 0.033787 |
| 0.005 | 1 | 75 | 1 | 1 |
| 0.005 | 2 | 75 | 1 | 1 |
| 0.010 | 1 | 300 | 1 | 1 |
| 0.010 | 2 | 300 | 1 | 1 |

At F=0.001, mean ARI was 0.188 for PopCluster and 0.032 for stockR. This is
descriptive evidence from two datasets, not a population-level superiority claim.
Both programs recovered all individuals correctly in the four stronger panels.
The result is compatible with S>1 not guaranteeing accurate individual recovery;
it does not independently validate or invalidate PPR.

## Start variability and implementation caveats

For weak replicate 1, PopCluster's three ARIs ranged from 0.0550 to 0.2150,
but its highest-likelihood run had ARI 0.0639. For replicate 2 the range was
0.2822–0.3116 and the highest-likelihood run also had the highest ARI. Thus simply
reporting the most accurate run would give an optimistic, non-deployable result.

Even with Model=1, this binary performed and printed an additional admixture
analysis. Only the explicit ML_C clustering labels and clustering LogL were used
here. Total execution times (approximately 79, 82, 44, 39, 23 and 24 seconds)
include that extra processing and output, so are not a fair optimizer speed
comparison. Whether this extra analysis is intentional is a focused question
for Jinliang if needed; no defect is asserted.

PopCluster's hard-partition likelihood and stockR's marginal mixture likelihood
are different objectives. Do not compare their numerical likelihood values
directly or attribute all performance differences to their optimizers.

## Extension: ten independent weak-differentiation datasets

Eight new datasets were generated using fresh seeds (911003–911010), bringing
the weak-differentiation comparison to ten independent replicates. Each has
300 training fish, 100 per true group, 10,000 retained independent SNPs,
generating F = 0.001 and generating PPR S = 3. The same grur simulation engine
was retained, including ancestral allele frequencies uniform between 0.05 and
0.95; this is not an exact replication of Wang's uniform-between-0-and-1 design.
Each replicate has newly drawn population allele frequencies and genotypes,
not merely another initialization of an existing dataset.

Both methods used three runs at supplied K = 3. The run with the highest
eligible training likelihood was selected before evaluating against truth.
These are native method comparisons, not equal-computation-budget comparisons;
stockR's native initialization includes its own multiple k-means starts.

| Replicate | Measured differentiation | Measured PPR S | PopCluster training ARI | stockR training ARI |
|---|---:|---:|---:|---:|
| 1 | 0.0009582 | 2.754 | 0.0639 | 0.0294 |
| 2 | 0.0009888 | 2.933 | 0.3116 | 0.0338 |
| 3 | 0.0010348 | 3.212 | 0.1113 | 0.0322 |
| 4 | 0.0008654 | 2.247 | 0.0511 | 0.0251 |
| 5 | 0.0010936 | 3.588 | 0.2244 | 0.0199 |
| 6 | 0.0010707 | 3.439 | 0.3456 | 0.0073 |
| 7 | 0.0010603 | 3.373 | 0.3353 | 0.0666 |
| 8 | 0.0009708 | 2.827 | 0.1382 | 0.0398 |
| 9 | 0.0009039 | 2.451 | 0.1365 | 0.0386 |
| 10 | 0.0010427 | 3.261 | 0.2276 | 0.0149 |

Measured differentiation is the Weir–Cockerham estimator calculated using true
simulated groups; it was not supplied to fitting or used to select runs.

| Method | Mean ARI | Median ARI | Range |
|---|---:|---:|---|
| PopCluster | 0.1945 | 0.1813 | 0.0511–0.3456 |
| stockR | 0.0307 | 0.0308 | 0.0073–0.0666 |

PopCluster scored higher in all ten paired datasets. The mean paired difference
was 0.1638; an exploratory 95% percentile bootstrap interval was 0.0993–0.2292
(10,000 resamples of whole paired datasets, seed 21092026). With ten replicates
at one condition, this interval should not be interpreted as evidence of
superiority across settings. The extension followed inspection of the pilot.

The direction of the initial comparison persists, but neither method achieved
reliable group recovery in this weak setting. This separates partial structure
recovery from accurate clustering. It is compatible with PPR S greater than one
not guaranteeing accurate assignment and does not validate the full PPR rule.
Automatic K selection, independent prediction and uncertain genotypes remain
separate, untested endpoints in this extension.

All eight new program executions exited successfully; all 24 new stockR fits
had no warnings or iteration-limit exits. Input/output sample ordering, genotype
counts and heterozygous counts passed checks. Native fitted training likelihoods
matched likelihoods recomputed by the prediction adapter. Tests also checked
simulation reproducibility, distinct replicate frequencies, genotype encoding
and agreement with an independent adjusted Rand index implementation.

## Next stage

### PCA–k-means baseline and native initialization audit

Completed 21 September 2026 on the exact same 80 weak/strong datasets.
Principal component analysis (PCA) used centred dosages standardized by
sqrt(2p(1-p)), with frequencies estimated from fitting genotypes. Prespecified
component counts: 2, 5, 10, 20, 50 and 100. Each setting used three batches of
25 k-means starts, selected by minimum within-cluster sum of squares. All six
settings are reported, never selected using simulation truth. No missing values
were imputed because these panels contain complete genotypes.

Weak-signal mean adjusted Rand index (ARI), for the two endpoint component
counts (intermediate settings are included in the artifact report):

| Total samples | Nominal markers | PCA–k-means: 2 components | PCA–k-means: 100 components | Native stockR final | PopCluster |
|---:|---:|---:|---:|---:|---:|
| 120 | 1,000 | 0.189 | 0.142 | 0.192 | 0.188 |
| 120 | 10,000 | 0.211 | 0.085 | 0.116 | 0.166 |
| 360 | 1,000 | 0.147 | 0.067 | 0.065 | 0.119 |
| 360 | 10,000 | 0.192 | 0.039 | 0.043 | 0.209 |

Retaining fewer components helped mean recovery in these weak settings, but
recovery remains poor and this does not establish a universal component default.
Generating F = 0.02 yielded perfect PCA–k-means recovery in 39/40 datasets at
2, 5 and 10 components (one dataset ARI = 0.974893), and 40/40 at 20, 50 and 100.

Native stockR uses different preprocessing: uncentred, unscaled rotation and
up to 100 components. Exact `init.grps` were read from its saved fits, not
reconstructed. For weak datasets, 33/40 likelihood-selected native fits had
identical initial/final partitions; across all starts, 101/120 were identical.
The seven changed selected fits improved truth ARI in six cases and worsened it
in one. For stronger datasets, all 120 native starts were already perfect and
all remained unchanged. Partition identity was checked directly, not inferred
from equal truth ARI. Membership probabilities may change even when hard
classifications do not.

This motivates a controlled initialization experiment, not a claim that stockR
is defective or that changing component count alone explains every difference.
Standalone and native initializers differ in centring/scaling as well as count.
All 480 standalone settings completed; 1,440 batches had no warnings or failed
convergence flags. Matched sample/marker checks, direct PCA decomposition,
seed reproducibility and objective-selection checks passed.

Artifacts: `pca-kmeans-20260921`, including `report.md`, component-sensitivity
figure, all fit objects and exact native partition audit. Runner, tests and
summary scripts are named `*-pca-kmeans-benchmark.R`; the direct partition
audit is `audit-native-pca-partitions.R`. No published package default or
previous PDF attachment was changed.

### Fixed generating F = 0.02: stronger-signal comparison

Completed 21 September 2026. Forty new independent datasets (seeds
924101–924140), ten per condition, true K = 3 supplied, complete independent
markers and balanced true groups. Same native runner, three starts per method,
and identical genotypes within each method pair. Unlike the weak equal-score
experiment, generating differentiation stays fixed as sample/marker counts vary.
WC84 denotes the Weir–Cockerham (1984) estimator.

| Total samples | Per group | Nominal markers | Retained markers | Mean measured WC84 FST | stockR mean ARI | PopCluster mean ARI |
|---:|---:|---:|---:|---:|---:|---:|
| 120 | 40 | 1,000 | 1,000 | 0.019935 | 1 | 1 |
| 120 | 40 | 10,000 | 9,999–10,000 | 0.019970 | 1 | 1 |
| 360 | 120 | 1,000 | 1,000 | 0.020097 | 1 | 1 |
| 360 | 120 | 10,000 | 10,000 | 0.020074 | 1 | 1 |

Every selected fit in both
methods recovered groups perfectly: 40/40 datasets, with zero across-dataset
standard deviation of ARI in all four cells. Two 120-sample/10,000-marker panels
lost one monomorphic marker; counts were matched between methods. All 40
PopCluster executions and 120 native stockR fits completed without recorded
warnings or iteration-limit/ineligible starts under the runner's checks.
Frozen software/source hashes match the weak experiment, and selected
objectives passed audit. Neither truth labels nor evaluation scores selected runs.

This stronger-signal result puts the weak-signal results in context: they are
not a general failure of either program. It is not proof of a universal FST
threshold, correct automatic K selection, independent prediction, or performance
on low-coverage/missing/linked genotypes. No experimental assigner model was
evaluated here. All results are training recovery with true K supplied.

Raw new artifacts: `fixed-f002-20260921`. Combined weak/strong report, figure,
per-panel results and source audit: `stockr-popcluster-weak-strong-20260921`.
The two-page attachment for Scott is
`output/pdf/stockr-popcluster-benchmark-for-scott.pdf` in the shared workspace.

### Equal-information extension: ten replicates per condition

Completed 21 September 2026: 32 new datasets (seeds 922101–922132) plus the
eight preserved pilot datasets. Forty datasets in total, ten per condition;
K = 3 supplied, generating S = 3, complete independent markers and three
native runs per method throughout. Identical software/source hashes were
verified across the pilot and extension. The extension followed pilot inspection.
Standard deviations (SD) below describe independent-dataset variability.

| Total samples | Markers | PopCluster mean ARI (SD) | stockR mean ARI (SD) | Mean paired difference | Exploratory 95% interval | PopCluster higher |
|---:|---:|---:|---:|---:|---:|---:|
| 120 | 1,000 | 0.188 (0.105) | 0.192 (0.124) | -0.004 | -0.062 to 0.052 | 6/10 |
| 360 | 1,000 | 0.119 (0.059) | 0.065 (0.049) | 0.054 | 0.031 to 0.075 | 9/10 |
| 120 | 10,000 | 0.166 (0.152) | 0.116 (0.102) | 0.050 | -0.051 to 0.148 | 5/10 |
| 360 | 10,000 | 0.209 (0.090) | 0.043 (0.024) | 0.166 | 0.123 to 0.213 | 10/10 |

Paired difference is PopCluster minus stockR. Intervals are percentile bootstrap
intervals from 10,000 resamples of whole paired datasets within each cell;
they are exploratory, not simultaneous intervals or broad superiority claims.
The apparent stockR advantage in the 120-sample/1,000-marker pilot did not
persist. PopCluster's clearest advantage is at 360 samples/10,000 markers,
but neither method provides reliably accurate group recovery at this setting.
The 120-sample/10,000-marker comparison is particularly variable.

Measured S spans 1.947–3.896 despite generating S = 3. These comparisons jointly
change sample size, marker count and generating differentiation: they do not
isolate the causal effect of any one factor or validate a universal PPR threshold.
They do not test automatic K selection or independent prediction.

All 40 PopCluster executions succeeded. All 120 stockR starts had no recorded
warnings or iteration-limit/ineligible exits under the runner's checks. Choices
were audited against the recorded objectives; no dataset was dropped.
The complete report, per-dataset measured differentiation, paired figure and
machine-readable summary are in workspace artifacts
`ppr-equal-information-ten-replicates-20260921`. Source extension outputs remain
in `ppr-equal-information-extension-20260921`. Reproduce the combined report
with `summarize-ppr-equal-information-extension.R`.

### Equal-information pilot completed 21 September 2026

Eight fresh datasets crossed 120 or 360 training samples with 1,000 or 10,000
independent markers, two replicates per combination. True K = 3 was supplied;
each true group contained 40 or 120 fish. Generating F was adjusted to keep
S = N L F² = 3. All nominal markers survived the variability filter.
These compare native stockR and PopCluster, not the experimental regularized
mixture implementation. Three runs per method were selected by each method's
own training objective before truth-based evaluation.

| Total samples | Markers | Generating F | Mean measured FST | Mean PopCluster ARI | Mean stockR ARI |
|---:|---:|---:|---:|---:|---:|
| 120 | 1,000 | 0.0050000 | 0.0053780 | 0.2027 | 0.3149 |
| 360 | 1,000 | 0.0028868 | 0.0027806 | 0.1275 | 0.0758 |
| 120 | 10,000 | 0.0015811 | 0.0015358 | 0.3024 | 0.1586 |
| 360 | 10,000 | 0.0009129 | 0.0009700 | 0.2675 | 0.0315 |

Equal generating PPR scores did not yield equal recovery in this pilot, and
the method ranking varied by condition. This does not establish a sample-size
or marker-count effect: generating differentiation also changed by design.
Measured scores ranged from 2.083 to 3.896, and each condition has only two
replicates. No method achieved consistently accurate recovery. This is not a
test of selecting K, held-out prediction, or a universal PPR transition.

All eight PopCluster executions succeeded. All 24 stockR fits had finite eligible
results, no captured warnings and no iteration-limit exits under the existing
runner's checks. Input/output identifiers, genotype counts and recomputed
stockR likelihoods passed checks. Full per-panel results, measured differentiation,
frozen sources and inputs are in workspace artifacts
`ppr-equal-information-20260921`; `report.md` and `summary.rds` are the summaries.

Next expand independent replication before making power claims, then vary S
above and below the current setting. Keep native K selection and independent
prediction as separate experiments, including null datasets for K selection.

The independent weak-differentiation extension is complete. Next vary sample
size, marker count and differentiation in matched PPR comparisons; do not infer
a transition curve from this single setting. Test native K selection and
held-out prediction separately, including unstructured null datasets. The wider
PPR protocol is in ppr-popcluster-protocol.md.

### Common reference for sample size, marker count and differentiation

Continue using the Patterson–Price–Reich rule throughout these experiments:

\[
S = N L F_{\mathrm{ST}}^2.
\]

N is the total training sample size, not the sample size per population and not
training plus held-out fish. L is the number of independent retained markers.
Vary sample size, marker count and differentiation separately, and include
combinations with equal S to test whether the same theoretical information
produces comparable recovery. Keep K and other design features fixed within
each matched comparison; examine their effects separately.

Each benchmark should report total and per-population sample sizes, nominal
and retained marker counts, generating differentiation, and measured
Weir–Cockerham differentiation separately. Label S calculated from the generating
parameter separately from S calculated using measured differentiation. Measure
differentiation with simulation truth for evaluation, not inferred groups.

Report three distinct outcomes: group recovery with K supplied, automatic
selection of K, and held-out prediction. Include null datasets and uncertainty
across independent replicates. PPR is a reference, not a pass/fail guarantee:
linked markers, uncertain genotypes, structured missingness and unequal sample
sizes can change performance. Do not treat raw correlated marker counts as
independent L or invent a coverage-to-effective-marker correction.

## Reproducibility

Workspace artifacts are under:
`/Users/gos040/.codex/.chatgpt-projects/g-p-6a8255ffbc948191b3592f138fe06cf5/artifacts/ppr-popcluster-matched-20260921`.
They include `run.R`, `summarize.R`, `audit.rds`, `comparison.rds`, parameter
files, genotype inputs, evaluation-only mappings, complete program outputs,
console/error logs and execution records with binary MD5. Genotype export code
and encoding tests are in this experimental directory.

The extension is stored beside the initial artifacts in
`artifacts/popcluster-weak-replicates-20260921`, with frozen design/source records,
simulation truth, native fitted models, complete PopCluster outputs, locked
choices, measured differentiation and `summary.rds`. Reproduction scripts are
`run-popcluster-weak-replicates.R` and `summarize-popcluster-weak-replicates.R`;
checks are in `test-popcluster-weak-design.R`. The runner requires a new output
directory and does not overwrite existing runs.
