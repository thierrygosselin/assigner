# Patterson–Price–Reich rule and matched PopCluster benchmark

The Patterson–Price–Reich (PPR) rule is the common reference for experiments
varying sample size, marker count and population differentiation (FST).
K denotes the number of populations.

## Frozen initial comparison

Reuse the existing independent single-nucleotide polymorphism (SNP),
complete-genotype benchmark. Never generate
a different dataset for each method. Retain training-variable markers using the
same rule as `run-stockr-native.R`; preserve marker order and keep test samples
out of fitting. The first PopCluster pilot comprises six panels: generating
drift 0.001, 0.005 and 0.01, two independent replicates, 300 training fish and
10,000 nominal markers. Use mixture-only mode, true K = 3, three runs,
no frequency scaling and the default pooled-frequency prior. Population truth
must not be supplied to initialization. This is a small recovery pilot, not
an estimate of power or a native K-selection experiment.

PopCluster's hard-partition objective differs from stockR's marginal mixture
likelihood. Compare recovery and runtime, not raw likelihoods across programs.
Select starts by the program's training objective, never truth or test accuracy.
Use the adjusted Rand index (ARI), a chance-adjusted measure of agreement with
true groups, to compare training recovery with training recovery. Independent prediction is a separate
endpoint and must not be silently compared with training recovery.

## Three endpoints

1. Known-K recovery: ARI, group sizes, start-to-start variability and failures.
2. K selection: native stockR Akaike information criterion (AIC), Bayesian
   information criterion (BIC), predictive cross-validation (CV) and audited native
   PopCluster criteria, each identified separately. Include K = 1 simulations;
   inspect whether criteria can select K = 1 and boundary K before comparison.
   PopCluster manual recommends at least ten runs per K for this stage.
3. Held-out prediction: fit frequencies and tune settings using training data
   only. Evaluate independent fish. Any custom PopCluster prediction extension
   must be identified as such, not described as native software output.

## PPR experiment

Use S = N L FST² as a reference axis: N is the total training sample size,
excluding held-out individuals; L is the independent retained marker count.
Vary N, L and differentiation separately, including combinations with equal S.
Keep K and other conditions fixed within each matched comparison. Report total
and per-population sample sizes and nominal and retained marker counts.
Report generating F separately from measured Weir–Cockerham (1984) FST evaluated using simulation
truth. Never estimate F from inferred groups and then claim independent evidence
of detectability. Existing data allow a descriptive overlay, not validation:
two replicates per condition are insufficient to estimate transition curves.
Distinguish S calculated from generating F from S calculated from measured FST.
PPR is a reference, not a guarantee of recovery, K selection or prediction.

Next simulation stage: true K = 2, 3, 4; training N = 120, 240, 360;
L = 1000, 3000, 5000, 10000; target S = 0.25, 1, 4, 16, 64, with generating
F = sqrt(S/(N*L)). Start with a costed subset rather than launching the whole
factorial design. Add independent null datasets, report replicate uncertainty,
and use independent seeds. Balanced complete genotypes come first. Subsequently vary
imbalance, linked markers, missingness and coverage for genotype likelihoods (GL)
or their phred-scaled representation (PL), one factor at a time.
Do not substitute a guessed coverage-based effective L into the PPR formula.

## Reading notes

### Equal-information pilot launched 21 September 2026

Before expanding the full factorial design, compare true K = 3 at generating
S = 3, with N = 120 or 360 and L = 1,000 or 10,000, two fresh datasets per
combination. F = sqrt(3/(N L)); seeds 921101–921108. The four combinations
therefore vary the components of the same generating information score.
Use three native runs per method and the same training genotypes. This is a
known-K native stockR/PopCluster pilot, not a test of selecting K or a comparison
with the experimental regularized model. Measure differentiation using true
groups only after run selection. Two replicates are insufficient for a power
curve. Frozen design and outputs are in workspace artifacts
`ppr-equal-information-20260921`.

`prepare-ppr-equal-information.R` creates the design, the optional second
argument of `run-popcluster-weak-replicates.R` accepts that design, and
`summarize-ppr-equal-information.R` checks completeness and generates the report.
The original runner without that argument retains its original design/seeds.

Replication extension: add replicates 3–10 for each of the same four cells,
using independent seeds 922101–922132 (32 new datasets, 40 in total). Keep
K = 3 supplied, three native runs per method, identical genotype inputs within
each pair, and unchanged objectives/settings. Preserve pilot artifacts rather
than rerunning them. This extension was chosen after inspecting the pilot.
`prepare-ppr-equal-information-extension.R` checks identifier/seed uniqueness
against that pilot. Outputs are in `ppr-equal-information-extension-20260921`.
The combined summary will check software/source hashes and selected objectives,
report standard deviations and ranges, and use 10,000 whole-pair bootstrap
resamples per cell for exploratory percentile intervals of mean paired ARI
differences (seeds 923001–923004). No simultaneous or general superiority claim
is intended; no K-selection or held-out endpoint is added by this extension.
Completed: all 32 new datasets passed, bringing the total to 40. Combined
report and figure: `ppr-equal-information-ten-replicates-20260921`. See
`popcluster-known-k-results.md` for the updated interpretation; preserve the
two-replicate pilot as historical rather than the current summary.

### Fixed generating F = 0.02 comparison

Authorized after the weak-signal extension: K = 3 supplied, N = 120 or 360,
L = 1,000 or 10,000, ten independent datasets per cell, fixed generating
F = 0.02 (seeds 924101–924140). Use the unchanged native runner and three
starts per method, identical genotypes within each pair. All groups balanced;
no missingness or linkage. Measured Weir–Cockerham differentiation remains an
evaluation output, not an exact imposed target. The generating PPR scores are
48, 144, 480 and 1,440, rather than held constant at 3. The earlier weak panels
are separate datasets, not paired with this new regime. This is a stronger-signal
context for communication with Scott, still not K selection or a benchmark of
the experimental genotype-likelihood model. Inputs and outputs are retained in
`fixed-f002-20260921`; design script `prepare-fixed-f002.R`, comparison report
script `summarize-fixed-f002.R`. No earlier results are overwritten.
Completed: 40/40 datasets, perfect recovery for both methods in every selected
fit. Mean measured differentiation ranges from 0.019935 to 0.020097 across
cells. Two 120-sample panels retained 9,999 rather than 10,000 markers; the
nominal PPR score of 480 is consequently 479.952 in those panels. Combined
weak/strong outputs: `stockr-popcluster-weak-strong-20260921`.

### Agreed diagnostic policy after the Privé paper review

Principal component analysis (PCA) diagnostics remain plot-first: inspect
sample coordinates, genome-wide marker contributions and technical metadata
before choosing thresholds. Flagging an outlier is not authorization to remove
it. Long-range linkage disequilibrium (LD) exclusions for population-structure
analysis must be recorded separately from retained regions for inversion work.

Balanced sampling sites or STRATA do not establish balanced true genetic
groups, or the true number of groups. Compare balanced and unbalanced fitting
subsets as sensitivity analyses without treating supplied strata as genetic truth.

For related pairs within a stratum, the preferred starting strategy is to retain
one representative, generally the individual with less missing data. For pairs
spanning sampling sites or strata, compare exclusion of both with reintroduction
of one representative and further investigation. These are reviewed sensitivity
analyses, not automatic deletion rules; larger relatedness networks need more
than independent pair-by-pair decisions. Preserve identifiers, exclusion reasons,
metadata and before/after results. Projection methods requiring unrelated
training and projected individuals need an additional relatedness check even
when only one member of a pair was retained for fitting.

Wang (2025), doi:10.1016/j.tpb.2025.03.001, evaluates accuracy with true K
supplied. A transition near S = 1 is not a guarantee of correct K or accurate
individual assignment. Wang (2022), doi:10.1038/s41437-022-00535-z, uses
L = 1000/FST in its weak-differentiation experiment: one million SNPs at 0.001.
Privé, Luu, Blum, McGrath and Vilhjálmsson (2020), doi:10.1093/bioinformatics/btaa520,
motivate linkage disequilibrium (LD), relatedness, imbalance and principal
component analysis (PCA) projection diagnostics. Their
supplement is not a separate benchmark paper. Low-missingness human PCA results
do not validate mean imputation for structured restriction-site-associated DNA
sequencing (RADseq) missingness.

## Platform check — resolved 21 September 2026

The supplied PopClusterMac binary is Intel x86_64; this Mac is arm64.
Initial attempts stalled or terminated without results. Rosetta was verified
working; macOS then displayed an explicit security block. After the user approved
the executable, the bundled ant example completed with exit status 0 and wrote
`ant377NoScale_K_10_R_1` and `ant377NoScale.K`. Reported version: 1.7.0.0,
build 20250715. Run the approved original binary through `arch -x86_64` from
the isolated example directory. No global security settings were disabled and
no email was sent. This confirms execution, not yet scientific benchmark validity.
