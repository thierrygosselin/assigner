# Experimental GL/PL stock mixture: review record

This internal prototype is not a public analysis function. Scott Foster's
statistical and code review is a prerequisite before public release of this
extension. His permission to use stockR code has been communicated by Thierry;
permission is not scientific endorsement or completed review.

## Attribution and source

Foster SD, Feutry P, Grewe PM, Berry O, Hui FKC, Davies CR (2018).
Reliably discriminating stock structure with genetic markers: Mixture models
with robust and fast computation. Molecular Ecology Resources 18:1310–1325.
https://doi.org/10.1111/1755-0998.12920

Inspected local stockR 1.0.76, R/stockR.R: calcCondProbs, calcCondMeans,
calcPis and calcTaus. Copyright attribution: Scott D. Foster; GPL (>= 2),
compatible with assigner's GPL-3. The new R code independently implements
the likelihood extension; it is not a verbatim copy of these functions.

## Model and deliberate differences

One discrete stock per individual; diploid biallelic genotypes follow
Hardy–Weinberg proportions within stocks; markers are treated as independent.
For each stock and marker, sum genotype likelihood times stock-specific
genotype probability over the three possible genotypes. Multiply over markers,
then mix across stocks. Membership probabilities are not ancestry proportions.

EM updates frequencies with posterior ALT counts conditional on both stock
membership and sequencing evidence. GL and PL are not genotype posteriors.
The engine normalizes each likelihood vector by its maximum: reported log
likelihoods omit these data-only constants and should only be compared using
the same observations and normalization.

Flat and wholly missing vectors contribute exactly zero log likelihood and are
excluded from frequency updates. Entirely uninformative individuals do not
estimate mixture proportions; their membership equals the fitted prior.
Markers with no information are rejected. Partially missing vectors are rejected.

stockR 1.0.76 floors frequencies at 1/(10*nFish) in calcCondProbs, but not
in calcCondMeans. The prototype instead uses a single explicitly bounded
frequency parameter in both steps. Exact comparisons must respect that
difference, especially for private/fixed alleles. The bound is numerical, not
a biological assertion that every stock carries every allele.

## Completed first stage

### Subsequent larger-panel proof of concept

See [fresh-replicate proof of concept](stock-proof-fresh.md) and
[native stockR comparison](stockr-native-comparison.md). These distinguish the
matched-engine experiments from native stockR Kmeans.EM/AIC/BIC and retain the
missed-K cases. Fresh fixed-strength results support further development at
FST near 0.005 without implying universal superiority or uncertain-GL validation.
These records should accompany Scott's eventual review.

The [sample-size/strength follow-up](strength-size-experiment.md) holds marker
count and initialization fixed. It documents improved recovery with larger
training cohorts, tuning's distinction between predictive density and hard
classification, and a separately labelled revisit of the earlier missed K.
The hard-call-limit accelerator is experimental and explicitly rejects uncertain
or missing GLs; it is not a substitute for the general GDS engine.

- Small in-memory reference EM engine with explicit starting values.
- Blockwise GDS engine, sequential multiple starts and per-start diagnostics.
- Stable log likelihoods and likelihood-decrease checks.
- Tests for hard-genotype limit, GL/PL equivalence, normalization, missing
  information and malformed inputs.
- Internal GDS adapter `.sm_fit_gds()` for existing paths or open SeqArray
  handles. It preserves active filters and aligns starting frequencies by
  variant.id and metadata by INDIVIDUALS. STRATA never enters estimation.
- Real VCF-to-GDS tests cover GL/PL agreement, missing floating-point storage,
  active filters, errors, metadata and borrowed handle ownership.

## Using the development GDS adapter

This is not exported. After loading the development source with devtools,
an explicit starting-frequency matrix can be supplied as follows:

```r
gds <- SeqArray::seqOpen("likelihoods.gds", readonly = TRUE)
# Apply any intended sample/variant filters before obtaining these IDs.
variant.ids <- SeqArray::seqGetData(gds, "variant.id")
# K=1 smoke test only: not a procedure for choosing the number of stocks.
start <- matrix(0.5, length(variant.ids), 1,
                dimnames = list(as.character(variant.ids), "STOCK_1"))
fit <- assigner:::.sm_fit_gds(gds, frequency = start, likelihood = "auto",
                            n.starts = 1, block.snps = 1000, seed = 42)
SeqArray::seqClose(gds)
fit$assignments
fit$strata
```

The auto setting prefers PL, then GL, without falling back to called genotypes
or posterior GP. The caller must verify diploid biallelic calling assumptions.
For K>1, supply distinct, defensible starting columns. The first start uses
these columns and supplied proportions. Additional starts use independent
random soft partitions of samples to initialize posterior allele counts, with
the mean supplied frequency as the initial genotype prior. These are not
stockR's K-means starts, hard genotype calls or fixed memberships. The default
is five starts; tune n.starts and assess likelihood repeatability for your data.

The GDS reader streams up to block.snps markers at a time (default 1000).
The max.genotypes guard now bounds sample-marker cells PER BLOCK, shrinking
blocks automatically rather than rejecting a large total marker count.
Frequency matrices (markers x stocks) and membership matrices (samples x
stocks) remain in RAM. Likelihood storage scales with samples x block.snps;
there is no full samples x markers x stocks posterior array, nor a full
sample-marker missingness matrix in the returned fit. Counts are returned as
sample.information and marker.information instead. This is not an absolute
RAM guarantee. Additional transient R copies, GDS buffers, and identifiers
also consume memory. Starts run sequentially and only the best fit is retained.

Each EM iteration makes two passes through the likelihood data. Membership
is calculated from ALL selected markers; the second pass updates frequencies
with those full-data responsibilities. This is not per-block clustering or
stochastic EM. Reduced memory costs additional I/O, especially on network
storage. Full-scale performance benchmarking is still required.

The selected fit is the highest-likelihood CONVERGED start. If no start
converges, return the highest-likelihood unfinished fit and warn. Failed and
unfinished starts remain visible in start.diagnostics; all failures stop the
call. LOGLIK_MINUS_SELECTED may be positive when an unfinished run exceeds
the selected converged run: increase iterations and investigate, rather than
ignoring that evidence. Component labels are arbitrary and are NOT aligned
across starts; equal likelihood alone does not establish identical partitions.
The seed, selected start and realized block size are returned. RNG state and
kind are restored even on error. Seed reproducibility applies within a fixed
software/data environment, not as a promise of cross-platform bitwise identity.
All-uninformative markers are rejected, not silently dropped. Samples with
no informative markers retain prior probabilities and an NA assigned stock.
Nonconvergence is warned about. Neither convergence nor high membership
probability proves a valid stock or a global optimum.

## Required before a public interface

- Realistic full-scale performance measurements; optional more informed starts.
- Cross-start partition agreement, beyond likelihoods (bootstrap refits now
  have exact frequency-based label alignment).
- Improved component-collapse recovery (currently failed starts are reported).
- Further stockR comparisons beyond the bundled hard-genotype example, which
  now passes a full 100-individual, 5000-marker matched-initialization test.
- Realistic simulations: weak differentiation, private alleles, uneven depth,
  unequal samples, relatedness, linkage, technical confounding and admixture.
- Calibration of the new fixed-K individual bootstrap and held-out K comparison;
  marker-block and non-independent-sampling extensions remain outstanding.
- Decide whether to implement explicit same-stock sampling-group constraints;
  never silently turn STRATA metadata into such constraints.
- Document that genetic clusters do not establish demographic independence.
- Scott's review of equations, updates, numerical conventions and validation.

No changes to assign_likelihood_dapc are part of this prototype.

See [stock-mixture-validation.md](stock-mixture-validation.md) for the current
benchmark results, equations, resampling interfaces and important limitations.

## Blockwise/multistart verification

The initial three focused test files passed 77 expectations, with no unexpected
warnings or failures. They cover dense-versus-streamed EM parity, different
block sizes (including a final short block), equivalent GL/PL files, active
filter restoration, seeded reproducibility, caller RNG preservation, retention
of failed starts, validation and explicit nonconvergence warnings. This is
numerical/interface verification, not evidence of biological calibration or
whole-genome performance. No package installation or public export is implied.

The subsequent uncertainty/K-comparison stage now passes 104 expectations
across four test files. The full bundled stockR comparison passes, but the
held-out K criterion fails to recover the generating K=3 on a 300-marker
subset and is sensitive to the frequency bound. See stock-mixture-validation.md;
automatic biological stock-count inference is explicitly not validated.
# Mixed initialization and regularization extension (2026-09-16)

The internal GDS fit now accepts opt-in mixed starts and symmetric or
training-pooled Beta MAP regularization. The default remains random starts
without regularization for a directly comparable baseline. See the vignette
and function help for controls, memory limits and objective definitions.

Scott's review should distinguish three questions: optimization of a fixed
objective; interpretation of a regularized estimator; and selection of K by
prediction. Neither the changed objective nor better simulation performance
establishes demographic stock independence. The implementation is not a copy
of GenoDive: its documentation motivated diverse repeated starts and explicit
missingness diagnostics, but no GenoDive code or random-filled genotype data
were incorporated. Simulated annealing is not implemented.

Points for methodological review:

- Is Beta(1+s*q, 1+s*(1-q)) MAP regularization appropriate for the intended use?
- Should q remain a training-pooled plug-in, or be part of a hierarchical model?
- How should prior strength be selected without optimistic held-out evaluation?
- Which private-allele and weak-differentiation regimes are harmed by shrinkage?
- Should frequency uncertainty be integrated in prediction instead of using MAP?
- What sampling unit and model checking are required for bootstrap coverage?
- What evidence should accompany a preferred predictive K, beyond convergence?

No automatic prior tuning, integrated Bayesian prediction, global optimum
certificate or interval-coverage claim is supplied by this extension.
