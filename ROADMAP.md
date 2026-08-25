# assigner development roadmap

This roadmap records methods that are scientifically useful but require
a careful data contract, validation, and documentation before they
become part of the public API. Inclusion here is not a promise of a
release date.

## 1. Absolute divergence (Dxy)

Develop pairwise Dxy for:

- called SNP genotypes;
- phased haplotypes;
- complete RAD loci or other aligned sequence blocks; and
- fixed, sliding, and Gaussian-kernel genomic summaries.

Important design requirements:

- distinguish absolute divergence from relative differentiation such as
  FST;
- define how missing genotypes and unequal sample sizes contribute to
  each population pair;
- retain callable invariant sites when Dxy is reported per base pair;
- clearly label SNP-only estimates that cannot use invariant sites;
- return the number of individuals, chromosomes, markers, callable
  sites, and sequence length contributing to every estimate; and
- validate results against established software and simulated datasets.

The first implementation should prioritize pairwise Dxy from
reference-guided GDS data. Haplotype-level Dxy should follow once phased
haplotypes have a stable representation.

## 2. AMOVA - initial implementation available

[`amova_genomic()`](https://thierrygosselin.github.io/assigner/reference/amova_genomic.md)
now provides the first Analysis of Molecular Variance framework with:

- GDS files or open GDS objects as the primary input;
- one or more hierarchical grouping columns from GDS metadata or strata
  data;
- SNP-dosage, identity, nucleotide, Manhattan, and custom distance
  models;
- locus-wise, filtered, and complete-case missing-data strategies;
- variance components and Phi-statistics at supported hierarchy levels;
- Euclidean validation and optional Lingoes correction;
- hierarchy-aware permutation tests with reproducible random seeds;
- locus and genomic-block bootstrap intervals and population jackknife
  sensitivity with a layered uncertainty report; and
- explicit reporting of retained loci, missing data, sample size,
  eigenvalues, corrections, and degrees of freedom.

AMOVA is a separate public function rather than an option hidden inside
[`fst_WC84()`](https://thierrygosselin.github.io/assigner/reference/fst_WC84.md).
Remaining work is validation against published examples, population and
full hierarchical bootstraps, plus experimental parametric and
wild-cluster bootstrap methods, simulations, and established programs
including `ade4`, `pegas`, `poppr`, Arlequin, GenoDive, and Stacks;
deeper testing of informative missingness; and profiling before moving
selected numerical kernels to C++.

## 3. Phased haplotypes

Coordinate the work across packages:

- genometranslator: import, validate, store, and export phased
  haplotypes, including STACKS 2 haplotypes and phased VCF genotypes;
- radr: haplotype-level quality control and filtering;
- assigner: population statistics and assignment models using validated
  haplotypes.

Candidate haplotype statistics include:

- haplotype counts and frequencies;
- haplotype and gene diversity;
- haplotype-level FST;
- PhiST and standardized F’ST;
- haplotype-level Dxy; and
- assignment likelihoods based on haplotype frequencies.

The representation must distinguish unresolved phase, partially observed
haplotypes, missing haplotypes, and true biological absence. STACKS 2
can shorten a locus to obtain complete haplotypes at a requested sample
threshold; assigner should never perform that filtering silently.

## Validation and interpretation

These methods will require tests covering:

- uneven population sample sizes;
- rare alleles and unsampled haplotypes;
- missingness and low coverage;
- linkage disequilibrium;
- genotyping and phasing errors;
- restriction-site polymorphism and allele dropout in RAD data;
- variable marker density across the genome; and
- reference-guided and de novo datasets.

Dxy, FST, PhiST, and related statistics answer different questions. The
documentation should encourage their joint interpretation rather than
present them as interchangeable measures.

## Foundational references

- Excoffier L, Smouse PE, Quattro JM (1992). Analysis of molecular
  variance inferred from metric distances among DNA haplotypes.
  Genetics, 131, 479-491.
- Nei M (1987). Molecular Evolutionary Genetics. Columbia University
  Press.
- Meirmans PG (2006). Using the AMOVA framework to estimate a
  standardized genetic differentiation measure. Evolution, 60,
  2399-2402.
- Cruickshank TE, Hahn MW (2014). Reanalysis suggests that genomic
  islands of speciation are due to reduced diversity, not reduced gene
  flow. Molecular Ecology, 23, 3133-3157.
- Davey JW, Cezard T, Fuentes-Utrilla P, Eland C, Gharbi K, Blaxter ML
  (2013). Special features of RAD sequencing data: implications for
  genotyping. Molecular Ecology, 22, 3151-3164.
- Rochette NC, Rivera-Colon AG, Catchen JM (2019). Stacks 2: Analytical
  methods for paired-end sequencing improve RADseq-based population
  genomics. Molecular Ecology, 28, 4737-4754.
