# assigner development roadmap

This roadmap records methods that are scientifically useful but require a
careful data contract, validation, and documentation before they become part of
the public API. Inclusion here is not a promise of a release date.

## 1. Absolute divergence (Dxy)

Develop `dxy()` for called genotypes in reference-guided Genomic Data Structure
(GDS) files. Absolute divergence measures how often two allele copies, one from
each group, differ at a nucleotide position. Nucleotide diversity (π) measures
the corresponding differences within a group. Examine both alongside the
fixation index ($F_{ST}$), which can rise where within-group diversity falls.
[Nei and Li (1979)](https://doi.org/10.1073/pnas.76.10.5269),
[Nei (1987, equation 10.20)](https://books.google.com/books?id=UhRSsLkjxDgC),
[Cruickshank and Hahn (2014)](https://doi.org/10.1111/mec.12796).

At site $i$, let $n_{Xi}$ and $n_{Yi}$ be the numbers of called allele copies in
groups $X$ and $Y$, and $n_{Xia}$ and $n_{Yia}$ the counts of nucleotide $a$.
There are $n_{Xi}n_{Yi}$ valid cross-group allele pairs and
$n_{Xi}n_{Yi} - \sum_a n_{Xia}n_{Yia}$ mismatching pairs. At a biallelic site,
the mismatch probability is $p_X(1-p_Y) + (1-p_X)p_Y$, where $p_X$ and $p_Y$
are the alternate allele frequencies among called copies.

For a window $W$, pool the mismatch and comparison counts across sites.

$$
\widehat{d}_{XY}(W) =
\frac{\sum_{i \in W} \left(n_{Xi}n_{Yi} - \sum_a n_{Xia}n_{Yia}\right)}
     {\sum_{i \in W} n_{Xi}n_{Yi}}.
$$

Only sites with called copies in both groups contribute valid pairs. A callable
site fixed for the same nucleotide in both groups contributes
valid pairs and zero mismatches. Report zero when comparisons exist and all
match, and `NA` when no valid pairs remain.
[Korunes and Samuk (2021)](https://doi.org/10.1111/1755-0998.13326),
[2022 corrigendum](https://doi.org/10.1111/1755-0998.13571).

Report **Dxy per callable base** from variant and callable invariant records.
An equivalent input would pair a complete set of variant calls at eligible
sites with matching per-sample callability data that recover the cross-group
pair counts at invariant sites. For a single-nucleotide polymorphism
(SNP) panel alone, report **mean divergence per analysed marker**, weighted by
the number of valid pairs. This estimate depends on which markers were
retained. Report sample sizes, callable sites, variable sites, valid pairs,
mismatching pairs and the physical span for each group pair and window.

Start with fixed, non-overlapping physical windows within chromosomes or
linkage groups. Validate the estimator against hand-worked cases and pixy on
the same callable sites, including unequal group sizes and missing calls. Add
sliding windows, Gaussian smoothing, phased haplotypes and direct input of
aligned sequence blocks after the initial implementation.

## 2. AMOVA - initial implementation available

`amova_genomic()` now provides the first Analysis of Molecular Variance
framework with:

- GDS files or open GDS objects as the primary input;
- one or more hierarchical grouping columns from GDS metadata or strata data;
- SNP-dosage, identity, nucleotide, Manhattan, and custom distance models;
- locus-wise, filtered, and complete-case missing-data strategies;
- variance components and Phi-statistics at supported hierarchy levels;
- Euclidean validation and optional Lingoes correction;
- hierarchy-aware permutation tests with reproducible random seeds;
- locus and genomic-block bootstrap intervals and population jackknife
  sensitivity with a layered uncertainty report; and
- explicit reporting of retained loci, missing data, sample size, eigenvalues,
  corrections, and degrees of freedom.

AMOVA is a separate public function rather than an option hidden inside
`fst_WC84()`. Remaining work is validation against published examples,
population and full hierarchical bootstraps, plus experimental parametric and
wild-cluster bootstrap methods,
simulations, and established programs including `ade4`, `pegas`, `poppr`,
Arlequin, GenoDive, and Stacks; deeper testing of informative missingness; and
profiling before moving selected numerical kernels to C++.

## 3. Phased haplotypes

Coordinate the work across packages:

- genometranslator: import, validate, store, and export phased haplotypes,
  including STACKS 2 haplotypes and phased VCF genotypes;
- radr: haplotype-level quality control and filtering;
- assigner: population statistics and assignment models using validated
  haplotypes.

Candidate haplotype statistics include:

- haplotype counts and frequencies;
- haplotype and gene diversity;
- haplotype-level FST;
- PhiST and standardized F'ST;
- haplotype-level Dxy; and
- assignment likelihoods based on haplotype frequencies.

The representation must distinguish unresolved phase, partially observed
haplotypes, missing haplotypes, and true biological absence. STACKS 2 can shorten
a locus to obtain complete haplotypes at a requested sample threshold; assigner
should never perform that filtering silently.

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
documentation should encourage their joint interpretation rather than present
them as interchangeable measures.

## Foundational references

- Excoffier L, Smouse PE, Quattro JM (1992). Analysis of molecular variance
  inferred from metric distances among DNA haplotypes. Genetics, 131, 479-491.
- Nei M (1987). Molecular Evolutionary Genetics. Columbia University Press.
- Meirmans PG (2006). Using the AMOVA framework to estimate a standardized
  genetic differentiation measure. Evolution, 60, 2399-2402.
- Cruickshank TE, Hahn MW (2014). Reanalysis suggests that genomic islands of
  speciation are due to reduced diversity, not reduced gene flow. Molecular
  Ecology, 23, 3133-3157.
- Davey JW, Cezard T, Fuentes-Utrilla P, Eland C, Gharbi K, Blaxter ML (2013).
  Special features of RAD sequencing data: implications for genotyping.
  Molecular Ecology, 22, 3151-3164.
- Rochette NC, Rivera-Colon AG, Catchen JM (2019). Stacks 2: Analytical methods
  for paired-end sequencing improve RADseq-based population genomics. Molecular
  Ecology, 28, 4737-4754.
