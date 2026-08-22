# assigner_dots

Extract and assign the dots-dots-dots Capture and assign the \`...\`
arguments in a robust and tidy way. This utility function: - Assigns
known \`...\` arguments (called "keepers") into the function's
environment. - Detects and warns about unknown or deprecated
arguments. - Assigns default values for expected \`...\` arguments not
supplied. - Logs everything in a tidy tibble grouped by type.

Especially useful for deeply nested or extensible functions (e.g.,
pipelines or packages).

Includes a \`dev.mode\` flag to allow testing in the global environment
interactively.

## Usage

``` r
assigner_dots(
  func.name = as.list(sys.call())[[1]],
  fd = NULL,
  args.list = NULL,
  dotslist = NULL,
  keepers = c("subsample.markers.stats", "force.stats", "id.stats", "subsample",
    "filter.reproducibility", "filter.individuals.missing",
    "filter.individuals.heterozygosity", "filter.individuals.coverage.total",
    "filter.individuals.coverage.median", "filter.individuals.coverage.iqr",
    "filter.common.markers", "filter.monomorphic", "filter.mac", "filter.coverage", "dp",
    "filter.genotyping", "filter.snp.position.read", "filter.snp.number",
    "filter.short.ld", "filter.long.ld", "long.ld.missing", "ld.method", "ld.figures", 
 
       "detect.mixed.genomes", "ind.heterozygosity.threshold",
    "detect.duplicate.genomes", "filter.hwe", "filter.strands", "random.seed",
    "path.folder", "filename", "parameters", "blacklist.genotypes", "erase.genotypes",
    "pop.levels", "pop.labels", "pop.select", "blacklist.id", "markers.info",
    "keep.allele.names", "keep.gds", "vcf.metadata", "vcf.stats", "wide",
    "whitelist.markers", "write.tidy", "missing.memory", "dart.sequence", "internal",
    "tidy.check", "tidy.vcf", "tidy.dart", "gt", "alt.dosage", "gt.vcf", 
     "gt.vcf.nuc",
    "calibrate.alleles", "forking"),
  deprecated = NULL,
  verbose = TRUE,
  dev.mode = FALSE
)
```

## Arguments

- func.name:

  The name of the function calling \`assigner_dots\`. Default:
  `as.list(sys.call())[[1]]`.

- fd:

  (optional) A character vector of formal argument names for the calling
  function. Default:
  [`rlang::fn_fmls_names()`](https://rlang.r-lib.org/reference/fn_fmls.html).

- args.list:

  (optional) A named list of arguments already present in the calling
  function's environment. Default:`args.list = as.list(environment())`.

- dotslist:

  The captured \`...\` using \`rlang::dots_list(...)\`. Default:
  ` dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE)`.

- keepers:

  (optional) A character vector of allowed dot arguments that should be
  assigned. Default:
  `keepers = c( "subsample.markers.stats", "subsample", "filter.reproducibility", "filter.individuals.missing", "filter.individuals.heterozygosity", "filter.individuals.coverage.total", "filter.individuals.coverage.median", "filter.individuals.coverage.iqr", "filter.common.markers", "filter.monomorphic", "filter.mac", "filter.coverage", "filter.genotyping", "filter.snp.position.read", "filter.snp.number", "filter.short.ld", "filter.long.ld", "long.ld.missing", "ld.method", "detect.mixed.genomes", "ind.heterozygosity.threshold", "detect.duplicate.genomes", "filter.hwe", "filter.strands", "random.seed", "path.folder", "filename", "blacklist.genotypes", "erase.genotypes", "gt", "alt.dosage", "gt.vcf", "gt.vcf.nuc", "pop.levels", "pop.labels", "pop.select", "blacklist.id", "markers.info", "keep.allele.names", "keep.gds", "vcf.metadata", "vcf.stats", "id.stats", "dp", "whitelist.markers", "write.tidy", "dart.sequence", "missing.memory", "internal", "tidy.check", "tidy.vcf", "tidy.dart", "calibrate.alleles", "forking" )`.

- deprecated:

  (optional) Vector of deprecated argument names to detect in \`...\`.
  Default: `deprecated = NULL`.

- verbose:

  (logical) Function will output more details. Will print messages
  describing assignments and issues. Default: `verbose = TRUE`.

- dev.mode:

  (logical) If \`TRUE\`, arguments will be assigned to the global
  environment (\`globalenv()\`) for interactive development and testing.
  Default: `dev.mode = FALSE`.

## Value

A \`tibble\` with columns \`ARGUMENTS\`, \`VALUES\`, and \`GROUPS\`,
indicating how each dot-argument was handled: - \`"fct.call.args"\`:
named arguments from the function call - \`"fct.call..."\`: valid dot
arguments passed - \`"default..."\`: valid dot arguments not passed
(defaults assigned) - \`"deprecated..."\`: deprecated arguments passed -
\`"unknowned..."\`: unknown arguments passed
