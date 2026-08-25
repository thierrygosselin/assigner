# Import and Standardize Genomic Data

Import and format genomic data into a tidy format compatible with
downstream analyses in genometranslator. This function detects the input
data type (tidy data frame, GDS file, or SeqVarGDSClass) and performs
the appropriate transformation.

## Usage

``` r
import_data(data, calibrate.alleles = FALSE, verbose = FALSE)
```

## Arguments

- data:

  Genomic dataset. Can be a tibble/data frame in wide format, a GDS file
  path, or a `SeqVarGDSClass` object.

- calibrate.alleles:

  Logical. If `TRUE`, force recalibration of REF/ALT allele genotypes
  using
  [`calibrate_alleles`](https://thierrygosselin.github.io/genometranslator/reference/calibrate_alleles.html).
  Default: `FALSE`.

- verbose:

  Logical. Show detailed messages during import. Default: `FALSE`.

## Value

A tidy data frame (tibble) with standardized genotype format and
required metadata (e.g., GT column).

## See also

[`detect_genomic_format`](https://thierrygosselin.github.io/genometranslator/reference/detect_genomic_format.html),
[`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html),
and
[`calibrate_alleles`](https://thierrygosselin.github.io/genometranslator/reference/calibrate_alleles.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# From tidy data
tidy_data <- import_data(data = my_tibble)

# From GDS file
tidy_data <- import_data(data = "myfile.gds")

# With calibration
tidy_data <- import_data(data = my_tibble, calibrate.alleles = TRUE)
} # }
```
