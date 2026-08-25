# Simulated dataset to test assigner

A dataset of genotypes containing 500 bi-allelic SNPs, simulated for 250
individuals, 5 populations.

## Usage

``` r
data_assigner_sim_01
```

## Format

A tibble with 125000 rows (genotypes) and 4 variables:

- MARKERS:

  SNPs markers

- POP_ID:

  Populations/strata for the samples

- INDIVIDUALS:

  Samples id)

- GT:

  Genotypes coded *a la genepop* format

## Source

The data was simulated with grur
<https://thierrygosselin.github.io/grur/reference/simulate_rad.html>

## Details

Dataset simulation caracteristics:

- num.pops: 5

- num.loci: 1000

- div.time: 25e3

- ne: 200

- nm: 0.5

- theta: 0.2

- mig.type: island

- mut.rate: 2.5e-4

- mig.rate: 0.0025

From this simulated dataset, 500 SNPs and 250 individuals (50 ind/pop)
were randomly selected.
