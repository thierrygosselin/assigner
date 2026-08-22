# Compare assigner and WGSassign unsampled-source scores

Join standardized unsampled-source scores produced by
\[assign_individuals()\] with independently generated WGSassign scores.
This is a validation aid; it does not call or reproduce WGSassign.

## Usage

``` r
compare_unsampled_scores(
  result,
  wgsassign,
  id.column = "INDIVIDUALS",
  score.column = "WGSASSIGN_Z"
)
```

## Arguments

- result:

  An object returned by \[assign_individuals()\] with \`unsampled.source
  = TRUE\`, or its \`assignment\` table.

- wgsassign:

  A data frame containing individual identifiers and WGSassign scores.

- id.column:

  Name of the WGSassign identifier column. Default:
  `id.column = "INDIVIDUALS"`.

- score.column:

  Name of the WGSassign score column. Default:
  `score.column = "WGSASSIGN_Z"`.

## Value

A list containing the joined scores, Pearson and Spearman correlations,
and a linear calibration model.
