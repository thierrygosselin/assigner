# Legacy name for assignment evaluation

\`assignment_ngs()\` is the former name of \[evaluate_assignment()\].
New code should use the more descriptive name.

## Usage

``` r
assignment_ngs(
  data,
  strata = NULL,
  pop.levels = NULL,
  assignment.analysis = c("gsi_sim", "adegenet"),
  markers.sampling = c("ranked", "random"),
  marker.number = "all",
  thl = 1,
  iteration.method = 10,
  subsample = NULL,
  iteration.subsample = 1,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1,
  ...
)
```

## Arguments

- data:

  Several input format are accepted. assigner uses genometranslator
  [`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.html)
  to import the data. See function documentation for more details.

- strata:

  Optional strata data or filename passed to readers that support it.
  Default: `strata = NULL`.

- pop.levels:

  (optional, string) This refers to the levels in a factor. In this
  case, the id of the pop. Use this argument to have the pop ordered
  your way instead of the default alphabetical or numerical order. e.g.
  `pop.levels = c("QUE", "ONT", "ALB")` instead of the default
  `pop.levels = c("ALB", "ONT", "QUE")`. White spaces in population
  names are replaced by underscore. Default: `pop.levels = NULL`.

- assignment.analysis:

  (character) Assignment analysis conducted with
  `assignment.analysis = "gsi_sim"` or
  `assignment.analysis = "adegenet"`. See **Details** section below for
  installing [gsi_sim](https://github.com/eriqande/gsi_sim).

- markers.sampling:

  (character) 2 options for markers selection:

  1.  `markers.sampling == "random"` the subset of markers are selected
      randomly, this is the classic *Leave-One-Out* (LOO) assignment.

  2.  `markers.sampling == "ranked"` the subset of markers are first
      ranked based on an overall *decreasing* Fst values. The Fst is
      computed with
      [`fst_WC84`](https://thierrygosselin.github.io/assigner/reference/fst_WC84.md)
      function, that uses a fast implementation of Weir and Cockerham
      1984 Fst/Theta equations. This selection method is used during
      *Training-Holdout-Leave One Out* (thl) assignment. How many
      markers are selected is controlled with argument `thl`.

- marker.number:

  (Integer or string of number or "all") The assignment analysis can use
  all your markers (default) or a subset of your markers. e.g. To test
  500, 1000, 2000 and all the markers:
  `marker.number = c(500, 1000, 2000, "all")`. To use only 500 makers
  `marker.number = 500`. How those markers are sampled is determined
  with the argument `markers.sampling`, next. Default =
  `marker.number = "all"`.

- thl:

  (character, integer, proportion) For `markers.sampling = "ranked"`
  only. Several options are available:

  1.  `thl = 1`: Only 1 individual sample is used as holdout. This
      individual is not participating in the markers ranking. For each
      marker number, the analysis will be repeated with all the
      indiviuals in the data set (e.g. 500 individuals, 500 times 500,
      1000, 2000 markers). This is the default.

  2.  `proportion`: e.g. `thl = 0.15`, 15 percent of the individuals, in
      each strata, are chosen randomly as holdout individuals.

  3.  `thl = "all"`: all individuals are used for ranking (not good) and
      the argument `iteration.method = 1` is set by default. This is for
      testing only.

  Different lists of holdout individuals are generated when the argument
  `iteration.method` (bootstrap) is used.

- iteration.method:

  (integer) With **random markers selection** method, the iterations
  argument = the number of iterations to repeat marker resampling.
  Default: `iteration.method = 10`.

  With `marker.number = c(500, 1000)` and default iterations setting,
  500 markers will be randomly chosen 10 times and 1000 markers will be
  randomly chosen 10 times.

  **Notes:** If all the markers are used, with `marker.number = "all"`
  or in a series of marker number groupings
  `marker.number = c(200, 500, "all")`, the number of iteration is
  automatically set to 1. The remaining groupings are unaffected.

  With **ranked makers selection** method, using `thl = 1`, the analysis
  will be repeated for each individuals in the data set for every
  `marker.number` selected. With a proportion argument `thl = 0.15`, 15
  percent of individuals in each populations are chosen randomly as
  holdout individuals and this process is reapeated the number of times
  chosen by the `iteration.method` value.

- subsample:

  (Integer or Character, optional) This argument subsample individuals.
  With `subsample = 36`, 36 individuals in each populations are chosen
  randomly to represent the dataset. This integer as to be smaller than
  the population with min sample size, if higher, the minimum sample
  size found in the data will be used as default. In doubt, use
  `subsample = "min"`, the function will use the smallest population
  sample size found in the data. The number of times this process is
  repeated is controlled by the argument `iteration.subsample`. Default:
  `subsample = NULL` (no subsampling).

- iteration.subsample:

  (Integer) The number of iterations to repeat subsampling of
  individuals. With `subsample = 20` and `iteration.subsample = 10`, 20
  individuals/populations will be randomly chosen 10 times. Default:
  `iteration.subsample = 1`.

- verbose:

  (optional, logical) `verbose = TRUE` to be chatty during execution.
  Default: `verbose = FALSE`.

- parallel.core:

  Number of processor cores passed to readers that support parallel
  processing. Default: `parallel.core = parallel::detectCores() - 1`.

- ...:

  Additional arguments passed to \[evaluate_assignment()\].

## Value

The value returned by \[evaluate_assignment()\].
