# Scott Foster review package

Research prototype - methodological review requested. Permission to use code is
not endorsement or completed review. No public release is implied.

Open `review.html` for the main report; `review.pdf` is the fixed review copy.
Edit `review.qmd`. From this directory render with:

```sh
quarto render review.qmd --to html
quarto render review.qmd --to pdf
```

Rendering needs R, knitr, rmarkdown and ggplot2; PDF additionally needs XeLaTeX.
It reads the included compact RDS and never reruns simulations.
To refresh data from the original development artifacts:

```sh
Rscript prepare-review.R /absolute/path/to/artifacts
```

`source-snapshot/` contains the implementation and shared GDS helpers at review
preparation time. `data/source-checksums.rds` records their checksums. Historical
benchmark versions may differ; the original experiment audit files remain in the
development artifact directories. This compact package is not the full raw-data
reproduction archive. Source remains subject to the assigner repository license.
