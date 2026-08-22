## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Prepare working environment, eval = FALSE--------------------------------
# rm(list = ls())
# if (!require("pak")) install.packages("pak")
# pak::pkg_install("thierrygosselin/assigner")
# library(assigner)
# assigner::install_gsi_sim(fromSource = TRUE)

## ----get whitelists, eval = FALSE---------------------------------------------
# all.whitelists <- list.files(path = "~/Desktop/whitelists_salmon", pattern = "whitelist")
# all.whitelists # to see if the whitelists are all there

## ----tailored assignment function, eval = FALSE-------------------------------
# whitelists_assigner <- function(x) {
#   res <- assigner::assignment_ngs(
#     data = "subset_melville_salmon.tsv", #change for what you want
#     whitelist.markers = x, # don't change this one,
#     assignment.analysis = "adegenet", #change for what you want
#     common.markers = TRUE,#change for what you want
#     marker.number = "all",# I think you get the idea...
#     pop.levels = c(1,2,3,4,5,6,7,8,9,10,11),
#     sampling.method = "random",
#     iteration.method = 1,
#     keep.gsi.files = FALSE
#   )
# }

## ----run the assignment, eval = FALSE-----------------------------------------
# salmon.assignment.res <- purrr::map(
# .x = all.whitelists, # your whitelist created above
# .f = whitelists_assigner # your new function
# )

