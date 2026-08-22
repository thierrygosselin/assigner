## ----echo = FALSE, message = FALSE--------------------------------------------
    knitr::opts_chunk$set(collapse = T, comment = "#>")
    options(tibble.print_min = 4L, tibble.print_max = 4L)

## ----eval= FALSE--------------------------------------------------------------
# rm(list=ls())

## ----eval=FALSE---------------------------------------------------------------
# if (!require("pak")) install.packages("pak")
# pak::pkg_install("thierrygosselin/assigner")
# library(assigner)

## ----eval= FALSE--------------------------------------------------------------
# data.fst <- readr::read_tsv(file = "https://raw.githubusercontent.com/thierrygosselin/package_data/master/assigner_data_fst.tsv.gz")

## ----eval= FALSE--------------------------------------------------------------
# fst <- assigner::fst_WC84(
#     data = data.fst,
#     pop.levels = c("pop1", "pop2", "pop3", "pop4", "pop5", "pop6", "pop7", "pop8", "pop9", "pop10", "pop11"),
#     pairwise = TRUE,
#     filename = "testing_fst", #will trigger the function to write the results in directory as well
#     verbose = TRUE
# )

## ----eval= FALSE--------------------------------------------------------------
# names(fst)

## ----eval= FALSE--------------------------------------------------------------
# df <- fst$pairwise.fst

## ----eval= FALSE--------------------------------------------------------------
# # to see as a tibble:
# fst.matrix <- tibble::as_tibble(fst$pairwise.fst.full.matrix, rownames = "POP")
# # to keep matrix:
# fst.matrix <- fst$pairwise.fst.full.matrix

## ----eval= FALSE--------------------------------------------------------------
# fst$pairwise.fst.ci.matrix # you will get this:
# [1] "confidence intervals not selected"

## ----eval= FALSE--------------------------------------------------------------
# fst.ci <- fst_WC84(data = data.fst,
#     pop.levels = c("pop1", "pop2", "pop3", "pop4", "pop5", "pop6", "pop7", "pop8", "pop9", "pop10", "pop11"),
#     pairwise = TRUE,
#     ci = TRUE,
#     iteration.ci = 100,
#     quantiles.ci = c(0.025, 0.975),
#     parallel.core = 12,
#     heatmap.fst = TRUE,
#     filename = "testing_fst",
#     verbose = TRUE
# )

## ----eval= FALSE--------------------------------------------------------------
# fst.ci$pairwise.fst.ci.matrix

## ----eval= FALSE--------------------------------------------------------------
# pairwise.fst.ci.df <- tibble::as_tibble(fst.ci$pairwise.fst.ci.matrix, rownames = "POP")
# # if you have numeric pop_id identifier you might have to do this to get proper column names:
# colnames(pairwise.fst.ci.df) <- colnames(fst.ci$pairwise.fst.ci.matrix)
# # to save:
# readr::write_tsv(x = pairwise.fst.ci.df, path = "pairwise.fst.ci.df.tsv")

## ----eval= FALSE--------------------------------------------------------------
# # build the tree:
# require(ape)
# tree <- ape::nj(X = fst.matrix) # fst.matrix as a matrix
# # for annotating bootstraps values on the tree:
# bootstrap.value <- ape::boot.phylo(
#     phy = tree,
#     x = fst.matrix,
#     FUN = function(x) ape::nj(x),
#     block = 1,
#     B = 10000,
#     trees = FALSE,
#     rooted = FALSE
#     )
#  # to get percentage values
# bootstrap.value <- round((bootstrap.value/10000)*100, 0)
# bootstrap.value
# # to include in the tree
# tree$node.label <- bootstrap.value

## ----eval= FALSE--------------------------------------------------------------
# require(stats)
# tree <- ape::as.phylo(stats::hclust(stats::dist(fst.matrix), method = "average"))
# bootstrap.value <- ape::boot.phylo(phy = tree, x = fst.matrix, FUN = function(xx) ape::as.phylo(stats::hclust(stats::dist(xx), method = "average")) , block = 1, B = 10000, trees = FALSE, rooted = TRUE)
# bootstrap.value <- round((bootstrap.value/10000)*100, 0)
# bootstrap.value
# tree$node.label <- bootstrap.value

## ----eval= FALSE--------------------------------------------------------------
# # get the latest development version of ggtree:
# if (!require("ggtree")) install_github("GuangchuangYu/ggtree")
# # If it's not working, use the Bioconductor version:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("ggtree")

## ----eval= FALSE--------------------------------------------------------------
# require(ggtree)
# require(ggplot2)
# tree.figure <- ggplot2::ggplot(tree, ggplot2::aes(x, y), ladderize = TRUE) +
#     ggtree::geom_tree() +
#     # geom_tiplab(size = 3, hjust = -0.05, vjust = 0.5)+ # for just the tip label
#     ggplot2::geom_text(ggplot2::aes(label = label), size = 3, hjust = -0.05, vjust = 0.5) + # for both tips and nodes
#     ggtree::theme_tree() +
#     ggplot2::xlim(0, 0.16) # to allocate more space for tip labels (trial/error)
# tree.figure
# ggplot2::ggsave(filename = "assigner.fst.tree.example.pdf", width = 15, height = 15, dpi = 600, units = "cm", useDingbats = FALSE)

