## ----Prepare your workspace, eval = FALSE-------------------------------------
#  library(assigner)

## ----load data, eval = FALSE--------------------------------------------------
#  data <- data_assigner_sim_01

## ----install gri_sim, eval = FALSE, include=TRUE------------------------------
#  assigner::install_gsi_sim(fromSource = TRUE)

## ----THL analysis with gsi_sim, eval = FALSE----------------------------------
#  test1 <- assigner::assignment_ngs(
#    data = data,
#    assignment.analysis = "gsi_sim",
#    markers.sampling = "ranked",
#    thl = 0.2,
#    iteration.method = 5
#  )
#  #>################################################################################
#  #>########################## assigner::assignment_ngs ############################
#  #>################################################################################
#  #>Execution date/time: 20190501@1104
#  #>Assignment analysis with gsi_sim
#  #>Folder created: assignment_analysis_method_ranked_20190501@1104
#  #>Calibrating REF/ALT alleles...
#  #>Subsampling: not selected
#  #>Conducting Assignment analysis using Training, Holdout, Leave-one-out
#  #>Using training samples to rank markers based on Fst
#  #>Holdout samples saved in your folder
#  #>Starting parallel computations, for progress monitor activity in folder...
#  #>
#  #>Computation time, overall: 7 sec
#  #>########################## assignment_ngs completed ############################

## ----returned in test1, eval = FALSE------------------------------------------
#  names(test1)
#  #>[1] "assignment"      "assignment.plot"

## ----test1 in directory, eval = FALSE-----------------------------------------
#  # 01_radiator_tidy_genomic: folder
#  # assigner_assignment_ngs_args_20190501@1102.tsv: tibble, file
#  # assignment_1: folder
#  # assignment_2: folder
#  # assignment_3: folder
#  # assignment_4: folder
#  # assignment_5: folder
#  # assignment.plot.pdf: figure
#  # assignment.ranked.results.iterations.raw.tsv: tibble, file
#  # assignment.ranked.results.iterations.summary.tsv: tibble, file
#  # assignment.results.summary.stats.tsv: tibble, file
#  # holdout.individuals.tsv: tibble, file

## ----thl figure, eval = FALSE-------------------------------------------------
#  test1$plot.assignment

## ----y axis full range, message=FALSE, warning=TRUE, eval = FALSE-------------
#  test1$plot.assignment + ggplot2::scale_y_continuous(limits = c(0,100))

## ----polymorphic markers in common, eval = FALSE------------------------------
#  data %<>%
#    radiator::filter_monomorphic(data = .) %>%
#    radiator::filter_common_markers(data = .)
#  #>Filter monomorphic markers
#  #>Number of individuals / strata / chrom / locus / SNP:
#  #>    Blacklisted: 0 / 0 / NA / NA / 3
#  #>
#  #>Filter common markers:
#  #>Number of individuals / strata / chrom / locus / SNP:
#  #>    Blacklisted: 0 / 0 / 0 / 0 / 0

## ----THL analysis with gsi_sim and subsample, eval = FALSE--------------------
#  test2 <- assigner::assignment_ngs(
#    data = data,
#    assignment.analysis = "gsi_sim",
#    markers.sampling = "ranked",
#    thl = 0.2,
#    iteration.method = 5,
#    marker.number = c(100, 200, 300, 400, "all"),
#    subsample = 30,
#    iteration.subsample = 3
#  )
#  #> ################################################################################
#  #> ########################## assigner::assignment_ngs ############################
#  #> ################################################################################
#  #> Execution date/time: 20190501@1158
#  #> Assignment analysis with gsi_sim
#  #> Folder created: assignment_analysis_method_ranked_20190501@1158
#  #> Calibrating REF/ALT alleles...
#  #> Subsampling: selected
#  #>     using subsample size of: 30
#  #>
#  #> Analyzing subsample: 1
#  #> Conducting Assignment analysis using Training, Holdout, Leave-one-out
#  #> Using training samples to rank markers based on Fst
#  #> Holdout samples saved in your folder
#  #> Starting parallel computations, for progress monitor activity in folder...
#  #>
#  #> Analyzing subsample: 2
#  #> Conducting Assignment analysis using Training, Holdout, Leave-one-out
#  #> Using training samples to rank markers based on Fst
#  #> Holdout samples saved in your folder
#  #> Starting parallel computations, for progress monitor activity in folder...
#  #>
#  #> Analyzing subsample: 3
#  #> Conducting Assignment analysis using Training, Holdout, Leave-one-out
#  #> Using training samples to rank markers based on Fst
#  #> Holdout samples saved in your folder
#  #> Starting parallel computations, for progress monitor activity in folder...
#  #>
#  #> Computation time, overall: 19 sec
#  #> ########################## assignment_ngs completed ############################

## ----output folder subsampling, eval = FALSE----------------------------------
#  # 01_radiator_tidy_genomic: folder
#  # assigner_assignment_ngs_args_20190501@1540.tsv: tibble, file
#  # assignment.plot.pdf: figure
#  # assignment.ranked.results.summary.stats.all.subsamples.tsv: tibble, file
#  # assignment.results.summary.stats.tsv: tibble, file
#  # subsample_1: folder
#  # subsample_2: folder
#  # subsample_3: folder
#  # subsampling_individuals.tsv: tibble, file

## ----y axis full range test2, message=FALSE, warning=TRUE, eval = FALSE-------
#  test2$plot.assignment + ggplot2::scale_y_continuous(limits = c(0,100))

## ----fst dataset1, eval = FALSE-----------------------------------------------
#  assigner::fst_WC84(data) %$% fst.overall$FST
#  #>[1] 0.39603

## ----load data low fst, eval = FALSE------------------------------------------
#  data <- data_assigner_sim_02

## ----THL analysis on low fst dataset using gsi_sim and subsample, eval = FALSE----
#  test3 <- assigner::assignment_ngs(
#    data = data,
#    assignment.analysis = "gsi_sim",
#    markers.sampling = "ranked",
#    thl = 0.2,
#    iteration.method = 5,
#    marker.number = c(100, 200, 300, 400, "all"),
#    subsample = 30,
#    iteration.subsample = 3
#  )

## ----y axis full range test3, message=FALSE, warning=TRUE, eval = FALSE-------
#  test3$plot.assignment + ggplot2::scale_y_continuous(limits = c(0,100))
#  # <img src="assignment_thl_test3.png">: works
#  #![](assignment_thl_test3.png): works
#  #knitr::include_graphics("assignment_thl_test3.png"):works

## ---- echo = FALSE------------------------------------------------------------
knitr::include_graphics("assignment_thl_test3.png")

## ----fst dataset2, eval = FALSE-----------------------------------------------
#  assigner::fst_WC84(data) %$% fst.overall$FST
#  #>[1] 0.001320833

## ----loo, eval = FALSE--------------------------------------------------------
#  test4 <- assigner::assignment_ngs(
#    data = data,
#    assignment.analysis = "gsi_sim",
#    markers.sampling = "random",
#    marker.number = "all"
#  )
#  #> ################################################################################
#  #> ########################## assigner::assignment_ngs ############################
#  #> ################################################################################
#  #> Execution date/time: 20190501@1317
#  #> Assignment analysis with gsi_sim
#  #> Folder created: assignment_analysis_method_random_20190501@1317
#  #> Calibrating REF/ALT alleles...
#  #> Subsampling: not selected
#  #> Conducting Assignment analysis with markers selected randomly
#  #> Making a list containing all the markers combinations
#  #> Starting parallel computations, for progress monitor activity in folder...
#  #> Summarizing the assignment analysis results by iterations and marker group
#  #> Compiling results
#  #> ########################## assignment_ngs completed ############################

## ----loo y axis full range test4, message=FALSE, warning=TRUE, eval = FALSE----
#  test4$plot.assignment + ggplot2::scale_y_continuous(limits = c(0,100))

## ---- echo = FALSE------------------------------------------------------------
knitr::include_graphics("assignment_loo_test4.png")

