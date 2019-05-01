#' Simulated dataset to test assigner
#'
#' A dataset of genotypes containing 500 bi-allelic SNPs, simulated 
#' for 250 individuals, 5 populations.
#'
#' @format A tibble with 125000 rows (genotypes) and 4 variables:
#' \describe{
#'   \item{MARKERS}{SNPs markers}
#'   \item{POP_ID}{Populations/strata for the samples}
#'   \item{INDIVIDUALS}{Samples id)}
#'   \item{GT}{Genotypes coded \emph{a la genepop} format}
#' }
#' @details 
#' Dataset simulation caracteristics:
#' \itemize{
#' \item num.pops: 5
#' \item num.loci: 1000
#' \item div.time: 25e3
#' \item ne: 200
#' \item nm: 0.5
#' \item theta: 0.2
#' \item mig.type: island
#' \item mut.rate: 2.5e-4
#' \item mig.rate: 0.0025
#' }
#' From this simulated dataset, 500 SNPs and 250 individuals (50 ind/pop) were
#' randomly selected.
#' @source The data was simulated with grur
#' \url{https://thierrygosselin.github.io/grur/reference/simulate_rad.html}

"data_assigner_sim_01"

#' Simulated dataset to test assigner
#'
#' A dataset of genotypes containing 500 bi-allelic SNPs, simulated 
#' for 250 individuals, 5 populations.
#'
#' @format A tibble with 125000 rows (genotypes) and 4 variables:
#' \describe{
#'   \item{MARKERS}{SNPs markers}
#'   \item{POP_ID}{Populations/strata for the samples}
#'   \item{INDIVIDUALS}{Samples id}
#'   \item{GT}{Genotypes coded \emph{a la genepop} format}
#' }
#' @details 
#' Dataset simulation caracteristics:
#' \itemize{
#' \item num.pops: 5
#' \item num.loci: 1000
#' \item div.time: 25e3
#' \item ne: 200
#' \item nm: 0.5
#' \item theta: 0.2
#' \item mig.type: island
#' \item mut.rate: 2.5e-4
#' \item mig.rate: 0.0025
#' }
#' From this simulated dataset, 500 SNPs and 250 individuals (50 ind/pop) were
#' randomly selected.
#' @source The data was simulated with grur
#' \url{https://thierrygosselin.github.io/grur/reference/simulate_rad.html}
"data_assigner_sim_02"
