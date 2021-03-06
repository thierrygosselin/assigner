% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{import_subsamples_fst}
\alias{import_subsamples_fst}
\title{Import the fst ranking from all the subsample runs inside
an assignment folder.}
\usage{
import_subsamples_fst(dir.path)
}
\arguments{
\item{dir.path}{The path to the directory containing the subsample folders.}
}
\value{
A data frame of all the Fst and ranking.
}
\description{
This function will import all the fst ranking from all the
subsample runs inside an assignment folder.
}
\examples{
\dontrun{
subsamples.data <- import_subsamples_fst(
dir.path = "assignment_analysis_method_ranked_no_imputations_20160425@2321"
)
}
}
\author{
Thierry Gosselin \email{thierrygosselin@icloud.com}
}
