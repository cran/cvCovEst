% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_functions.R
\name{matrixMetrics}
\alias{matrixMetrics}
\title{General Matrix Metrics}
\usage{
matrixMetrics(estimate)
}
\arguments{
\item{estimate}{A \code{matrix} corresponding to a covariance matrix
estimate.}
}
\value{
A named \code{list} containing the three values.
}
\description{
\code{matrixMetrics} computes the condition number, sparsity,
and sign of a covariance matrix estimate.
}
\keyword{internal}
