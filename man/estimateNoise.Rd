% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spikedCovarianceHelpers.R
\name{estimateNoise}
\alias{estimateNoise}
\title{Estimate Noise in Spiked Covariance Matrix Model}
\usage{
estimateNoise(eig_vals, p_n_ratio)
}
\arguments{
\item{eig_vals}{A \code{numeric} vector of estimated eigenvalues.}

\item{p_n_ratio}{A \code{numeric} indicating the asymptotic ratio of the
number of features, p, and the number of observations, n. This ratio is
assumed to be between 0 and 1.}
}
\value{
A \code{numeric} estimate of the noise term in a spiked covariance
matrix model.
}
\description{
\code{estimateNoise()} estimates the unknown noise term in a
Gaussian spiked covariance matrix model, where the covariance matrix is
assumed to be the identity matrix multiplied by the unknown noise, save for
a few "spiked" entries. This procedures is described in
\insertCite{donoho2018;textual}{cvCovEst}.
}
\references{
\insertAllCited{}
}
\keyword{internal}
