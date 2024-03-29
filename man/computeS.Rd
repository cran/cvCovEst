% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spikedCovarianceHelpers.R
\name{computeS}
\alias{computeS}
\title{Estimate S of Spiked Covariance Matrix Estimator}
\usage{
computeS(c_donoho)
}
\arguments{
\item{c_donoho}{A \code{numeric} vector output by \code{\link{computeC}()}.}
}
\value{
A \code{numeric} vector.
}
\description{
\code{computeS()} computes the s(ell) value described in
\insertCite{donoho2018;textual}{cvCovEst}.
}
\references{
\insertAllCited{}
}
\keyword{internal}
