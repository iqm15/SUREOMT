% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{saffron_OnlineSuperUnif}
\alias{saffron_OnlineSuperUnif}
\title{saffron_OnlineSuperUnif}
\usage{
saffron_OnlineSuperUnif(alpha, w0, raw_pvalues, gamma, lambda, capping = FALSE)
}
\arguments{
\item{alpha}{A numeric for the desired level of type I error control.}

\item{w0}{A numeric in \link{0, alpha} representing the amount of earning back one wants.}

\item{gamma}{A vector: the gamma spending sequence.}

\item{lambda}{A numeric in \link{0, 1} : the threshold for adaptivity.}

\item{raw.pvalues}{A vector containing raw p-values.}

\item{pCDFlist}{A list containing the support of the discrete CDF of the p-values.
Each support is represented by a vector in increasing order.}

\item{gamma_prime}{A vector: the gamma prime smoothing sequence (it can be the same as gamma). Default to NULL.}
}
\value{
A list containing a vector of the sequence of critical values and
a vector of the indices of rejected hypothesis.
}
\description{
Function that computes SAFFRON of Ramdas et al. (2018).
}
