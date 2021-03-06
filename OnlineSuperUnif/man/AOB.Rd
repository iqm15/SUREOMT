% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/procedures.R
\name{AOB}
\alias{AOB}
\title{AOB}
\usage{
AOB(alpha, raw.pvalues, gamma, lambda)
}
\arguments{
\item{alpha}{A numeric for the desired level of type I error control.}

\item{raw.pvalues}{A vector containing the raw p-values.}

\item{gamma}{A vector: the gamma spending sequence.}

\item{lambda}{A numeric in \link{0, 1} : the adaptivity threshold.}

\item{pCDFlist}{A list containing the support of the discrete CDF of the p-values.
Each support is represented by a vector in increasing order.}
}
\value{
A list containing a vector of the sequence of critical values and
a vector of the indices of rejected hypothesis.
}
\description{
Function that can compute AOB (Adaptive-Spending of Tian & Ramdas (2021)) procedure.
}
