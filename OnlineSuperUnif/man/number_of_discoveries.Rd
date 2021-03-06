% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_functions.R
\name{number_of_discoveries}
\alias{number_of_discoveries}
\title{number_of_discoveries}
\usage{
number_of_discoveries(rej_index, alternative_index, error_metric)
}
\arguments{
\item{rej_index}{A vector containing the indices of the rejected hypothesis.}

\item{alternative_index}{A vector containing the indices (in the stream of hypothesis) of the signal.}

\item{error_metric}{A string, either "FWER" or "mFDR" to indicate the error metric the user is studying.}
}
\value{
A list containing
ratio_true_discoveries : Ratio between the nb of true discoveries
and the number of non-nulls (= signals),
Nb of true discoveries,
error_quantity : depending on the error metric;
either a boolean stating the presence of a false discovery (FWER),
or the number of false discoveries (mFDR).
}
\description{
Function that allows to get the necessary quantities to estimate the
error (power, FWER or mFDR).
}
