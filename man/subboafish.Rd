% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{subboafish}
\alias{subboafish}
\title{Returns the Fisher Information matrix and its inverse for the AEP}
\usage{
subboafish(size = 1L, bl = 2, br = 2, m = 0, al = 1, ar = 1, O_munknown = 0L)
}
\arguments{
\item{size}{(numeric) - number of observations (Default: 01)}

\item{bl}{(numeric) - set the left exponent (Default: 2.0)}

\item{br}{(numeric) - set the right exponent (Default: 2.0)}

\item{m}{(numeric) - the location parameter (Default: 0.0)}

\item{al}{(numeric) - the left scale parameter (Default: 1.0)}

\item{ar}{(numeric) - the right scale parameter (Default: 1.0)}

\item{O_munknown}{(numeric) - if true assumes \eqn{m} is known}
}
\value{
a list containing three elements:
\itemize{
\item std_error - the standard error for the parameters
\item infmatrix - the Fisher Information Matrix
\item inv_infmatrix - the Inverse Fisher Information Matrix
}
}
\description{
Returns the Fisher Information matrix and its inverse for the Asymmetric
Power Exponential distribution for the given parameters.
}
