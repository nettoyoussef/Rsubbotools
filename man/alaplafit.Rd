% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{alaplafit}
\alias{alaplafit}
\title{Fit an Asymmetric Laplace Distribution via maximum likelihood}
\usage{
alaplafit(data, verb = 0L, interv_step = 10L, provided_m_ = NULL)
}
\arguments{
\item{data}{(NumericVector) - the sample used to fit the distribution.}

\item{verb}{(int) - the level of verbosity. Select one of:
\itemize{
\item 0  just the final result
\item 1  details of optim. routine
}}

\item{interv_step}{int - the number of intervals to be explored after
the last minimum was found in the interval optimization. Default is 10.}

\item{provided_m_}{NumericVector - if NULL, the m parameter is estimated
by the routine. If numeric, the estimation fixes m to the given value.}
}
\value{
a list containing the following items:
\itemize{
\item "dt" - dataset containing parameters estimations and standard deviations.
\item "log-likelihood" - negative log-likelihood value.
\item "matrix" - the covariance matrix for the parameters.
}
}
\description{
\code{alaplafit} returns the parameters, standard errors. negative
log-likelihood and covariance matrix of the Asymmetric Laplace Distribution
for a sample. See details below.
}
\details{
The Asymmetric Laplace distribution is a distribution controlled
by three parameters, with formula:
\deqn{f(x;a_l,a_r,m) =
\begin{cases}
\frac{1}{A} e^{-|\frac{x-m}{a_l}| }, & x < m \\
\frac{1}{A} e^{-|\frac{x-m}{a_r}| }, & x > m
\end{cases}}
with:
\deqn{A = a_l + a_r}
where \eqn{a*} are scale parameters, and \eqn{m} is a location parameter.
It is basically derived from the Asymmetric Exponential Power distribution
by setting \eqn{b_l = b_r = b}.
The estimations are produced by maximum likelihood, where
analytical formulas are available for the \eqn{a*} parameters.
The \eqn{m} parameter is found by an iterative method, using the median as
the initial guess. The method explore intervals around the last minimum
found, similar to the \code{subboafit} routine.
Details on this method can be found on the package vignette.
}
\examples{
sample_subbo <- rpower(1000, 1, 1)
alaplafit(sample_subbo)
}
