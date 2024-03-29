% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{sepfit}
\alias{sepfit}
\title{Fit a Skewed Exponential Power density via maximum likelihood}
\usage{
sepfit(
  data,
  verb = 0L,
  par = as.numeric(c(0, 1, 0, 2)),
  g_opt_par = as.numeric(c(0.1, 0.01, 100, 0.001, 1e-05, 2))
)
}
\arguments{
\item{data}{(NumericVector) - the sample used to fit the distribution.}

\item{verb}{(int) - the level of verbosity. Select one of:
\itemize{
\item 0  just the final result
\item 1  headings and summary table
\item 2  intermediate steps results
\item 3  intermediate steps internals
\item 4+  details of optim. routine
}}

\item{par}{NumericVector - vector containing the initial guess for
parameters mu, sigma, lambda and alpha, respectively. Default values of are
c(0, 1, 0, 2).}

\item{g_opt_par}{NumericVector - vector containing the global optimization
parameters.
The optimization parameters are:
\itemize{
\item step  - (num) initial step size of the searching algorithm.
\item tol   - (num) line search tolerance.
\item iter  - (int) maximum number of iterations.
\item eps   - (num) gradient tolerance. The stopping criteria is \eqn{||\text{gradient}||<\text{eps}}.
\item msize - (num) simplex max size. stopping criteria given by \eqn{||\text{max edge}||<\text{msize}}
\item algo  - (int) algorithm. the optimization method used:
\itemize{
\item 0 Fletcher-Reeves
\item 1 Polak-Ribiere
\item 2 Broyden-Fletcher-Goldfarb-Shanno
\item 3 Steepest descent
\item 4 Nelder-Mead simplex
\item 5 Broyden-Fletcher-Goldfarb-Shanno ver.2
}
}

Details for each algorithm are available on the \href{https://www.gnu.org/software/gsl/doc/html/multimin.html}{GSL Manual}.
Default values are c(.1, 1e-2, 100, 1e-3, 1e-5, 2).}
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
\code{sepfit} returns the parameters, standard errors. negative
log-likelihood and covariance matrix of the skewed power exponential
for a sample. The process performs a global minimization over the negative
log-likelihood function. See details below.
}
\details{
The  SEP is a exponential power distribution controlled
by four parameters, with formula:
\deqn{ f(x; \mu, \alpha, \lambda, \sigma) =
2 \Phi(w) e^{-|z|^\alpha/\alpha}/ ( \sigma C)}
where:
\deqn{z = (x-\mu)/\sigma}
\deqn{w = sign(z) |z|^{(\alpha/2)} \lambda \sqrt{2/\alpha}}
\deqn{C = 2 \alpha^{(1/\alpha-1)} \Gamma(1/\alpha)}
with \eqn{\Phi} the cumulative normal distribution with mean zero and variance
one.
Details on this method are available on the package vignette.
}
\examples{
sample_subbo <- rpower(1000, 1, 2)
sepfit(sample_subbo)
}
