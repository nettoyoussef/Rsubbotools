% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{subboafit}
\alias{subboafit}
\title{Fit an Asymmetric Power Exponential density via maximum likelihood}
\usage{
subboafit(
  data,
  verb = 0L,
  method = 6L,
  interv_step = 10L,
  provided_m_ = NULL,
  par = as.numeric(c(2, 2, 1, 1, 0)),
  g_opt_par = as.numeric(c(0.1, 0.01, 100, 0.001, 1e-05, 2)),
  itv_opt_par = as.numeric(c(0.01, 0.001, 200, 0.001, 1e-05, 5))
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

\item{method}{int - the steps that should be used to estimate the
parameters.
\itemize{
\item 0 no optimization perform - just return the log-likelihood from initial guess.
\item 1 global optimization not considering lack of smoothness in m
\item 2 interval optimization taking non-smoothness in m into consideration
}}

\item{interv_step}{int - the number of intervals to be explored after
the last minimum was found in the interval optimization. Default is 10.}

\item{provided_m_}{NumericVector - if NULL, the m parameter is estimated
by the routine. If numeric, the estimation fixes m to the given value.}

\item{par}{NumericVector - vector containing the initial guess for
parameters bl, br, al, ar and m, respectively. Default values of are
c(2, 2, 1, 1, 0).}

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

\item{itv_opt_par}{NumericVector - interval optimization parameters. Fields
are the same as the ones for the global optimization. Default values
are c(.01, 1e-3, 200, 1e-3, 1e-5, 5).}
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
\code{subboafit} returns the parameters, standard errors. negative
log-likelihood and covariance matrix of the asymmetric power exponential for
a sample. The process can execute two steps, dependending on the level of
accuracy required. See details below.
}
\details{
The AEP is a exponential power distribution controlled
by five parameters, with formula:
\deqn{ f(x;a_l,a_r,b_l,b_r,m) =
\begin{cases}
\frac{1}{A} e^{- \frac{1}{b_l} |\frac{x-m}{a_l}|^{b_l} }, & x < m \\
\frac{1}{A} e^{- \frac{1}{b_r} |\frac{x-m}{a_r}|^{b_r} }, & x > m
\end{cases} }
with:
\deqn{A = a_lb_l^{1/b_l}\Gamma(1+1/b_l) + a_rb_r^{1/b_r}\Gamma(1+1/b_r)}
where \eqn{l} and \eqn{r} represent left and right tails, \eqn{a*} are
scale parameters, \eqn{b*} control the tails (lower values represent
fatter tails), and \eqn{m} is a location parameter. Due to its lack of
simmetry, and differently from the Subbotin, there is no simple equations
available to use the method of moments, so we start directly by minimizing
the negative log-likelihood. This global optimization is executed without
restricting any parameters. If required (default), after the global
optimization is finished, the method proceeds to iterate over the intervals
between several two observations, iterating the same algorithm of the
global optimization. The last method happens because of the lack of
smoothness on the \eqn{m} parameter, and intervals must be used since the
likelihood function doesn't have a derivative whenever \eqn{m} equals a
sample observation. Due to the cost, these iterations are capped at most
\emph{interv_step} (default 10) from the last minimum observed.
Details on this method are available on the package vignette.
}
\examples{
sample_subbo <- rpower(1000, 1, 2)
subboafit(sample_subbo)
}
