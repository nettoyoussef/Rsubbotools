% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{psep}
\alias{psep}
\title{Returns CDF from the Skewed Exponential Power distribution}
\usage{
psep(x, m = 0, a = 2, b = 1, lambda = 0)
}
\arguments{
\item{x}{\itemize{
\item vector with values to evaluate CDF.
}}

\item{m}{\itemize{
\item the location parameter.
}}

\item{a}{\itemize{
\item the scale parameter.
}}

\item{b}{\itemize{
\item the shape parameter
}}

\item{lambda}{\itemize{
\item the skewness parameter.
}}
}
\description{
The \code{psep} returns the Cumulative Distribution Function at point x for
the Skewed Exponential Power distribution with parameters \eqn{a}, \eqn{b}.
}
\details{
The  SEP is a exponential power distribution controlled
by four parameters, with formula:
\deqn{ f(x; m, b, a, \lambda) = 2 \Phi(w) e^{-|z|^b/b}/(c)}
where:
\deqn{z = (x-m)/a}
\deqn{w = sign(z) |z|^{(b/2)} \lambda \sqrt{2/b}}
\deqn{c = 2 ab^{(1/b)-1} \Gamma(1/b)}
with \eqn{\Phi} the cumulative normal distribution with mean zero and variance
one. The CDF is calculated through numerical integration using the GSL suite.
}
