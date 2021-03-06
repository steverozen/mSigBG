% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SeparateBackground.R
\name{SeparateSignatureFromBackground}
\alias{SeparateSignatureFromBackground}
\title{Estimate a signature from experimentally exposed spectra minus a background signature.}
\usage{
SeparateSignatureFromBackground(
  spectra,
  bg.sig.info,
  m.opts = NULL,
  start.b.fraction = 0.1
)
}
\arguments{
\item{spectra}{The spectra from which to subtract the background,
as a matrix or \code{\link[ICAMS]{ICAMS}} catalog.}

\item{bg.sig.info}{Information about the background signature. See
\code{\link{background.info}}.}

\item{m.opts}{Options to pass to \code{\link[nloptr]{nloptr}}.}

\item{start.b.fraction}{The estimated fraction of the mutations in
\code{spectra} due to the background signature.}
}
\value{
A list with the elements \describe{

\item{\code{inferred.target.sig}}{The estimated target signature as a numerical
vector.}

\item{exposures.to.target.sig}{The estimated total number of mutations due
to the target signature in each input spectrum.}

\item{\code{exposures.to.bg.sig}}{The estimated total number of mutations due
to the background in each input spectrum.}

\item{\code{message}}{The \code{message} element of \code{all.opt.ret}.}

\item{\code{all.opt.ret}}{The entire return value from the optimization.
See \code{\link[nloptr]{nloptr}}}
}
}
\description{
We index mutation channels (e.g. \code{ACA > AAA, ACC > AAC, ...}) by \eqn{j}, \eqn{j \in 1...96}.

We index input mutational spectra by \eqn{i}.

Let

\eqn{g = g_1, g_2,\ldots , g_{96}}, with \eqn{\Sigma g_j = 1},
be the previously determined, input background signature profile,

\eqn{s^i, i \in 1, 2,\ldots} be the input spectra,
from exposed samples, usually only 2 or 3,

\eqn{b^i, i \in 1, 2,\ldots} be the (to-be-estimated)
numbers of mutations due to the background signature in each
\eqn{s^i}, and

\eqn{t = t_1, t_2,\ldots , t_{96}}, with \eqn{\Sigma t_j = 1},
be the (to-be-estimated) target signature due to an exposure.

We want to maximize \eqn{\Pi^iP(s^i|b^i,t)P(b^i)} over
\eqn{b^1, b^2,\dots} and \eqn{t}.  (Note that the code
actually minimizes the additive inverse of this.)

\eqn{P(b^i)} is estimated from the distribution of previously observed
numbers of mutations in untreated samples, with the additional
constraint that \eqn{b^i \le |s^i|}), where \eqn{|s^i|} is defined as
the total number of mutations in spectrum \eqn{s^i}, i.e.
\eqn{|s^i| = \Sigma_j s^i_j}, \eqn{j \in 1...96}.

\eqn{P(s^i|b^i,t)} is estimated as follows:

The expected number of mutations in each
mutation category, \eqn{j}, is estimated as

\eqn{e^i_j = g_jb^i + t_j(|s^i| - b^i)}.

Then \eqn{P(s^i|e^i)} is estimated as \eqn{\Pi_jP(s^i_j|e^i_j)}.

\eqn{P(s^i_j|e^i_j)} is estimated from a negative binomial distribution
centered on each \eqn{e^i_j}; see the \code{sig.nbinom.size} elements of
the \code{\link{background.info}} package variables.
}
\details{
See \code{\link{ObjFn1}}.
}
