\name{bounds}
\alias{[.bounds}
\alias{bounds}
\alias{bounds.MRC}

\title{Bounds based on MRC}
\description{Computes two-sided bounds for a collection of intervals based on a multiresolution criterion (MRC).}

\usage{

bounds(y, type = "MRC", \dots)
bounds.MRC(y, q, alpha = 0.05, r = ceiling(50 / min(alpha, 1 - alpha)),
  lengths = if(family == "gaussKern")
    2^(floor(log2(length(y))):ceiling(log2(length(param$kern)))) else
    2^(floor(log2(length(y))):0), penalty = c("none", "len", "var", "sqrt"),
  name = if(family == "gaussKern") ".MRC.ktable" else ".MRC.table", pos = .GlobalEnv,
  family = c("gauss", "gaussvar", "poisson", "binomial","gaussKern"), param = NULL,
  subset, max.iter = 1e2, eps = 1e-3)
\method{[}{bounds}(x, subset)
}

\arguments{
\item{y}{a numeric vector containing the serial data}
\item{type}{so far only bounds of type \code{"MRC"} are implemented}
\item{\dots}{further arguments to be passed on to \code{bounds.MRC}}
\item{q}{quantile of the MRC; if specified, \code{alpha} and \code{r} will be ignored}
\item{alpha}{level of significance}
\item{r}{number of simulations to use to obtain quantile of MRC for specified \code{alpha}}
\item{lengths}{vector of interval lengths to use, dyadic intervals by default}
\item{penalty}{penalty term in the multiresolution statistic: \code{"none"} for no penalty, \code{"len"} for penalizing the length of an interval, \code{"var"} for penalizing the variance over an interval, and \code{"sqrt"} for penalizing the square root of the MRC}
\item{family, param}{specifies distribution of data, see \link{family}}
\item{subset}{a subset of indices of \code{y} for which bounds should be aggregated}
\item{name,pos}{under which name and where precomputed results are stored, or retrieved, see \code{\link{assign}}}
\item{max.iter}{maximal iterations in Newton's method to compute non-Gaussian MRC bounds}
\item{eps}{tolerance in Newton's method}
\item{x}{an object of class \code{bounds}}
}

\value{
Returns an object of class \code{bounds}, i.e. a list whose entry \code{bounds} contains two-sided bounds (\code{lower} and \code{upper}) of the considered intervals (with left index \code{li} and right index \code{ri}) in a \code{\link{data.frame}}, along with a vector \code{start} specifying in which row of entry \code{bounds} intervals with corresponding \code{li} start (if any; specified as a C-style index), and a \code{\link{logical}} \code{feasible} telling whether a feasible solution exists for these bounds (always \code{TRUE} for MRC bounds which are not restricted to a \code{subset}).
}

\seealso{\code{\link{stepbound}}, \link{family}}

\examples{
# simulate signal of 100 data points
Y <- rpois(100, 1:100 / 10)
# compute bounds for intervals of dyadic lengths
b <- bounds(Y, penalty="len", family="poisson", q=4)
# compute bounds for all intervals
b <- bounds(Y, penalty="len", family="poisson", q=4, lengths=1:100)
}

\keyword{nonparametric}
