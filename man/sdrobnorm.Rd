\name{sdrobnorm}
\alias{sdrobnorm}

\title{Robust standard deviation estimate}
\description{Robust estimation of the standard deviation of Gaussian data.}

\usage{
sdrobnorm(x, p = c(0.25, 0.75), lag = 1)
}

\arguments{
\item{x}{a vector of numerical observations}
\item{p}{vector of two distinct probabilities}
\item{lag}{the lag of the difference used, see \code{\link{diff}}}
}

\details{
Compares the difference between the estimated sample quantile corresponding to \code{p} after taking (\code{lag}ged) differences) with the corresponding theoretical quantiles of Gaussian white noise to determine the standard deviation under a Gaussian assumption; if the data contain (few) jumps, this will (on average) be a slight overestimate of the true standard deviation.
}

\value{
Returns the estimate of the sample's standard deviation.
}

\seealso{\code{\link{sd}}, \code{\link{diff}}, \link{family}}

\references{
Davies, P.L., Kovac, A. (2001). Local extremes, runs, strings and multiresolution. \emph{The Annals of Statistics} 29, 1-65.
}

\examples{
# simulate data sample
y <- rnorm(100, c(rep(1, 50), rep(10, 50)), 2)
# estimate standard deviation
sdrobnorm(y)
}

\keyword{nonparametric}
