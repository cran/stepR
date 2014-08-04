\name{stepbound}
\alias{stepbound}
\alias{stepbound.default}
\alias{stepbound.stepcand}

\title{Jump estimation under restrictions}
\description{Computes piecewise constant maximum likelihood estimators with minimal number of jumps under given restrictions on subintervals.}

\usage{
stepbound(y, bounds, \dots)
\method{stepbound}{default}(y, bounds, x = 1:length(y), x0 = 2 * x[1] - x[2],
  max.cand = NULL, family = c("gauss", "gaussvar", "poisson", "binomial", "gaussKern"),
  param = NULL, weights = rep(1, length(y)), refit = y,
  jumpint = confband, confband = FALSE, \dots)
\method{stepbound}{stepcand}(y, bounds, refit = TRUE, \dots)
}

\arguments{
\item{y}{a vector of numerical observations}
\item{bounds}{bounds on the value allowed on intervals; typically computed with \code{\link{bounds}}}
\item{x}{a numeric vector of the same length as \code{y} containing the corresponding sample points}
\item{x0}{a single numeric giving the last unobserved sample point directly before sampling started}
\item{max.cand, weights}{see \code{\link{stepcand}}}
\item{family, param}{specifies distribution of data, see \link{family}}
\item{refit}{\code{\link{logical}}, for \code{family = "gaussKern"}; determines whether a fit taken the filter kernel into account will be computed at the end}
\item{jumpint}{\code{\link{logical}} (\code{FALSE} by default), indicates if confidence sets for jumps should be computed}
\item{confband}{\code{\link{logical}}, indicates if a confidence band for the piecewise-continuous function should be computed}
\item{\dots}{arguments to be passed to generic methods}
}

\value{
An object of class \code{\link{stepfit}} that contains the fit; if \code{jumpint = TRUE} function \code{\link{jumpint}} allows to extract the confidence interval for the jumps, if \code{confband = TRUE} function \code{\link{confband}} allows to extract the confidence band.
}

\seealso{\code{\link{bounds}}, \code{\link{smuceR}}, \code{\link{jsmurf}}, \code{\link{stepsel}}, \code{\link{stepfit}}, \code{\link{jumpint}}, \code{\link{confband}}}

\encoding{utf8}
\references{
Frick, K., Munk, A., and Sieling, H. (2013). Multiscale Change-Point Inference. \emph{eprint arXiv}:1301.7212 [stat.ME], \url{http://arxiv.org/abs/1301.7212}.
Hotz, T., \enc{Schütte}{Schuette}, O., Sieling, H., Polupanow, T., Diederichsen, U., Steinem, C., and Munk, A. (2013). Idealizing Ion Channel Recordings by a Jump Segmentation Multiresolution Filter. (submitted)
}

\examples{
# simulate poisson data with two levels
y <- rpois(100, c(rep(1, 50), rep(4, 50)))
# compute bounds
b <- bounds(y, penalty="len", family="poisson", q=4)
# fit step function to bounds
sb <- stepbound(y, b, family="poisson", confband=TRUE)
plot(y)
lines(sb)
# plot confidence intervals for jumps on axis
points(jumpint(sb), col="blue")
# confidence band
lines(confband(sb), lty=2, col="blue")
}

\keyword{nonparametric}
