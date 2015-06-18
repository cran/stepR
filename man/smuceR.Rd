\name{smuceR}
\alias{smuceR}
\alias{thresh.smuceR}

\title{Piecewise constant regression with SMUCE}
\description{Computes the SMUCE estimator for one-dimensional data.}

\usage{
smuceR(y, x = 1:length(y), x0 = 2 * x[1] - x[2], q = thresh.smuceR(length(y)), alpha, r,
  lengths, family = c("gauss", "gaussvar", "poisson", "binomial"), param,
  jumpint = confband, confband = FALSE)
thresh.smuceR(v)
}

\arguments{
\item{y}{a numeric vector containing the serial data}
\item{x}{a numeric vector of the same length as \code{y} containing the corresponding sample points}
\item{x0}{a single numeric giving the last unobserved sample point directly before sampling started}
\item{q}{threshold value, by default chosen automatically according to Frick et al.~(2013)}
\item{alpha}{significance level; if set to a value in (0,1), \code{q} is chosen as the corresponding quantile of the asymptotic (if \code{r} is not given) null distribution (and any value specified for \code{q} is silently ignored)}
\item{r}{numer of simulations; if specified along \code{alpha}, \code{q} is chosen as the corresponding quantile of the simulated null distribution}
\item{lengths}{length of intervals considered; by default up to a sample size of 1000 all lengths, otherwise only dyadic lengths}
\item{family, param}{specifies distribution of data, see \link{family}}
\item{jumpint}{\code{\link{logical}} (\code{FALSE} by default), indicates if confidence sets for change-points should be computed}
\item{confband}{\code{\link{logical}}, indicates if a confidence band for the piecewise-continuous function should be computed}
\item{v}{number of data points}
}

\value{
For \code{smuceR}, an object of class \code{\link{stepfit}} that contains the fit; if \code{jumpint = TRUE} function \code{\link{jumpint}} allows to extract the \code{1- alpha} confidence interval for the jumps, if \code{confband = TRUE} function \code{\link{confband}} allows to extract the \code{1- alpha} confidence band.

For \code{thresh.smuceR}, a precomputed threshhold value, see reference.
}

\seealso{\code{\link{stepbound}}, \code{\link{bounds}}, \link{family}, \code{\link{MRC.asymptotic}}, \code{\link{sdrobnorm}}, \code{\link{stepfit}}}

\references{
Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point Inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} 76(3), 495-580.

Futschik, A., Hotz, T., Munk, A. Sieling, H. (2014). Multiresolution DNA partitioning: statistical evidence for segments. \emph{Bioinformatics},  30(16), 2255-2262.
} 

\examples{
# simulate poisson data with two levels
y <- rpois(100, c(rep(1, 50), rep(4, 50)))
# compute fit, q is chosen automatically
fit <- smuceR(y, family="poisson", confband = TRUE)
# plot result
plot(y)
lines(fit)
# plot confidence intervals for jumps on axis
points(jumpint(fit), col="blue")
# confidence band
lines(confband(fit), lty=2, col="blue")

# simulate binomial data with two levels
y <- rbinom(200,3,rep(c(0.1,0.7),c(110,90)))
# compute fit, q is the 0.9-quantile of the (asymptotic) null distribution
fit <- smuceR(y, alpha=0.1, family="binomial", param=3, confband = TRUE)
# plot result
plot(y)
lines(fit)
# plot confidence intervals for jumps on axis
points(jumpint(fit), col="blue")
# confidence band
lines(confband(fit), lty=2, col="blue")
}

\keyword{nonparametric}
