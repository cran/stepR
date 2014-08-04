\name{stepR-package}
\alias{stepR-package}
\alias{stepR}
\docType{package}

\title{Fit step-functions}

\description{Allows to fit step-functions to univariate serial data where the number of jumps is used as a penalty.}

\details{
\tabular{ll}{
Package: \tab stepR\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-17\cr
Authors@R: \tab c(person("Thomas", "Hotz", role = c("aut", "cre"), email = "thomas.hotz@tu-ilmenau.de"),
person("Hannes", "Sieling", role = "aut", email = "hsielin@uni-goettingen.de"))\cr
Depends: \tab R (>= 2.15.0), stats\cr
License: \tab GPL-3\cr
LazyData: \tab yes\cr
Built: \tab R 2.15.3; i686-pc-linux-gnu; 2013-09-17 07:05:08 UTC; unix\cr
}

Index:
\preformatted{
BesselPolynomial        Bessel Polynomials
MRC                     Compute Multiresolution Criterion
MRC.1000                Values of the MRC statistic for 1,000
                        observations (all intervals)
MRC.asymptotic          "Asymptotic" values of the MRC statistic (all
                        intervals)
MRC.asymptotic.dyadic   "Asymptotic" values of the MRC statistic
                        (dyadic intervals)
bounds                  Bounds based on MRC
compareBlocks           Compare fit blockwise with ground truth
contMC                  Continuous time Markov chain
dfilter                 Digital filters
family                  Family of distributions
jsmurf                  Reconstruct filtered piecewise constant
                        functions with noise
neighbours              Neighbouring integers
sdrobnorm               Robust standard deviation estimate
smuceR                  Piecewise constant regression with SMUCE
stepR-package           Fit step-functions
stepblock               Step function
stepbound               Jump estimation under restrictions
stepcand                Forward selection of candidate jumps
stepfit                 Fitted step function
steppath                Solution path of step-functions
stepsel                 Automatic selection of number of jumps
transit                 TRANSIT algorithm for detecting jumps
}

%Further information is available in the following vignettes:
%\tabular{ll}{
%\code{stepRusage} \tab Using stepR
%}
}

\author{
Thomas Hotz and Hannes Sieling

Maintainer: Thomas Hotz <thomas.hotz@tu-ilmenau.de>
}

\encoding{utf8}
\references{
Davies, P.L., Kovac, A. (2001). Local extremes, runs, strings and multiresolution. \emph{The Annals of Statistics} 29, 1-65.
Frick, K., Munk, A., and Sieling, H. (2014). Multiscale Change-Point Inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} 76(3), 495-580.
Friedrich, F., Kempe, A., Liebscher, V., Winkler, G. (2008). Complexity penalized M-estimation: fast computation. \emph{Journal of Computational and Graphical Statistics} 17(1), 201-224.
Futschik, A., Hotz, T., Munk, A. Sieling, H. (2013). Multiresolution DNA partitioning: statistical evidence for segments. \emph{Bioinformatics}, online first: \url{http://bioinformatics.oxfordjournals.org/content/early/2014/05/09/bioinformatics.btu180.full}.
Hotz, T., \enc{Schütte}{Schuette}, O., Sieling, H., Polupanow, T., Diederichsen, U., Steinem, C., and Munk, A. (2013). Idealizing Ion Channel Recordings by a Jump Segmentation Multiresolution Filter. \emph{IEEE Transactions on NanoBioscience} 12(4), 376-386.
VanDongen, A.M.J. (1996). A New Algorithm for Idealizing Single Ion Channel Data Containing Multiple Unknown Conductance Levels. \emph{Biophysical Journal} 70(3), 1303-1315.
}

\keyword{package,nonparametric}

\seealso{
\code{\link{smuceR}}, \code{\link{jsmurf}}, \code{\link{stepbound}}, \code{\link{steppath}}, \code{\link{stepsel}}, \code{\link{transit}}, \code{\link{compareBlocks}}, \link{family}
}

\examples{
# estimating step-functions with Gaussian white noise added
# simulate a Gaussian hidden Markov model of length 1000 with 2 states
# with identical transition rates 0.01, and signal-to-noise ratio 2
sim <- contMC(1e3, 0:1, matrix(c(0, 0.01, 0.01, 0), 2), param=1/2)
plot(sim$data, cex = 0.1)
lines(sim$cont, col="red")
# maximum-likelihood estimation under multiresolution constraints
fit.MRC <- smuceR(sim$data$y, sim$data$x)
lines(fit.MRC, col="blue")
# choose number of jumps using BIC
path <- steppath(sim$data$y, sim$data$x, max.blocks=1e2)
fit.BIC <- path[[stepsel.BIC(path)]]
lines(fit.BIC, col="green3", lty = 2)

# estimate after filtering
# simulate filtered ion channel recording with two states
set.seed(9)
# sampling rate 10 kHz
sampling <- 1e4
# tenfold oversampling
over <- 10
# 1 kHz 4-pole Bessel-filter, adjusted for oversampling
cutoff <- 1e3
df.over <- dfilter("bessel", list(pole=4, cutoff=cutoff / sampling / over))
# two states, leaving state 1 at 10 Hz, state 2 at 20 Hz
rates <- rbind(c(0, 10), c(20, 0))
# simulate 0.5 s, level 0 corresponds to state 1, level 1 to state 2
# noise level is 0.3 after filtering
Sim <- contMC(0.5 * sampling, 0:1, rates, sampling=sampling, family="gaussKern",
  param = list(df=df.over, over=over, sd=0.3))
plot(Sim$data, pch = ".")
lines(Sim$discr, col = "red")
# fit under multiresolution constraints using filter corresponding to sample rate
df <- dfilter("bessel", list(pole=4, cutoff=cutoff / sampling))
Fit.MRC <- jsmurf(Sim$data$y, Sim$data$x, param=df, r=1e2)
lines(Fit.MRC, col = "blue")
# fit using TRANSIT
Fit.trans <- transit(Sim$data$y, Sim$data$x)
lines(Fit.trans, col = "green3", lty=2)
}
% created with cat(unlist(promptPackage("stepR", filename=NA)), sep = "\n")
