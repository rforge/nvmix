\name{pmultinorm}
\alias{pmultinorm}
\title{Distribution Function of the Multivariate Normal Distribution}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t} for non-integer degrees of freedom).
}
\usage{
pmultinorm(upper, lower = rep(-Inf, length(upper)), mean = rep(0, length(upper)), 
       scale, standardized = FALSE, gam = 3.3, abserr = 0.001, Nmax = 1e8, 
       B = 12, n_init = 2^6, precond = TRUE, method = "sobol")
}
\arguments{
  \item{upper}{vector of length \eqn{d}.}
  \item{lower}{vector of length \eqn{d}.}
  \item{mean}{mean vector of length \eqn{d}}
  \item{scale}{positive definite \eqn{(d,d)}-covariance matrix.}
  \item{standardized}{\code{logical}. If \code{TRUE}, \code{scale} is assumed to be a correlation matrix; if \code{FALSE} (default), lower, upper and scale will be normalized.}
  \item{abserr}{numeric and non-negative. Absolute precision required. If \code{abserr = 0}, algorithm will run 
        until total number of function evaluations exceeds \code{Nmax}.}
  \item{gam}{Monte Carlo confidence multiplier. Algorithm runs until  \eqn{estimated standard error < gam * abserr}.
      \code{gam = 3.3} (the default) means that one can expect that in 99.9 percent of the cases the actual absolute error is less than \code{abserr}.}
  \item{Nmax}{maximum number of function evaluations, can be used to
    control run time.}
  \item{B}{number of repetitions to get an error estimate in the
    randomized quasi-Monte Carlo approach.}
  \item{n_init}{size of the first point set being used to estimate
    the probability. Any positive integer allowed, powers or at least multiples of 2 are recommended for \code{method = sobol}.}
  \item{precond}{\code{\link{logical}} indicating if preconditioning
    is applied, that is, reordering the variables.}
  \item{method}{Character string indicating method to be used. Allowed are \code{"sobol"}, \code{"ghalton"} and \code{"prng"}.}
}
\value{
  \code{pmultinorm()} returns a list of length five, containing the
  the estimated probabilities, the number of iterations, the total
  number of function evaluations, an error estimate and the estimated variance of the randomized Quasi Monte Carlo estimator. 
}
\details{
  \code{pmultinorm()} is a user-friendly wrapper and calls \code{pnvmix(..., mix = "constant")}, see \code{\code{\link{pnvmix}()}. 
  In the univariate case, this function calls \code{\link{pnorm}()}.
  
  Note that this procedure calls underlying C code. Currently, the
  dimensions \eqn{d\ge 16510}{d >= 16510} are not supported for the default method sobol.
  
  Care should be taken when changing the algorithm-specific parameters, notably \code{B}, \code{Nmax}, \code{method} and \code{precond}. Error estimates will not be reliable for too small \code{B} and the performance of the algorithm depends heavily on the (Quasi-) Monte Carlo point-set used. 
  
  If the absolute error tolerance \code{abserr} cannot be achieved with \code{Nmax} function evaluations, an additional warning will be returned. 
  
}
\author{Marius Hofert, Erik Hintz and Christiane Lemieux}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\examples{
## Generate a random correlation matrix in three dimensions
d <- 5
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
## Evaluate N(0,P) distribution function
a <- runif(d) * sqrt(d) * (-3) # random lower limit
b <- runif(d) * sqrt(d) * 3 # random upper limit
pn <- pmultinorm(upper = b, lower = a, scale = P)
stopifnot(all.equal(pn$Prob, 0.46884, tol = 5e-4))
}
\keyword{distribution}