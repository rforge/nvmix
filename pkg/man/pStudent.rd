\name{pStudent}
\alias{pStudent}
\title{Distribution Function of the Multivariate t Distribution}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t} for non-integer degrees of freedom).
}
\usage{
pStudent(upper, lower = rep(-Inf, length(upper)), shift = rep(0, length(upper)), scale, df, standardized = FALSE,  
       gam = 3.3, abserr = 0.001, Nmax = 1e8, N = 12, n_init = 2^5, precond = TRUE, 
       method = "sobol")
}
\arguments{
  \item{upper}{vector of length \eqn{d}.}
  \item{lower}{vector of length \eqn{d}.}
  \item{shift}{shift vector of length \eqn{d}. If \code{df > 1}, this is the mean of the distribution.}
  \item{scale}{positive definite \eqn{(d,d)}-covariance matrix.}
  \item{df}{degress of freedom (positive real or Inf in which case the corresponding probability for a N(0, scale) is returned.}
  \item{standardized}{\code{logical}. If \code{TRUE}, \code{scale} is assumed to be a correlation matrix; if \code{FALSE} (default), lower, upper and scale will be normalized.}
  \item{abserr}{numeric and non-negative. Absolute precision required. If \code{abserr = 0}, algorithm will run 
        until total number of function evaluations exceeds \code{Nmax}.}
  \item{gam}{Monte Carlo confidence multiplier. Algorithm runs until  \eqn{estimated standard error < gam * abserr}.
      \code{gam = 3.3} (the default) means that one can expect that in 99.9 percent of the cases the actual absolute error is less than \code{abserr}.}
  \item{Nmax}{maximum number of function evaluations, can be used to
    control run time.}
  \item{N}{number of repetitions to get an error estimate in the
    randomized quasi-Monte Carlo approach.}
  \item{n_init}{size of the first point set being used to estimate
    the probability. Any positive integer allowed, powers or at least multiples of 2 are recommended for \code{method = sobol}.}
  \item{precond}{\code{\link{logical}} indicating if preconditioning
    is applied, that is, reordering the variables.}
  \item{method}{Character string indicating method to be used. Allowed are \code{"sobol"}, \code{"ghalton"} and \code{"prng"}.}
}
\value{
  \code{pStudent()} returns a list of length five, containing the
  the estimated probabilities, the number of iterations, the total
  number of function evaluations, an error estimate and the estimated variance of the randomized Quasi Monte Carlo estimator. 
}
\details{
  Note that this procedure calls underlying C code. Currently, the
  dimensions \eqn{d\ge 16510}{d >= 16510} are not supported for the default method sobol.
  
  Care should be taken when changing the algorithm-specific parameters, notably \code{N}, \code{Nmax}, \code{method} and \code{precond}. Error estimates will not be reliable for too small \code{N} and the performance of the algorithm depends heavily on the (Quasi-) Monte Carlo point-set used. 
  
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
d <- 3
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
## Evaluate t_{3.5} distribution function
a <- runif(d) * sqrt(d) * (-3) # random lower limit
b <- runif(d) * sqrt(d) * 3 # random upper limit
pt <- pStudent(upper = b, lower = a, scale = P, df = 3.5)
stopifnot(all.equal(pt$Prob, 0.8061, tol = 5e-4))
}
\keyword{distribution}