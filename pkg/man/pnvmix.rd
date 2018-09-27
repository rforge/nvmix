\name{pnvmix}
\alias{pnvmix}
\title{Distribution Function of Multivariate Normal Variance Mixtures}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t} for non-integer degrees of freedom).
}
\usage{
pnvmix(upper, lower = rep(-Inf, length(upper)), loc = rep(0, length(upper)),
      scale, mix, mean.sqrt.mix = NA, standardized = FALSE, gam = 3.3, abstol = 0.001,
      Nmax = 1e8, B = 12, n_init = 2^6, precond = TRUE, method = "sobol", ...)
}
\arguments{
  \item{upper}{vector of length \eqn{d}.}
  \item{lower}{vector of length \eqn{d}.}
  \item{loc}{location vector of length \eqn{d}. If \code{mix} has a mean, this is the mean of the normal variance mixture distribution.}
  \item{scale}{positive definite \eqn{(d,d)}-covariance matrix.}
    \item{mix}{specification of the mixing variable \eqn{W}; see McNeil et
    al. (2015). Supported are the following types of specification (see
    also the examples below):
    \describe{
      \item{\code{\link{character}}:}{a \code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1} and thus the cdf of
        the multivariate normal distribution with mean vector
	 \code{loc} and covariance matrix \code{scale} results) and
	 \code{"inverse.gamma"} (in which case \eqn{W} is an inverse gamma distribution with shape and rate parameters
	 \code{df}/2 resulting in the cdf the multivariate
	 Student \emph{t} distribution with degrees of freedom
	 \code{df}; note that \code{df} needs to be provided via
	 the ellipsis argument then; see the examples below).}
      \item{\code{\link{list}}:}{a \code{\link{list}} of length at least
	one, where the first component is a \code{\link{character}}
	string specifying the base name of a distribution which has a quantile function
	accessible via prefix \code{"q"}; an example is \code{"exp"}
        for which \code{"qexp"} exists. If the list is
        of length larger than one, the remaining elements contain
        additional parameters of the distribution; for \code{"exp"},
        this can be the parameter \code{rate}.}
      \item{\code{\link{function}}:}{a \code{\link{function}}
	interpreted as the quantile function of the mixing
	variable \eqn{W}; internally, sampling is then done with the
	inversion method by applying this function to U(0,1) random variates via \code{mix(u, ...)}. Additional arguments for \code{mix} can be passed via the ellipsis argument.}
    }
  }
  \item{mean.sqrt.mix}{Mean of \code{sqrt(W)}, where \eqn{W} is the mixing variable. If not provided, it will be estimated. This is only needed for reordering, hence a rather
  crude approximation is fine.}
  \item{standardized}{\code{logical}. If \code{TRUE}, \code{scale} is assumed to be a correlation matrix; if \code{FALSE} (default), lower, upper and scale will be normalized.}
  \item{abstol}{numeric and non-negative. Absolute precision required. If abstol = 0, algorithm will run
        until total number of function evaluations exceeds Nmax (see also Nmax).}
  \item{gam}{Monte Carlo confidence multiplier. Algorithm runs until  \eqn{estimated standard error < gam * abstol}.
      gam = 3.3 (the default) means that one can expect that in 99.9 percent of the cases the actual absolute error is less than \eqn{abstol}.}
  \item{Nmax}{maximum number of function evaluations, can be used to
    control run time.}
  \item{B}{number of repetitions to get an error estimate in the
    randomized quasi-Monte Carlo approach.}
  \item{n_init}{size of the first point set being used to estimate
    the probability. Any positive integer allowed, powers or at least multiples of 2 are recommended for method = sobol}
  \item{precond}{\code{\link{logical}} indicating if preconditioning
    is applied, that is, reordering the variables.}
  \item{method}{Character string indicating method to be used. Allowed are "sobol", "ghalton" and "prng".}
  \item{\dots}{additional arguments containing parameters of
    mixing distributions when \code{mix} is a \code{\link{character}}
    string or \code{\link{function}}.}
}
\value{
  \code{pnvmix()} returns a list of length five, containing the
  the estimated probabilities, the number of iterations, the total
  number of function evaluations, an error estimate and the estimated variance of the randomized Quasi Monte Carlo estimator.
}
\details{
  Note that this procedure calls underlying C code. Currently, the
  dimensions \eqn{d\ge 16510}{d >= 16510} are not supported for the default method sobol.

  Care should be taken when changing the algorithm-specific parameters, notably \code{B}, \code{Nmax}, \code{method} and \code{precond}. Error estimates will not be reliable for too small \code{B} and the performance of the algorithm depends heavily on the (Quasi-) Monte Carlo point-set used.

  If the absolute error tolerance \code{abstol} cannot be achieved with \code{Nmax} function evaluations, an additional warning will be returned.

  For the case of a multivariate normal or multivariate t distributions, user-friendly wrappers are provided in \code{\link{pNorm}} and \code{\link{pStudent}}.

}
\author{Marius Hofert, Erik Hintz and Christiane Lemieux}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{pNorm}}(), \code{\link{pStudent}()}
}
\examples{
## Example 1: Multivariate t distribution
## Generate a random correlation matrix in three dimensions
d <- 3
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
## Evaluate t_{0.5} distribution function
df <- 0.5
a <- runif(d) * sqrt(d) * (-3) # random lower limit
b <- runif(d) * sqrt(d) * 3 # random upper limit
## Using "inverse.gamma"
set.seed(1) # result depends slightly on .random.seed
pt1 <- pnvmix(upper = b, lower = a, scale = P, mix = "inverse.gamma", df = df)
## The same can be achieved by defining a function:
mix. <- function(u, df){
  df2 <- df/2
  return(1 / qgamma(u, shape = df2, rate = df2))
}
## In this case meansqrtmix is known:
c <- sqrt(df) * gamma(df/2) / ( sqrt(2) * gamma( (df+1) / 2 ) )
set.seed(1)
pt2 <- pnvmix(upper = b, lower = a, scale = P, mix = mix., mean.sqrt.mix = c, df  = df)
stopifnot(all.equal(pt1, pt2, tol = 5e-4, check.attributes = FALSE))
## mean.sqrt.mix will be approximated internally if not provided.
## This leads to slightly different results
set.seed(1)
pt3 <- pnvmix(upper = b, lower = a, scale = P, mix = mix., df  = df)
print(abs(pt3 - pt2))
}
\keyword{distribution}