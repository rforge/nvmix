\name{dnvmix}
\alias{dnvmix}
\title{Density of Multivariate Normal Variance Mixtures}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t}).
}
\usage{
dnvmix(x, loc = rep(0, d), scale, mix,
       factor = tryCatch(factorize(scale), error = function(e) e),
       abstol = 0.001, gam = 3.3, Nmax = 1e8, B = 12, n_init = 2^6, method = "sobol",
       log = FALSE, ...)
}
\arguments{
  \item{x}{\eqn{(n, d)}-\code{\link{matrix}} of evaluation points.}
  \item{loc}{location vector of dimension \eqn{d} (the number of columns
    of \code{x}).}
  \item{scale}{covariance matrix of dimension \eqn{(d, d)}.}
  \item{mix}{specification of the mixing variable \eqn{W}; see McNeil et
    al. (2015). Supported are the following types of specification (see
    also the examples below):
    \describe{
      \item{\code{\link{character}}:}{a \code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1} and thus the pdf of
        the multivariate normal distribution with mean vector
	 \code{loc} and covariance matrix \code{scale} results) and
	 \code{"inverse.gamma"} (in which case \eqn{W} is an inverse gamma distribution with shape and rate parameters
	 \code{df}/2 resulting in the pdf the multivariate
	 Student \emph{t} distribution with degrees of freedom
	 \code{df}; note that \code{df} needs to be provided via
	 the ellipsis argument then; see the examples below). In both cases, the closed formula for the density of a multivariate
	 t / normal distribution is used.}
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
  \item{factor}{factorization matrix of the covariance matrix
    \code{scale}. A \code{\link{matrix}} \eqn{R} with \eqn{d} rows
    such that \eqn{R^T R} equals \code{scale}. Although not tested
    (for efficiency reasons), \eqn{R} (thus \code{factor}) has
    to be \strong{upper triangular} (otherwise the determinant
        of \code{scale} involved in the density is not computed correctly).
    Also note that \code{factor} can be \code{NULL} or inherit from
    \code{"error"} in which case the degenerate density is computed
    as \code{Inf} at \code{loc} and 0 everywhere else.}
  \item{abstol}{numeric and non-negative. Absolute precision required. If \code{abstol = 0}, algorithm will run
        until total number of function evaluations exceeds Nmax (see also \code{Nmax}). If \eqn{n > 1}, algorithm runs until precision is reached for all \eqn{n} estimates for the         density.}
  \item{gam}{Monte Carlo confidence multiplier. Algorithm runs until  \eqn{estimated standard error < gam * abstol}.
      \code{gam = 3.3} (the default) means that one can expect that in 99.9 percent of the cases the actual absolute error is less than \code{abstol}.}
  \item{Nmax}{maximum number of function evaluations, can be used to
    control run time.}
  \item{B}{number of repetitions to get an error estimate in the
    randomized quasi-Monte Carlo approach.}
  \item{n_init}{size of the first point set being used to estimate
    the densities. Any positive integer allowed, powers or at least multiples of 2 are recommended for \code{method = sobol}.}
  \item{method}{Character string indicating method to be used. Allowed are "sobol", "ghalton" and "prng".}
  \item{log}{\code{\link{logical}} indicating whether the logarithm
    of the density is to be computed.}
  \item{\dots}{additional arguments containing parameters of
    mixing distributions when \code{mix} is a \code{\link{character}}
    string or \code{\link{function}}.}
}
\value{
  \code{dnvmix()} returns a list of length five: The first element is a \eqn{n}-\code{\link{vector}} with the (estimated)
  density values (default) or log-density values (if \code{log})
  of the specified multivariate normal variance mixture with
  location vector \code{loc} and scale matrix \code{scale}
  (a covariance matrix). The other elements are the number of iterations, the total
  number of function evaluations, an error estimate and the estimated variance of the randomized Quasi Monte Carlo estimator. If \eqn{n > 1}, the error estimate and the variance    reported is the largest error estimate and variance among the \eqn{n} estimates.

}
\details{
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given.

  The default factorization used is the Cholesky decomposition.
  To this end, \code{scale} needs to have full rank.

  Care should be taken when changing the algorithm-specific parameters, notably \code{B} and \code{Nmax}. Error estimates will not be reliable for too small \code{B}.

  If the absolute error tolerance \code{abstol} cannot be achieved with \code{Nmax} function evaluations, an additional warning will be returned.

  For the case of a multivariate normal or multivariate t distributions, user-friendly wrappers are provided in \code{\link{dNorm}} and \code{\link{dStudent}}. Note that in these cases, \code{dnvmix} uses the closed formulas for the density of a multivariate normal / t distribution.
}
\author{Marius Hofert and Erik Hintz}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{pnvmix}()}, \code{\link{rnvmix}()}, \code{\link{dNorm}()}, \code{\link{dStudent}()}
}
\examples{
## Generate a random correlation matrix in three dimensions
d <- 3
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
## Evaluate t_{3.5} density
df <- 3.5
x <- matrix(1:12/12, ncol = d) # evaluation points
dt <- dnvmix(x, mix = "inverse.gamma", df = df, scale = P)$Density
stopifnot(all.equal(dt, c(0.013266542, 0.011967156, 0.010760575, 0.009648682),
                    tol = 1e-7))

## The same can be achieved by defining a function for mix:
mix. <- function(u, df){
  df2 <- df/2
  return(1 / qgamma(u, shape = df2, rate = df2))
}
dt2 <-dnvmix(x, mix = mix., df = df, scale = P)$Density
stopifnot(all.equal(dt, dt2, tol = 5e-4))

## Evaluate normal density
dn <- dnvmix(x, mix = "constant", scale = P)$Density
stopifnot(all.equal(dn, c(0.013083858, 0.011141923, 0.009389987, 0.007831596),
                    tol = 1e-7))

## Missing data
x[3,2] <- NA
x[4,3] <- NA
dt <- dnvmix(x, mix = "inverse.gamma", df = df, scale = P)$Density
stopifnot(is.na(dt) == rep(c(FALSE, TRUE), each = 2))

## Univariate case
x <- matrix(1:10/10, ncol = 1)
dt <- dnvmix(x, mix = "inverse.gamma", df = df, factor = 1)$Density
dt. <- dt(as.vector(x), df = df)
stopifnot(all.equal(dt, dt.))
}
\keyword{distribution}