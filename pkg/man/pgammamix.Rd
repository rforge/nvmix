\name{pgammamix}
\alias{pgammamix}
\title{Distribution Function of Gamma - mixture distribution}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including Student \emph{t} and normal distributions).
}
\usage{
pgammamix(m, qmix, d, lower.tail = TRUE, control = list(), verbose = TRUE, ...)
}
\arguments{
  \item{m}{\eqn{n}-\code{\link{vector}} of evaluation points.}
  \item{qmix}{specification of the mixing variable \eqn{W}; see McNeil et
    al. (2015, Chapter 6). Supported are the following types of
    specification (see also the examples below):
    \describe{\item{\code{\link{character}}:}{\code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1}) and
	\code{"inverse.gamma"} (in which case \eqn{W} is
	inverse gamma distributed with shape and rate parameters
	\code{df}/2).}
      \item{\code{\link{list}}:}{\code{\link{list}} of length at least
	one, where the first component is a \code{\link{character}}
	string specifying the base name of a distribution whose
	quantile function can be accessed via the prefix \code{"q"};
	an example is \code{"exp"} for which \code{"qexp"} exists.
	If the list is of length larger than one, the remaining elements
	contain additional parameters of the distribution; for \code{"exp"},
        for example, this can be the parameter \code{rate}.}
      \item{\code{\link{function}}:}{\code{\link{function}}
	interpreted as the quantile function of the mixing
	variable \eqn{W}.}
    }
  }
  \item{d}{degree of freedom parameter.}
  \item{lower.tail}{\code{\link{logical}}; if \code{TRUE} (default), probabilities are
  \eqn{P(X\le m)}, otherwise \eqn{P(X > m)}.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see ?pnvmix.}   
  \item{verbose}{% *not* the same as in ?dnvmix()
    \code{\link{logical}}
    indicating whether a warning is given if the required precision
    has not been reached.}
  \item{\dots}{additional arguments (for example, parameters) passed to the
    underlying mixing distribution when \code{qmix} is a
    \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  \code{pgammamix()} returns
  a \code{\link{numeric}} \eqn{n}-vector with the computed probabilities
  and corresponding attributes \code{"error"} (error estimates of the
  RQMC estimator) and \code{"numiter"} (number of iterations).
}
\details{
  TO COME - Not changed yet!!

  \code{pgammamix} evaluates the distribtution function of the mahalanobis distance
  of a multivariate normal variance mixture distribution. 
  
  
  One should highlight that evaluating normal variance mixtures
  is a non-trivial tasks which, at the time of development of
  \pkg{nvmix}, was not available in \R before, not even
  the special case of a multivariate Student \emph{t} distribution
  for non-integer degrees of freedom, which frequently appears in
  applications in finance, insurance and risk management after
  estimating such distributions.


  Internally, an iterative randomized Quasi-Monte Carlo (RQMC) approach is used
  to estimate the probabilities. It is an iterative algorithm that evaluates
  the integrand at a point-set (with size as specified by \code{increment})
  in each iteration until the pre-specified error tolerance \code{abstol}
  is reached. The attribute \code{"numiter"} gives the number of such
  iterations needed.
  
  Algorithm specific parameters (such as the above mentioned \code{'increment'}
  or \code{'abstol'}) can be passed as a list via \code{control}. 
  It can contain any of the following:
  
  \describe{
  \item{\code{method}}{\code{\link{character}} string indicating the method
    to be used to compute the integral. Available are:
    \describe{
      \item{\code{"sobol"}:}{Sobol' sequence (default).}
      \item{\code{"ghalton"}:}{generalized Halton sequence.}
      \item{\code{"PRNG"}:}{plain Monte Carlo based on a pseudo-random
         number generator.}
    }
  }
  \item{\code{mean.sqrt.mix}}{expectation of the square root \eqn{\sqrt(W)} of
    the mixing variable \eqn{W}. If \code{NULL}, it will be
    estimated via QMC; this is only needed for determining the reordering
    of the integration bounds, so a rather crude approximation
    is fine.}
  \item{\code{abstol}}{non-negative \code{\link{numeric}} providing the
    absolute precision required, defaults to \code{1e-3}. If \code{abstol = 0},
    the algorithm will typically run until the total number of function
    evaluations exceeds \code{fun.eval[2]}. If \eqn{n > 1} (so
    \code{x} has more than one row), the algorithm
    runs until the precision requirement is reached for all
    \eqn{n} density estimates.}
  \item{\code{precond}}{\code{\link{logical}} indicating whether preconditioning
    is applied, that is, reordering of the integration variables. If \code{TRUE},
    integration limits as well as \code{scale} are internally re-ordered in
    a way such that the overall variance of the integrand is usually smaller than
    with the original ordering; this usually leads smaller run-times.}
  \item{\code{increment}}{\code{\link{character}} string indicating how the
    sample size should be increased in each iteration. Available are:
    \describe{
      \item{\code{"doubling"}:}{next iteration has as many sample points
	as all the previous iterations combined.}
      \item{\code{"num.init"}:}{all iterations use an additional
	\code{fun.eval[1]}-many points.}
    }
  }  
  \item{\code{CI.factor}}{multiplier of the Monte Carlo confidence interval
    bounds.  The algorithm runs until \code{CI.factor} times the estimated
    standard error is less than \code{abstol}. If \code{CI.factor = 3.3}
    (the default), one can expect the actual absolute error to be less
    than \code{abstol} in 99.9\% of the cases.}
  \item{\code{fun.eval}}{\code{\link{numeric}(2)} providing the size of the
    first point set to be used to estimate the densities
    (typically a power of 2) and the maximal number of function
    evaluations. \code{fun.eval} defaults to \code{c(2^7, 1e8)}.}  
  \item{\code{max.iter.rqmc}}{\code{\link{numeric}}, providing the maximum number of 
    iterations allowed in the RQMC approach; the default is 15.}  
  \item{\code{B}}{number of randomizations for obtaining an error estimate in the
    randomized quasi-Monte Carlo (RQMC) approach; the default is 12.}
  }  
  
  
  Care should be taken when changing the algorithm-specific parameters,
  notably \code{method}, \code{precond}, \code{fun.eval[2]} and \code{B}.
  Error estimates will not be reliable for too small \code{B} and the
  performance of the algorithm depends heavily on the (quasi-)Monte
  Carlo point-set used.

  If the absolute error tolerance \code{abstol} cannot be achieved with
  \code{fun.eval[2]} function evaluations, an additional warning is
  thrown.

  \code{pStudent()} and \code{pNorm()} are wrappers of
  \code{pnvmix(, qmix = "inverse.gamma", df = df)} and
  \code{pnvmix(, qmix = "constant")}, respectively.
  In the univariate case, the functions
  \code{\link{pt}()} and \code{\link{pnorm}()} are used.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  McNeil, A. J., Frey, R. and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.

  Genz, A. and Bretz, F. (1999).
  Numerical computation of multivariate t-probabilities with
  application to power calculation of multiple contrasts.
  \emph{Journal of Statistical Computation and Simulation} 63(4),
  103--117.

  Genz, A. and Bretz, F. (2002).
  Comparison of methods for the computation of multivariate \emph{t}
  probabilities.
  \emph{Journal of Computational and Graphical Statistics} 11(4),
  950--971.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{rnvmix}()}
}
\examples{
### Examples for pnvmix() ######################################################
}
\keyword{distribution}