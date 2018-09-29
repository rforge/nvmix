\name{pnvmix}
\alias{pnvmix}
\alias{pStudent}
\alias{pNorm}
\title{Distribution Function of Multivariate Normal Variance Mixtures}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including Student \emph{t} and normal distributions).
}
\usage{
pnvmix(upper, lower = matrix(-Inf, nrow = n, ncol = d), mix, mean.sqrt.mix = NULL,
       loc = rep(0, d), scale = diag(d), standardized = FALSE,
       method = c("sobol", "ghalton", "PRNG"), precond = TRUE,
       abstol = 1e-3, CI.factor = 3.3, fun.eval = c(2^6, 1e8),
       increment = c("doubling", "num.init"), B = 12, ...)

pStudent(upper, lower = rep(-Inf, d),
         df, loc = rep(0, d), scale = diag(d), standardized = FALSE,
         method = c("sobol", "ghalton", "PRNG"), precond = TRUE,
         abstol = 1e-3, CI.factor = 3.3, fun.eval = c(2^6, 1e8), B = 12)
pNorm(upper, lower = rep(-Inf, d),
      loc = rep(0, d), scale = diag(d), standardized = FALSE,
      method = c("sobol", "ghalton", "PRNG"), precond = TRUE,
      abstol = 1e-3, CI.factor = 3.3, fun.eval = c(2^6, 1e8), B = 12)
}
\arguments{
  \item{upper}{\eqn{(n, d)}-\code{\link{matrix}} of upper integration limits;
    each row represents a \eqn{d}-vector of upper integration limits.}
  \item{lower}{\eqn{(n, d)}-\code{\link{matrix}} of lower integration limits
    (componentwise less than or equal to \code{upper});
    each row represents a \eqn{d}-vector of lower integration limits.}
  \item{mix}{specification of the mixing variable \eqn{W}; see McNeil et
    al. (2015, Chapter 6). Supported are the following types of
    specification (see also the examples below):
    \describe{
      \item{\code{\link{character}}:}{\code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1} and thus
	the multivariate normal distribution with mean vector
	\code{loc} and covariance matrix \code{scale} results) and
	\code{"inverse.gamma"} (in which case \eqn{W} is
	inverse gamma distributed with shape and rate parameters
	\code{df}/2 and thus the multivariate
	Student \emph{t} distribution with degrees of freedom
	\code{df} (required to be provided via
	the ellipsis argument) results); additional arguments
	(such as parameters) can be passed via the ellipsis argument.}
      \item{\code{\link{list}}:}{\code{\link{list}} of length at least
	one, where the first component is a \code{\link{character}}
	string specifying the base name of a distribution whose
	quantile function can be accessed via the prefix \code{"q"};
	an example is \code{"exp"} for which \code{"qexp"} exists.
	If the list is of length larger than one, the remaining elements
	contain additional parameters of the distribution; for \code{"exp"},
        this can be the parameter \code{rate}.}
      \item{\code{\link{function}}:}{\code{\link{function}}
	interpreted as the quantile function of the mixing
	variable \eqn{W}; additional arguments (such as parameters)
        can be passed via the ellipsis argument.}
    }
  }
  \item{mean.sqrt.mix}{expectation of the square root \eqn{\sqrt(W)} of
    the mixing variable \eqn{W}. If \code{NULL}, it will be
    estimated via QMC; this is only needed for determining the reordering
    of the integration bounds, so a rather crude approximation
    is fine.}
  \item{df}{degress of freedom (positive real or \code{Inf} in which case
    the distribution is interpreted as multivariate normal with mean
    vector \code{loc} and covariance matrix \code{scale}).}
  \item{loc}{location vector of dimension \eqn{d}; this equals the mean
    vector of a random vector following the specified normal variance
    mixture distribution if (and only if) the latter exists.}
  \item{scale}{scale matrix (a covariance matrix entering the
    distribution as a parameter) of dimension \eqn{(d, d)};
    this equals the covariance matrix of a random vector following
    the specified normal variance mixture distribution divided by
    the expecation of the mixing variable \eqn{W} if (and only if) the
    former exists.}
  \item{standardized}{\code{logical} indicating whether \code{scale} is
    assumed to be a correlation matrix.}
  \item{method}{\code{\link{character}} string indicating the method
    to be used to compute the integral. Available are:
    \describe{
      \item{\code{"sobol"}:}{}
      \item{\code{"ghalton"}:}{}
      \item{\code{"PRNG"}:}{}
    }
  }
  \item{precond}{\code{\link{logical}} indicating whether preconditioning
    is applied, that is, reordering of the integration variables.}
  \item{abstol}{non-negative \code{\link{numeric}} providing the
    absolute precision required; if \code{abstol} is zero,
    the algorithm will run until the total number of function
    evaluations exceeds \code{fun.eval[2]}.}
  \item{CI.factor}{multiplier of the Monte Carlo confidence interval
    bounds.  The algorithm runs until \code{CI.factor} times the estimated
    standard error is less than \code{abstol}. If \code{CI.factor = 3.3}
    (the default), one can expect the actual absolute error to be less
    than \code{abstol} in 99.9\% of the cases.}
  \item{fun.eval}{\code{\link{numeric}(2)} providing the size of the
    first point set to be used to estimate the probabilities
    (typically a power of 2) and the maximal number of function
    evaluations.}
  \item{increment}{\code{\link{character}} string indicating how the
    sample size should be increased in each iteration. Available are:
    \describe{
      \item{\code{"doubling"}:}{next iteration has as many sample points
	as all the previous iterations combined.}
      \item{\code{"num.init"}:}{all iterations use an additional
	\code{fun.eval[1]} many points.}
    }
  }
  \item{B}{number of randomizations to get an error estimate in the
    randomized quasi-Monte Carlo approach.}
  \item{\dots}{additional arguments (e.g., parameters) passed to the
    underlying mixing distribution when \code{mix} is a
    \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  \code{pnvmix()}, \code{pStudent()} and \code{pNorm()} return
  a \code{\link{numeric}} vector with the computed probabilities
  and attributes \code{"error"} (error estimate of the RQMC estimator)
  and  \code{"numiter"} (number of iterations).
}
\details{
  One should highlight that evaluating normal variance mixtures
  in general and the multivariate Student \emph{t} distribution
  for non-integer degrees of freedom in particular are non-trivial
  tasks which, at the time of development of \pkg{nvmix}, was
  not available in \R before.

  Note that the procedures call underlying C code. Currently,
  dimensions \eqn{d\ge 16510}{d >= 16510} are not supported for the
  default method \code{sobol}.

  Care should be taken when changing the algorithm-specific parameters,
  notably \code{method}, \code{precond}, \code{fun.eval[2]} and \code{B}.
  Error estimates will not be reliable for too small \code{B} and the
  performance of the algorithm depends heavily on the (quasi-)Monte
  Carlo point-set used.

  If the absolute error tolerance \code{abstol} cannot be achieved with
  \code{fun.eval[2]} function evaluations, an additional warning will be
  returned.

  \code{pStudent()} and \code{pNorm()} are wrappers of
  \code{pnvmix(, mix = "inverse.gamma", df = df)} and
  \code{pnvmix(, mix = "constant")}, respectively.
  In the univariate case, the functions
  \code{\link{pt}()} and \code{\link{pnorm}()} are used.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{rnvmix}()}
}
\examples{
### Examples for pnvmix() ######################################################

## Generate a random correlation matrix in d dimensions
d <- 3
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))

## Evaluate a t_{1/2} distribution function
a <- -3 * runif(d) * sqrt(d) # random lower limit
b <-  3 * runif(d) * sqrt(d) # random upper limit
df <- 0.5 # note that this is *non-integer*
set.seed(1)
pt1 <- pnvmix(b, lower = a, mix = "inverse.gamma", df = df, scale = P)

## Here is a version providing the quantile function of the mixing distribution
qW <- function(u, df) 1 / qgamma(u, shape = df/2, rate = df/2)
mean.sqrt.mix <- sqrt(df) * gamma(df/2) / (sqrt(2) * gamma((df+1) / 2))
set.seed(1)
pt2 <- pnvmix(b, lower = a, mix = qW, mean.sqrt.mix = mean.sqrt.mix, df  = df,
              scale = P)

## Compare
stopifnot(all.equal(pt1, pt2, tol = 7e-4, check.attributes = FALSE))

## mean.sqrt.mix will be approximated by QMC internally if not provided,
## so the results will differ slightly.
set.seed(1)
pt3 <- pnvmix(b, lower = a, mix = qW, df = df, scale = P)
stopifnot(all.equal(pt3, pt1, tol = 7e-4, check.attributes = FALSE))

## Case with missing data and a matrix of lower and upper bounds
a. <- matrix(rep(a, each = 4), ncol = d)
b. <- matrix(rep(b, each = 4), ncol = d)
a.[2,1] <- NA
b.[3,2] <- NA
pt <- pnvmix(b., lower = a., mix = "inverse.gamma", df = df, scale = P)
stopifnot(is.na(pt) == c(FALSE, TRUE, TRUE, FALSE))


### Examples for pStudent() and pNorm() ########################################

## Evaluate a t_{3.5} distribution function
set.seed(271)
pt <- pStudent(b, lower = a, df = 3.5, scale = P)
stopifnot(all.equal(pt, 0.6180, tol = 5e-5, check.attributes = FALSE))

## Evaluate a normal distribution function
set.seed(271)
pn <- pNorm(b, lower = a, scale = P)
stopifnot(all.equal(pn, 0.7001, tol = 1e-4, check.attributes = FALSE))
}
\keyword{distribution}