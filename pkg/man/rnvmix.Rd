\name{rnvmix}
\alias{rnvmix}
\title{Random Number Generator for Multivariate Normal Variance Mixtures}
\description{
  Generate vectors of random variates from multivariate normal variance
  mixtures (including normal and Student \emph{t}).
}
\usage{
rnvmix(n, mix, loc = rep(0, d), scale, factor = factorize(scale), ...)
}
\arguments{
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{mix}{specification of the mixing variable \eqn{W}; see McNeil et
    al. (2015). Supported are the following types of specification (see
    also the examples below):
    \describe{
      \item{\code{\link{character}}:}{a \code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1} and thus a sample
	 from the multivariate normal distribution with mean vector
	 \code{loc} and covariance matrix \code{scale} results) and
	 \code{"inverse.gamma"} (in which case \eqn{W} is drawn from
	 and inverse gamma distribution with shape and rate parameters
	 \code{df}/2 resulting in a sample from the multivariate
	 Student \emph{t} distribution with degrees of freedom
	 \code{df}; note that \code{df} needs to be provided via
	 the ellipsis argument then; see the examples below).}
      \item{\code{\link{list}}:}{a \code{\link{list}} of length at least
	one, where the first component is a \code{\link{character}}
	string specifying the base name of a distribution which can be
	sampled via prefix \code{"r"}; an example is \code{"exp"}
        for which \code{"rexp"} exists for sampling. If the list is
        of length larger than one, the remaining elements contain
        additional parameters of the distribution; for \code{"exp"},
        this can be the parameter \code{rate}.}
      \item{\code{\link{function}}:}{a \code{\link{function}}
	interpreted as the quantile function of the mixing
	variable \eqn{W}; internally, sampling is then done with the
	inversion method by applying this function to U(0,1) random variates.}
      \item{\code{\link{numeric}}:}{a \code{\link{numeric}} vector
	of length \code{n} providing a random sample of the mixing
	variable \eqn{W}.}
    }
  }
  \item{loc}{location vector of dimension \eqn{d}.}
  \item{scale}{covariance matrix of dimension \eqn{(d, d)}.}
  \item{factor}{factorization matrix of the covariance matrix
    \code{scale}. A \code{\link{matrix}} \eqn{R} with \eqn{d}
    rows such that \eqn{R^T R} equals \code{scale}. \eqn{R} is internally
    multiplied from the right to an \eqn{(n, d)}-matrix of independent
    standard normals to obtain a sample from a multivariate normal
    with zero mean vector and covariance matrix \code{scale}.}
  \item{\dots}{additional arguments containing parameters of
    mixing distributions when \code{mix} is a \code{\link{character}}
    string.}
}
\value{
  \code{rnvmix()} returns an \eqn{(n, d)}-\code{\link{matrix}}
  containing \eqn{n} samples of the specified (via \code{mix})
  \eqn{d}-dimensional multivariate normal variance mixture with location
  vector \code{loc} and scale matrix \code{scale} (a covariance matrix).
}
\details{
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given.

  The default factorization used to obtain \code{factor} is the Cholesky
  decomposition via \code{\link{chol}()}. To this end, \code{scale}
  needs to have full rank. See \code{\link{factorize}()} for other options.

  The number of rows of \code{factor} are the dimension \eqn{d} of
  the sample. Typically (but not necessarily), \code{factor} is square.
}
\author{Marius Hofert}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{pnvmix}()}, \code{\link{factorize}()}
}
\examples{
## Generate a random correlation matrix in three dimensions
d <- 3
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))

## Draw random variates and compare
df <- 3.5
n <- 1000
set.seed(271)
X <- rnvmix(n, mix = "inverse.gamma", df = df, scale = P) # providing scale
set.seed(271)
X. <- rnvmix(n, mix = "inverse.gamma", df = df, factor = chol(P)) # providing the factor
stopifnot(all.equal(X, X.))

## Checking 'df = Inf'
set.seed(271)
X <- rnvmix(n, mix = "constant", scale = P)
set.seed(271)
X. <- rnvmix(n, mix = "inverse.gamma", df = Inf, scale = P)
stopifnot(all.equal(X, X.))

## Univariate case (dimension = number of rows of 'factor' = 1 here)
set.seed(271)
X.1d <- rnvmix(n, mix = "inverse.gamma", df = df, factor = 1/2)
set.seed(271)
X.1d. <- rnvmix(n, mix = "inverse.gamma", df = df, factor = 1)/2 # manual scaling
stopifnot(all.equal(X.1d, X.1d.))

## Checking different ways of providing 'mix'
## 1) By providing a character string (and corresponding '...' arguments)
set.seed(271)
X.mix1 <- rnvmix(n, mix = "inverse.gamma", df = df, factor = P)
## 2) By providing a list; the first element has to be an existing distribution
##    with random number generator available with prefix "r"
rinverse.gamma <- function(n, df) 1 / rgamma(n, shape = df/2, rate = df/2)
set.seed(271)
X.mix2 <- rnvmix(n, mix = list("inverse.gamma", df = df), factor = P)
## 3) The same without extra arguments (need the extra list() here to
##    distinguish from Case 1))
rinverseGamma <- function(n) 1 / rgamma(n, shape = df/2, rate = df/2)
set.seed(271)
X.mix3 <- rnvmix(n, mix = list("inverseGamma"), factor = P)
## 4) By providing a quantile function
##    Note: P(1/Y <= x) = P(Y >= 1/x) = 1-F_Y(1/x) = y <=> x = 1/F_Y^-(1-y)
set.seed(271)
X.mix4 <- rnvmix(n, mix = function(p) 1 / qgamma(1-p, shape = df/2, rate = df/2),
                 factor = P)
## 5) By providing random variates
set.seed(271) # if seed is set here, results are comparable to the above methods
W <- rinverse.gamma(n, df = df)
X.mix5 <- rnvmix(n, mix = W, factor = P)
## Compare
## Note that X.mix4 is not 'all equal' with X.mix1 (or the other samples)
## since rgamma() != qgamma(runif()) (or qgamma(1-runif()))
stopifnot(all.equal(X.mix2, X.mix1),
          all.equal(X.mix3, X.mix1),
          all.equal(X.mix5, X.mix1))
}
\keyword{distribution}