\name{qnvmix1}
\alias{qnvmix1}
\title{Quantile Function of univariate Normal Variance Mixture Distribution}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t} for non-integer degrees of freedom).
}
\usage{
qnvmix1(u, shift = 0, scale = 1, mix, N = 1e5, ... )
}
\arguments{
  \item{u}{Probability.}
  \item{shift}{shift vector of length \eqn{d}. If \code{mix} has a mean, this is the mean of the normal variance mixture distribution.}
  \item{scale}{positive scale (variance of the normal distribution.}
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
	inversion method by applying this function to U(0,1) random variates.}
    }
  }
  \item{N}{number of Sobol points used to approximate the univariate cdf.}
  \item{\dots}{additional arguments containing parameters of
    mixing distributions when \code{mix} is a \code{\link{character}}
    string.}
}
\value{
  \code{qnvmix1()} returns a number \eqn{q} satisfying \eqn{q = inf_x { F(x) \geq u}} where \eqn{F(x)} denotes the univariate cdf of the normal variance mixture.  
}
\details{
  Just a workhorse, needs to be improved. 
}
\author{Marius Hofert, Erik Hintz and Christiane Lemieux}
\references{
  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\examples{
## Consider a t distribution
u <- 0.4
df <- 1
stopifnot(all.equal(qt(u, df = df), qnvmix1(u, mix = "inverse.gamma", df = df), tol = 5e-4))
}
\keyword{distribution}