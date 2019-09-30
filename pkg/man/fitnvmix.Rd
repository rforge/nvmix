\name{fitnvmix}
\alias{fitnvmix}
\title{Fitting Multivariate Normal Variance Mixtures}
\description{
  Functionalities for fitting multivariate normal variance mixtures via an
  EMCE algorithm.
}
\usage{
fitnvmix(x, qmix, mix.param.bounds, nu.init = NA, 
         init.size.subsample = min(n, 100), size.subsample = n,
         control = list(), verbose = 2)
}
\arguments{
  \item{x}{\eqn{(n, d)}-data \code{\link{matrix}}.}
  \item{qmix}{specification of the mixing variable \eqn{W}; see McNeil et
    al. (2015, Chapter 6). Supported are the following types of
    specification (see also the examples below):
    \describe{\item{\code{\link{character}}:}{\code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1} and thus
	the multivariate normal distribution with mean vector
	\code{loc} and covariance matrix \code{scale} results), 
	\code{"inverse.gamma"} (in which case \eqn{W} is
	inverse gamma distributed with shape and rate parameters
	\code{df}/2 and thus the multivariate
	Student \emph{t} distribution with \code{df} degrees of freedom
	results) and \code{"pareto"} (in which case \eqn{W} is
	Pareto distributed with scale equal to unity and shape equal to
	\code{alpha}).}
      \item{\code{\link{function}}:}{\code{\link{function}}
	interpreted as the quantile function of the mixing
	variable \eqn{W}. In this case, \code{qmix} \emph{must} have the form \code{qmix = function(u, nu)}, where the argument \code{nu} corresponds to the parameter (vector) specifying the distribution
	of the mixing variable.}
    }
  }
  \item{mix.param.bounds}{eiher \code{\link{numeric}(2)} or a
  \code{\link{matrix}} with two columns. The first/second column
  corresponds to the lower/upper bound of \eqn{nu_i}, the ith
  component of the parameter vector \eqn{nu} of the mixing variable \eqn{W}. All elements need to be finite, numeric values. 
  Note: The algorithm tends to converge quicker if the parameter
  ranges supplied are smaller.}
  \item{nu.init}{either \code{NA} or an initial guess for the parameter (vector) \eqn{nu}. In the former case an initial
  estimate is calculated by the algorithm. If \code{nu.init} is provided, the algorithm often converges faster; the better the starting value, the faster convergence.}
  \item{init.size.subsample}{\code{\link{numeric}}, non-negative, giving the sub-samplesize used to get an initial estimate for \eqn{nu}. Only used if \code{is.na(nu.init)}, otherwise ignored.}
  \item{size.subsample}{\code{\link{numeric}}, non-negative, specifying the size of the subsample that is being used in the EMCE iterations to optimize the log-likelihood over \eqn{nu}. Defaults to \code{n}, so that the full sample is being used. Decreasing this number can lead to faster run-times (as fewer log-densities need to be estimated) but also to an increase in bias and variance.}
  \item{control}{\code{\link{list}} specifying algorithm specific parameters; see below under 'Details' and \code{?get.set.parameters}.}
  \item{verbose}{\code{\link{numeric}} or \code{\link{logical}} (in which case it is converted to \code{numeric}) specifying the
  amount of tracing to be done. If 0 or \code{FALSE}, neither tracing nor warnigns are communicated; if 1, only warnigns are communicated, if 2 or 3, warnings and (shorter or longer) tracing information is provided.}
}
\value{
  \code{fitnvmix()} by returns a \code{\link{list}} containing \code{nu}, 
  \code{loc}, \code{scale} as MLEs for \eqn{nu}, the location vector and the scale 
  matrix, respectively; the list also contains \code{iter} (number of EMCE iterations performed)
  and \code{max.ll} (log-likelihood at the MLEs).
}
\details{
  The function \code{fitnvmix} uses an EMCE algorithm to approximate the MLEs of
  the parameters \code{nu}, \code{loc} and \code{scale} of a normal variance mixture specified by
  \code{qmix}. The underlying procedure successively estimates \code{nu} (with 
  given \code{loc} and \code{scale}) by maximizing the likelihood which is approximated
  by \code{\link{dnvmix}()} (unless \code{qmix} is a \code{character} string, in which
  case analytical formulas for the log-densities are used) and \code{scale} and \code{loc}
  (given \code{nu}) using weights (which again need to be approximated) related
  to the posterior distribution, details can be found in (...).
  
  It should be highlighted that (unless unless \code{qmix} is a \code{character} string),
  every log-likelihood and every weight needed in the estimation is numerically
  approximated via RQMC methods. For large dimensions and sample sizes this procedure
  can therefore be slow. 
  
  Various tolerances and convergence criteria can be changed by the user via the
  \code{control} argument. For more details, see \code{\link{get.set.parameters}()}.
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
  \code{\link{dnvmix}()}, \code{\link{rnvmix}()}, \code{\link{pnvmix}()}
}
\examples{
set.seed(274) # for reproducibility

## For sampling:
nu                <- 3 # parameter used to sample data
d                 <- 5 # dimension
n                 <- 100 # small sample size to have examples run fast 
loc               <- rep(0, d) # location vector
A                 <- matrix(runif(d * d), ncol = d)
diag_vars         <- diag(runif(d, min = 2, max = 5))
scale             <- diag_vars \%*\% cov2cor(A \%*\% t(A)) \%*\% diag_vars # scale matrix
mix.param.bounds  <- c(1, 5) # nu in [1, 5]

### Example 1: Fitting a multivariate t distribution ###########################

## Define 'qmix' as the quantile function of an IG(nu/2, nu/2) distribution
qmix <- function(u, nu){
   1 / qgamma(1 - u, shape = nu/2, rate = nu/2)
}
## Sample data using 'rnvmix':
x <- rnvmix(n, qmix = qmix, nu = nu, loc = loc, scale = scale)
## Call 'fitvnmix' with 'qmix' as a function (so all densities/weights are estimated)
(MyFit11 <- fitnvmix(x, qmix = qmix, mix.param.bounds = mix.param.bounds))
## Call 'fitnvmix' with 'qmix = "inverse.gamma"' in which case analytical formulas
## for weights and densities are used:
(MyFit12 <- fitnvmix(x, qmix = "inverse.gamma", mix.param.bounds = mix.param.bounds))
## Check
stopifnot(all.equal(MyFit11$nu, MyFit12$nu, tol = 5e-2))


## Visual goodness-of-fit test: QQ Plot of mahalanobis distances obtained
## from fitted parameters using 'qqplot.maha()':
qqplot.maha(x, qmix = "inverse.gamma", loc = MyFit11$loc, scale = MyFit11$scale, 
            df = MyFit11$nu)

\dontrun{
### Example 2: Fitting a Pareto mixture ########################################
## Define 'qmix' as the quantile function of a Par(nu, 1) distribution
qmix <- function(u, nu){
   (1-u)^(-1/nu)
}
## Sample data using 'rnvmix':
x <- rnvmix(n, qmix = qmix, nu = nu, loc = loc, scale = scale)
## Call 'fitvnmix' with 'qmix' as function (so all densities/weights are estimated)
(MyFit21 <- fitnvmix(x, qmix = qmix, mix.param.bounds = mix.param.bounds))
## Call 'fitnvmix' with 'qmix = "pareto"' in which case an analytical formula
## for the density is used:
(MyFit22 <- fitnvmix(x, qmix = "pareto", mix.param.bounds = mix.param.bounds))
stopifnot(all.equal(MyFit21$nu, MyFit22$nu, tol = 5e-2))}
}
\keyword{distribution}