### pStudent() ###################################################################

##' @title Re-order Variables According to their Expected Integration Limits
##'        (Precondition) for the case of a multivariate t distribution. See [genzbretz2002, p. 957].
##' @param a vector, lower integration limit 
##' @param b vector, upper integration limit    
##' @param C Cholesky (lower triangular) factor of the sccale matrix R
##' @param R scale matrix 
##' @param q dimension of the problem
##' @param nu degrees of freedom 
##' @return list a, b, R, C of integration limits/scale matrix and Cholesky factor after reordering has been performed
##' @author Erik Hintz

precond_t <- function(a, b, R, C, q, nu)
{
  y <- rep(0, q-1)

  Rbar <- sqrt(nu) * gamma(nu / 2 ) / ( sqrt(2) * gamma( (nu+1) / 2 ) ) 
  
  ## Find i = argmin_j { <expected length of interval> }
  i <- which.min( apply( pnorm( cbind(a, b) / (Rbar * sqrt(diag(R)) ) ), 1, diff) )
  
  if(i != 1) {
    ## Swap 1 and i
    tmp <- swap(a = a, b = b, R = R, i = 1, j = i)
    a <- tmp$a
    b <- tmp$b
    R <- tmp$R
  }
  
  ## Store y1
  y[1] <- - ( dnorm(b[1]/Rbar) - dnorm(a[1]/Rbar)) / (pnorm(b[1]/Rbar) - pnorm(a[1]/Rbar))
  
  ## Update Cholesky
  C[1,1] <- sqrt( R[1,1] )
  C[2:q,1] <- as.matrix( R[2:q,1] / C[1,1] )
  
  for(j in 2:(q-1)){

    denom <- sqrt( diag(R)[j:q] - rowSums( as.matrix( C[j:q,1:(j-1)] )^2 ) )
    c <- as.matrix( C[j:q,1:j-1] ) %*% y[1:(j-1)]
    
    # Find i = argmin { <expected length of interval j> }
    i <- which.min( pnorm( (b[j:q]/Rbar - c) / denom ) - pnorm( (a[j:q]/Rbar - c) / denom ) ) + j - 1
    
    if(i != j){
      ## Swap i and j
      tmp <- swap(a = a, b = b,R = R,i = i, j =j)
      a <- tmp$a
      b <- tmp$b
      R <- tmp$R
      
      C[c(i,j),]    <- as.matrix(C[c(j,i),])
      C[j, (j+1):i] <- as.matrix(0, ncol = i-j, nrow = 1)
    }
    
    ## Update Cholesky
    C[j,j] <- sqrt(R[j,j] - sum( C[j,1:(j-1)]^2 ) )
    if(j< (q-1)) C[(j+1):q,j] <- ( R[(j+1):q,j] - as.matrix( C[(j+1):q,1:(j-1)] ) %*% C[j,1:(j-1)] ) / C[j,j]
    else C[(j+1):q,j] <- ( R[(j+1):q,j] - C[(j+1):q,1:(j-1)] %*% C[j,1:(j-1)] ) / C[j,j]
    
    ## Get yj
    ajbj <-  c( (a[j] / Rbar - C[j,1:(j-1)] %*% y[1:(j-1)] ) , (b[j] / Rbar - C[j,1:(j-1)] %*% y[1:(j-1)]) ) / C[j,j]
    y[j] <- ( dnorm(ajbj[1]) - dnorm(ajbj[2]) ) / ( pnorm(ajbj[2]) - pnorm(ajbj[1]) )
  }
  
  C[q,q] <- sqrt(R[q,q] - sum(C[q,1:(q-1)]^2))
  list(a = a, b = b, R = R, C = C)
}




##' @title Evaluatiing the tntegrand for pStudent. 
##' @param n sample size
##' @param skip number of points to be skipped (ignored unless method == "sobol")
##' @param a vector of lower limits
##' @param b vector of upper limits
##' @param C cholesky factor of scale
##' @param df degrees of freedom
##' @param ONE largest number x<1 such that x != 1
##' @param ZERO smallest number x>0 such that x!= 0
##' @param method character string indicating method to be used. Allowed are
##'         - "sobol" for a Sobol sequence.
##'         - "ghalton" for a generalized Halton sequence.
##'         - "prng" for a pure Monte Carlo approach. 
##' @return mean( f(U) ) where f is the integrand and U is the point-set specified by method. 
##' @author Erik Hintz

int_pStudent <- function(n, skip, a, b, C, q, df, ONE, ZERO, method){
  
  U <- switch(method,
              "sobol"   = {
                qrng::sobol(n = n, d = q, randomize = TRUE, skip = skip)
              },
              "gHalton" = {
                qrng::ghalton(n = n, d = q, method = "generalized")
              },
              "prng"    = {
                matrix(runif( n * q), ncol = q)
              })
  
  m <- .Call("eval_int_t_",
             n    = as.integer(n),
             q    = as.integer(q),
             U    = as.double(U),
             a    = as.double(a),
             b    = as.double(b),
             C    = as.double(C),
             nu   = as.double(df),
             ONE  = as.double(ONE),
             ZERO = as.double(ZERO))
  
  return(m)
}




##' @title Distribution Function of the Multivariate t Distribution
##' @param upper vector of upper limits
##' @param lower vector of lower limits
##' @param scale covariance matrix of dimension (q, q)
##' @param df degrees of freedom (positive real or Inf in which case the density
##'        of a N(loc, scale) is evaluated)
##' @param abserr numeric and non-negative. Absolute precision required. If abserr = 0, algorithm will run until total number of function evaluations exceeds Nmax (see also Nmax)
##' @param gam Monte Carlo confidence multiplier. Algorithm runs until gam * (estimated standard error) < abserr. gam = 3.3 means that one can expect 
##'        that in 99.9% of the cases the actual absolute error is less than abserr.
##' @param Nmax Total number of function evaluations allowed. 
##' @param N Number of randomizations to get error estimate. 
##' @param n_init First loop uses n_init function evaluations. Any positive integer allowed, powers or at least multiples of 2 are recommended. 
##' @param precond Logical. If TRUE (recommended), variable reordering as described in [genzbretz2002] pp. 955-956 is performed. Variable reordering can lead to a significant variance
##'        reduction and decrease in computational time. 
##' @param method Character string indicating method to be used. Allowed are
##'         - "sobol" for a Sobol sequence.
##'         - "ghalton" for a generalized Halton sequence.
##'         - "prng" for a pure Monte Carlo approach. 
##' @author Erik Hintz
pStudent <- function(upper, lower = rep(-Inf, length(upper)), scale, df, gam = 3.3, abserr = 0.001, Nmax = 1e8, N = 12, n_init = 2^5, precond = TRUE, method = "sobol")
{

  if(!is.matrix(scale)) scale <- as.matrix(scale)
  q <- dim(scale)[1] # dimension of the problem
  
  ## Checks
  if(method != "sobol" && method != "ghalton" && method != "prng") stop("Only sobol, ghalton and prng are allowed as methods")
  if( length(lower) != length(upper) ) stop("Lenghts of lower and upper differ")
  if( any(lower > upper) ) stop("lower needs to be smaller than upper (componentwise)")
  if( q != length(lower) ) stop("Dimension of scale does not match dimension of lower")
  
  ## Find infinite limits
  infina  <-  (lower == -Inf)
  infinb  <-  (upper == Inf)
  infinab <-  infina * infinb
  
  ## Remove double infinities
  if( sum(infinab) >0 )
  {
    whichdoubleinf <- which( infinab == TRUE)
    lower <- lower[ -whichdoubleinf ]
    lower <- lower[ -whichdoubleinf ]
    scale <- scale[ -whichdoubleinf, -whichdoubleinf ]
    ## Update dimension
    q <- dim(scale)[1]
  }
  
  ## Deal with the univariate case
  if(q == 1) return( pt(upper, df = df) - pt(lower, df = df) )
  
  ## Get Cholesky factor (lower triangular)
  C <- t(chol(scale))
  
  ## precondtioning (resorting the limits (cf precond_t); only for q > 2)
  if(precond && q>2) {
    temp <- precond_t(a = lower, b = upper, R = scale, C = C, q = q, nu = df)
    lower <- temp$a
    upper <- temp$b
    C <- temp$C
  }
  
  gam <- gam / sqrt(N) # instead of dividing sigma by sqrt(N) each time
  n. <- n_init # initial n
  T. <- rep(NA, N) # vector to safe RQMC estimates
  
  ONE <- 1-.Machine$double.neg.eps
  ZERO <- .Machine$double.eps
  
  if(method == "sobol") seed <- .Random.seed # need to reset to the seed later when sobol is being used.
  
  for(l in 1:N){
    T.[l] <- int_pStudent(n = n., skip = 0, a = lower, b = upper, C = C, q = q, df = df, ONE = ONE, ZERO = ZERO, method = method)
  }
  
  N. <- 2 * N * n. # N. will count the total number of function evaluations
  
  sig <- sd(T.)
  err <- gam * sig # note that this gam is actually gamma/sqrt(N)
  i. <- 0
  
  while(err > abserr && N. < Nmax)
  {
    if(method == "sobol") .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )
    
    for(l in 1:N){
      T.[l] <- ( T.[l] + int_pStudent(n = n., skip = n.,  a = lower, b = upper, C = C, q = q, df = df, ONE = ONE, ZERO = ZERO, method = method) )/2
      # note that T.[l] and int_pStudent(...) both depend on n. evaluations; hence we can just average them
    }

    N. <- N. + 2 * N * n. # update number of function evaluations. 
    n. <- 2 * n.
    sig <- sd(T.)
    err <- gam * sig # note that this gam is actually gamma/sqrt(N)
    i. <- i. + 1
  }
  T <- mean(T.)
  var <- (sig/sqrt(N))^2
  list(Prob = T, N = N., i = i., ErrEst = err, Var = var)
}




