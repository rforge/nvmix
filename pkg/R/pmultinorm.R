### pmultinorm() ###################################################################

##' @title Re-order Variables According to their Expected Integration Limits
##'        (Precondition) for the case of a multivariate normal distribution.
##' @param a vector, lower integration limit 
##' @param b vector, upper integration limit    
##' @param C Cholesky (lower triangular) factor of the sccale matrix R
##' @param R scale matrix 
##' @param q dimension of the problem
##' @return list a, b, R, C of integration limits/scale matrix and Cholesky factor after reordering has been performed
##' @author Erik Hintz


precond_norm <- function(a, b, R, C, q, nu)
{
  y <- rep(0, q - 1)
  
  ## Find i = argmin_j { <expected length of interval> }
  i <- which.min( apply( pnorm( cbind(a, b) / sqrt(diag(R)) ), 1, diff) )
  
  if(i != 1) {
    ## Swap 1 and i
    tmp <- swap(a = a, b = b, R = R, i = 1, j = i)
    a <- tmp$a
    b <- tmp$b
    R <- tmp$R
  }
  
  ## Store y1
  y[1] <- - ( dnorm( b[1] ) - dnorm( a[1] ) ) / ( pnorm( b[1] ) - pnorm( a[1] ))
  
  ## Update Cholesky
  C[1,1] <- sqrt( R[1,1] )
  C[2:q,1] <- as.matrix( R[2:q,1] / C[1,1] )
  
  for(j in 2:(q-1)){
    
    denom <- sqrt( diag(R)[j:q] - rowSums( as.matrix( C[j:q,1:(j-1)] )^2 ) )
    c <- as.matrix( C[j:q,1:j-1] ) %*% y[1:(j-1)]
    
    # Find i = argmin { <expected length of interval j> }
    i <- which.min( pnorm( (b[j:q] - c) / denom ) - pnorm( (a[j:q] - c) / denom ) ) + j - 1
    
    if(i != j){
      ## Swap i and j
      tmp <- swap(a = a, b = b, R = R, i = i, j =j)
      a <- tmp$a
      b <- tmp$b
      R <- tmp$R
      
      C[c(i,j),]    <- as.matrix(C[c(j,i),])
      C[j, (j+1):i] <- as.matrix(0, ncol = i-j, nrow = 1)
    }
    
    ## Update Cholesky
    C[j,j] <- sqrt( R[j,j] - sum( C[j,1:(j-1)]^2 ) )
    if(j < (q-1) ) C[(j+1):q,j] <- ( R[(j+1):q,j] - as.matrix( C[(j+1):q,1:(j-1)] ) %*% C[j,1:(j-1)] ) / C[j,j]
    else C[(j+1):q,j] <- ( R[(j+1):q,j] - C[(j+1):q,1:(j-1)] %*% C[j,1:(j-1)] ) / C[j,j]
    
    ## Get yj
    ajbj <-  c( (a[j] - C[j,1:(j-1)] %*% y[1:(j-1)] ) , (b[j] - C[j,1:(j-1)] %*% y[1:(j-1)]) ) / C[j,j]
    y[j] <- ( dnorm(ajbj[1]) - dnorm(ajbj[2]) ) / ( pnorm(ajbj[2]) - pnorm(ajbj[1]) )
  }
  
  C[q,q] <- sqrt(R[q,q] - sum(C[q,1:(q-1)]^2))
  list(a = a, b = b, R = R, C = C)
}






##' @title Distribution Function of the Multivariate Normal Distribution
##' @param upper vector of upper limits
##' @param lower vector of lower limits
##' @param scale covariance matrix of dimension (q, q)
##' @param standardized logical. If TRUE, scale is assumed to be a correlation matrix; if FALSE (default), lower, upper and scale will be normalized.
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
pmultinorm <- function(upper, lower = rep(-Inf, length(upper)), mean = rep(0, length(upper)), scale, standardized = FALSE, gam = 3.3, abserr = 0.001, Nmax = 1e8, N = 12, n_init = 2^5, precond = TRUE, method = "sobol")
{
  
  if(!is.matrix(scale)) scale <- as.matrix(scale)
  q <- dim(scale)[1] # dimension of the problem
  
  ## Checks
  if( method != "sobol" && method != "ghalton" && method != "prng") stop("Only sobol, ghalton and prng are allowed as methods")
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
  
  ## Subtract mean if necessary: 
  if(any(mean != 0)){
    lower <- lower - mean
    upper <- upper - mean
  }
  
  ## Standardize if necessary:
  if(!standardized){
    Dinv <- diag(1/sqrt(diag(scale)))
    scale <- Dinv %*% scale %*% Dinv
    lower <- as.vector(Dinv %*% lower)
    upper <- as.vector(Dinv %*% upper)
  }
  
  
  ## Deal with the univariate case
  if(q == 1) return( pnorm(upper) - pnorm(lower) )
  
  ## Get Cholesky factor (lower triangular)
  C <- t(chol(scale))
  
  ## precondtioning (resorting the limits (cf precond_t); only for q > 2)
  if(precond && q>2) {
    temp <- precond_norm(a = lower, b = upper, R = scale, C = C, q = q)
    lower <- temp$a
    upper <- temp$b
    C <- temp$C
  }
  
  gam <- gam / sqrt(N) # instead of dividing sigma by sqrt(N) each time
  n. <- n_init # initial n
  T. <- rep(NA, N) # vector to store RQMC estimates
  
  ONE <- 1-.Machine$double.neg.eps
  ZERO <- .Machine$double.eps
  
  if(method == "sobol") seed <- .Random.seed # need to reset to the seed later if a Sobol sequence is being used.
  
  err <- abserr + 42 # initialize err to something bigger than abserr so that we can enter the while loop
  N. <- 0 # N. will count the total number of function evaluations
  i. <- 0 # initialize counter; this will count the number of iterations in the while loop
  
  while(err > abserr && N. < Nmax)
  {
    
    ## First iteration seperate since it does not have a skip and no averaging is being done here:
    if(i. == 0){
      
      ## Get N RQCM estimates
      for(l in 1:N){
        
        ## Get the pointset
        U <- switch(method,
                    "sobol"   = {
                      qrng::sobol(n = n., d = q - 1, randomize = TRUE, skip = 0)
                    },
                    "gHalton" = {
                      qrng::ghalton(n = n., d = q - 1, method = "generalized")
                    },
                    "prng"    = {
                      matrix(runif( n. * (q - 1)), ncol = q - 1)
                    })
        
        ## Evaluate the integrand at the point-set and save it in T[] (calls C Code)
        T.[l] <- .Call("eval_int_normal_",
                       n    = as.integer(n.),
                       q    = as.integer(q),
                       U    = as.double(U),
                       a    = as.double(lower),
                       b    = as.double(upper),
                       C    = as.double(C),
                       ONE  = as.double(ONE),
                       ZERO = as.double(ZERO))
      }
      
      sig <- sd(T.) # get standard deviation of the estimator 
      err <- gam * sig # update error. Note that this gam is actually gamma/sqrt(N)
      N. <- 2 * N * n. # update the total number of function evaluations; mutliplied by 2 since antithetic variates are being used in eval_int_t_ 
      i. <- i. + 1 # update counter
      
    } else {
      
      ## In this case, we already had one iteration:
      
      if(method == "sobol") .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )
      
      for(l in 1:N){
        
        ## Get the next pointset (observe that we skip the first n. points). Since we did reset the seed, we will get the correct skipped pointset in each iteration, 
        ## since N did not change. 
        
        U <- switch(method,
                    "sobol"   = {
                      qrng::sobol(n = n., d = q - 1, randomize = TRUE, skip = n.)
                    },
                    "gHalton" = {
                      qrng::ghalton(n = n., d = q - 1, method = "generalized")
                    },
                    "prng"    = {
                      matrix(runif( n. * q), ncol = q - 1)
                    })
        
        ## Evaluate the integrand at the point-set and update T[l]: The estimator saved in T[l] and the new estimator based on U both depend on the same number
        ## of points, namely n., so we can just average them:
        T.[l] <- (T.[l] + .Call("eval_int_normal_",
                                n    = as.integer(n.),
                                q    = as.integer(q),
                                U    = as.double(U),
                                a    = as.double(lower),
                                b    = as.double(upper),
                                C    = as.double(C),
                                ONE  = as.double(ONE),
                                ZERO = as.double(ZERO)) )/2
      }
      
      N. <- N. + 2 * N * n. # update number of function evaluations. 
      n. <- 2 * n. # double sample size 
      sig <- sd(T.) # get standard deviation of the estimator 
      err <- gam * sig # note that this gam is actually gamma/sqrt(N)
      i. <- i. + 1 # update counter 
    }
    
  }
  
  ## Calculate the RQMC estimator:
  T <- mean(T.)
  
  ## Get the variance 
  var <- (sig/sqrt(N))^2
  
  ## And return
  list(Prob = T, N = N., i = i., ErrEst = err, Var = var)
}