### pnvmix() ###################################################################

##' @title Re-order Variables According to their Expected Integration Limits
##'        (Precondition). See [genzbretz2002, p. 957]
##' @param C Cholesky (lower triangular) factor of the covariance matrix R
##' @param q dimension of the problem
##' @param ... all the other parameters same as in pnvmix
##' @return list a, b, R, C of integration limits/covariance matrix/...
##'         after reordering has been performed
##' @return d vector indicating the new ordering
##' @author Erik Hintz
##' @note the outermost integral has the smallest limits, etc.
precond_t <- function(a, b, R, C, q, nu)
{
    y <- rep(0, q-1)
    sumysq <- 0
    d <- 1:q

    ## Find i = argmin_j { <expected length of interval> }
    i <- which.min( apply( pt( cbind(a, b) / sqrt(diag(R) ), nu), 1, diff) )

    if(i != 1) {
        ## Swap 1 and i
        tmp <- swap(a,b,R,1,i)
        a <- tmp$a
        b <- tmp$b
        R <- tmp$R
        d[c(1,i)] <- d[c(i,1)]
    }

    ## Store y1
    y[1] <- gamma( (nu+1)/2 )/(gamma( nu/2 ) * sqrt( nu*pi ) ) * nu / (nu+3) *
        ( (1+ a[1]^2/2)^(-(nu+3)/2 ) - (1+b[1]^2/2)^(-(nu+3)/2)) / ( pt(b[1],nu) - pt(a[1],nu) )
    ## Initialize sum of y^2
    sumysq <- y[1]^2

    ## Update Cholesky
    C[1,1] <- sqrt(R[1,1])
    C[2:q,1] <- as.matrix(R[2:q,1]/C[1,1])

    for(j in 2:(q-1)){
        ## c0 <- rowSums(as.matrix(C[j:q,1:(j-1)])^2)
        c01<- sqrt( (nu+j-1)/(nu+sumysq))
        c1 <- c01/sqrt(diag(R)[j:q]-rowSums(as.matrix(C[j:q,1:(j-1)])^2))
        c2 <- as.matrix(C[j:q,1:j-1]) %*%y [1:(j-1)]

        ## Find i = argmin { <expected length of interval j> }
        i <- which.min(pt(c1 * (b[j:q] - c2), nu+j-1) - pt(c1 * (a[j:q] - c2), nu+j-1))+j-1

        if(i!=j){
            ## Swap i and j
            tmp <- swap(a,b,R,i,j)
            a <- tmp$a
            b <- tmp$b
            R <- tmp$R

            C[c(i,j),]    <- as.matrix(C[c(j,i),])
            C[j, (j+1):i] <- as.matrix(0, ncol = i-j, nrow = 1)

            d[c(i,j)] <- d[c(j,i)]
        }
        ## Update Cholesky
        C[j,j] <- sqrt(R[j,j]-sum(C[j,1:(j-1)]^2))
        if(j< (q-1)) C[(j+1):q,j] <- (R[(j+1):q,j]-as.matrix(C[(j+1):q,1:(j-1)])%*%C[j,1:(j-1)])/C[j,j]
        else C[(j+1):q,j] <- (R[(j+1):q,j]-C[(j+1):q,1:(j-1)]%*%C[j,1:(j-1)])/C[j,j]

        ## Get yj
        ajbj <- c01 * c( (a[j]-C[j,1:(j-1)]%*%y[1:(j-1)]),(b[j]-C[j,1:(j-1)]%*%y[1:(j-1)]))/C[j,j]
        y[j] <- exp(lgamma((nu+j)/2)-lgamma((nu+j-1)/2)) / (sqrt((nu+j-1)*pi))*(nu+j-1)/(nu+j+2)*
            ( (1+ajbj[1]^2/2)^(-(nu+j+2)/2) - (1+ajbj[2]^2/2)^(-(nu+j+2)/2))/(pt(ajbj[2],nu+j-1)-pt(ajbj[1],nu+j-1))*
            sqrt( (nu+sumysq)/(nu+j-1))

        sumysq <- sumysq + y[j]^2
    }
    C[q,q] <- sqrt(R[q,q] - sum(C[q,1:(q-1)]^2))
    list(a = a, b = b, R = R, C = C, d = d)
}

##' @title Evaluating the Integrand
##' @param n sample size
##' @param skip skip for sobol
##' @param C Cholesky factor
##' @param q dimension of the problem
##' @param ONE largest number x < 1 such that x != 1
##' @param ZERO smallest number x>0 such that x != 0
##' @param ... all the other parameters same as pnvmix
##' @return mean((f(U)+f(1-U))/2) (scalar)
##' @author Erik Hintz
##' @note Switches between f_ and g_ depending on formt
func <- function(n, skip, a, b, C, q, nu, ONE, ZERO, allinfina)
{

    U <- sobol(n = n, d = (q - 1), randomize = TRUE, skip = skip)
    if(allinfina){
        m <- .Call("evalfbonly_",
                   n    = as.integer(n),
                   q    = as.integer(q),
                   U    = as.double(U),
                   b    = as.double(b),
                   C    = as.double(C),
                   nu   = as.double(nu),
                   ONE  = as.double(ONE),
                   ZERO = as.double(ZERO))
    } else {
        m <- .Call("evalf_",
                   n    = as.integer(n),
                   q    = as.integer(q),
                   U    = as.double(U),
                   a    = as.double(a),
                   b    = as.double(b),
                   C    = as.double(C),
                   nu   = as.double(nu),
                   ONE  = as.double(ONE),
                   ZERO = as.double(ZERO))
    }
    return(m)
}

##' @title Distribution Function of Multivariate Normal Variance Mixtures
##' @param a vector of lower limits
##' @param b vector of upper limits
##' @param R covariance matrix
##' @param nu degrees of freedom
##' @param gam,eps,Nmax,N,n_init tuning parameters. N=number of rep'ns to get sigmahat, algorithm runs unitl gam*sighat < eps or total number evaluation >= Nmax. First loop uses n_init samples in each of the N runs.
##' @param precond logical. If TRUE, preconditioning as described in [genzbretz2002] pp. 955-956 is performed.
##' @author Erik Hintz
pnvmix <- function(a, b, R, nu, gam = 3.3, eps = 0.001, Nmax = 1e8, N = 10, n_init = 2^5, precond = TRUE)
{
    if(!is.matrix(R)) R <- as.matrix(R)
    q <- dim(R)[1] # dimension of the problem
    ## Checks
    if( length(a) != length(b) ) stop("Lenghts of a and b differ")
    if( any(a > b) ) stop("a needs to be smaller than b")
    if( q != length(a) ) stop("Dimension of R does not match dimension of a")

    ## Find infinite limits
    infina  <-  (a == -Inf)
    infinb  <-  (b == Inf)
    infinab <-  infina * infinb

    ## Check if all a_i are -Inf
    allinfina <- (sum(infina)==q)

    ## Remove double infinities
    if( sum(infinab) >0 )
    {
        whichdoubleinf <- which( infinab == TRUE)
        a <- a[ -whichdoubleinf ]
        b <- b[ -whichdoubleinf ]
        R <- R[ -whichdoubleinf, -whichdoubleinf ]
        ## Update dimension
        q <- dim(R)[1]
    }

    ## Deal with the univariate case
    if(q == 1) return( pt(b, df = nu) - pt(a, df = nu) )

    ## Get Cholesky factor (lower triangular)
    C <- t(chol(R))

    ## precondtioning (resorting the limits (cf precond_t); only for q > 2)
    if(precond && q>2) {
        temp <- precond_t(a, b, R, C, q, nu)
        a <- temp$a
        b <- temp$b
        C <- temp$C
    }

    gam <- gam / sqrt(N) # instead of dividing sigma by sqrt(N) each time
    n. <- n_init # initial n
    T. <- rep(NA, N) # vector to safe RQMC estimates

    ONE <- 1-.Machine$double.neg.eps
    ZERO <- .Machine$double.eps
    seed <- .Random.seed # need to reset to the seed later when sobol is being used.

    for(l in 1:N)
        T.[l] <- func(n = n., skip = 0, a = a, b = b, C = C, q = q, nu = nu, ONE = ONE, ZERO = ZERO, allinfina = allinfina)
        ## func returns the average of f(u) and f(1-u)
    N. <- 2 * N * n. # N. will count the total number of f. evaluations

    sig <- sd(T.)
    err <- gam * sig # note that this gam is actually gamma/sqrt(N)

    i. <- 0
    while(err > eps && N. < Nmax)
    {
        .Random.seed <- seed # reset seed to have the same shifts in sobol( ... )
        for(l in 1:N)
            T.[l] <- ( T.[l] + func(n = n., skip = n.,  a = a, b = b, C = C, q = q, nu = nu, ONE = ONE, ZERO = ZERO, allinfina = allinfina) )/2
        ## Note that T.[l] and func(...) both depend on n. evaluations; hence we can just average them
        N. <- N. + 2 * N * n.
        n. <- 2 * n.
        sig <- sd(T.)
        err <- gam * sig # note that this gam is actually gamma/sqrt(N)
        i. <- i. + 1
    }
    T <- mean(T.)
    list(Prob = T, N = N., i = i., ErrEst = err)
}

## ##' @title Evaluate t Copulas
## ##' @param u (n, d) matrix of evaluation points
## ##' @param R covariance matrix
## ##' @param nu degrees of freedom
## ##' @return list of four: estimated prob, total number of function evaluations needed, number of iterations
## ##'         in the while loop, estimated error
## ##' @author Erik Hintz and Marius Hofert
## t_cop_prob <- function(u, R, nu)
## {
##     if(!is.matrix(u)) u <- rbind(u)
##     t_prob(rep(-Inf, ncol(u)), b = qt(u, nu), R = R, nu = nu)
## }
