## qnvmix1()

qnvmix1 <- function(u, shift = 0, scale = 1, mix, N = 1e5, ... )
{
    if(u < 0.5){
        u <- 1 - u
        lower <- TRUE
    } else {
        lower <- FALSE
    }

    ## Define the quantile function of the mixing variable:
    ## (Copy and paste from pnvmix())
    W <- if(is.character(mix)) { # 'mix' is a character vector specifying supported mixture distributions (utilizing '...')
             mix <- match.arg(mix, choices = c("constant", "inverse.gamma"))
             switch(mix,
                    "constant" = {
                        const <- TRUE
                        function(u) 1
                    },
                    "inverse.gamma" = {
                        if(hasArg(df)) df <- list(...)$df else stop("'mix = \"inverse.gamma\"' requires 'df' to be provided.")
                        ## Still allow df = Inf (normal distribution)
                        stopifnot(is.numeric(df), length(df) == 1, df > 0)
                        if(is.finite(df)) {
                            inv.gam <- TRUE
                            df2 <- df / 2
                            meansqrtmix <- sqrt(df) * gamma(df2) / ( sqrt(2) * gamma( (df+1) / 2 ) ) # mean of sqrt(W) in this case, will be used for preconditioning
                            function(u)
                                1 / qgamma(u, shape = df2, rate = df2)
                        } else {
                            const <- TRUE
                            meansqrtmix <- 1 # mean of sqrt(W) in this case, will be used for preconditioning
                            function(u) 1
                        }
                    },
                    stop("Currently unsupported 'mix'"))
         } else if(is.list(mix)) { # 'mix' is a list of the form (<character string>, <parameters>)
             stopifnot(length(mix) >= 1, is.character(distr <- mix[[1]]))
             qmix <- paste0("q", distr)
             if(!existsFunction(qmix))
                 stop("No function named '", qmix, "'.")
             function(u) do.call(qmix, c(u, mix[-1]))
         } else if(is.function(mix)) { # 'mix' is interpreted as the quantile function F_W^- of the mixture distribution F_W of W
             function(u) mix(u, ...)
         }

    ## Get realizations of sqrt(W)
    ## Mutliply by sqrt(scale) for standardization
    mixings <- sqrt(W(qrng::sobol(n = N, d = 1, randomize = TRUE)) * scale)

    ## Evaluates the cdf of the mixture distribution using QMC. Can use the same realizations of mixings in each run -> faster (like common random numbers)
    cdf <- function(b)
        mean( pnorm( (b - shift) / mixings))

    q <- uniroot( function(x){cdf(x)-u}, lower = 0, upper = 10^8, extendInt = "upX")$root # TODO: return whole object
    if(lower) -q else q
}
