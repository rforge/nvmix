### get.set.parameters() #######################################################

##' @title  Retrieve algorithm specific default parameters and overwrite them
##' @return list with default values for all functions in the 'nvmix' package
##' @note newton.logdens.abstol not very small as it is only needed to get the
##' the next iteration in the Newton procedure *unless* it is being called 
##' from dnvmixcop() in which case this tolerance is changed there. 
##' @author Erik Hintz and Marius Hofert

get.set.parameters <- function(control = list()){
  ## Set up default controls:
  ctrl <- list(
    ## For pnvmix(): 
    mean.sqrt.mix = NULL, 
    precond = TRUE, 
    pnvmix.abstol = 1e-3, 
    ## For dnvmix():
    dnvmix.abstol = 1e-3, 
    dnvmix.reltol = 0.025, # If !NA, 'reltol' is used instead of 'abstol'
    dnvmix.max.iter.rqmc.pilot = 4,
    dnvmix.tol.int.lower = 1e-30,
    dnvmix.tol.bisec.w = 0.1,
    dnvmix.tol.stratlength = 1e-20,
    dnvmix.max.iter.bisec.w = 55,
    ## For pgammamix:
    pgammamix.reltol = NA,
    pgammamix.abstol = 1e-3,
    ## For qnvmix():
    max.iter.newton = 40, 
    newton.conv.abstol = 1e-4,
    newton.df.abstol = 1e-4,
    newton.logdens.abstol = 1e-2, 
    ## For fitnvmix():
    weights.abstol = 1e-1, 
    ECMEstep.do.nu = TRUE,
    laststep.do.nu = TRUE,
    ECME.maxiter = 20,
    ECME.conv.tol = c(rep(1e-1, 2), 1e-2), # [1] => 'loc'; [2] => 'scale'; [3] => 'nu'
    ## For all (randomized) algorithms:
    method = "sobol", 
    increment = "doubling", # "doubling" or "num.init" 
    max.iter.rqmc = NA, # defined below, depending on 'increment'
    CI.factor = 3.3,
    fun.eval = c(2^7, 1e12), 
    B = 15)
  if(length(control) > 0){
    ## If input provided, grab input controls and overwrite:
    names.control <- names(ctrl)
    ctrl[(names.provided <- names(control))] <- control
    ## Did the user provide something that is not used?
    if (length(unmatched <- names.provided[!names.provided %in% names.control])) 
      warning("unknown names in control: ", paste(unmatched, collapse = ", "))
    ## Check if 'method' and 'increment' were provided correctly
    ctrl$method     <- match.arg(ctrl$method, 
                                 choices = c("sobol", "ghalton", "PRNG"))
    ctrl$increment  <- match.arg(ctrl$increment, 
                                 choices = c("doubling", "num.init"))
    ## Now some more checkings: ('max.iter.rqmc' checked at the end)
    stopifnot(is.logical(ctrl$precond),
              ctrl$pnvmix.abstol >= 0,
              ctrl$dnvmix.abstol >= 0,
              ctrl$dnvmix.reltol >= 0,
              ctrl$dnvmix.max.iter.rqmc.pilot >= 1,
              ctrl$dnvmix.tol.int.lower > 0,
              ctrl$dnvmix.tol.bisec.w >0,
              ctrl$dnvmix.tol.stratlength > 0,
              ctrl$dnvmix.max.iter.bisec.w > 0, 
              ctrl$max.iter.newton >= 0,
              ctrl$newton.conv.abstol >= 0,
              ctrl$newton.df.abstol >= 0,
              ctrl$newton.logdens.abstol >= 0,
              ctrl$weights.abstol >= 0,
              is.logical(ctrl$ECMEstep.do.nu),
              is.logical(ctrl$laststep.do.nu),
              ctrl$ECME.maxiter >= 0,
              length(ctrl$ECME.conv.tol) == 3, ctrl$ECME.conv.tol >= 0, 
              ctrl$CI.factor >= 0,
              length(ctrl$fun.eval) == 2, ctrl$fun.eval >= 0,
              ctrl$B > 1) # If B = 1 error estimates are NA => need B > 1
  }
  ## Define 'max.iter.rqmc': If it was not provided (=> NA), set defaults
  if(is.na(ctrl$max.iter.rqmc)){
    ctrl$max.iter.rqmc <- if(ctrl$increment == "doubling") 15 else 200
  } else {
    ## If it was provided (=> not NA), check if it's reasonable 
    stopifnot(ctrl$max.iter.rqmc > 1)
  }
  ## Return
  ctrl
}

