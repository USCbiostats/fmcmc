#' Deprecated methods in fmcmc
#' 
#' These functions will no longer be included starting version 0.6-0. Instead,
#' use the functions in [mcmc-output].
#' 
#' @export
#' @name fmcmc-deprecated
#' @return The function `last_elapsed` returns the elapsed time of the last call
#' to [MCMC]. In particular, the `MCMC` function records the running time of R
#' at the beginning and end of the function using [proc.time()]. So this function
#' returns the difference between the two (`time_end - time_start`).
last_elapsed <- function() {
  last_("time_end") - last_("time_start")
}


#' @export
#' @rdname fmcmc-deprecated
LAST_MCMC <- structure(list(), class = "fmcmc_last_mcmc")

#' @export
`[[.fmcmc_last_mcmc` <- function(i, j, ..., exact=TRUE) {
  .Deprecated("MCMC_OUTPUT", old = "LAST_MCMC")
  MCMC_OUTPUT[[i]]
}

#' @export
`$.fmcmc_last_mcmc` <- function(x, name) {
  .Deprecated("MCMC_OUTPUT", old = "LAST_MCMC")
  MCMC_OUTPUT[[name]]
}

#' @export
#' @rdname fmcmc-deprecated
last_nsteps <- function() last_("nsteps")

#' @export
#' @rdname fmcmc-deprecated
last_nchains <- function() last_("nchains")

#' @export
#' @rdname fmcmc-deprecated
last_kernel <- function() last_("kernel")

#' @export
#' @rdname fmcmc-deprecated
last_conv_checker <- function() last_("conv_checker")

#' @export
#' @rdname fmcmc-deprecated
#' @param x Name of the object to retrieve.
last_ <- function(x) {
  .Deprecated("get_", "The -last_*- methods will be deprecated in the next version of -fmcmc-. Use -get_*- instead.")
  get_(x)
}
