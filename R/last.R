#' Information about the last `MCMC` call
#' 
#' This environment holds a copy of the last call to [MCMC], including the start
#' and end time (to compute total elapsed time) of the call. Since the resulting
#' object of `MCMC` is an object of class [coda::mcmc], this is a way to capture
#' more information in case the user needs it.
#' 
#' @name last-mcmc
#' @return `last_*` returns the corresponding variable.
#' @export
LAST_MCMC <- structure(new.env(), class = c("fmcmc_last", "environment"))

#' @export
print.fmcmc_last <- function(x, ...) {
  
  if (length(x) == 0L)
    cat("-MCMC- has not been called yet. Nothing to show.\n")
  else {
    cat("Last call to MCMC holds the following elements:\n")
    print(utils::ls.str(x))
  }
  
  invisible(x)
}

LAST_MCMC_add <- function(...) {
  
  if (any(is.na(...names())))
    stop(
      "All elements passed to -LAST_MCMC_add- should have a name.",
      call. = FALSE
    )
  
  for (i in 1:...length()) {
    
    assign(
      x = ...names()[i],
      value = ...elt(i),
      envir = LAST_MCMC
      )
    
    LAST_MCMC$.internal_elements <- c(
      LAST_MCMC$.internal_elements,
      ...names()[i]
    )
    
  }
  
  invisible()
}

MCMC_init <- function(...) {
  
  # Getting the caller environment
  # LAST_MCMC$time_start   <- proc.time()
  LAST_MCMC_add(time_start = proc.time())
  env <- parent.frame()
  
  # Initializing the variables
  for (n in names(env))
    if (n != "...") 
      
      assign(n, get(n, envir = env), envir = LAST_MCMC)
    
  # Assigning dots
  dots <- list(...)
  if (length(dots))
    for (n in names(dots)) {
      assign(n, dots[[n]], envir = LAST_MCMC)
    }
    
  invisible(NULL)
  
}

MCMC_finalize <- function() {
  LAST_MCMC_add(time_end = proc.time())
  # LAST_MCMC$time_end <- proc.time()
}



#' @export
#' @rdname last-mcmc
#' @return The function `last_elapsed` returns the elapsed time of the last call
#' to [MCMC]. In particular, the `MCMC` function records the running time of R
#' at the beginning and end of the function using [proc.time()]. So this function
#' returns the difference between the two (`time_end - time_start`).
last_elapsed <- function() {
  last_("time_end") - last_("time_start")
}

#' @export
#' @rdname last-mcmc
last_nsteps <- function() last_("nsteps")

#' @export
#' @rdname last-mcmc
last_nchains <- function() last_("nchains")

#' @export
#' @rdname last-mcmc
last_kernel <- function() last_("kernel")

#' @export
#' @rdname last-mcmc
last_conv_checker <- function() last_("conv_checker")

#' @export
#' @rdname last-mcmc
#' @param x Character scalar. Name of an argument to retrieve. If `x` was not
#' passed to the last call, the function returns with an error.
last_ <- function(x) {
  
  if (length(LAST_MCMC) == 0L)
    stop("-MCMC- has not been called yet.", call. = FALSE)
  
  if (!exists(x, envir = LAST_MCMC))
    stop(
      "The object -", x, "- was not found in the last MCMC call.",
      call. = FALSE
      )
  
  LAST_MCMC[[x]]
}

#' @export
#' @rdname last-mcmc
#' @details The function `get_logpost` returns the `logposterior` value at each
#' iteration. The values correspond to a named numeric vector.
get_logpost <- function() {
  
  last_("logpost")
  
}

#' @noRd
#' Sets the the value for LAST_MCMC prob
set_last_mcmc_ <- function(x, value, nchain) {
  
  # Making space
  if (!exists(x, envir = LAST_MCMC))
    assign(x, rep(NULL, last_nchains()), envir = LAST_MCMC)
  
  if (last_nchains() > 1L) {
    
    LAST_MCMC[[x]][[nchain]] <- c(LAST_MCMC[[x]][[nchain]], value)
    
  } else {
    
    LAST_MCMC[[x]] <- c(LAST_MCMC[[x]], value)
    
  }
  
  invisible()
  
}