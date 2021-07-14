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
LAST_RUN <- structure(
  list2env(list(
    data.   = list(),
    ptr     = NULL,
    i       = NA_integer_,
    nchains = 0L
  ), envir = new.env()),
  class = c("fmcmc_run", "environment"))

#' @noRd Clears the LAST_RUN environment and sets the number of chains
#' @param nchains an integer
#' @param env an optional environment (this is useful in the context of parallel chains.)
LAST_RUN$clear <- function(nchains, env) {
  
  lapply(LAST_RUN$data., function(d) rm(list = ls(all.names = TRUE, envir = d), envir = d))
  
  if (missing(env))
    LAST_RUN$data. <- replicate(
      nchains, list2env({
        list(logpost = numeric())
      }
      ), simplify = FALSE)
  else {
    
    if (nchains != 1L)
      stop("When clearing with an environment, nchains should be 1.")
    
    LAST_RUN$data. <- list(env)
    
  }
    
  LAST_RUN$nchains <- nchains
  LAST_RUN$info    <- new.env()
  invisible()
  
}

#' @noRd This function sets the current pointer (useful when running in
#' serial fashion.)
LAST_RUN$set_ptr <- function(i) {
  
  if (i > LAST_RUN$nchains)
    stop("fmcmc_run pointer out of range.", call. = FALSE)
  
  LAST_RUN$ptr <- LAST_RUN$data.[[i]]
  LAST_RUN$i   <- i
  
  invisible()
  
}

#' @noRd Combine
LAST_RUN$c_ <- function(x, val) {
  
  assign(x = x, value = c(LAST_RUN$ptr[[x]], val), envir = LAST_RUN$ptr)
  
}

#' @export
print.fmcmc_run <- function(x, ...) {
  
  if (x$nchains == 0L)
    cat("-MCMC- has not been called yet. Nothing to show.\n")
  else {
    cat("Last call to MCMC holds the following elements:\n")
    lapply(x$data., function(d) print(utils::ls.str(d)))
  }
  
  invisible(x)
}

MCMC_init <- function(...) {
  
  # Getting the caller environment
  # LAST_RUN$time_start   <- proc.time()
  LAST_RUN$time_start <- proc.time()
  env <- parent.frame()
  
  # Initializing the variables
  for (n in names(env))
    if (n != "...") 
      assign(n, get(n, envir = env), envir = LAST_RUN$info)
    
  # Assigning dots
  for (i in seq_len(...length())) {
    assign(...names()[i], ...elt(i), envir = LAST_RUN$info)
  }
    
  invisible(NULL)
  
}

MCMC_finalize <- function() {
  LAST_RUN$time_end <- proc.time()
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
  
  if (LAST_RUN$nchains == 0L)
    stop("-MCMC- has not been called yet.", call. = FALSE)
  
  # Checking if it exists outside
  if (exists(x, envir = LAST_RUN)) {
    return(LAST_RUN[[x]])
  }
  
  if (exists(x, envir = LAST_RUN$info)) {
    return(LAST_RUN$info[[x]])
  }
  
  # Otherwise, it should be part of the run
  lapply(LAST_RUN$data., function(env) {
    
    if (!exists(x, envir = env))
      stop(
        "The object -", x, "- was not found in the last MCMC call.",
        call. = FALSE
      )
    
    env[[x]]
    
  })
  
}

#' @export
#' @rdname last-mcmc
#' @details The function `get_logpost` returns the `logposterior` value at each
#' iteration. The values correspond to a named numeric vector.
get_logpost <- function() {
  
  last_("logpost")
  
}

