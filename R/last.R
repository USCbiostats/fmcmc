#' Information about the last `MCMC` call
#' 
#' This environment holds a copy of the last call to [MCMC], including the start
#' and end time (to compute total elapsed time) of the call. Since the resulting
#' object of `MCMC` is an object of class [coda::mcmc], this is a way to capture
#' more information in case the user needs it.
#' 
#' @name fmcmc-info
#' @return `get_*` returns the corresponding variable passed to the last call
#' of [MCMC].
#' @export
MCMC_INFO <- structure(
  list2env(list(
    data.   = list(),
    ptr     = NULL,
    i       = NA_integer_,
    nchains = 0L
  ), envir = new.env()),
  class = c("fmcmc_info", "environment"))

#' Clears the MCMC_INFO environment and sets the number of chains
#' @param nchains an integer
#' @param env an optional environment (this is useful in the context of parallel chains.)
#' @noRd 
MCMC_INFO$clear <- function(nchains, env) {
  
  lapply(MCMC_INFO$data., function(d) rm(list = ls(all.names = TRUE, envir = d), envir = d))
  
  if (missing(env))
    MCMC_INFO$data. <- replicate(
      nchains, list2env({
        list(logpost = numeric(), draws = NULL)
      }
      ), simplify = FALSE)
  else {
    
    if (nchains != 1L)
      stop("When clearing with an environment, nchains should be 1.")
    
    MCMC_INFO$data. <- list(env)
    
  }
    
  MCMC_INFO$nchains <- nchains
  MCMC_INFO$info    <- new.env()
  invisible()
  
}

#' This function sets the current pointer (useful when running in
#' serial fashion.)
#' @noRd 
MCMC_INFO$set_ptr <- function(i) {
  
  if (i > MCMC_INFO$nchains)
    stop("fmcmc_info pointer out of range.", call. = FALSE)
  
  MCMC_INFO$ptr <- MCMC_INFO$data.[[i]]
  MCMC_INFO$i   <- i
  
  invisible()
  
}

#' Combine
#' @noRd
MCMC_INFO$c_ <- function(x, val) {
  
  assign(x = x, value = c(MCMC_INFO$ptr[[x]], val), envir = MCMC_INFO$ptr)
  
}

#' Combine
#' @noRd
MCMC_INFO$rbind_ <- function(x, val) {
  
  assign(x = x, value = rbind(MCMC_INFO$ptr[[x]], val), envir = MCMC_INFO$ptr)
  
}

#' @export
print.fmcmc_info <- function(x, ...) {
  
  if (x$nchains == 0L)
    cat("-MCMC- has not been called yet. Nothing to show.\n")
  else {
    cat("Last call to MCMC holds the following elements:\n")
    lapply(x$data., function(d) print(utils::ls.str(d)))
  }
  
  invisible(x)
}

#' @export
print.fmcmc_last_mcmc <- function(x, ...) {
  
  .Deprecated("MCMC_INFO")
  print(MCMC_INFO)
  invisible(x)
  
}

MCMC_init <- function(...) {
  
  # Getting the caller environment
  # MCMC_INFO$time_start   <- proc.time()
  MCMC_INFO$time_start <- proc.time()
  env <- parent.frame()
  
  # Initializing the variables
  for (n in names(env))
    if (n != "...") 
      assign(n, get(n, envir = env), envir = MCMC_INFO$info)
    
  # Assigning dots
  for (i in seq_len(...length())) {
    assign(...names()[i], ...elt(i), envir = MCMC_INFO$info)
  }
    
  invisible(NULL)
  
}

MCMC_finalize <- function() {
  MCMC_INFO$time_end <- proc.time()
}


#' @export
#' @rdname fmcmc-info
#' @param x Character scalar. Name of an argument to retrieve. If `x` was not
#' passed to the last call, the function returns with an error.
get_ <- function(x) {
  
  if (MCMC_INFO$nchains == 0L)
    stop("-MCMC- has not been called yet.", call. = FALSE)
  
  # Checking if it exists outside
  if (exists(x, envir = MCMC_INFO, inherits = FALSE)) {
    return(MCMC_INFO[[x]])
  }
  
  if (exists(x, envir = MCMC_INFO$info, inherits = FALSE)) {
    return(MCMC_INFO$info[[x]])
  }
  
  # Otherwise, it should be part of the run
  lapply(MCMC_INFO$data., function(env) {
    
    if (!exists(x, envir = env, inherits = FALSE))
      stop(
        "The object -", x, "- was not found in the last MCMC call.",
        call. = FALSE
      )
    
    env[[x]]
    
  })
  
}

#' @export
#' @rdname fmcmc-info
#' @details The function `get_logpost` returns the `logposterior` value at each
#' iteration. The values correspond to a named numeric vector. If `nchains > 1`
#' then it will return a list of length `nchains` with the corresponding logpost
#' values for each chain.
#' @examples 
#' set.seed(23133)
#' x <- rnorm(200)
#' y <- x*2 + rnorm(200)
#' f <- function(p) {
#'   sum(dnorm(y - x*p, log = TRUE))
#' }
#' 
#' ans <- MCMC(fun = f, initial = c(0), nsteps=2000)
#' plot(get_logpost(), type = "l") # Plotting the logpost from the last run
get_logpost <- function() {
  
  if (get_nchains() == 1L)
    return(get_("logpost")[[1L]])
  else
    return(get_("logpost"))
  
}

#' @export
#' @details The function `get_draws()` retrieves the proposed samples from the
#' kernel function.
#' @rdname fmcmc-info
get_draws <- function() {
  
  if (get_nchains() == 1L)
    return(get_("draws")[[1L]])
  else
    return(get_("draws"))
  
}


#' @export
#' @rdname fmcmc-info
get_elapsed <- function() {
  get_("time_end") - get_("time_start")
}

#' @export
#' @rdname fmcmc-info
get_nsteps <- function() get_("nsteps")

#' @export
#' @rdname fmcmc-info
get_nchains <- function() get_("nchains")

#' @export
#' @rdname fmcmc-info
get_kernel <- function() get_("kernel")

#' @export
#' @rdname fmcmc-info
get_conv_checker <- function() get_("conv_checker")

#' @export
#' @rdname fmcmc-info
get_seed <- function() get_("seed")

#' @export
#' @rdname fmcmc-info
get_burnin <- function() get_("seed")

#' @export
#' @rdname fmcmc-info
get_thin <- function() get_("seed")


#' @export
#' @rdname fmcmc-deprecated
LAST_MCMC <- structure(list(), class = "fmcmc_last_mcmc")

#' @export
`[[.fmcmc_last_mcmc` <- function(i, j, ..., exact=TRUE) {
  .Deprecated("MCMC_INFO", old = "LAST_MCMC")
  MCMC_INFO[[i]]
}

#' @export
`$.fmcmc_last_mcmc` <- function(x, name) {
  .Deprecated("MCMC_INFO", old = "LAST_MCMC")
  MCMC_INFO[[name]]
}

#' Deprecated methods in fmcmc
#' 
#' These functions will no longer be included starting version 0.6-0. Instead,
#' use the functions in [fmcmc-info].
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
