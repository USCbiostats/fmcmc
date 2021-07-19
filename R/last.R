#' Information about the last `MCMC` call
#' 
#' This environment holds a copy of the last call to [MCMC], including the start
#' and end time (to compute total elapsed time) of the call. Since the resulting
#' object of `MCMC` is an object of class [coda::mcmc], this is a way to capture
#' more information in case the user needs it.
#' 
#' 
#' @name fmcmc-info
#' @return The `MCMC_INFO` object is an environment of class, 
#' `c("fmcmc_info", "environment")` that has the following structure:
#' 
#' - `time_start`, `time_end` Objects of class [proc.time]. Mark the start and
#'   end of the `MCMC` call.
#' 
#' - `data.` A list of environments of length `get_nchains()`. Each environment
#'   will hold information about the particular chain. By default, each environment
#'   holds the elements `logpost` (named numeric vector) and `draws` (numeric matrix).
#'   
#'   The `draws` matrix contains the draws from the proposal kernel function. Both 
#'   `logpost` and `draws` have indices that match those of the chain. (see details).
#'   
#'   `data.` can also be accessed by the user to store information if needed.
#'   
#' - `ptr` An environment. This is used as a pointer that is defined at the beginning
#'   of the MCMC process. The environment will be pointing to the current chain, thus,
#'   if `MCMC` is running chain 2 of 4, `ptr = data.[[2]]`.
#'   
#' - `i` Integer. Index of the current chain, so if `MCMC` is running chain 3 of 4,
#'   then `i = 3` (and `ptr = data.[[3]]`, get it?).
#'   
#' - `nchains` Integer. The number of chains specified in `MCMC`.
#' 
#' - `...` further arguments passed to `MCMC`, e.g., `initial`, `fun`, `nsteps`,
#'   `kernel`, `thin`, etc.
#'   
#' It also contains the following **helper functions**:
#' 
#' - `c_(x, val)` Combine elements. It will access the current value of `x` in
#'   `ptr`, and will combine it with `val` using [c()]. 
#' 
#' - `rbind_(x, val)` Row-combine elements. It will access the current value of `x` in
#'   `ptr`, and will combine it with `val` using [rbind()].
#' 
#' - `cbind_(x, val)` Column-combine elements. It will access the current value of `x` in
#'   `ptr`, and will combine it with `val` using [cbind()].
#' 
#' 
#' `get_*` returns the corresponding variable passed to the last call
#' of [MCMC].
#' @export
MCMC_INFO <- structure(
  list2env(list(
    time_start = NULL,
    time_end   = NULL,
    data.      = list(),
    ptr        = NULL,
    i          = NA_integer_,
    nchains    = 0L,
    loop_envir = NULL   
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

#' Combine
#' @noRd
MCMC_INFO$cbind_ <- function(x, val) {
  
  assign(x = x, value = cbind(MCMC_INFO$ptr[[x]], val), envir = MCMC_INFO$ptr)
  
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
#' # Getting the logpost -------------------------------------------------------
#' set.seed(23133)
#' x <- rnorm(200)
#' y <- x*2 + rnorm(200)
#' f <- function(p) {
#'   sum(dnorm(y - x*p, log = TRUE))
#' }
#' 
#' ans <- MCMC(fun = f, initial = c(0), nsteps=2000)
#' plot(get_logpost(), type = "l") # Plotting the logpost from the last run
#' 
#' 
#' # Printing information every 500 step ---------------------------------------
#' # for this we use ith_step()
#' 
#' f <- function(p) {
#' 
#'   # Capturing info from within the loop
#'   i      <- ith_step("i")
#'   nsteps <- ith_step("nsteps")
#'   
#'   if (!(i %% 500)) {
#'   
#'     cat(
#'       "////////////////////////////////////////////////////\n",
#'       "Step ", i, " of ", nsteps,". Values in the loop:\n",
#'       "theta0: ", ith_step("theta0"), "\n",
#'       "theta1: ", ith_step()$theta1, "\n",
#'       sep = ""
#'     )

#'   }
#'     
#' 
#'   sum(dnorm(y - x*p, log = TRUE))
#' }
#' 
#' MCMC(fun = f, initial = c(0), nsteps=2000, progress = FALSE, seed = 22)
#' # ////////////////////////////////////////////////////
#' # Step 500 of 2000. Values in the loop:
#' # theta0: 2.025379
#' # theta1: 1.04524
#' # ////////////////////////////////////////////////////
#' # Step 1000 of 2000. Values in the loop:
#' # theta0: 2.145967
#' # theta1: 0.2054037
#' # ////////////////////////////////////////////////////
#' # Step 1500 of 2000. Values in the loop:
#' # theta0: 2.211691
#' # theta1: 2.515361
#' # ////////////////////////////////////////////////////
#' # Step 2000 of 2000. Values in the loop:
#' # theta0: 1.998789
#' # theta1: 1.33034
#' 
#' 
#' # Printing information if the current logpost is greater than max -----------
#' f <- function(p) {
#' 
#'   i            <- ith_step("i")
#'   logpost_prev <- max(ith_step("logpost")[1:(i-1)])
#'   logpost_curr <- sum(dnorm(y - x*p, log = TRUE))
#'   
#'   # Only worthwhile after the first step
#'   if ((i > 1L) && logpost_prev < logpost_curr)
#'     cat("At a higher point!:", logpost_curr, ", step:", i,"\n")
#'     
#'   return(logpost_curr)
#' 
#' }
#' MCMC(fun = f, initial = c(0), nsteps=1000, progress = FALSE, seed = 22)
#' # At a higher point!: -357.3584 , step: 2 
#' # At a higher point!: -272.6816 , step: 6 
#' # At a higher point!: -270.9969 , step: 7 
#' # At a higher point!: -269.8128 , step: 24 
#' # At a higher point!: -269.7435 , step: 46 
#' # At a higher point!: -269.7422 , step: 543 
#' # At a higher point!: -269.7421 , step: 788 
get_logpost <- function() {
  
  if (get_nchains() == 1L)
    return(get_("logpost")[[1L]])
  else
    return(get_("logpost"))
  
}

#' @export
#' @details The function `get_draws()` retrieves the proposed states from the
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

#' @export
#' @section Advanced usage:
#' The function [ith_step(x)] is a convenience function that provides
#' access to the environment within which the main loop of the MCMC call is
#' being evaluated. This is a wrapper of `MCMC_INFO$loop_envir` that will
#' either return the value `x` of the entire environment if `x` is missing. If
#' `ith_step()` is called outside of the `MCMC` call, then it will return with
#' an error.
#' 
#' For example, if you wanted to print information if the current value
#' of logpost is greater than the previous value of logpost, you could define
#' the objective function as follows:
#' 
#' ```
#' f <- function(p) {
#' 
#'   i            <- ith_step("i")
#'   logpost_prev <- ith_step("logpost")[i - 1L]
#'   logpost_curr <- sum(dnorm(y - x*p, log = TRUE))
#'   
#'   if (logpost_prev < logpost_curr)
#'     cat("At a higher point!\n")
#'     
#'   return(logpost_curr)
#' 
#' }
#' ```
#' 
#' More examples below.
#' 
#' @rdname fmcmc-info
ith_step <- function(x) {
  
  if (is.null(MCMC_INFO$loop_envir))
    stop("-ith_step()- should only be used within an -MCMC()- call.", call. = FALSE)
  
  if (missing(x))
    return(MCMC_INFO$loop_envir)
  else
    return(get(x, envir = MCMC_INFO$loop_envir))
  
}