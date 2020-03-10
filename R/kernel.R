#' Checks the dimensions of the parameter according to k.
#' @noRd
check_dimensions <- function(x, k) {
  
  # Getting the name of the variable
  name <- deparse(match.call()$x)
  
  if (length(x) > 1L &&  length(x) != k)
    stop("Incorrect length of -", name, "-.", call. = FALSE)
  
  if (length(x) == 1L && k > 1L)
    return(rep(x, k))
  
  return(x)
  
  
}

#' Parameters' update sequence
#' @param k Integer. Number of parameters
#' @param nsteps Integer. Number of steps.
#' @param fixed Logical scalar or vector of length `k`. Indicates which parameters
#' will be treated as fixed or not. Single values are recycled.
#' @param scheme Scheme in which the proposals are made (see details).
#' 
#' @details 
#' The parameter `scheme` present on the currently available kernels sets the way
#' in which proposals are made. By default, `scheme = "joint"`, proposals are done
#' jointly, this is, at each step of the chain we are proposing new states for
#' each parameter of the model. When `scheme = "ordered"`, a sequential update schema
#' is followed, in which, at each step of the chain, proposals are made one
#' variable at a time, If `scheme = "random"`, proposals are also made one
#' variable at a time but in a random scheme.
#' 
#' Finally, users can specify their own sequence of proposals for the variables
#' by passing a numeric vector to `scheme`, for example, if the user wants to make
#' sequential proposals following the scheme 2, 1, 3, then scheme must be set to
#' be `scheme = c(2, 1, 3)`.
#' 
#' @return A logical vector of size `nsteps` x `k`.
#' @export
plan_update_sequence <- function(k, nsteps, fixed, scheme) {
  
  # Setting the scheme in which the variables will be updated
  if (length(scheme) > 1L && is.numeric(scheme)) {
    
    # Is it the right length?
    if (length(scheme) != sum(!fixed))
      stop(
        "When setting the update scheme, it should have the same length ",
        "as the number of variables that will not be fixed. ",
        "Right now length(scheme) = ", length(scheme), " while sum(!fixed) = ",
        sum(!fixed),".", call. = FALSE
        )
    
    # Is the full sequence included?
    test <- which(!(which(!fixed) %in% scheme))
    if (length(test))
      stop(
        "One or more variables was not included in the ordering sequence. ",
        "The full list follows: ", paste(test, collapse=", "), ". ",
        "Only variables that are not fixed can be included in this list.",
        call. = FALSE
        )
    
    update_sequence <- matrix(FALSE, nrow = nsteps, ncol = k)
    suppressWarnings(update_sequence[cbind(1:nsteps, scheme)] <- TRUE)
    
    
  } else if (scheme == "joint") {
    
    update_sequence <- matrix(TRUE, nrow = nsteps, ncol = k)
    
    for (j in which(fixed))
      update_sequence[, j] <- FALSE
  
  } else if (scheme == "ordered") {
    
    update_sequence <- matrix(FALSE, nrow = nsteps, ncol = k)
    suppressWarnings(update_sequence[cbind(1:nsteps, which(!fixed))] <- TRUE)
    
  } else if (scheme == "random") {
    
    update_sequence <- matrix(FALSE, nrow = nsteps, ncol = k)
    update_sequence[cbind(
      1:nsteps,
      sample(which(!fixed), nsteps, TRUE))
      ] <- TRUE
    update_sequence <- update_sequence
    
  } else {
    
    stop(
      "-scheme- update must be either an integer sequence, 'joint', 'ordered', ",
      "or 'random'.",
      call. = FALSE
      )
    
  }
  
  if (sum(update_sequence[1,]) == 0L)
    stop("The number of parameters to update, i.e. not fixed, cannot be ",
         "zero. Check the value -fixed- in the kernel initialization.", 
         call. = FALSE)
    
  
  return(update_sequence)
  
}

#' Transition Kernels for MCMC
#' 
#' The function `kernel_new` is a helper function that allows creating
#' `fmcmc_kernel` objects which are passed to the [MCMC()] function.
#' 
#' @param proposal,logratio Functions. Both receive a single argument, an environment.
#' This functions are called later within [MCMC] (see details).
#' @param kernel_env Environment. This will be used as the main container of the
#' kernel's components. It is returned as an object of class `c("environment", "fmcmc_kernel")`.
#' @param ... In the case of `kernel_new`, further arguments to be stored with
#' the kernel.
#' 
#' @details
#' 
#' The objects `fmcmc_kernels` are environments that in general contain the 
#' following objects:
#' 
#' -  `proposal`: The function used to propose changes in the chain based
#'    on the current state. The function must return a vector of length equal
#'    to the number of parameters in the model.
#'    
#' -  `logratio`: This function is called after a new state has been proposed,
#'    and is used to compute the log of the Hastings ratio.
#'    
#'    In the case that the `logratio` function is not specified, then it is assumed
#'    that the transition kernel is symmetric, this is, log-ratio is then implemented
#'    as `function(env) {env$f1 - env$f0}`
#'    
#' -  `...`: Further objects that are used within those functions.
#' 
#' Both functions, `proposal` and `logratio`, receive a single argument, an
#' environment, which is passed by the [MCMC()] function during each step using
#' the function [environment()].
#' 
#' The passed environment is actually the environment in which the `MCMC`
#' function is running, in particular, this environment contains the following
#' objects:
#' 
#' \tabular{lcl}{
#' \strong{Object} \tab \tab \strong{Description} \cr
#' `i` \tab \tab Integer. The current iteration. \cr
#' `theta1` \tab \tab Numeric vector. The last proposed state. \cr
#' `theta0` \tab \tab Numeric vector. The current state \cr
#' `f`\tab \tab The log-unnormalized posterior function (a wrapper of `fun` passed 
#' to [MCMC]). \cr
#' `f1` \tab \tab The last value of `f(theta1)` \cr
#' `f0` \tab \tab The last value of `f(theta0)` \cr
#' `kernel` \tab \tab The actual `fmcmc_kernel` object. \cr
#' `ans` \tab \tab The matrix of samples defined up to `i - 1`.
#' }
#' 
#' These are the core component of the `MCMC` function. The following block
#' of code is how this is actually implemented in the package:
#' 
#' ```
#' for (i in 1L:nsteps) {
#'   # Step 1. Propose
#'   theta1[] <- kernel$proposal(environment())
#'   f1       <- f(theta1)
#'   
#'   # Checking f(theta1) (it must be a number, can be Inf)
#'   if (is.nan(f1) | is.na(f1) | is.null(f1)) 
#'     stop(
#'       "fun(par) is undefined (", f1, ")",
#'       "Check either -fun- or the -lb- and -ub- parameters.",
#'       call. = FALSE
#'     )
#'   
#'   # Step 2. Hastings ratio
#'   if (R[i] < kernel$logratio(environment())) {
#'     theta0 <- theta1
#'     f0     <- f1
#'   }
#'   
#'   # Step 3. Saving the state
#'   ans[i,] <- theta0
#'   
#' }
#' ```
#' 
#' For an extensive example on how to create new kernel objects see the vignette
#' `vignette("user-defined-kernels", "fmcmc")`.
#' 
#' @section Behavior:
#' 
#' In some cases, calls to the `proposal()` and `logratio()` functions in
#' `fmcmc_kernels` can trigger changes or updates of variables stored within them.
#' A concrete example is with adaptive kernels.
#' 
#' Adaptive Metropolis and Robust Adaptive Metropolis implemented in the functions
#' [kernel_adapt()] and [kernel_ram()] both update a covariance matrix used
#' during the proposal stage, and furthermore, have a `warmup` stage that sets
#' the point at which both will start adapting. Because of this, both kernels
#' have internal counters of the absolute step count which allows them activating,
#' scaling, etc. the proposals correctly.
#' 
#' 1. When running multiple chains, `MCMC` will create independent copies of a
#'    baseline passed `fmcmc_kernel` object. These are managed together in a
#'    `fmcmc_kernel_list` object.
#'    
#' 2. Even if the chains are run in parallel, if a predefined kernel object is
#'    passed it will be updated to reflect the last state of the kernels
#'    before the `MCMC` call returns.
#'    
#' @family kernels
#' @aliases fmcmc_kernel kernels
#' @examples 
#' 
#' # Example creating a multivariate normal kernel using the mvtnorm R package
#' # for a bivariate normal distribution
#' library(mvtnorm)
#' 
#' # Define your Sigma
#' sigma <- matrix(c(1, .2, .2, 1), ncol = 2)
#' 
#' # How does it looks like?
#' sigma
#' #      [,1] [,2]
#' # [1,]  1.0  0.2
#' # [2,]  0.2  1.0
#' 
#' # Create the kernel
#' kernel_mvn <- kernel_new(
#'   proposal = function(env) {
#'   env$theta0 + as.vector(mvtnorm::rmvnorm(1, mean = 0, sigma = sigma.))
#'   },
#'   sigma. = sigma
#' )
#' 
#' # As you can see, in the previous call we passed sigma as it will be used by
#' # the proposal function
#' # The logaratio function was not necesary to be passed since this kernel is
#' # symmetric.
#' 
#' @return An environment of class `fmcmc_kernel` which contains the following:
#' 
#' - `proposal` A function that receives a single argument, an environment. This
#'   is the proposal function used within [MCMC()].
#'   
#' - `logratio` A function to compute log ratios of the current vs the proposed
#'   step of the chain. Also used within [MCMC()].
#' 
#' - `...` Further arguments passed to `kernel_new`.
#' 
#' @references 
#' Brooks, S., Gelman, A., Jones, G. L., & Meng, X. L. (2011). Handbook of
#' Markov Chain Monte Carlo. Handbook of Markov Chain Monte Carlo.
#' @export
kernel_new <- function(
  proposal,
  ... ,
  logratio   = NULL,
  kernel_env = new.env(hash = TRUE)
  ) {
  
  # Checks
  if (length(formals(proposal)) != 1L)
    stop(
      "The `proposal` function should receive a single arguments (an environment).",
      call. = FALSE
      )
  if (!is.null(logratio) && length(formals(logratio)) != 1L)
    stop(
      "The `logratio` function should receive a single argument (an environment) ",
      call. = FALSE
      )
  
  if (is.null(logratio))
    logratio <- function(env) {env$f1 - env$f0}
  
  environment(proposal) <- kernel_env
  environment(logratio) <- kernel_env
  
  kernel_env$proposal <- proposal
  kernel_env$logratio <- logratio
  
  dots <- list(...)
  for (n in names(dots))
    kernel_env[[n]] <- dots[[n]]
  
  structure(kernel_env, class = c("environment", "fmcmc_kernel"))
  
}

#' @export
print.fmcmc_kernel <- function(x, ...) {
  
  cat("\nAn environment of class fmcmc_kernel:\n\n")
  print(utils::ls.str(x))
  cat("\n")
  invisible(x)
  
}

#' @export
c.fmcmc_kernel <- function(...) {
  
  dots <- list(...)
  if (any(!sapply(dots, inherits, what = "fmcmc_kernel")))
    stop("One or more terms passed in ... are not of class fmcmc_kernel.", call. = FALSE)
  
  kernels_env <- new.env(hash = TRUE)
  for (i in seq_along(dots))
    assign(x = paste0("kernel", i), value = dots[[i]], envir = kernels_env)
  
  structure(kernels_env, class = c("environment", "fmcmc_kernel_list"))
  
}

#' Creates independent duplicates of an fmcmc_kernel object
#' 
#' This is intended to be used in the case of having multiple chains running.s
#' @noRd
rep_kernel <- function(x, nchains, ...) {
  
  # If it is already duplicated, then we don't need to duplicate it again...
  # unless the number of kernels doesn't matches the number of chains!
  if (is_kernel_list(x)) 
    stop("The kernel is already a list (fmcmc_kernel_list).", call. = FALSE)

  if (!inherits(x, "fmcmc_kernel"))
    stop("x is not of class fmcmc_kernel.", call. = FALSE)
  
  # This is a list of environments, we need to add it to the original
  # environment
  ans <- replicate(n = nchains, expr = copy_kernel(x))
  
  # Cleaning and replacing
  rm(list = ls(envir = x), envir = x)
  
  ans <- lapply(
    X = seq_along(ans),
    FUN = function(i) assign(paste0("kernel", i), ans[[i]], envir = x)
    )
  
  class(x) <- c("environment", "fmcmc_kernel_list")
  
  invisible(NULL)
  
}

#' Copy a kernel by going env2list2env and wrap in fmcmc_kernel class
#' @noRd
#' 
copy_kernel <- function(x) {
  
  if (!inherits(x, "fmcmc_kernel"))
    stop("x must be an object of class fmcmc_kernel.", call. = FALSE)
  
  kernel_env <- new.env(hash = TRUE)
  lapply(names(x), function(obj_name) {
    assign(obj_name, get(obj_name, envir = x), kernel_env)
  })
  
  environment(kernel_env$proposal) <- kernel_env
  environment(kernel_env$logratio) <- kernel_env
  
  structure(kernel_env, class = c("environment", "fmcmc_kernel"))
  
}

#' @export
`[[.fmcmc_kernel_list` <- function(x, i) {
  
  get(names(x)[i], envir = x)
  
}

#' This function copies the contents of k1 into k0
#' @noRd
update_kernel <- function(k0, k1, positions = 1:length(k0)) {
  
  # Basic checks (to make sure we are not busting it)
  if (!is_kernel_list(k0))
    stop("k0 must be an object of class kernel_list.", call. = FALSE)
  
  if (!is_kernel_list(k1))
    stop("k1 must be an object of class kernel_list.", call. = FALSE)
  
  if (length(k0) != length(k1))
    stop("k0 and k1 must be of the same length.", call. = FALSE)
  
  # Assigning the values in the new environment
  for (i in positions) 
    assign(names(k1)[i], k1[[i]], envir = k0)
  
  invisible(NULL)
}

is_kernel_list <- function(x) {
  inherits(x, what = "fmcmc_kernel_list")
}

#' @export
print.fmcmc_kernel_list <- function(x, ...) {
  
  cat("A list of", length(x), "fmcmc_kernels.\n")
  invisible(x)
  
}



#' Reflective Boundaries
#' 
#' Adjust a proposal according to its support by reflecting it. This is the workhorse
#' of [kernel_normal_reflective] and [kernel_unif_reflective]. It is intended
#' for internal use only. 
#' 
#' @param x A numeric vector. The proposal
#' @param lb,ub Numeric vectors of length `length(x)`. Lower and upper bounds.
#' @param which Integer vector. Index of variables to be updated.
#' 
#' @return An adjusted proposal vector.
#' @export
reflect_on_boundaries <- function(
  x,
  lb,
  ub,
  which
) {
  
  test_above <- which(x[which] > ub[which])
  test_below <- which(x[which] < lb[which])
  
  d <- ub - lb
  
  if (length(test_above)) {
    
    # Putting in context
    test_above <- which[test_above]
    
    # Direction of update
    d_above <- x[test_above] - ub[test_above]
    odd     <- (d_above %/% d[test_above]) %% 2
    d_above <- d_above %% d[test_above]
    
    x[test_above] <- (lb[test_above] + d_above)*odd + 
      (ub[test_above] - d_above)*(1-odd)
  } 
  
  if (length(test_below)) {
    
    # Putting in context
    test_below <- which[test_below]
    
    # Direction of update
    d_below <- lb[test_below] - x[test_below]
    odd     <- (d_below %/% d[test_below]) %% 2
    d_below <- d_below %% d[test_below]
    
    x[test_below] <- (ub[test_below] - d_below)*odd + 
      (lb[test_below] + d_below)*(1-odd)
    
  } 
  
  x
  
}
