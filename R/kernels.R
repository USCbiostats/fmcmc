#' Various kernel functions for MCMC
#' 
#' @param mu,scale Either a numeric vector or a scalar. Proposal mean and scale.
#' @param lb,ub Either a numeric vector or a scalar. Lower and upper bounds for
#' bounded kernels (currently the `kernel_reflective` only).
#' @param fixed Logical scalar or vector. When `TRUE` fixes the corresponding
#' parameter, avoiding new proposals.
#' @details
#' The objects `fmcmc_kernels` are environments that in general contain the 
#' following objects:
#' 
#' -  `proposal`: The function used to propose changes in the chain based
#'    on the current state. The function must return a vector of length equal
#'    to the number of parameters in the model.
#'    
#' -  `logration`: This function is called after a new state has been proposed,
#'    and is used to compute the log of the hastings ratio.
#'    
#'    In the case that the `logratio` function is not specified, then it is assumed
#'    that the transition kernel is symmetric, this is, logratio is then implemented
#'    as `function(env) {env$f1 - env$f0}`
#'    
#' -  `...`: Further objects that are used within those functions.
#' 
#' Both functions, `proposal` and `logratio`, receive a single argument, an
#' environment, which is passed by the [MCMC] function during each step using
#' the function [environment]. The passed
#' environment is actually the environemnt in which the `MCMC` function is running,
#' in particular, this environment contains the following objects:
#' 
#' \tabular{ll}{
#' `i` \tab Integer. The current iteration. \cr
#' `theta1` \tab Numeric vector. The last proposed state. \cr
#' `theta0` \tab Numeric vector. The current state \cr
#' `f` \tab The log-unnormalized posterior function (a wrapper of `fun` passed 
#' to [MCMC]). \cr
#' `f1` \tab The last value of `f(theta1)` \cr
#' `f0` \tab The last value of `f(theta0)` \cr
#' `kernel` \tab The actual `fmcmc_kernel` object.
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
#' For an extended example see the vignette "Personalized kernel functions".
#' 
#' 
#' @name kernels
#' @aliases fmcmc_kernel fmcmc-kernel
#' @examples 
#' \dontrun{
#' # Example creating a multivariate normal kernel using the mvtnorm R package
#' # for a bivariate normal distribution
#' library(mvtnorm)
#' 
#' # Define your Sigma
# sigma <- matrix(c(1, .2, .2, 1), ncol = 2)
#' 
#' # How does it looks like?
# sigma
#'      [,1] [,2]
#' [1,]  1.0  0.2
#' [2,]  0.2  1.0
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
#' }
NULL

#' @export
#' @rdname kernels
#' @param proposal,logratio Functions. The function receives a single argument, an environment.
#' This functions are called later within [MCMC] (see details).
#' @param ... In the case of `kernel_new`, further arguments to be stored with
#' the kernel.
#' @section Creating your own kernels:
#' The function `kernel_new` is a helper function that allows creating
#' `fmcmc_kernel` which is used with the `MCMC` function. The `fmcmc_kernel`
#' are the backbone of the [MCMC] function.
kernel_new <- function(proposal, logratio = NULL, ...) {
  
  # Checks
  if (length(formals(proposal)) != 1L)
    stop("The `proposal` function should receive a single argument.", call. = FALSE)
  if (!is.null(logratio) && length(formals(logratio)) != 1L)
    stop("The `logratio` function should receive a single argument.", call. = FALSE)
  
  if (is.null(logratio))
    logratio <- function(env) {env$f1 - env$f0}
  
  env <- new.env()
  environment(proposal) <- env
  environment(logratio) <- env
  
  env$proposal <- proposal
  env$logratio <- logratio
  
  dots <- list(...)
  for (n in names(dots))
    env[[n]] <- dots[[n]]
  
  structure(env, class = c("environment", "fmcmc_kernel"))
  
}

#' @export
#' @rdname kernels
#' @param x An object of class `fmcmc_kernel`.
print.fmcmc_kernel <- function(x, ...) {
  
  cat("\nAn object of class fmcmc_kernel. The following two functions are available:\n\n")
  print(utils::ls.str(x))
  cat("\n")
  invisible(x)
  
}

#' @export
#' @rdname kernels
#' @section Kernels:
#' The `kernel_unif` function provides a uniform transition kernel. This (symmetric)
#' kernel function by default adds the current status values between \[0,1\].
kernel_unif <- function(lb = -1.0, ub = 1.0) {
  kernel_new(
    proposal = function(env) env$theta0 + runif(1, lb, ub),
    lb = lb,
    ub = ub
  )
}

#' @export
#' @rdname kernels
#' @section Kernels:
#' The `kernel_normal` function provides the cannonical normal kernel
#' with symmetric transition probabilities.
kernel_normal <- function(mu = 0, scale = 1) {
  
  kernel_new(
    proposal = function(env) 
      env$theta0 + stats::rnorm(length(env$theta0), mean = mu, sd = scale),
    logratio = function(env) env$f1 - env$f0,
    mu       = mu,
    scale    = scale
    )
}

#' @export 
#' @rdname kernels
#' @section Kernels:
#' The `kernel_reflective` implements the normal kernel with reflective
#' boundaries. Lower and upper bounds are treated using reflecting boundaries, this is, 
#' if the proposed \eqn{\theta'} is greater than the \code{ub}, then \eqn{\theta' - ub}
#' is substracted from \eqn{ub}. At the same time, if it is less than \code{lb}, then
#' \eqn{lb - \theta'} is added to \code{lb} iterating until \eqn{\theta} is within
#' \code{[lb, ub]}.
#' 
#' In this case, the transition probability is symmetric (just like the normal
#' kernel).
kernel_reflective <- function(
  mu    = 0,
  scale = 1,
  lb    = -.Machine$double.xmax,
  ub    = .Machine$double.xmax,
  fixed = FALSE
  ) {
  
  k <- NULL
  
  kernel_new(
    proposal = function(env) {
      
      # Checking whether k exists. We should do this only once. When doing this
      # we have to make sure that the length of lb and ub are according to the
      # length of model parameter.
      #
      # We also want to restart k in the first run. Since it is based on
      # environments, the user may move this around... so it is better to just
      # restart this every time that the MCMC function starts from scratch.
      if (env$i == 1L | is.null(k)) {
        
        k <<- length(env$theta0)
        
        # Checking boundaries
        if (length(ub) > 1 && (k != length(ub)))
          stop("Incorrect length of -ub-", call. = FALSE)
        
        if (length(lb) > 1 && (k != length(lb)))
          stop("Incorrect length of -lb-", call. = FALSE)
        
        # Repeating boundaries
        if (length(ub) == 1)
          ub <<- rep(ub, k)
        
        if (length(lb) == 1)
          lb <<- rep(lb, k)
        
        if (any(ub <= lb))
          stop("-ub- cannot be <= than -lb-.", call. = FALSE)
        
        # mu, scale and fixed according to the number of parameters in the
        # model if needed. 
        if (length(mu) == 1)
          mu <<- rep(mu, k)
        if (length(scale) == 1)
          scale <<- rep(scale, k)
        if (length(fixed) == 1)
          fixed <<- rep(fixed, k)
        
      }
      
      normal_reflective(env$theta0, lb, ub, mu, scale, fixed)
      },
    logratio = function(env) env$f1 - env$f0,
    mu       = mu,
    scale    = scale, 
    ub       = ub, 
    lb       = lb, 
    fixed    = fixed,
    k        = k
  )
  
}

normal_reflective <- function(
  x,
  lb, ub,
  mu    = rep(0, length(x)),
  scale = rep(1, length(x)),
  fixed = rep(FALSE, length(x))
) {
  
  notfixed <- which(!fixed)
  x[notfixed] <- x[notfixed] + stats::rnorm(length(notfixed), mu) * scale[notfixed]
  
  test_above <- which(x[notfixed] > ub[notfixed])
  test_below <- which(x[notfixed] < lb[notfixed])
  
  d <- ub - lb
  
  if (length(test_above)) {
    
    # Direction of update
    d_above <- x[test_above] - ub[test_above]
    odd     <- (d_above %/% d[test_above]) %% 2
    d_above <- d_above %% d[test_above]
    
    x[test_above] <- (lb[test_above] + d_above)*odd + 
      (ub[test_above] - d_above)*(1-odd)
  } 
  
  if (length(test_below)) {
    
    # Direction of update
    d_below <- lb[test_below] - x[test_below]
    odd     <- (d_below %/% d[test_below]) %% 2
    d_below <- d_below %% d[test_below]
    
    x[test_below] <- (ub[test_below] - d_below)*odd + 
      (lb[test_below] + d_below)*(1-odd)
    
  } 
  
  x
  
}


