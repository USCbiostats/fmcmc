#' Various kernel functions for MCMC
#' 
#' @param mean,scale Either a numeric vector or a scalar. Proposal scale.
#' @param lb,ub Either a numeric vector or a scalar. Lower and upper bounds for
#' bounded kernels (currently the `kernel_reflective` only).
#' @param fixed Logical scalar or vector. When `TRUE` fixes the corresponding
#' parameter, avoiding new proposals.
#' @name kernels
#' @aliases amcmc_kernel amcmc-kernel
#' @examples 
#' \dontrun{
#' # Example creating a multivariate normal kernel using the mvtnorm R package
#' library(mvtnorm)
#' 
#' # Define your Sigma
#' sigma <- ...
#' 
#' # Create the kernel
#' mvn_kernel <- kernel_new(
#'   proposal = function(env) {
#'   env$theta0 + as.vector(mvtnorm::rmvnorm(1, mean = 0, sigma = sigma, ...))
#'   },
#'   logratio = function(env) {
#'     env$f1 - env$f0 +
#'       log(mvtnorm::pmvnorm(upper = env$theta0 - env$theta1)) -
#'       log(mvtnorm::pmvnorm(upper = env$theta1 - env$theta0))
#'   },
#'   sigma = sigma
#' )
#' 
#' }
NULL

#' @export
#' @rdname kernels
#' @param proposal Function. The function receives a single argument, an environment.
#' This function is later called in [MCMC] to propose changes in the chain based
#' on the current state. The function must return a vector of length equal
#' to the number of parameters in the model.
#' @param logratio Function. This function receives a single argument, an environment.
#' This function is called after a new state has been proposed, and is used to
#' compute the log of the hastings ratio.
#' @param ... In the case of `kernel_new`, further arguments to be stored with
#' the kernel.
#' 
#' @details The function `kernel_new` is a helper function that allows creating
#' `amcmc_kernel` which is used with the `MCMC` function.
kernel_new <- function(proposal, logratio, ...) {
  
  env <- new.env()
  environment(proposal) <- env
  environment(logratio) <- env
  
  env$proposal <- proposal
  env$logratio <- logratio
  
  dots <- list(...)
  for (n in names(dots))
    env[[n]] <- dots[[n]]
  
  structure(env, class = c("environment", "amcmc_kernel"))
  
}

#' @export
#' @rdname kernels
#' @param x An object of class `amcmc_kernel`.
print.amcmc_kernel <- function(x, ...) {
  
  cat("\nAn object of class amcmc_kernel. The following two functions are available:\n\n")
  print(utils::ls.str(x))
  cat("\n")
  invisible(x)
  
}

#' @export
#' @rdname kernels
#' @details The `kernel_normal` function provides the cannonical normal kernel
#' with symmetric transition probabilities.
kernel_normal <- function(mean = 0, scale = 1) {
  
  kernel_new(
    proposal = function(env) 
      env$theta0 + stats::rnorm(length(env$theta0), mean = mean, sd = scale),
    logratio = function(env) env$f1 - env$f0,
    mean  = mean,
    scale = scale
    )
}

#' @export 
#' @rdname kernels
#' @details The `kernel_reflective` implements the normal kernel with reflective
#' boundaries. Lower and upper bounds are treated using reflecting boundaries, this is, 
#' if the proposed \eqn{\theta'} is greater than the \code{ub}, then \eqn{\theta' - ub}
#' is substracted from \eqn{ub}. At the same time, if it is less than \code{lb}, then
#' \eqn{lb - \theta'} is added to \code{lb} iterating until \eqn{\theta} is within
#' \code{[lb, ub]}.
#' 
#' In this case, the transition probability is symmetric (just like the normal
#' kernel).
kernel_reflective <- function(
  scale = 1,
  lb    = -.Machine$double.xmax,
  ub    = .Machine$double.xmax,
  fixed = FALSE
  ) {

  kernel_new(
    proposal = function(env) {
      
      # Checking whether k exists. We should do this only once. When doing this
      # we have to make sure that the length of lb and ub are according to the
      # length of model parameter.
      #
      # We also want to restart k in the first run. Since it is based on
      # environments, the user may move this around... so it is better to just
      # restart this every time that the MCMC function starts from scratch.
      if (env$i == 1L | !length(environment(proposal)[["k"]])) {
        
        assign("k", length(env$theta0), envir = environment(proposal))
        
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
        
        # Repeating 
        if (length(fixed) == 1)
          fixed <<- rep(fixed, k)
        
        # Repeating scale
        if (length(scale) == 1)
          scale <<- rep(scale, k)
        
      }
      
      normal_prop(env$theta0, lb, ub, scale, fixed)
      },
    logratio = function(env) env$f1 - env$f0,
    scale = scale, ub = ub, lb = lb, fixed = fixed
  )
  
}


