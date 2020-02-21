#' @rdname kernels
#' @export
#' @section Kernels:
#' `kernel_nmirror` and `kernel_umirror` implemented in Thawornwattana et al.
#' (2018). Provides simple symmetric transition kernels that pivot around
#' an approximation of the asymptotic mean. The `warmup` parameter here sets
#' the number of steps before replacing `mu` and `scale` (sd) with those obtained
#' from the sample until the `warmup` number of iteration.
#' 
#' In the multidimensional case, this implementation just draws a vector of
#' independent draws from the proposal kernel, instead of using, for example,
#' a multivariate distribution of some kind. This will be implemented in the
#' next update of the package.
#' @references 
#' Thawornwattana, Y., Dalquen, D., & Yang, Z. (2018). Designing Simple and
#' Efficient Markov Chain Monte Carlo Proposal Kernels. Bayesian Analysis, 13(4),
#' 1037–1063. \url{https://doi.org/10.1214/17-BA1084}
kernel_nmirror <- function(
  mu     = 0,
  scale  = 1,
  warmup = 100,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  fixed  = FALSE
) {
  
  k      <- NULL
  which. <- NULL
  
  kernel_new(
    proposal = function(env) {
      
      # In the case of the first iteration
      if (env$i == 1L | is.null(k)) {
        
        k      <<- length(env$theta0)
        mu     <<- check_dimensions(mu, k)
        scale  <<- check_dimensions(scale, k)
        ub     <<- check_dimensions(ub, k)
        lb     <<- check_dimensions(lb, k)
        fixed  <<- check_dimensions(fixed, k)
        which. <<- which(!fixed)
      }
      
      # Computing the mu and scale
      if ((warmup > 0L) && (env$i > warmup) && (env$i > 2)) {
        
        # Once passed, we can set the warmup to a negative number so that
        # way, when wrapped around convergence checker, it won't re-update.
        warmup <<- structure(-999L, original = c(warmup = warmup))
        
        mu    <<- colMeans(env$ans[1L:(env$i - 1L), , drop=FALSE])
        scale <<- apply(env$ans[1L:(env$i - 1L), , drop=FALSE], 2L, stats::sd)
      }
      
      # Making the proposal
      theta1 <- stats::rnorm(
        k,
        mean = 2*mu - env$theta0,
        sd   = scale
        )
      
      reflect_on_boundaries(theta1, lb, ub, which.)
      
    },
    mu     = mu,
    scale  = scale,
    lb     = lb,
    ub     = ub,
    warmup = warmup,
    k      = k,
    fixed  = fixed,
    which. = which.
    )
  
}

#' @export
#' @rdname kernels
kernel_umirror <- function(
  mu     = 0,
  scale  = 1,
  warmup = 100,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  fixed  = FALSE
) {
  
  k      <- NULL
  which. <- NULL
  sqrt3  <- sqrt(3)
  
  kernel_new(
    proposal = function(env) {
      
      # In the case of the first iteration
      if (env$i == 1L | is.null(k)) {
        
        k      <<- length(env$theta0)
        mu     <<- check_dimensions(mu, k)
        scale  <<- check_dimensions(scale, k)
        ub     <<- check_dimensions(ub, k)
        lb     <<- check_dimensions(lb, k)
        fixed  <<- check_dimensions(fixed, k)
        which. <<- which(!fixed)
      }
      
      # Computing the mu and scale
      if ((warmup > 0L) && (env$i > warmup) && (env$i > 2)) {
        
        # Once passed, we can set the warmup to a negative number so that
        # way, when wrapped around convergence checker, it won't re-update.
        warmup <<- structure(-999L, original = c(warmup = warmup))
        
        mu    <<- colMeans(env$ans[1L:(env$i - 1L), , drop=FALSE])
        scale <<- apply(env$ans[1L:(env$i - 1L), , drop=FALSE], 2L, stats::sd)
        
      }
      
      # Making the proposal
      theta1 <- stats::runif(
        k,
        min = 2*mu - env$theta0 - sqrt3 * scale,
        max = 2*mu - env$theta0 + sqrt3 * scale
        )
      
      reflect_on_boundaries(theta1, lb, ub, which.)
      
    },
    mu     = mu,
    scale  = scale,
    lb     = lb,
    ub     = ub,
    warmup = warmup,
    k      = k,
    fixed  = fixed,
    which. = which.,
    sqrt3  = sqrt3
  )
  
}
