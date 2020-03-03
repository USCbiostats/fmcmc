#' Mirror Transition Kernels
#' @template mu-scale
#' @template lb-ub
#' @param warmup Integer. Number of steps required before starting adapting the
#' chains.
#' @param tadapt Integer. Number of steps during which the adaptation will take
#' place.
#' @param arate Double. Target acceptance rate used as a reference during the
#' adaptation process.
#' @family kernels
#' @details
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
#' @name kernel_mirror
NULL

#' @export
#' @rdname kernel_mirror
#' @family kernels
kernel_nmirror <- function(
  mu     = 0,
  scale  = 1,
  warmup = 100,
  tadapt = 100,
  arate  = .4,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  fixed  = FALSE,
  scheme = "joint",
  eps    = 1e-5
) {
  
  k <- NULL
  update_sequence <- NULL
  
  # We create copies of this b/c each chain has its own values
  mu0       <- mu
  mu        <- list()
  scale0    <- scale
  scale     <- list()
  abs_iter  <- list()
  obs_arate <- list()
  
  kernel_new(
    proposal = function(env) {
      
      # In the case of the first iteration
      if (env$i == 1L | is.null(k)) {
        
        k      <<- length(env$theta0)
        ub     <<- check_dimensions(ub, k)
        lb     <<- check_dimensions(lb, k)
        fixed  <<- check_dimensions(fixed, k)
        
        # Looking at this chain in particular, we grow the thing
        if ( is.null(abs_iter[env$chain_id][[1L]]) ) {
          
          mu[env$chain_id]        <<- list(check_dimensions(mu0, k)) 
          scale[env$chain_id]     <<- list(check_dimensions(scale0, k))
          abs_iter[env$chain_id]  <<- list(0L)

        }
        
        # Setting the scheme in which the variables will be updated
        update_sequence <<- plan_update_sequence(
          k      = k,
          nsteps = env$nsteps,
          fixed  = fixed,
          scheme = scheme
        )
        
        k <<- sum(update_sequence[1L, ])
        
      }
      
      # Computing the mu and scale
      if (abs_iter[[env$chain_id]] >= warmup && abs_iter[[env$chain_id]] <= tadapt) {

        # Updating the mean
        if (abs_iter[[env$chain_id]] == warmup) {
          
          mu[env$chain_id] <<- list(colMeans(env$ans[1:(env$i - 1), , drop = FALSE]))
          
        } else {
          
          mu[[env$chain_id]] <<- mean_recursive(
            X_t         = env$ans[env$i - 1L, , drop = FALSE],
            Mean_t_prev = mu[[env$chain_id]],
            t.          = abs_iter[[env$chain_id]]
          )
          
        }
        
        # Updating acceptance rate
        if (abs_iter[[env$chain_id]] == warmup) {
          obs_arate[env$chain_id] <<- list(
            1.0 - mean(rowSums(diff(env$ans[1L:(env$i - 1L), , drop = FALSE])^2) == 0.0)
            )
        } else {
          obs_arate[[env$chain_id]] <<- mean_recursive(
            X_t         = as.double(env$ans[env$i - 1L, ] != env$ans[env$i - 2L, ]),
            Mean_t_prev = obs_arate[[env$chain_id]],
            t.          = abs_iter[[env$chain_id]]
          )
        }
        
        # Updating the scale
        scale[[env$chain_id]] <<- scale[[env$chain_id]] *
          tan(pi/2.0 * obs_arate[[env$chain_id]]) /
          tan(pi/2.0 * arate) + eps
        
        # if (env$i>1100) {
        #   warning("upsi!")
        # }
        
      }
      
      # Which to update (only whichever is to be updated in the sequence)
      # of updates.
      which. <- which(update_sequence[env$i, ])
      
      # Making the proposal
      theta1 <- env$theta0
      
      theta1[which.] <- stats::rnorm(
        k,
        mean = 2*mu[[env$chain_id]][which.] - env$theta0[which.],
        sd   = scale[[env$chain_id]][which.]
        )
      
      # Increasing the absolute number of iteration
      abs_iter[[env$chain_id]] <<- abs_iter[[env$chain_id]] + 1L
      
      reflect_on_boundaries(x = theta1, lb = lb, ub = ub, which = which.)
      
    },
    mu0       = mu0,
    mu        = mu,
    scale0    = scale0,
    scale     = scale,
    abs_iter  = abs_iter,
    obs_arate = obs_arate,
    lb        = lb,
    ub        = ub,
    warmup    = warmup,
    tadapt    = tadapt,
    arate     = arate,
    k         = k,
    fixed     = fixed,
    scheme    = scheme,
    update_sequence = update_sequence,
    eps       = eps
    )
  
}

#' @export
#' @rdname kernel_mirror
#' @family kernels
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
