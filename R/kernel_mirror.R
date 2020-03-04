#' Mirror Transition Kernels
#' 
#' NMirror and UMirror transition kernels described in Thawornwattana et al.
#' (2018).
#' 
#' @template mu-scale
#' @template lb-ub
#' @template scheme
#' @param warmup Integer. Number of steps required before starting adapting the
#' chains.
#' @param nadapt Integer. Number of times the scale is adjusted for adaptation
#' during the warmup (burnin) period.
#' @param arate Double. Target acceptance rate used as a reference during the
#' adaptation process.
#' @family kernels
#' 
#' @details
#' 
#' The `kernel_nmirror` and `kernel_umirror` functions implement simple symmetric
#' transition kernels that pivot around an approximation of the asymptotic mean.
#' 
#' In the multidimensional case, this implementation just draws a vector of
#' independent draws from the proposal kernel, instead of using, for example,
#' a multivariate distribution of some kind. This will be implemented in the
#' next update of the package.
#' 
#' During the warmup period (or burnin as described in the paper), the algorithm
#' adapts both the scale and the reference mean of the proposal distribution.
#' While the mean is adapted continuously, the scale is updated only a handful
#' of times, in particular, `nadapt` times during the warmup time. The adaptation
#' is done as proposed by Yang and Rodriguez (2013) in which the
#' scale is adapted four times.
#' 
#' @return An object of class [fmcmc_kernel].
#' 
#' @references 
#' Thawornwattana, Y., Dalquen, D., & Yang, Z. (2018). Designing Simple and
#' Efficient Markov Chain Monte Carlo Proposal Kernels. Bayesian Analysis, 13(4),
#' 1037–1063. \url{https://doi.org/10.1214/17-BA1084}
#' 
#' Yang, Z., & Rodriguez, C. E. (2013). Searching for efficient Markov chain
#' Monte Carlo proposal kernels. Proceedings of the National Academy of Sciences,
#' 110(48), 19307–19312. \url{https://doi.org/10.1073/pnas.1311790110}
#' @name kernel_mirror
NULL

#' @export
#' @rdname kernel_mirror
kernel_nmirror <- function(
  mu     = 0,
  scale  = 1,
  warmup = 500L,
  nadapt = 4L,
  arate  = .4,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  fixed  = FALSE,
  scheme = "joint"
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
      
      # Updating the mean
      if ((abs_iter[[env$chain_id]] >= 1L) && (abs_iter[[env$chain_id]] <= warmup)) {
        
        mu[[env$chain_id]] <<- mean_recursive(
          X_t         = env$ans[env$i - 1L, , drop = FALSE],
          Mean_t_prev = mu[[env$chain_id]],
          t.          = abs_iter[[env$chain_id]]
        )
        
      }

      # Updating acceptance rate
      if (abs_iter[[env$chain_id]] == nadapt[1]) {
        
        obs_arate[env$chain_id] <<- list(
          1.0 - mean(rowSums(diff(env$ans[1L:(env$i - 1L), , drop = FALSE])^2) == 0.0)
        )
        
      } else if (abs_iter[[env$chain_id]] > nadapt[1] && abs_iter[[env$chain_id]] <= warmup) {
        
        obs_arate[[env$chain_id]] <<- mean_recursive(
          X_t         = as.double(env$ans[env$i - 1L, ] != env$ans[env$i - 2L, ]),
          Mean_t_prev = obs_arate[[env$chain_id]],
          t.          = abs_iter[[env$chain_id]]
        )
        
      }
      
      # Is this the iteration when we adapt?
      if (abs_iter[[env$chain_id]] %in% nadapt) {
        
        scale[[env$chain_id]] <<- scale[[env$chain_id]] *
          tan(pi/2.0 * obs_arate[[env$chain_id]]) /
          tan(pi/2.0 * arate)
        
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
    nadapt    = floor(seq(1L, warmup, length.out = nadapt + 1)[-1]),
    arate     = arate,
    k         = k,
    fixed     = fixed,
    scheme    = scheme,
    update_sequence = update_sequence
    )
  
}

#' @export
#' @rdname kernel_mirror
kernel_umirror <- function(
  mu     = 0,
  scale  = 1,
  warmup = 500L,
  nadapt = 4L,
  arate  = .4,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  fixed  = FALSE,
  scheme = "joint"
) {
  
  k      <- NULL
  sqrt3  <- sqrt(3)
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
      
      # Updating the mean
      if ((abs_iter[[env$chain_id]] >= 1L) && (abs_iter[[env$chain_id]] <= warmup)) {
        
        mu[[env$chain_id]] <<- mean_recursive(
          X_t         = env$ans[env$i - 1L, , drop = FALSE],
          Mean_t_prev = mu[[env$chain_id]],
          t.          = abs_iter[[env$chain_id]]
        )
        
      }
      
      # Updating acceptance rate
      if (abs_iter[[env$chain_id]] == nadapt[1]) {
        
        obs_arate[env$chain_id] <<- list(
          1.0 - mean(rowSums(diff(env$ans[1L:(env$i - 1L), , drop = FALSE])^2) == 0.0)
        )
        
      } else if (abs_iter[[env$chain_id]] > nadapt[1] && abs_iter[[env$chain_id]] <= warmup) {
        
        obs_arate[[env$chain_id]] <<- mean_recursive(
          X_t         = as.double(env$ans[env$i - 1L, ] != env$ans[env$i - 2L, ]),
          Mean_t_prev = obs_arate[[env$chain_id]],
          t.          = abs_iter[[env$chain_id]]
        )
        
      }
      
      # Is this the iteration when we adapt?
      if (abs_iter[[env$chain_id]] %in% nadapt) {
        
        scale[[env$chain_id]] <<- scale[[env$chain_id]] *
          tan(pi/2.0 * obs_arate[[env$chain_id]]) /
          tan(pi/2.0 * arate)
        
      }
      
      # Which to update (only whichever is to be updated in the sequence)
      # of updates.
      which. <- which(update_sequence[env$i, ])
      
      # Making the proposal
      theta1 <- env$theta0
      
      # Making the proposal
      theta1[which.] <- stats::runif(
        k,
        min = 2*mu[[env$chain_id]] - theta1[which.] - sqrt3 * scale[[env$chain_id]],
        max = 2*mu[[env$chain_id]] - theta1[which.] + sqrt3 * scale[[env$chain_id]]
        )
      
      # Increasing the absolute number of iteration
      abs_iter[[env$chain_id]] <<- abs_iter[[env$chain_id]] + 1L
      
      reflect_on_boundaries(theta1, lb, ub, which.)
      
    },
    mu        = mu,
    mu0       = mu0,
    scale     = scale,
    scale0    = scale0,
    abs_iter  = abs_iter,
    obs_arate = obs_arate,
    lb        = lb,
    ub        = ub,
    warmup    = warmup,
    nadapt    = floor(seq(1L, warmup, length.out = nadapt + 1)[-1]),
    arate     = arate,
    k         = k,
    fixed     = fixed,
    scheme    = scheme,
    update_sequence = update_sequence,
    sqrt3     = sqrt3
  )
  
}
