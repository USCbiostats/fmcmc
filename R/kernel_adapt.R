#' Adaptive Metropolis (AM) Transition Kernel
#' 
#' Implementation of Haario et al. (2001)'s Adaptive Metropolis.
#' 
#' @param fixed Logical scalar or vector of length `k`. Indicates which parameters
#' will be treated as fixed or not. Single values are recycled.
#' @template lb-ub
#' @template mu-Sigma
#' @template until
#' @param bw Integer scalar. The bandwidth, is the number of observations to
#' include in the computation of the variance-covariance matrix.
#' @param freq Integer scalar. Frequency of updates. How often the
#' variance-covariance matrix is updated. The implementation is different from that
#' described in the original paper (see details).
#' @param Sd Overall scale for the algorithm. By default, the variance-covariance
#' is scaled to \eqn{2.4^2/d}, with \eqn{d} the number of dimensions.
#' 
#' @details
#' 
#' `kernel_adapt` Implements the adaptive Metropolis (AM) algorithm of Haario
#' et al. (2001). If the value of bw is greater than zero, then the algorithm
#' folds back AP, a  previous version which is known to have ergodicity problems.
#' 
#' The parameter `eps` has two functions. The first one is to set the initial
#' scale for the multivariate normal kernel, which is replaced after `warmup`
#' steps with the actual variance-covariance computed by the main algorithm.
#' The second usage is in the equation that ensures that the variance-covariance
#' is greater than zero, this is, the \eqn{\varepsilon}{epsilon} parameter in the
#' original paper.
#' 
#' The update of the covariance matrix is done using [cov_recursive()] function,
#' which makes the updates faster. The `freq` parameter, besides of indicating the
#' frequency with which the updates are done, it specifies what are the samples
#' included in each update, in other words, like a thinning parameter, only every
#' `freq` samples will be used to compute the covariance matrix. Since this
#' implementation uses the recursive formula for updating the covariance, there is
#' no practical need to set `freq != 1`.
#' 
#' @references 
#' Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive Metropolis algorithm.
#' Bernoulli, 7(2), 223–242.
#' \url{https://projecteuclid.org/euclid.bj/1080222083}
#' 
#' @template kernel
#' 
#' @export 
#' @examples 
#' # Update every-step and wait 1,000 steps before starting to adapt
#' kern <- kernel_adapt(freq = 1, warmup = 1000)
#' 
#' # Two parameters model, the second parameter with a restricted range, i.e.
#' # a lower bound of 1
#' kern <- kernel_adapt(lb = c(-.Machine$double.xmax, 0))
kernel_adapt <- function(
  mu     = 0,
  bw     = 0L,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  freq   = 1L,
  warmup = 500L,
  Sigma  = NULL,
  Sd     = NULL,
  eps    = 1e-4,
  fixed  = FALSE,
  until  = Inf
) {
  
  k      <- NULL
  sd     <- NULL
  Ik     <- NULL
  which. <- NULL
  
  # Variables for fast cov (see cov_recursive)
  Mean_t_prev <- NULL
  # t.          <- NULL
  
  # We create copies of this b/c each chain has its own values
  abs_iter <- 0L
  
  if (bw > 0L && bw > warmup)
    stop("The `warmup` parameter must be greater than `bw`.", call. = FALSE)
  
  kernel_new(
    proposal = function(env) {
      
      # In the case of the first iteration
      if (env$i == 1L | is.null(k)) {
        
        k     <<- length(env$theta0)

        mu     <<- check_dimensions(mu, k)
        ub     <<- check_dimensions(ub, k)
        lb     <<- check_dimensions(lb, k)
        fixed  <<- check_dimensions(fixed, k)
        which. <<- which(!fixed)
        
        # Need to adapt this afterwards
        k  <<- sum(!fixed)
        Ik <<- diag(k) * eps
        mu <<- mu[which.]
        
        if (is.null(Sigma))
          Sigma <<- Ik
        
        if (any(ub <= lb))
          stop("-ub- cannot be <= than -lb-.", call. = FALSE)
        
        if (is.null(Sd))
          Sd <<- 5.76 / k # 2.38^2
      }
      
      # Updating the scheme
      if ((until > abs_iter) && (abs_iter > warmup) && (env$i > 2L) && !(env$i %% freq)) {
        
        ran <- if (bw <= 0L) 1L:(env$i - 1L) 
        else (env$i - bw + 1L):(env$i - 1L)
        
        if (bw > 0L) {
          
          Sigma <<- Sd * (stats::cov(env$ans[ran, which., drop = FALSE]) + Ik)
          
        } else {
          
          # Update mean
          if (is.null(Mean_t_prev))
            Mean_t_prev <<- colMeans(env$ans[1:(env$i - 1), which., drop = FALSE])
          
          # Figuring out range (this is greater than a single observation if the
          # updates happen with freq > 1).
          update_range <- (env$i - freq):(env$i - 1L)
          
          # if (is.null(t.))
          #   t. <<- update_range[1L]

          # If we are updating every freq != 1, then this is a matrix
          Mean_t <- mean_recursive(
            X_t         = env$ans[update_range, which., drop = FALSE],
            Mean_t_prev = Mean_t_prev,
            t.          = abs_iter - freq
            )
          
          # Update sigma
          Sigma <<- cov_recursive(
            X_t         = env$ans[update_range, which., drop = FALSE],
            Cov_t       = Sigma,
            Mean_t      = Mean_t,
            Mean_t_prev = Mean_t_prev,
            t.          = abs_iter - freq,
            eps         = 1e-5,
            Ik          = Ik
          )[, , freq, drop = TRUE]
          
          if (freq > 1L) 
            Mean_t_prev <<- Mean_t[freq, ]
          else
            Mean_t_prev <<- Mean_t
          
          # Add one to the counter, we need this for the recursive update of the
          # Covariance matrix.
          # t. <<- t. + freq
        }
      }
      
      # Increasing the absolute number of iteration
      abs_iter <<- abs_iter + 1L
      
      # Making the proposal
      theta1 <- env$theta0
      theta1[which.] <- env$theta0[which.] +
        MASS::mvrnorm(
          mu    = mu,
          Sigma = Sigma
          )
      
      reflect_on_boundaries(theta1, lb, ub, which.)
      
    },
    mu          = mu,
    bw          = bw,
    lb          = lb,
    ub          = ub,
    freq        = freq,
    warmup      = warmup,
    Sigma       = Sigma,
    Sd          = Sd,
    eps         = eps,
    fixed       = fixed,
    k           = k,
    Ik          = Ik,
    which.      = which.,
    Mean_t_prev = Mean_t_prev,
    # t.          = t.,
    abs_iter    = abs_iter,
    until       = until
  )
  
}

#' @export
#' @rdname kernel_adapt
#' @details 
#' `kernel_am` is just an alias for `kernel_adapt`.
kernel_am <- kernel_adapt


