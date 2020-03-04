#' @export 
#' @rdname kernels
#' @param bw Integer scalar. The bandwidth, is the number of observations to
#' include in the computation of the variance-covariance matrix.
#' @param freq Integer scalar. Frequency of updates. How often the
#' variance-covariance matrix is updated.
#' @param warmup Integer scalar. The number of iterations that the algorithm has
#' to wait before starting to do the updates.
#' @param Sigma The variance-covariance matrix. By default this will be an
#' identity matrix during the warmup period.
#' @param eps Double scalar. Default size of the initial step (see details).
#' @param Sd Overall scale for the algorithm. By default, the variance-covariance
#' is scaled to \eqn{2.4^2/d}, with \eqn{d} the number of dimensions.
#' @section Kernels:
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
#' which makes the updates faster.
#' 
#' @references 
#' Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive Metropolis algorithm.
#' Bernoulli, 7(2), 223–242.
#' \url{https://projecteuclid.org/euclid.bj/1080222083}
#' 
kernel_adapt <- function(
  mu     = 0,
  bw     = 0L,
  lb     = -.Machine$double.xmax,
  ub     = .Machine$double.xmax,
  freq   = 50L,
  warmup = 500L,
  Sigma  = NULL,
  Sd     = NULL,
  eps    = 1e-4,
  fixed  = FALSE
) {
  
  k      <- NULL
  sd     <- NULL
  Ik     <- NULL
  which. <- NULL
  
  # Variables for fast cov (see cov_recursive)
  Mean_t_prev <- list()
  t.          <- list()
  
  # We create copies of this b/c each chain has its own values
  Sigma0  <- Sigma
  Sigma   <- list()
  
  abs_iter <- list()
  
  if (bw > 0L && bw > warmup)
    stop("The `warmup` parameter must be greater than `bw`.", call. = FALSE)
  
  kernel_new(
    proposal = function(env) {
      
      # In the case of the first iteration
      if (env$i == 1L | is.null(k)) {
        
        k     <<- length(env$theta0)
        Ik    <<- diag(k) * eps
        
        if (is.null(Sigma0))
          Sigma0 <<- Ik
        
        # Looking at this chain in particular, we grow the thing
        if ( is.null(Sigma[env$chain_id][[1L]]) ) {
          
          Sigma[env$chain_id]    <<- list(Sigma0) # Ik * eps
          t.[env$chain_id]       <<- list(1L) # Just starting
          abs_iter[env$chain_id] <<- list(0L)
          
        }
        
        mu     <<- check_dimensions(mu, k)
        ub     <<- check_dimensions(ub, k)
        lb     <<- check_dimensions(lb, k)
        fixed  <<- check_dimensions(fixed, k)
        which. <<- which(!fixed)
        
        if (any(ub <= lb))
          stop("-ub- cannot be <= than -lb-.", call. = FALSE)
        
        if (is.null(Sd))
          Sd <<- 5.76 / length(which.) # 2.38^2
      }
      
      # Updating the scheme
      if (abs_iter[[env$chain_id]] > warmup && env$i > 2L && !(env$i %% freq)) {
        
        ran <- if (bw <= 0L) 1L:(env$i - 1L) 
        else (env$i - bw + 1L):(env$i - 1L)
        
        if (bw > 0L) {
          
          Sigma[[env$chain_id]] <<- Sd * (stats::cov(env$ans[ran, , drop = FALSE]) + Ik)
          
        } else {
          
          # Update mean
          if (is.null(Mean_t_prev[env$chain_id][[1L]]))
            Mean_t_prev[env$chain_id] <<- list(colMeans(env$ans[1:(env$i - 1), ]))
          
          Mean_t <- mean_recursive(
            X_t         = env$ans[env$i - 1L, ],
            Mean_t_prev = Mean_t_prev[[env$chain_id]],
            t.          = t.[[env$chain_id]]
            )
          
          # Update sigma
          Sigma[[env$chain_id]] <<- cov_recursive(
            X_t         = env$ans[env$i - 1, ],
            Cov_t       = Sigma[[env$chain_id]],
            Mean_t      = Mean_t,
            Mean_t_prev = Mean_t_prev[[env$chain_id]],
            t.          = t.[[env$chain_id]],
            eps         = 1,
            Ik          = Ik
            )
          
          Mean_t_prev[[env$chain_id]] <<- Mean_t
          
          # Add one to the counter, we need this for the recursive update of the
          # Covariance matrix.
          t.[[env$chain_id]] <<- t.[[env$chain_id]] + 1L
        }
      }
      
      # Increasing the absolute number of iteration
      abs_iter[[env$chain_id]] <<- abs_iter[[env$chain_id]] + 1L
      
      # Making the proposal
      theta1 <- env$theta0
      theta1[which.] <- env$theta0[which.] +
        MASS::mvrnorm(
          mu    = mu[which.],
          Sigma = Sigma[[env$chain_id]][which., , drop = FALSE][,which. , drop = FALSE]
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
    Sigma0      = Sigma0,
    Sd          = Sd,
    eps         = eps,
    fixed       = fixed,
    k           = k,
    Ik          = Ik,
    which.      = which.,
    Mean_t_prev = Mean_t_prev,
    t.          = t.,
    abs_iter    = abs_iter
  )
  
}


#' Recursive algorithms for computing variance and mean
#' 
#' These algorithms are used in [kernel_adapt()] to simplify variance-covariance
#' recalculation at every step of the algorithm.
#' 
#' @param X_t Last value of the sample
#' @param Cov_t Covariance in t
#' @param Mean_t,Mean_t_prev Vectors of averages in time `t` and `t-1` respectively.
#' @param t. Sample size up to `t-1`.
#' @param Sd,eps,Ik See [kernel_adapt()].
#' @export
#' @details
#' The variance covariance algorithm was described in Haario, Saksman and
#' Tamminen (2002). 
#' 
#' @references 
#' Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive Metropolis algorithm.
#' Bernoulli, 7(2), 223–242.
#' \url{https://projecteuclid.org/euclid.bj/1080222083}
#' @examples 
#' # Generating random data (only four points to see the difference)
#' set.seed(1231)
#' n <- 3
#' X <- matrix(rnorm(n*4), ncol = 4)
#' 
#' # These two should be equal
#' mean_recursive(
#'   X_t         = X[1, ],
#'   Mean_t_prev = colMeans(X[-1, ]),
#'   t.          = n - 1
#' )
#' colMeans(X)
#' 
#' # These two should be equal
#' cov_recursive(
#'   X_t         = X[1, ], 
#'   Cov_t       = cov(X[-1, ]), 
#'   Mean_t      = colMeans(X),
#'   Mean_t_prev = colMeans(X[-1, ]),
#'   t           = n - 1
#' )
#' cov(X)
#' 
#' # Speed example -------------------------------------------------------------
#' set.seed(13155511)
#' X <- matrix(rnorm(1e3*100), ncol = 100)
#' 
#' ans0 <- cov(X[-1, ])
#' t0 <- system.time({
#'   ans1 <- cov(X)
#' })
#' 
#' t1 <- system.time(ans2 <- cov_recursive(
#'   X[1, ], ans0,
#'   Mean_t      = colMeans(X),
#'   Mean_t_prev = colMeans(X[-1, ]),
#'   t. = 1e3 - 1
#' ))
#' 
#' # Comparing accuracy and speed
#' range(ans1 - ans2)
#' t0/t1
#' 
cov_recursive <- function(
  X_t,
  Cov_t,
  Mean_t,
  Mean_t_prev,
  t.,
  eps = 0, Sd = 1, Ik = diag(length(X_t))
  ) {
  
  (t. - 1)/t. * Cov_t + 
    Sd/t. * (
      t. * tcrossprod(Mean_t_prev) -
        (t. + 1) * tcrossprod(Mean_t) +
        tcrossprod(X_t) + 
        eps * Ik
      )
  
}

#' @export
#' @rdname cov_recursive
mean_recursive <- function(X_t, Mean_t_prev, t.) {
  
  (Mean_t_prev * t. + X_t)/ (t. + 1)
  
}


