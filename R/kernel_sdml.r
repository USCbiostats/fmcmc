
#' State-Dependent Memory-Less Adaptive Transition Kernel
#' 
#' @param x Current state
#' @param f Objective function
#' @param gr (optional) Gradient of the objective function.
#' @param rz Random number function.
#' @param dz Fensity function of `z`.
#' @param rz.args List of parameters passed to `rz`.
#' @param ... Further arguments passed to `f` and, if specified, `gr`
#' @param alpha Number in `[0, 1)`.
#' @param eps If no `gr` is specified, epsilon used to approximate the gradient.
#' @export
sdml_adaptive_kernel <- function(
  x,
  f,
  gr,
  rz      = stats::rnorm,
  dz      = stats::dnorm,
  rz.args = list(mean=0, sd=1),
  ...,
  alpha = .9,
  eps = rep(1e-5, length(x))
) {
  
  # Number of parameters
  k <- length(x)
  
  # Evaluating functions
  fx  <- f(x, ...)

  # Gradient function
  gx <- if (missing(gr))
    (f(x + eps, ...) - fx)/
        (2 * eps)
  else
    gr(x, ...)
  
  
  # Proposal
  beta_x <- 1 - stats::plogis(gx)*alpha
  x_new  <- x + do.call(rz, c(list(k), rz.args))*beta_x
  fx_new <- f(x_new, ...)
  
  # New gradient
  gx_new <- if (missing(gr))
    (f(x_new + eps, ...) - fx_new)/
    (2 * eps)
  else
    gr(x_new, ...)
  
  # Hastings ratio
  beta_x_new <- 1 - stats::plogis(gx_new)*alpha
  fz_x_x_new <- do.call(dz, c(list((x_new - x)/beta_x), rz.args))/beta_x
  fz_x_new_x <- do.call(dz, c(list((x - x_new)/beta_x_new), rz.args))/beta_x_new
  
  list(
    x = x_new,
    h = min(1, fx_new * fz_x_x_new/(fx * fz_x_new_x + 1e-15))
  )
  
}
