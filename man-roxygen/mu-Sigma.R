#' @param mu Either a numeric vector or a scalar. Proposal mean.
#' If scalar, values are recycled to match the number of parameters in the
#' objective function.
#' @param Sigma The variance-covariance matrix. By default this will be an
#' identity matrix during the warmup period.
#' @param eps Double scalar. Default size of the initial step (see details).
#' @param warmup Integer scalar. The number of iterations that the algorithm has
#' to wait before starting to do the updates.
NULL

