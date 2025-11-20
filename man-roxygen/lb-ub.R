#' @param lb,ub Either a numeric vector or a scalar. Lower and upper bounds for
#' bounded kernels. When of length 1, the values are recycled to match the number
#' of parameters in the objective function. Use `NA` to indicate unbounded parameters
#' (equivalent to `-.Machine$double.xmax` for lower bounds and `.Machine$double.xmax`
#' for upper bounds). Named vectors are supported for clarity.
NULL

