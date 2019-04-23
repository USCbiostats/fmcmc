#' @importFrom Rcpp evalCpp
#' @importFrom coda mcmc mcmc.list
#' @importFrom parallel makePSOCKcluster stopCluster clusterExport clusterEvalQ
#'   clusterApply detectCores clusterSetRNGStream
#' @useDynLib fmcmc
#' @importFrom methods formalArgs
#' @importFrom stats runif plogis
NULL

#' A friendly MCMC framework
#' 
#' @docType package
#' @name fmcmc
#' 
NULL