#' @importFrom Rcpp evalCpp
#' @importFrom coda mcmc mcmc.list
#' @importFrom parallel makePSOCKcluster stopCluster clusterExport clusterEvalQ
#'   clusterApply detectCores clusterSetRNGStream
#' @useDynLib amcmc
#' @importFrom methods formalArgs
#' @importFrom stats runif
NULL

#' Adaptative MCMC
#' 
#' @docType package
#' @name amcmc
#' 
NULL