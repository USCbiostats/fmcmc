#' @importFrom coda mcmc mcmc.list
#' @importFrom parallel makePSOCKcluster stopCluster clusterExport clusterEvalQ
#'   clusterApply detectCores clusterSetRNGStream
#' @importFrom methods formalArgs
#' @importFrom stats runif plogis
#' @importFrom MASS mvrnorm
#' @importFrom utils head
NULL

#' A friendly MCMC framework
#' 
#' The `fmcmc` package provides a flexible framework for implementing MCMC models
#' using a lightweight in terms of dependencies. Among its main features, `fmcmc`
#' allows:
#' 
#' - Implementing arbitrary transition kernels.
#' 
#' - Incorporating convergence monitors for automatic stop.
#' 
#' - Out-of-the-box parallel computing implementation for running multiple chains
#' simultaneously.
#' 
#' For more information see the packages vignettes:
#' 
#' ```
#' vignette("workflow-with-fmcmc", "fmcmc")
#' 
#' vignette("user-defined-kernels", "fmcmc")
#' ```
#' 
#' @references
#' Vega Yon et al., (2019). fmcmc: A friendly MCMC framework. Journal of Open
#' Source Software, 4(39), 1427, \doi{10.21105/joss.01427}
#' 
#' @docType package
#' @name fmcmc
#' 
NULL