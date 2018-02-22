#' Markov Chain Monte Carlo
#' 
#' Metropolis-Hastings algorithm using a random walk kernel with reflecting boundaries.
#' 
#' @param fun A function. Returns the log-likelihood
#' @param initial A numeric vector. Initial values of the parameters.
#' @param nbatch Integer scalar. Number of MCMC runs.
#' @param nchains Integer scalar. Number of chains to run (in parallel).
#' @param cl A cluster object passed to \code{\link[parallel:clusterApply]{clusterApply}}.
#' @param thin Integer scalar. Passed to \code{\link[coda:mcmc]{coda::mcmc}}.
#' @param scale Either a numeric vector of length \code{length(initial)} or a
#' scalar. Step size for the transition kernel (see details).
#' @param burnin Integer scalar. Number of burn-in samples. Passed to 
#' \code{\link[coda:mcmc]{coda::mcmc}} as \code{init}.
#' @param lb Numeric vector of length \code{length(initial)}. Lower bounds
#' @param ub Numeric vector of length \code{length(initial)}. Upper bounds
#' @param useCpp Logical scalar. When \code{TRUE}, loops using a Rcpp implementation.
#' @param fixed Logical vector. If the kth position is \code{TRUE}, then that
#' value will be fixed.
#' @param ... Further arguments passed to \code{fun}.
#' 
#' @details This function implements MCMC using the Metropolis-Hastings ratio with
#' scaled standard normal propositions for each parameter. For each parameter
#' the transition function is 
#' 
#' \deqn{
#' \theta' = \theta + scale*z
#' }
#' 
#' Where \eqn{z} has standard normal distribution. The MCMC follows a block
#' sampling scheme, i.e. proposed states are either accepted or rejected
#' altogether. If \code{length(initial) > 1} and \code{length(scale) == 1},
#' the value will be recycled so that \code{length(initial) == length(scale)}.
#' 
#' Lower and upper bounds are treated using reflecting boundaries, this is, 
#' if the proposed \eqn{\theta'} is greater than the \code{ub}, then \eqn{\theta' - ub}
#' is substracted from \eqn{ub}. At the same time, if it is less than \code{lb}, then
#' \eqn{lb - \theta'} is added to \code{lb} iterating until \eqn{\theta} is within
#' \code{[lb, ub]}.
#' 
#' If \code{name(initial) == NULL}, then a names in the form of \code{par1, par2, ...}
#' will be assigned to the variables.
#' 
#' When \code{nchains > 1}, the function will run multiple chains. Furthermore,
#' if \code{cl} is not passed, \code{MCMC} will create a \code{PSOCK} cluster
#' using \code{\link[parallel:makePSOCKcluster]{makePSOCKcluster}} with
#' \code{\link[parallel:detectCores]{detectCores}}
#' clusters and try to run it using multiple cores. Internally, the function does
#' the following:
#' 
#' \preformatted{
#'   # Creating the cluster
#'   ncores <- parallel::detectCores()
#'   ncores <- ifelse(nchains < ncores, nchains, ncores)
#'   cl     <- parallel::makePSOCKcluster(ncores)
#'   
#'   # Loading the package and setting the seed using clusterRNGStream
#'   invisible(parallel::clusterEvalQ(cl, library(amcmc)))
#'   parallel::clusterSetRNGStream(cl, .Random.seed)
#' }
#' 
#' In such case, when running in parallel, objects that are
#' used within \code{fun} must be passed throught \code{...}, otherwise the cluster
#' will return with an error.
#' 
#' 
#' 
#' @return An object of class \code{\link[coda:mcmc]{mcmc}} from the \CRANpkg{coda}
#' package. The \code{mcmc} object is a matrix with one column per parameter,
#' and \code{nbatch} rows. If \code{nchains > 1}, then it returns a \code{\link[coda:mcmc]{mcmc.list}}.
#' 
#' 
#' @export
#' @examples 
#' # Univariate distributed data with multiple parameters ----------------------
#' # Parameters
#' set.seed(1231)
#' n <- 1e3
#' pars <- c(mean = 2.6, sd = 3)
#' 
#' # Generating data and writing the log likelihood function
#' D <- rnorm(n, pars[1], pars[2])
#' fun <- function(x) {
#'   x <- log(dnorm(D, x[1], x[2]))
#'   sum(x)
#' }
#' 
#' # Calling MCMC, but first, loading the coda R package for
#' # diagnostics
#' library(coda)
#' ans <- MCMC(fun, c(mu=1, sigma=1), nbatch = 2e3, scale = .1, ub = 10, lb = 0)
#' 
#' # Ploting the output
#' oldpar <- par(no.readonly = TRUE)
#' par(mfrow = c(1,2))
#' boxplot(as.matrix(ans), 
#'         main = expression("Posterior distribution of"~mu~and~sigma),
#'         names =  expression(mu, sigma), horizontal = TRUE,
#'         col  = blues9[c(4,9)],
#'         sub = bquote(mu == .(pars[1])~", and"~sigma == .(pars[2]))
#' )
#' abline(v = pars, col  = blues9[c(4,9)], lwd = 2, lty = 2)
#' 
#' plot(apply(as.matrix(ans), 1, fun), type = "l",
#'      main = "LogLikelihood",
#'      ylab = expression(L("{"~mu,sigma~"}"~"|"~D)) 
#' )
#' par(oldpar)
#' 
#' \dontrun{
#' # In this example we estimate the parameter for a dataset with ----------------
#' # With 5,000 draws from a MVN() with parameters M and S.
#' 
#' # Loading the required packages
#' library(mvtnorm)
#' library(coda)
#' 
#' # Parameters and data simulation
#' S <- cbind(c(.8, .2), c(.2, 1))
#' M <- c(0, 1)
#' 
#' set.seed(123)
#' D <- rmvnorm(5e3, mean = M, sigma = S)
#' 
#' # Function to pass to MCMC
#' fun <- function(pars) {
#'   # Putting the parameters in a sensible way
#'   m <- pars[1:2]
#'   s <- cbind( c(pars[3], pars[4]), c(pars[4], pars[5]) )
#'   
#'   # Computing the unnormalized log likelihood
#'   sum(log(dmvnorm(D, m, s)))
#' }
#' 
#' # Calling MCMC
#' ans <- MCMC(
#'   fun,
#'   initial = c(mu0=5, mu1=5, s0=5, s01=0, s2=5), 
#'   lb      = c(-10, -10, .01, -5, .01),
#'   ub      = 5,
#'   nbatch  = 1e5,
#'   thin    = 20,
#'   scale   = .01,
#'   burnin  = 5e3,
#'   useCpp  = TRUE
#' )
#' 
#' # Checking out the outcomes
#' plot(ans)
#' summary(ans)
#' 
#' # Multiple chains -----------------------------------------------------------
#' 
#' # As we want to run -fun- in multiple cores, we have to
#' # pass -D- explicitly (unless using Fork Clusters)
#' # just like specifying that we are calling a function from the
#' # -mvtnorm- package.
#'   
#' fun <- function(pars, D) {
#'   # Putting the parameters in a sensible way
#'   m <- pars[1:2]
#'   s <- cbind( c(pars[3], pars[4]), c(pars[4], pars[5]) )
#'   
#'   # Computing the unnormalized log likelihood
#'   sum(log(mvtnorm::dmvnorm(D, m, s)))
#' }
#' 
#' # Two chains
#' ans <- MCMC(
#'   fun,
#'   initial = c(mu0=5, mu1=5, s0=5, s01=0, s2=5), 
#'   nchains = 2,
#'   lb      = c(-10, -10, .01, -5, .01),
#'   ub      = 5,
#'   nbatch  = 1e5,
#'   thin    = 20,
#'   scale   = .01,
#'   burnin  = 5e3,
#'   useCpp  = TRUE,
#'   D       = D
#' )
#' 
#' summary(ans)
#' }
#' 
#' @aliases Metropolis-Hastings
MCMC <- function(
  fun,
  initial, 
  nbatch,
  nchains = 1L,
  thin    = 1L,
  scale   = rep(1, length(initial)),
  burnin  = 1e3L,
  ub      = rep(.Machine$double.xmax, length(initial)),
  lb      = rep(-.Machine$double.xmax, length(initial)),
  useCpp  = FALSE,
  cl      = NULL,
  fixed   = rep(FALSE, length(initial)),
  ...
  ) {
  
  # Filling the gap on parallel
  if ((nchains > 1L) && !length(cl)) {
    
    # Creating the cluster
    ncores <- parallel::detectCores()
    ncores <- ifelse(nchains < ncores, nchains, ncores)
    cl     <- parallel::makePSOCKcluster(ncores)
  
    # Loading the package and setting the seed using clusterRNGStream
    invisible(parallel::clusterEvalQ(cl, library(amcmc)))
    parallel::clusterSetRNGStream(cl, .Random.seed)
    
    on.exit(parallel::stopCluster(cl))
    
  }
  
  if (nchains > 1L) {

    # Running the cluster
    ans <- parallel::clusterApply(
      cl, 1:nchains, fun=
        function(i, Fun, initial, nbatch, thin, scale, burnin, ub, lb, useCpp, fixed, ...) {
          MCMC(
            fun     = Fun,
            initial = initial,
            nbatch  = nbatch,
            nchains = 1L,
            thin    = thin,
            scale   = scale,
            burnin  = burnin,
            ub      = ub,
            lb      = lb,
            useCpp  = useCpp,
            fixed   = fixed,
            ...
            )
          }, Fun = fun, nbatch=nbatch, initial = initial, thin = thin, scale = scale,
      burnin = burnin, ub = ub, lb = lb, useCpp = useCpp, fixed = fixed, ...)
    
    return(coda::mcmc.list(ans))
    
  } else {
  
    # Adding names
    cnames <- names(initial)
    if (!length(cnames))
      cnames <- paste0("par",1:length(initial))
    
    # Wrapping function. If ellipsis is there, it will wrap it
    # so that the MCMC call only uses a single argument
    passedargs <- names(list(...))
    funargs    <- methods::formalArgs(fun)
    
    # ... has extra args
    if (length(passedargs)) {
      # ... has stuff that fun doesnt
      if (any(!(passedargs %in% funargs))) {
        
        stop("Some arguments passed via -...- are not present in -fun-.")
      
      # fun has stuff that ... doesnt
      } else if (length(funargs) > 1 && any(!(funargs[-1] %in% passedargs))) {
        
        stop("-fun- requires more arguments to be passed via -...-.")
      
      # Everything OK
      } else {
        f <- function(z) {
          fun(z, ...)
        }
      }
    # ... doesnt have extra args, but funargs does!
    } else if (length(funargs) > 1) {
      
      stop("-fun- has extra arguments not passed by -...-.")
      
    # Everything OK
    } else {
      
      f <- function(z) fun(z)
      
    }
     
    # Checking boundaries
    if (length(ub) > 1 && (length(initial) != length(ub)))
      stop("Incorrect length of -ub-")
    
    if (length(lb) > 1 && (length(initial) != length(lb)))
      stop("Incorrect length of -lb-")
    
    # Repeating boundaries
    if (length(ub) == 1)
      ub <- rep(ub, length(initial))
  
    if (length(lb) == 1)
      lb <- rep(lb, length(initial))
    
    if (any(ub <= lb))
      stop("-ub- cannot be <= than -lb-.")
    
    # Repeating scale
    if (length(scale) == 1)
      scale <- rep(scale, length(initial))
    
    # Checkihg burnins
    if (burnin >= nbatch)
      stop("-burnin- (",burnin,") cannot be >= than -nbatch- (",nbatch,").")
  
    # Checking thin
    if (thin > nbatch)
      stop("-thin- (",thin,") cannot be > than -nbatch- (",nbatch,").")
    
    if (thin < 1L)
      stop("-thin- should be >= 1.")
    
    if (useCpp) {
      ans <- MCMCcpp(f, initial, nbatch, lb, ub, scale, fixed)
      dimnames(ans) <- list(1:nbatch, cnames)
      
    } else {
      R <- stats::runif(nbatch)
      ans <- matrix(ncol = length(initial), nrow = nbatch,
                    dimnames = list(1:nbatch, cnames))
      
      
      theta0 <- initial
      f0     <- f(theta0)
      for (i in 1:nbatch) {
        # Step 1. Propose
        theta1 <- normal_prop(theta0, lb, ub, scale, fixed)
        f1     <- f(theta1)
        
        # Checking f(theta1) (it must be a number, can be Inf)
        if (is.nan(f1) | is.na(f1) | is.null(f1)) 
          stop(
            "fun(par) is undefined (", f1, ")",
            "Check either -fun- or the -lb- and -ub- parameters."
          )
        
        # Step 2. Hastings ratio
        r <- min(1, exp(f1 - f0))
        
        # Updating the value
        if (R[i] < r) {
          theta0 <- theta1
          f0     <- f(theta0)
        }
        
        # Storing
        ans[i,] <- theta0
        
      }
    }
    
    
    # Thinning the data
    if (burnin) ans <- ans[-c(1:burnin), , drop=FALSE]
    if (thin)   ans <- ans[(1:nrow(ans) %% thin) == 0, , drop=FALSE]
    
    # Returning an mcmc object from the coda package
    # if the coda package hasn't been loaded, then return a warning
    if (!("package:coda" %in% search()))
      warning("The -coda- package has not been loaded.")
    
    return(
      coda::mcmc(
        ans,
        start = as.integer(rownames(ans)[1]),
        end   = as.integer(rownames(ans)[nrow(ans)]),
        thin = thin
        )
      )
  }
}

