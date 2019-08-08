#' Markov Chain Monte Carlo
#' 
#' A flexible implementation of the Metropolis-Hastings MCMC algorithm.
#' 
#' @param fun A function. Returns the log-likelihood.
#' @param initial Either a numeric matrix or vector, or an object of class [coda::mcmc]
#' or [coda::mcmc.list] (see details).
#' initial values of the parameters for each chain (See details).
#' @param nsteps Integer scalar. Length of each chain.
#' @param nchains Integer scalar. Number of chains to run (in parallel).
#' @param cl A `cluster` object passed to [parallel::clusterApply].
#' @param thin Integer scalar. Passed to [coda::mcmc].
#' @param kernel An object of class [fmcmc_kernel].
#' @param burnin Integer scalar. Length of burn-in. Passed to 
#' [coda::mcmc] as \code{start}.
#' @param multicore Logical. If `FALSE` then chains will be executed in serial.
#' @param ... Further arguments passed to \code{fun}.
#' @param conv_checker A function that receives an object of class [coda::mcmc.list],
#' and returns a logical value with `TRUE` indicating convergence.
#' 
#' @details This function implements MCMC using the Metropolis-Hastings ratio with
#' flexible transition kernels. Users can specify either one of the available
#' transition kernels or define one of their own (see [kernels]). Furthermore,
#' it allows easy parallel implementation running multiple chains in parallel. In
#' addition, we incorporate a variety of convergence diagnostics, alternatively
#' the user can specify their own (see [convergence-checker]).
#' 
#' We now give details of the various options included in the function.
#' 
#' @section Starting point:
#' 
#' By default, if `initial` is of class `mcmc`, `MCMC` will take the last `nchains`
#' points from the chain as starting point for the new sequence. If `initial` is
#' of class `mcmc.list`, the number of chains in `initial` must match the `nchains`
#' parameter. 
#' 
#' If `initial` is a vector, then it must be of length equal to the number of
#' parameters used in the model. When using multiple chains, if `initial` is not
#' an object of class `mcmc` or `mcmc.list`, then it must be a numeric matrix
#' with as many rows as chains, and as many columns as parameters in the model.
#' 
#' @section Multiple chains:
#' 
#' When \code{nchains > 1}, the function will run multiple chains. Furthermore,
#' if \code{cl} is not passed, \code{MCMC} will create a \code{PSOCK} cluster
#' using [parallel::makePSOCKcluster] with
#' [parallel::detectCores]
#' clusters and attempt to execute using multiple cores. Internally, the function does
#' the following:
#' 
#' \preformatted{
#'   # Creating the cluster
#'   ncores <- parallel::detectCores()
#'   ncores <- ifelse(nchains < ncores, nchains, ncores)
#'   cl     <- parallel::makePSOCKcluster(ncores)
#'   
#'   # Loading the package and setting the seed using clusterRNGStream
#'   invisible(parallel::clusterEvalQ(cl, library(fmcmc)))
#'   parallel::clusterSetRNGStream(cl, .Random.seed)
#' }
#' 
#' When running in parallel, objects that are
#' used within \code{fun} must be passed through \code{...}, otherwise the cluster
#' will return with an error.
#' 
#' The user controls the initial value of the parameters of the MCMC algorithm
#' using the argument `initial`. When using multiple chains, i.e., `nchains > 1`,
#' the user can specify multiple starting points, which is recommended. In such a
#' case, each row of `initial` is use as a starting point for each of the
#' chains. If `initial` is a vector and `nchains > 1`, the value is recycled, so
#' all chains start from the same point (not recommended, the function throws a
#' warning message).
#' 
#' @section Automatic stop:
#' 
#' By default, no automatic stop is implemented. If one of the functions in 
#' [convergence-checker] is used, then the MCMC is done by bulks as specified
#' by the convergence checker function, and thus the algorithm will stop if,
#' the `conv_checker` returns `TRUE`. For more information see [convergence-checker].
#' 
#' @return An object of class [coda::mcmc] from the \CRANpkg{coda}
#' package. The \code{mcmc} object is a matrix with one column per parameter,
#' and \code{nsteps} rows. If \code{nchains > 1}, then it returns a [coda::mcmc.list].
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
#' ans <- MCMC(
#'   fun, initial = c(mu=1, sigma=1), nsteps = 2e3,
#'   kernel = kernel_normal_reflective(scale = .1, ub = 10, lb = 0)
#'   )
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
#'   initial = c(mu0=5, mu1=5, s0=5, s01=0, s2=5), 
#'   fun,
#'   kernel  = kernel_normal_reflective(
#'     lb    = c(-10, -10, .01, -5, .01),
#'     ub    = 5
#'     scale = 0.01
#'   ),
#'   nsteps  = 1e5,
#'   thin    = 20,
#'   burnin  = 5e3
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
#'   initial = c(mu0=5, mu1=5, s0=5, s01=0, s2=5), 
#'   fun,
#'   nchains = 2,
#'   kernel  = kernel_normal_reflective(
#'     lb    = c(-10, -10, .01, -5, .01),
#'     ub    = 5
#'     scale = 0.01
#'   ),
#'   nsteps  = 1e5,
#'   thin    = 20,
#'   burnin  = 5e3,
#'   D       = D
#' )
#' 
#' summary(ans)
#' }
#' 
#' @aliases Metropolis-Hastings
MCMC <- function(
  initial,
  fun,
  nsteps,
  ...,
  nchains      = 1L,
  burnin       = 0L,
  thin         = 1L,
  kernel       = kernel_normal(),
  multicore    = FALSE,
  conv_checker = NULL, 
  cl           = NULL
) UseMethod("MCMC")

#' @export
#' @rdname MCMC
MCMC.mcmc <- function(
  initial,
  fun,
  nsteps,
  ...,
  nchains      = 1L,
  burnin       = 0L,
  thin         = 1L,
  kernel       = kernel_normal(),
  multicore    = FALSE,
  conv_checker = NULL, 
  cl           = NULL
) {
  
  MCMC.default(
    initial      = utils::tail(initial, nchains - 1L),
    fun          = fun,
    nsteps       = nsteps,
    ...,     
    nchains      = nchains,
    burnin       = burnin,
    thin         = thin,
    kernel       = kernel,
    multicore    = multicore,
    conv_checker = conv_checker,
    cl           = cl
  )
  
}

#' @export
#' @rdname MCMC
MCMC.mcmc.list <- function(
  initial,
  fun,
  nsteps,
  ...,
  nchains      = 1L,
  burnin       = 0L,
  thin         = 1L,
  kernel       = kernel_normal(),
  multicore    = FALSE,
  conv_checker = NULL, 
  cl           = NULL
) {
  
  if (nchains != length(initial))
    stop(
      "The parameter `nchains` must equal the number of chains passed by ",
      "`initial`.", call. = FALSE
      )
  
  MCMC.default(
    initial      = do.call(rbind, utils::tail(initial, 0)),
    fun          = fun,
    nsteps       = nsteps,
    ...,     
    nchains      = nchains,
    burnin       = burnin,
    thin         = thin,
    kernel       = kernel,
    multicore    = multicore,
    conv_checker = conv_checker,
    cl           = cl
  )
  
}



#' @export
#' @rdname MCMC
MCMC.default <- function(
  initial,
  fun,
  nsteps,
  ...,
  nchains      = 1L,
  burnin       = 0L,
  thin         = 1L,
  kernel       = kernel_normal(),
  multicore    = FALSE,
  conv_checker = NULL, 
  cl           = NULL
  ) {
  
  # # if the coda package hasn't been loaded, then return a warning
  # if (!("package:coda" %in% search()))
  #   warning("The -coda- package has not been loaded.", call. = FALSE, )

  # Checking initial argument
  initial <- check_initial(initial, nchains)
  
  if (multicore && nchains == 1L) 
    stop("When `multicore = TRUE`, `nchains` should be greater than 1.",
         call. = FALSE)
  
  if (nchains < 1L)
    stop("`nchains` must be an integer greater than 1.", call. = FALSE)
  
    # Checkihg burnins
  if (burnin >= nsteps)
    stop("-burnin- (",burnin,") cannot be >= than -nsteps- (",nsteps,").", call. = FALSE)
  
  # Checking thin
  if (thin >= nsteps)
    stop("-thin- (",thin,") cannot be > than -nsteps- (",nsteps,").", call. = FALSE)
  
  if (thin < 1L)
    stop("-thin- should be >= 1.", call. = FALSE)
  
  # Filling the gap on parallel
  if (multicore && !length(cl)) {
    
    # Creating the cluster
    ncores <- parallel::detectCores()
    ncores <- ifelse(nchains < ncores, nchains, ncores)
    cl     <- parallel::makePSOCKcluster(ncores)
    
    # Loading the package and setting the seed using clusterRNGStream
    invisible(parallel::clusterEvalQ(cl, library(fmcmc)))
    parallel::clusterSetRNGStream(cl, .Random.seed)
    
    on.exit(parallel::stopCluster(cl))
  }
  
  if (nchains > 1L) {
    
    # Preparing the call for multicore
    fmcmc_call <- as.call(
      c(
        if (multicore) 
          list(quote(parallel::clusterApply), cl=quote(cl), x = quote(1L:nchains)) 
        else 
          list(quote(lapply), X = quote(1L:nchains)),
        list(
          FUN = quote(function(
            i, fun., initial., nsteps., thin., kernel., burnin., ...) {
          
          MCMC(
            fun          = fun.,
            ...,
            initial      = initial.[i, , drop = FALSE],
            nsteps       = nsteps.,
            burnin       = burnin.,
            thin         = thin.,
            kernel       = kernel.,
            nchains      = 1L,
            multicore    = FALSE,
            cl           = NULL,
            conv_checker = NULL
            )
          
        }),
        fun.     = quote(fun),
        initial. = quote(initial),
        nsteps.  = quote(nsteps),
        burnin.  = quote(burnin),
        thin.    = quote(thin),
        kernel.  = quote(kernel),
        quote(...)
        )
      )
    )
    
    # updating names
    if (multicore)
      names(fmcmc_call)[names(fmcmc_call) == "FUN"] <- "fun"
    
    fmcmc_call <- as.call(list(quote(do.call),quote(coda::mcmc.list), fmcmc_call))
    
  } else if (!is.null(conv_checker)) {
    
    # If not multicore, still we need to make sure that we are passing some
    # variables as symbols and not as constants. As in the convergence checker
    # function we modify the current environment in order to adapt the algorithm.
    fmcmc_call              <- match.call()
    fmcmc_call$fun          <- quote(fun)
    fmcmc_call$nsteps       <- quote(nsteps)
    fmcmc_call$thin         <- quote(thin)
    fmcmc_call$burnin       <- quote(burnin)
    fmcmc_call$conv_checker <- enquote(NULL)
    
  }
  
  # If conv_checker, we run it with the conv checker and return. Notice that
  # we already adapted the code for the case in which we are runing multiple
  # chains
  if (!is.null(conv_checker)) {
    
    fmcmc_call <- call(
      "with_autostop",
      fmcmc_call,
      conv_checker = quote(conv_checker)
      )
    
    ans <- eval(fmcmc_call)
    return(ans)
    
  # If we are not using conv_checker, but still have multiple chains, then
  # we still have to run this somewhat recursively.
  } else if (nchains > 1L) {
    ans <- eval(fmcmc_call)
    return(ans)
  }
  
  # Adding names
  initial <- initial[1,,drop=TRUE]
  cnames  <- names(initial)
  
  # Wrapping function. If ellipsis is there, it will wrap it
  # so that the MCMC call only uses a single argument
  passedargs <- names(list(...))
  # print(match.call())
  funargs    <- methods::formalArgs(eval(fun))
  
  # Compiling
  # cfun <- compiler::cmpfun(fun)
  
  # ... has extra args
  if (length(passedargs)) {
    # ... has stuff that fun doesnt
    if (any(!(passedargs %in% funargs))) {
      
      stop("The following arguments passed via -...- are not present in -fun-:\n - ",
           paste(setdiff(passedargs, funargs), collapse=",\n - "),".\nThe function",
           "was expecting:\n - ", paste0(funargs, collapse=",\n - "), ".", call. = FALSE)
    
    # fun has stuff that ... doesnt
    } else if (length(funargs) > 1 && any(!(funargs[-1] %in% passedargs))) {
      
      stop("-fun- requires more arguments to be passed via -...-.", call. = FALSE)
    
    # Everything OK
    } else {
      
      f <- function(z) {
        fun(z, ...)
      }
      
    }
  # ... doesnt have extra args, but funargs does!
  } else if (length(funargs) > 1) {
    
    stop("-fun- has extra arguments not passed by -...-.", call. = FALSE)
    
  # Everything OK
  } else {
    
    f <- function(z) fun(z)
    
  }
  
  # MCMC algorithm -----------------------------------------------------------
    
  theta0 <- initial
  theta1 <- theta0
  f0     <- f(theta0)
  f1     <- f(theta1)
  
  # The updates can be done jointly or sequentially
  klogratio <- kernel$logratio(environment())
  if (length(klogratio) > 1L) {
    joint_rate <- FALSE
    R   <- matrix(stats::runif(nsteps * length(initial)), nrow = nsteps)
    R[] <- log(R)
  } else {
    joint_rate <- TRUE
    R <- matrix(log(stats::runif(nsteps)), nrow = nsteps)
  }
  
  ans <- matrix(ncol = length(initial), nrow = nsteps,
                dimnames = list(1:nsteps, cnames))
  
  for (i in 1L:nsteps) {
    # Step 1. Propose
    theta1[] <- kernel$proposal(environment())
    f1       <- f(theta1)
    
    # Checking f(theta1) (it must be a number, can be Inf)
    if (is.nan(f1) | is.na(f1) | is.null(f1)) 
      stop(
        "fun(par) is undefined (", f1, "). ",
        "Check either -fun- or the -lb- and -ub- parameters. ",
        "This error ocurred during step i = ", i, " and proposal parameters ",
        "theta1 = \n", sprintf(" %10s: %.4f\n", names(theta1), theta1), "\n",
        call. = FALSE
      )
    
    # Step 2. Hastings ratio
    klogratio <- kernel$logratio(environment())
    if (joint_rate) {
      
      if (R[i] < klogratio) {
        theta0 <- theta1
        f0     <- f1
      }
      
    } else {
      
      klogratio <- (R[i, ] < klogratio)
      theta0[klogratio] <- theta1[klogratio]
      
    }
    
    
    # Step 3. Saving the state
    ans[i,] <- theta0
    
  }
  
  # Thinning the data
  if (burnin) ans <- ans[-c(1:burnin), , drop = FALSE]
  if (thin)   ans <- ans[(1:nrow(ans) %% thin) == 0, , drop = FALSE]
  
  # Returning an mcmc object from the coda package
  return(
    coda::mcmc(
      ans,
      start = as.integer(rownames(ans)[1]),
      end   = as.integer(rownames(ans)[nrow(ans)]),
      thin  = thin
      )
    )
  
}
