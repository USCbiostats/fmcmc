#' Markov Chain Monte Carlo
#' 
#' A flexible implementation of the Metropolis-Hastings MCMC algorithm, users
#' can utilize arbitrary transition kernels as well as set-up an automatic 
#' stop criterion using a convergence check test.
#' 
#' @param fun A function. Returns the log-likelihood.
#' @param initial Either a numeric matrix or vector, or an object of class [coda::mcmc]
#' or [coda::mcmc.list] (see details).
#' initial values of the parameters for each chain (See details).
#' @param seed If not null, passed to [set.seed].
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
#' and returns a logical value with `TRUE` indicating convergence. See the
#' "Automatic stop" section and the [convergence-checker] manual.
#' @param progress Logical scalar. When set to `TRUE` shows a progress bar. A new
#' bar will be show every time that the convergence checker is called.
#' @param chain_id Integer scalar (internal use only). This is an argument
#' passed to the kernel function and it allows it identify in which of the
#' chains the process is taking place. This could be relevant for some kernels
#' (see [kernel_new()]).
#' 
#' @details This function implements Markov Chain Monte Carlo (MCMC) using the
#' Metropolis-Hastings ratio with
#' flexible transition kernels. Users can specify either one of the available
#' transition kernels or define one of their own (see [kernels]). Furthermore,
#' it allows easy parallel implementation running multiple chains in parallel.
#' The function also allows using convergence diagnostics tests to set-up a
#' criterion for automatically stopping the algorithm  (see [convergence-checker]).
#' 
#' The canonical form of the Metropolis Hastings algorithm consists on accepting
#' a move from state \eqn{x} to state \eqn{y} based on the Hastings ratio \eqn{r(x,y)}:
#' 
#' \deqn{%
#' r(x,y) = \frac{h(y)q(y,x)}{h(x)q(x,y)},%
#' }{%
#' r(x,y) = [h(y)q(y,x)]/[h(x)q(x,y)],%
#' }
#' 
#' where \eqn{h} is the unnormalized density of the specified distribution (
#' the posterior probability), and \eqn{q} has the conditional probability of
#' moving from state \eqn{x} to \eqn{y} (the proposal density). The move
#' \eqn{x \to y}{x->y} is then accepted with probability
#' 
#' \deqn{%
#' \alpha(x,y) = \min\left(1, r(x,y)\right)%
#' }{%
#' alpha(x, y) = min {1, r(x, y)}%
#' }
#' 
#' Observe that, in the case that \eqn{q()} is symmetric, meaning \eqn{q(x, y) = q(y, x)},
#' the Hastings ration reduces to \eqn{h(y)/h(x)}. Starting version 0.5-0, the value 
#' of the log unnormalized density and the proposed states `y` can be accessed using
#' the functions [get_logpost()] and [get_draws()].
#' 
#' We now give details of the
#' various options included in the function.
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
#' using \code{\link[parallel:makeCluster]{makePSOCKcluster}} with
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
#' @return 
#' `MCMC` returns an object of class [coda::mcmc] from the \CRANpkg{coda}
#' package. The \code{mcmc} object is a matrix with one column per parameter,
#' and \code{nsteps} rows. If \code{nchains > 1}, then it returns a [coda::mcmc.list].
#'    
#' @references 
#' Brooks, S., Gelman, A., Jones, G. L., & Meng, X. L. (2011). Handbook of
#' Markov Chain Monte Carlo. Handbook of Markov Chain Monte Carlo.
#' 
#' Vega Yon, G., & Marjoram, P. (2019). fmcmc: A friendly MCMC framework.
#' Journal of Open Source Software, 4(39), 1427. \doi{10.21105/joss.01427}
#' 
#' @export
#'  
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
#' 
#' # In this example we estimate the parameter for a dataset with ----------------
#' # With 5,000 draws from a MVN() with parameters M and S.
#' \donttest{
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
#'     ub    = 5,
#'     scale = 0.01
#'   ),
#'   nsteps  = 1e3,
#'   thin    = 10,
#'   burnin  = 5e2
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
#'     ub    = 5,
#'     scale = 0.01
#'   ),
#'   nsteps  = 1e3,
#'   thin    = 10,
#'   burnin  = 5e2,
#'   D       = D
#' )
#' 
#' summary(ans)
#' }
#' 
#' # Example using a user-defined cl object -------------------------------------
#' 
#' # A silly function that gets the first two parameters of the
#' # vector. Using the default multicore example will cause and error
#' get_m <- function(pars) {
#'   pars[1:2]
#' }
#' 
#' fun <- function(pars, D) {
#'   # Putting the parameters in a sensible way
#'   m <- get_m(pars)
#'   s <- cbind( c(pars[3], pars[4]), c(pars[4], pars[5]) )
#'   
#'   # Computing the unnormalized log likelihood
#'   sum(log(mvtnorm::dmvnorm(D, m, s)))
#' }
#' 
#' if (FALSE) {
#' 
#'   # Thi will fail with the error
#'   # Error in checkForRemoteErrors(val) :
#'   #   4 nodes produced errors; first error: could not find function "get_m"
#'   ans <- MCMC(
#'     initial = c(mu0=5, mu1=5, s0=5, s01=0, s2=5), 
#'     fun,
#'     kernel  = kernel_normal_reflective(
#'       lb    = c(-10, -10, .01, -5, .01),
#'       ub    = 5,
#'       scale = 0.01
#'     ),
#'     nsteps  = 1e3,
#'     thin    = 10,
#'     burnin  = 5e2,
#'     D       = D,
#'     # Multiple chains using parallel computing
#'     multicore = TRUE,
#'     nchains = 4
#'   )
#'   
#'   summary(ans)
#' 
#' }
#' 
#' # To solve the error, we need to build the cluster object
#' library(parallel)
#' cl <- makePSOCKcluster(4)
#' 
#' # Export the function `get_m` to the cluster. The function `fun` and data `D`
#' # are automatically exported by the `MCMC` function.
#' clusterExport(cl, "get_m")
#' 
#' # Run the MCMC
#' ans <- MCMC(
#'   initial = c(mu0=5, mu1=5, s0=5, s01=0, s2=5), 
#'   fun,
#'   kernel  = kernel_normal_reflective(
#'     lb    = c(-10, -10, .01, -5, .01),
#'     ub    = 5,
#'     scale = 0.01
#'   ),
#'   nsteps  = 1e3,
#'   thin    = 10,
#'   burnin  = 5e2,
#'   D       = D,
#'   # Multiple chains using parallel computing
#'   multicore = TRUE,
#'   nchains = 4,
#'   # Use the cluster object
#'   cl = cl
#' )
#' 
#' summary(ans)
#' 
#' 
#' @aliases Metropolis-Hastings
#' @seealso [get_logpost()], [get_logpost()] ([mcmc-output]) for post execution of `MCMC`, and
#' [ith_step()] for accessing objects within an `MCMC` call.
MCMC <- function(
  initial,
  fun,
  nsteps,
  ...,
  seed         = NULL,
  nchains      = 1L,
  burnin       = 0L,
  thin         = 1L,
  kernel       = kernel_normal(),
  multicore    = FALSE,
  conv_checker = NULL, 
  cl           = NULL,
  progress     = interactive() && !multicore,
  chain_id     = 1L
) UseMethod("MCMC")

#' @export
# @rdname MCMC
MCMC.mcmc <- function(
  initial,
  fun,
  nsteps,
  ...,
  seed         = NULL,
  nchains      = 1L,
  burnin       = 0L,
  thin         = 1L,
  kernel       = kernel_normal(),
  multicore    = FALSE,
  conv_checker = NULL, 
  cl           = NULL,
  progress     = interactive() && !multicore,
  chain_id     = 1L
) {
  
  MCMC.default(
    initial      = utils::tail(initial, nchains - 1L),
    fun          = fun,
    nsteps       = nsteps,
    ...,     
    seed         = seed,
    nchains      = nchains,
    burnin       = burnin,
    thin         = thin,
    kernel       = kernel,
    multicore    = multicore,
    conv_checker = conv_checker,
    cl           = cl,
    progress     = progress,
    chain_id     = chain_id
  )
  
}

#' @export
# @rdname MCMC
MCMC.mcmc.list <- function(
  initial,
  fun,
  nsteps,
  ...,
  seed         = NULL,
  nchains      = 1L,
  burnin       = 0L,
  thin         = 1L,
  kernel       = kernel_normal(),
  multicore    = FALSE,
  conv_checker = NULL, 
  cl           = NULL,
  progress     = interactive() && !multicore,
  chain_id     = 1L
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
    seed         = seed,
    nchains      = nchains,
    burnin       = burnin,
    thin         = thin,
    kernel       = kernel,
    multicore    = multicore,
    conv_checker = conv_checker,
    cl           = cl,
    progress     = progress,
    chain_id     = chain_id
  )
  
}

#' @export
# @rdname MCMC
MCMC.default <- function(
  initial,
  fun,
  nsteps,
  ...,
  seed         = NULL,
  nchains      = 1L,
  burnin       = 0L,
  thin         = 1L,
  kernel       = kernel_normal(),
  multicore    = FALSE,
  conv_checker = NULL, 
  cl           = NULL,
  progress     = interactive() && !multicore,
  chain_id     = 1L
) {
  
  # Initializing recording of variables
  MCMC_OUTPUT$clear(nchains = nchains)
  MCMC_OUTPUT$set_ptr(1L)
  
  MCMC_init(...)
  
  MCMC_CALL <- if (!is.null(conv_checker))
    MCMC_with_conv_checker
  else
    MCMC_without_conv_checker
  
  # Checking the seed
  if (!is.null(seed))
    set.seed(seed)
  
  ans <- MCMC_CALL(
    initial      = initial,
    fun          = fun,
    nsteps       = nsteps,
    nchains      = nchains,
    burnin       = burnin,
    thin         = thin,
    kernel       = kernel,
    multicore    = multicore,
    conv_checker = conv_checker,
    cl           = cl,
    progress     = progress,
    chain_id     = chain_id,
    ...
  )
  
  # Stopping the timer
  MCMC_finalize()
  
  return(ans)

}

#' @export
#' @details The function `MCMC_without_conv_checker` is for internal use
#' only.
#' @rdname MCMC
MCMC_without_conv_checker <- function(
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
  cl           = NULL,
  progress     = interactive() && !multicore,
  chain_id     = 1L
  ) {
  
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
  
  # If we are running in multiple chains, and this is not a kernel_list
  # object, then we need to replicate it. This will modify the original kernel
  # object creating a new environment with `nchains` copies of the original
  # kernel.
  if (nchains > 1L && !is_kernel_list(kernel))
    rep_kernel(kernel, nchains = nchains)
  else if (nchains == 1L && is_kernel_list(kernel))
    stop(
      "The passed kernel is for MCMC with more than one chain. ",
      "Right now, -kernel- is of length ", length(kernel),
      call. = FALSE
      )
  
  # Filling the gap on parallel
  if (multicore && !length(cl)) {
    
    # Creating the cluster
    ncores <- parallel::detectCores()
    ncores <- ifelse(nchains < ncores, nchains, ncores)
    cl     <- parallel::makePSOCKcluster(ncores)
    
    parallel::clusterSetRNGStream(cl, .Random.seed)
    
    on.exit(parallel::stopCluster(cl))
    
  }
  
  # We need to actively export these so we can capture it later
  if (multicore) {
    
    # Loading the package and setting the seed using clusterRNGStream
    invisible(parallel::clusterEvalQ(cl, {
      if (!"fmcmc" %in% .packages())
      library(fmcmc)
      }))
    
    parallel::clusterExport(cl, "kernel", envir = environment())
    parallel::clusterEvalQ(cl, {
      .FMCMC_PLL_KERNEL  <- new.env()
      .FMCMC_MCMC_OUTPUT <- new.env()
      .FMCMC_MCMC_OUTPUT_usr <- new.env()
      })
    
    # Exporting data
    invisible({
      parallel::parLapply(
        cl, MCMC_OUTPUT$data., function(i) assign("ptr__", i, envir = .FMCMC_MCMC_OUTPUT)
        )
      
      parallel::parLapply(
        cl, MCMC_OUTPUT$data_usr., function(i) assign("ptr_usr__", i, envir = .FMCMC_MCMC_OUTPUT_usr)
      )
    })
    
    # Setting up data
    invisible(parallel::clusterEvalQ(cl, {
      
      MCMC_OUTPUT$clear(
        1L,
        env_data     = get("ptr__", envir = .FMCMC_MCMC_OUTPUT),
        env_data_usr = get("ptr_usr__", envir = .FMCMC_MCMC_OUTPUT_usr)
        )
      
      MCMC_OUTPUT$set_ptr(1L)
    }))
    
  }
    
  if (nchains > 1L && multicore) {
    
    # Recursively calling the function
    ans <- parallel::parLapply(
      cl = cl, X = seq_len(nchains),
      fun = function(
        i, initial, fun., nsteps, nchains, burnin, thin, 
        progress
        ) {
        res. <- fmcmc::MCMC_without_conv_checker(
          initial      = initial[i, , drop=FALSE],
          fun          = fun.,
          nsteps       = nsteps,
          nchains      = 1L,
          burnin       = burnin,
          thin         = thin,
          kernel       = kernel[[i]],
          multicore    = FALSE,
          conv_checker = NULL,
          cl           = NULL,
          progress     = progress,
          chain_id     = i,
          ...
        )
        
        # Making sure we are returning the right kernel
        assign("kernel", kernel[[i]], envir = .FMCMC_PLL_KERNEL)
        
        res.
      },
      initial      = initial,
      fun.         = fun,
      nsteps       = nsteps,
      nchains      = nchains,
      burnin       = burnin,
      thin         = thin,
      progress     = progress
    )
    
    # Updating the kernel
    kernel_list <- parallel::clusterEvalQ(cl, get("kernel", envir = .FMCMC_PLL_KERNEL, inherits = FALSE))
    update_kernel(kernel, do.call(c, kernel_list))
    
    # Updating run
    MCMC_OUTPUT$data.     <- parallel::clusterEvalQ(cl, get("ptr", envir = MCMC_OUTPUT, inherits = FALSE))
    MCMC_OUTPUT$data_usr. <- parallel::clusterEvalQ(cl, get("ptr_usr", envir = MCMC_OUTPUT, inherits = FALSE))
    
    # So things are pointing to the right place
    MCMC_OUTPUT$set_ptr(1L)
    
    # Appending the chains
    return(coda::as.mcmc.list(ans))
    
  } else if (nchains > 1L && !multicore) {
    
    # Recursively calling the function
    ans <- vector("list", nchains)
    for (i in seq_len(nchains)) {
      
      # Where to record logpost and other stuff
      MCMC_OUTPUT$set_ptr(i)
      
      ans[[i]] <- MCMC_without_conv_checker(
        initial      = initial[i,,drop=FALSE],
        fun          = fun,
        nsteps       = nsteps,
        nchains      = 1L,
        burnin       = burnin,
        thin         = thin,
        kernel       = kernel[[i]],
        multicore    = multicore,
        conv_checker = NULL,
        cl           = NULL,
        progress     = progress,
        chain_id     = i,
        ...
        )
      
    }
    
    # Appending the chains
    return(coda::as.mcmc.list(ans))
    
  }
  
  # Last checks before we lunch the MCMC ---------------------------------------
  
  # Adding names
  initial <- initial[1,,drop=TRUE]
  cnames  <- names(initial)
  
  # Wrapping function. If ellipsis is there, it will wrap it
  # so that the MCMC call only uses a single argument
  passedargs <- names(list(...))
  funargs    <- methods::formalArgs(eval(fun))

  # ... has extra args
  if (length(passedargs)) {
    # ... has stuff that fun doesn't
    if (any(!(passedargs %in% funargs))) {
      
      stop("The following arguments passed via -...- are not present in -fun-:\n - ",
           paste(setdiff(passedargs, funargs), collapse=",\n - "),".\nThe function",
           "was expecting:\n - ", paste0(funargs, collapse=",\n - "), ".", call. = FALSE)
    
    # fun has stuff that ... doesn't
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
  
  MCMC_OUTPUT$loop_envir <- environment()
  on.exit(MCMC_OUTPUT$loop_envir <- NULL, add = TRUE)
  
  # The updates can be done jointly or sequentially
  R <- matrix(log(stats::runif(nsteps)), nrow = nsteps)
  
  ans <- matrix(ncol = length(initial), nrow = nsteps,
                dimnames = list(1:nsteps, cnames))
  logpost <- vector("numeric", nsteps)
  draws   <- matrix(
    NA_real_, nrow = nsteps, ncol = length(initial),
    dimnames = list(1:nsteps, cnames)
    )
  
  # We start assuming the first call is accepted
  i           <- 1L # Need that in case we have it
  theta0      <- initial
  theta1      <- theta0
  draws[1L,]  <- initial
  ans[1L, ]   <- theta0
  logpost[1L] <- f(theta0)
  f0          <- logpost[1L]
  
  if (progress)
    progress_bar <- new_progress_bar(nsteps)
  
  # We start from the second, since we already completed the first
  for (i in 2L:nsteps) {

    # Step 1. Propose
    draws[i,] <- kernel$proposal(environment())
    theta1[]       <- draws[i,]
    logpost[i]     <- f(theta1)
    f1             <- logpost[i]

    # Checking f(theta1) (it must be a number, can be Inf)
    if (is.nan(f1) || is.na(f1) || is.null(f1)) 
      stop(
        "fun(par) is undefined (", f1, "). ",
        "Check either -fun- or the -lb- and -ub- parameters. ",
        "This error ocurred during step i = ", i, " and proposal parameters ",
        "theta1 = \n", sprintf(" %10s: %.4f\n", names(theta1), theta1), "\n",
        call. = FALSE
      )
    
    # Step 2. Hastings ratio
    klogratio <- kernel$logratio(environment())

    if (R[i] < klogratio) {
      
      theta0 <- theta1
      f0     <- logpost[i]
      
    }
      
    # Step 3. Saving the state
    ans[i,] <- theta0
    
    if (progress)
      progress_bar(i)
    
  }
  
  # Thinning the data
  if (burnin) {
    
    to_burn <- -c(1L:burnin)
    
    ans     <- ans[to_burn, , drop = FALSE]
    logpost <- logpost[to_burn]
    draws   <- draws[to_burn, , drop=FALSE]
    
    if (length(MCMC_OUTPUT$ptr_usr))
      MCMC_OUTPUT$ptr_usr[["userdata"]] <- MCMC_OUTPUT$ptr_usr[["userdata"]][
        to_burn, , drop = FALSE
      ]
    
  }
  
  if (thin) {
    
    to_include <- which((1L:nrow(ans) %% thin) == 0L)
    logpost    <- logpost[to_include]
    ans        <- ans[to_include, , drop = FALSE]
    draws      <- draws[to_include, , drop=FALSE]
    
    if (length(MCMC_OUTPUT$ptr_usr))
      MCMC_OUTPUT$ptr_usr[["userdata"]] <- MCMC_OUTPUT$ptr_usr[["userdata"]][
        to_include, , drop = FALSE
      ]
    
  }
  
  names(logpost)  <- rownames(ans)
  rownames(draws) <- names(logpost)
  
  if (length(MCMC_OUTPUT$ptr_usr))
    rownames(MCMC_OUTPUT$ptr_usr[["userdata"]]) <- names(logpost)
  
  # Storing the logpost and draws
  MCMC_OUTPUT$ptr[["logpost"]] <- logpost
  MCMC_OUTPUT$ptr[["draws"]]   <- draws
  
  # Cleaning the space
  on.exit(rm(list = ls(envir = environment())), add = TRUE)
  
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

#' @noRd
MCMC_with_conv_checker <- function(
  initial,
  fun,
  nsteps,
  ...,
  nchains      ,
  burnin       ,
  thin         ,
  kernel       ,
  multicore    ,
  conv_checker , 
  cl           ,
  progress     ,
  chain_id
){
  
  if (is.null(conv_checker))
    stop("The convergence checker for this call cannot be null.", call. = FALSE)


  # Getting the parent environment
  freq   <- attr(conv_checker, "freq")
  if (is.null(freq)) {
    
    freq <- floor(nsteps/2)
    
    warning(
      "The -conv_checker- function has no freq attribute. ",
      "Default value set to be ", freq, call. = FALSE, immediate. = TRUE
      )
    
  }

  # Correcting the freq
  if (freq * 2L > nsteps) 
    freq <- 0L
  
  # Calculating lengths. The bulk vector sets what will be
  # nsteps in each call. This excludes burnin
  bulks <- if (freq > 0L)
    rep(freq, (nsteps - burnin) %/% freq)
  else
    nsteps
  
  if (freq > 0 && (nsteps - burnin) %% freq)
    bulks <- c(bulks, (nsteps - burnin) - sum(bulks))
  
  # We need to add the burnin to the first
  bulks[1] <- bulks[1] + burnin
  
  # Do while no convergence
  converged   <- FALSE
  i           <- 0L
  ans         <- NULL
  logpost     <- NULL
  free_params <- NULL
  
  # Cleaning the convergence environment
  convergence_data_flush()
  
  for (i in seq_along(bulks)) {
    
    # Updating the nsteps argument
    nsteps <- bulks[i]
    
    if (i > 1) {
      
      burnin  <- 0L
      initial <- ans[coda::niter(ans),]
      if (coda::is.mcmc.list(ans))
        initial <- do.call(rbind, initial)
      
      # Retrieving information from the previous run (we need to combine it)
      prev_logpost       <- get_logpost()
      prev_draws         <- get_draws()
      
      # Userdata not always exists, so we only store if there's something
      # to store!
      prev_userdata <- tryCatch(get_userdata(), error = function(e) e)
      if (inherits(prev_userdata, "error"))
        prev_userdata <- NULL
    
    }
    
    # Running the MCMC and adding it to the tail
    tmp <- MCMC_without_conv_checker(
      initial      = initial,
      fun          = fun,
      nsteps       = nsteps,
      nchains      = nchains,
      burnin       = burnin,
      thin         = thin,
      kernel       = kernel,
      multicore    = multicore,
      conv_checker = NULL,
      cl           = cl,
      progress     = progress,
      chain_id     = chain_id,
      ...
    )
    
    # A friendly call to the garbage collector
    if (multicore && !is.null(cl))
      parallel::clusterEvalQ(cl, gc())
    
    # Appending retults
    ans <- append_chains(ans, tmp)
    
    # Checking the set of free parameters
    if (is.null(free_params)) {
      
      free_params <- if (inherits(kernel, "fmcmc_kernel_list"))
        kernel[[1]]$fixed
      else
        kernel$fixed
      
      if (is.null(free_params))
        free_params <- seq_along(initial)
      else
        free_params <- which(!free_params)
      
    }
    
    # Resetting the convergence message
    convergence_msg_set()
    
    # Checking convergence on the gree parameters only
    converged <- conv_checker(ans[, free_params, drop = FALSE])
    
    # Anything to say?
    msg <- convergence_msg_get()
    
    # Wrapping MCMC_DATA
    rnames <- if (nchains > 1L)
      rownames(ans[[1L]])
    else
      rownames(ans)
    
    if (i > 1L) {
      MCMC_OUTPUT$append_("logpost", prev_logpost, "data.", "c", rnames)
      MCMC_OUTPUT$append_("draws", prev_draws, "data.", "rbind", rnames)
      
      # Userdata not always exists, so we only append if there's something
      # to append!
      if (!is.null(prev_userdata))
        MCMC_OUTPUT$append_("userdata", prev_userdata, "data_user.", "rbind", rnames)
    }
    
    
    if (converged) {
      
      message(
        "Convergence has been reached with ", sum(bulks[1:i]), " steps. ",
        ifelse(is.na(msg), "", paste0(msg, " ")), "(",
        coda::niter(ans), " final count of samples)."
      )
      break
      
    } else {
      
      message(
        "No convergence yet (steps count: ", sum(bulks[1:i]), "). ",
        ifelse(is.na(msg), "", paste0(msg, " ")),
        "Trying with the next bulk."
      )
    }
    
  }
  
  # Did it converged?
  if (!is.null(conv_checker) && (i == length(bulks) & !converged))
    message("No convergence reached after ", sum(bulks[1:i]), " steps (",
            coda::niter(ans), " final count of samples).")
  
  # Returning
  return(ans)
  
  
}
