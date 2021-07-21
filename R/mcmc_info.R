#' Information about the last `MCMC` call
#' 
#' This environment holds a copy of the last call to [MCMC], including the start
#' and end time (to compute total elapsed time) of the call. Since the resulting
#' object of `MCMC` is an object of class [coda::mcmc], this is a way to capture
#' more information in case the user needs it.
#' 
#' 
#' @name mcmc-output
NULL

#' @return The `MCMC_OUTPUT` object is an environment of class, 
#' `c("fmcmc_output", "environment")` that has the following structure:
#' 
#' - `time_start`, `time_end` Objects of class [proc.time]. Mark the start and
#'   end of the `MCMC` call.
#' 
#' - `data.` A list of environments of length `get_nchains()`. Each environment
#'   will hold information about the particular chain. By default, each environment
#'   holds the elements `logpost` (named numeric vector) and `draws` (numeric matrix).
#'   
#'   The `draws` matrix contains the draws from the proposal kernel function. Both 
#'   `logpost` and `draws` have indices that match those of the chain. (see details).
#'   
#'   `data.` can also be accessed by the user to store information if needed.
#'   
#' - `ptr` An environment. This is used as a pointer that is defined at the beginning
#'   of the MCMC process. The environment will be pointing to the current chain, thus,
#'   if `MCMC` is running chain 2 of 4, `ptr = data.[[2]]`.
#'   
#' - `i` Integer. Index of the current chain, so if `MCMC` is running chain 3 of 4,
#'   then `i = 3` (and `ptr = data.[[3]]`, get it?).
#'   
#' - `nchains` Integer. The number of chains specified in `MCMC`.
#' 
#' - `...` further arguments passed to `MCMC`, e.g., `initial`, `fun`, `nsteps`,
#'   `kernel`, `thin`, etc.
#'   
#' It also contains the following **helper functions**:
#' 
#' - `c_(x, val)` Combine elements. It will access the current value of `x` in
#'   `ptr`, and will combine it with `val` using [c()]. 
#' 
#' - `rbind_(x, val)` Row-combine elements. It will access the current value of `x` in
#'   `ptr`, and will combine it with `val` using [rbind()].
#' 
#' - `cbind_(x, val)` Column-combine elements. It will access the current value of `x` in
#'   `ptr`, and will combine it with `val` using [cbind()].
#' 
#' 
#' `get_*` returns the corresponding variable passed to the last call
#' of [MCMC].
#' @noRd
NULL

#' @rdname MCMC
#' @return While the main output of `MCMC` is the `mcmc`(`.list`) object, other information
#' and intermediate outputs of the process are stored in `MCMC_OUTPUT`. For information
#' about how to retrieve/set data, see [mcmc-output].
#' @format NULL
#' @export
MCMC_OUTPUT <- structure(
  list2env(list(
    time_start = NULL,
    time_end   = NULL,
    data.      = list(),
    data_usr.  = list(),
    ptr        = NULL,
    ptr_usr    = NULL,
    i          = NA_integer_,
    nchains    = 0L,
    loop_envir = NULL   
  ), envir = new.env()),
  class = c("fmcmc_output", "environment"))

#' Clears the MCMC_OUTPUT environment and sets the number of chains
#' @param nchains an integer
#' @param env an optional environment (this is useful in the context of parallel chains.)
#' @noRd 
MCMC_OUTPUT$clear <- function(nchains, env_data, env_data_usr) {
  
  lapply(MCMC_OUTPUT$data., function(d) rm(list = ls(all.names = TRUE, envir = d), envir = d))
  
  if (missing(env_data)) {
    
    # Making some room
    MCMC_OUTPUT$data. <- replicate(
      nchains, list2env({
        list(logpost = numeric(), draws = NULL)
      }
      ), simplify = FALSE)
    
    # Making more space for the user
    MCMC_OUTPUT$data_usr. <- replicate(
      nchains, list2env({list()}), simplify = FALSE
      )
    
  } else {
    
    if (nchains != 1L)
      stop("When clearing with an environment, nchains should be 1.")
    
    MCMC_OUTPUT$data.     <- list(env_data)
    MCMC_OUTPUT$data_usr. <- list(env_data_usr)
    
  }
    
  MCMC_OUTPUT$nchains <- nchains
  MCMC_OUTPUT$info    <- new.env()
  invisible()
  
}

#' This function sets the current pointer (useful when running in
#' serial fashion.)
#' @noRd 
MCMC_OUTPUT$set_ptr <- function(i) {
  
  if (i > MCMC_OUTPUT$nchains)
    stop("fmcmc_output pointer out of range.", call. = FALSE)
  
  MCMC_OUTPUT$ptr     <- MCMC_OUTPUT$data.[[i]]
  MCMC_OUTPUT$ptr_usr <- MCMC_OUTPUT$data_usr.[[i]]
  
  MCMC_OUTPUT$i   <- i
  
  invisible()
  
}

# set_user_data <- function()
# get_user_data

#' Append data
#' @param x Name of the object.
#' @param val Value to append.
#' @param what Either "data." or "data_usr.".
#' @param call. Function for combining.
#' @noRd
MCMC_OUTPUT$append_ <- function(x, val, what, call., names.) {
  
  if (get_nchains() > 1L) {
    
    env <- MCMC_OUTPUT[[what]][[i]]
    
    for (i in 1L:get_nchains()) {
      
      assign(x = x[[i]], value = do.call(call., list(val, env[[x]])), envir = env)
      
      if (call. == "c")
        names(env[[x]]) <- names.
      else if (call. == "rbind")
        rownames(env[[x]]) <- names.
      
    }
  } else {
    
    assign(
      x = x, value = do.call(call., list(
        val,
        MCMC_OUTPUT[[what]][[1L]][[x]]
        )),
      envir = MCMC_OUTPUT[[what]][[1L]]
      )
    
    if (call. == "c")
      names(MCMC_OUTPUT[[what]][[1L]][[x]]) <- names.
    else if (call. == "rbind")
      rownames(MCMC_OUTPUT[[what]][[1L]][[x]]) <- names.
    
  }

}


#' @export
print.fmcmc_output <- function(x, ...) {
  
  if (get_nchains() == 0L)
    cat("-MCMC- has not been called yet. Nothing to show.\n")
  else {
    cat("Last call to MCMC holds the following elements:\n")
    
    for (i in 1:get_nchains()) {
      cat("Chain N: ", i, "\n", sep = "")
      print(utils::ls.str(x[["data."]][[i]]))
    }
    
    if (length(x$ptr_usr)) {
      
      cat("\nIncluding the following userdata (use -get_userdata()- to access it):\n")
      tmp <- get_userdata()
      for (i in 1:get_nchains()) {
        cat("Chain N ", i, "\n", sep = "")
        print(head(tmp[[i]], n = 6L))
        if (nrow(tmp[[i]]) > 6L)
          cat("...", nrow(tmp[[i]]) - 6L, " more...\n")
      }
      
    }
  }
  
  invisible(x)
  
}

#' @export
print.fmcmc_last_mcmc <- function(x, ...) {
  
  .Deprecated("MCMC_OUTPUT")
  print(MCMC_OUTPUT)
  invisible(x)
  
}

MCMC_init <- function(...) {
  
  # Getting the caller environment
  MCMC_OUTPUT$time_start <- proc.time()
  env <- parent.frame()
  
  # Initializing the variables
  for (n in names(env))
    if (n != "...") 
      assign(n, get(n, envir = env), envir = MCMC_OUTPUT$info)
    
  # Assigning dots
  for (i in seq_len(...length())) {
    assign(...names()[i], ...elt(i), envir = MCMC_OUTPUT$info)
  }
    
  invisible(NULL)
  
}

MCMC_finalize <- function() {
  MCMC_OUTPUT$time_end <- proc.time()
}


#' @export
#' @rdname mcmc-output
#' @param x Character scalar. Name of an argument to retrieve. If `x` was not
#' passed to the last call, the function returns with an error.
get_ <- function(x) {
  
  if (MCMC_OUTPUT$nchains == 0L)
    stop("-MCMC- has not been called yet.", call. = FALSE)
  
  # Checking if it exists outside
  if (exists(x, envir = MCMC_OUTPUT, inherits = FALSE)) {
    return(MCMC_OUTPUT[[x]])
  }
  
  if (exists(x, envir = MCMC_OUTPUT$info, inherits = FALSE)) {
    return(MCMC_OUTPUT$info[[x]])
  }
  
  # Otherwise, it should be part of the run
  res <- tryCatch(lapply(MCMC_OUTPUT$data., function(env) {
    
    if (!exists(x, envir = env, inherits = FALSE))
      stop()
    
    env[[x]]
    
  }), error = function(e) e)
  
  if (!inherits(res, "error")) {
    
    # Simplifying the output
    if (get_nchains() == 1L)
      res <- res[[1L]]
    
    return(res)
    
  }
  
  
  res <- lapply(MCMC_OUTPUT$data_usr., function(env) {
    
    if (!exists(x, envir = env, inherits = FALSE))
      stop(
        "The object -", x, "- was not found in the last MCMC call.",
        call. = FALSE
      )
    
    env[[x]]
    
  })
  
  # Simplifying the output
  if (get_nchains() == 1L)
    res <- res[[1L]]
  
  return(res)
  
}

#' @export
#' @rdname mcmc-output
#' @details The function `get_logpost` returns the `logposterior` value at each
#' iteration. The values correspond to a named numeric vector. If `nchains > 1`
#' then it will return a list of length `nchains` with the corresponding logpost
#' values for each chain.
#' @examples 
#' # Getting the logpost -------------------------------------------------------
#' set.seed(23133)
#' x <- rnorm(200)
#' y <- x*2 + rnorm(200)
#' f <- function(p) {
#'   sum(dnorm(y - x*p, log = TRUE))
#' }
#' 
#' ans <- MCMC(fun = f, initial = c(0), nsteps=2000)
#' plot(
#'   # Plotting the logpost from the last run
#'   get_logpost(), 
#'   
#'   # Getting the number of chains
#'   main = paste0("nchains: ", get_nchains()),
#'   
#'   # And the elapsed time
#'   sub  = sprintf("Run time: %.4f(s)", get_elapsed()[3]),
#'   type = "l"
#' ) 
get_logpost <- function() {
  
  get_("logpost")
  
}

#' @export
#' @details The function `get_draws()` retrieves the proposed states from the
#' kernel function.
#' @rdname mcmc-output
get_draws <- function() {
  
  get_("draws")
  
}


#' @export
#' @rdname mcmc-output
get_elapsed <- function() {
  get_("time_end") - get_("time_start")
}

#' @export
#' @rdname mcmc-output
get_initial <- function() get_("initial")

#' @export
#' @rdname mcmc-output
get_fun <- function() get_("fun")

#' @export
#' @rdname mcmc-output
get_nsteps <- function() get_("nsteps")

#' @export
#' @rdname mcmc-output
get_seed <- function() get_("seed")

#' @export
#' @rdname mcmc-output
get_nchains <- function() get_("nchains")

#' @export
#' @rdname mcmc-output
get_burnin <- function() get_("burnin")

#' @export
#' @rdname mcmc-output
get_thin <- function() get_("thin")

#' @export
#' @rdname mcmc-output
get_kernel <- function() get_("kernel")

#' @export
#' @rdname mcmc-output
get_multicore <- function() get_("multicore")

#' @export
#' @rdname mcmc-output
get_conv_checker <- function() get_("conv_checker")

#' @export
#' @rdname mcmc-output
get_cl <- function() get_("cl")

#' @export
#' @rdname mcmc-output
get_progress <- function() get_("progress")

#' @export
#' @rdname mcmc-output
get_chain_id <- function() get_("chain_id")


#' Functions to interact with the main loop
#' 
#' You can use these functions to read variables, store, and retrieve data
#' during the MCMC process.
#' 
#' @export
#' @section Advanced usage:
#' The function [ith_step()] is a convenience function that provides
#' access to the environment within which the main loop of the MCMC call is
#' being evaluated. This is a wrapper of `MCMC_OUTPUT$loop_envir` that will
#' either return the value `x` or, if missing, the entire environment. If
#' `ith_step()` is called outside of the `MCMC` call, then it will return with
#' an error.
#' 
#' For example, if you wanted to print information if the current value
#' of logpost is greater than the previous value of logpost, you could define
#' the objective function as follows:
#' 
#' ```
#' f <- function(p) {
#' 
#'   i            <- ith_step("i")
#'   logpost_prev <- ith_step("logpost")[i - 1L]
#'   logpost_curr <- sum(dnorm(y - x*p, log = TRUE))
#'   
#'   if (logpost_prev < logpost_curr)
#'     cat("At a higher point!\n")
#'     
#'   return(logpost_curr)
#' 
#' }
#' ```
#' 
#' More examples below.
#' @examples 
#' #' # Getting the logpost -------------------------------------------------------
#' set.seed(23133)
#' x <- rnorm(200)
#' y <- x*2 + rnorm(200)
#' f <- function(p) {
#'   sum(dnorm(y - x*p, log = TRUE))
#' }
#' 
#' ans <- MCMC(fun = f, initial = c(0), nsteps=2000)
#' plot(get_logpost(), type = "l") # Plotting the logpost from the last run
#' 
#' 
#' # Printing information every 500 step ---------------------------------------
#' # for this we use ith_step()
#' 
#' f <- function(p) {
#' 
#'   # Capturing info from within the loop
#'   i      <- ith_step("i")
#'   nsteps <- ith_step("nsteps")
#'   
#'   if (!(i %% 500)) {
#'   
#'     cat(
#'       "////////////////////////////////////////////////////\n",
#'       "Step ", i, " of ", nsteps,". Values in the loop:\n",
#'       "theta0: ", ith_step("theta0"), "\n",
#'       "theta1: ", ith_step()$theta1, "\n",
#'       sep = ""
#'     )

#'   }
#'     
#' 
#'   sum(dnorm(y - x*p, log = TRUE))
#' }
#' 
#' MCMC(fun = f, initial = c(0), nsteps=2000, progress = FALSE, seed = 22)
#' # ////////////////////////////////////////////////////
#' # Step 500 of 2000. Values in the loop:
#' # theta0: 2.025379
#' # theta1: 1.04524
#' # ////////////////////////////////////////////////////
#' # Step 1000 of 2000. Values in the loop:
#' # theta0: 2.145967
#' # theta1: 0.2054037
#' # ////////////////////////////////////////////////////
#' # Step 1500 of 2000. Values in the loop:
#' # theta0: 2.211691
#' # theta1: 2.515361
#' # ////////////////////////////////////////////////////
#' # Step 2000 of 2000. Values in the loop:
#' # theta0: 1.998789
#' # theta1: 1.33034
#' 
#' 
#' # Printing information if the current logpost is greater than max -----------
#' f <- function(p) {
#' 
#'   i            <- ith_step("i")
#'   logpost_prev <- max(ith_step("logpost")[1:(i-1)])
#'   logpost_curr <- sum(dnorm(y - x*p, log = TRUE))
#'   
#'   # Only worthwhile after the first step
#'   if ((i > 1L) && logpost_prev < logpost_curr)
#'     cat("At a higher point!:", logpost_curr, ", step:", i,"\n")
#'     
#'   return(logpost_curr)
#' 
#' }
#' MCMC(fun = f, initial = c(0), nsteps=1000, progress = FALSE, seed = 22)
#' # At a higher point!: -357.3584 , step: 2 
#' # At a higher point!: -272.6816 , step: 6 
#' # At a higher point!: -270.9969 , step: 7 
#' # At a higher point!: -269.8128 , step: 24 
#' # At a higher point!: -269.7435 , step: 46 
#' # At a higher point!: -269.7422 , step: 543 
#' # At a higher point!: -269.7421 , step: 788 
#' @name mcmc-loop
#' @param x Name of the element to retrieve. If missing, it will return the entire
#' environment in which the main MCMC loop is running.
ith_step <- function(x) {
  
  if (is.null(MCMC_OUTPUT$loop_envir))
    stop("-ith_step()- should only be used within an -MCMC()- call.", call. = FALSE)
  
  if (missing(x))
    return(MCMC_OUTPUT$loop_envir)
  else
    return(get(x, envir = MCMC_OUTPUT$loop_envir))
  
}

#' @export
#' @rdname mcmc-loop
#' @param ... Named values to be appended to the user data.
#' @return The function [set_userdata()] returns [invisible()]. The only side
#' effect is appending the information by row.
#' @examples 
#' # Saving extra information --------------------------------------------------
#' data("lifeexpect")
#' 
#' # Defining the logposterior
#' logpost <- function(p) {
#' 
#'   set_userdata(i = ith_step("i"), sum_of_p = sum(p))
#' 
#'   with(lifeexpect, {
#'     sum(dnorm(age - (p[1] + p[2]*smoke + p[3]*female), sd = p[4], log = TRUE))
#'   })
#'   
#' }
#' 
#' # A kernel where sd is positive, the first is average age, so we 
#' # make it positive too
#' kern <- kernel_ram(lb = c(10, -20, -20, .0001), eps = .01)
#' ans <- MCMC(
#'   initial = c(70, -2, 2, 1), fun = logpost, kernel = kern, nsteps = 1000, seed = 1
#'   )
#' 
set_userdata <- function(...) {
  
  # In the first step, setup the data.frame
  i <- ith_step("i")
  if (i == 1L) {
    
    # Building the data
    MCMC_OUTPUT$ptr_usr[["userdata"]] <- do.call(data.frame, list(...))
    
    # Resizing to make sure we don't resize every time
    MCMC_OUTPUT$ptr_usr[["userdata"]] <- do.call(rbind, replicate(
      ith_step("nsteps"),
      MCMC_OUTPUT$ptr_usr[["userdata"]],
      simplify = FALSE
      ))
    
    return()
    
  }
  
  MCMC_OUTPUT$ptr_usr[["userdata"]][i, ] <- list(...)
  return()
  
}

#' @export
#' @rdname mcmc-loop
get_userdata <- function() {
  get_("userdata")
}