#' Append MCMC chains (objects of class [coda::mcmc])
#' @param ... A list of `mcmc` or `mcmc.list` class objects.
#' @export
append_chains <- function(...) UseMethod("append_chains")

#' @export
#' @rdname append_chains
append_chains.mcmc.list <- function(...) {
  
  dots <- list(...)
  
  nchains <- sapply(dots, coda::nchain)
  if (length(unique) != 1)
    stop()
  
  ans <- vector("list", nchains[1])
  for (i in 1:nchains[1])
    ans[[i]] <- do.call(append_chains.mcmc, lapply(dots, "[[", i))
  
  coda::as.mcmc.list(ans)
  
}

#' @export
#' @rdname append_chains
append_chains.mcmc <- function(...) {
  
  # retrieving the list of objects
  dots <- list(...)
  
  # Thinning
  thin <- sapply(dots, coda::thin)
  
  if (length(unique(thin)) != 1L)
    stop("All `mcmc` objects have to have the same `thin` parameter.",
         "Observed: ", paste(thin, collapse=,", "), ".", call.=FALSE)
  
  # Number of parameters
  nvar <- sapply(dots, coda::nvar)
  
  if (length(unique(nvar)) != 1L)
    stop("All `mcmc` objects have to have the same number of parameters.",
         "Observed: ", paste(nvar, collapse=,", "), ".", call.=FALSE)
  
  # Measuring lengths
  mcpar <- lapply(dots, coda::mcpar)
  start <- sapply(mcpar, "[[", 1)
  end   <- sapply(mcpar, "[[", 2)
  
  coda::mcmc(
    data  = unname(do.call(rbind, unclass(dots))),
    start = start[1],
    end   = sum(end),
    thin  = thin[1]
  )
  
  
}

