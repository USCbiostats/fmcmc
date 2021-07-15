#' Append MCMC chains (objects of class [coda::mcmc])
#' 
#' Combines two or more MCMC runs into a single run. If runs have
#' multiple chains, it will check that all have the same number of chains, and
#' it will join chains using the [rbind] function.
#' 
#' @param ... A list of `mcmc` or `mcmc.list` class objects.
#' @return If `mcmc.list`, an object of class `mcmc.list`, otherwise,
#' an object of class `mcmc`.
#' @export
append_chains <- function(...) UseMethod("append_chains")

#' @export
# @rdname append_chains
append_chains.default <- function(...) {
  
  dots <- list(...)
  dots <- dots[sapply(dots, length) > 0L]
  
  if (coda::is.mcmc(dots[[1]]))
    do.call(append_chains.mcmc, dots)
  else if (coda::is.mcmc.list(dots[[1]]))
    do.call(append_chains.mcmc, dots)
  else
    stop("No method available to append these chains.", call. = FALSE)
  
}

#' @export
# @rdname append_chains
append_chains.mcmc.list <- function(...) {
  
  dots <- list(...)
  
  if (length(dots) == 1L)
    return(dots[[1]])
  
  nchains <- sapply(dots, coda::nchain)
  if (length(unique(nchains)) != 1)
    stop(
      "All mcmc.list objects must have the same number of chains. ",
      "The passed objects have ", paste(nchains, collapse = ", "),
      " respectively.",
      call. = FALSE
      )
  
  ans <- vector("list", nchains[1])
  for (i in 1:nchains[1])
    ans[[i]] <- do.call(append_chains.mcmc, lapply(dots, "[[", i))
  
  coda::as.mcmc.list(ans)
  
}

#' @export
# @rdname append_chains
append_chains.mcmc <- function(...) {
  
  # retrieving the list of objects
  dots <- list(...)
  
  # A single one
  if (length(dots) == 1L)
    return(dots[[1]])
  
  # Thinning
  thin <- sapply(dots, coda::thin)
  
  if (length(unique(thin)) != 1L)
    stop("All `mcmc` objects have to have the same `thin` parameter.",
         "Observed: ", paste(thin, collapse=", "), " respectively.", call.=FALSE)
  
  # Number of parameters
  nvar <- sapply(dots, coda::nvar)
  
  if (length(unique(nvar)) != 1L)
    stop("All `mcmc` objects have to have the same number of parameters.",
         "Observed: ", paste(nvar, collapse=", "), " respectively.", call.=FALSE)
  
  # Measuring lengths
  mcpar <- lapply(dots, coda::mcpar)
  start <- sapply(mcpar, "[[", 1)
  end   <- sapply(mcpar, "[[", 2)
  
  # checking rownames
  rn <- lapply(dots, function(d) {
    as.integer(rownames(unclass(d)))
  })
  
  midnames <- Map(function(rn., add) {
    as.character(rn. - rn.[1] + add + thin[1])
  }, rn = rn[-1], add=cumsum(end[-length(end)])
  )
  
  rn <- c(as.character(rn[[1]]), unlist(midnames))
  
  # Correcting endings
  end[-1] <- end[-1] + thin[-1] - start[-1]
  
  
  dat <- do.call(rbind, unclass(dots))
  rownames(dat) <- rn
  
  coda::mcmc(
    data  = dat,
    start = start[1],
    end   = sum(end),
    thin  = thin[1]
  )
  
  
}

