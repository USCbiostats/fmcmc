z <- function(x = 1, w = 4, multicore = FALSE, nchains = 1L, ...) {
  
  if (nchains > 1) {
    
    # Recovering arguments of the call
    used_args <- formals()
    used_args <- used_args[setdiff(names(used_args), "...")]
    
    # Baseline call
    cll <- if (multicore) {
      as.call(
        c(list(
        quote(parallel::clusterApply),
        cl = str2lang("cl"),
        x  = str2lang(paste0("1:", nchains))
      ), used_args)
      )
    } else {
      as.call(c(list(
        quote(lapply),
        X  = str2lang(paste0("1:", nchains))
      ), used_args)
      )
    }
    
    # Adjusting call for the parallel case
    call. <- match.call(expand.dots = TRUE)
    call. <- lapply
    call.$multicore <- FALSE
    call.$nchains   <- 1L
    
    fun <- sprintf(
      "function(i, %s) {\n%s\n}",
      paste(names(formals()), collapse="., "),
      deparse(call.)
    )
    
    if (multicore)
      cll$fun <- str2lang(fun)
    else
      cll$FUN <- str2lang(fun)
    
    print(cll)
    
    return(eval(cll))
  }
  
  x + w
  
}

cl <- parallel::makePSOCKcluster(4)
z(nchains = 4L, multicore = TRUE, omega = 1)