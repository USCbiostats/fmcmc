#' Checks the initial values of the MCMC
#' 
#' This function is for internal use only.
#' 
#' @param initial Either a vector or matrix,.
#' @param nchains Integer scalar. Number of chains.
#' @details When `initial` is a vector, the values are recycled to form a matrix of
#' size `nchains * length(initial)`.
#' @return
#'  A named matrix.
#' @examples
#' init <- c(.4, .1)
#' check_initial(init, 1)
#' check_initial(init, 2)
#' 
#' init <- matrix(1:9, ncol=3)
#' check_initial(init, 3)
#' 
#' # check_initial(init, 2) # Returns an error
#' 
#' @export
check_initial <- function(initial, nchains) {
  
  # If function, then return function
  if (is.vector(initial, mode = "numeric")) {
    
    if (nchains > 1)
      warning("While using multiple chains, a single initial point has been ",
              "passed via `initial`: c(", paste(initial, collapse=", "),
              "). The values will be recycled. Ideally you would want to ",
              "start each chain from different locations.", call. = FALSE)
  
    # The number of dimensions cannot be null
    if (!length(initial))
      stop("The `initial` vector is of length zero.", call. = FALSE)
      
    initial <- matrix(initial, ncol=length(initial), nrow=nchains,
                      dimnames = list(NULL, names(initial)), byrow=TRUE)
    
  } else if (inherits(initial, "mcmc.list")) {
    
    
    
  } else if (!is.matrix(initial))
    stop("When `initial` is not a numeric vector, it should be a matrix. Right now it is ",
         "an object of class `", class(initial), "`.", call. = FALSE)
  else if (is.matrix(initial) && nrow(initial) != nchains)
    stop("The number of rows of `initial` (", nrow(initial),") must coincide with the ",
         "number of chains (", nchains, ").", call. = FALSE)
  
  # Checking parnames
  if (!length(colnames(initial)))
    colnames(initial) <- paste0("par", 1:ncol(initial))
  
  # Returns the corresponding seed
  initial
  
}

