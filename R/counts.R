sample_size <- function(size, burnin = 0, thin = 1) {
  
  list(
    niter  = size*thin + burnin,
    burnin = burnin,
    thin   = thin
  )
  
}

sample_burn <- function(x, burnin) {
  x[-1:burnin, , drop=FALSE]
}

sample_thin <- function(x, thin) {
  
}