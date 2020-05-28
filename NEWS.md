# fmcmc 0.3-0-99

* `kernel_am` and `kernel_ram` no longer fail when at least one parameter is
  an offset (`fixed = TRUE` for some parameter).

* Now `kernel_ram` tries first to find the cholesky decomp. If it fails, then
  it uses `Matrix::nearPD` and re-tries. This is following what is done in the
  `adaptMCMC` package.
  

# fmcmc 0.3-0

* Adding Vihola (2012)'s Robust Adaptive Metropolis, Haario et al. (2001)'s
  Adaptive Metropolis, and Thawornwattana et al. (2018)'s mirror kernel
  algorithms.

* The argument `progress` is no longer ignored. When set to `TRUE`, the function
  will print the progress of the MCMC algorithm.

* Improved coverage and fixed minor bugs.

* When running with convergence check, fixed parameters (offset), as tagged in
  the `fmcmc_kernel` object will be excluded from the call to `conv_checker`.


# fmcmc 0.2-0

* Added a `NEWS.md` file to track changes to the package.
