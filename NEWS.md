# fmcmc 0.6-0

* Added support for using `NA` to specify unbounded parameters in kernel
  functions (closes [#22](https://github.com/USCbiostats/fmcmc/issues/22), as suggested by [@dmi3kno](https://github.com/dmi3kno)). This provides a more
  intuitive alternative to `.Machine$double.xmax`. For example:
  `kernel_ram(lb = c(alpha = NA, beta = 0), ub = c(alpha = NA, beta = 1))`
  where `alpha` is unbounded and `beta` is bounded on [0, 1].

* Added an example using cluster objects.


# fmcmc 0.5-2

* Adressing [`roxygen2` issue #1491](https://github.com/r-lib/roxygen2/issues/1491)


# fmcmc 0.5-1

* Fixed issue #21: Restricting search scope in MCMC temp environment data.

* Removed annoying warning when using convergence checker.


# fmcmc 0.5-0

* The function `fun` passed to `MCMC` is now called two times less. It shouldn't
  significantly affect any previous results.
  
* `convergence_gelman` now stores the Gelman and Rubin's statistics in the correct
  order, i.e., the most recent at the end of the array in `convergence_data_get("val")`.
  
* Users can now pass `seed` to `MCMC`. If `is.null(seed) != TRUE`, then `seed` is
  passed to `set.seed()`.
  
* The function `convergence_auto()` now behaves as expected. Before, it was not checking
  convergence.
  
* The set of functions `last_*` and `LAST_MCMC` will be deprecated in favor of
  `get_*` and `MCMC_OUTPUT`.
  
* The new function `get_logpost()` returns the computed values of the objective
  function from the last `MCMC` run.
  
* The new function `get_draws()` returns the MCMC draws from the kernel's 
  proposal function (proposed states).
  
* The new function `set_userdata(...)` allows storing information into a data.frame
  as the MCMC process runs. Users can retrieve the data with the function
  `get_userdata()`.

* The new function `ith_step()` provides access to objects within the MCMC
  loop during the run. The new function comes with a vignette that illustrates
  its usage.
  
* The function `append_chains()` was randomly dropping one sample of the final
  set.
  
* A new artificial dataset `lifeexpect` is shipped with the package. This simulates
  1,000 observations of `age` at death using US's statistics.
  

# fmcmc 0.4-0

* `kernel_am` and `kernel_ram` no longer fail when at least one parameter is
  an offset (`fixed = TRUE` for some parameter).

* Now `kernel_ram` tries first to find the cholesky decomp. If it fails, then
  it uses `Matrix::nearPD` and re-tries. This is following what is done in the
  `adaptMCMC` package.

* Workflow for running MCMC with `conv_checker` re-designed (less error prone).

* Environments `LAST_RUN` and `LAST_CONV_CHECK` provide information about the
  last call to `MCMC` and the corresponding convergence checker. Users can
  access these environments via getter and setter functions.
  
* `MCMC` with convergence checker now reports the status of the convergence
  statistic using the `LAST_CONV_CHECK` environment and corresponding
  functions.
  
* The functions to compute mean and variance recursively now allow us to do so
  using windows.
  

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
