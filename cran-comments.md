## Test environments

* Ubuntu 16.04.6 LTS (local), R 3.6.1
* Ubuntu 16.04.6 LTS (on travis-ci), R 3.6.1
* Mac OS X (on travis-ci), R 3.6.1
* Windows Server 2012 R2 x64 (AppVeyor), R 3.6.1

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## CRAN comments:

This re-submission addresses the following CRAN comments:

* Adding `\donttest` to sections of the manual where examples take more than 5s.

* "[...] what you mean with 'A friendly' in the title." It means flexible. Just
  added that to the Description.

* "If there are references describing the methods in your package, please add [...]"
  Added a general reference for the package (The Handbook of MCMC).

* Fixed the "Unexcutable code" on kernels.Rd and MCMC.Rd (thanks for spotting
  those!).

Thanks!


