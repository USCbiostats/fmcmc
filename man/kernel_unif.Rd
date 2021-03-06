% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_unif.R
\name{kernel_unif}
\alias{kernel_unif}
\alias{kernel_unif_reflective}
\title{Uniform Transition Kernel}
\usage{
kernel_unif(min. = -1, max. = 1, fixed = FALSE, scheme = "joint")

kernel_unif_reflective(
  min. = -1,
  max. = 1,
  lb = min.,
  ub = max.,
  fixed = FALSE,
  scheme = "joint"
)
}
\arguments{
\item{min., max.}{Passed to \code{\link[stats:Uniform]{runif}}.}

\item{fixed, scheme}{For multivariate functions, sets the update plan.
See \code{\link[=plan_update_sequence]{plan_update_sequence()}}.}

\item{lb, ub}{Either a numeric vector or a scalar. Lower and upper bounds for
bounded kernels. When of length 1, the values are recycled to match the number
of parameters in the objective function.}
}
\value{
An object of class \link{fmcmc_kernel}. \code{fmcmc_kernel} objects are intended
to be used with the \code{\link[=MCMC]{MCMC()}} function.
}
\description{
Uniform Transition Kernel
}
\details{
The \code{kernel_unif} function provides a uniform transition kernel. This (symmetric)
kernel function by default adds the current status values between [-1,1].

The \code{kernel_unif_reflective} is similar to \code{kernel_unif} with the
main difference that proposals are bounded to be within \verb{[lb, ub]}.
}
\examples{
# Multivariate setting with 4 parameters in which we set the kernel to make
# proposals one parameter at-a-time in a random ordering.
kern <- kernel_unif(scheme = "random")
}
\seealso{
Other kernels: 
\code{\link{kernel_adapt}()},
\code{\link{kernel_mirror}},
\code{\link{kernel_new}()},
\code{\link{kernel_normal}()},
\code{\link{kernel_ram}()}
}
\concept{kernels}
