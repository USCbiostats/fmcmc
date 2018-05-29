
plot((function(x) dgamma(x, shape = 1, scale=1/50))(seq(0, 1, length.out = 2000)), type="l")

dat <- rgamma(10, 1, 1/50)

f <- function(p, ...) prod(dgamma(dat, p[1], p[2]))
x <- c(1, 1)

# Adaptive kernel
prob <- replicate(1e3, sdml_adaptive_kernel(
  x = c(1,1),
  f = f, alpha=.9
)$h)

plot(log(prob[!is.nan(prob)]), type="l")

# Non adaptive kernel
# normal_prop()