# Test get_userdata and set_userdata -------------------------------------------

# Simple test with single chain: counter and predefined vector -----------------
set.seed(123)
n <- 100
predefined_vec <- rnorm(n)

# Objective function that sets user data at each step
fun_single <- function(p) {
  i <- ith_step("i")
  # Store the step counter and the predefined vector value for this step
  set_userdata(step = i, value = predefined_vec[i])
  dnorm(p, log = TRUE)
}

ans_single <- MCMC(
  fun     = fun_single,
  initial = 0,
  nsteps  = n,
  kernel  = kernel_normal(scale = 1),
  seed    = 42
)

# Retrieve user data and verify
userdata_single <- get_userdata()

# Test: should be a data.frame
expect_true(is.data.frame(userdata_single))

# Test: should have n rows (one per step)
expect_equal(nrow(userdata_single), n)

# Test: step counter should be 1:n
expect_equal(userdata_single$step, 1:n)

# Test: values should match the predefined vector
expect_equal(userdata_single$value, predefined_vec)


# Test with multiple chains ----------------------------------------------------
set.seed(456)
n_multi <- 50
nchains_multi <- 2
predefined_vec_multi <- rnorm(n_multi)

fun_multi <- function(p) {
  i <- ith_step("i")
  # Store the step counter and the predefined vector value for this step
  set_userdata(step = i, value = predefined_vec_multi[i])
  dnorm(p, log = TRUE)
}

ans_multi <- MCMC(
  fun     = fun_multi,
  initial = 0,
  nsteps  = n_multi,
  nchains = nchains_multi,
  kernel  = kernel_normal(scale = 1),
  seed    = 42
)

# Retrieve user data and verify
userdata_multi <- get_userdata()

# Test: should be a list of length nchains_multi (one per chain)
expect_true(is.list(userdata_multi))
expect_equal(length(userdata_multi), nchains_multi)

# Test: each element should be a data.frame with n_multi rows
for (chain in seq_len(nchains_multi)) {
  expect_true(is.data.frame(userdata_multi[[chain]]))
  expect_equal(nrow(userdata_multi[[chain]]), n_multi)
  expect_equal(userdata_multi[[chain]]$step, 1:n_multi)
  expect_equal(userdata_multi[[chain]]$value, predefined_vec_multi)
}
