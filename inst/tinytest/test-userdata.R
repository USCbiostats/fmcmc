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


# Test add_userdata() function -------------------------------------------------

# Test 1: Basic functionality with single chain --------------------------------
set.seed(789)
n_add <- 50

fun_add_single <- function(p) {
  i <- ith_step("i")
  set_userdata(iteration = i, param_squared = p^2)
  dnorm(p, log = TRUE)
}

ans_add_single <- MCMC(
  fun     = fun_add_single,
  initial = 0,
  nsteps  = n_add,
  kernel  = kernel_normal(scale = 1),
  seed    = 123
)

# Add userdata to the mcmc object
combined_single <- add_userdata(ans_add_single)

# Test: should be mcmc object (single chain returns mcmc, not mcmc.list)
expect_true(inherits(combined_single, "mcmc"))

# Test: combined object should have more columns than original
expect_true(ncol(combined_single) > ncol(ans_add_single))

# Test: should have the original columns plus userdata columns
orig_cols <- colnames(ans_add_single)
combined_cols <- colnames(combined_single)
expect_true(all(orig_cols %in% combined_cols))
expect_true("iteration" %in% combined_cols)
expect_true("param_squared" %in% combined_cols)

# Test: number of iterations should match
expect_equal(nrow(combined_single), nrow(ans_add_single))


# Test 2: Basic functionality with multiple chains -----------------------------
set.seed(321)
n_add_multi <- 60
nchains_add <- 3

fun_add_multi <- function(p) {
  i <- ith_step("i")
  set_userdata(iter = i, log_p = log(abs(p) + 1))
  dnorm(p, log = TRUE)
}

ans_add_multi <- MCMC(
  fun     = fun_add_multi,
  initial = 0,
  nsteps  = n_add_multi,
  nchains = nchains_add,
  kernel  = kernel_normal(scale = 1),
  seed    = 456
)

# Add userdata to the mcmc object
combined_multi <- add_userdata(ans_add_multi)

# Test: should be mcmc.list
expect_true(inherits(combined_multi, "mcmc.list"))

# Test: should have the correct number of chains
expect_equal(length(combined_multi), nchains_add)

# Test: all chains should have userdata columns
for (chain_idx in seq_len(nchains_add)) {
  combined_cols <- colnames(combined_multi[[chain_idx]])
  expect_true("iter" %in% combined_cols)
  expect_true("log_p" %in% combined_cols)
  expect_equal(nrow(combined_multi[[chain_idx]]), n_add_multi)
}


# Test 3: Verify data integrity ------------------------------------------------
# Check that userdata values are correctly combined with mcmc data
set.seed(654)
n_integrity <- 30

fun_integrity <- function(p) {
  i <- ith_step("i")
  set_userdata(step_id = i, double_p = 2 * p)
  dnorm(p, log = TRUE)
}

ans_integrity <- MCMC(
  fun     = fun_integrity,
  initial = 1,
  nsteps  = n_integrity,
  kernel  = kernel_normal(scale = 0.5),
  seed    = 789
)

combined_integrity <- add_userdata(ans_integrity)

# Test: step_id should match row indices
userdata_integrity <- get_userdata()
expect_equal(userdata_integrity$step_id, 1:n_integrity)

# Test: double_p should be approximately 2 * parameter values
# (allowing for small numerical differences)
param_values <- as.numeric(ans_integrity)
expected_double <- 2 * param_values
actual_double <- as.numeric(combined_integrity[, "double_p"])
expect_equal(actual_double, expected_double)


# Test 4: Edge case - check error when number of chains mismatch --------------
# Create a scenario where get_userdata() doesn't match the mcmc.list
# This tests the validation in add_userdata()

set.seed(111)
n_edge1 <- 20

fun_edge1 <- function(p) {
  i <- ith_step("i")
  set_userdata(idx = i)
  dnorm(p, log = TRUE)
}

ans_edge1 <- MCMC(
  fun     = fun_edge1,
  initial = 0,
  nsteps  = n_edge1,
  nchains = 2,
  kernel  = kernel_normal(scale = 1),
  seed    = 111
)

# Now run a different MCMC with different number of chains
# This will overwrite the userdata with different chain count
ans_edge1_wrong <- MCMC(
  fun     = fun_edge1,
  initial = 0,
  nsteps  = n_edge1,
  nchains = 3,  # Different number of chains
  kernel  = kernel_normal(scale = 1),
  seed    = 222
)

# Test: should error because number of chains don't match
# ans_edge1 has 2 chains but get_userdata() now returns 3 chains
expect_error(
  add_userdata(ans_edge1),
  "Number of chains does not match"
)


# Test 5: Verify mcpar attributes are preserved --------------------------------
set.seed(555)
n_mcpar <- 40

fun_mcpar <- function(p) {
  i <- ith_step("i")
  set_userdata(iteration = i)
  dnorm(p, log = TRUE)
}

ans_mcpar <- MCMC(
  fun     = fun_mcpar,
  initial = 0,
  nsteps  = n_mcpar,
  kernel  = kernel_normal(scale = 1),
  seed    = 333
)

combined_mcpar <- add_userdata(ans_mcpar)

# Test: mcpar attributes should match
orig_mcpar <- coda::mcpar(ans_mcpar)
combined_mcpar_attr <- coda::mcpar(combined_mcpar)
expect_equal(orig_mcpar, combined_mcpar_attr)
