# generate vals for tests -------------------------------------------------
set.seed(1)



total_dp <- 10 # total number of decision points

# expected availability E(I_t) for t = 1,...,total_dp
tau_t <- rep(0.8, total_dp) 

p_t <- rep(0.4, total_dp) # randomization probability over time
gamma <- 0.05 # type I error
b <- 0.2 # type II error; power = 1 - b

### specify g_t and alpha ###
g_t <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1,g_t)
alpha <- as.matrix(c(-0.2, -0.1), ncol = 1)

# check that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 0) is between 0 and 1.
# E(Y_{t+1} = 1 | I_t = 1, A_t = 0) for t = 1,...,total_dp
mu0_t <- exp(g_t %*% alpha) 

### specify f_t and beta ###
f_t <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1, t)
beta <- as.matrix(c(0.15, - 0.01), ncol = 1)

# check that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 1) is between 0 and 1.
# MEE(t) for t = 1,...,total_dp
mee_t <- f_t %*% beta 

# E(Y_{t+1} = 1 | I_t = 1, A_t = 1) for t = 1,...,total_dp
mu1_t <- mu0_t * exp(mee_t) 


p_t <- rep(0.4, total_dp) # randomization probability over time
gamma <- 0.05 # type I error
b <- 0.2 # type II error; power = 1 - b

### specify g_t and alpha ###
g_t <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1,g_t)
alpha <- as.matrix(c(-0.2, -0.1), ncol = 1)

# check that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 0) is between 0 and 1.
# E(Y_{t+1} = 1 | I_t = 1, A_t = 0) for t = 1,...,total_dp
mu0_t <- exp(g_t %*% alpha) 

### specify f_t and beta ###
f_t <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1, t)
beta <- as.matrix(c(0.15, - 0.01), ncol = 1)

# check that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 1) is 0 and 1.
mee_t <- f_t %*% beta # MEE(t) for t = 1,...,total_dp
# E(Y_{t+1} = 1 | I_t = 1, A_t = 1) for t = 1,...,total_dp
mu1_t <- mu0_t * exp(mee_t) 


# make matrices for quadratic tests
g_new <- cbind(rep(1, total_dp), 1:total_dp, (1:total_dp)^2)
alpha_new <- as.matrix(c(-0.2, -0.1, .01), ncol = 1)

f_new <- cbind(rep(1, total_dp), 1:total_dp, (1:total_dp)^2)
beta_new <- as.matrix(c(0.15, - 0.01, -.1), ncol = 1)


# mrt_binary_power tests -----------------------------------------

size1 <- mrt_binary_ss(tau_t, f_t, g_t, beta, 
                       alpha, p_t, gamma, .3, exact=TRUE)

size2 <- mrt_binary_ss(tau_t, f_t, g_t, beta, 
                       alpha, p_t, gamma, .1, exact=TRUE)
test_that(
  "check TQ's sample",
  {
    expect_warning(
      mrt_binary_power(tau_t, f_t, g_t, beta, 
                       alpha, p_t, gamma, size1),
      "n should be an integer")
  }
)

test_that(
  "check TQ's sample at different sample size",
  {
    expect_equal(
      mrt_binary_power(tau_t, f_t, g_t, beta, 
                       alpha, p_t, gamma, round(size2)),
      .9, tol=0.0001)
  }
)

test_that(
  "check that it works as an 'inverse' of mrt_binary_ss",
  {
    expect_equal(
      mrt_binary_power(tau_t, f_t, g_new, beta, 
                       alpha_new, p_t, gamma,
                       mrt_binary_ss(
                         tau_t, f_t, g_new, beta,
                         alpha_new, p_t, gamma, 1-1/pi, FALSE)),
      1/pi, tolerance=.01)
  }
)

# test warnings
test_that(
  "check example with invalid dimension f, g",
  {
    expect_warning(
      mrt_binary_power(tau_t, f_new, g_t, beta_new, 
                       alpha, p_t, gamma, 20))
  }
)

f_warn <- cbind(rep(1, 10), rep(c(1,0), times=5))

test_that(
  "check for warning about p_t*f_t not being in span of g_t",
  {
    expect_warning(
      mrt_binary_power(tau_t, f_warn, g_t, beta, 
                       alpha, p_t, gamma, round(size2)))
  }
)

# test errors
test_that(
  "check that min sample size stops function",
  {
    expect_error(
      mrt_binary_power(tau_t, f_t, g_new, beta, 
                       alpha_new, p_t, gamma,
                       1),
      "n is too small"
    )
  }
)

test_that(
  "check that min sample size stops function",
  {
    expect_error(
      mrt_binary_power(tau_t, f_t, g_new, beta, 
                       alpha_new, p_t, gamma,
                       -1),
      "n is too small"
    )
  }
)



test_that(
  "check example with invalid dimension f and beta",
  {
    expect_error(
      mrt_binary_power(tau_t, f_t, g_t, beta_new, 
                       alpha, p_t, gamma, 1000),
      "Dimensions of f_t and beta do not agree.")
  }
)


test_that(
  "check example with invalid dimension g and alpha",
  {
    expect_error(
      mrt_binary_power(tau_t, f_t, g_t, beta, 
                       alpha_new, p_t, gamma, 99),
      "Dimensions of g_t and alpha do not agree.")
  }
)

test_that(
  "test for incorrect number of time points",
  {
    expect_error(
      mrt_binary_power(rep(.4, times=2), f_t, g_t, beta, 
                       alpha_new, p_t, gamma, 55),
      "All arguments must agree on number of time points.")
  }
)

test_that(
  "test for invalid type I error",
  {
    expect_error(
      mrt_binary_power(tau_t, f_t, g_t, beta, 
                       alpha, p_t, -3, 44),
      "gamma, type I error, should be between 0 and 1")
  }
)


test_that(
  "test for invalid type I error",
  {
    expect_error(
      mrt_binary_power(tau_t, f_t, g_t, beta, 
                       alpha, p_t, 'calc', 44),
      "gamma, type I error, should be between 0 and 1")
  }
)

test_that(
  "test for incorrect type of f_t",
  {
    expect_error(
      mrt_binary_power(tau_t, "test", g_t, beta, 
                       alpha, p_t, 'calc', .44),
      "f_t and g_t should be matrices")
  }
)

test_that(
  "test for invalid type I error",
  {
    expect_error(
      mrt_binary_power(tau_t, f_t, pi, beta, 
                       alpha, p_t, 'calc', .44),
      "f_t and g_t should be matrices")
  }
)
