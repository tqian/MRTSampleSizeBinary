# generate data for tests -------------------------------------------------
set.seed(1)



total_dp <- 10 # total number of decision points
# p and q are not directly used in the function, so I commented them out.
# p <- 2
# q <- 2

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



g_new <- cbind(rep(1, total_dp), 1:total_dp, (1:total_dp)^2)
alpha_new <- as.matrix(c(-0.2, -0.1, .01), ncol = 1)

f_new <- cbind(rep(1, total_dp), 1:total_dp, (1:total_dp)^2) # f_t = (1, t)
beta_new <- as.matrix(c(0.15, - 0.01, -.1), ncol = 1)

# mrt_binary_ss tests ------------------------------------



# check outputs for valid inputs
context("testing that outputs match original function")
test_that(
  "check TQ's sample",
  {
    expect_equal(
      mrt_binary_ss(tau_t, f_t, g_t, beta, 
                    alpha, p_t, gamma, b, TRUE),
      274.0055127)
  }
)

test_that(
  "check that the round up feature is working",
  {
    expect_equal(
      mrt_binary_ss(tau_t, f_t, g_t, beta, 
                    alpha, p_t, gamma, b, FALSE),
      275)
  }
)


test_that(
  "check quadratic example",
  {
    expect_equal(
      mrt_binary_ss(tau_t, f_new, g_new, beta_new, 
                    alpha_new, p_t, gamma, b, TRUE),
      32.003286527546599)
  }
)

test_that(
  "check example with different dimension f, g",
  {
    expect_equal(
      mrt_binary_ss(tau_t, f_t, g_new, beta, 
                    alpha_new, p_t, gamma, b, TRUE),
      184.43823325060882)
  }
)

test_that(
  "check that it works as an 'inverse' of mrt_binary_power",
  {
    expect_equal(
      mrt_binary_ss(tau_t, f_t, g_new, beta, 
                    alpha_new, p_t, gamma,
                    1-mrt_binary_power(
                      tau_t, f_t, g_new, beta,
                      alpha_new, p_t, gamma, 50), TRUE),
      50)
  }
)

# test warning
test_that(
  "check example with invalid dimension f, g",
  {
    expect_warning(
      mrt_binary_ss(tau_t, f_new, g_t, beta_new, 
                    alpha, p_t, gamma, b, TRUE))
  }
)

f_warn <- cbind(rep(1, 10), rep(c(1,0), times=5))

test_that(
  "check for warning about p_t*f_t not being in span of g_t",
  {
    expect_warning(
      mrt_binary_ss(tau_t, f_warn, g_t, beta, 
                    alpha, p_t, gamma, 0.4))
  }
)
# test errors
test_that(
  "check example with invalid dimension f and beta",
  {
    expect_error(
      mrt_binary_ss(tau_t, f_t, g_t, beta_new, 
                    alpha, p_t, gamma, b, TRUE),
      "Dimensions of f_t and beta do not agree.")
  }
)


test_that(
  "check example with invalid dimension g and alpha",
  {
    expect_error(
      mrt_binary_ss(tau_t, f_t, g_t, beta, 
                    alpha_new, p_t, gamma, b, TRUE),
      "Dimensions of g_t and alpha do not agree.")
  }
)

test_that(
  "test for incorrect number of time points",
  {
    expect_error(
      mrt_binary_ss(rep(.4, times=2), f_t, g_t, beta, 
                    alpha_new, p_t, gamma, b, TRUE),
      "All arguments must agree on number of time points.")
  }
)

test_that(
  "test for invalid type I error",
  {
    expect_error(
      mrt_binary_ss(tau_t, f_t, g_t, beta, 
                    alpha, p_t, -3, b, TRUE),
      "gamma, type I error, should be between 0 and 1")
  }
)

test_that(
  "test for invalid type II error",
  {
    expect_error(
      mrt_binary_ss(tau_t, f_t, g_t, beta, 
                    alpha, p_t, gamma, 100*b, TRUE),
      "b, type II error, should be between 0 and 1")
  }
)

test_that(
  "test for incorrect type of f_t",
  {
    expect_error(
      mrt_binary_ss(tau_t, c(1,3,0), g_t, beta, 
                    alpha, p_t, gamma, 100*b, TRUE),
      "f_t and g_t should be matrices")
  }
)

test_that(
  "test for incorrect type of g_t",
  {
    expect_error(
      mrt_binary_ss(tau_t, f_t, 0, beta, 
                    alpha, p_t, gamma, 100*b, TRUE),
      "f_t and g_t should be matrices")
  }
)
