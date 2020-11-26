library(testthat)


# simulate data for tests -------------------------------------------------
# TQ added
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
mu0_t # look at the mu0_t values
stopifnot(all(mu0_t < 1) & all(mu0_t > 0))

### specify f_t and beta ###
f_t <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1, t)
beta <- as.matrix(c(0.15, - 0.01), ncol = 1)

# check that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 1) is between 0 and 1.
# MEE(t) for t = 1,...,total_dp
mee_t <- f_t %*% beta 

# E(Y_{t+1} = 1 | I_t = 1, A_t = 1) for t = 1,...,total_dp
mu1_t <- mu0_t * exp(mee_t) 
mu1_t # look at the mu1_t values
stopifnot(all(mu1_t < 1) & all(mu1_t > 0))


p_t <- rep(0.4, total_dp) # randomization probability over time
gamma <- 0.05 # type I error
b <- 0.2 # type II error; power = 1 - b

### specify g_t and alpha ###
g_t <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1,g_t)
alpha <- as.matrix(c(-0.2, -0.1), ncol = 1)

# check that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 0) is between 0 and 1.
# E(Y_{t+1} = 1 | I_t = 1, A_t = 0) for t = 1,...,total_dp
mu0_t <- exp(g_t %*% alpha) 
mu0_t # look at the mu0_t values
stopifnot(all(mu0_t < 1) & all(mu0_t > 0))

### specify f_t and beta ###
f_t <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1, t)
beta <- as.matrix(c(0.15, - 0.01), ncol = 1)

# check that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 1) is 0 and 1.
mee_t <- f_t %*% beta # MEE(t) for t = 1,...,total_dp
# E(Y_{t+1} = 1 | I_t = 1, A_t = 1) for t = 1,...,total_dp
mu1_t <- mu0_t * exp(mee_t) 
mu1_t # look at the mu1_t values
stopifnot(all(mu1_t < 1) & all(mu1_t > 0))



g_new <- cbind(rep(1, total_dp), 1:total_dp, (1:total_dp)^2)
alpha_new <- as.matrix(c(-0.2, -0.1, .01), ncol = 1)

f_new <- cbind(rep(1, total_dp), 1:total_dp, (1:total_dp)^2) # f_t = (1, t)
beta_new <- as.matrix(c(0.15, - 0.01, -.1), ncol = 1)

# calculate_mrt_bin_samplesize_f tests ------------------------------------

# check outputs for valid inputs
context("testing that outputs match original function")
test_that(
  "check TQ's sample",
  {
    expect_equal(
      calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta, 
                                     alpha, p_t, gamma, b, TRUE),
    274.0055127)
  }
)

test_that(
  "check that the round up feature is working",
  {
    expect_equal(
      calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta, 
                                     alpha, p_t, gamma, b, FALSE),
      275)
  }
)


test_that(
  "check quadratic example",
  {
    expect_equal(
      calculate_mrt_bin_samplesize_f(tau_t, f_new, g_new, beta_new, 
                                     alpha_new, p_t, gamma, b, TRUE),
      32.003286527546599)
  }
)

test_that(
  "check example with different dimension f, g",
  {
    expect_equal(
      calculate_mrt_bin_samplesize_f(tau_t, f_t, g_new, beta, 
                                     alpha_new, p_t, gamma, b, TRUE),
      184.43823325060882)
  }
)

test_that(
  "check that it works as an 'inverse' of calculate_mrt_bin_power_f",
  {
    expect_equal(
      calculate_mrt_bin_samplesize_f(tau_t, f_t, g_new, beta, 
                                alpha_new, p_t, gamma,
                                1-calculate_mrt_bin_power_f(
                                  tau_t, f_t, g_new, beta,
                                  alpha_new, p_t, gamma, 10), FALSE),
      10)
  }
)

# test warning
test_that(
  "check example with invalid dimension f, g",
  {
    expect_warning(
      calculate_mrt_bin_samplesize_f(tau_t, f_new, g_t, beta_new, 
                                     alpha, p_t, gamma, b, TRUE),
      "f should lie in span of g")
  }
)

# test errors
test_that(
  "check example with invalid dimension f and beta",
  {
    expect_error(
      calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta_new, 
                                     alpha, p_t, gamma, b, TRUE),
      "Dimensions of f_t and beta do not agree.")
  }
)


test_that(
  "check example with invalid dimension g and alpha",
  {
    expect_error(
      calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta, 
                                     alpha_new, p_t, gamma, b, TRUE),
      "Dimensions of g_t and alpha do not agree.")
  }
)

test_that(
  "test for incorrect number of time points",
  {
    expect_error(
      calculate_mrt_bin_samplesize_f(rep(.4, times=2), f_t, g_t, beta, 
                                     alpha_new, p_t, gamma, b, TRUE),
      "All arguments must agree on number of time points.")
  }
)

test_that(
  "test for invalid type I error",
  {
    expect_error(
      calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta, 
                                     alpha, p_t, -3, b, TRUE),
      "gamma, type I error, should be between 0 and 1")
  }
)

test_that(
  "test for invalid type II error",
  {
    expect_error(
      calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta, 
                                     alpha, p_t, gamma, 100*b, TRUE),
      "b, type II error, should be between 0 and 1")
  }
)
calculate_mrt_bin_samplesize_f(tau_t, f_new, g_t, beta_new, 
                              alpha, p_t, gamma, b, TRUE)

# compute_m_sigma tests ---------------------------------------------------
test_that(
  "check if first error check of invalid probabilities catches error",
  {
    expect_error(
      compute_m_sigma(tau_t, f_t, -g_t, beta, 
                                     alpha, p_t),
      message="g_t and alpha values led to invalid probabilities")
  }
)

test_that(
  "check if second error check of invalid probabilities catches error",
  {
    expect_error(
      compute_m_sigma(tau_t, f_t, g_t, beta+1, 
                                     alpha, p_t, gamma, b),
      message="f_t and beta values led to invalid probabilities")
  }
)


# calculate_mrt_bin_power_f tests -----------------------------------------

size1 <- calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta, 
                          alpha, p_t, gamma, .3, exact=TRUE)

size2 <- calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta, 
                                        alpha, p_t, gamma, .1, exact=TRUE)
test_that(
  "check TQ's sample",
  {
    expect_equal(
      calculate_mrt_bin_power_f(tau_t, f_t, g_t, beta, 
                                     alpha, p_t, gamma, size1),
      .7)
  }
)

test_that(
  "check TQ's sample",
  {
    expect_equal(
      calculate_mrt_bin_power_f(tau_t, f_t, g_t, beta, 
                                alpha, p_t, gamma, size2),
      .9)
  }
)

test_that(
  "check that it works as an 'inverse' of calculate_mrt_bin_samplesize_f",
  {
    expect_equal(
      calculate_mrt_bin_power_f(tau_t, f_t, g_new, beta, 
                                alpha_new, p_t, gamma,
                                calculate_mrt_bin_samplesize_f(
                                  tau_t, f_t, g_new, beta,
                                  alpha_new, p_t, gamma, 1-1/pi, TRUE)),
      1/pi, tolerance=.001)
  }
)

