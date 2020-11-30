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


# power_summary tests -----------------------------------------------------

test_that(
  "error check for low power",
  {
    expect_error(
      power_summary(tau_t, f_t, g_t, beta, alpha, p_t, 0.05, c(0.000001)))
  }
)

test_that(
  "error check for negative power",
  {
    expect_error(
      power_summary(tau_t, f_t, g_t, beta, alpha, p_t, 0.05, c(-0.000001)),
      "Power should be between 0 and 1")
  }
)

test_that(
  "error check for too large power",
  {
    expect_error(
      power_summary(tau_t, f_t, g_t, beta, alpha, p_t, 0.05, c(100, 100000)),
      "Power should be between 0 and 1")
  }
)
