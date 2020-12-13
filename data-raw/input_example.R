## code to prepare `input_example` dataset goes here


set.seed(1)

total_dp <- 10 # total number of decision points
# p and q are not directly used in the function; they are defined when we define alpha_1 and beta_1.
# p <- 2
# q <- 2

tau_t_1 <- rep(0.8, total_dp) # expected availability E(I_t) for t = 1,...,total_dp

p_t_1 <- rep(0.4, total_dp) # randomization probability over time
gamma_1 <- 0.05 # type I error
b_1 <- 0.2 # type II error; power = 1 - b

### specify g_t and alpha ###
g_t_1 <- cbind(rep(1, total_dp), 1:total_dp)  # g_t = (1, t)
alpha_1 <- c(-0.2, -0.1)

# checking that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 0) should always be between 0 and 1.
mu0_t_1 <- exp(g_t_1 %*% as.matrix(alpha_1, ncol = 1)) # E(Y_{t+1} = 1 | I_t = 1, A_t = 0) for t = 1,...,total_dp
mu0_t_1 # look at the mu0_t values
stopifnot(all(mu0_t_1 < 1) & all(mu0_t_1 > 0))

### specify f_t and beta ###
f_t_1 <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1, t)
beta_1 <- c(0.15, - 0.01)

# checking that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 1) should always be between 0 and 1.
mee_t <- f_t_1 %*% as.matrix(beta_1, ncol = 1) # MEE(t) for t = 1,...,total_dp
mu1_t_1 <- mu0_t_1 * exp(mee_t) # E(Y_{t+1} = 1 | I_t = 1, A_t = 1) for t = 1,...,total_dp
mu1_t_1 # look at the mu1_t values
stopifnot(all(mu1_t_1 < 1) & all(mu1_t_1 > 0))



sigma_matrix_1 <- matrix(c(0.4166547, 2.371095, 2.371095, 16.608106), nrow = 2, byrow = TRUE)
m_matrix_1 <- matrix(c(0.9846607, 4.585791, 4.5857913, 29.055175), nrow = 2, byrow = TRUE)



usethis::use_data(tau_t_1, overwrite = TRUE)
usethis::use_data(f_t_1, overwrite = TRUE)
usethis::use_data(g_t_1, overwrite = TRUE)
usethis::use_data(p_t_1, overwrite = TRUE)
usethis::use_data(beta_1, overwrite = TRUE)
usethis::use_data(alpha_1, overwrite = TRUE)
usethis::use_data(m_matrix_1, overwrite = TRUE)
usethis::use_data(sigma_matrix_1, overwrite = TRUE)


