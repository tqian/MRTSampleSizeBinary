## code to prepare `input_example` dataset goes here


set.seed(1)

total_dp <- 10 # total number of decision points
# p and q are not directly used in the function, so I commented them out.
# p <- 2
# q <- 2

tau_t_1 <- rep(0.8, total_dp) # expected availability E(I_t) for t = 1,...,total_dp

p_t_1 <- rep(0.4, total_dp) # randomization probability over time
gamma_1 <- 0.05 # type I error
b_1 <- 0.2 # type II error; power = 1 - b

### specify g_t and alpha ###
g_t_1 <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1,g_t)
alpha_1 <- as.matrix(c(-0.2, -0.1), ncol = 1)

# checking that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 0) should always be between 0 and 1.
mu0_t_1 <- exp(g_t %*% alpha) # E(Y_{t+1} = 1 | I_t = 1, A_t = 0) for t = 1,...,total_dp
mu0_t_1 # look at the mu0_t values
stopifnot(all(mu0_t < 1) & all(mu0_t > 0))

### specify f_t and beta ###
f_t_1 <- cbind(rep(1, total_dp), 1:total_dp)  # f_t = (1, t)
beta_1 <- as.matrix(c(0.15, - 0.01), ncol = 1)

# checking that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 1) should always be between 0 and 1.
mee_t <- f_t_1 %*% beta_1 # MEE(t) for t = 1,...,total_dp
mu1_t <- mu0_t * exp(mee_t) # E(Y_{t+1} = 1 | I_t = 1, A_t = 1) for t = 1,...,total_dp
mu1_t # look at the mu1_t values
stopifnot(all(mu1_t < 1) & all(mu1_t > 0))

### saving the inputs in a list bc they have different kinds of data, vectors, matrices and scalars 
input_example <- list(total_dp = total_dp,
                      tau_t_1 = tau_t_1, 
                      f_t_1 = f_t_1, 
                      g_t_1 = g_t_1, 
                      p_t_1 = p_t_1, 
                      beta_1 = beta_1,
                      alpha_1 = alpha_1,
                      b_1 = b_1, 
                      gamma_1 = gamma_1
                      )

usethis::use_data(input_example, overwrite = TRUE)
