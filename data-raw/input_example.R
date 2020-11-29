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
mu0_t_1 <- exp(g_t_1 %*% alpha_1) # E(Y_{t+1} = 1 | I_t = 1, A_t = 0) for t = 1,...,total_dp
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


############ to creat an example for compute_ncp, we use the compute m-sigma function 

compute_m_sigma <- function(avail_pattern,  
                            f_t,             
                            g_t,            
                            beta,            
                            alpha,           
                            p_t              
) {
  p <- length(beta)
  q <- length(alpha)
  
  ## The M and Sigma matrices (needed to compute lambda)
  m_matrix <- matrix(data=0, nrow=length(beta), ncol=length(beta))
  sigma_matrix <- matrix(data=0, nrow=length(beta), ncol=length(beta))
  
  # E(Y_{t+1} = 1 | I_t = 1, A_t = 0) for t = 1,...,total_dp
  mu0_t <- exp(g_t %*% alpha) 
  
  if(!(all(mu0_t < 1) & all(mu0_t > 0))) {
    stop("g_t and alpha values led to invalid probabilities")
  }
  
  # check that probability E(Y_{t+1} = 1 | I_t = 1, A_t = 1) is 0 and 1.
  mee_t <- f_t %*% beta # MEE(t) for t = 1,...,total_dp
  # E(Y_{t+1} = 1 | I_t = 1, A_t = 1) for t = 1,...,total_dp
  mu1_t <- mu0_t * exp(mee_t) 
  if(!(all(mu1_t < 1) & all(mu1_t > 0))){
    stop("f_t and beta values led to invalid probabilities.")
  }
  
  # For each decision point 
  # (T = total # of decision points is taken as the length of p_t, for now)
  for (i in 1:length(p_t)){
    
    # breaking down steps to identify bug and improve robustness of code
    this_f_t <- as.numeric(f_t[i, ])
    this_g_t <- as.numeric(g_t[i, ])
    
    if(length(this_f_t) != p){
      stop("Incorrect dimensions for f_t.")
    }
    
    if(length(this_g_t) != q){
      stop("Incorrect dimsions for g_t.")
    }
    
    # both are vectors, so use this way to calculate inner product
    this_f_t_times_beta <- sum(this_f_t * beta) 
    this_g_t_times_alpha <- sum(this_g_t * alpha)
    
    this_f_t_f_t <- outer(this_f_t, this_f_t) # this is f_t %*% f_t^T
    
    this_m <- as.numeric(avail_pattern[i] * 
                           (exp(p_t[i] * this_f_t_times_beta)) * 
                           exp(this_g_t_times_alpha) *
                           (1-p_t[i]) *
                           p_t[i]) * this_f_t_f_t
    
    # added a check of dimension
    if(!("matrix" %in% class(this_m) & all(dim(this_m) == c(p, p)))){
      stop("Incorrect dimensions for M matrix, or improper type.")
    } 
    
    # A running sum so that we end up with each entry being the sum of that
    # entry across all time points
    m_matrix <- m_matrix + this_m 
    
    this_sigma <- as.numeric(avail_pattern[i] * 
                               (exp(2 * p_t[i] * this_f_t_times_beta)) *
                               exp(this_g_t_times_alpha) * 
                               (1-p_t[i]) * 
                               p_t[i] * 
                               ((1-p_t[i]) * 
                                  exp(-1 * this_f_t_times_beta) +
                                  p_t[i] - exp(this_g_t_times_alpha))) * 
      this_f_t_f_t
    
    
    sigma_matrix <- sigma_matrix + this_sigma
    if(!("matrix" %in% class(this_sigma) & 
         all(dim(this_sigma) == c(p, p)))){
      stop("Incorrect dimensions for M matrix, or improper type.")
    }
  }
  
  return(list(m = m_matrix, sigma = sigma_matrix))
}


m_matrix_1 <- compute_m_sigma(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1, p_t_1)$m

sigma_matrix_1 <-compute_m_sigma(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1, p_t_1)$sigma
x <- calculate_mrt_bin_samplesize_f(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1, p_t_1, gamma_1, b_1)

usethis::use_data(x, overwrite = TRUE)
usethis::use_data(total_dp, overwrite = TRUE)
usethis::use_data(tau_t_1, overwrite = TRUE)
usethis::use_data(f_t_1, overwrite = TRUE)
usethis::use_data(g_t_1, overwrite = TRUE)
usethis::use_data(p_t_1, overwrite = TRUE)
usethis::use_data(beta_1, overwrite = TRUE)
usethis::use_data(alpha_1, overwrite = TRUE)
usethis::use_data(b_1, overwrite = TRUE)
usethis::use_data(gamma_1, overwrite = TRUE)
usethis::use_data(m_matrix_1, overwrite = TRUE)
usethis::use_data(sigma_matrix_1, overwrite = TRUE)


