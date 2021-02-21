#' Computes "M" and "Sigma" matrices for the sandwich estimator of
#' variance-covariance matrix. 
#' 
#' A helper function for
#' mrt_binary_power() and mrt_binary_ss().
#'
#' @param avail_pattern A vector of length T that is the average availability at
#'   each time point
#' @param f_t           Defines marginal excursion effect MEE(t) under
#'   alternative together with beta. Assumed to be matrix of size T*p.
#' @param g_t           Defines success probability null curve together with
#'   alpha. Assumed to be matrix of size T*q.
#' @param beta          Length p vector that defines marginal excursion effect
#'   MEE(t) under alternative together with f_t.
#' @param alpha         Length q vector that defines success probability null
#'   curve together with g_t.
#' @param p_t           Length T vector of randomization probabilities at each time point
#' @return              List containing two matrices. The first is the M matrix
#'   and the second is the Sigma matrix.
#' @export
#'
#' @examples            compute_m_sigma(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1,
#'                                       p_t_1)
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
