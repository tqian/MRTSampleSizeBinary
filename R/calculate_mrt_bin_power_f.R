#' Returns power of a test given a prespecified sample size in the context of a
#' binary MRT.
#'
#' @param avail_pattern A vector of length T that is the average availability at
#'   each time point
#' @param f_t           Defines marginal excursion effect MEE(t) under
#'   alternative together with beta. Assumed to be matrix of size T*p.
#' @param g_t           Defines success probability null curve together with
#'   alpha. Assumed to be matrix of size T*q.
#' @param beta          Length p vector that defines marginal excursion effect
#'   MEE(t) under alternative together with g_t.
#' @param alpha         Length q vector that defines success probability null
#'   curve together with f_t.
#' @param p_t           Length T vector of Randomization probabilities at each
#'   time point.
#' @param gamma         Desired Type I error
#' @param n             Number of participants
#'
#' @return              Power of the test for fixed null/alternative and sample
#'   size.
#' @importFrom          stats qf pf
#'
#' @export
#'
#' @examples calculate_mrt_bin_power_f(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1, p_t_1, gamma_1, 300)
calculate_mrt_bin_power_f <- function(avail_pattern,  
                                           f_t,             
                                           g_t,             
                                           beta,            
                                           alpha,           
                                           p_t,             
                                           gamma,          
                                           n)               
{
  pts <- length(avail_pattern)
  
  dim_vec <- c(dim(f_t)[1], dim(g_t)[1], length(p_t))
  
  if(!all(dim_vec == pts)){
    stop("All arguments must agree on number of time points.")
  }
  
  if(dim(f_t)[2] > dim(g_t)[2]){
    warning("f should lie in span of g")
  }
  
  if(dim(f_t)[2] != length(beta)) {
    stop("Dimensions of f_t and beta do not agree.")
  }
  
  if(dim(g_t)[2] != length(alpha)) {
    stop("Dimensions of g_t and alpha do not agree.")
  }
  
  if(gamma < 0 | gamma > 1) {
    stop("gamma, type I error, should be between 0 and 1")
  }
  
  if(!(all(avail_pattern <= 1 & avail_pattern > 0))){
    stop("avail_pattern must have values between 0 and 1")
  }
  

  p <- length(beta)
  q <- length(alpha)
  
  if(n < (p+q)) {
    stop("n is too small")
  }
  
  if(n %% 1 != 0){
    warning("n should be an integer")
  }
  
  
  m_and_sigma <- compute_m_sigma(avail_pattern, f_t, g_t, beta, alpha, p_t)
  m_matrix <- m_and_sigma$m
  sigma_matrix <- m_and_sigma$sigma
  
  return(1-pf(q=qf(p=(1-gamma), df1=p, df2=n-q-p), 
     df1=p, 
     df2=n-q-p, 
     ncp=compute_ncp(n, beta, m_matrix, sigma_matrix)))
}
