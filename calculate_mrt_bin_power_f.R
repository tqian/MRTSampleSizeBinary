#' Returns power of a test given a prespecified sample size in the context of a
#' binary MRT.
#'
#' @param avail_pattern A vector of length T that is the average availability at
#'   each time point
#' @param f_t           Defines marginal excursion effect MEE(t) under
#'   alternative together with beta
#' @param g_t           Defines success probability null curve together with
#'   alpha
#' @param beta          Defines marginal excursion effect MEE(t) under
#'   alternative together with g_t
#' @param alpha         Defines success probability null curve together with f_t
#' @param p_t           Randomization probability at each time point
#' @param gamma         Desired Type I error
#' @param n             Number of participants
#'
#' @return              Power of the test.
#' @export
#'
#' @examples
calculate_mrt_bin_power_f <- function(avail_pattern,  
                                           f_t,             
                                           g_t,             
                                           beta,            
                                           alpha,           
                                           p_t,             
                                           gamma,          
                                           n)               
{
  p <- length(beta)
  q <- length(alpha)
  
  m_and_sigma <- compute_m_sigma(avail_pattern, f_t, g_t, beta, alpha, p_t)
  m_matrix <- m_and_sigma$m
  sigma_matrix <- m_and_sigma$sigma
  
  return(1-pf(q=qf(p=(1-gamma), df1=p, df2=n-q-p), 
     df1=p, 
     df2=n-q-p, 
     ncp=compute_ncp(n, beta, m_matrix, sigma_matrix)))
}
