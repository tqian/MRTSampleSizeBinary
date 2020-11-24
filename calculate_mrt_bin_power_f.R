
# returns power of test given a fixed sample size and specified 
# null/alternative hypotheses
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
  
  return(pf(q=qf(p=(1-gamma), df1=p, df2=n-q-p), 
     df1=p, 
     df2=n-q-p, 
     ncp=compute_ncp(n, beta, m_matrix, sigma_matrix)))
}
