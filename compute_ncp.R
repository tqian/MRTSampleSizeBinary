# helper function for calculate_mrt_bin_samplesize_f
# computes the non-centrality parameter for the F distribution when solving for
# sample size
compute_ncp <- function(x, beta, m_matrix, sigma_matrix){
  as.numeric(x * t(beta) %*%
               solve(solve(m_matrix) %*%
                       sigma_matrix %*%
                       t(solve(m_matrix))) %*% 
               beta)
}