#' Computes the non-centrality parameter for a an F distributed random variable
#' in the context of a binary MRT. This is primarily a helper function for
#' calculate_mrt_bin_power_f and calculate_mrt_bin_samplesize_f.
#'
#' @param x            Sample size
#' @param beta         Marginal excursion effect, assumed dimension p
#' @param m_matrix     "Bread" of sandwich estimator for variance
#' @param sigma_matrix "Meat" of sandwich estimator for variance
#'
#' @return Returns non-centrality parameter for an F distributed random
#'   variable.
#' @export
#'
#' @examples
compute_ncp <- function(x, beta, m_matrix, sigma_matrix){
  return(as.numeric(x * t(beta) %*%
               solve(solve(m_matrix) %*%
                       sigma_matrix %*%
                       t(solve(m_matrix))) %*% 
               beta))
}
