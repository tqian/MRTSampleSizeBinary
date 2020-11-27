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
  if(det(m_matrix) == 0){
    stop("m_matrix must be nonsingular")
  }
  
  if(length(beta) != dim(m_matrix)[1]){
    stop("Dimensions of beta and m_matrix do not agree")
  }
  
  if(dim(m_matrix)[1] != dim(sigma_matrix)[2] | 
     dim(m_matrix)[2] != dim(sigma_matrix)[1]){
    stop("Dimensions of m_matrix and sigma_matrix do not agree")
  }
  
  if(det(solve(m_matrix) %*%
         sigma_matrix %*%
         t(solve(m_matrix))) == 0){
    stop("sigma_matrix must be nonsingular")
  }
  
  return(as.numeric(x * t(beta) %*%
               solve(solve(m_matrix) %*%
                       sigma_matrix %*%
                       t(solve(m_matrix))) %*% 
               beta))
}