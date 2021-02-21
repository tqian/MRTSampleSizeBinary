#' Calculate power for binary outcome MRT
#'
#' Returns power of the hypothesis test of marginal excursion effect (see Details) given 
#' a specified sample size in the context of an MRT with binary outcomes
#' with small sample correction using F-distribution. See the vignette for
#' more details.
#'
#' @param avail_pattern A vector of length m that is the average availability at
#'   each time point
#' @param f_t           Defines marginal excursion effect MEE(t) under
#'   alternative together with beta. Assumed to be matrix of size m*p.
#' @param g_t           Defines success probability null curve together with
#'   alpha. Assumed to be matrix of size m*q.
#' @param beta          Length p vector that defines marginal excursion effect
#'   MEE(t) under alternative together with f_t.
#' @param alpha         Length q vector that defines success probability null
#'   curve together with g_t.
#' @param p_t           Length m vector of Randomization probabilities at each
#'   time point.
#' @param gamma         Desired Type I error
#' @param n             Sample size
#'
#' @importFrom          stats qf pf
#'
#' @return              Power of the test for fixed null/alternative and sample
#'   size.
#' @export
#'
#' @examples            mrt_binary_power(tau_t_1, f_t_1, g_t_1, beta_1,
#'                                               alpha_1, p_t_1, 0.05, 100)
mrt_binary_power <- function(avail_pattern,  
                             f_t,             
                             g_t,             
                             beta,            
                             alpha,           
                             p_t,             
                             gamma,          
                             n)               
{
  pts <- length(avail_pattern)
  
  if(!(is.matrix(f_t)) | !(is.matrix(g_t))){
    stop("f_t and g_t should be matrices")
  }
  
  dim_vec <- c(dim(f_t)[1], dim(g_t)[1], length(p_t))
  
  if(!all(dim_vec == pts)){
    stop("All arguments must agree on number of time points.")
  }
  
  if(dim(f_t)[2] > dim(g_t)[2]){
    warning("p_t *f_t should lie in span of g_t")
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
  
  if(n <= 10) {
    warning("n <= 10 may result in inaccurate power calculation, because the sample size formula is based on an asymptotic result.")
  }
  
  # check that f_t is of full column rank
  if(!is_full_column_rank(f_t)) {
    stop("f_t has linearly dependent columns.")
  }
  
  # check that g_t is of full column rank
  if(!is_full_column_rank(g_t)) {
    stop("g_t has linearly dependent columns.")
  }
  
  # check that p_t * f_t is in the linear span of g_t
  lincombo_flag <- FALSE
  for (icol in 1:ncol(f_t)) {
    if (is_full_column_rank(cbind(p_t * f_t[, icol], g_t))) {
      # p_t * f_t[, icol] is not in the linear span of g_t
      warning(paste0("p_t * f_t[, ", icol, 
                     "] is not in the linear span of g_t."))
      lincombo_flag <- TRUE
    }
  }
  
  # check that p_t * f_t is in the linear span of g_t
  if (lincombo_flag) {
    warning("p_t * f_t is not in the linear span of g_t,
          so the sample size result can be inaccurate.\n
          Consider revising g_t.")
  }
  
  
  m_and_sigma <- compute_m_sigma(avail_pattern, f_t, g_t, beta, alpha, p_t)
  m_matrix <- m_and_sigma$m
  sigma_matrix <- m_and_sigma$sigma
  
  return(1-pf(q=qf(p=(1-gamma), df1=p, df2=n-q-p), 
              df1=p, 
              df2=n-q-p, 
              ncp=compute_ncp(n, beta, m_matrix, sigma_matrix)))
}
