#' Calculate sample size for binary outcome MRT
#'
#' Returns sample size needed to achieve a specified power for the hypothesis test
#' of marginal excursion effect (see Details) in the context of an MRT with binary outcomes
#' with small sample correction using F-distribution. See the vignette for
#' more details.
#'
#' When the calculator finds out that a sample size less than
#' or equal to 10 is sufficient to attain the desired power, the calculator does
#' not output the exact sample size but produces an error message. This is because 
#' the sample size calculator is based on an asymptotic result, and in
#' this situation the sample size result may not be as accurate.
#' (A small sample correction is built in the calculator, but even with the correction
#' the sample size result may still be inaccurate when it is <= 10.) In general,
#' when the output sample size is small, one might reconsider the following: (1)
#' whether you are correctly or conservatively guessing the average of expected
#' availability, (2) whether the duration of study is too long, (3) whether the
#' treatment effect is overestimated, and (4) whether the power is set too low.
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
#' @param b             Desired Type II error
#' @param exact         Determines if exact n or ceiling will be returned
#'
#' @return              Sample size to achieve desired power.
#' @importFrom          stats uniroot qf pf
#' @export
#'
#' @examples mrt_binary_ss(tau_t_1, f_t_1, g_t_1, 
#'                                          beta_1, alpha_1, p_t_1, 
#'                                          0.05, .2, FALSE)
mrt_binary_ss <- function(avail_pattern,  
                          f_t,             
                          g_t,             
                          beta,            
                          alpha,           
                          p_t,             
                          gamma,          
                          b,
                          exact=FALSE)               
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
    warning("p_t * f_t should lie in span of g_t")
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
  
  if(b <= 0 | b > 1) {
    stop("b, type II error, should be between 0 and 1")
  }
  
  if(!(all(avail_pattern <= 1 & avail_pattern > 0))){
    stop("avail_pattern must have values between 0 and 1")
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
  
  ## towards the end of the main function execution
  if (lincombo_flag) {
    warning("p_t * f_t is not in the linear span of g_t,
          so the sample size result can be inaccurate.\n
          Consider revising g_t.")
  }
  
  p <- length(beta)
  q <- length(alpha)
  
  m_and_sigma <- compute_m_sigma(avail_pattern, f_t, g_t, beta, alpha, p_t)
  m_matrix <- m_and_sigma$m
  sigma_matrix <- m_and_sigma$sigma
  
  suppressWarnings(
    ten_power <- mrt_binary_power(avail_pattern, f_t, g_t, beta, alpha, 
                                  p_t, gamma, 10))
  
  if(1-b <= ten_power) {
    stop(strwrap(paste0("The required sample size is <=10 to attain ", 
                        1-b, 
                        " power for this setting. See help(mrt_binary_ss) for
                        more details"), exdent=1))
  }
  
  # Set up the function that we will ultimately solve to get the sample size
  power_f <- function(n){
    
    right_hand_side <- pf(q=qf(p=(1-gamma), df1=p, df2=n-q-p), 
                          df1=p, 
                          df2=n-q-p, 
                          ncp=compute_ncp(n, beta, m_matrix, sigma_matrix))
    left_hand_side <- b
    return(right_hand_side - left_hand_side)
  }
  
  # find min n to achieve desired power
  sample_size <- uniroot(power_f, lower=p+q+1, upper=1000000)$root
  
  # round up if non-exact size is requested
  if(exact == FALSE){
    sample_size <- ceiling(sample_size)
  }
  
  return(sample_size)
}



