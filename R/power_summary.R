#' Reports the power of a binary outcome MRT at a range of power levels.
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
#' @param power_levels  Vector of powers to find sample size for.
#'
#' @return              Table containing needed sample size to achieve some
#'   user-specified power values.
#' @export
#' @importFrom          knitr kable
#'
#' @examples            power_summary(tau_t_1, f_t_1, g_t_1, 
#' beta_1, alpha_1, p_t_1, gamma_1)
power_summary <- function(avail_pattern,  
                          f_t,             
                          g_t,             
                          beta,            
                          alpha,           
                          p_t,             
                          gamma,
                          power_levels = seq(from=0.6, to=0.95, by=0.05)) {
  
  if(!(all(power_levels > 0.5 & power_levels < 1))){
    stop("Power should be between 0.5 and 1")
  }
  
  power_size <- cbind(power=power_levels, sample_size=0)
  
  for(r in 1:length(power_levels)){
    power_size[r,2] <-  calculate_mrt_bin_samplesize_f(avail_pattern, f_t, g_t, 
                                                        beta, alpha, p_t, gamma, 
                                                        1-power_levels[r])
  }
  
  power_size
}
