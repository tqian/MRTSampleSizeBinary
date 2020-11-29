#' Title This function returns a table of selected sample sizes based on power values
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
#'
#' @return              Table of the input values, and sample sizes with some selected power values
#' @export
#' @import dplyr
#'
#' @examples   power_summary(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1, p_t_1, gamma_1)


library(mrtbincalc)
power_summary <- function(avail_pattern,
                          f_t,
                          g_t,
                          beta,
                          alpha,
                          p_t,
                          gamma
)
{
  
  min_n <- length(beta) + length(alpha)
  n_vec <- c(seq(min_n +1, min_n + 100, by = 1))
  l <- length(n_vec)
  power_vec <- c(rep(NA, l))
  
  for (n_i in n_vec){
    indx <- which(n_vec == n_i)
    power_vec[indx] <- calculate_mrt_bin_power_f(avail_pattern, f_t, g_t, beta, alpha, p_t, gamma,n_i)
    
  }
  #power_vec
  #n_vec
  
  
  power_n_data <- data.frame(power_vec = round(power_vec, 3),  n_vec)
  selected_power <- dplyr::filter(power_n_data, power_vec %in% seq(0.6, 0.95, by = 0.05))
  input <- list("Available Pattern" = avail_pattern,
                "f(t)" = f_t,
                "g(t)" = g_t,
                "beta" = beta,
                "alpha" = alpha,
                "Randomization prob" = p_t,
                "Desired Type I error" = gamma)
  #input
  (selected_power)
  
  #power_n_data
  
  
}
power_summary(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1, p_t_1, gamma_1)


(z = dplyr::filter(data.frame(x =c(1,2,3), y= c(4,5,6)), x %in% c(1,10)) )
