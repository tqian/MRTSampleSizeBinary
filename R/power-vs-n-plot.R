
#' Returns a plot of power vs sample size in the context of MRT
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
#' @param n             Number of participants, sample size
#'
#' @return              Plot of power and sample size
#' @export
#' @import ggplot2
#'
#' @examples plot_power_vs_samplesize(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1, p_t_1, gamma_1, 300)
plot_power_vs_samplesize <- function(avail_pattern,
                                     f_t,
                                     g_t,
                                     beta,
                                     alpha,
                                     p_t,
                                     gamma,
                                     n)
{
  
  min_n <- length(beta) + length(alpha)
  n_vec <- c(seq(min_n +1, n + 200, by = 1))
  l <- length(n_vec)
  power_vec <- c(rep(NA, l))
  
  for (n_i in n_vec){
    indx <- which(n_vec == n_i)
   power_vec[indx] <- calculate_mrt_bin_power_f(avail_pattern, f_t, g_t, beta, alpha, p_t, gamma,n_i)
 
  }
  #power_vec
  #n_vec
  power_n_data <- data.frame(power_vec, n_vec)
  #power_n_data
  # length(power_vec)
  
  ggplot(power_n_data)+
    geom_line(aes(x = n_vec, y = power_vec))+
    ggtitle("Power vs. Sample Size") +
    xlab("Sample size") + ylab("Power")
  
  
  
  
  # plot(n_vec, 
  #      power_vec,
  #      main = "Power and Sample Size", 
  #      xlab = "Sample Size", 
  #      ylab = "Power", 
  #      ylim = c(0 , 1))
}

