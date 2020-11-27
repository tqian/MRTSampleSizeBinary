
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
#' @return              Plot of power and sample sice
#' @export
#'
#' @examples
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
  
  power_vec
  length(power_vec)
  plot(n_vec, power_vec)
}
