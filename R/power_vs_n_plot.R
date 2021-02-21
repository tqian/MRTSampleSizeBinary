#' Returns a plot of power vs sample size in the context of a binary outcome
#' MRT. See the vignette for more details.
#'
#' @param avail_pattern A vector of length T that is the average availability at
#'   each time point
#' @param f_t           Defines marginal excursion effect MEE(t) under
#'   alternative together with beta. Assumed to be matrix of size T*p.
#' @param g_t           Defines success probability null curve together with
#'   alpha. Assumed to be matrix of size T*q.
#' @param beta          Length p vector that defines marginal excursion effect
#'   MEE(t) under alternative together with f_t.
#' @param alpha         Length q vector that defines success probability null
#'   curve together with g_t.
#' @param p_t           Length T vector of Randomization probabilities at each
#'   time point.
#' @param gamma         Desired Type I error
#' @param min_n         Minimum of range of sample sizes to plot. Should be
#'   greater than the sum of the dimensions of alpha and beta.
#' @param max_n         Maximum of range of sample sizes to plot. Should be
#'   greater than min_n.
#'
#' @return              Plot of power and sample size
#' @export
#' @import ggplot2
#'
#' @examples            power_vs_n_plot(tau_t_1, f_t_1, g_t_1, beta_1, alpha_1,
#'                         p_t_1, 0.05, 15, 700)
power_vs_n_plot <- function(avail_pattern,
                            f_t,
                            g_t,
                            beta,
                            alpha,
                            p_t,
                            gamma,
                            min_n=max(min_samp(alpha, beta), 11),
                            max_n=max_samp(min_n)
)
{
  
  if(min_n >= max_n) {
    stop("max_n should be greater than min_n")
  }
  
  if(min_n < length(beta) + length(alpha)){
    stop(strwrap("min_n is too small. min_n must be greater than the sum of 
                 the dimensions of alpha and beta", exdent=1))
  }
  
  #min_n <- length(beta) + length(alpha)
  n_vec <- c(seq(min_n, max_n, by = 1))
  l <- length(n_vec)
  power_vec <- c(rep(NA, l))
  
  for (n_i in n_vec){
    indx <- which(n_vec == n_i)
    power_vec[indx] <- mrt_binary_power(avail_pattern, f_t, g_t, 
                                        beta, alpha, p_t, gamma,n_i)
    
  }
  
  power_n_data <- data.frame(power_vec, n_vec)
  
  
  ggplot(power_n_data)+
    geom_line(aes(y = power_vec, x = n_vec), color = "deepskyblue3")+
    ggtitle("Power vs. Sample Size") +
    xlab("Sample Size") + ylab("Power")+
    ylim(0,1)+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14))
  
  
  
}

