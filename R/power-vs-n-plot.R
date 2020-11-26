
# we want n to vary 
library(ggplot2)



#' Title
#'
#' @param n Sample size
#'
#' @return plots power vs. sample size, based on all other inputs used in the main function
#' @export
#'
#' @examples
plot_power <- function(n) {
  if( n<0 ) {
    stop(" n has to be greater than zero")
  }
  
  n_vec <- c(seq(0,n, by = 1))  
  
  power <- function(...){
    calculate_mrt_bin_power_f(tau_t, f_t, g_t  , beta, alpha, p_t, gamma,...)
    
  }
  plot(n_vec, power(n_vec))
  
}


plot_power(n = 699)

ggplot(data, aes(y = power, x= n_vec)) +  
  geom_point()
coord_cartesian( xlim = c(0, 1000), ylim=c(0, 1) )+ 
  labs(title="Power vs. Sample Size", y="Power", x="Sample size")



##### the other way around   

b_vect <-  c(seq(0.01,0.999991, by = 0.015))
length(b_vect)

sampe_size_vect <- function(...){
  calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta, alpha, p_t, gamma, ...)
}


sampe_size_vect(b_vect)


calculate_mrt_bin_samplesize_f(tau_t, f_t, g_t, beta, alpha, p_t, gamma,0.100)
