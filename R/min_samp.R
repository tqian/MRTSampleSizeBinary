#' Compute minimum sample size.
#' 
#' Returns a default minimum sample size to start power_vs_n_plot() at.
#'
#' @param alph Vector to describe the MEE under the alternative.
#' @param bet  Vector to describe the MEE under the null.
#'
#' @return A default minimum sample size to start power_vs_n_plot() at.
#' @export
#'
#' @examples min_samp(alpha_1, beta_1)
min_samp <- function(alph, bet){
  return(length(alph) + length(bet) + 1)
}
