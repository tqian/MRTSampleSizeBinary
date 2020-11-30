#' Returns default maximum sample size to end power_vs_n_plot().
#'
#' @param min_samp The starting sample size of the plot.
#'
#' @return A default maximum sample size to end power_vs_n_plot().
#' @export
#'
#' @examples max_samp(100)
max_samp <- function(min_samp){
  return(min_samp + 1000)
}
