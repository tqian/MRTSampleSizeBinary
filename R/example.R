#'
#' @format A list of generated values to use as an example for MRT sample size and power calculations
#' \describe{
#'   \item {tau_t_1} {A vector of length T that is the average availability at
#'   each time point}
#'   \item {f_t_1} {Defines marginal excursion effect MEE(t) under
#'   alternative together with beta
#'   \item {g_t_1} {Defines success probability null curve together with}
#'   alpha
#'   \item{beta_1} {Defines marginal excursion effect MEE(t) under}
#'   alternative together with g_t
#'   \item{alpha_1} {Defines success probability null curve together with f_t_1}
#'   \item{p_t_1} {Randomization probability at each time point}
#'   \item{gamma_1} {Desired Type I error}
#'   \item{m_matrix_1} { "Bread" of sandwich estimator for variance}
#'   \item{sigma_matrix_1} {"Meat" of the sandwich estimator}
#' }
#' @format 
#'   "total_dp",
#'   "tau_t_1",
#'   "f_t_1",
#'   "g_t_1",
#'   "p_t_1",
#'   "beta_1", 
#'   "alpha_1", 
#'   "b_1", 
#'   "gamma_1",
#'   "m_matrix_1",
#'   "sigma_matrix_1",
#'   "x"
#'   
#'   
#'   
#' 
#' 
#' 
#' 
#'  

#' @format numeric value
#' \describe{
#' 10
#' }
"total_dp"


#' @format vector of length T
#' \describe{A vector of length T that is the aver availability at each time point.}
"tau_t_1"

#' @format a 2 by 10 matrix defining a linear effect
#' \describe{
#' vector for the null
#' }
"f_t_1"
   
   
#' @format a 2 by 10 matrix defining a linear effect
#' \describe{
#' vector for alternative
#' }
"g_t_1"
   
#' @format a length T vector of probability for each time
#' \describe{
#' vectro of probs
#' }

"p_t_1"

#' @format beta
#' \describe{
#' adsf}
"beta_1"

#' @format beta
#' \describe{
#' adsf}
"alpha_1"


#' @format beta
#' \describe{
#' adsf}
"gamma_1"

#' @format beta
#' \describe{
#' adsf}
"b_1"
