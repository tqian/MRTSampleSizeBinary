#' Vector that holds the average availability at each time point.
#' @format vector of length T
#' \describe{A vector of length T that is the average availability
#'  at each decision point.}
"tau_t_1"

#' A matrix defining the MEE under the alternative hypothesis.
#' @format a 10 by 2 matrix
#' \describe{
#' In this example it is a log-linear trend.
#' }
"f_t_1"

#' A matrix defining the success probability null curve.
#' @format a 10 by 2 matrix
#' \describe{
#' In this example it is a log-linear trend.
#' }
"g_t_1"


#' A vector of randomization probabilities for each time point.  
#' @format a length T vector
#' \describe{
#' Vector of randomization probabilities. 
#' }
"p_t_1"

#' Vector that defines the MEE under the alternative hypothesis.
#' @format a length 2 vector
#' \describe{
#' The matrix multiplication of this vector with f_t_1 defines the MEE under the
#' alternative hypothesis.
#' }
"beta_1"

#' Vector that defines the success probability null curve.
#' @format a length 2 vector
#' \describe{
#' The matrix multiplication of this vector with g_t_1 defines the MEE under the
#' null hypothesis.
#' }
"alpha_1"


#' An example matrix for "meat" of sandwich estimator of variance.
#' @format A 2 by 2 matrix
#' \describe{
#' Generated from a toy example.
#' }
"sigma_matrix_1"

#' An example matrix for "bread" of sandwich estimator of variance.
#' @format A 2 by 2 matrix
#' \describe{
#' Generated from a toy example.
#' }
"m_matrix_1"


