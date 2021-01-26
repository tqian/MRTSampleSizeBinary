#' Check if a matrix is full column rank.
#' 
#' Used in checking if p_t*f_t is in the linear span of g_t.
#'
#' @param mat A matrix.
#'
#' @return Boolean TRUE/FALSE for if matrix is full column rank.
#' @importFrom Matrix rankMatrix
#' @export
#' @examples   is_full_column_rank(diag(4))
is_full_column_rank <- function(mat) {
  # Full column rank if only if all columns are linearly independent
  # See: https://en.wikipedia.org/wiki/Rank_(linear_algebra)
  return(rankMatrix(mat)[[1]] == ncol(mat))
}
