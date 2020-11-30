test_that(
  "check true",
  {
    expect_equal(
      is_full_column_rank(diag(9)),
      TRUE)
  }
)

test_that(
  "check false",
  {
    expect_equal(
      is_full_column_rank(diag(rep(0, times=4))),
      FALSE)
  }
)

test_that(
  "check false",
  {
    expect_equal(
      is_full_column_rank(diag(c(rep(0, times=4), 5))),
      FALSE)
  }
)

t_mat <- cbind(rep(1, 4), 1:4)

test_that(
  "check false",
  {
    expect_equal(
      is_full_column_rank(t_mat),
      TRUE)
  }
)
