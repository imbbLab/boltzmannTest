test_that("matrix_hasFullRowRank: `G must be a numeric matrix", {
  ## G is character
  G <- "A"
  expect_error(
    matrix_hasFullRowRank(G),
    regexp = "`G must be a numeric matrix"
  )
  ## G is a list
  G <- list(1)
  expect_error(
    matrix_hasFullRowRank(G),
    regexp = "`G must be a numeric matrix"
  )
  ## G is a data.frame
  G <- data.frame(1)
  expect_error(
    matrix_hasFullRowRank(G),
    regexp = "`G must be a numeric matrix"
  )
})

test_that("eta_isFeasible: `G must be a numeric matrix",{
  ## G is character
  eta <- 1
  G <- "A"
  expect_error(
    eta_isFeasible(G, eta),
    regexp = "`G must be a numeric matrix"
  )
  ## G is a list
  G <- list(1)
  expect_error(
    eta_isFeasible(G, eta),
    regexp = "`G must be a numeric matrix"
  )
  ## G is a data.frame
  G <- data.frame(1)
  expect_error(
    eta_isFeasible(G, eta),
    regexp = "`G must be a numeric matrix"
  )
})

test_that("eta_isFeasible: `eta` must be a numeric vector",{
  ## eta is character
  eta <- "A"
  G <- rbind(
    norm = 1,
    x = c(-1, 1, 0)
  )
  expect_error(
    eta_isFeasible(G, eta),
    regexp = "`eta` must be a numeric vector"
  )
  ## eta is a list
  eta <- list(1)
  expect_error(
    eta_isFeasible(G, eta),
    regexp = "`eta` must be a numeric vector"
  )
  ## eta is a data.frame
  eta <- data.frame(1)
  expect_error(
    eta_isFeasible(G, eta),
    regexp = "`eta` must be a numeric vector"
  )
})

test_that("eta_isFeasible: `eta` must have one value per row of `G`",{
  ## eta has less entries than G rows
  eta <- 1
  G <- rbind(
    norm = 1,
    x = c(-1, 1, 0)
  )
  expect_error(
    eta_isFeasible(G, eta),
    regexp = "`eta` must have one value per row of `G`"
  )
  ## eta has more entries than G rows
  eta <- c(1,2,3)
  expect_error(
    eta_isFeasible(G, eta),
    regexp = "`eta` must have one value per row of `G`"
  )

})

test_that("which_approx_equal: `a` is not a numeric vector",{
  ## a is a character
  a <- "A"
  b <- rep(1/10, 10)
  expect_error(
    which_approx_equal(a, b),
    regexp = "`a` is not a numeric vector"
  )
  ## a is a list
  a <- list(1)

  expect_error(
    which_approx_equal(a, b),
    regexp = "`a` is not a numeric vector"
  )
  ## a is a data.frame
  a <- data.frame(1)
  expect_error(
    which_approx_equal(a, b),
    regexp = "`a` is not a numeric vector"
  )
})
test_that("which_approx_equal: `b` is not a numeric vector",{
  ## b is a character
  a <- rep(1/10, 10)
  b <- "A"
  expect_error(
    which_approx_equal(a, b),
    regexp = "`b` is not a numeric vector"
  )
  ## b is a list
  b <- list(1)

  expect_error(
    which_approx_equal(a, b),
    regexp = "`b` is not a numeric vector"
  )
  ## b is a data.frame
  b <- data.frame(1)
  expect_error(
    which_approx_equal(a, b),
    regexp = "`b` is not a numeric vector"
  )
})
