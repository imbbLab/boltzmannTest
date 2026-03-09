test_that("`p` and `q` must have the same length",{

  p <- rep(1/10, 10)
  q <- rep(1/5, 5)

  expect_error(
    iDivergence(p, q),
    regexp = "`p` and `q` must have the same length"
  )

})
test_that("`p` must be a numeric vector",{
  p <- rep("A", 10)
  q <- rep(1/10, 10)

  expect_error(
    iDivergence(p, q),
    regexp = "`p` must be a numeric vector"
  )
})
test_that("`q` must be a numeric vector",{
  p <- rep(1/10, 10)
  q <- rep("A", 10)

  expect_error(
    iDivergence(p, q),
    regexp = "`q` must be a numeric vector"
  )
})

test_that("`p` must be non-negative",{
  p <- rep(1/10, 10)
  q <- rep(1/10, 10)

  p[1] = -0.1
  expect_error(
    iDivergence(p, q),
    regexp = "`p` must be non-negative"
  )
})

test_that("`q` must be positive",{
  p <- rep(1/10, 10)
  q <- rep(1/10, 10)

  ## one q is zero
  q[1] = 0
  expect_error(
    iDivergence(p, q),
    regexp = "`q` must be positive"
  )
  ## one q is negative
  q[1] = -1
  expect_error(
    iDivergence(p, q),
    regexp = "`q` must be positive"
  )
})

test_that("`p` should sum to 1. Normalizing to 1",{
  p <- rep(1/10, 10)
  q <- rep(1/10, 10)

  p[1] = 0.2
  expect_warning(
    iDivergence(p, q),
    regexp = "`p` should sum to 1. Normalizing to 1"
  )

})

test_that("`q` should sum to 1. Normalizing to 1",{
  p <- rep(1/10, 10)
  q <- rep(1/10, 10)

  q[1] = 0.2
  expect_warning(
    iDivergence(p, q),
    regexp = "`q` should sum to 1. Normalizing to 1"
  )

})
