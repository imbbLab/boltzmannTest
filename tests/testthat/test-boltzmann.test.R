test_that("different number of expectations and rows in G", {
  data("kidneyStones")
  outcomes <- outcomes_tibble(kidneyStones)
  ## coefficient matrix
  G <- with(
    outcomes,
    rbind(
      norm = 1,
      # structural expectation fraction treatment A
      treatment_A = treatment == "A",
      # ambient expectation fraction small stone size given treatment B
      stoneSize_small.treatment_B =
        (stoneSize == "small" & treatment == "B") / 0.5,
      # ambient expectation fraction small stone size given treatment A
      stoneSize_small.treatment_A =
        (stoneSize == "small" & treatment == "A") / 0.5,
      # tested expectation: success rate treatment A versus B
      "success_yes:B_vs_A" = c(0, -2, 0, 2, 0, -2, 0, 2)
    )
  )
  ## hypothesized values of the expectations
  eta <- c(
    norm = 1,
    treatment_A = 0.5,
    # overall frequency of small stones is 0.51
    stoneSize_small.treatment_B = 0.51,
    stoneSize_small.treatment_A = 0.51,
    "success_yes:B_vs_A" = 0.0
  )

  expect_error(
    boltzmann.test(
      outcomes = outcomes,
      G = G[-2, ],
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4)
    ),
    regexp = "the number of generalized moments given by"
  )


})
test_that("testedExpectations out of range", {
  data("kidneyStones")
  outcomes <- outcomes_tibble(kidneyStones)
  ## coefficient matrix
  G <- with(
    outcomes,
    rbind(
      norm = 1,
      # structural expectation fraction treatment A
      treatment_A = treatment == "A",
      # ambient expectation fraction small stone size given treatment B
      stoneSize_small.treatment_B =
        (stoneSize == "small" & treatment == "B") / 0.5,
      # ambient expectation fraction small stone size given treatment A
      stoneSize_small.treatment_A =
        (stoneSize == "small" & treatment == "A") / 0.5,
      # tested expectation: success rate treatment A versus B
      "success_yes:B_vs_A" = c(0, -2, 0, 2, 0, -2, 0, 2)
    )
  )
  ## hypothesized values of the expectations
  eta <- c(
    norm = 1,
    treatment_A = 0.5,
    # overall frequency of small stones is 0.51
    stoneSize_small.treatment_B = 0.51,
    stoneSize_small.treatment_A = 0.51,
    "success_yes:B_vs_A" = 0.0
  )

  expect_error(
    boltzmann.test(
      outcomes = outcomes,
      G = G,
      eta = eta,
      testedExpectations = 6,
      ambientExpectations = c(3,4)
    ),
    regexp = "the indices of the tested moments are out of range"
  )


})
test_that("ambientExpectations out of range", {
  data("kidneyStones")
  outcomes <- outcomes_tibble(kidneyStones)
  ## coefficient matrix
  G <- with(
    outcomes,
    rbind(
      norm = 1,
      # structural expectation fraction treatment A
      treatment_A = treatment == "A",
      # ambient expectation fraction small stone size given treatment B
      stoneSize_small.treatment_B =
        (stoneSize == "small" & treatment == "B") / 0.5,
      # ambient expectation fraction small stone size given treatment A
      stoneSize_small.treatment_A =
        (stoneSize == "small" & treatment == "A") / 0.5,
      # tested expectation: success rate treatment A versus B
      "success_yes:B_vs_A" = c(0, -2, 0, 2, 0, -2, 0, 2)
    )
  )
  ## hypothesized values of the expectations
  eta <- c(
    norm = 1,
    treatment_A = 0.5,
    # overall frequency of small stones is 0.51
    stoneSize_small.treatment_B = 0.51,
    stoneSize_small.treatment_A = 0.51,
    "success_yes:B_vs_A" = 0.0
  )

  expect_error(
    boltzmann.test(
      outcomes = outcomes,
      G = G,
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4,-8)
    ),
    regexp = "the indices of the ambient expectations are out of range"
  )


})
test_that("ambientExpectations overlap with testedExpectations", {
  data("kidneyStones")
  outcomes <- outcomes_tibble(kidneyStones)
  ## coefficient matrix
  G <- with(
    outcomes,
    rbind(
      norm = 1,
      # structural expectation fraction treatment A
      treatment_A = treatment == "A",
      # ambient expectation fraction small stone size given treatment B
      stoneSize_small.treatment_B =
        (stoneSize == "small" & treatment == "B") / 0.5,
      # ambient expectation fraction small stone size given treatment A
      stoneSize_small.treatment_A =
        (stoneSize == "small" & treatment == "A") / 0.5,
      # tested expectation: success rate treatment A versus B
      "success_yes:B_vs_A" = c(0, -2, 0, 2, 0, -2, 0, 2)
    )
  )
  ## hypothesized values of the expectations
  eta <- c(
    norm = 1,
    treatment_A = 0.5,
    # overall frequency of small stones is 0.51
    stoneSize_small.treatment_B = 0.51,
    stoneSize_small.treatment_A = 0.51,
    "success_yes:B_vs_A" = 0.0
  )

  expect_error(
    boltzmann.test(
      outcomes = outcomes,
      G = G,
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4,5)
    ),
    regexp = "the ambient expectations must not include the tested expectations"
  )


})

test_that("hypothesized expectations contain NA or non-finite number", {
  data("kidneyStones")
  outcomes <- outcomes_tibble(kidneyStones)
  ## coefficient matrix
  G <- with(
    outcomes,
    rbind(
      norm = 1,
      # structural expectation fraction treatment A
      treatment_A = treatment == "A",
      # ambient expectation fraction small stone size given treatment B
      stoneSize_small.treatment_B =
        (stoneSize == "small" & treatment == "B") / 0.5,
      # ambient expectation fraction small stone size given treatment A
      stoneSize_small.treatment_A =
        (stoneSize == "small" & treatment == "A") / 0.5,
      # tested expectation: success rate treatment A versus B
      "success_yes:B_vs_A" = c(0, -2, 0, 2, 0, -2, 0, 2)
    )
  )
  ## hypothesized values of the expectations
  eta <- c(
    norm = 1,
    treatment_A = 0.5,
    # overall frequency of small stones is 0.51
    stoneSize_small.treatment_B = NA,
    stoneSize_small.treatment_A = 0.51,
    "success_yes:B_vs_A" = 0.0
  )

  expect_error(
    boltzmann.test(
      outcomes = outcomes,
      G = G,
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4)
    ),
    regexp = "`eta` contains non-finite or missing values"
  )
  ## hypothesized values of the expectations
  eta <- c(
    norm = 1,
    treatment_A = 0.5,
    # overall frequency of small stones is 0.51
    stoneSize_small.treatment_B = 0.51,
    stoneSize_small.treatment_A = 0.51,
    "success_yes:B_vs_A" = Inf
  )
  expect_error(
    boltzmann.test(
      outcomes = outcomes,
      G = G,
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4)
    ),
    regexp = "`eta` contains non-finite or missing values"
  )
  ## hypothesized values of the expectations
  eta <- c(
    norm = 1,
    treatment_A = NaN,
    # overall frequency of small stones is 0.51
    stoneSize_small.treatment_B = 0.51,
    stoneSize_small.treatment_A = 0.51,
    "success_yes:B_vs_A" = 0.0
  )
  expect_error(
    boltzmann.test(
      outcomes = outcomes,
      G = G,
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4)
    ),
    regexp = "`eta` contains non-finite or missing values"
  )
})

test_that("nested hypothesis testing", {
  data("kidneyStones")
  outcomes <- outcomes_tibble(kidneyStones)
  ## coefficient matrix
  G <- with(
    outcomes,
    rbind(
      norm = 1,
      # structural expectation fraction treatment A
      treatment_A = treatment == "A",
      # ambient expectation fraction small stone size given treatment B
      stoneSize_small.treatment_B =
        (stoneSize == "small" & treatment == "B") / 0.5,
      # ambient expectation fraction small stone size given treatment A
      stoneSize_small.treatment_A =
        (stoneSize == "small" & treatment == "A") / 0.5,
      # tested expectation: success rate treatment A versus B
      "success_yes:B_vs_A" = c(0, -2, 0, 2, 0, -2, 0, 2)
    )
  )
  ## hypothesized values of the expectations
  eta <- c(
    norm = 1,
    treatment_A = 0.5,
    # overall frequency of small stones is 0.51
    stoneSize_small.treatment_B = 0.51,
    stoneSize_small.treatment_A = 0.51,
    "success_yes:B_vs_A" = 0.0
  )

  expect_no_error(
    bt <- boltzmann.test(
      outcomes = outcomes,
      G = G,
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4)
    )
  )
  expect_equal(round(bt$pValue, 4), 0.0592)

})

test_that("hypothesis distribution is not feasible", {
  outcomes <- outcomes_tibble(
    data.frame(
      x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
    )
  )

  G <- rbind(
    norm = 1,
    x = outcomes$x
  )
  eta <- c(
    norm = 1,
    x = 2
  )
  expect_warning(
    bt <- boltzmann.test(outcomes, G, eta, testedExpectations = 2),
    regexp = "`eta` is not a feasible"
  )
  expect_true(is.na(bt$statistic))
  expect_true(is.na(bt$iDivergence))
  expect_true(is.na(bt$pValue))
  expect_true(all(is.na(bt$hypothesisDistribution)))

})

test_that("hypothesis distribution is has zero probabilities", {
  outcomes <- outcomes_tibble(
    data.frame(
      x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
    )
  )

  G <- rbind(
    norm = 1,
    x = outcomes$x
  )
  eta <- c(
    norm = 1,
    x = 1
  )
  expect_warning(
    bt <- boltzmann.test(outcomes, G, eta, testedExpectations = 2),
    regexp = "hypothesis conditions led to zero probabilities"
  )
  expect_true(is.infinite(bt$statistic))
  expect_true(is.infinite(bt$iDivergence))
  expect_true(bt$pValue == 0)
  expect_true(all(any(bt$hypothesisDistribution < .Machine$double.eps)))

})

test_that("checks for mu",{

  expect_error(
    boltzmann.test(x = c(rep(-1, 4), rep(0, 4), rep(1, 2)), mu = NULL),
    regexp = "`mu` must be a single number"
  )
  expect_error(
    boltzmann.test(x = c(rep(-1, 4), rep(0, 4), rep(1, 2)), mu = c(1,5)),
    regexp = "`mu` must be a single number"
  )
  expect_error(
    boltzmann.test(x = c(rep(-1, 4), rep(0, 4), rep(1, 2)), mu = "A"),
    regexp = "`mu` must be a single number"
  )
})

test_that("checks for numeric y", {
  expect_error(
    boltzmann.test(x = c(rep(-1, 4), rep(0, 4), rep(1, 2)), y = rep("A", 10)),
    regexp = "`y` must a be a numeric vector"
  )
})
test_that("check for same length for x and y in paired test",{
  expect_error(
    boltzmann.test(x = rep(0, 10), y = rep(1, 12), paired = TRUE),
    regexp = "must have the same length"
  )
})
test_that("check for same length for x and no y in paired test",{
  expect_error(
    boltzmann.test(x = rep(0, 10), paired = TRUE),
    regexp = "`y` is missing for paired test"
  )
})

test_that("one sample test", {
  data(nhanes)

  expect_no_error(
    with(
      nhanes,
      boltzmann.test(
        x = BMXWT,
        mu = 75
      )
    )
  )
})

test_that("two sample test", {
  data(nhanes)

  expect_no_error(
    with(
      nhanes,
      boltzmann.test(
        x = BMXWT[RIAGENDR == "male"],
        y = BMXWT[RIAGENDR == "female"]
      )
    )
  )
  expect_no_error(
    boltzmann.test(
      x =stats::rnorm(50, mean = 1),
      y = exp(stats::rnorm(50, mean = -0.5))
    )
  )
})

test_that("paired two sample test", {
  x <- stats::rpois(50, lambda = 1)
  y <- stats::rpois(50, lambda = 1)
  expect_no_error(
    boltzmann.test(
      x,y,
      paired = TRUE
    )
  )

})

