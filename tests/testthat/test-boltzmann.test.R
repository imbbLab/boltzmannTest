## boltzmann.test.outcomes_tibble
test_that("non-numeric or non-matrix G", {
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
      G = as.list(G),
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4)
    ),
    regexp = "`G` is not a numeric matrix"
  )
  expect_error(
    boltzmann.test(
      outcomes = outcomes,
      G = as.numeric(G),
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4)
    ),
    regexp = "`G` is not a numeric matrix"
  )
  expect_error(
    boltzmann.test(
      outcomes = outcomes,
      G = as.character(G),
      eta = eta,
      testedExpectations = 5,
      ambientExpectations = c(3,4)
    ),
    regexp = "`G` is not a numeric matrix"
  )
})


test_that("non-numeric eta", {
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
    stoneSize_small.treatment_B = "A",
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
    regexp = "eta` must be numeric"
  )
})

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
## boltzmann.test.numeric

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

test_that("one sample test with paired = TRUE", {
  x <- stats::rpois(50, lambda = 1)
  expect_error(
    boltzmann.test(
      x,
      paired = TRUE
    ),
    regexp = "`y` is missing for paired test"
  )

})

## boltzmann.test.formula

test_that("variables in the formula are not in the data", {
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = factor(rep(c("A", "B"),5))
  )
  expect_error(
    boltzmann.test(y ~ z, data = data),
    regexp = "all variables in the formula must be columns in data"
  )
})

test_that("variables in formula must not be character vectors", {
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = rep(c("A", "B"),5)
  )
  expect_error(
    boltzmann.test(y ~ x, data = data),
    regexp = "variables in `formula` are character vectors."
  )
})

test_that("warning for non-finite and/or missing data", {
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(Inf, 2)),
    x = factor(rep(c("A", "B"),5))
  )
  expect_warning(
    boltzmann.test(y ~ x, data = data),
    regexp = "selected columns in data contain unobserved values. "
  )
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = factor(rep(c("A", "B"),5))
  )
  data[1,2] <- NA
  expect_warning(
    boltzmann.test(y ~ x, data = data),
    regexp = "selected columns in data contain unobserved values. "
  )
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = factor(rep(c("A", "B"),5))
  )
  data[1,1] <- NaN
  expect_warning(
    boltzmann.test(y ~ x, data = data),
    regexp = "selected columns in data contain unobserved values. "
  )
})

test_that("formula is misspecified", {
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = factor(rep(c("A", "B"),5))
  )
  expect_error(
    boltzmann.test(y ~ x ~ 1, data = data),
    regexp = "formula is misspecified."
  )
})

test_that("targets are misspecified", {
  data <- data.frame(
    z = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = factor(rep(c("A", "B"),5)),
    y = rep(1, 10)
  )
  expect_error(
    boltzmann.test(z + x ~ y, data = data),
    regexp = "targets are misspecified."
  )
})

test_that("lhs is a matrix",{
  data <- data.frame(
    x = stats::rnorm(10),
    y = stats::rnorm(10),
    a = factor(rep(c("A", "B"), 5))
  )
  bt <- expect_no_error(boltzmann.test(cbind(x, y) ~ a, data = data))
  expect_equal(rownames(bt$coefficientMatrix), c("norm","B","x:B_vs_A","y:B_vs_A"))
})

test_that("lhs is only one variable",{
  data(nhanes)
  bt <- expect_no_error(
    boltzmann.test(BMXWT / (BMXHT / 100)^2 ~ RIAGENDR, data = nhanes)
  )
  expect_equal(
    rownames(bt$coefficientMatrix),
    c("norm","female", "BMXWT/(BMXHT/100)^2:female_vs_male")
  )
})

test_that("multifactorial groups", {
  data <- data.frame(
    x = stats::rnorm(100),
    a = factor(rep(c("A", "B"), 50)),
    b = factor(rep(c("1", "2"), each = 50))
  )
  bt <- expect_no_error(boltzmann.test(x ~ a + b, data = data))
  expect_equal(
    rownames(bt$coefficientMatrix),
    c(
      "norm", "A.2", "B.1", "B.2",
      "x:A.2_vs_A.1", "x:B.1_vs_A.1", "x:B.2_vs_A.1"
    )
  )
})

test_that("all grouping variables must be factors",{
  data <- data.frame(
    z = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = factor(rep(c("A", "B"),5)),
    y = rep(1, 10)
  )
  expect_error(
    boltzmann.test(z ~ x + y, data = data),
    regexp = "all grouping variables must be factors."
  )
})

test_that("not enough groups to compare",{
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = factor(rep(c("A"),10))
  )
  expect_error(
    boltzmann.test(y ~ x, data = data),
    regexp = "not enough groups to compare"
  )
})

test_that("categorical target .+ has only one level",{
  data <- data.frame(
    y = factor(rep(c("A"),10)),
    x = factor(rep(c("X", "Y"), each = 5))
  )
  expect_error(
    boltzmann.test(y ~ x, data = data),
    regexp = "categorical target .+ has only one level"
  )
  expect_error(
    boltzmann.test(~ y + x, data = data),
    regexp = "categorical target .+ has only one level"
  )
})

test_that("the number of expectations in `nu` does not match",{
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = factor(rep(c("A", "B"),5))
  )
  expect_error(
    boltzmann.test(y ~ x, data = data, nu = c(0,0)),
    regexp = "the number of expectations in `nu` does not match"
  )
  expect_error(
    boltzmann.test(~ y + x, data = data, nu = c(0,0, 10)),
    regexp = "the number of expectations in `nu` does not match"
  )
})

test_that("given expectations are tested", {
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    a = factor(rep(c("A", "B"), 5))
  )
  nu = 0.4
  names(nu) = "y:B_vs_A"
  bt <- expect_no_error(
    boltzmann.test(y ~ a, data = data, nu = nu)
  )
  expect_equal(bt$hypothesisExpectations[3], nu)
})

test_that("hypothesized prevalences outside (0, 1)",{
  data <- data.frame(
    y = c(rep(-1, 4), rep(0, 4), rep(1, 2)),
    x = factor(rep(c("A", "B"),5))
  )
  expect_error(
    boltzmann.test(~ y + x, data = data),
    regexp = "hypothesized prevalences outside \\(0, 1\\)"
  )
})

