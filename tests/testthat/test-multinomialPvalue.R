test_that("exact multininomial p-value", {
  data <- data.frame(
    x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
  )
  bt <- boltzmann.test(~x, data = data, nu = 0.5)
  outcomes <- outcomes_tibble(data)
  exactPvalue <- multinomialPvalue(
    G = bt$coefficientMatrix[-1, , drop = FALSE],
    eta = bt$hypothesisExpectations[-1],
    mu = bt$alternativeExpectations[-1],
    v = bt$hypothesisDistribution,
    N = bt$sampleSize
  )

  ## get all possible data sets
  dataSets <- starsAndBars(
    N = bt$sampleSize,
    k = nrow(outcomes)
  )
  ## get multinomial probabilities for the hypothesis distribution
  dmult <- apply(
    dataSets, 2,
    function(x)
      stats::dmultinom(x, prob = bt$hypothesisDistribution)
  )
  ## get the average for x for each data sets
  sampleMeans <- apply(
    dataSets, 2,
    function(x){
      sum(outcomes$x * x) / sum(x)
    }
  )
  ## summarize probabilities per sample mean
  dSampleMean <- tapply(
    dmult, sampleMeans, sum
  )

  expectedExactPvalue <- sum(dSampleMean[dSampleMean <= dSampleMean["-0.2"]])
  expect_equal(exactPvalue, expectedExactPvalue)

})

test_that("`G must be a numeric matrix",{
  data <- data.frame(
    x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
  )
  bt <- boltzmann.test(~x, data = data, nu = 0.5)
  G <- bt$coefficientMatrix[-1, , drop = FALSE]
  expect_no_error(
    exactPvalue <- multinomialPvalue(
      G = G,
      eta = bt$hypothesisExpectations[-1],
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    )
  )
  ## G is a character
  G[1, 1] <- "A"
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = G,
      eta = bt$hypothesisExpectations[-1],
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`G must be a numeric matrix"
  )
  ## G is a list
  G <- list(1)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = G,
      eta = bt$hypothesisExpectations[-1],
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`G must be a numeric matrix"
  )
  ## G is a data.frame
  G <- data.frame(1)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = G,
      eta = bt$hypothesisExpectations[-1],
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`G must be a numeric matrix"
  )

})

test_that("`eta` must be a numeric vector",{
  data <- data.frame(
    x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
  )
  bt <- boltzmann.test(~x, data = data, nu = 0.5)
  eta <- bt$hypothesisExpectations[-1]
  expect_no_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = eta,
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    )
  )
  ## eta is a character
  eta <- "A"
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = eta,
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`eta` must be a numeric vector"
  )
  ## eta is a list
  eta <- list(1)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = eta,
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`eta` must be a numeric vector"
  )
  ## eta is a data.frame
  eta <- data.frame(1)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = eta,
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`eta` must be a numeric vector"
  )

})

test_that("`mu` must be a numeric vector",{
  data <- data.frame(
    x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
  )
  bt <- boltzmann.test(~x, data = data, nu = 0.5)
  mu <- bt$alternativeExpectations[-1]
  expect_no_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = bt$hypothesisExpectations[-1],
      mu = mu,
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    )
  )
  ## mu is a character
  mu <- "A"
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = bt$hypothesisExpectations[-1],
      mu = mu,
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`mu` must be a numeric vector"
  )
  ## eta is a list
  mu <- list(1)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = bt$hypothesisExpectations[-1],
      mu = mu,
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`mu` must be a numeric vector"
  )
  ## eta is a data.frame
  mu<- data.frame(1)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = bt$hypothesisExpectations[-1],
      mu = mu,
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`mu` must be a numeric vector"
  )

})

test_that("`v` must be a numeric vector", {
  data <- data.frame(
    x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
  )
  bt <- boltzmann.test(~x, data = data, nu = 0.5)
  v <- bt$hypothesisDistribution

  ## v is a character
  v <- "A"
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = bt$hypothesisExpectations[-1],
      mu = bt$alternativeExpectations[-1],
      v = v,
      N = bt$sampleSize
    ),
    regexp = "`v` must be a numeric vector"
  )
  ## v is a list
  v <- list(1)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = bt$hypothesisExpectations[-1],
      mu = bt$alternativeExpectations[-1],
      v = v,
      N = bt$sampleSize
    ),
    regexp = "`v` must be a numeric vector"
  )
  ## v is a data.frame
  v <- data.frame(1)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = bt$hypothesisExpectations[-1],
      mu = bt$alternativeExpectations[-1],
      v = v,
      N = bt$sampleSize
    ),
    regexp = "`v` must be a numeric vector"
  )
})

test_that("`eta` must have one value per row of `G`",{
  data <- data.frame(
    x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
  )
  bt <- boltzmann.test(~x, data = data, nu = 0.5)
  ## one more
  eta <- c(1,2)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = eta,
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`eta` must have one value per row of `G`"
  )
  ## one less
  eta <- bt$hypothesisExpectations[-1]
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix,
      eta = eta,
      mu = bt$alternativeExpectations[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`eta` must have one value per row of `G`"
  )

})

test_that("`eta` and `mu` must be the same length",{
  data <- data.frame(
    x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
  )
  bt <- boltzmann.test(~x, data = data, nu = 0.5)
  ## one more
  mu <- c(1,2)
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = bt$hypothesisExpectations[-1],
      mu = mu,
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`eta` and `mu` must be the same length"
  )

})

test_that("`v` must have one value per column of `G`",{
  data <- data.frame(
    x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
  )
  bt <- boltzmann.test(~x, data = data, nu = 0.5)
  ## one less
  v <- bt$hypothesisDistribution[-1]
  expect_error(
    exactPvalue <- multinomialPvalue(
      G = bt$coefficientMatrix[-1, , drop = FALSE],
      eta = bt$hypothesisExpectations[-1],
      mu = bt$alternativeDistribution[-1],
      v = bt$hypothesisDistribution,
      N = bt$sampleSize
    ),
    regexp = "`v` must have one value per column of `G`"
  )
})
