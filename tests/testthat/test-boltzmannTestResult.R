test_that("missing fields",{
  fields <- c(
    "statistic", "iDivergence", "degreesOfFreedom", "sampleSize",
    "pValue", "alternativeDistribution", "alternativeExpectations",
    "hypothesisDistribution", "hypothesisExpectations", "testedExpectations",
    "dataName", "coefficientMatrix", "method", "ambientDistribution",
    "ambientExpectations"
  )
  for (i in seq_len(length(fields))){
    object <- as.list(
      fields[-i]
    )
    expect_error(
      validate_boltzmannTestResult(object),
      regexp = "missing fields"
    )
  }


})

nestedBoltzmannTest <- function(){
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
  ## perform nested Boltzmann Test
  boltzmannTest(
    outcomes = outcomes,
    G = G,
    eta = eta,
    testedExpectations = 5,
    ambientExpectations = c(3,4)
  )
}

test_that("`statistic` must be a single numeric value",{
  bt <- nestedBoltzmannTest()

  ## statistic is not a number
  bt$statistic <- "A"
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`statistic` must be a single numeric value"
  )
  ## statistic is a vector
  bt$statistic <- c(1,1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`statistic` must be a single numeric value"
  )
  ## statistic is a list
  bt$statistic  <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`statistic` must be a single numeric value"
  )
  ## statistic is a data.frame
  bt$statistic  <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`statistic` must be a single numeric value"
  )

})

test_that("`iDivergence` must be a single numeric value", {
  bt <- nestedBoltzmannTest()

  ## iDivergence is not a number
  bt$iDivergence <- "A"
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`iDivergence` must be a single numeric value"
  )
  ## iDivergence is a vector
  bt$iDivergence <- c(1,1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`iDivergence` must be a single numeric value"
  )
  ## iDivergence is a list
  bt$iDivergence <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`iDivergence` must be a single numeric value"
  )
  ## iDivergence is a data.frame
  bt$iDivergence <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`iDivergence` must be a single numeric value"
  )
})

test_that("`degreesOfFreedom` must be a single integer value",{
  bt <- nestedBoltzmannTest()


  ## degreesOfFreedom is not a number
  bt$degreesOfFreedom <- "A"
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`degreesOfFreedom` must be a single integer value"
  )
  ## degreesOfFreedom is a vector
  bt$degreesOfFreedom <- c(1,1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`degreesOfFreedom` must be a single integer value"
  )
  ## degreesOfFreedom is a list
  bt$degreesOfFreedom <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`degreesOfFreedom` must be a single integer value"
  )
  ## degreesOfFreedom is a data.frame
  bt$degreesOfFreedom <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`degreesOfFreedom` must be a single integer value"
  )
  ## degreesOfFreedom is not an integer
  bt$degreesOfFreedom <- 1.1
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`degreesOfFreedom` must be a single integer value"
  )
})

test_that("`degreesOfFreedom` must be at least 1",{
  bt <- nestedBoltzmannTest()
  bt$degreesOfFreedom <- 0L
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`degreesOfFreedom` must be at least 1"
  )
})

test_that("`sampleSize` must be a single integer value",{
  bt <- nestedBoltzmannTest()


  ## sampleSize is not a number
  bt$sampleSize <- "A"
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`sampleSize` must be a single integer value"
  )
  ## sampleSize is a vector
  bt$sampleSize <- c(1,1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`sampleSize` must be a single integer value"
  )
  ## sampleSize is a list
  bt$sampleSize <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`sampleSize` must be a single integer value"
  )
  ## sampleSize is a data.frame
  bt$sampleSize <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`sampleSize` must be a single integer value"
  )
  ## sampleSize is not an integer
  bt$sampleSize <- 1.1
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`sampleSize` must be a single integer value"
  )
})

test_that("`sampleSize` must be at least 2",{
  bt <- nestedBoltzmannTest()
  bt$sampleSize <- 1L
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`sampleSize` must be at least 2"
  )
})

test_that("`pValue` must be a single numeric value",{
  bt <- nestedBoltzmannTest()

  ## pValue is not a number
  bt$pValue <- "A"
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`pValue` must be a single numeric value"
  )
  ## pValue is a vector
  bt$pValue <- c(1,1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`pValue` must be a single numeric value"
  )
  ## pValue is a list
  bt$pValue <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`pValue` must be a single numeric value"
  )
  ## pValue is a data.frame
  bt$pValue <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`pValue` must be a single numeric value"
  )
})

test_that("`pValue` must be between 0 and 1",{
  bt <- nestedBoltzmannTest()

  ## pValue is smaller than 0
  bt$pValue <- -0.1
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`pValue` must be between 0 and 1"
  )
  ## pValue is larger than 1
  bt$pValue <- 1.1
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`pValue` must be between 0 and 1"
  )
})
test_that("`alternativeDistribution` must be a numeric vector",{
  bt <- nestedBoltzmannTest()
  ## alternativeDistribution is a character
  bt$alternativeDistribution <- rep("A", length(bt$alternativeDistribution))
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`alternativeDistribution` must be a numeric vector"
  )
  ## alternativeDistribution is a list
  bt$alternativeDistribution <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`alternativeDistribution` must be a numeric vector"
  )
  ## alternativeDistribution is a data.frame
  bt$alternativeDistribution <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`alternativeDistribution` must be a numeric vector"
  )
})

test_that("`alternativeExpectations` must be a numeric vector",{
  bt <- nestedBoltzmannTest()
  ## alternativeExpectations is a character
  bt$alternativeExpectations <- rep("A", length(bt$alternativeExpectations))
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`alternativeExpectations` must be a numeric vector"
  )
  ## alternativeExpectations is a list
  bt$alternativeExpectations <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`alternativeExpectations` must be a numeric vector"
  )
  ## alternativeExpectations is a data.frame
  bt$alternativeExpectations <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`alternativeExpectations` must be a numeric vector"
  )
})


test_that("`hypothesisDistribution` must be a numeric vector",{
  bt <- nestedBoltzmannTest()
  ## hypothesisDistribution is a character
  bt$hypothesisDistribution <- rep("A", length(bt$hypothesisDistribution))
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`hypothesisDistribution` must be a numeric vector"
  )
  ## hypothesisDistribution is a list
  bt$hypothesisDistribution <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`hypothesisDistribution` must be a numeric vector"
  )
  ## hypothesisDistribution is a data.frame
  bt$hypothesisDistribution <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`hypothesisDistribution` must be a numeric vector"
  )
})

test_that("`hypothesisExpectations` must be a numeric vector",{
  bt <- nestedBoltzmannTest()
  ## hypothesisExpectations is a character
  bt$hypothesisExpectations <- rep("A", length(bt$hypothesisExpectations))
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`hypothesisExpectations` must be a numeric vector"
  )
  ## hypothesisExpectations is a list
  bt$hypothesisExpectations <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`hypothesisExpectations` must be a numeric vector"
  )
  ## hypothesisExpectations is a data.frame
  bt$hypothesisExpectations <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`hypothesisExpectations` must be a numeric vector"
  )
})

test_that("`ambientDistribution` must be a numeric vector",{
  bt <- nestedBoltzmannTest()
  ## ambientDistribution is a character
  bt$ambientDistribution <- rep("A", length(bt$ambientDistribution))
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`ambientDistribution` must be a numeric vector"
  )
  ## ambientDistribution is a list
  bt$ambientDistribution <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`ambientDistribution` must be a numeric vector"
  )
  ## ambientDistribution is a data.frame
  bt$ambientDistribution <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`ambientDistribution` must be a numeric vector"
  )
})

test_that("`ambientExpectations` must be a numeric vector",{
  bt <- nestedBoltzmannTest()
  ## ambientExpectations is a character
  bt$ambientExpectations <- rep("A", length(bt$ambientExpectations))
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`ambientExpectations` must be a numeric vector"
  )
  ## ambientExpectations is a list
  bt$ambientExpectations <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`ambientExpectations` must be a numeric vector"
  )
  ## ambientExpectations is a data.frame
  bt$ambientExpectations <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`ambientExpectations` must be a numeric vector"
  )
})

test_that(
"`hypothesisDistribution` and `alternativeDistribution` must
have the same length",{
  bt <- nestedBoltzmannTest()
  ## hypothesis Distribution is shorter than alternative Distribution
  bt$hypothesisDistribution = rep(0.1, 5)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "hypothesisDistribution` and `alternativeDistribution` ",
      "must have the same length"
    )
  )
  ## hypothesis Distribution is longer than alternative Distribution
  bt$hypothesisDistribution = rep(0.1, 9)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "hypothesisDistribution` and `alternativeDistribution` ",
      "must have the same length"
    )
  )
})
test_that(
  "`hypothesisExpectations` and `alternativeExpectations` must
have the same length",{
  bt <- nestedBoltzmannTest()
  ## hypothesis Expectations is shorter than alternative Expectations
  bt$hypothesisExpectations = rep(0.1, 4)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "hypothesisExpectations` and `alternativeExpectations` ",
      "must have the same length"
    )
  )
  ## hypothesis Expectations is longer than alternative Expectations
  bt$hypothesisExpectations = rep(0.1, 6)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "hypothesisExpectations` and `alternativeExpectations` ",
      "must have the same length"
    )
  )
})

test_that("`testedExpectations` must be a integer vector",{
  bt <- nestedBoltzmannTest()
  ## testExpectations is a character
  bt$testedExpectations <- "A"
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`testedExpectations` must be a integer vector"
  )
  ## testExpectations is a list
  bt$testedExpectations <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`testedExpectations` must be a integer vector"
  )
  ## testExpectations is a data.frame
  bt$testedExpectations <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`testedExpectations` must be a integer vector"
  )
  ## testExpectations is a numeric
  bt$testedExpectations <- 1.1
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`testedExpectations` must be a integer vector"
  )

})
test_that(
  "`testedExpectations` must be between 1
  and the number of hypothesis expectations",{
  bt <- nestedBoltzmannTest()
  bt$testedExpectations <- 0L
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "`testedExpectations` must be between 1 ",
      "and the number of hypothesis expectations"
    )
  )
  bt$testedExpectations <- 6L
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "`testedExpectations` must be between 1 ",
      "and the number of hypothesis expectations"
    )
  )

})

test_that("`dataName` must be a string",{
  bt <- nestedBoltzmannTest()
  ## dataName is a number
  bt$dataName <- 1
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`dataName` must be a string"
  )
  ## dataName is a list
  bt$dataName <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`dataName` must be a string"
  )
  ## dataName is a data.frame
  bt$dataName <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`dataName` must be a string"
  )
})

test_that("`coefficientMatrix` must be a numeric matrix", {
  bt <- nestedBoltzmannTest()
  ## coefficientMatrix is a character matrix
  bt$coefficientMatrix[1,1] <- "A"
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`coefficientMatrix` must be a numeric matrix"
  )
  ## coefficientMatrix is a single number
  bt$coefficientMatrix <- 1
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`coefficientMatrix` must be a numeric matrix"
  )
  ## coefficientMatrix is a list
  bt$coefficientMatrix <- list(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`coefficientMatrix` must be a numeric matrix"
  )
  ## coefficientMatrix is a data.frame
  bt$coefficientMatrix <- data.frame(1)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = "`coefficientMatrix` must be a numeric matrix"
  )
})

test_that(
  "`hypothesisDistribution` and `ambientDistribution` must
have the same length",{
  bt <- nestedBoltzmannTest()
  ## hypothesis Distribution is shorter than ambient Distribution
  bt$ambientDistribution = rep(0.1, 5)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "hypothesisDistribution` and `ambientDistribution` ",
      "must have the same length"
    )
  )
  ## hypothesis Distribution is longer than ambient Distribution
  bt$ambientDistribution = rep(0.1, 9)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "hypothesisDistribution` and `ambientDistribution` ",
      "must have the same length"
    )
  )
})
test_that(
  "`hypothesisExpectations` and `ambientExpectations` must
have the same length",{
  bt <- nestedBoltzmannTest()
  ## ambient Expectations is shorter than hypothesis Expectations
  bt$ambientExpectations = rep(0.1, 4)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "hypothesisExpectations` and `ambientExpectations` ",
      "must have the same length"
    )
  )
  ## ambient Expectations is longer than hypothesis Expectations
  bt$ambientExpectations = rep(0.1, 6)
  expect_error(
    validate_boltzmannTestResult(bt),
    regexp = paste0(
      "hypothesisExpectations` and `ambientExpectations` ",
      "must have the same length"
    )
  )
})

test_that("print boltzmannTestResult",{
  bt <- nestedBoltzmannTest()

  expect_output(
    print(bt),
    regexp = "Boltzmann Test"
  )

  expect_output(
    print(bt),
    regexp = "data:  outcomes"
  )

  expect_output(
    print(bt),
    regexp = "statistic = 3.559549, df = 1, N = 700, p-value = 0.0592"
  )
  expect_output(
    print(bt),
    regexp =
      "                            hypothesis   ambient alternative tested"
  )
  expect_output(
    print(bt),
    regexp =
      "norm                              1.00 1.0000000  1.00000000       "
  )

  expect_output(
    print(bt),
    regexp =
      "treatment_A                       0.50 0.5000000  0.50000000       "
  )
  expect_output(
    print(bt),
    regexp =
      "stoneSize_small.treatment_B       0.51 0.7714286  0.77142857       "
  )
  expect_output(
    print(bt),
    regexp =
      "stoneSize_small.treatment_A       0.51 0.2485714  0.24857143       "
  )
  expect_output(
    print(bt),
    regexp =
      "success_yes:B_vs_A                0.00 0.1011363  0.04571429      *"
  )


})
