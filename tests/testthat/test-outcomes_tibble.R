test_that("`empirical` attribute is missing",{
  object <- data.frame(
    A = rep(c(1,2), 5)
  )
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "`empirical` attribute is missing"
  )

})

test_that("`sampleSize` attribute is missing",{
  object <- data.frame(
    A = rep(c(1,2), 5)
  )
  attr(object, "empirical") = rep(1/10, 10)
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "`sampleSize` attribute is missing"
  )
})

test_that("`empirical` must be numeric",{
  object <- data.frame(
    A = rep(c(1,2), 5)
  )
  attr(object, "empirical") = rep("A", 10)
  attr(object, "sampleSize") = 10
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "empirical` must be numeric"
  )
})


test_that("`empirical` must have one value per row of `data`.",{
  object <- data.frame(
    A = rep(c(1,2), 5)
  )
  attr(object, "empirical") = rep(0.1, 5)
  attr(object, "sampleSize") = 10
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "`empirical` must have one value per row of `data`."
  )
})

test_that("`empirical` contains negative entries", {
  object <- data.frame(
    A = rep(c(1,2), 5)
  )
  attr(object, "empirical") = c(-0.1, rep(0.9 / 9, 9))
  attr(object, "sampleSize") = 10
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "`empirical` contains negative entries"
  )

})

test_that("`empirical` does not sum to 1",{
  object <- data.frame(
    A = rep(c(1,2), 5)
  )
  attr(object, "empirical") = c(0.2, rep(1 / 9, 9))
  attr(object, "sampleSize") = 10
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "`empirical` does not sum to 1"
  )
})
test_that("sampleSize` must be a single number.",{
  object <- data.frame(
    A = rep(c(1,2), 5)
  )
  attr(object, "empirical") = rep(1/10, 10)
  attr(object, "sampleSize") = "A"
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "sampleSize` must be a single number."
  )
  attr(object, "sampleSize") = c(5, 5)
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "sampleSize` must be a single number."
  )

})

test_that("`sampleSize` must be positive.",{
  object <- data.frame(
    A = rep(c(1,2), 5)
  )
  attr(object, "empirical") = rep(1/10, 10)
  attr(object, "sampleSize") = -1
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "`sampleSize` must be positive."
  )
  attr(object, "sampleSize") = 0
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "`sampleSize` must be positive."
  )
})

test_that("fewer observations than outcomes",{
  object <- data.frame(
    A = rep(c(1,2), 5)
  )
  attr(object, "empirical") = rep(1/10, 10)
  attr(object, "sampleSize") = 9
  expect_error(
    validate_outcomes_tibble(object),
    regexp = "fewer observations (9) than outcomes (10)"
  )
})





test_that("input has to be a data.frame", {
  expect_error(outcomes_tibble(list(a = 1, b = 1:10)))
  expect_no_error(outcomes_tibble(data.frame(a = 1, b = 1:10)))
})

test_that("input contains NA and non-finite values", {
  withInvalidData <- data.frame(
    a = c(1, Inf, NaN, 4, 5),
    b = c("A", "A", "B", NA, "B")
  )
  onlyValidData <- data.frame(
    a = c(1, 5),
    b = factor(c("A", "B"))
  )
  expect_warning(outcomes_tibble(withInvalidData))
  expect_no_warning(outcomes_tibble(onlyValidData))
  expect_equal(
    outcomes_tibble(withInvalidData),
    outcomes_tibble(onlyValidData)
  )
})
test_that("dplyr row slicing works on outcomes_tibble",{
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  outcomes <- dplyr::filter(outcomes, RIAGENDR == "male")
  expect_equal(sum(empirical(outcomes)), 1)
  expect_equal(length(empirical(outcomes)), 2551L)
  expect_equal(sampleSize(outcomes), 2551L)
})
test_that("subsetting with [] works on outcomes_tibble", {
  data(nhanes)
  outcomes <- outcomes_tibble(nhanes)
  outcomes <- outcomes[1:100, c("RIAGENDR", "BMXWT")]
  expect_equal(sum(empirical(outcomes)), 1)
  expect_equal(length(empirical(outcomes)), 96L)
  expect_equal(sampleSize(outcomes), 100L)
})
