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
